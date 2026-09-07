#!/usr/bin/env python3
"""Read schema-v1 RIXS state-character shards and analyze local multipoles.

The persisted density is the unnormalized physical ket-bra density in the
native representation,

    rho(m,s,m',s') = C(m,s) * conj(C(m',s')),

with m=-2,...,+2 and interleaved spin order s=1 (+1/2), s=2 (-1/2).
The HDF matrix row/column index is ``2*(m+2) + (s-1)`` in zero-based Python
indexing.  No transpose or conjugation is applied when reading it.

Structural multipoles are formed only after applying the stored common-frame
map ``T = kron(orbital_global_to_local, spin_native_to_local)`` as
``T @ rho_native @ T.conj().T``.  The Python transform is a post-processing
mirror of FLEUR's canonical Fortran Nordstrom/Bultmark implementation.
"""

from __future__ import annotations

from dataclasses import dataclass
from math import factorial, sqrt
from pathlib import Path
from typing import Iterable, Mapping, Sequence

import h5py
import numpy as np


SCHEMA_VERSION = 1
RHO_CONVENTION = "unnormalized native local MT d-spin reduced density matrix"
NATIVE_BASIS = "complex Y_2m m=-2..2, each with native-MT spin components 1,2"
IDENTITY_COLUMNS = "ikpt band physical_atom atom_type iatom_l role_mask"
SCALAR_COLUMNS = "k_weight energy_Ha occupation d_weight t2g_weight eg_weight t2g_check"
L = 2
T_OFFSET = 2 * L + 1
Channel = tuple[int, int, int]
StateKey = tuple[int, int, int]
# Schema-v1 stores the first-variation spinor occupation used by RIXS's
# production empty-state test as 1.0-occupation; one is therefore full.
SCHEMA_V1_MAX_OCCUPATION = 1.0


def valid_components() -> tuple[tuple[int, int, int, int], ...]:
    """Return the 100 valid l=2 components in canonical (k,p,r,t) order."""
    return tuple(
        (k, p, r, t)
        for k in range(2 * L + 1)
        for p in range(2)
        for r in range(abs(k - p), k + p + 1)
        for t in range(-r, r + 1)
    )


VALID_COMPONENTS = valid_components()
if len(VALID_COMPONENTS) != 100:
    raise RuntimeError("internal l=2 multipole component enumeration is invalid")


def valid_channels() -> tuple[Channel, ...]:
    """Return the 18 valid l=2 channels in canonical (k,p,r) order."""
    return tuple(dict.fromkeys(component[:3] for component in VALID_COMPONENTS))


VALID_CHANNELS = valid_channels()
if len(VALID_CHANNELS) != 18 or sum(2 * r + 1 for _, _, r in VALID_CHANNELS) != 100:
    raise RuntimeError("internal l=2 multipole channel enumeration is invalid")

_VALID_COMPONENT_SET = frozenset(VALID_COMPONENTS)
_VALID_CHANNEL_SET = frozenset(VALID_CHANNELS)
_CHANNEL_ORDER = {channel: index for index, channel in enumerate(VALID_CHANNELS)}


@dataclass(frozen=True)
class DominantChannel:
    """Largest and next-largest finite character values in a candidate set."""

    channel: Channel
    value: float
    second_channel: Channel | None
    second_value: float | None
    gap: float | None


@dataclass(frozen=True)
class SiteCharacter:
    """Site-level frame data shared by all state records for one physical atom."""

    physical_atom: int
    atom_type: int
    iatom_l: int
    atomic_number: int
    orbital_global_to_local: np.ndarray
    spin_native_to_local: np.ndarray
    datasets: Mapping[str, np.ndarray]
    source_shards: tuple[Path, ...]

    @property
    def combined_native_to_local(self) -> np.ndarray:
        """Return T in interleaved (m,spin) ordering."""
        return np.kron(self.orbital_global_to_local, self.spin_native_to_local)


@dataclass(frozen=True)
class StateCharacter:
    """One unweighted, unnormalized schema-v1 band/site state record."""

    ikpt: int
    band: int
    physical_atom: int
    atom_type: int
    iatom_l: int
    role_mask: int
    k_vector: np.ndarray
    k_weight: float
    band_energy: float
    occupation: float
    rho_native: np.ndarray
    source_shard: Path

    @property
    def key(self) -> tuple[int, int, int]:
        return self.ikpt, self.band, self.physical_atom


@dataclass(frozen=True)
class StateCharacterData:
    """Merged records and consistent site metadata from one or more shards."""

    records: tuple[StateCharacter, ...]
    sites: Mapping[int, SiteCharacter]

    def select(
        self,
        *,
        physical_atom: int | None = None,
        ikpt: int | None = None,
        band: int | None = None,
        role_mask: int | None = None,
    ) -> tuple[StateCharacter, ...]:
        """Select records without applying weights, occupations, or normalization.

        When ``role_mask`` is supplied, every requested role bit must be present.
        """
        selected = []
        for record in self.records:
            if physical_atom is not None and record.physical_atom != physical_atom:
                continue
            if ikpt is not None and record.ikpt != ikpt:
                continue
            if band is not None and record.band != band:
                continue
            if role_mask is not None and (record.role_mask & role_mask) != role_mask:
                continue
            selected.append(record)
        return tuple(selected)

    def structural_density(self, record: StateCharacter) -> np.ndarray:
        """Transform one physical native density into its stored structural frame."""
        try:
            site = self.sites[record.physical_atom]
        except KeyError as exc:
            raise ValueError(
                f"state {record.key} has no matching site metadata"
            ) from exc
        if record.atom_type != site.atom_type or record.iatom_l != site.iatom_l:
            raise ValueError(f"state/site identity mismatch for state {record.key}")
        transform = site.combined_native_to_local
        return transform @ record.rho_native @ transform.conj().T

    def multipoles(self, record: StateCharacter) -> np.ndarray:
        """Return the complete dense structural-frame l=2 tensor for one state."""
        return density_to_multipoles(self.structural_density(record))

    def sum_structural_density(
        self,
        records: Iterable[StateCharacter] | None = None,
        *,
        weights: Iterable[float] | None = None,
    ) -> np.ndarray:
        """Sum selected structural densities with only explicit real weights.

        Unit weights are used when ``weights`` is omitted.  Neither occupations
        nor k-point weights stored in the records are applied implicitly.
        """
        chosen = tuple(self.records if records is None else records)
        coefficients = _explicit_weights(weights, len(chosen))
        total = np.zeros((10, 10), dtype=np.complex128)
        for coefficient, record in zip(coefficients, chosen, strict=True):
            total += coefficient * self.structural_density(record)
        return total

    def sum_multipoles(
        self,
        records: Iterable[StateCharacter] | None = None,
        *,
        weights: Iterable[float] | None = None,
    ) -> np.ndarray:
        """Transform a density sum; this is not a sum of quadratic strengths."""
        return density_to_multipoles(
            self.sum_structural_density(records, weights=weights)
        )


@dataclass(frozen=True)
class KResolvedSelection:
    """Inclusive energy-window and/or explicit-band state selection.

    When both forms are supplied, a record must satisfy both.  Energies are
    selected by ``energy_min <= band_energy-energy_reference <= energy_max``;
    all three values use the stored schema-v1 energy unit (Hartree).  The
    reference is caller supplied and has no implicit Fermi-level meaning.
    """

    energy_min: float | None = None
    energy_max: float | None = None
    energy_reference: float | None = None
    bands: tuple[int, ...] | None = None

    def __post_init__(self) -> None:
        energy_values = (self.energy_min, self.energy_max, self.energy_reference)
        has_energy_selection = any(value is not None for value in energy_values)
        if has_energy_selection and not all(
            value is not None for value in energy_values
        ):
            raise ValueError(
                "energy_min, energy_max, and energy_reference must be supplied together"
            )
        if has_energy_selection:
            assert self.energy_min is not None
            assert self.energy_max is not None
            assert self.energy_reference is not None
            energy_min = float(self.energy_min)
            energy_max = float(self.energy_max)
            energy_reference = float(self.energy_reference)
            if not np.all(np.isfinite((energy_min, energy_max, energy_reference))):
                raise ValueError("energy selection values must be finite")
            if energy_min > energy_max:
                raise ValueError("energy_min must not exceed energy_max")
            object.__setattr__(self, "energy_min", energy_min)
            object.__setattr__(self, "energy_max", energy_max)
            object.__setattr__(self, "energy_reference", energy_reference)

        if self.bands is not None:
            supplied_bands = tuple(self.bands)
            if any(
                isinstance(band, bool) or not isinstance(band, (int, np.integer))
                for band in supplied_bands
            ):
                raise ValueError("bands must contain integer physical band numbers")
            canonical_bands = tuple(sorted({int(band) for band in supplied_bands}))
            if not canonical_bands or any(band < 1 for band in canonical_bands):
                raise ValueError("bands must contain positive physical band numbers")
            object.__setattr__(self, "bands", canonical_bands)

        if not has_energy_selection and self.bands is None:
            raise ValueError(
                "an energy window and/or explicit band selection is required"
            )

    def matches(self, record: StateCharacter) -> bool:
        """Return whether one record satisfies this selection."""
        if self.bands is not None and record.band not in self.bands:
            return False
        if self.energy_min is not None:
            assert self.energy_reference is not None
            assert self.energy_max is not None
            relative_energy = record.band_energy - self.energy_reference
            if relative_energy < self.energy_min or relative_energy > self.energy_max:
                return False
        return True


@dataclass(frozen=True)
class KResolvedPoint:
    """Density-first local multipole analysis for one sampled k point."""

    ikpt: int
    k_vector: np.ndarray
    k_weight: float
    selected_bands: tuple[int, ...]
    selected_energies: tuple[float, ...]
    selected_occupations: tuple[float, ...]
    applied_weights: tuple[float, ...]
    d_weight: float
    trace_imaginary: float
    rho_hs_norm: float
    rho_structural: np.ndarray
    multipoles: np.ndarray
    raw_strengths: Mapping[Channel, float]
    hs_weights: Mapping[Channel, float]
    character: Mapping[Channel, float]
    norm_fractions: Mapping[Channel, float]
    dominant: DominantChannel | None
    normalized_masked: bool
    local_purity: float

    @property
    def n_selected(self) -> int:
        return len(self.selected_bands)


@dataclass(frozen=True)
class KResolvedAnalysis:
    """Ordered map-ready results for one explicit physical atom."""

    physical_atom: int
    atom_type: int
    iatom_l: int
    selection: KResolvedSelection
    weighting_mode: str
    min_weight: float | None
    candidate_channels: tuple[Channel, ...]
    exclude_scalar: bool
    include_empty: bool
    points: tuple[KResolvedPoint, ...]

    def rows(
        self, *, channels: Iterable[Channel] | None = None
    ) -> tuple[dict[str, object], ...]:
        """Return scalar map rows, optionally adding four columns per channel."""
        requested = _ordered_channel_subset(channels, allow_empty=True)
        output = []
        for point in self.points:
            dominant = point.dominant
            row: dict[str, object] = {
                "ikpt": point.ikpt,
                "kx": float(point.k_vector[0]),
                "ky": float(point.k_vector[1]),
                "kz": float(point.k_vector[2]),
                "k_weight": point.k_weight,
                "n_selected": point.n_selected,
                "selected_bands": point.selected_bands,
                "selected_energies": point.selected_energies,
                "selected_occupations": point.selected_occupations,
                "applied_weights": point.applied_weights,
                "D": point.d_weight,
                "trace_imaginary": point.trace_imaginary,
                "rho_hs_norm": point.rho_hs_norm,
                "local_purity": point.local_purity,
                "normalized_masked": point.normalized_masked,
                "dominant_k": None if dominant is None else dominant.channel[0],
                "dominant_p": None if dominant is None else dominant.channel[1],
                "dominant_r": None if dominant is None else dominant.channel[2],
                "dominant_value": float("nan") if dominant is None else dominant.value,
                "second_k": None
                if dominant is None or dominant.second_channel is None
                else dominant.second_channel[0],
                "second_p": None
                if dominant is None or dominant.second_channel is None
                else dominant.second_channel[1],
                "second_r": None
                if dominant is None or dominant.second_channel is None
                else dominant.second_channel[2],
                "second_value": float("nan")
                if dominant is None or dominant.second_value is None
                else dominant.second_value,
                "dominance_gap": float("nan")
                if dominant is None or dominant.gap is None
                else dominant.gap,
            }
            for channel in requested:
                label = "".join(str(index) for index in channel)
                row[f"S_raw_{label}"] = point.raw_strengths[channel]
                row[f"Q_{label}"] = point.hs_weights[channel]
                row[f"C_orth_{label}"] = point.character[channel]
                row[f"F_{label}"] = point.norm_fractions[channel]
            output.append(row)
        return tuple(output)


def analyze_k_resolved(
    data: StateCharacterData,
    *,
    physical_atom: int,
    selection: KResolvedSelection,
    weighting_mode: str,
    custom_weights: Mapping[StateKey, float] | None = None,
    min_weight: float | None = None,
    candidate_channels: Iterable[Channel] | None = None,
    exclude_scalar: bool = False,
    include_empty: bool = False,
    metadata_tolerance: float = 1.0e-13,
    trace_tolerance: float = 1.0e-12,
) -> KResolvedAnalysis:
    """Analyze density-first local multipoles at each observed k point.

    The selected structural densities are summed before the multipole transform.
    This manifold quantity is not an incoherent sum of per-state strengths.
    ``k_weight`` is retained only as metadata.  Empty points can only be included
    for k-point indices represented by at least one schema-v1 record for the atom.
    """
    atom = _physical_atom(physical_atom)
    if not isinstance(selection, KResolvedSelection):
        raise TypeError("selection must be a KResolvedSelection")
    mode = _weighting_mode(weighting_mode)
    tolerance = _nonnegative_finite(metadata_tolerance, "metadata_tolerance")
    trace_limit = _nonnegative_finite(trace_tolerance, "trace_tolerance")
    threshold = (
        None if min_weight is None else _nonnegative_finite(min_weight, "min_weight")
    )
    candidates = _candidate_channels(candidate_channels, exclude_scalar)

    try:
        site = data.sites[atom]
    except KeyError as exc:
        raise ValueError(f"physical atom {atom} has no site metadata") from exc
    atom_records = tuple(
        record for record in data.records if record.physical_atom == atom
    )
    if not atom_records:
        raise ValueError(f"physical atom {atom} has no state records")

    groups: dict[int, list[StateCharacter]] = {}
    for record in atom_records:
        groups.setdefault(record.ikpt, []).append(record)
    ordered_groups = []
    for ikpt in sorted(groups):
        records = tuple(sorted(groups[ikpt], key=lambda record: record.band))
        _validate_k_group(records, tolerance)
        ordered_groups.append((ikpt, records))

    selected_by_k = {
        ikpt: tuple(record for record in records if selection.matches(record))
        for ikpt, records in ordered_groups
    }
    selected_records = tuple(
        record for ikpt, _ in ordered_groups for record in selected_by_k[ikpt]
    )
    custom = _validated_custom_weights(custom_weights, selected_records, mode)

    points = []
    for ikpt, records in ordered_groups:
        selected = selected_by_k[ikpt]
        if not selected and not include_empty:
            continue
        coefficients = _k_resolved_weights(selected, mode, custom)
        density = data.sum_structural_density(selected, weights=coefficients)
        points.append(
            _analyze_k_point(
                ikpt,
                records[0].k_vector,
                records[0].k_weight,
                selected,
                coefficients,
                density,
                threshold,
                candidates,
                trace_limit,
            )
        )

    return KResolvedAnalysis(
        physical_atom=atom,
        atom_type=site.atom_type,
        iatom_l=site.iatom_l,
        selection=selection,
        weighting_mode=mode,
        min_weight=threshold,
        candidate_channels=candidates,
        exclude_scalar=exclude_scalar,
        include_empty=include_empty,
        points=tuple(points),
    )


def read_state_character_shards(
    paths: str | Path | Sequence[str | Path],
) -> StateCharacterData:
    """Read and merge schema-v1 state-character shards without modifying them.

    Duplicate ``(ikpt,band,physical_atom)`` keys are rejected.  Site metadata
    repeated across shards must agree to double-precision roundoff.
    """
    if isinstance(paths, (str, Path)):
        shard_paths = (Path(paths),)
    else:
        shard_paths = tuple(Path(path) for path in paths)
    if not shard_paths:
        raise ValueError("at least one state-character shard is required")

    records: list[StateCharacter] = []
    keys: set[tuple[int, int, int]] = set()
    sites: dict[int, SiteCharacter] = {}
    ligand_atomic_number: int | None = None

    for path in shard_paths:
        with h5py.File(path, "r") as handle:
            _validate_root(handle, path)
            current_ligand_z = _scalar_integer_attribute(handle, "ligand_atomic_number")
            if ligand_atomic_number is None:
                ligand_atomic_number = current_ligand_z
            elif current_ligand_z != ligand_atomic_number:
                raise ValueError(f"inconsistent ligand atomic number in {path}")

            shard_sites = _read_sites(handle["sites"], path)
            for physical_atom, site in shard_sites.items():
                if physical_atom in sites:
                    sites[physical_atom] = _merge_site(sites[physical_atom], site)
                else:
                    sites[physical_atom] = site

            for record in _read_states(handle["states"], path):
                if record.key in keys:
                    raise ValueError(f"duplicate state key {record.key} in {path}")
                keys.add(record.key)
                records.append(record)

    for record in records:
        if record.physical_atom not in sites:
            raise ValueError(f"state {record.key} has no matching site metadata")
        site = sites[record.physical_atom]
        if record.atom_type != site.atom_type or record.iatom_l != site.iatom_l:
            raise ValueError(f"state/site identity mismatch for state {record.key}")

    records.sort(key=lambda record: record.key)
    return StateCharacterData(tuple(records), sites)


def unflatten_density(rho_flat: np.ndarray) -> np.ndarray:
    """Return rho[m+2,s-1,m'+2,s'-1] from an interleaved 10x10 matrix."""
    density = np.asarray(rho_flat, dtype=np.complex128)
    if density.shape != (10, 10):
        raise ValueError("physical d-spin density must have shape (10,10)")
    return density.reshape(5, 2, 5, 2)


def flatten_density(rho: np.ndarray) -> np.ndarray:
    """Flatten a (5,2,5,2) density using the schema-v1 interleaved ordering."""
    density = np.asarray(rho, dtype=np.complex128)
    if density.shape != (5, 2, 5, 2):
        raise ValueError("physical d-spin density must have shape (5,2,5,2)")
    return density.reshape(10, 10)


def density_to_multipoles(rho_flat: np.ndarray) -> np.ndarray:
    """Compute the canonical complex Nordstrom/Bultmark l=2 tensor.

    The input is the unnormalized physical ket-bra density in schema-v1
    interleaved ordering.  The returned dense array has shape ``(5,2,6,11)``;
    component ``(k,p,r,t)`` is stored at ``w[k,p,r,t+5]``.  Invalid entries
    are exactly zero.
    """
    rho = unflatten_density(rho_flat)
    chi = _spherical_pauli_matrices()
    wxy = np.zeros((5, 2, 9, 3), dtype=np.complex128)
    result = np.zeros((5, 2, 6, 11), dtype=np.complex128)

    for k in range(5):
        normalization = _orbital_normalization(k)
        for p in range(2):
            for x in range(-k, k + 1):
                for y in range(-p, p + 1):
                    value = 0.0j
                    for ma in range(-L, L + 1):
                        for mb in range(-L, L + 1):
                            orbital_factor = (
                                _phase(L - mb)
                                * _integer_wigner_3j(L, k, L, -mb, x, ma)
                                / normalization
                            )
                            for sa in range(2):
                                for sb in range(2):
                                    value += (
                                        orbital_factor
                                        * chi[sb, sa, p, y + 1]
                                        * rho[ma + L, sa, mb + L, sb]
                                    )
                    wxy[k, p, x + 4, y + 1] = value

    for k in range(5):
        for p in range(2):
            for r in range(abs(k - p), k + p + 1):
                bracket = _generalized_three_bracket(k, p, r)
                for t in range(-r, r + 1):
                    value = 0.0j
                    for x in range(-k, k + 1):
                        for y in range(-p, p + 1):
                            value += (
                                wxy[k, p, x + 4, y + 1]
                                * _integer_wigner_3j(k, r, p, -x, t, -y)
                                * _phase(-x - y)
                            )
                    result[k, p, r, t + T_OFFSET] = _phase(k + p) * value / bracket
    return result


def multipole_component(w: np.ndarray, k: int, p: int, r: int, t: int) -> complex:
    """Return one valid component from the dense t-offset representation."""
    tensor = _dense_multipoles(w)
    if (k, p, r, t) not in _VALID_COMPONENT_SET:
        raise ValueError(f"invalid l=2 multipole component {(k, p, r, t)}")
    return complex(tensor[k, p, r, t + T_OFFSET])


def compact_multipoles(w: np.ndarray) -> dict[tuple[int, int, int, int], complex]:
    """Return the 100 valid components as a canonical-order mapping."""
    return {
        component: multipole_component(w, *component) for component in VALID_COMPONENTS
    }


def channel_raw_strengths(w: np.ndarray) -> dict[Channel, float]:
    """Return Nordstrom/Bultmark raw strengths without normalization.

    For each channel this is ``sum_t abs(w_t)**2``.  Structurally invalid
    entries in the dense tensor are ignored.  Raw strengths are rotationally
    invariant within a channel but are not norm fractions across channels.
    """
    tensor = _dense_multipoles(w)
    return {
        (k, p, r): float(
            np.sum(np.abs(tensor[k, p, r, T_OFFSET - r : T_OFFSET + r + 1]) ** 2)
        )
        for k, p, r in VALID_CHANNELS
    }


def channel_metric_factor(k: int, p: int, r: int) -> float:
    """Return A_(kpr) converting raw strength to Hilbert-Schmidt weight."""
    _validate_channel(k, p, r)
    normalization = _orbital_normalization(k)
    bracket = _generalized_three_bracket(k, p, r)
    return float((2 * k + 1) * normalization**2 * (2 * r + 1) * abs(bracket) ** 2 / 2.0)


def channel_hilbert_schmidt_weights(w: np.ndarray) -> dict[Channel, float]:
    """Return the orthogonal Hilbert-Schmidt density weight of each channel."""
    raw = channel_raw_strengths(w)
    return {
        channel: channel_metric_factor(*channel) * value
        for channel, value in raw.items()
    }


def channel_character(
    w: np.ndarray, d_weight: float, *, min_weight: float | None = None
) -> dict[Channel, float]:
    """Return Q_(kpr)/d_weight**2, or NaN values when the weight is unavailable.

    ``min_weight`` is an optional caller-selected scientific threshold.  With
    no threshold, only nonpositive d weight is masked.  No density, occupation,
    or k-point normalization is applied implicitly.
    """
    weights = channel_hilbert_schmidt_weights(w)
    denominator = _normalization_denominator(d_weight, min_weight, "d_weight")
    if denominator is None:
        return _unavailable_channels()
    return {channel: value / denominator**2 for channel, value in weights.items()}


def channel_norm_fractions(
    w: np.ndarray, rho_hs_norm: float, *, min_norm: float | None = None
) -> dict[Channel, float]:
    """Return Q_(kpr)/Tr(rho^dagger rho), or NaN values if unavailable."""
    weights = channel_hilbert_schmidt_weights(w)
    denominator = _normalization_denominator(rho_hs_norm, min_norm, "rho_hs_norm")
    if denominator is None:
        return _unavailable_channels()
    return {channel: value / denominator for channel, value in weights.items()}


def dominant_channel(
    character: Mapping[Channel, float],
    *,
    channels: Iterable[Channel] | None = None,
    exclude_scalar: bool = False,
) -> DominantChannel | None:
    """Rank C_orth or Q values using canonical-order tie breaking.

    All 18 channels are candidates unless ``channels`` is supplied.  Setting
    ``exclude_scalar`` removes only ``(0,0,0)``; no other channel is silently
    excluded.  Raw strengths are not cross-channel comparable and therefore
    are not valid input by contract.  Fully masked candidate sets return
    ``None``.
    """
    candidates = _candidate_channels(channels, exclude_scalar)
    missing = [channel for channel in candidates if channel not in character]
    if missing:
        raise ValueError(f"character mapping is missing channel {missing[0]}")
    finite = [
        channel for channel in candidates if np.isfinite(float(character[channel]))
    ]
    if not finite:
        return None
    ranked = sorted(
        finite,
        key=lambda channel: (-float(character[channel]), _CHANNEL_ORDER[channel]),
    )
    first = ranked[0]
    if len(ranked) == 1:
        return DominantChannel(first, float(character[first]), None, None, None)
    second = ranked[1]
    first_value = float(character[first])
    second_value = float(character[second])
    return DominantChannel(
        first, first_value, second, second_value, first_value - second_value
    )


def channel_family(k: int, p: int, r: int) -> str:
    """Return the broad parity/spin family of a valid l=2 channel."""
    _validate_channel(k, p, r)
    if p == 0:
        return "charge" if k % 2 == 0 else "current/orbital-current"
    return "magnetization" if k % 2 == 0 else "spin-current/spin-orbital"


def channel_subset_sum(
    values: Mapping[Channel, float], channels: Iterable[Channel]
) -> float:
    """Sum an explicit canonical channel subset.

    This is directly meaningful for orthogonal ``Q``, ``F``, and ``C_orth``
    values.  A sum of differently normalized raw strengths must not be labeled
    as a Hilbert-Schmidt family weight.
    """
    requested = _ordered_channel_subset(channels, allow_empty=False)
    missing = [channel for channel in requested if channel not in values]
    if missing:
        raise ValueError(f"channel mapping is missing channel {missing[0]}")
    return float(sum(float(values[channel]) for channel in requested))


def _dense_multipoles(w: np.ndarray) -> np.ndarray:
    tensor = np.asarray(w)
    if tensor.shape != (5, 2, 6, 11):
        raise ValueError("dense l=2 multipole tensor must have shape (5,2,6,11)")
    return tensor


def _validate_channel(k: int, p: int, r: int) -> None:
    if (k, p, r) not in _VALID_CHANNEL_SET:
        raise ValueError(f"invalid l=2 multipole channel {(k, p, r)}")


def _normalization_denominator(
    value: float, minimum: float | None, name: str
) -> float | None:
    denominator = float(value)
    if not np.isfinite(denominator):
        raise ValueError(f"{name} must be finite")
    threshold = 0.0 if minimum is None else float(minimum)
    if not np.isfinite(threshold) or threshold < 0.0:
        raise ValueError(f"minimum {name} must be finite and nonnegative")
    return denominator if denominator > threshold else None


def _unavailable_channels() -> dict[Channel, float]:
    return {channel: float("nan") for channel in VALID_CHANNELS}


def _candidate_channels(
    channels: Iterable[Channel] | None, exclude_scalar: bool
) -> tuple[Channel, ...]:
    if channels is None:
        requested = set(VALID_CHANNELS)
    else:
        requested = set()
        for channel in channels:
            if len(channel) != 3:
                raise ValueError(f"invalid l=2 multipole channel {channel}")
            canonical = (int(channel[0]), int(channel[1]), int(channel[2]))
            _validate_channel(*canonical)
            requested.add(canonical)
    if exclude_scalar:
        requested.discard((0, 0, 0))
    candidates = tuple(channel for channel in VALID_CHANNELS if channel in requested)
    if not candidates:
        raise ValueError("dominant-channel candidate set is empty")
    return candidates


def _explicit_weights(weights: Iterable[float] | None, count: int) -> np.ndarray:
    if weights is None:
        return np.ones(count, dtype=np.float64)
    supplied = np.asarray(tuple(weights))
    if supplied.ndim != 1 or supplied.size != count:
        raise ValueError(f"weights must contain exactly {count} values")
    if np.iscomplexobj(supplied):
        raise ValueError("density weights must be real")
    try:
        coefficients = supplied.astype(np.float64)
    except (TypeError, ValueError) as exc:
        raise ValueError("density weights must be real numbers") from exc
    if not np.all(np.isfinite(coefficients)):
        raise ValueError("density weights must be finite")
    return coefficients


def _physical_atom(value: int) -> int:
    if isinstance(value, bool) or not isinstance(value, (int, np.integer)):
        raise ValueError("physical_atom must be one explicit integer atom index")
    atom = int(value)
    if atom < 1:
        raise ValueError("physical_atom must be positive")
    return atom


def _weighting_mode(value: str) -> str:
    if value not in ("unit", "occupation", "hole", "custom"):
        raise ValueError("weighting_mode must be unit, occupation, hole, or custom")
    return value


def _nonnegative_finite(value: float, name: str) -> float:
    result = float(value)
    if not np.isfinite(result) or result < 0.0:
        raise ValueError(f"{name} must be finite and nonnegative")
    return result


def _validate_k_group(records: tuple[StateCharacter, ...], tolerance: float) -> None:
    if not records:
        raise ValueError("internal empty k-point group")
    reference = records[0]
    if reference.k_vector.shape != (3,) or not np.all(np.isfinite(reference.k_vector)):
        raise ValueError(f"invalid k-vector metadata for ikpt {reference.ikpt}")
    if not np.isfinite(reference.k_weight):
        raise ValueError(f"invalid k-weight metadata for ikpt {reference.ikpt}")
    seen_bands = set()
    for record in records:
        if record.band in seen_bands:
            raise ValueError(
                f"duplicate band {record.band} for ikpt {record.ikpt} and atom {record.physical_atom}"
            )
        seen_bands.add(record.band)
        if record.k_vector.shape != (3,) or not np.all(np.isfinite(record.k_vector)):
            raise ValueError(f"invalid k-vector metadata for ikpt {record.ikpt}")
        if not np.allclose(
            record.k_vector, reference.k_vector, rtol=0.0, atol=tolerance
        ):
            raise ValueError(f"inconsistent k-vector metadata for ikpt {record.ikpt}")
        if not np.isfinite(record.k_weight) or not np.isclose(
            record.k_weight, reference.k_weight, rtol=0.0, atol=tolerance
        ):
            raise ValueError(f"inconsistent k-weight metadata for ikpt {record.ikpt}")
        if not np.isfinite(record.band_energy) or not np.isfinite(record.occupation):
            raise ValueError(f"nonfinite state metadata for state {record.key}")


def _validated_custom_weights(
    weights: Mapping[StateKey, float] | None,
    records: tuple[StateCharacter, ...],
    mode: str,
) -> dict[StateKey, float] | None:
    if mode != "custom":
        if weights is not None:
            raise ValueError("custom_weights are only valid with custom weighting")
        return None
    if weights is None:
        raise ValueError(
            "custom weighting requires custom_weights keyed by state identity"
        )

    expected = {record.key for record in records}
    supplied: dict[StateKey, float] = {}
    for key, value in weights.items():
        if not isinstance(key, tuple) or len(key) != 3:
            raise ValueError("custom weight keys must be (ikpt,band,physical_atom)")
        if any(
            isinstance(index, bool) or not isinstance(index, (int, np.integer))
            for index in key
        ):
            raise ValueError("custom weight keys must contain integer state indices")
        canonical = (int(key[0]), int(key[1]), int(key[2]))
        if canonical in supplied:
            raise ValueError(f"duplicate custom weight key {canonical}")
        if np.iscomplexobj(value):
            raise ValueError("custom state weights must be real")
        coefficient = float(value)
        if not np.isfinite(coefficient):
            raise ValueError("custom state weights must be finite")
        supplied[canonical] = coefficient
    missing = expected - supplied.keys()
    extra = supplied.keys() - expected
    if missing:
        raise ValueError(f"custom weight is missing state key {min(missing)}")
    if extra:
        raise ValueError(f"custom weight has unselected state key {min(extra)}")
    return supplied


def _k_resolved_weights(
    records: tuple[StateCharacter, ...],
    mode: str,
    custom_weights: Mapping[StateKey, float] | None,
) -> np.ndarray:
    if mode == "unit":
        return np.ones(len(records), dtype=np.float64)
    if mode in ("occupation", "hole"):
        occupations = np.asarray(
            [record.occupation for record in records], dtype=np.float64
        )
        if np.any(occupations < 0.0) or np.any(occupations > SCHEMA_V1_MAX_OCCUPATION):
            raise ValueError(
                "schema-v1 occupation/hole weighting requires occupations in [0,1]"
            )
        if mode == "occupation":
            return occupations
        return SCHEMA_V1_MAX_OCCUPATION - occupations
    if custom_weights is None:
        raise ValueError("internal missing custom weight mapping")
    return np.asarray(
        [custom_weights[record.key] for record in records], dtype=np.float64
    )


def _analyze_k_point(
    ikpt: int,
    k_vector: np.ndarray,
    k_weight: float,
    records: tuple[StateCharacter, ...],
    coefficients: np.ndarray,
    density: np.ndarray,
    min_weight: float | None,
    candidate_channels: tuple[Channel, ...],
    trace_tolerance: float,
) -> KResolvedPoint:
    trace = complex(np.trace(density))
    trace_scale = max(1.0, abs(trace.real))
    if abs(trace.imag) > trace_tolerance * trace_scale:
        raise ValueError(
            f"aggregated density at ikpt {ikpt} has non-negligible imaginary trace {trace.imag}"
        )
    d_weight = float(trace.real)
    rho_hs_norm = float(np.vdot(density, density).real)
    multipoles = density_to_multipoles(density)
    raw = channel_raw_strengths(multipoles)
    hs_weights = channel_hilbert_schmidt_weights(multipoles)
    character = channel_character(multipoles, d_weight, min_weight=min_weight)
    normalized_masked = not np.isfinite(character[VALID_CHANNELS[0]])
    if normalized_masked:
        fractions = _unavailable_channels()
        local_purity = float("nan")
    else:
        fractions = channel_norm_fractions(multipoles, rho_hs_norm)
        local_purity = rho_hs_norm / d_weight**2
    dominant = dominant_channel(character, channels=candidate_channels)
    return KResolvedPoint(
        ikpt=ikpt,
        k_vector=np.asarray(k_vector, dtype=np.float64).copy(),
        k_weight=float(k_weight),
        selected_bands=tuple(record.band for record in records),
        selected_energies=tuple(record.band_energy for record in records),
        selected_occupations=tuple(record.occupation for record in records),
        applied_weights=tuple(float(value) for value in coefficients),
        d_weight=d_weight,
        trace_imaginary=float(trace.imag),
        rho_hs_norm=rho_hs_norm,
        rho_structural=density,
        multipoles=multipoles,
        raw_strengths=raw,
        hs_weights=hs_weights,
        character=character,
        norm_fractions=fractions,
        dominant=dominant,
        normalized_masked=normalized_masked,
        local_purity=local_purity,
    )


def _ordered_channel_subset(
    channels: Iterable[Channel] | None, *, allow_empty: bool
) -> tuple[Channel, ...]:
    if channels is None:
        return () if allow_empty else VALID_CHANNELS
    requested = set()
    for channel in channels:
        if len(channel) != 3:
            raise ValueError(f"invalid l=2 multipole channel {channel}")
        canonical = (int(channel[0]), int(channel[1]), int(channel[2]))
        _validate_channel(*canonical)
        requested.add(canonical)
    result = tuple(channel for channel in VALID_CHANNELS if channel in requested)
    if not result and not allow_empty:
        raise ValueError("channel subset is empty")
    return result


def _validate_root(handle: h5py.File, path: Path) -> None:
    version = _scalar_integer_attribute(handle, "schema_version")
    if version != SCHEMA_VERSION:
        raise ValueError(
            f"unsupported state-character schema version {version} in {path}"
        )
    if _text_attribute(handle, "rho_convention") != RHO_CONVENTION:
        raise ValueError(f"unsupported rho convention in {path}")
    if _text_attribute(handle, "native_basis") != NATIVE_BASIS:
        raise ValueError(f"unsupported native basis in {path}")
    if "states" not in handle or "sites" not in handle:
        raise ValueError(f"missing states/sites group in {path}")


def _read_states(group: h5py.Group, path: Path) -> tuple[StateCharacter, ...]:
    required = ("identity", "k_vector", "scalars", "rho_native_real", "rho_native_imag")
    _require_datasets(group, required, path)
    if _text_attribute(group, "identity_columns") != IDENTITY_COLUMNS:
        raise ValueError(f"unsupported state identity columns in {path}")
    if _text_attribute(group, "scalar_rows") != SCALAR_COLUMNS:
        raise ValueError(f"unsupported state scalar columns in {path}")

    identity = np.asarray(group["identity"])
    k_vector = np.asarray(group["k_vector"], dtype=np.float64)
    scalars = np.asarray(group["scalars"], dtype=np.float64)
    rho_real = np.asarray(group["rho_native_real"], dtype=np.float64)
    rho_imag = np.asarray(group["rho_native_imag"], dtype=np.float64)
    count = identity.shape[0] if identity.ndim == 2 else -1
    expected_shapes = {
        "identity": (count, 6),
        "k_vector": (count, 3),
        "scalars": (count, 7),
        "rho_native_real": (count, 10, 10),
        "rho_native_imag": (count, 10, 10),
    }
    actual = {
        "identity": identity.shape,
        "k_vector": k_vector.shape,
        "scalars": scalars.shape,
        "rho_native_real": rho_real.shape,
        "rho_native_imag": rho_imag.shape,
    }
    for name, expected in expected_shapes.items():
        if actual[name] != expected:
            raise ValueError(
                f"dataset states/{name} has shape {actual[name]}, expected {expected} in {path}"
            )

    output = []
    for index in range(count):
        values = identity[index]
        rho_native = rho_real[index] + 1j * rho_imag[index]
        output.append(
            StateCharacter(
                ikpt=int(values[0]),
                band=int(values[1]),
                physical_atom=int(values[2]),
                atom_type=int(values[3]),
                iatom_l=int(values[4]),
                role_mask=int(values[5]),
                k_vector=k_vector[index].copy(),
                k_weight=float(scalars[index, 0]),
                band_energy=float(scalars[index, 1]),
                occupation=float(scalars[index, 2]),
                rho_native=rho_native,
                source_shard=path,
            )
        )
    return tuple(output)


def _read_sites(group: h5py.Group, path: Path) -> dict[int, SiteCharacter]:
    required = (
        "identity",
        "orbital_global_to_local_real",
        "orbital_global_to_local_imag",
        "spin_native_to_local_real",
        "spin_native_to_local_imag",
    )
    _require_datasets(group, required, path)
    identity = np.asarray(group["identity"])
    if identity.ndim != 2 or identity.shape[1] != 4:
        raise ValueError(f"dataset sites/identity must have shape (nsite,4) in {path}")
    count = identity.shape[0]
    for name, shape in (
        ("orbital_global_to_local_real", (count, 5, 5)),
        ("orbital_global_to_local_imag", (count, 5, 5)),
        ("spin_native_to_local_real", (count, 2, 2)),
        ("spin_native_to_local_imag", (count, 2, 2)),
    ):
        if group[name].shape != shape:
            raise ValueError(
                f"dataset sites/{name} has shape {group[name].shape}, expected {shape} in {path}"
            )

    all_datasets = {name: np.asarray(dataset) for name, dataset in group.items()}
    output = {}
    for index in range(count):
        values = identity[index]
        physical_atom = int(values[0])
        if physical_atom in output:
            raise ValueError(f"duplicate site identity {physical_atom} in {path}")
        datasets = {}
        for name, data in all_datasets.items():
            if data.ndim < 1 or data.shape[0] != count:
                raise ValueError(
                    f"site dataset {name} has inconsistent leading extent in {path}"
                )
            datasets[name] = data[index].copy()
        orbital = (
            datasets["orbital_global_to_local_real"]
            + 1j * datasets["orbital_global_to_local_imag"]
        )
        spin = (
            datasets["spin_native_to_local_real"]
            + 1j * datasets["spin_native_to_local_imag"]
        )
        output[physical_atom] = SiteCharacter(
            physical_atom=physical_atom,
            atom_type=int(values[1]),
            iatom_l=int(values[2]),
            atomic_number=int(values[3]),
            orbital_global_to_local=orbital,
            spin_native_to_local=spin,
            datasets=datasets,
            source_shards=(path,),
        )
    return output


def _merge_site(existing: SiteCharacter, current: SiteCharacter) -> SiteCharacter:
    if (existing.atom_type, existing.iatom_l, existing.atomic_number) != (
        current.atom_type,
        current.iatom_l,
        current.atomic_number,
    ):
        raise ValueError(
            f"inconsistent identity for repeated site {existing.physical_atom}"
        )
    if existing.datasets.keys() != current.datasets.keys():
        raise ValueError(
            f"inconsistent datasets for repeated site {existing.physical_atom}"
        )
    for name in existing.datasets:
        left = existing.datasets[name]
        right = current.datasets[name]
        if np.issubdtype(left.dtype, np.integer):
            consistent = np.array_equal(left, right)
        else:
            consistent = np.allclose(
                left, right, rtol=0.0, atol=1.0e-13, equal_nan=True
            )
        if not consistent:
            raise ValueError(
                f"inconsistent site {existing.physical_atom} dataset {name}"
            )
    return SiteCharacter(
        physical_atom=existing.physical_atom,
        atom_type=existing.atom_type,
        iatom_l=existing.iatom_l,
        atomic_number=existing.atomic_number,
        orbital_global_to_local=existing.orbital_global_to_local,
        spin_native_to_local=existing.spin_native_to_local,
        datasets=existing.datasets,
        source_shards=existing.source_shards + current.source_shards,
    )


def _require_datasets(group: h5py.Group, names: Sequence[str], path: Path) -> None:
    for name in names:
        if name not in group:
            raise ValueError(f"missing dataset {group.name}/{name} in {path}")


def _scalar_integer_attribute(handle: h5py.Group | h5py.File, name: str) -> int:
    if name not in handle.attrs:
        raise ValueError(f"missing attribute {name} on {handle.name}")
    values = np.asarray(handle.attrs[name]).reshape(-1)
    if values.size != 1:
        raise ValueError(f"attribute {name} on {handle.name} is not scalar")
    return int(values[0])


def _text_attribute(handle: h5py.Group | h5py.File, name: str) -> str:
    if name not in handle.attrs:
        raise ValueError(f"missing attribute {name} on {handle.name}")
    values = np.asarray(handle.attrs[name]).reshape(-1)
    if values.size == 1:
        value = values[0]
        return value.decode("utf-8") if isinstance(value, bytes) else str(value)
    if values.dtype.kind in ("S", "U"):
        # FLEUR's HDF helper stores character attributes as rank-one arrays;
        # h5py exposes blank Fortran characters as empty one-byte strings.
        characters = []
        for value in values:
            character = (
                value.decode("utf-8") if isinstance(value, bytes) else str(value)
            )
            characters.append(character if character else " ")
        return "".join(characters).rstrip()
    raise ValueError(f"attribute {name} on {handle.name} is not textual")


def _phase(number: int) -> float:
    return 1.0 if abs(number) % 2 == 0 else -1.0


def _integer_wigner_3j(j1: int, j2: int, j3: int, m1: int, m2: int, m3: int) -> float:
    if m1 + m2 + m3 != 0:
        return 0.0
    if abs(m1) > j1 or abs(m2) > j2 or abs(m3) > j3:
        return 0.0
    if j3 < abs(j1 - j2) or j3 > j1 + j2:
        return 0.0
    delta = sqrt(
        factorial(j1 + j2 - j3)
        * factorial(j1 - j2 + j3)
        * factorial(-j1 + j2 + j3)
        / factorial(j1 + j2 + j3 + 1)
    )
    magnetic = sqrt(
        factorial(j1 + m1)
        * factorial(j1 - m1)
        * factorial(j2 + m2)
        * factorial(j2 - m2)
        * factorial(j3 + m3)
        * factorial(j3 - m3)
    )
    lower = max(0, j2 - j3 - m1, j1 - j3 + m2)
    upper = min(j1 + j2 - j3, j1 - m1, j2 + m2)
    racah_sum = 0.0
    for z in range(lower, upper + 1):
        racah_sum += _phase(z) / (
            factorial(z)
            * factorial(j1 + j2 - j3 - z)
            * factorial(j1 - m1 - z)
            * factorial(j2 + m2 - z)
            * factorial(j3 - j2 + m1 + z)
            * factorial(j3 - j1 - m2 + z)
        )
    return _phase(j1 - j2 - m3) * delta * magnetic * racah_sum


def _orbital_normalization(k: int) -> float:
    return factorial(2 * L) / sqrt(factorial(2 * L - k) * factorial(2 * L + k + 1))


def _double_factorial(number: int) -> int:
    result = 1
    for value in range(number, 0, -2):
        result *= value
    return result


def _generalized_three_bracket(k: int, p: int, r: int) -> complex:
    g = k + p + r
    if g % 2 == 0:
        return complex(_integer_wigner_3j(k, p, r, 0, 0, 0))
    magnitude = sqrt(
        factorial(g - 2 * k)
        * factorial(g - 2 * p)
        * factorial(g - 2 * r)
        / factorial(g + 1)
    )
    magnitude *= _double_factorial(g) / (
        _double_factorial(g - 2 * k)
        * _double_factorial(g - 2 * p)
        * _double_factorial(g - 2 * r)
    )
    return 1j**g * magnitude


def _spherical_pauli_matrices() -> np.ndarray:
    chi = np.zeros((2, 2, 2, 3), dtype=np.complex128)
    chi[0, 0, 0, 1] = 1.0
    chi[1, 1, 0, 1] = 1.0
    chi[0, 0, 1, 1] = 1.0
    chi[1, 1, 1, 1] = -1.0
    chi[0, 1, 1, 2] = -sqrt(2.0)
    chi[1, 0, 1, 0] = sqrt(2.0)
    return chi


__all__ = [
    "SCHEMA_VERSION",
    "SCHEMA_V1_MAX_OCCUPATION",
    "DominantChannel",
    "KResolvedAnalysis",
    "KResolvedPoint",
    "KResolvedSelection",
    "StateCharacter",
    "StateCharacterData",
    "SiteCharacter",
    "VALID_CHANNELS",
    "VALID_COMPONENTS",
    "channel_character",
    "channel_family",
    "channel_hilbert_schmidt_weights",
    "channel_metric_factor",
    "channel_norm_fractions",
    "channel_raw_strengths",
    "channel_subset_sum",
    "compact_multipoles",
    "density_to_multipoles",
    "dominant_channel",
    "flatten_density",
    "multipole_component",
    "read_state_character_shards",
    "analyze_k_resolved",
    "unflatten_density",
    "valid_channels",
    "valid_components",
]
