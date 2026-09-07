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
            raise ValueError(f"state {record.key} has no matching site metadata") from exc
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
        return density_to_multipoles(self.sum_structural_density(records, weights=weights))


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
    return {component: multipole_component(w, *component) for component in VALID_COMPONENTS}


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
    return {channel: channel_metric_factor(*channel) * value for channel, value in raw.items()}


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
    finite = [channel for channel in candidates if np.isfinite(float(character[channel]))]
    if not finite:
        return None
    ranked = sorted(finite, key=lambda channel: (-float(character[channel]), _CHANNEL_ORDER[channel]))
    first = ranked[0]
    if len(ranked) == 1:
        return DominantChannel(first, float(character[first]), None, None, None)
    second = ranked[1]
    first_value = float(character[first])
    second_value = float(character[second])
    return DominantChannel(first, first_value, second, second_value, first_value - second_value)


def channel_family(k: int, p: int, r: int) -> str:
    """Return the broad parity/spin family of a valid l=2 channel."""
    _validate_channel(k, p, r)
    if p == 0:
        return "charge" if k % 2 == 0 else "current/orbital-current"
    return "magnetization" if k % 2 == 0 else "spin-current/spin-orbital"


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


def _validate_root(handle: h5py.File, path: Path) -> None:
    version = _scalar_integer_attribute(handle, "schema_version")
    if version != SCHEMA_VERSION:
        raise ValueError(f"unsupported state-character schema version {version} in {path}")
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
            raise ValueError(f"dataset states/{name} has shape {actual[name]}, expected {expected} in {path}")

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
            raise ValueError(f"dataset sites/{name} has shape {group[name].shape}, expected {shape} in {path}")

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
                raise ValueError(f"site dataset {name} has inconsistent leading extent in {path}")
            datasets[name] = data[index].copy()
        orbital = datasets["orbital_global_to_local_real"] + 1j * datasets[
            "orbital_global_to_local_imag"
        ]
        spin = datasets["spin_native_to_local_real"] + 1j * datasets["spin_native_to_local_imag"]
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
        raise ValueError(f"inconsistent identity for repeated site {existing.physical_atom}")
    if existing.datasets.keys() != current.datasets.keys():
        raise ValueError(f"inconsistent datasets for repeated site {existing.physical_atom}")
    for name in existing.datasets:
        left = existing.datasets[name]
        right = current.datasets[name]
        if np.issubdtype(left.dtype, np.integer):
            consistent = np.array_equal(left, right)
        else:
            consistent = np.allclose(left, right, rtol=0.0, atol=1.0e-13, equal_nan=True)
        if not consistent:
            raise ValueError(f"inconsistent site {existing.physical_atom} dataset {name}")
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
            character = value.decode("utf-8") if isinstance(value, bytes) else str(value)
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
    "DominantChannel",
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
    "compact_multipoles",
    "density_to_multipoles",
    "dominant_channel",
    "flatten_density",
    "multipole_component",
    "read_state_character_shards",
    "unflatten_density",
    "valid_channels",
    "valid_components",
]
