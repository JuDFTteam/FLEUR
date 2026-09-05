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
        self, records: Iterable[StateCharacter] | None = None
    ) -> np.ndarray:
        """Sum explicitly selected structural densities without implicit weights."""
        chosen = self.records if records is None else tuple(records)
        total = np.zeros((10, 10), dtype=np.complex128)
        for record in chosen:
            total += self.structural_density(record)
        return total

    def sum_multipoles(
        self, records: Iterable[StateCharacter] | None = None
    ) -> np.ndarray:
        """Transform the unweighted sum of explicitly selected densities."""
        return density_to_multipoles(self.sum_structural_density(records))


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
    tensor = np.asarray(w)
    if tensor.shape != (5, 2, 6, 11):
        raise ValueError("dense l=2 multipole tensor must have shape (5,2,6,11)")
    if (k, p, r, t) not in _VALID_COMPONENT_SET:
        raise ValueError(f"invalid l=2 multipole component {(k, p, r, t)}")
    return complex(tensor[k, p, r, t + T_OFFSET])


def compact_multipoles(w: np.ndarray) -> dict[tuple[int, int, int, int], complex]:
    """Return the 100 valid components as a canonical-order mapping."""
    return {component: multipole_component(w, *component) for component in VALID_COMPONENTS}


_VALID_COMPONENT_SET = frozenset(VALID_COMPONENTS)


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
    "StateCharacter",
    "StateCharacterData",
    "SiteCharacter",
    "VALID_COMPONENTS",
    "compact_multipoles",
    "density_to_multipoles",
    "flatten_density",
    "multipole_component",
    "read_state_character_shards",
    "unflatten_density",
    "valid_components",
]
