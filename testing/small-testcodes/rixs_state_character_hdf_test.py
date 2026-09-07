#!/usr/bin/env python3
"""Validate the standalone RIXS state-character HDF5 round trip."""

from __future__ import annotations

import argparse
from dataclasses import replace
import hashlib
import importlib.util
from pathlib import Path
import re
import sys
import tempfile

import h5py
import numpy as np


REPOSITORY_ROOT = Path(__file__).resolve().parents[2]
ANALYZER_PATH = REPOSITORY_ROOT / "src" / "tools" / "rixs_state_multipoles.py"
SPEC = importlib.util.spec_from_file_location("rixs_state_multipoles", ANALYZER_PATH)
if SPEC is None or SPEC.loader is None:
    raise RuntimeError(f"cannot import analyzer from {ANALYZER_PATH}")
ANALYZER = importlib.util.module_from_spec(SPEC)
sys.modules[SPEC.name] = ANALYZER
SPEC.loader.exec_module(ANALYZER)


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def golden_density() -> np.ndarray:
    """Reproduce the deterministic physical density used by the Fortran oracle."""
    density = np.zeros((10, 10), dtype=np.complex128)
    for row in range(1, 11):
        density[row - 1, row - 1] = 0.17 * row + 0.011 * row * row
        for column in range(row + 1, 11):
            density[row - 1, column - 1] = complex(
                0.013 * (row + 2 * column) / (1 + column - row),
                0.009 * (2 * row - column) / (1 + column - row),
            )
            density[column - 1, row - 1] = density[row - 1, column - 1].conjugate()
    return density


def golden_multipoles() -> dict[tuple[int, int, int, int], complex]:
    """Read the independently frozen oracle already maintained by the Fortran test."""
    source_path = REPOSITORY_ROOT / "testing" / "small-testcodes" / "spin_orbital_multipoles_test.F90"
    source = source_path.read_text(encoding="utf-8")
    start = source.index("COMPLEX, PARAMETER :: golden_l2(100) = [")
    end = source.index("] ! (4,1,5,5)", start)
    block = source[start : end + len("] ! (4,1,5,5)")]
    pattern = re.compile(
        r"CMPLX\(\s*([-+0-9.eE]+),\s*([-+0-9.eE]+)\s*\)[^!]*!\s*"
        r"\((\d+),(\d+),(\d+),(-?\d+)\)"
    )
    fixture = {
        (int(k), int(p), int(r), int(t)): complex(float(real), float(imaginary))
        for real, imaginary, k, p, r, t in pattern.findall(block)
    }
    assert tuple(fixture) == ANALYZER.VALID_COMPONENTS
    return fixture


def valid_values(tensor: np.ndarray) -> np.ndarray:
    return np.asarray(
        [ANALYZER.multipole_component(tensor, *component) for component in ANALYZER.VALID_COMPONENTS]
    )


def assert_invalid_entries_zero(tensor: np.ndarray) -> None:
    valid = set(ANALYZER.VALID_COMPONENTS)
    for k in range(5):
        for p in range(2):
            for r in range(6):
                for t in range(-5, 6):
                    value = tensor[k, p, r, t + 5]
                    if (k, p, r, t) not in valid:
                        assert value == 0.0j


def rotation_factors() -> tuple[np.ndarray, np.ndarray]:
    """Return deterministic nonidentity l=2 and spin-1/2 frame maps."""
    orbital_angle = 0.37
    orbital = np.diag(np.exp(1j * orbital_angle * np.arange(-2, 3)))
    spin_angle = -0.41
    axis = np.asarray([1.0, 2.0, 3.0])
    axis /= np.linalg.norm(axis)
    sigma_x = np.asarray([[0.0, 1.0], [1.0, 0.0]], dtype=np.complex128)
    sigma_y = np.asarray([[0.0, -1j], [1j, 0.0]], dtype=np.complex128)
    sigma_z = np.asarray([[1.0, 0.0], [0.0, -1.0]], dtype=np.complex128)
    generator = axis[0] * sigma_x + axis[1] * sigma_y + axis[2] * sigma_z
    spin = np.cos(spin_angle / 2.0) * np.eye(2) - 1j * np.sin(spin_angle / 2.0) * generator
    np.testing.assert_allclose(orbital.conj().T @ orbital, np.eye(5), rtol=0.0, atol=1e-14)
    np.testing.assert_allclose(spin.conj().T @ spin, np.eye(2), rtol=0.0, atol=1e-14)
    return orbital, spin


def common_rotation_factors() -> tuple[np.ndarray, np.ndarray]:
    """Return l=2 and spin-1/2 representations of one common SO(3) rotation."""
    axis = np.asarray([1.0, 2.0, -0.7])
    axis /= np.linalg.norm(axis)
    angle = 0.731

    m_values = np.arange(-2, 3, dtype=np.float64)
    orbital_z = np.diag(m_values)
    orbital_plus = np.zeros((5, 5), dtype=np.complex128)
    for column, m_value in enumerate(m_values[:-1]):
        orbital_plus[column + 1, column] = np.sqrt(6.0 - m_value * (m_value + 1.0))
    orbital_x = 0.5 * (orbital_plus + orbital_plus.conj().T)
    orbital_y = (orbital_plus - orbital_plus.conj().T) / (2.0j)

    spin_x = 0.5 * np.asarray([[0.0, 1.0], [1.0, 0.0]], dtype=np.complex128)
    spin_y = 0.5 * np.asarray([[0.0, -1j], [1j, 0.0]], dtype=np.complex128)
    spin_z = 0.5 * np.asarray([[1.0, 0.0], [0.0, -1.0]], dtype=np.complex128)

    def rotate(generator: np.ndarray) -> np.ndarray:
        eigenvalues, eigenvectors = np.linalg.eigh(generator)
        return (eigenvectors * np.exp(-1j * angle * eigenvalues)) @ eigenvectors.conj().T

    orbital = rotate(axis[0] * orbital_x + axis[1] * orbital_y + axis[2] * orbital_z)
    spin = rotate(axis[0] * spin_x + axis[1] * spin_y + axis[2] * spin_z)
    return orbital, spin


def deterministic_psd_density() -> np.ndarray:
    """Return a fixed trace-one positive-semidefinite spin-orbital density."""
    generator = np.random.default_rng(20260907)
    matrix = generator.normal(size=(10, 5)) + 1j * generator.normal(size=(10, 5))
    density = matrix @ matrix.conj().T
    return density / np.trace(density).real


def write_synthetic_shard(
    path: Path,
    identity: tuple[int, int, int, int, int, int],
    rho_native: np.ndarray,
    orbital: np.ndarray,
    spin: np.ndarray,
    *,
    schema_version: int = 1,
) -> None:
    with h5py.File(path, "w") as handle:
        handle.attrs["schema_version"] = schema_version
        handle.attrs["rho_convention"] = ANALYZER.RHO_CONVENTION
        handle.attrs["native_basis"] = ANALYZER.NATIVE_BASIS
        handle.attrs["ligand_atomic_number"] = 8
        states = handle.create_group("states")
        states.attrs["identity_columns"] = ANALYZER.IDENTITY_COLUMNS
        states.attrs["scalar_rows"] = ANALYZER.SCALAR_COLUMNS
        states.create_dataset("identity", data=np.asarray([identity], dtype=np.int32))
        states.create_dataset("k_vector", data=np.asarray([[0.1, 0.2, 0.3]]))
        states.create_dataset("scalars", data=np.asarray([[0.04, 1.25, 0.75, 13.585, 8.0, 5.585, 8.0]]))
        states.create_dataset("rho_native_real", data=np.real(rho_native)[None, :, :])
        states.create_dataset("rho_native_imag", data=np.imag(rho_native)[None, :, :])
        sites = handle.create_group("sites")
        sites.create_dataset("identity", data=np.asarray([[identity[2], identity[3], identity[4], 77]], dtype=np.int32))
        sites.create_dataset("orbital_global_to_local_real", data=np.real(orbital)[None, :, :])
        sites.create_dataset("orbital_global_to_local_imag", data=np.imag(orbital)[None, :, :])
        sites.create_dataset("spin_native_to_local_real", data=np.real(spin)[None, :, :])
        sites.create_dataset("spin_native_to_local_imag", data=np.imag(spin)[None, :, :])


def expect_value_error(action, text: str) -> None:
    try:
        action()
    except ValueError as error:
        assert text in str(error)
    else:
        raise AssertionError(f"expected ValueError containing {text!r}")


def check_channel_table_and_factors() -> tuple[float, tuple[int, int, int]]:
    expected = {
        (0, 0, 0): 1 / 10,
        (0, 1, 1): 1 / 10,
        (1, 0, 1): 1 / 5,
        (1, 1, 0): 1 / 15,
        (1, 1, 1): 9 / 40,
        (1, 1, 2): 2 / 15,
        (2, 0, 2): 1 / 7,
        (2, 1, 1): 2 / 35,
        (2, 1, 2): 25 / 168,
        (2, 1, 3): 3 / 35,
        (3, 0, 3): 1 / 20,
        (3, 1, 2): 3 / 140,
        (3, 1, 3): 49 / 960,
        (3, 1, 4): 1 / 35,
        (4, 0, 4): 1 / 140,
        (4, 1, 3): 1 / 315,
        (4, 1, 4): 81 / 11200,
        (4, 1, 5): 1 / 252,
    }
    assert tuple(expected) == ANALYZER.VALID_CHANNELS
    assert len(ANALYZER.VALID_CHANNELS) == 18
    assert sum(2 * r + 1 for _, _, r in ANALYZER.VALID_CHANNELS) == 100
    errors = {
        channel: abs(ANALYZER.channel_metric_factor(*channel) - value)
        for channel, value in expected.items()
    }
    location = max(errors, key=errors.get)
    maximum = errors[location]
    assert maximum < 2.0e-15
    expect_value_error(lambda: ANALYZER.channel_metric_factor(0, 0, 1), "invalid l=2")
    assert ANALYZER.channel_family(0, 0, 0) == "charge"
    assert ANALYZER.channel_family(0, 1, 1) == "magnetization"
    assert ANALYZER.channel_family(1, 0, 1) == "current/orbital-current"
    assert ANALYZER.channel_family(1, 1, 1) == "spin-current/spin-orbital"
    return maximum, location


def check_channel_strength_identities() -> dict[str, float]:
    identity = np.eye(10, dtype=np.complex128) / 10.0
    densities = (identity, golden_density(), deterministic_psd_density())
    parseval_error = 0.0
    baseline_error = 0.0
    fraction_error = 0.0
    purity_error = 0.0
    for density in densities:
        multipoles = ANALYZER.density_to_multipoles(density)
        raw = ANALYZER.channel_raw_strengths(multipoles)
        weights = ANALYZER.channel_hilbert_schmidt_weights(multipoles)
        hs_norm = float(np.vdot(density, density).real)
        d_weight = float(np.trace(density).real)
        character = ANALYZER.channel_character(multipoles, d_weight)
        fractions = ANALYZER.channel_norm_fractions(multipoles, hs_norm)
        assert tuple(raw) == ANALYZER.VALID_CHANNELS
        assert tuple(weights) == ANALYZER.VALID_CHANNELS
        parseval_error = max(parseval_error, abs(sum(weights.values()) - hs_norm))
        baseline_error = max(baseline_error, abs(character[(0, 0, 0)] - 0.1))
        fraction_error = max(fraction_error, abs(sum(fractions.values()) - 1.0))
        purity_error = max(
            purity_error,
            abs(sum(character.values()) - hs_norm / d_weight**2),
        )
    assert parseval_error < 2.0e-12
    assert baseline_error < 2.0e-14
    assert fraction_error < 2.0e-14
    assert purity_error < 2.0e-14

    scalar_w = ANALYZER.density_to_multipoles(identity)
    scalar_raw = ANALYZER.channel_raw_strengths(scalar_w)
    assert abs(scalar_raw[(0, 0, 0)] - 1.0) < 2.0e-14
    assert max(value for channel, value in scalar_raw.items() if channel != (0, 0, 0)) < 2.0e-28
    invalid_changed = scalar_w.copy()
    invalid_changed[0, 0, 1, ANALYZER.T_OFFSET] = 1.0e20
    assert ANALYZER.channel_raw_strengths(invalid_changed) == scalar_raw
    return {
        "parseval": parseval_error,
        "baseline": baseline_error,
        "fraction": fraction_error,
        "purity": purity_error,
    }


def check_channel_rotation_and_scale() -> dict[str, float]:
    density = deterministic_psd_density()
    orbital, spin = common_rotation_factors()
    transform = np.kron(orbital, spin)
    rotated = transform @ density @ transform.conj().T
    original_w = ANALYZER.density_to_multipoles(density)
    rotated_w = ANALYZER.density_to_multipoles(rotated)
    original_raw = ANALYZER.channel_raw_strengths(original_w)
    rotated_raw = ANALYZER.channel_raw_strengths(rotated_w)
    original_q = ANALYZER.channel_hilbert_schmidt_weights(original_w)
    rotated_q = ANALYZER.channel_hilbert_schmidt_weights(rotated_w)

    raw_absolute = max(abs(original_raw[c] - rotated_raw[c]) for c in ANALYZER.VALID_CHANNELS)
    q_absolute = max(abs(original_q[c] - rotated_q[c]) for c in ANALYZER.VALID_CHANNELS)
    raw_relative = max(
        abs(original_raw[c] - rotated_raw[c]) / max(1.0, original_raw[c], rotated_raw[c])
        for c in ANALYZER.VALID_CHANNELS
    )
    q_relative = max(
        abs(original_q[c] - rotated_q[c]) / max(1.0, original_q[c], rotated_q[c])
        for c in ANALYZER.VALID_CHANNELS
    )
    assert raw_relative < 2.0e-13
    assert q_relative < 2.0e-13

    scale = 0.37
    scaled_density = scale * density
    scaled_w = ANALYZER.density_to_multipoles(scaled_density)
    scaled_q = ANALYZER.channel_hilbert_schmidt_weights(scaled_w)
    q_scale_error = max(
        abs(scaled_q[channel] - scale**2 * original_q[channel])
        for channel in ANALYZER.VALID_CHANNELS
    )
    d_weight = float(np.trace(density).real)
    scaled_d_weight = float(np.trace(scaled_density).real)
    d_scale_error = abs(scaled_d_weight - scale * d_weight)
    original_character = ANALYZER.channel_character(original_w, d_weight)
    scaled_character = ANALYZER.channel_character(scaled_w, scaled_d_weight)
    character_scale_error = max(
        abs(original_character[channel] - scaled_character[channel])
        for channel in ANALYZER.VALID_CHANNELS
    )
    original_norm = float(np.vdot(density, density).real)
    scaled_norm = float(np.vdot(scaled_density, scaled_density).real)
    original_fractions = ANALYZER.channel_norm_fractions(original_w, original_norm)
    scaled_fractions = ANALYZER.channel_norm_fractions(scaled_w, scaled_norm)
    fraction_scale_error = max(
        abs(original_fractions[channel] - scaled_fractions[channel])
        for channel in ANALYZER.VALID_CHANNELS
    )
    assert max(q_scale_error, d_scale_error, character_scale_error, fraction_scale_error) < 2.0e-13
    return {
        "raw_rotation_absolute": raw_absolute,
        "raw_rotation_relative": raw_relative,
        "q_rotation_absolute": q_absolute,
        "q_rotation_relative": q_relative,
        "q_scale": q_scale_error,
        "d_scale": d_scale_error,
        "character_scale": character_scale_error,
        "fraction_scale": fraction_scale_error,
    }


def check_masking_and_dominance() -> None:
    density = np.eye(10, dtype=np.complex128) / 10.0
    multipoles = ANALYZER.density_to_multipoles(density)
    assert all(np.isfinite(value) for value in ANALYZER.channel_character(multipoles, 1.0).values())
    just_above = ANALYZER.channel_character(multipoles, 0.5000000000001, min_weight=0.5)
    assert all(np.isfinite(value) for value in just_above.values())
    below = ANALYZER.channel_character(multipoles, 0.49, min_weight=0.5)
    zero = ANALYZER.channel_character(multipoles, 0.0)
    zero_norm = ANALYZER.channel_norm_fractions(multipoles, 0.0)
    assert all(np.isnan(value) for value in below.values())
    assert all(np.isnan(value) for value in zero.values())
    assert all(np.isnan(value) for value in zero_norm.values())
    assert ANALYZER.dominant_channel(below) is None
    expect_value_error(
        lambda: ANALYZER.channel_character(multipoles, 1.0, min_weight=-1.0),
        "nonnegative",
    )

    values = {channel: 0.0 for channel in ANALYZER.VALID_CHANNELS}
    values[(0, 0, 0)] = 5.0
    values[(0, 1, 1)] = 1.0
    values[(1, 0, 1)] = 2.0
    values[(1, 1, 0)] = 3.0
    values[(4, 1, 5)] = 4.0
    overall = ANALYZER.dominant_channel(values)
    anisotropic = ANALYZER.dominant_channel(values, exclude_scalar=True)
    restricted = ANALYZER.dominant_channel(
        values,
        channels=((0, 1, 1), (1, 0, 1), (1, 1, 0), (4, 1, 5)),
    )
    assert overall is not None and overall.channel == (0, 0, 0) and overall.second_channel == (4, 1, 5)
    assert anisotropic is not None and anisotropic.channel == (4, 1, 5)
    assert restricted is not None and restricted.channel == (4, 1, 5)

    values[(0, 1, 1)] = 4.0
    tied = ANALYZER.dominant_channel(
        values,
        channels=((4, 1, 5), (0, 1, 1)),
    )
    assert tied is not None and tied.channel == (0, 1, 1)
    assert tied.second_channel == (4, 1, 5) and tied.gap == 0.0


def check_complete_golden() -> tuple[float, tuple[int, int, int, int], float, float]:
    density = golden_density()
    expected = golden_multipoles()
    tensor = ANALYZER.density_to_multipoles(density)
    errors = {component: abs(ANALYZER.multipole_component(tensor, *component) - value) for component, value in expected.items()}
    location = max(errors, key=errors.get)
    maximum = errors[location]
    assert maximum < 2.0e-12
    assert_invalid_entries_zero(tensor)
    trace_error = abs(ANALYZER.multipole_component(tensor, 0, 0, 0, 0) - np.trace(density))
    assert trace_error < 2.0e-12
    hermiticity = max(
        abs(
            ANALYZER.multipole_component(tensor, k, p, r, t).conjugate()
            - (-1) ** t * ANALYZER.multipole_component(tensor, k, p, r, -t)
        )
        for k, p, r, t in ANALYZER.VALID_COMPONENTS
    )
    assert hermiticity < 2.0e-12
    return maximum, location, trace_error, hermiticity


def check_negative_conventions() -> dict[str, float]:
    correct_density = golden_density()
    correct = valid_values(ANALYZER.density_to_multipoles(correct_density))
    rho4 = ANALYZER.unflatten_density(correct_density)
    wrong = {
        "orbital_partial_transpose": ANALYZER.flatten_density(rho4.transpose(2, 1, 0, 3)),
        "swapped_spin_indices": ANALYZER.flatten_density(rho4[:, ::-1, :, ::-1]),
        "complex_conjugation": correct_density.conjugate(),
        "wrong_flattening": rho4.transpose(1, 0, 3, 2).reshape(10, 10),
    }
    separations = {
        name: float(np.max(np.abs(valid_values(ANALYZER.density_to_multipoles(value)) - correct)))
        for name, value in wrong.items()
    }
    assert all(value > 1.0e-8 for value in separations.values())
    return separations


def check_frame_and_shards() -> dict[str, float]:
    rho_structural = golden_density()
    orbital, spin = rotation_factors()
    transform = np.kron(orbital, spin)
    rho_native = transform.conj().T @ rho_structural @ transform
    second_structural = 0.37 * rho_structural
    second_native = transform.conj().T @ second_structural @ transform

    with tempfile.TemporaryDirectory(prefix="rixs_state_multipoles_") as directory:
        root = Path(directory)
        shard_one = root / "test_state_character_rank0000.hdf"
        shard_two = root / "test_state_character_rank0001.hdf"
        write_synthetic_shard(shard_one, (3, 5, 4, 2, 2, 3), rho_native, orbital, spin)
        write_synthetic_shard(shard_two, (4, 6, 4, 2, 2, 2), second_native, orbital, spin)
        data = ANALYZER.read_state_character_shards((shard_two, shard_one))
        assert [record.key for record in data.records] == [(3, 5, 4), (4, 6, 4)]
        assert data.select(physical_atom=4, ikpt=3, band=5, role_mask=1) == (data.records[0],)
        assert data.sites[4].source_shards == (shard_two, shard_one)

        reconstructed = data.structural_density(data.records[0])
        frame_error = float(np.max(np.abs(reconstructed - rho_structural)))
        assert frame_error < 2.0e-12
        expected = valid_values(ANALYZER.density_to_multipoles(rho_structural))
        reconstructed_w = valid_values(data.multipoles(data.records[0]))
        frame_tensor_error = float(np.max(np.abs(reconstructed_w - expected)))
        assert frame_tensor_error < 2.0e-12
        omitted_frame_separation = float(
            np.max(np.abs(valid_values(ANALYZER.density_to_multipoles(rho_native)) - expected))
        )
        assert omitted_frame_separation > 1.0e-8

        summed_once = data.sum_multipoles()
        summed_separately = data.multipoles(data.records[0]) + data.multipoles(data.records[1])
        linearity_error = float(np.max(np.abs(summed_once - summed_separately)))
        assert linearity_error < 2.0e-12
        total_q = ANALYZER.channel_hilbert_schmidt_weights(summed_once)
        first_q = ANALYZER.channel_hilbert_schmidt_weights(data.multipoles(data.records[0]))
        second_q = ANALYZER.channel_hilbert_schmidt_weights(data.multipoles(data.records[1]))
        strength_nonadditivity = max(
            abs(total_q[channel] - first_q[channel] - second_q[channel])
            for channel in ANALYZER.VALID_CHANNELS
        )
        assert strength_nonadditivity > 1.0e-8

        coefficients = (0.25, -0.4)
        weighted_density = data.sum_structural_density(weights=coefficients)
        expected_weighted_density = coefficients[0] * rho_structural + coefficients[1] * second_structural
        weighted_density_error = float(np.max(np.abs(weighted_density - expected_weighted_density)))
        assert weighted_density_error < 2.0e-12
        weighted_multipoles = data.sum_multipoles(weights=coefficients)
        expected_weighted_multipoles = (
            coefficients[0] * data.multipoles(data.records[0])
            + coefficients[1] * data.multipoles(data.records[1])
        )
        weighted_linearity_error = float(
            np.max(np.abs(weighted_multipoles - expected_weighted_multipoles))
        )
        assert weighted_linearity_error < 2.0e-12
        expect_value_error(lambda: data.sum_structural_density(weights=(1.0,)), "exactly 2")
        expect_value_error(lambda: data.sum_structural_density(weights=(1.0, 1.0j)), "must be real")

        changed_records = tuple(
            replace(record, k_weight=17.0 + index, occupation=-4.0 - index)
            for index, record in enumerate(data.records)
        )
        changed_metadata = ANALYZER.StateCharacterData(changed_records, data.sites)
        metadata_weighting_error = float(
            np.max(np.abs(changed_metadata.sum_structural_density() - data.sum_structural_density()))
        )
        assert metadata_weighting_error == 0.0

        duplicate = root / "duplicate_state_character_rank0002.hdf"
        write_synthetic_shard(duplicate, (3, 5, 4, 2, 2, 3), rho_native, orbital, spin)
        expect_value_error(
            lambda: ANALYZER.read_state_character_shards((shard_one, duplicate)), "duplicate state key"
        )

        inconsistent = root / "inconsistent_state_character_rank0003.hdf"
        changed_orbital = orbital.copy()
        changed_orbital[0, 0] += 1.0e-4
        write_synthetic_shard(inconsistent, (5, 7, 4, 2, 2, 1), rho_native, changed_orbital, spin)
        expect_value_error(
            lambda: ANALYZER.read_state_character_shards((shard_one, inconsistent)), "inconsistent site"
        )

        unknown = root / "unknown_state_character_rank0004.hdf"
        write_synthetic_shard(unknown, (5, 7, 4, 2, 2, 1), rho_native, orbital, spin, schema_version=2)
        expect_value_error(lambda: ANALYZER.read_state_character_shards(unknown), "unsupported state-character schema")

    return {
        "frame_density": frame_error,
        "frame_tensor": frame_tensor_error,
        "omitted_frame": omitted_frame_separation,
        "linearity": linearity_error,
        "strength_nonadditivity": strength_nonadditivity,
        "weighted_density": weighted_density_error,
        "weighted_linearity": weighted_linearity_error,
        "metadata_weighting": metadata_weighting_error,
    }


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("file", type=Path)
    args = parser.parse_args()

    before_hash = sha256(args.file)
    with h5py.File(args.file, "r") as handle:
        def scalar_attr(name: str) -> int:
            return int(np.asarray(handle.attrs[name]).reshape(-1)[0])

        assert scalar_attr("schema_version") == 1
        assert scalar_attr("global_rank") == 7
        assert scalar_attr("ligand_atomic_number") == 8
        assert scalar_attr("number_of_states") == 1
        assert scalar_attr("number_of_sites") == 1

        states = handle["states"]
        np.testing.assert_array_equal(states["identity"][0], [3, 5, 4, 2, 2, 3])
        np.testing.assert_allclose(states["k_vector"][0], [0.1, 0.2, 0.3], rtol=0.0, atol=1e-14)
        np.testing.assert_allclose(
            states["scalars"][0], [0.04, 1.25, 0.75, 0.9, 0.45, 0.45, 0.45], rtol=0.0, atol=1e-14
        )
        np.testing.assert_allclose(
            states["real_orbital_weights"][0], [0.1, 0.2, 0.15, 0.25, 0.2], rtol=0.0, atol=1e-14
        )
        np.testing.assert_allclose(states["jeff_weights"][0], [0.3, 0.15], rtol=0.0, atol=1e-14)
        np.testing.assert_allclose(
            states["jeff_mj_weights"][0], [0.1, 0.2, 0.01, 0.02, 0.03, 0.09], rtol=0.0, atol=1e-14
        )
        expected_rho = np.empty((10, 10), dtype=np.complex128)
        for column in range(1, 11):
            for row in range(1, 11):
                expected_rho[row - 1, column - 1] = complex(0.01 * (row + 10 * column), 0.001 * (row - column))
        stored_rho = states["rho_native_real"][0] + 1j * states["rho_native_imag"][0]
        np.testing.assert_allclose(stored_rho, expected_rho, rtol=0.0, atol=1e-14)

        sites = handle["sites"]
        np.testing.assert_array_equal(sites["identity"][0], [4, 2, 2, 77])
        np.testing.assert_array_equal(sites["ligand_atom_ids"][0], [11, 12, 13, 14, 15, 16])
        np.testing.assert_allclose(sites["local_to_global"][0], np.eye(3), rtol=0.0, atol=1e-14)
        np.testing.assert_allclose(sites["reference_frame"][0], np.eye(3), rtol=0.0, atol=1e-14)
        np.testing.assert_allclose(sites["bond_distances"][0], [1.0, 1.1, 1.2, 1.3, 1.4, 1.5])
        np.testing.assert_allclose(sites["shell_diagnostics"][0], [1.8, 0.3, 1.2], rtol=0.0, atol=1e-14)
        np.testing.assert_array_equal(sites["opposite_pairs"][0], [[1, 2], [3, 4], [5, 6]])

    loaded = ANALYZER.read_state_character_shards(args.file)
    assert len(loaded.records) == 1
    assert loaded.records[0].key == (3, 5, 4)
    np.testing.assert_allclose(loaded.records[0].rho_native, expected_rho, rtol=0.0, atol=1e-14)
    expected_orbital, expected_spin = rotation_factors()
    np.testing.assert_allclose(loaded.sites[4].orbital_global_to_local, expected_orbital, rtol=0.0, atol=1e-14)
    np.testing.assert_allclose(loaded.sites[4].spin_native_to_local, expected_spin, rtol=0.0, atol=1e-14)
    expected_transform = np.kron(expected_orbital, expected_spin)
    expected_structural = expected_transform @ expected_rho @ expected_transform.conj().T
    np.testing.assert_allclose(
        loaded.structural_density(loaded.records[0]), expected_structural, rtol=0.0, atol=1e-13
    )
    rho4 = ANALYZER.unflatten_density(expected_rho)
    for m in range(-2, 3):
        for spin in range(1, 3):
            for mp in range(-2, 3):
                for spinp in range(1, 3):
                    row = 2 * (m + 2) + spin - 1
                    column = 2 * (mp + 2) + spinp - 1
                    assert rho4[m + 2, spin - 1, mp + 2, spinp - 1] == expected_rho[row, column]
    assert sha256(args.file) == before_hash

    factor_error, factor_location = check_channel_table_and_factors()
    strength_checks = check_channel_strength_identities()
    rotation_scale_checks = check_channel_rotation_and_scale()
    check_masking_and_dominance()
    golden_error, golden_location, trace_error, hermiticity_error = check_complete_golden()
    shard_checks = check_frame_and_shards()
    negative = check_negative_conventions()
    print(f"Channel metric-factor error: {factor_error:.6e} at {factor_location}")
    print(f"Channel Parseval error: {strength_checks['parseval']:.6e}")
    print(f"Universal C_orth(000) error: {strength_checks['baseline']:.6e}")
    print(f"Channel F normalization error: {strength_checks['fraction']:.6e}")
    print(f"Channel C_orth purity error: {strength_checks['purity']:.6e}")
    print(
        "Raw-strength rotation error: "
        f"abs={rotation_scale_checks['raw_rotation_absolute']:.6e} "
        f"rel={rotation_scale_checks['raw_rotation_relative']:.6e}"
    )
    print(
        "Hilbert-Schmidt rotation error: "
        f"abs={rotation_scale_checks['q_rotation_absolute']:.6e} "
        f"rel={rotation_scale_checks['q_rotation_relative']:.6e}"
    )
    print(
        "Scale-invariance errors: "
        f"D={rotation_scale_checks['d_scale']:.6e} "
        f"Q={rotation_scale_checks['q_scale']:.6e} "
        f"C_orth={rotation_scale_checks['character_scale']:.6e} "
        f"F={rotation_scale_checks['fraction_scale']:.6e}"
    )
    print(f"Complete 100-component golden error: {golden_error:.6e} at {golden_location}")
    print(f"Trace identity error: {trace_error:.6e}")
    print(f"Hermiticity relation error: {hermiticity_error:.6e}")
    print(f"Structural frame-chain density error: {shard_checks['frame_density']:.6e}")
    print(f"Structural frame-chain tensor error: {shard_checks['frame_tensor']:.6e}")
    print(f"Omitted-frame separation: {shard_checks['omitted_frame']:.6e}")
    print(f"Density-sum versus multipole-sum error: {shard_checks['linearity']:.6e}")
    print(f"Density-first strength nonadditivity: {shard_checks['strength_nonadditivity']:.6e}")
    print(f"Explicitly weighted density error: {shard_checks['weighted_density']:.6e}")
    print(f"Explicitly weighted multipole error: {shard_checks['weighted_linearity']:.6e}")
    print(f"Implicit metadata-weighting error: {shard_checks['metadata_weighting']:.6e}")
    for name, value in negative.items():
        print(f"Negative-control separation {name}: {value:.6e}")
    print("RIXS STATE CHARACTER HDF TEST: PASS")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
