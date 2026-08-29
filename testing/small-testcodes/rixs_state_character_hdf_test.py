#!/usr/bin/env python3
"""Validate the standalone RIXS state-character HDF5 round trip."""

from __future__ import annotations

import argparse
from pathlib import Path

import h5py
import numpy as np


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("file", type=Path)
    args = parser.parse_args()

    with h5py.File(args.file, "r") as handle:
        scalar_attr = lambda name: int(np.asarray(handle.attrs[name]).reshape(-1)[0])
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

    print("RIXS STATE CHARACTER HDF TEST: PASS")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
