import re
import pytest


@pytest.mark.fleur
@pytest.mark.bulk
@pytest.mark.lo
@pytest.mark.hdf
@pytest.mark.fast
def test_lapwfn_DiamondLO(default_fleur_test):
    """
    CASINO export (output/lapwfn, src/contributed/CASINO): a non-relativistic SCF of diamond with
    a 3p local orbital (warp_factor), then the writer. Checks the header numbers of lapwfn.data,
    that every band at the k = G/2 points is written real, and the inversion relation between the
    two atoms on the LAPW and on the local-orbital coefficients.
    """
    res = default_fleur_test("lapwfn/DiamondLO", mpi_procs=1)
    for name in ["lapwfn.data", "pwfn.data", "lapwfn_measurements.out", "lapwfn_ylm.out"]:
        assert name in res

    with open(res["lapwfn.data"]) as f:
        lines = f.read().splitlines()

    def after(label):
        return lines[lines.index(label) + 1].strip()

    assert lines[0] == "LAPWFN.DATA: LAPW orbitals, format version 1"
    assert int(after("Number of electrons per primitive cell (core included)")) == 12
    assert int(after("Number of valence electrons per primitive cell (in the WAVE FUNCTION section)")) == 8
    assert float(after("Total energy (au per primitive cell)")) == pytest.approx(-75.5613677331, abs=1e-8)
    assert float(after("Fermi energy (au)")) == pytest.approx(0.44050331, abs=2e-5)
    assert int(after("Number of k-points")) == 4
    assert after("Number of radial functions for l = 0 to L").split() == ["2", "3", "2", "2", "2", "2", "2"]

    flags = [cur.split()[-1] for prev, cur in zip(lines, lines[1:]) if prev.startswith("Band, spin, eigenvalue")]
    assert len(flags) >= 4 * 4
    assert all(flag == "T" for flag in flags)

    with open(res["lapwfn_measurements.out"]) as f:
        inversion = [line for line in f if line.startswith("inversion")]
    assert len(inversion) == 4
    for line in inversion:
        res_lapw, res_lo = [float(x) for x in re.search(r"\(-1\)\^m\s+(\S+)\s+(\S+)", line).groups()]
        assert res_lapw < 1e-10 and res_lo < 1e-10
