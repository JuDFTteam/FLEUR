
from pathlib import Path
from xml.etree import ElementTree

import pytest

"""
Differential regression tests for the noncollinear correctness fixes.

Unlike the reference-out.xml rows in tests.md these assert a physical invariant
rather than stored numbers, so they fail against the pre-fix code.
See issues #797, #798 and #800.

All runs share the session work dir and are chained deliberately: a run that
follows another picks up its cdn.hdf, which is how the seeded constraint cases
below get their tilted start density. Pass rm_files=['.'] to start clean.
"""

@pytest.mark.fleur
@pytest.mark.noco
@pytest.mark.soc
@pytest.mark.ldau
@pytest.mark.hdf
def test_full_mt_noco_zero_u_invariant(execute_fleur):
    """A skipped zero-U term must not change the first eigenproblem (#798).

    The LDA+U branch is gated on the *count* n_u+n_hia+n_opc, not on the value of
    U, so U=0.0 still enters the spin off-diagonal path and still calls rad_ovlp.
    Before the fix rad_ovlp handed the caller's shared t_usdus to radfun, which
    overwrote the radial boundary data the Hamiltonian is built from -- so adding
    a numerically inert U=0.0 term shifted the eigenvalues.
    """
    input_dir = "./inputfiles/noco/FFNZeroUInvariant"
    common = ["kpts.xml", "sym.xml"]

    no_u = execute_fleur(input_dir, rm_files=['.'],
                         only_copy=[["inp_no_u.xml", "inp.xml"], *common], mpi_procs=1)
    no_u_eig, no_u_mom = _eigen_and_moments(no_u["out.xml"])

    zero_u = execute_fleur(input_dir, rm_files=['.'],
                           only_copy=[["inp_zero_u.xml", "inp.xml"], *common], mpi_procs=1)
    zero_u_eig, zero_u_mom = _eigen_and_moments(zero_u["out.xml"])
    zero_u_stdout = Path(zero_u["stdout"]).read_text()

    assert no_u_eig
    assert len(zero_u_eig) == len(no_u_eig)
    assert len(zero_u_mom) == len(no_u_mom)
    assert max(abs(a - b) for a, b in zip(zero_u_eig, no_u_eig)) < 1.0e-9
    assert max(abs(a - b) for a, b in zip(zero_u_mom, no_u_mom)) < 1.0e-9
    assert "no density matrix found ... skipping LDA+U" in zero_u_stdout


def _eigen_and_moments(outxml):
    root = ElementTree.parse(outxml).getroot()
    eigenvalues = [float(v) for e in root.findall(".//eigenvaluesAt") for v in e.text.split()]
    moments = [float(v)
               for e in root.findall(".//magneticMomentsInMTSpheres/globalMagMoment")
               for v in e.attrib["vec"].split()]
    return eigenvalues, moments
