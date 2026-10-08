
import math
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

BASE = "./inputfiles/noco/FFNConstraint"


def _moments(outxml, kind="localMagMoment"):
    root = ElementTree.parse(outxml).getroot()
    return [tuple(float(v) for v in e.attrib["vec"].split())
            for e in root.findall(f".//magneticMomentsInMTSpheres/{kind}")]


def _misalignments(outxml):
    """Signed |m_perp|/|m| in the local frame, one entry per SCF iteration."""
    out = []
    for mx, my, mz in _moments(outxml):
        norm = math.sqrt(mx * mx + my * my + mz * mz)
        assert norm > 0.0, "zero-magnitude local moment"
        out.append(math.copysign(math.hypot(mx, my), mx + my) / norm)
    return out


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


@pytest.mark.fleur
@pytest.mark.noco
@pytest.mark.hdf
@pytest.mark.slow
def test_ffn_transverse_constraint_reaches_hamiltonian(execute_fleur):
    """The transverse constraining field must survive vgen_finalize in FFN (#800).

    For l_mtNocoPot=T, vgen_finalize calls rotate_mt_den_from_local, which zeroes
    vTot%mt and rebuilds all four spin components from the local-frame diagonal
    pair plus theta_mt/phi_mt. vgen_constraint therefore has to run *after* that
    step. When it ran before it, l_constrained=T produced output bit-identical to
    l_constrained=F -- which is what the second assertion below catches.
    """
    def seeded(seed_dir, seed_file, run_dir, run_file):
        """Flip the moment to build a tilted cdn.hdf, then run the SCF from it."""
        execute_fleur(seed_dir, rm_files=['.'],
                      only_copy=[[seed_file, "inp.xml"]], mpi_procs=1)
        # no rm_files here on purpose: this run inherits the seed's cdn.hdf
        res = execute_fleur(run_dir, only_copy=[[run_file, "inp.xml"]], mpi_procs=1)
        return _misalignments(res["out.xml"])

    mis_on = seeded(f"{BASE}/stage1", "inp.xml", f"{BASE}/stage2", "inp.xml")
    mis_off = seeded(f"{BASE}/stage1", "inp.xml", BASE, "inp_ffn_off.xml")
    mis_on_m = seeded(BASE, "inp_seed_m.xml", f"{BASE}/stage2", "inp.xml")

    assert len(mis_on) >= 3 and len(mis_on) == len(mis_off) == len(mis_on_m)

    # both runs start from the same density, so iteration 1 must agree
    assert abs(mis_on[0] - mis_off[0]) < 1.0e-9

    # the constraint must actually act -- this is what was false before the fix
    assert abs(mis_on[-1] - mis_off[-1]) > 1.0e-4, \
        "l_constrained=T is indistinguishable from l_constrained=F in FFN mode"

    # ... and it must pull the moment back towards the local axis
    assert abs(mis_on[-1]) < abs(mis_on[0]), f"constraint is not restoring: {mis_on}"
    assert abs(mis_on[-1]) < 0.75 * abs(mis_off[-1])

    # opposite initial deviation -> mirrored response of the same magnitude
    assert mis_on[0] * mis_on_m[0] < 0.0
    assert max(abs(abs(a) - abs(b)) for a, b in zip(mis_on, mis_on_m)) < 1.0e-6


@pytest.mark.fleur
@pytest.mark.noco
@pytest.mark.hdf
def test_mt_collinear_constraint_path_unchanged(execute_fleur):
    """The non-FFN (l_mtNocoPot=F) constrained path must be unaffected (#800).

    vgen_constraint moved from before to after vgen_finalize and became additive.
    For l_mtNocoPot=F components 3/4 are still zero at the new insertion point, so
    nothing may change. Here the local moment is aligned with the local axis, so
    the constraining field is exactly zero and switching l_constrained on must
    change nothing at all -- this catches a relocation that injects a spurious
    transverse field into the collinear-MT path.
    """
    on = execute_fleur(BASE, rm_files=['.'],
                       only_copy=[["inp_mtc_on.xml", "inp.xml"]], mpi_procs=1)

    moments = _moments(on["out.xml"])
    assert moments, "no moments written by the non-FFN constrained run"
    assert all(v == v for m in moments for v in m), "NaN in the non-FFN constrained moments"

    root = ElementTree.parse(on["out.xml"]).getroot()
    assert root.findall(".//constrainingField"), \
        "constrainingField block missing from the non-FFN constrained run"
    # prefix match on purpose: the emitted tag is truncated to 15 characters by a
    # fixed-length label in update_b_cons, so it reads DeltaBConstrain.
    deltas = [float(e.attrib["DBConX"]) for e in root.iter()
              if e.tag.startswith("DeltaBConstrain")]
    assert deltas, "no DeltaBConstrain entries recorded"
    assert max(abs(d) for d in deltas) < 1.0e-8, \
        f"spurious transverse constraining field in the collinear-MT path: {deltas}"
