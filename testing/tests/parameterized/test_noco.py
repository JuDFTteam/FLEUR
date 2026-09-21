
import pytest
from pathlib import Path
from xml.etree import ElementTree
"""
Boiler-plate code for testset
"""
from read_tests import read_tests
all_tests = read_tests("noco")

@pytest.mark.fleur
@pytest.mark.noco
@pytest.mark.parametrize(("dir","desc","cmdline","mpi_procs"), all_tests)
def test_noco(dir,desc,cmdline,mpi_procs,default_fleur_test):
    """
    """

    assert default_fleur_test(dir, cmdline_args=cmdline, mpi_procs=mpi_procs)


@pytest.mark.fleur
@pytest.mark.noco
@pytest.mark.soc
@pytest.mark.ldau
@pytest.mark.hdf
def test_full_mt_noco_zero_u_invariant(execute_fleur, work_dir):
    """A skipped zero-U term must not change the first eigenproblem."""

    input_dir = "./inputfiles/noco/FFNZeroUInvariant"
    common_files = ["kpts.xml", "sym.xml"]
    Path(work_dir, "no_u").mkdir()
    Path(work_dir, "zero_u").mkdir()

    no_u = execute_fleur(
        input_dir,
        only_copy=[["inp_no_u.xml", "inp.xml"], *common_files],
        sub_dir="no_u",
        mpi_procs=1,
    )
    zero_u = execute_fleur(
        input_dir,
        only_copy=[["inp_zero_u.xml", "inp.xml"], *common_files],
        sub_dir="zero_u",
        mpi_procs=1,
    )

    def first_iteration_values(outxml):
        root = ElementTree.parse(outxml).getroot()
        eigenvalues = [
            float(value)
            for element in root.findall(".//eigenvaluesAt")
            for value in element.text.split()
        ]
        moments = [
            float(value)
            for element in root.findall(".//magneticMomentsInMTSpheres/globalMagMoment")
            for value in element.attrib["vec"].split()
        ]
        return eigenvalues, moments

    no_u_eigenvalues, no_u_moments = first_iteration_values(no_u["out.xml"])
    zero_u_eigenvalues, zero_u_moments = first_iteration_values(zero_u["out.xml"])
    assert no_u_eigenvalues
    assert len(zero_u_eigenvalues) == len(no_u_eigenvalues)
    assert len(zero_u_moments) == len(no_u_moments)
    assert max(abs(a - b) for a, b in zip(zero_u_eigenvalues, no_u_eigenvalues)) < 1.0e-9
    assert max(abs(a - b) for a, b in zip(zero_u_moments, no_u_moments)) < 1.0e-9
    assert "no density matrix found ... skipping LDA+U" in Path(zero_u["stdout"]).read_text()


@pytest.mark.fleur
@pytest.mark.noco
@pytest.mark.hdf
def test_ffn_directional_constraint_reaches_hamiltonian(execute_fleur, work_dir):
    """The transverse constraining field must survive vgen_finalize in FullyFullyNoco.

    For l_mtNocoPot=T, vgen_finalize calls rotate_mt_den_from_local, which zeroes
    vTot%mt and rebuilds all four spin components from the local-frame diagonal pair
    plus theta_mt/phi_mt. vgen_constraint therefore has to be applied *after* that
    step; when it ran before it, l_constrained=T produced output bit-identical to
    l_constrained=F. This test fails (on/off indistinguishable) without that fix.
    """
    import math
    import shutil

    input_dir = "./inputfiles/noco/FFNConstraint"

    def misalignments(outxml):
        """|m_perp|/|m| in the local frame, one entry per SCF iteration."""
        root = ElementTree.parse(outxml).getroot()
        out = []
        for element in root.findall(".//magneticMomentsInMTSpheres/localMagMoment"):
            mx, my, mz = (float(v) for v in element.attrib["vec"].split())
            out.append(math.copysign(math.hypot(mx, my), mx + my) / math.sqrt(mx*mx + my*my + mz*mz))
        return out

    def seeded_run(tag, seed_inp, run_inp):
        """Build the tilted start density, then run the SCF from it."""
        Path(work_dir, f"{tag}_seed").mkdir()
        Path(work_dir, tag).mkdir()
        seed = execute_fleur(input_dir, only_copy=[[seed_inp, "inp.xml"]],
                             sub_dir=f"{tag}_seed", mpi_procs=1)
        shutil.copy(Path(seed["cdn.hdf"]), Path(work_dir, tag, "cdn.hdf"))
        return execute_fleur(input_dir, only_copy=[[run_inp, "inp.xml"]],
                             sub_dir=tag, mpi_procs=1)

    on_p = seeded_run("on_p", "inp_seed_p.xml", "inp_ffn_on.xml")
    off_p = seeded_run("off_p", "inp_seed_p.xml", "inp_ffn_off.xml")
    on_m = seeded_run("on_m", "inp_seed_m.xml", "inp_ffn_on.xml")

    mis_on, mis_off, mis_on_m = (misalignments(r["out.xml"]) for r in (on_p, off_p, on_m))

    assert len(mis_on) >= 3 and len(mis_on) == len(mis_off) == len(mis_on_m)

    # both runs start from the same density, so iteration 1 must agree
    assert abs(mis_on[0] - mis_off[0]) < 1.0e-9

    # the constraint must actually act: this is what failed before the fix
    assert abs(mis_on[-1] - mis_off[-1]) > 1.0e-4, \
        "l_constrained=T is indistinguishable from l_constrained=F in FFN mode"

    # ... and it must be restoring, monotonically
    assert all(abs(b) < abs(a) for a, b in zip(mis_on, mis_on[1:])), \
        f"misalignment is not decreasing under the constraint: {mis_on}"
    assert abs(mis_on[-1]) < 0.75 * abs(mis_off[-1])

    # opposite initial deviation -> mirrored response, same magnitude
    assert mis_on[0] * mis_on_m[0] < 0.0
    assert max(abs(abs(a) - abs(b)) for a, b in zip(mis_on, mis_on_m)) < 1.0e-6


@pytest.mark.fleur
@pytest.mark.noco
@pytest.mark.hdf
def test_mt_collinear_constraint_path_unchanged(execute_fleur, work_dir):
    """The non-FFN (l_mtNocoPot=F) constrained path must be unaffected by the fix.

    vgen_constraint was moved from before to after vgen_finalize and made additive.
    For l_mtNocoPot=F, components 3/4 are still zero at the new insertion point, so
    the result must be unchanged. Here the local moment is aligned with the local
    axis, so the constraining field is exactly zero and switching l_constrained on
    must change nothing at all -- this catches a relocation that injects a spurious
    transverse field into the collinear-MT path.
    """
    input_dir = "./inputfiles/noco/FFNConstraint"

    # no seed density here: the self-generated start density is collinear along the
    # local axis, which is exactly the aligned reference state this test needs.
    on = execute_fleur(input_dir, only_copy=[["inp_mtc_on.xml", "inp.xml"]],
                       mpi_procs=1)

    def moments(outxml):
        root = ElementTree.parse(outxml).getroot()
        return [float(v)
                for e in root.findall(".//magneticMomentsInMTSpheres/localMagMoment")
                for v in e.attrib["vec"].split()]

    on_m = moments(on["out.xml"])
    assert on_m, "no moments written by the non-FFN constrained run"
    assert all(v == v for v in on_m), "NaN in the non-FFN constrained moments"

    # the constraint bookkeeping must run (element present) and be identically zero
    root = ElementTree.parse(on["out.xml"]).getroot()
    fields = root.findall(".//constrainingField")
    assert fields, "constrainingField block missing from the non-FFN constrained run"
    deltas = [float(e.attrib["DBConX"]) for e in root.iter()
              if e.tag.startswith("DeltaBConstrain")]
    assert deltas, "no DeltaBConstrain entries recorded"
    assert max(abs(d) for d in deltas) < 1.0e-8, \
        f"spurious transverse constraining field in the collinear-MT path: {deltas}"
