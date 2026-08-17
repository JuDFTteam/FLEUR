
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
