
import pytest
"""
Boiler-plate code for testset
"""
from read_tests import read_tests
all_tests = read_tests("basic")

@pytest.mark.fleur
@pytest.mark.parametrize(("dir","desc","cmdline","mpi_procs"), all_tests)
def test_basic(dir,desc,cmdline,mpi_procs,default_fleur_test):
    """
    """

    assert default_fleur_test(dir, cmdline_args=cmdline, mpi_procs=mpi_procs)

@pytest.mark.fleur
@pytest.mark.bulk
@pytest.mark.dos
@pytest.mark.fast
def test_CuDM(default_fleur_test, grep_number):
    """
    Band-resolved density matrix (symmetrized, s, p, d and all l-l' blocks):
    the traces of the diagonal blocks must reproduce the l-resolved DOS weights
    """
    res_files = default_fleur_test("basic/CuDM", mpi_procs=2)
    deviation = grep_number(res_files['out'], "from the l-resolved DOS weights:")
    assert deviation < 1e-10
