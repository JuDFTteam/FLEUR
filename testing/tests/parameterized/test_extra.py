
import pytest
"""
Boiler-plate code for testset
"""
from read_tests import read_tests
all_tests = read_tests("extra")

@pytest.mark.fleur
@pytest.mark.extra
@pytest.mark.parametrize(("dir","desc","cmdline","mpi_procs"), all_tests)
def test_extra(dir,desc,cmdline,mpi_procs,default_fleur_test):
    """
    """

    assert default_fleur_test(dir, cmdline_args=cmdline, mpi_procs=mpi_procs)
    