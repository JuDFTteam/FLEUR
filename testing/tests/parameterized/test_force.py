
import pytest
"""
Boiler-plate code for testset
"""
from read_tests import read_tests
all_tests = read_tests("forces")

@pytest.mark.fleur
@pytest.mark.forces
@pytest.mark.parametrize(("dir","desc"), all_tests)
def test_forces(dir,desc,default_fleur_test):
    """
    """
    
    assert default_fleur_test(dir)
    