
import pytest

#Here noco tests can be simply added by:
# creating a new entry in the list
# first string is the directory of the input
# second string is the short descriptions
# in addition marks can be set.
# it is also a good idea to add an explicit id to avoid complex directory names
all_tests = [
   
]

@pytest.mark.serial
@pytest.mark.fleur
@pytest.mark.bulk
@pytest.mark.noco
@pytest.mark.parametrize(("dir","desc"), all_tests)
def test_noco(dir,desc,default_fleur_test):
    """
    """
    
    assert default_fleur_test(dir)
    