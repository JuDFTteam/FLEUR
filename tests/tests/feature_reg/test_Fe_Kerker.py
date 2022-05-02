
import pytest

@pytest.mark.serial
def test_Fe_Kerker(execute_fleur, grep_exists, grep_number):
    """Fleur Fe Kerker XML

    Test of the Kerker preconditioner with one step:
    1.Generate a starting density and run 3 iterations and compare distance
    """
    test_file_folder = './inputfiles/Fe_Kerker/'
    res_files = execute_fleur(test_file_folder)
    res_file_names = list(res_files.keys())
    should_files = ['out']
    for file1 in should_files:
        assert file1 in res_file_names

    assert grep_exists(res_files['out'], "it=  3  is completed")
    dist = grep_number(res_files['out'], "distance of charge densities for it=    3", ":")

    assert abs(dist - 11.745805) <= 0.001
