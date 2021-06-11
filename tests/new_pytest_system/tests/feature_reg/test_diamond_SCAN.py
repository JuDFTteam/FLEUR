
import pytest

@pytest.mark.libxc
@pytest.mark.skip('We have to investigate LibXC Meta-GGAs before we reenable this.')
def test_diamond_SCAN(execute_fleur, grep_number, grep_exists):
    """Fleur Diamond SCAN
    Converge diamond for one k-point with scan
    
    # Does this even Fail if scan is missing?
    """
    test_file_folder = './inputfiles/Diamond_SCAN/'
    
    res_files = execute_fleur(test_file_folder)
    res_file_names = list(res_files.keys())
    should_files = ['out']
    for file1 in should_files:
        assert (file1 in res_file_names), file1

    assert grep_exists(res_files['out'], "it= 12  is completed")
    tenergy = grep_number(res_files['out'], "total energy=", "=")
    eff_pot = grep_number(res_files['out'], "density-effective potential i", "l i")
    dist = grep_number(res_files['out'], "charge density-ex.-corr.energ", "energ ")
    
    assert abs(tenergy - -75.1367097131) <= 0.001
    assert abs(eff_pot - -111.6343435080) <= 0.001
    assert abs(dist - -11.0673799988) <= 0.001