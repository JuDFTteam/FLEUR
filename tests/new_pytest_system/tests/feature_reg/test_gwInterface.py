import pytest

@pytest.mark.noci
@pytest.mark.gw
@pytest.mark.hdf
def test_gw1Interface(execute_fleur, grep_number, grep_exists):
    """GW=1 switch for Si

    Simple test of Fleur with one step:
    1.Generate a starting density and run 1 iteration with gw=1 and compare several quantities written out
    """
    test_file_folder = './inputfiles/gw1Interface/'
    res_files = execute_fleur(test_file_folder)
    res_file_names = list(res_files.keys())
    should_files = ['out', 'basis.hdf', 'pot.hdf', 'ecore', 'sym.out']
    for file1 in should_files:
        assert (file1 in res_file_names), f'{file1} missing'
    
    assert grep_exists(res_files['out'], "it=  1  is completed")


@pytest.mark.noci
@pytest.mark.gw
@pytest.mark.hdf
def test_gw2Interface(execute_fleur, grep_number, grep_exists):
    """GW=2 switch for Si

    Simple test of Fleur with one step:
    1.Generate a starting density and run 1 iteration with gw=1 and compare several quantities written out
    """
    test_file_folder = './inputfiles/gw2Interface/'
    res_files = execute_fleur(test_file_folder)
    res_file_names = list(res_files.keys())
    should_files = ['out', 'eig_gw.hdf']
    for file1 in should_files:
        assert (file1 in res_file_names), f'{file1} missing'
    
    assert grep_exists(res_files['out'], "0.50000.*0.37500")

