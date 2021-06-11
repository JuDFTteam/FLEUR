
import pytest

@pytest.mark.fleur
@pytest.mark.bulk
def test_CrystalFieldOutput(execute_fleur, grep_exists, grep_number):
    """Fleur Crystalfield Output
    Simple test of Fleur crystalfiedl output:
    1.Generate a starting density, run 1 iteration, generate both potential and charge density for
    crystalfield calculations and check that the CFdata.hdf file was created
    """
    test_file_folder = './inputfiles/CrystalFieldOutput/'
    res_files = execute_fleur(test_file_folder,
                              only_copy=['inp.xml', 'sym.xml'])
    res_file_names = list(res_files.keys())
    should_files = ['out']

    if 'cdn.hdf' in res_file_names:
        should_files.append('CFdata.hdf')
    else:
        should_files = should_files + ['n4f.1.dat', 'V_20.1.dat', 'V_40.1.dat', 'V_60.1.dat', 
        'V_66.1.dat']

    for file1 in should_files:
        assert file1 in res_file_names, f'{file1} missing'

    #Test that there was a second vgen call with modified potential
    #Here the big multipole moments are grepped (Maybe not the best way)
    # The python grep has interprets the E+01 as some expression
    assert grep_exists(res_files['out'], "0.14060E")#+01")
    assert grep_exists(res_files['out'], "0.82207E")#+00")
    assert grep_exists(res_files['out'], "0.89945E")#+00")
    assert grep_exists(res_files['out'], "0.65772E")#+01")