
import pytest

@pytest.mark.lo
@pytest.mark.xml
@pytest.mark.serial
def test_SiLOXML(execute_fleur, grep_number, grep_exists):
    """Fleur Si with LOs - XML

    Simple test of Fleur with XML input with one step:
    1.Generate a starting density and run a few iteration and compare fermi-energy & total energy
    """
    test_file_folder = './inputfiles/SiLOXML/'
    res_files = execute_fleur(test_file_folder)
    res_file_names = list(res_files.keys())
    should_files = ['out']
    for file1 in should_files:
        assert (file1 in res_file_names), f'{file1} missing'

    assert grep_exists(res_files['out'], "it=  6  is completed")
    efermi = grep_number(res_files['out'], "first approx. to ef", ":")
    bandgap = grep_number(res_files['out'], "bandgap", ":")
    tenergy = grep_number(res_files['out'], "    total energy=", "=")
    dist = grep_number(res_files['out'], "distance of charge densities for spin  1                 it=    6", ":")

    assert abs(efermi - 0.18481) <= 0.0001
    assert abs(bandgap - 0.08353) <= 0.0001
    assert abs(tenergy - -580.064565) <= 0.0001
    assert abs(dist - 0.0261) <= 0.0001
