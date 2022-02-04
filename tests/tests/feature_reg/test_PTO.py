
import pytest

@pytest.mark.serial
@pytest.mark.lo
@pytest.mark.soc
@pytest.mark.xml
@pytest.mark.magnetism
def test_PTO_SOCXML(execute_fleur, grep_number, grep_exists):
    """PTO one iteration with SOC (XML codepath)

    Simple test of Fleur with one step:
    1.Generate a starting density and run a single iteration and test for ef, total energy
    """
    test_file_folder = './inputfiles/PTO-SOCXML/'
    res_files = execute_fleur(test_file_folder, only_copy=['inp.xml'])
    # there are other files in the folder, which are not used, but would change the result of the
    # bandgap
    res_file_names = list(res_files.keys())
    should_files = ['out']
    for file1 in should_files:
        assert (file1 in res_file_names), f'{file1} missing'

    assert grep_exists(res_files['out'], "it=  1  is completed")
    bandgap = grep_number(res_files['out'], "bandgap", ":")
    tenergy = grep_number(res_files['out'], "total energy=", "=")

    assert abs(bandgap - 0.0470) <= 0.001
    assert abs(tenergy - -22006.622) <= 0.01

@pytest.mark.serial
@pytest.mark.lo
@pytest.mark.xml
@pytest.mark.magnetism
def test_PTOXML(execute_fleur, grep_number, grep_exists):
    """PTO one iteration with LOs (XML codepath)

    Simple test of Fleur with one step:
    1.Generate a starting density and run a single iteration and test for ef, total energy
    """
    test_file_folder = './inputfiles/PTOXML/'
    res_files = execute_fleur(test_file_folder)
    res_file_names = list(res_files.keys())
    should_files = ['out']
    for file1 in should_files:
        assert (file1 in res_file_names), f'{file1} missing'

    assert grep_exists(res_files['out'], "it=  1  is completed")
    efermi = grep_number(res_files['out'], "new fermi", ":")
    tenergy = grep_number(res_files['out'], "total energy=", "=")

    assert abs(efermi - 0.211258) <= 0.001
    assert abs(tenergy - -22006.6222524) <= 0.0001
