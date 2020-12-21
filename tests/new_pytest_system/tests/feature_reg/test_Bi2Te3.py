
import pytest

@pytest.mark.fleur
@pytest.mark.bulk
@pytest.mark.soc
@pytest.mark.lo
def test_Bi2Te3(execute_fleur, grep_exists, grep_number):
    """Bi2Te3 with magnetisation
    Simple test of Fleur with two steps:
    1.Generate a starting density
    2.Run a single iteration and test for ef, total energy and orbital moment
    Uses: SOC,LOs
    """
    test_file_folder = './inputfiles/Bi2Te3/'
    cmd_params = []

    # Stage 1
    res_files = execute_fleur(cmd_params, test_file_folder, only_copy=['inp', 'kpts'])
    should_files = ['out']
    res_file_names = list(res_files.keys())
    if 'cdn.hdf' in res_file_names:
        with_hdf = True
        should_files.append('cdn.hdf')
    else:
        with_hdf = False
        should_files.append('cdn1')
    for file1 in should_files:
        assert file1 in res_file_names
    
    assert grep_exists(res_files['out'], "total charge")
    qfix = grep_number(res_files['out'], "qfix=", "qfix= ")
    assert abs(qfix - 1.0) <= 0.00001
    
    # Stage 2
    res_files = execute_fleur(cmd_params, test_file_folder, only_copy=[{'inp2': 'inp'}, 'kpts'])
    res_file_names = list(res_files.keys())
    should_files = ['out']
    for file1 in should_files:
        assert file1 in res_file_names
    
    assert grep_exists(res_files['out'], "it=  1  is completed")
    fermi = grep_number(res_files['out'], "first approx. to ef", ":")
    tenergy = grep_number(res_files['out'], "total energy=", "=")
    mm = grep_number(res_files['out'], "mm      15    ", "mm      15 ")

    assert abs(fermi - 0.188) <= 0.001
    assert abs(tenergy - -190662.67147) <= 0.01
    assert abs(mm - -0.022) <= 0.01

@pytest.mark.fleur
@pytest.mark.bulk
@pytest.mark.soc
@pytest.mark.lo
@pytest.mark.xml
def test_Bi2Te3XML(execute_fleur, grep_exists, grep_number):
    """Bi2Te3 with magnetisation (XML codepath)
    Simple test of Fleur with one step:
    1.Generate a starting density and run a single iteration and test for ef, total energy and orbital moment
    Uses: SOC,LOs
    """
    test_file_folder = './inputfiles/Bi2Te3XML/'
    cmd_params = []

    # Stage 1
    res_files = execute_fleur(cmd_params, test_file_folder)
    res_file_names = list(res_files.keys())
    should_files = ['out']
    for file1 in should_files:
        assert file1 in res_file_names
    
    assert grep_exists(res_files['out'], "it=  1  is completed")
    fermi = grep_number(res_files['out'], "first approx. to ef", ":")
    tenergy = grep_number(res_files['out'], "total energy=", "=")
    mm = grep_number(res_files['out'], "mm      15    ", "mm      15 ")

    assert abs(fermi - 0.188) <= 0.001
    assert abs(tenergy - -190662.68039) <= 0.01
    assert abs(mm - -0.02) <= 0.01