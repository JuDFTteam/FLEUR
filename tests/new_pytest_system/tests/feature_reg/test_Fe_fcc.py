
import pytest

@pytest.mark.serial
@pytest.mark.bulk
@pytest.mark.xml
@pytest.mark.magnetism
def test_Fe_fccXML(execute_fleur, grep_number, grep_exists):
    """Fleur Fe fcc spin-spiral XML

    Simple test of Fleur with two steps:
    1.Generate a starting density
    2.Run 20 iterations and compare convergence, fermi-energy & total energy

    "SOC",0,"complex",1,"EVP",0
    """
    test_file_folder = './inputfiles/Fe_fccXML/'
    res_files = execute_fleur(test_file_folder)
    should_files = ['out']
    res_file_names = list(res_files.keys())
    if 'cdn.hdf' not in res_file_names:
        should_files.append('cdn1')
    for file1 in should_files:
        assert file1 in res_file_names

    assert grep_exists(res_files['out'], "total charge")
    qfix = grep_number(res_files['out'], "qfix=", "qfix=")
#    assert grep_exists(res_files['out'], "it= 20  is completed")

    fermi = grep_number(res_files['out'], "to ef", ":")
    tenergy = grep_number(res_files['out'], "total energy=", "=")
    dist = grep_number(res_files['out'], "distance of charge densities for it=", ":")
    mm = grep_number(res_files['out'], "mm       1", "mm       1")

    assert abs(qfix - 1.0) <= 0.00001
    assert abs(fermi - 0.322013) <= 0.01
    assert abs(tenergy - -1273.0878841822) <= 0.00002
    assert abs(dist - 0.0001) <= 0.01
    assert abs(mm - 3.93483) <= 0.001

@pytest.mark.magnetism
@pytest.mark.soc
@pytest.mark.lo
@pytest.mark.hdf
def test_FeFFNLOsSOC(execute_fleur, grep_number, grep_exists):
    """Iron LO's and SOC test in FFN

    Simple test of Fleur with one step:
    1.Generate a starting density and run 1 iteration and check magnetization of Fe. 
    """
    test_file_folder = './inputfiles/FeFFNLOsSOC/'
    res_files = execute_fleur(test_file_folder)
    should_files = ['out']
    res_file_names = list(res_files.keys())
    for file1 in should_files:
        assert file1 in res_file_names

    assert grep_exists(res_files['out'], r"mx= [\s\-]0\.00000 my= [\s\-]0\.00000 mz=  2\.293")
