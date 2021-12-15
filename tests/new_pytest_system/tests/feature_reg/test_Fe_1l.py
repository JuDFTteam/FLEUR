
import pytest

@pytest.mark.serial
@pytest.mark.soc
@pytest.mark.film
@pytest.mark.xml
def test_Fe_1l_SOCXML(execute_fleur, grep_number, grep_exists):
    """Fleur Fe Monolayer SOC XML
    Simple test of Fleur with one steps:
    1.Generate a starting density and run 1 iteration and compare convergence, fermi-energy & total energy
    """
    test_file_folder = './inputfiles/Fe_1l_SOCXML/'

    res_files = execute_fleur(test_file_folder, only_copy=['inp.xml'])
    # there are also other files in the test folder, which are somehow not needed
    # and the use of these leads to a different result for mm
    should_files = ['out']
    res_file_names = list(res_files.keys())
    for file1 in should_files:
        assert file1 in res_file_names

    assert grep_exists(res_files['out'], "it=  1  is completed")
    mm = grep_number(res_files['out'], "mm       1", " 1 ")
    qfix = grep_number(res_files['out'], "qfix=", "x=")
    fermi = grep_number(res_files['out'], "new fermi energy", ":")
    tenergy = grep_number(res_files['out'], "total energy=", "=")

    assert abs(mm - 0.279) <= 0.001
    assert abs(qfix - 1.0) <= 0.0001
    assert abs(fermi - -0.2450) <= 0.0001
    assert abs(tenergy - -1272.6885) <= 0.001

@pytest.mark.serial
@pytest.mark.film
def test_Fe_1l_Tria(execute_fleur, grep_number, grep_exists):
    """Fleur Fe Monolayer Triangular method
    Simple test of Fleur with XML input with one step:
    1.Generate a starting density and run a single iteration and compare convergence, fermi-energy & total energy
    (with linear triangular method for fermi energy evaluation)
    """
    test_file_folder = './inputfiles/Fe_1l_Tria/'

    res_files = execute_fleur(test_file_folder)
    should_files = ['out']
    res_file_names = list(res_files.keys())
    for file1 in should_files:
        assert file1 in res_file_names

    assert grep_exists(res_files['out'], "it=  5  is completed")
    fermi = grep_number(res_files['out'], "fermi energy=", "=")
    tenergy = grep_number(res_files['out'], "total energy=", "=")
    dist = grep_number(res_files['out'], "distance of charge densities for it=", ":")


    assert abs(fermi - -0.16919) <= 0.005
    assert abs(tenergy - -1272.6376271846) <= 0.01
    assert abs(dist - 9.75) <= 0.5


@pytest.mark.serial
@pytest.mark.film
@pytest.mark.xml
def test_Fe_1lXML(execute_fleur, grep_number, grep_exists):
    """Fleur Fe Monolayer
    Simple test of Fleur with XML input with one step:
    1.Generate a starting density and run a single iteration and compare convergence, fermi-energy & total energy
    """
    test_file_folder = './inputfiles/Fe_1lXML/'

    res_files = execute_fleur(test_file_folder)
    should_files = ['out']
    res_file_names = list(res_files.keys())
    for file1 in should_files:
        assert file1 in res_file_names

    assert grep_exists(res_files['out'], "it=  1  is completed")
    fermi = grep_number(res_files['out'], "new fermi energy", ":")
    tenergy = grep_number(res_files['out'], "total energy=", "=")
    dist = grep_number(res_files['out'], "distance of charge densitie", "1:")

    assert abs(fermi - -0.242) <= 0.005
    assert abs(tenergy - -1272.68) <= 0.01
    assert abs(dist - 19.5) <= 0.5

