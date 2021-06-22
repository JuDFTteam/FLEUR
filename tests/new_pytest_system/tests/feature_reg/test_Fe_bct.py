
import pytest

@pytest.mark.bulk
@pytest.mark.magnetism
@pytest.mark.libxc
def test_Fe_bct_LibXC(execute_fleur, grep_number, grep_exists):
    """Fleur Fe bct non-collinear XML

    Simple test of Fleur with two steps:
    1.Generate a starting density
    2.Run 20 iterations and compare convergence, fermi-energy & total energy
    """
    test_file_folder = './inputfiles/Fe_bct_LibXC/'
    res_files = execute_fleur(test_file_folder)
    should_files = ['out']
    res_file_names = list(res_files.keys())
    if 'cdn.hdf' not in res_file_names:
        should_files.append('cdn1')
    for file1 in should_files:
        assert file1 in res_file_names, f'{file1} missing'

    assert grep_exists(res_files['out'], "total charge")
    qfix = grep_number(res_files['out'], "qfix=", "qfix=")
    assert grep_exists(res_files['out'], "it=  3  is completed")
    fermi = grep_number(res_files['out'], "new fermi energy", ":")
    tenergy = grep_number(res_files['out'], "total energy=", "=")
    dist = grep_number(res_files['out'], "distance of charge densities for it=    3", ":")
    mm = grep_number(res_files['out'], "mm       2", "mm       2")

    assert abs(qfix - 1.0) <= 0.00001
    assert abs(fermi - 0.338) <= 0.005
    assert abs(tenergy - -2545.600) <= 0.005
    assert abs(dist - 12.700) <= 0.09
    assert abs(mm - 1.96) <= 0.03

@pytest.mark.serial
@pytest.mark.bulk
@pytest.mark.magnetism
@pytest.mark.xml
def test_Fe_bct_LOXML(execute_fleur, grep_number, grep_exists):
    """Fleur Fe bct noco LOs XML

    Simple test of Fleur with two steps:
    1.Generate a starting density (+get some errors afterwards)
    2.Run 4 noco iterations and compare convergence, fermi-energy & total energy
    """
    test_file_folder = './inputfiles/Fe_bct_LOXML/'
    res_files = execute_fleur(test_file_folder)
    should_files = ['out']
    res_file_names = list(res_files.keys())
    if 'cdn.hdf' not in res_file_names:
        should_files.append('cdn1')
    for file1 in should_files:
        assert file1 in res_file_names

    assert grep_exists(res_files['out'], "total charge")
    qfix = grep_number(res_files['out'], "qfix=", "qfix=")
    assert grep_exists(res_files['out'], "it=  4  is completed")
    fermi = grep_number(res_files['out'], "new fermi energy", ":")
    tenergy = grep_number(res_files['out'], "total energy=", "=")
    dist = grep_number(res_files['out'], "distance of charge densities for it=", ":")
    mm = grep_number(res_files['out'], "mm       2", "mm       2")

    assert abs(qfix - 1.0) <= 0.00001
    assert abs(fermi - 0.341) <= 0.01
    assert abs(tenergy - -2545.5968) <= 0.01
    assert abs(dist - 1.593616) <= 0.3
    assert abs(mm - 1.90) <= 0.01

@pytest.mark.serial
@pytest.mark.bulk
@pytest.mark.magnetism
@pytest.mark.soc
def test_Fe_bct_SOCXML(execute_fleur, grep_number, grep_exists):
    """Fleur Fe bct noco SOC XML

    Simple test of Fleur with two steps:
    1.Generate a starting density
    2.Run 20 iterations and compare convergence, fermi-energy & total energy
    """
    test_file_folder = './inputfiles/Fe_bct_SOCXML/'
    res_files = execute_fleur(test_file_folder)
    should_files = ['out']
    res_file_names = list(res_files.keys())
    if 'cdn.hdf' not in res_file_names:
        should_files.append('cdn1')
    for file1 in should_files:
        assert file1 in res_file_names

    assert grep_exists(res_files['out'], "total charge")
    qfix = grep_number(res_files['out'], "qfix=", "qfix=", res_index=-1)
    assert grep_exists(res_files['out'], "it= 20  is completed")
    fermi = grep_number(res_files['out'], "new fermi energy", ":")
    tenergy = grep_number(res_files['out'], "total energy=", "=")
    dist = grep_number(res_files['out'], "distance of charge densities for it=   20", ":")
    mm = grep_number(res_files['out'], "mm       2", "mm       2")

    assert abs(qfix - 1.0) <= 0.00001
    assert abs(fermi - 0.326) <= 0.005
    assert abs(tenergy - -2545.607623611) <= 0.005
    assert abs(dist - 0.00003) <= 0.1
    assert abs(mm - 0.05411) <= 0.01

@pytest.mark.serial
@pytest.mark.bulk
@pytest.mark.magnetism
@pytest.mark.xml
def test_Fe_bctXML(execute_fleur, grep_number, grep_exists):
    """Fleur Fe bct non-collinear XML

    Simple test of Fleur with two steps:
    1.Generate a starting density
    2.Run 20 iterations and compare convergence, fermi-energy & total energy
    """
    test_file_folder = './inputfiles/Fe_bctXML/'

    res_files = execute_fleur(test_file_folder)
    should_files = ['out']
    res_file_names = list(res_files.keys())
    if 'cdn.hdf' not in res_file_names:
        should_files.append('cdn1')
    for file1 in should_files:
        assert file1 in res_file_names

    assert grep_exists(res_files['out'], "total charge")
    qfix = grep_number(res_files['out'], "qfix=", "qfix=")
    assert grep_exists(res_files['out'], "it= 20  is completed")
    fermi = grep_number(res_files['out'], "new fermi energy", ":")
    tenergy = grep_number(res_files['out'], "total energy=", "=")
    dist = grep_number(res_files['out'], "distance of charge densities for it=   20", ":")
    mm = grep_number(res_files['out'], "mm       2", "mm       2")

    assert abs(qfix - 1.0) <= 0.00001 #, 'qfix'
    assert abs(fermi - 0.326) <= 0.005 #, 'Fermi energy'
    assert abs(tenergy - -2545.607) <= 0.005 #, 'Total energy'
    assert abs(dist - 0.00003) <= 0.09 #, 'Distance'
    assert abs(mm - 1.73) <= 0.03 #, 'Magnetic moment 2'
