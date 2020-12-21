
import pytest

@pytest.mark.fleur
@pytest.mark.bulk
@pytest.mark.xml
@pytest.mark.fast
def test_CuBulkXML(execute_fleur, grep_exists, grep_number, stage_for_parser_test):
    """
    Simple test of Fleur with XML input with one step:
    1.Generate a starting density and run 1 iteration and compare fermi-energy & total energy
    """
    test_file_folder = './inputfiles/CuBulkXML/'

    cmd_params = []
    res_files = execute_fleur(cmd_params, test_file_folder)
    should_files = ['out.xml', 'out', 'usage.json']
    res_file_names = list(res_files.keys())
    
    # Test if all files are there
    for file1 in should_files:
        assert file1 in res_file_names
    
    assert grep_exists(res_files['out'], "it=  1  is completed")
    fermi = grep_number(res_files['out'], "new fermi energy", ":")
    tenergy = grep_number(res_files['out'], "total energy=", "=")
    dist = grep_number(res_files['out'], "distance of charge densitie", "1:")
    
    assert abs(fermi - 0.4233) <= 0.001
    assert abs(tenergy - -3305.016) <= 0.001
    assert abs(dist - 45.8259) <= 0.001


@pytest.mark.fleur
@pytest.mark.bulk
@pytest.mark.xml
@pytest.mark.band
@pytest.mark.fast
def test_CuBandXML(execute_fleur, grep_exists, stage_for_parser_test):
    """
    Simple test of Fleur band structure calculation with XML input with one step:
    1.Generate a starting density, run 1 iteration, and generate band structure. Ensure that the files are created.
    """
    test_file_folder = './inputfiles/CuBandXML/'
    cmd_params = []
    res_files = execute_fleur(cmd_params, test_file_folder)
    res_file_names = list(res_files.keys())
    should_files = ['out.xml', 'bands.1']
    if 'cdn.hdf' in res_file_names:
        with_hdf = True
    else:
        with_hdf = False

    # Test if all files are there
    for file1 in should_files:
        assert file1 in res_file_names

    # Note: For the test the correction of the Fermi energy with HDF5 is wrong because
    #       the band structure is calculated directly on top of the starting
    #       density (previous eFermi = 0.0).
    bands_file = res_files['bands.1']
    if with_hdf:
        assert grep_exists(bands_file, "0.033629.*-0.400779") #line 2
        assert grep_exists(bands_file, "0.269033.*0.614193") #line 9
        assert grep_exists(bands_file, "0.336291.*32.833232") #line 351
        assert grep_exists(bands_file, "0.638954.*33.485766") #line 360
    else:
        assert grep_exists(bands_file, "0.033629.*-11.565510") #line 2
        assert grep_exists(bands_file, "0.269033.*-10.550538") #line 9
        assert grep_exists(bands_file, "0.336291.*21.668500") #line 351
        assert grep_exists(bands_file, "0.638954.*22.321034") #line 360


@pytest.mark.fleur
@pytest.mark.bulk
@pytest.mark.libxc
@pytest.mark.fast
def test_CuBulkLibXC(execute_fleur, grep_exists, grep_number):
    """Simple test of Fleur with XML input with one step:
    1.Generate a starting density and run 1 iteration and compare fermi-energy & total energy
    
    If fleur is not complied and linked to libXC this test will fail
    """
    test_file_folder = './inputfiles/CuBulkLibXC/'
    cmd_params = []
    res_files = execute_fleur(cmd_params, test_file_folder)
    res_file_names = list(res_files.keys())
    should_files = ['out']
    for file1 in should_files:
        assert file1 in res_file_names
    
    assert grep_exists(res_files['out'], "it=  1  is completed")
    fermi = grep_number(res_files['out'], "new fermi energy", ":")
    tenergy = grep_number(res_files['out'], "total energy=", "=")
    dist = grep_number(res_files['out'], "distance of charge densitie", "1:")

    assert abs(fermi - 0.4250) <= 0.001
    assert abs(tenergy - -3305.007) <= 0.001
    assert abs(dist - 48.986892) <= 0.001


@pytest.mark.fleur
@pytest.mark.bulk
@pytest.mark.xml
@pytest.mark.dos
@pytest.mark.fast
def test_CuDOSXML(execute_fleur, grep_exists, grep_number, stage_for_parser_test):
    """Simple test of Fleur DOS calculation with XML input with one step:
    1.Generate a starting density, run 1 iteration, and generate DOS. Ensure that the files are created.

    """
    test_file_folder = './inputfiles/CuDOSXML/'
    cmd_params = []
    res_files = execute_fleur(cmd_params, test_file_folder)
    res_file_names = list(res_files.keys())
    should_files = ['out']
    if 'cdn.hdf' in res_file_names:
        with_hdf = True
        should_files.append('banddos.hdf')
    else:
        with_hdf = False
        should_files.append('Local.1')
    for file1 in should_files:
        assert file1 in res_file_names

    fermi = grep_number(res_files['out'], "fermi energy =", "=")
    assert abs(fermi - 0.4239) <= 0.001
    
    if not with_hdf:
        assert grep_exists(res_files['Local.1'], "0.207025")
        assert grep_exists(res_files['Local.1'], "0.326256")


@pytest.mark.fleur
@pytest.mark.bulk
@pytest.mark.xml
@pytest.mark.dos
@pytest.mark.fast
def test_CuOrb(execute_fleur, grep_exists, grep_number):
    """Simple test of Fleur DOS calculation with XML input with one step:
    1.Generate a starting density, run 1 iteration, and generate DOS. Ensure that the files are created.
    """
    # TODO these are the same checks as in the CuDOSXML test, is this correct?
    # we do not check orbital decom...
    test_file_folder = './inputfiles/CuOrb/'
    cmd_params = []
    res_files = execute_fleur(cmd_params, test_file_folder)
    res_file_names = list(res_files.keys())
    should_files = ['out']
    if 'cdn.hdf' in res_file_names:
        with_hdf = True
        should_files.append('banddos.hdf')
    else:
        with_hdf = False
        should_files.append('Local.1')
    for file1 in should_files:
        assert file1 in res_file_names

    fermi = grep_number(res_files['out'], "fermi energy =", "=")
    assert abs(fermi - 0.4239) <= 0.001
    
    if not with_hdf:
        assert grep_exists(res_files['Local.1'], "0.207025")
        assert grep_exists(res_files['Local.1'], "0.326256")
