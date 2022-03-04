
import pytest

@pytest.mark.serial
@pytest.mark.fleur
@pytest.mark.bulk
@pytest.mark.xml
@pytest.mark.dos
@pytest.mark.mcd
def test_CoMCDXML(execute_fleur, grep_exists, grep_number):
    """Fleur Co MCD XML
    Simple test of Fleur with two steps:
    1.Generate a starting density and run 2 iterations
    2.Calculate and verify MCD spectra
    """
    test_file_folder = './inputfiles/CoMCDXML/'

    # Stage 1
    res_files = execute_fleur(test_file_folder, only_copy=['inp.xml', 'sym.out'])
    res_file_names = list(res_files.keys())
    should_files = ['out']
    with_hdf = False
    if 'cdn.hdf' in res_file_names:
        with_hdf = True
    else:
        should_files.append('cdn1')
    for file1 in should_files:
        assert file1 in res_file_names

    tenergy = grep_number(res_files['out'], "total energy=", "=")
    dist = grep_number(res_files['out'], "distance of charge densities for it=    2", ": ")
    dist_spin = grep_number(res_files['out'], "distance of spin densities for it=    2", ": ")
    
    assert abs(tenergy - -2786.95650275) <= 0.001
    assert abs(dist - 8.376682) <= 0.001
    assert abs(dist_spin - 17.14035) <= 0.001


    # Stage 2
    res_files = execute_fleur(test_file_folder, only_copy=[['inp2.xml', 'inp.xml']])
    res_file_names = list(res_files.keys())
    should_files = []
    if with_hdf:
        should_files.append('banddos.hdf')
    else:
        should_files.append('MCD.1')
        
    for file1 in should_files:
        assert file1 in res_file_names
    
    if not with_hdf:
        assert grep_exists(res_files['MCD.1'], "0.95976")
        assert grep_exists(res_files['MCD.1'], "0.36559")


@pytest.mark.fleur
@pytest.mark.bulk
@pytest.mark.band
def test_CoUnfold(execute_fleur, grep_exists, grep_number):
    """Co band unfolding test
    Simple test of Fleur with one step:
    1.Generate a starting density and run 1 iteration with band unfolding and compare several quantities
    """
    test_file_folder = './inputfiles/CoUnfold/'

    res_files = execute_fleur(test_file_folder, only_copy=['inp.xml', 'sym.xml'], mpi_procs=3)
    res_file_names = list(res_files.keys())
    should_files = ['out', 'bands_sc.1', 'bands_sc.2']
    for file1 in should_files:
        assert file1 in res_file_names

    if 'cdn.hdf' in res_file_names:
        with_hdf = True
    else:
        with_hdf = False

    if with_hdf:
		#    assert grep_exists(res_files['out'], "it=  1  is completed")
        assert grep_exists(res_files['bands_sc.1'], "0.91625.*0.94323")
        assert grep_exists(res_files['bands_sc.1'], "14.0343.*0.03976")
        assert grep_exists(res_files['bands_sc.2'], "18.1958.*0.62231")
        assert grep_exists(res_files['bands_sc.2'], "27.1348.*0.00942")
    else:
        assert grep_exists(res_files['bands_sc.1'], "-8.92164.*0.94323")
        assert grep_exists(res_files['bands_sc.1'], "6.028964.*0.039764")
        assert grep_exists(res_files['bands_sc.2'], "10.19046.*0.622318")
        assert grep_exists(res_files['bands_sc.2'], "19.12943.*0.009426")

