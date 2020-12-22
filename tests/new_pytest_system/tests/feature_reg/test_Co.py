

import pytest

@pytest.mark.fleur
@pytest.mark.bulk
@pytest.mark.hybrid
def test_CoHybridPBE0(execute_fleur, grep_exists, grep_number):
    """Fleur Co Hybrid PBE0
    Simple test of Fleur with XML input with one step:
    1.Generate a starting density and run 1 HF iterations and compare total energies and other quantities
    """
    test_file_folder = './inputfiles/CoHybridPBE0/'
    cmd_params = []
    res_files = execute_fleur(test_file_folder)
    res_file_names = list(res_files.keys())
    should_files = ['out']
    for file1 in should_files:
        assert file1 in res_file_names
    
    tenergy = grep_number(res_files['out'], "HF total energy=", "= ")
    mm = grep_number(res_files['out'], "mm       1", "mm       1")

    assert abs(tenergy - -2786.7235930101) <= 0.000001
    assert abs(mm - 1.64918) <= 0.001


@pytest.mark.fleur
@pytest.mark.bulk
@pytest.mark.xml
@pytest.mark.band
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
        assert grep_exists(res_files['MCD.1'], "0.36559906")


@pytest.mark.fleur
@pytest.mark.bulk
@pytest.mark.band
def test_CoUnfold(execute_fleur, grep_exists, grep_number):
    """Co band unfolding test
    Simple test of Fleur with one step:
    1.Generate a starting density and run 1 iteration with band unfolding and compare several quantities
    """
    test_file_folder = './inputfiles/CoUnfold/'

    res_files = execute_fleur(test_file_folder, only_copy=['inp.xml', 'sym.xml'])
    res_file_names = list(res_files.keys())
    should_files = ['out', 'bands_sc.1', 'bands_sc.2']
    for file1 in should_files:
        assert file1 in res_file_names

    assert grep_exists(res_files['out'], "it=  1  is completed")
    assert grep_exists(res_files['bands_sc.1'], "0.91625.*0.94323")
    assert grep_exists(res_files['bands_sc.1'], "14.16776.*0.034036")
    assert grep_exists(res_files['bands_sc.1'], "18.195862.*0.622318")
    assert grep_exists(res_files['bands_sc.1'], "27.134829.*0.009426")
