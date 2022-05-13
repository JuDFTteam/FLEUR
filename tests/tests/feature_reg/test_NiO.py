import pytest


@pytest.mark.serial
@pytest.mark.ldau
@pytest.mark.slow
@pytest.mark.magnetism
def test_NiOldaUAMF(execute_fleur, grep_number, grep_exists):
    """NiO: test with LDA+U with Around mean field double counting and magnetism

    Simple test of Fleur with one steps:
    1.Generate a starting density and run 25 iterations with LDA+U with AMF double counting
    """
    test_file_folder = './inputfiles/NiOldaUAMF/'

    # Stage 1
    res_files = execute_fleur(test_file_folder, only_copy=['inp.xml'])
    res_file_names = list(res_files.keys())
    should_files = ['out']
    if 'cdn.hdf' not in res_file_names:
        should_files.append('cdn1')
    for file1 in should_files:
        assert (file1 in res_file_names), f'{file1} missing'

    assert grep_exists(res_files['out'], "it= 12  is completed")
    tenergy = grep_number(res_files['out'], "total energy=", "=")
    dist = grep_number(res_files['out'], "distance of charge densities for it=   12", ":")
    efermi = grep_number(res_files['out'], "first approx. to ef", ":")
    bandgap = grep_number(res_files['out'], "bandgap", ":")
    mm_atom1 = grep_number(res_files['out'], "mm       1", " 1 ")
    mm_atom3 = grep_number(res_files['out'], "mm       3", " 3 ")

    assert abs(tenergy - -3189.615) <= 0.001
    assert abs(dist - 0.33) <= 0.01
    assert abs(efermi - 0.246) <= 0.001
    assert abs(bandgap - 0.085) <= 0.001
    assert abs(mm_atom1 - 1.85) <= 0.01
    assert abs(mm_atom3 + 1.85) <= 0.01

@pytest.mark.serial
@pytest.mark.ldau
@pytest.mark.slow
@pytest.mark.magnetism
def test_NiOldaUFLL(execute_fleur, grep_number, grep_exists):
    """NiO: test with LDA+U with fully localized double counting and magnetism

    Simple test of Fleur with one steps:
    1.Generate a starting density and run 12 iterations with LDA+U with FLL double counting
    """
    test_file_folder = './inputfiles/NiOldaUFLL/'

    # Stage 1
    res_files = execute_fleur(test_file_folder, only_copy=['inp.xml'])
    res_file_names = list(res_files.keys())
    should_files = ['out']
    if 'cdn.hdf' not in res_file_names:
        should_files.append('cdn1')
    for file1 in should_files:
        assert (file1 in res_file_names), f'{file1} missing'

    assert grep_exists(res_files['out'], "it= 12  is completed")
    tenergy = grep_number(res_files['out'], "total energy=", "=")
    dist = grep_number(res_files['out'], "distance of charge densities for it=   12", ":")
    efermi = grep_number(res_files['out'], "first approx. to ef", ":")
    bandgap = grep_number(res_files['out'], "bandgap", ":")
    mm_atom1 = grep_number(res_files['out'], "mm       1", " 1 ")
    mm_atom3 = grep_number(res_files['out'], "mm       3", " 3 ")

    assert abs(tenergy - -3191.458) <= 0.001
    assert abs(dist - 0.58) <= 0.01
    assert abs(efermi - 0.207) <= 0.001
    assert abs(bandgap - 0.124) <= 0.001
    assert abs(mm_atom1 - 1.85) <= 0.01
    assert abs(mm_atom3 - -1.85) <= 0.01