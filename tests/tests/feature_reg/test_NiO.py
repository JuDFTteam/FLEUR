import pytest


@pytest.mark.serial
@pytest.mark.ldau
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

    assert grep_exists(res_files['out'], "it= 25  is completed")
    tenergy = grep_number(res_files['out'], "total energy=", "=")
    dist = grep_number(res_files['out'], "distance of charge densities for it=   25", ":")
    efermi = grep_number(res_files['out'], "first approx. to ef", ":")

    assert abs(tenergy - -3189.616) <= 0.001
    assert abs(dist - 0.002572) <= 0.0001
    assert abs(efermi - 0.240) <= 0.001