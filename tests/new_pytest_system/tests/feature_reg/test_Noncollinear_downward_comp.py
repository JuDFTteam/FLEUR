
import pytest

@pytest.mark.magnetism
@pytest.mark.non_collinear
@pytest.mark.hdf
def test_Noncollinear_downward_compatible(execute_fleur, grep_number, grep_exists):
    """Tests downward compatibility of collinear, noncollinear and FFN calclations.

    Simple test of Fleur with four steps:
    1. Generate a starting density and run 1 iteration and check Energy of Fe in collinear steup. Run iteration with noco and check again.
    2. Run iteration with FFN and check again.
    3. Generate a starting density and run 1 iteration and check Energy of Fe in collinear setup. Run iteration with FFN and check again.
    4. Invoke LDA+U on Density from step three and check for distance.
    """
    test_file_folder = './inputfiles/Noncollinear_downward_compatible/'

    # Stage 1
    res_files = execute_fleur(test_file_folder, only_copy=['inp.xml', 'kpts.xml', 'sym.xml', 'JUDFT_WARN_ONLY'])
    res_files = execute_fleur(test_file_folder, only_copy=[['inp-2.xml', 'inp.xml']], rm_files=['mixing_history'])
    res_file_names = list(res_files.keys())
    should_files = ['out', 'cdn.hdf']
    for file1 in should_files:
        assert (file1 in res_file_names), f'{file1} missing'
    

    tenergy = grep_number(res_files['out'], "total energy=", "=")
    assert abs(tenergy - -1270.4886) <= 0.0002

    # Stage 2
    res_files = execute_fleur(test_file_folder, only_copy=[['inp-3.xml', 'inp.xml']], rm_files=['mixing_history'])
    res_file_names = list(res_files.keys())
    should_files = ['out']
    for file1 in should_files:
        assert (file1 in res_file_names), f'{file1} missing'
    

    tenergy = grep_number(res_files['out'], "total energy=", "=")
    assert abs(tenergy - -1270.4873) <= 0.0002

    # Stage 3
    res_files = execute_fleur(test_file_folder, only_copy=['inp.xml'], rm_files=['cdn.hdf', 'mixing_history'])
    res_file_names = list(res_files.keys())
    should_files = ['out', 'cdn.hdf']
    for file1 in should_files:
        assert (file1 in res_file_names), f'{file1} missing'
    
    tenergy = grep_number(res_files['out'], "total energy=", "=")
    assert abs(tenergy - -1270.491) <= 0.0002

    res_files = execute_fleur(test_file_folder, only_copy=[['inp-3.xml', 'inp.xml']], rm_files=['mixing_history'])
    res_file_names = list(res_files.keys())
    should_files = ['out']
    for file1 in should_files:
        assert (file1 in res_file_names), f'{file1} missing'

    tenergy = grep_number(res_files['out'], "total energy=", "=")
    assert abs(tenergy - -1270.4886) <= 0.0002

    # Stage 4

    res_files = execute_fleur(test_file_folder, only_copy=[['inp-4.xml', 'inp.xml'], ['sym-2.xml', 'sym.xml']], rm_files=['mixing_history'])
    res_file_names = list(res_files.keys())
    should_files = ['out']
    for file1 in should_files:
        assert (file1 in res_file_names), f'{file1} missing'

    dist = grep_number(res_files['out'], "distance of charge densities for it=    4", ":")

    assert abs(dist - 130.55) <= 0.01
