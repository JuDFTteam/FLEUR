
import pytest

@pytest.mark.relaxation
@pytest.mark.magnetism
@pytest.mark.hdf
def DEACTIVATED_test_RelaxMTFeature(execute_fleur, grep_number, grep_exists):
    """Tests relaxation feature of FFN in the MT

    Simple test of Fleur with one step:
    1.Generate a starting density and run 5 iterations and check direction of magnetization. 
    """
    test_file_folder = './inputfiles/RelaxMTFeature/'
    res_files = execute_fleur(test_file_folder)
    res_file_names = list(res_files.keys())
    should_files = ['out']
    for file1 in should_files:
        assert (file1 in res_file_names), f'{file1} missing'

    noco_beta = grep_number(res_files['out'], "nococonv%beta=", "nococonv%beta=")
    noco_alpha = grep_number(res_files['out'], "nococonv%alpha=", "nococonv%alpha=")

    assert abs(noco_beta - 1.54839) <= 0.0005
    assert abs(noco_alpha - 1.57079) <= 0.0005

