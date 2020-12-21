
import pytest

@pytest.mark.wannier
def test_Cwann(execute_fleur, grep_number, grep_exists):
    """C: simple test for the Wannier code
    Simple test of Fleur with three steps:
    1.Generate a starting density and run 1 iteration
    2.Generate projections for Wannier functions
    3.Generate input for Wannier90 code
    """
    test_file_folder = './inputfiles/Cwann/'
    cmd_params = []
    
    # Stage 1
    res_files = execute_fleur(cmd_params, test_file_folder, only_copy=['inp', 'sym.out', 'enpara', 'kpts'])
    res_file_names = list(res_files.keys())
    should_files = ['out']
    for file1 in should_files:
        assert file1 in res_file_names
    
    # Stage 2
    res_files = execute_fleur(cmd_params, test_file_folder, only_copy=[['wann_inp2', 'wann_inp'], 'projgen_inp', ['kpts-2', 'kpts']])
    res_file_names = list(res_files.keys())
    should_files = ['proj']
    for file1 in should_files:
        assert file1 in res_file_names
    assert grep_exists(res_files['proj'],"           8           8")
    assert grep_exists(res_files['proj'],"  1 -3  1  0")
    assert grep_exists(res_files['proj'],"  2 -3  4  0")

    # Stage 3
    res_files = execute_fleur(cmd_params, test_file_folder, only_copy=[['wann_inp2', 'wann_inp']])
    res_file_names = list(res_files.keys())
    should_files = ['WF1.amn', 'WF1.mmn', 'WF1.eig', 'WF1.win', 'WF1.wout']
    for file1 in should_files:
        assert file1 in res_file_names
    
    assert grep_exists(res_files['WF1.eig'], "           1           1   -9.96004957")
    assert grep_exists(res_files['WF1.eig'], "           8           1   24.36259922")
    assert grep_exists(res_files['WF1.eig'], "           8           8   27.60170066")
    assert grep_exists(res_files['WF1.eig'], "           5           6   16.51363508")
    #Note:WF1.amn and WF1.mmn seem to differ strongly from the reference result.
    #But this does not seem to be relevant as invoking wannier90.x WF1 yields high precision results.
    assert grep_exists(res_files['WF1.amn'], "    8    8    8")
    assert grep_exists(res_files['WF1.amn'], "    1    1    1      -0.21709150")
    assert grep_exists(res_files['WF1.amn'], "    8    2    1       0.000000000000    0.29174760")
    assert grep_exists(res_files['WF1.amn'], "    2    6    2       0.000000000000   -0.01714566")
    assert grep_exists(res_files['WF1.amn'], "    3    5    2      -0.18082932")
    assert grep_exists(res_files['WF1.amn'], "    8    8    8       0.000000000000   -0.11210751")
    assert grep_exists(res_files['WF1.amn'], "    6    7    8      -0.33695808")
    assert grep_exists(res_files['WF1.mmn'], "    8    8    8")
    assert grep_exists(res_files['WF1.mmn'], "    1    2      0   0   0")
    assert grep_exists(res_files['WF1.mmn'], "   -0.757938912603")
    assert grep_exists(res_files['WF1.mmn'], "    0.257799685127")
    assert grep_exists(res_files['WF1.mmn'], "    8    7      0   0   1")


@pytest.mark.wannier
@pytest.mark.skip('test not implemented')
def test_CwannXML(execute_fleur, grep_number, grep_exists):
    """
    """
    assert False