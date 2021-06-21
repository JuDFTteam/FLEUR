
import pytest

@pytest.mark.wannier
def test_CwannXML(execute_fleur, grep_number, grep_exists):
    """C: simple test for the Wannier code with inp.xml
    Simple test of Fleur with three steps:
    1.Generate a starting density and run 1 iteration
    2.Generate projections for Wannier functions
    3.Generate input for Wannier90 code
    """
    test_file_folder = './inputfiles/CwannXML/'
    
    # Stage 1
    res_files = execute_fleur(test_file_folder, only_copy=['inp.xml', 'sym.xml', 'kpts.xml'])
    res_file_names = list(res_files.keys())
    should_files = ['out']
    for file1 in should_files:
        assert (file1 in res_file_names), file1

    # Stage 2
    res_files = execute_fleur(test_file_folder, only_copy=[['wann_inp-1.xml', 'inp.xml'], 'projgen_inp'])
    res_file_names = list(res_files.keys())
    should_files = ['proj']
    for file1 in should_files:
        assert (file1 in res_file_names), file1
    assert grep_exists(res_files['proj'],"           8           8")
    assert grep_exists(res_files['proj'],"  1 -3  1  0")
    assert grep_exists(res_files['proj'],"  2 -3  4  0")
    assert grep_exists(res_files['proj'],"     0.000000  0.000000  0.000000  0.000000 1.00")

    # Stage 3
    res_files = execute_fleur(test_file_folder, only_copy=[['wann_inp-2.xml', 'inp.xml']])
    res_file_names = list(res_files.keys())
    should_files = ['WF1.amn', 'WF1.mmn', 'WF1.eig', 'WF1.win', 'WF1.wout']
    for file1 in should_files:
        assert (file1 in res_file_names), file1

    # These eigenvalues differ from the wann test without xml above
    assert grep_exists(res_files['WF1.eig'], "           1           1  -10.236815")
    assert grep_exists(res_files['WF1.eig'], "           7           1   16.700987")
    assert grep_exists(res_files['WF1.eig'], "           7           8   19.833238")
    assert grep_exists(res_files['WF1.eig'], "           5           6   16.235248")
    assert grep_exists(res_files['WF1.mmn0'], "    8    8      8       1.000000000000")
    assert grep_exists(res_files['WF1.mmn0'], "    1    1      1       1.000000000000")
