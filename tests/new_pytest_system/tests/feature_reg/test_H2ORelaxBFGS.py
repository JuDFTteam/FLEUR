
import pytest

@pytest.mark.serial
@pytest.mark.slow
@pytest.mark.relaxation
@pytest.mark.forces
def test_H2ORelaxBFGS(execute_fleur, grep_number, grep_exists):
    """H2O: structural relaxation with BFGS

    Simple test of Fleur with a single steps:
    1.Run Fleur 3 times and check the agreement of the relaxation with the reference.
    Uses: forces, BFGS
    #TODo: Maybe change this into a full regression test of the relax.xml file
    """
    test_file_folder = './inputfiles/H2ORelaxBFGS/'
    # Run fleur 3 times
    res_files = execute_fleur(test_file_folder)
    res_files = execute_fleur(test_file_folder)
    res_files = execute_fleur(test_file_folder)

    res_file_names = list(res_files.keys())
    should_files = ['out', 'relax.xml']
    for file1 in should_files:
        assert (file1 in res_file_names), f'{file1} missing'

    #grep for displacements
    assert grep_exists(res_files['relax.xml'], "0.267210")
    assert grep_exists(res_files['relax.xml'], "0.117072")
    assert grep_exists(res_files['relax.xml'], "0.138861")


    #grep for forces in different calls
    assert grep_exists(res_files['relax.xml'], "0.155704")
    assert grep_exists(res_files['relax.xml'], "0.005702")
    assert grep_exists(res_files['relax.xml'], "0.078131")

    assert grep_exists(res_files['relax.xml'], "0.060642")
    assert grep_exists(res_files['relax.xml'], "0.024061")
    assert grep_exists(res_files['relax.xml'], "0.031505")

    assert grep_exists(res_files['relax.xml'], "0.016925")
    assert grep_exists(res_files['relax.xml'], "0.032709")
    assert grep_exists(res_files['relax.xml'], "0.009706")

    tenergy = grep_number(res_files['out'], "total energy=", "=")
    assert abs(tenergy - -75.961947) <= 0.000002
