import pytest

@pytest.mark.serial
@pytest.mark.magnetism
@pytest.mark.dos
@pytest.mark.jdos
def test_SmAtomjDOS(execute_fleur, grep_number, grep_exists):
    """Fleur Sm Atom jDOS

    Test the decomposition into total angular momentum states in three stages:
    1. Run 15 iterations in a non-spin-polarized system with SOC and test for the fermi energy &total energy
    2. Generate a spin-polarized density with no magnetic moment using swsp (no checks here except for crash)
    3. Generate the jDOS and check the height of the peak around the fermi level
    """
    test_file_folder = './inputfiles/SmAtomjDOS/'

    # Stage 1
    res_files = execute_fleur(test_file_folder, only_copy=[['inp1.xml', 'inp.xml']])
    res_file_names = list(res_files.keys())
    should_files = ['out']
    for file1 in should_files:
        assert (file1 in res_file_names), f'{file1} missing'

    assert grep_exists(res_files['out'], "it= 15  is completed")
    efermi = grep_number(res_files['out'], "new fermi energy", ":")
    tenergy = grep_number(res_files['out'], "    total energy=", "=")
    dist = grep_number(res_files['out'], "distance of charge densities for spin  1                 it=   15", ":")

    assert abs(efermi - -0.0957) <= 0.005
    assert abs(tenergy - -10434.5474) <= 0.005
    assert abs(dist - 0.0000) <= 0.001


    # Stage 2
    res_files = execute_fleur(test_file_folder, only_copy=[['inp2.xml', 'inp.xml']])
    res_file_names = list(res_files.keys())
    should_files = ['out']
    for file1 in should_files:
        assert (file1 in res_file_names), f'{file1} missing'
    
    # Stage 3
    res_files = execute_fleur(test_file_folder, only_copy=[['inp3.xml', 'inp.xml']])
    res_file_names = list(res_files.keys())
    should_files = ['out']
    for file1 in should_files:
        assert (file1 in res_file_names), f'{file1} missing'

    assert ('banddos.hdf' in res_file_names) or ('jDOS.1' in res_file_names)

    assert grep_exists(res_files['out'], "0.2073")
    assert grep_exists(res_files['out'], "5.835")
    assert grep_exists(res_files['out'], "0.020")

    if 'jDOS.1' in res_file_names:
        assert grep_exists(res_files['jDOS.1'], "0.850633")
        assert grep_exists(res_files['jDOS.1'], "0.300025")
        assert grep_exists(res_files['jDOS.1'], "0.980324")
        assert grep_exists(res_files['jDOS.1'], "0.263514")
