import pytest

@pytest.mark.serial
@pytest.mark.libxc
def test_AlLibxcPbe(execute_fleur, stage_workdir, grep_number, grep_exists):
    """Fleur Al Bulk libxc PBE

    Simple test of Fleur with XML input in two steps:
    1.Generate a starting density and run 1 iteration for the inbuild PBE-GGA functional.
    2.Run a second iteration with inbuild and libxc PBE respectively to compare the results.
    """
    test_file_folder = './inputfiles/Al_libxc_PBE/'

    # Stage 1
    res_files = execute_fleur(test_file_folder)
    #res_files = execute_fleur(test_file_folder, only_copy=[['inp.xml', 'inp.xml']])
    res_file_names = list(res_files.keys())
    should_files = ['out']
    for file1 in should_files:
        assert (file1 in res_file_names), f'{file1} missing'

    assert grep_exists(res_files['out'], "it=  1  is completed")
    efermi = grep_number(res_files['out'], "new fermi energy", ":")
    tenergy = grep_number(res_files['out'], "    total energy=", "=")
    dist = grep_number(res_files['out'], "distance of charge densities for spin  1                 it=    1", ":")

    assert abs(efermi - 0.2856) <= 0.001
    assert abs(tenergy - -242.820) <= 0.001
    assert abs(dist - 3.869) <= 0.001

    stage_workdir(foldername='inb')
    stage_workdir(foldername='lib')

    # Stage 2.1 [inbuilt, 2nd iteration]
    res_files = execute_fleur(test_file_folder, sub_dir='inb')
    res_file_names = list(res_files.keys())
    should_files = ['out']
    for file1 in should_files:
        assert (file1 in res_file_names), f'{file1} missing'

    assert grep_exists(res_files['out'], "it=  1  is completed")
    efermi = grep_number(res_files['out'], "new fermi energy", ":")
    tenergy = grep_number(res_files['out'], "    total energy=", "=")
    dist = grep_number(res_files['out'], "distance of charge densities for spin  1                 it=    2", ":")

    assert abs(efermi - 0.2856) <= 0.001
    assert abs(tenergy - -242.820) <= 0.001
    assert abs(dist - 3.660) <= 0.001

    # Stage 2.2 [libxc, 2nd iteration]
    res_files = execute_fleur(test_file_folder, only_copy=[['inp_libxc.xml', 'inp.xml']], sub_dir='lib')
    res_file_names = list(res_files.keys())
    should_files = ['out']
    for file1 in should_files:
        assert (file1 in res_file_names), f'{file1} missing'

    assert grep_exists(res_files['out'], "it=  1  is completed")
    efermi = grep_number(res_files['out'], "new fermi energy", ":")
    tenergy = grep_number(res_files['out'], "    total energy=", "=")
    dist = grep_number(res_files['out'], "distance of charge densities for spin  1                 it=    2", ":")

    assert abs(efermi - 0.2856) <= 0.001
    assert abs(tenergy - -242.820) <= 0.001
    assert abs(dist - 3.660) <= 0.001
