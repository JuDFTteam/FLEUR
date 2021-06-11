"""
Collection of tests running only inpgen
"""

import pytest


@pytest.mark.inpgen
@pytest.mark.bulk
@pytest.mark.fast
def test_inpgen_Si_simple_interface_files(execute_inpgen, clean_workdir):
    """Simple test of Inpgen:
    1.Generate input files where no namelists are given & compare with expected files
    2.Generate input files with explicit option & compare with expected files
    """
    test_file_folder = './inputfiles/inpgen/Si_no_para/'

    # Test if inpgen runs and required files are there
    cmd_params = ['-f', 'inp_simple']
    res_files = execute_inpgen(test_file_folder, cmd_params)
    should_files = ['inp.xml', 'struct.xsf', 'default.econfig', 'out', 'usage.json', 'sym.xml', 'kpts.xml']
    res_file_names = list(res_files.keys())
    
    for file1 in should_files:
        assert file1 in res_file_names

    clean_workdir()

    # Test if inpgen can run with --explicit and different includes are right in inp.xml
    cmd_params = ['-explicit', '-inc', 'all', '-f', 'inp_simple']
    res_files = execute_inpgen(test_file_folder, cmd_params)
    should_files = ['inp.xml', 'struct.xsf', 'default.econfig', 'out', 'usage.json']
    not_there = ['sym.xml', 'kpts.xml']
    res_file_names = list(res_files.keys())
    for file1 in should_files:
        assert file1 in res_file_names
    for file1 in not_there:
        assert file1 not in res_file_names


@pytest.mark.inpgen
@pytest.mark.bulk
@pytest.mark.fast
def test_inpgen_Si_inputfile_content(execute_inpgen, clean_workdir):
    """Simple test of Inpgen:
    1.Generate input files where no namelists are given & compare with expected files
    Be very careful when changing the inpgen input file, it is parsed by various other tools:
    
    Masci-tools (used by aiida-fleur): https://github.com/JuDFTteam/masci-tools/tree/develop/masci_tools/io/io_fleur_inpgen.py
    Phononpy: https://github.com/phonopy/phonopy/blob/develop/phonopy/interface/fleur.py
    # ASE
    """
    test_file_folder = './inputfiles/inpgen/Si_full_para/'

    # Test if inpgen runs and required files are there
    cmd_params = ['-f', 'inp_simple']
    res_files = execute_inpgen(test_file_folder, cmd_params)
    should_files = ['inp.xml', 'struct.xsf', 'default.econfig', 'out', 'usage.json', 'sym.xml', 'kpts.xml']
    res_file_names = list(res_files.keys())
    
    # if it would complain about a list the inp.xml would not be there
    for file1 in should_files:
        assert file1 in res_file_names

    # Now we test if the namelists are really set in the inp.xml
