import pytest


@pytest.mark.inpgen
@pytest.mark.bulk
@pytest.mark.fast
def test_inpgen_Si_simple_interface(execute_inpgen, clean_workdir):
    """Simple test of Inpgen with one steps:
    1.Generate input files where most namelists are given & compare with expected files
    """
    test_file_folder = './inputfiles/inpgen/Si_full_para/files/'

    
    cmd_params = ['-f', 'inp_simple']
    res_files = execute_inpgen(cmd_params, test_file_folder)
    should_files = ['inp.xml', 'struct.xsf', 'default.econfig', 'out', 'usage.json', 'sym.xml', 'kpts.xml']
    res_file_names = list(res_files.keys())
    
    # Test if all files are there
    for file1 in should_files:
        assert file1 in res_file_names

    clean_workdir()

    # do some more testing
