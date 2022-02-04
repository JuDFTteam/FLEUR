import pytest

@pytest.mark.hybrid
def test_FeHybridPBE0(execute_fleur, check_value_outfile, fleur_binary):
    """Fleur Fe Hybrid PBE0

    Simple test of Fleur with XML input with one step:
    1.Generate a starting density and run 2 HF iterations and compare total energies and other quantities
    """
    test_file_folder = './inputfiles/FeHybridPBE0/'
    cmd_params = ['-trace']

    res_files = execute_fleur(test_file_folder, cmdline_param=cmd_params, mpi_procs=3)
    res_file_names = list(res_files.keys())
    should_files = ['out']
    for file1 in should_files:
        assert file1 in res_file_names
    
    assert check_value_outfile(res_files['out'], "HF total energy=", "htr", [-1272.7289531184, -1272.7324870224], 0.000001)
    exp_mm1 = 27*[None]
    exp_mm1[14] = 3.40037
    exp_mm1[-1] = 3.42074
    check_value_outfile(res_files['out'], "--> mm       1", " ", exp_mm1, 0.00001)


@pytest.mark.hybrid
@pytest.mark.eigpara
def test_FeHybridPBE0_eigpar(execute_fleur, check_value_outfile, fleur_binary):
    """Fleur Fe Hybrid PBE0

    Simple test of Fleur with XML input with one step:
    1.Generate a starting density and run 2 HF iterations and compare total energies and other quantities

    # How do you know eigpar really happened and that it worked well?
    # Also if it is a serial binary it does not do eigpar, right?
    """
    test_file_folder = './inputfiles/FeHybridPBE0_eigpar/'
    cmd_params = ['-trace']

    res_files = execute_fleur(test_file_folder, cmdline_param=cmd_params, mpi_procs=6)
    res_file_names = list(res_files.keys())
    should_files = ['out']
    for file1 in should_files:
        assert file1 in res_file_names
    
    assert check_value_outfile(res_files['out'], "HF total energy=", "htr", [-1272.7289531184, -1272.7324870224], 0.000001)
    exp_mm1 = 27*[None]
    exp_mm1[14] = 3.40037
    exp_mm1[-1] = 3.42074
    check_value_outfile(res_files['out'], "--> mm       1", " ", exp_mm1, 0.00001)