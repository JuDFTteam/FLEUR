import pytest

@pytest.mark.hybrid
def test_MnHybridNoinv(execute_fleur, check_value_outfile, fleur_binary):
    """Fleur Mn Hybrid Noinv

    Simple test of Fleur with XML input with one step:
    1.Generate a starting density and run 2 HF iterations and compare total energies and other quantities
    """
    test_file_folder = './inputfiles/MnHybridNoinv/'
    cmd_params = ['-trace']

    res_files = execute_fleur(test_file_folder, cmdline_param=cmd_params, mpi_procs=3)
    res_file_names = list(res_files.keys())
    should_files = ['out']
    for file1 in should_files:
        assert file1 in res_file_names
    
    assert check_value_outfile(res_files['out'], "HF total energy=", "htr", [-2317.1670995071, -2317.1757499769], 0.000001)
    exp_mm1 = 36*[None]
    exp_mm1[-1] = 2.4744
    assert check_value_outfile(res_files['out'], "--> mm       1", " ", exp_mm1, 0.0001)
 
    exp_mm2 = 36*[None]
    exp_mm2[-1] =-2.4744
    assert check_value_outfile(res_files['out'], "--> mm       2", " ", exp_mm2, 0.0001)
 
    # only check the last bandgap
    exp_bandgap = 36 * [None]
    exp_bandgap[-1] = 0.086591
    assert check_value_outfile(res_files['out'], "bandgap                     :", "htr", exp_bandgap, 0.0001)

@pytest.mark.hybrid
@pytest.mark.eigpara
def test_MnHybridNoinv_eigpar(execute_fleur, check_value_outfile, fleur_binary):
    """Fleur Mn Hybrid Noinv

    Simple test of Fleur with XML input with one step:
    1.Generate a starting density and run 2 HF iterations and compare total energies and other quantities
    """
    test_file_folder = './inputfiles/MnHybridNoinv_eigpar/'
    cmd_params = ['-trace']

    res_files = execute_fleur(test_file_folder, cmdline_param=cmd_params, mpi_procs=6)
    res_file_names = list(res_files.keys())
    should_files = ['out']
    for file1 in should_files:
        assert file1 in res_file_names
    
    assert check_value_outfile(res_files['out'], "HF total energy=", "htr", [-2317.1670995071, -2317.1757499769], 0.000001)
    exp_mm1 = 36*[None]
    exp_mm1[-1] = 2.4744
    assert check_value_outfile(res_files['out'], "--> mm       1", " ", exp_mm1, 0.0001)

    exp_mm2 = 36*[None]
    exp_mm2[-1] =-2.4744
    assert check_value_outfile(res_files['out'], "--> mm       2", " ", exp_mm2, 0.0001)

    # only check the last bandgap
    exp_bandgap = 36 * [None]
    exp_bandgap[-1] = 0.086591
    assert check_value_outfile(res_files['out'], "bandgap                     :", "htr", exp_bandgap, 0.0001)