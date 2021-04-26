
import pytest

@pytest.mark.hybrid
def test_KClHybridPBE0(execute_fleur, check_value_outfile, fleur_binary):
    """Fleur KCl Hybrid PBE0

    Simple test of Fleur with XML input with one step:
    1.Generate a starting density and run 2 HF iterations and compare total energies and other quantities
    """
    # TODO: rewirte execute fleur in confests such, that it gets cmd param before the executable
    # In this case the judf_MPI environment variable should be ignored
    test_file_folder = './inputfiles/KClHybridPBE0/'
    cmd_params = ['-trace']

    # special for this hybrid test:
    fleur, parallel = fleur_binary
    if parallel:
        env = {'juDFT_MPI': "mpirun -n 3 --allow-run-as-root --mca btl vader,self"}
    else:
        env = {}

    res_files = execute_fleur(test_file_folder, cmdline_param=cmd_params, env=env)
    res_file_names = list(res_files.keys())
    should_files = ['out']
    for file1 in should_files:
        assert file1 in res_file_names
    
    assert check_value_outfile(res_files['out'], "HF total energy=", "htr", [-1063.8587731477, -1063.8383730939], 0.000001)
    # only check the last bandgap
    exp_bandgap = 28 * [None]
    exp_bandgap[27] = 0.2771
    assert check_value_outfile(res_files['out'], "bandgap                     :", "htr", exp_bandgap, 0.0001)

@pytest.mark.eigpara
@pytest.mark.hybrid
def test_KClHybridPBE0_eigpar(execute_fleur, check_value_outfile, fleur_binary):
    """Fleur KCl Hybrid PBE0

    Simple test of Fleur with XML input with one step:
    1.Generate a starting density and run 2 HF iterations and compare total energies and other quantities
    """
    test_file_folder = './inputfiles/KClHybridPBE0_eigpar/'
    cmd_params = ['-trace']

    # special for this hybrid test:
    fleur, parallel = fleur_binary
    if parallel:
        env = {'juDFT_MPI': "mpirun -n 6 --allow-run-as-root --mca btl vader,self"}
    else:
        env = {}

    res_files = execute_fleur(test_file_folder, cmdline_param=cmd_params, env=env)
    res_file_names = list(res_files.keys())
    should_files = ['out']
    for file1 in should_files:
        assert file1 in res_file_names
    
    assert check_value_outfile(res_files['out'], "HF total energy=", "htr", [-1063.8587731477, -1063.8383730939], 0.000001)
    # only check the last bandgap
    exp_bandgap = 28 * [None]
    exp_bandgap[27] = 0.2771
    assert check_value_outfile(res_files['out'], "bandgap                     :", "htr", exp_bandgap, 0.0001)