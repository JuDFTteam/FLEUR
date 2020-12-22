import pytest

@pytest.mark.hybrid
def test_GaAsHybridPBEO(execute_fleur, grep_number, grep_exists):
    """Fleur GaAs Hybrid PBE0

    Simple test of Fleur with XML input with one step:
    1.Generate a starting density and run 2 HF iterations and compare total energies and other quantities

    """
    test_file_folder = './inputfiles/GaAsHybridPBEO/'
    cmd_params = ['-trace']

    # special for this hybrid test:
    fleur, parallel = fleur_binary
    if parallel:
        env = {'judft_MPI': "mpirun -n 3 --allow-run-as-root --mca btl vader,self"}
    else:
        env = {}

    res_files = execute_fleur(test_file_folder, cmdline_param=cmd_params, env=env)
    res_file_names = list(res_files.keys())
    should_files = ['out']
    for file1 in should_files:
        assert file1 in res_file_names
    
    assert check_value_outfile(res_files['out'], "HF total energy=", "htr", [-4205.2193168, -4204.9174214647], 0.000001)
    # only check the last bandgap
    exp_bandgap = 38 * [None]
    exp_bandgap[37] = 0.09105
    check_value_outfile(res_files['out'], "bandgap                     :", "htr", exp_bandgap, 0.0001)

@pytest.mark.eigpar
@pytest.mark.hybrid
def test_GaAsHybridPBEO_eigpar(execute_fleur, grep_number, grep_exists):
    """Fleur GaAs Hybrid PBE0

    Simple test of Fleur with XML input with one step:
    1.Generate a starting density and run 2 HF iterations and compare total energies and other quantities
    """
    test_file_folder = './inputfiles/GaAsHybridPBEO_eigpar/'
    cmd_params = ['-trace']

    # special for this hybrid test:
    fleur, parallel = fleur_binary
    if parallel:
        env = {'judft_MPI': "mpirun -n 6 --allow-run-as-root --mca btl vader,self"}
    else:
        env = {}

    res_files = execute_fleur(test_file_folder, cmdline_param=cmd_params, env=env)
    res_file_names = list(res_files.keys())
    should_files = ['out']
    for file1 in should_files:
        assert file1 in res_file_names
    
    assert check_value_outfile(res_files['out'], "HF total energy=", "htr", [-4205.2193168, -4204.9174214647], 0.000001)
    # only check the last bandgap
    exp_bandgap = 38 * [None]
    exp_bandgap[37] = 0.09105
    check_value_outfile(res_files['out'], "bandgap                     :", "htr", exp_bandgap, 0.0001)


@pytest.mark.ldau
@pytest.mark.forces
def test_GaAsMulitUForceXML(execute_fleur, grep_number, grep_exists):
    """GaAs: displaced atoms with U parameters and forces

    Simple test of Fleur with three steps:
    1.Generate a starting density and run 1 iteration
    2.Remove Broyden files, last line of n_mmp_mat and run 12 iterations
    3.Calculate forces
    """
    
    
    # Stage 1
    # Stage 2
    # Stage 3

@pytest.mark.soc
@pytest.mark.wannier
def test_GaAsWannSOC(execute_fleur, grep_number, grep_exists):
    """GaAs: simple test for the Wannier code with inp.xml and SOC=T

    Simple test of Fleur with three steps:
    1.Generate a starting density and run 1 iteration
    2.Generate projections for Wannier functions
    3.Generate input for Wannier90 code
    """
    
    test_file_folder = './inputfiles/GaAsWannSOC/'
    
    # Stage 1
    res_files = execute_fleur(test_file_folder, only_copy=['inp', 'sym.out', 'enpara', 'kpts'])
    res_file_names = list(res_files.keys())
    should_files = ['out']
    for file1 in should_files:
        assert (file1 in res_file_names), file1
    
    # Stage 2
    res_files = execute_fleur(test_file_folder, only_copy=[['wann_inp2', 'wann_inp'], 'projgen_inp', ['kpts-2', 'kpts']])
    res_file_names = list(res_files.keys())
    should_files = ['proj']
    for file1 in should_files:
        assert (file1 in res_file_names), file1
    assert grep_exists(res_files['proj'],"           8           8")
    assert grep_exists(res_files['proj'],"  1 -3  1  0")
    assert grep_exists(res_files['proj'],"  2 -3  4  0")
    assert grep_exists(res_files['proj'],"     0.000000  0.000000  0.000000  0.000000 1.00")

    # Stage 3
    res_files = execute_fleur(test_file_folder, only_copy=[['wann_inp2', 'wann_inp']])
    res_file_names = list(res_files.keys())
    should_files = ['WF1.amn', 'WF1.mmn', 'WF1.eig', 'WF1.win', 'WF1.wout']
    for file1 in should_files:
        assert (file1 in res_file_names), file1
    
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