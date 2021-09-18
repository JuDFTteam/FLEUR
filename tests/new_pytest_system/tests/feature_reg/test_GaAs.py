import pytest

@pytest.mark.hybrid
def test_GaAsHybridPBE0(execute_fleur, check_value_outfile, fleur_binary):
    """Fleur GaAs Hybrid PBE0

    Simple test of Fleur with XML input with one step:
    1.Generate a starting density and run 2 HF iterations and compare total energies and other quantities

    """
    test_file_folder = './inputfiles/GaAsHybridPBE0/'
    cmd_params = ['-trace']

    res_files = execute_fleur(test_file_folder, cmdline_param=cmd_params, mpi_procs=3)
    res_file_names = list(res_files.keys())
    should_files = ['out']
    for file1 in should_files:
        assert (file1 in res_file_names), f'{file1} missing'
    
    assert check_value_outfile(res_files['out'], "HF total energy=", "htr", [-4205.2193168, -4204.9174214647], 0.000001)
    # only check the last bandgap
    exp_bandgap = 38 * [None]
    exp_bandgap[37] = 0.09105
    check_value_outfile(res_files['out'], "bandgap                     :", "htr", exp_bandgap, 0.0001)

@pytest.mark.eigpara
@pytest.mark.hybrid
def test_GaAsHybridPBE0_eigpar(execute_fleur, check_value_outfile, fleur_binary):
    """Fleur GaAs Hybrid PBE0

    Simple test of Fleur with XML input with one step:
    1.Generate a starting density and run 2 HF iterations and compare total energies and other quantities
    """
    test_file_folder = './inputfiles/GaAsHybridPBE0_eigpar/'
    cmd_params = ['-trace']

    res_files = execute_fleur(test_file_folder, cmdline_param=cmd_params, mpi_procs=6)
    res_file_names = list(res_files.keys())
    should_files = ['out']
    for file1 in should_files:
        assert (file1 in res_file_names), f'{file1} missing'
    
    assert check_value_outfile(res_files['out'], "HF total energy=", "htr", [-4205.2193168, -4204.9174214647], 0.000001)
    # only check the last bandgap
    exp_bandgap = 38 * [None]
    exp_bandgap[37] = 0.09105
    check_value_outfile(res_files['out'], "bandgap                     :", "htr", exp_bandgap, 0.0001)

@pytest.mark.serial
@pytest.mark.ldau
@pytest.mark.forces
def test_GaAsMultiUForceXML(execute_fleur, grep_number, grep_exists):
    """GaAs: displaced atoms with U parameters and forces

    Simple test of Fleur with three steps:
    1.Generate a starting density and run 1 iteration
    2.Remove Broyden files, last line of n_mmp_mat and run 12 iterations
    3.Calculate forces
    """
    test_file_folder = './inputfiles/GaAsMultiUForceXML/'

    # Stage 1
    res_files = execute_fleur(test_file_folder, only_copy=['inp.xml'])
    res_file_names = list(res_files.keys())
    should_files = ['out']
    if not 'cdn.hdf' in res_file_names:
        should_files.append('cdn1')
    for file1 in should_files:
        assert (file1 in res_file_names), f'{file1} missing'

    tenergy = grep_number(res_files['out'], "total energy=", "=")
    dist = grep_number(res_files['out'], "distance of charge densitie", "1:")

    assert abs(tenergy - -4205.455) <= 0.001
    assert abs(dist - 6.656070) <= 0.0001

    # Stage 2
    res_files = execute_fleur(test_file_folder, only_copy=[['inp-2.xml', 'inp.xml'], 'n_mmp_mat'], rm_files=['mixing_history'])
    res_file_names = list(res_files.keys())
    should_files = ['out']
    for file1 in should_files:
        assert (file1 in res_file_names), f'{file1} missing'

    assert grep_exists(res_files['out'], "it= 12  is completed")
    efermi = grep_number(res_files['out'], "first approx. to ef", ":")
    tenergy = grep_number(res_files['out'], "total energy=", "=")

    assert abs(efermi - 0.1917) <= 0.001
    assert abs(tenergy - -4205.4319) <= 0.0001

    # Stage 3
    res_files = execute_fleur(test_file_folder, only_copy=[['inp-3.xml', 'inp.xml']], rm_files=['mixing_history'])
    res_file_names = list(res_files.keys())
    should_files = ['out', 'relax.xml']
    for file1 in should_files:
        assert (file1 in res_file_names), f'{file1} missing'

    assert grep_exists(res_files['relax.xml'], "4205.431")
    assert grep_exists(res_files['relax.xml'], "1.3806000000   -0.0094")


@pytest.mark.soc
@pytest.mark.wannier
@pytest.mark.skip('Test was also deactivated in the old test system. This has to be investigated')
def test_GaAsWannSOC(execute_fleur, grep_number, grep_exists):
    """GaAs: simple test for the Wannier code with inp.xml and SOC=T

    Simple test of Fleur with three steps:
    1.Generate a starting density and run 1 iteration
    2.Generate projections for Wannier functions
    3.Generate input for Wannier90 code
    """
    
    test_file_folder = './inputfiles/GaAsWannSOC/'
    
    # Stage 1
    res_files = execute_fleur(test_file_folder, only_copy=['inp.xml', 'sym.xml', 'kpts.xml'])
    res_file_names = list(res_files.keys())
    should_files = ['out']
    for file1 in should_files:
        assert (file1 in res_file_names), f'{file1} missing'
    
    # Stage 2
    res_files = execute_fleur(test_file_folder, only_copy=[['wann_inp-1.xml', 'inp.xml'], 'projgen_inp'])
    res_file_names = list(res_files.keys())
    should_files = ['proj']
    for file1 in should_files:
        assert (file1 in res_file_names), f'{file1} missing'
    assert grep_exists(res_files['proj'],"          16          36   t ")
    assert grep_exists(res_files['proj'],"  1 -3  1  0  1")
    assert grep_exists(res_files['proj'],"  2 -3  4  0 -1")
    assert grep_exists(res_files['proj'],"     0.000000  0.000000  0.000000  0.000000 1.00")

    # Stage 3
    res_files = execute_fleur(test_file_folder, only_copy=[['wann_inp-2.xml', 'inp.xml']])
    res_file_names = list(res_files.keys())
    should_files = ['WF1.amn', 'WF1.mmn', 'WF1.mmn0','WF1.eig', 'WF1.win', 'WF1.wout']
    for file1 in should_files:
        assert (file1 in res_file_names), f'{file1} missing'
    
    assert grep_exists(res_files['WF1.eig'], "           1           1   -7.869")
    assert grep_exists(res_files['WF1.eig'], "           8           1    4.840")
    assert grep_exists(res_files['WF1.eig'], "           8           8    3.769")
    assert grep_exists(res_files['WF1.eig'], "           5           6    2.167")
    assert grep_exists(res_files['WF1.mmn0'], "    8    8      8       1.000000000000")
    assert grep_exists(res_files['WF1.mmn0'], "    1    1      1       1.000000000000")
    assert grep_exists(res_files['WF1.wout'], "WF centre and spread    1  ( -2.041")
    assert grep_exists(res_files['WF1.wout'], "WF centre and spread    2  ( -2.041")
    assert grep_exists(res_files['WF1.wout'], "dis")
    assert grep_exists(res_files['WF1.wout'], "WF centre and spread    4  ( -0.642")
