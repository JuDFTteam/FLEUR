import pytest



@pytest.mark.bulk
@pytest.mark.greensfunction
@pytest.mark.magnetism
@pytest.mark.serial
def test_GreensFunction_Sphavg(execute_fleur, grep_number, grep_exists):
    """Fleur Fe bcc Green's function
    Simple test of the green's function calculation in FLEUR with one step:
    1. Generate starting density, run 1 Iteration and calculate Green's function
       for d-orbitals. Ensure that the occupations from the Green's function are
       close to the MT-charges obtained
    """
    test_file_folder = './inputfiles/Fe_bcc_GreensFunction/'

    res_files = execute_fleur(test_file_folder)
    should_files = ['out']
    res_file_names = list(res_files.keys())
    if 'cdn.hdf' in res_file_names:
        should_files.append('greensf.hdf')
    for file1 in should_files:
        assert file1 in res_file_names

    assert grep_exists(res_files['out'], "it=  1  is completed")
    spinup_trace = grep_number(res_files['out'], "Spin-Up trace:", ": ")
    spindown_trace = grep_number(res_files['out'], "Spin-Down trace:", ": ")

    assert abs(spinup_trace - 4.0508) <= 0.0005
    assert abs(spindown_trace - 1.8636) <= 0.0005


@pytest.mark.serial
@pytest.mark.film
@pytest.mark.greensfunction
def test_GreensFunction_SphavgFilm(execute_fleur, grep_number, grep_exists):
    """Fleur Fe Monolayer Green's function
    Simple test of the green's function calculation in FLEUR for films with one step:
    1. Generate starting density, run 1 Iteration and calculate Green's function
    for d-orbitals. Ensure that the occupations from the Green's function are
    close to the MT-charges obtained
    """
    test_file_folder = './inputfiles/Fe_1l_GreensFunction/'

    res_files = execute_fleur(test_file_folder)
    should_files = ['out']
    res_file_names = list(res_files.keys())
    if 'cdn.hdf' in res_file_names:
        should_files.append('greensf.hdf')
    for file1 in should_files:
        assert file1 in res_file_names

    assert grep_exists(res_files['out'], "it=  1  is completed")
    spinup_trace = grep_number(res_files['out'], "Spin-Up trace:", ":")
    spindown_trace = grep_number(res_files['out'], "Spin-Down trace:", ":")

    assert abs(spinup_trace - 4.8586) <= 0.0005
    assert abs(spindown_trace - 2.6652) <= 0.0005


@pytest.mark.serial
@pytest.mark.greensfunction
def test_GreensFunction_MultiContour(execute_fleur, grep_number, grep_exists):
    """Greens Function MultiContour

    Simple test of the green's function calculation in FLEUR with one step and multiple defined energy contours:
    1. Generate starting density, run 1 Iteration and calculate Green's function
       for d-orbitals. Ensure that the occupations are like expected for all contours
    """

    test_file_folder = './inputfiles/GreensFunction_MultiContour/'
    res_files = execute_fleur(test_file_folder)
    res_file_names = list(res_files.keys())
    should_files = ['out']
    if 'cdn.hdf' in res_file_names:
        should_files.append('greensf.hdf')
    for file1 in should_files:
        assert (file1 in res_file_names), f'{file1} missing'

    assert grep_exists(res_files['out'], "it=  1  is completed")
    # Test for the occupations (may need adjustments)
    # Rectangle
    cont_reg_sup = grep_number(res_files['out'], r"Contour\(RectangleExample\)    Spin-Up trace:", ":")
    cont_reg_sdown = grep_number(res_files['out'], r"Contour\(RectangleExample\)    Spin-Down trace:", ":")
    assert abs(cont_reg_sup - 3.9404) <= 0.0005
    assert abs(cont_reg_sdown - 1.8541) <= 0.0005

    # Semicircle
    cont_sem_sup = grep_number(res_files['out'], r"Contour\(default\)    Spin-Up trace:", ":")
    cont_sem_sdown = grep_number(res_files['out'], r"Contour\(default\)    Spin-Down trace:", ":")
    assert abs(cont_sem_sup - 4.0508) <= 0.0005
    assert abs(cont_sem_sdown - 1.8636) <= 0.0005

    # DOS (not weighted with fermi function)
    cont_dos_sup = grep_number(res_files['out'], r"Contour\(DOSExample\)    Spin-Up trace:", ":")
    cont_dos_sdown = grep_number(res_files['out'], r"Contour\(DOSExample\)    Spin-Down trace:", ":")
    assert abs(cont_dos_sup - 4.9905) <= 0.0005
    assert abs(cont_dos_sdown - 4.9905) <= 0.0005

@pytest.mark.serial
@pytest.mark.greensfunction
def test_GreensFunctionRadial(execute_fleur, grep_number, grep_exists):
    """Fleur Fe bcc Green's function Radial

    Simple test of the green's function calculation with radial dependence in FLEUR with one step:
    1. Generate starting density, run 1 Iteration and calculate Green's function
       for d-orbitals. Ensure that the occupations from the Green's function are
       close to the MT-charges obtained
    """
    test_file_folder = './inputfiles/GreensFunctionRadial/'
    res_files = execute_fleur(test_file_folder, only_copy=['inp.xml', 'JUDFT_WARN_ONLY'])
    res_file_names = list(res_files.keys())
    should_files = ['out']
    if 'cdn.hdf' in res_file_names:
        should_files.append('greensf.hdf')
    for file1 in should_files:
        assert (file1 in res_file_names), f'{file1} missing'

    assert grep_exists(res_files['out'], "it=  1  is completed")
    cont_reg_sup2 = grep_number(res_files['out'], r"l--> 2 Contour\(default\)    Spin-Up trace:", ":")
    cont_reg_sdown2 = grep_number(res_files['out'], r"l--> 2 Contour\(default\)    Spin-Down trace:", ":")
    cont_reg_sup1 = grep_number(res_files['out'], r"l--> 1 Contour\(default\)    Spin-Up trace:", ":")
    cont_reg_sdown1 = grep_number(res_files['out'], r"l--> 1 Contour\(default\)    Spin-Down trace:", ":")

    assert abs(cont_reg_sup2 - 4.0540) <= 0.0005
    assert abs(cont_reg_sdown2 - 1.8645) <= 0.0005
    assert abs(cont_reg_sup1 - 3.1967) <= 0.0005
    assert abs(cont_reg_sdown1 - 3.2171) <= 0.0005

@pytest.mark.serial
@pytest.mark.greensfunction
@pytest.mark.lo
def test_GreensFunctionRadial_LO(execute_fleur, grep_number, grep_exists):
    """Fleur Fe bcc Green's function Radial with local orbitals

    Simple test of the green's function calculation with radial dependence and local orbitals in FLEUR with two steps:
    1. Generate starting density, run 1 Iteration and calculate Green's function
       for d and p-orbitals (One SCLO on the p-orbitals). Ensure that the occupations from the Green's function are
       close to the MT-charges obtained
    1. Generate starting density, run 1 Iteration and calculate Green's function
       for d and p-orbitals (One SCLO on the p-orbitals,Each one HELO on d and p-orbitals). Ensure that the occupations from the Green's function are
       close to the MT-charges obtained
    """
    test_file_folder = './inputfiles/GreensFunctionRadial_LO/'

    # Stage 1
    res_files = execute_fleur(test_file_folder, only_copy=[['inp-1.xml', 'inp.xml'],
'JUDFT_WARN_ONLY'])
    res_file_names = list(res_files.keys())
    should_files = ['out']
    if 'cdn.hdf' in res_file_names:
        should_files.append('greensf.hdf')
    for file1 in should_files:
        assert (file1 in res_file_names), f'{file1} missing'

    assert grep_exists(res_files['out'], "it=  1  is completed")
    cont_reg_sup2 = grep_number(res_files['out'], r"2 Contour\(default\)    Spin-Up trace", ":")
    cont_reg_sdown2 = grep_number(res_files['out'], r"2 Contour\(default\)    Spin-Down trace", ":")
    cont_reg_sup1 = grep_number(res_files['out'], r"1 Contour\(default\)    Spin-Up trace", ":")
    cont_reg_sdown1 = grep_number(res_files['out'], r"1 Contour\(default\)    Spin-Down trace", ":")

    assert abs(cont_reg_sup2 - 4.0540) <= 0.0005
    assert abs(cont_reg_sdown2 - 1.8645) <= 0.0005
    assert abs(cont_reg_sup1 - 3.1967) <= 0.0005
    assert abs(cont_reg_sdown1 - 3.2171) <= 0.0005

    # Stage 2

    res_files = execute_fleur(test_file_folder, only_copy=[['inp-2.xml', 'inp.xml'], 'JUDFT_WARN_ONLY'])
    res_file_names = list(res_files.keys())
    should_files = ['out']
    if 'cdn.hdf' in res_file_names:
        should_files.append('greensf.hdf')
    for file1 in should_files:
        assert (file1 in res_file_names), f'{file1} missing'

    # now test output (Occupations are different because of different LO setup and it is the second iteration but consistent with mt charges)

    assert grep_exists(res_files['out'], "it=  1  is completed")
    cont_reg_sup2 = grep_number(res_files['out'], r"l--> 2 Contour\(default\)    Spin-Up trace:", ":")
    cont_reg_sdown2 = grep_number(res_files['out'], r"l--> 2 Contour\(default\)    Spin-Down trace:", ":")
    cont_reg_sup1 = grep_number(res_files['out'], r"l--> 1 Contour\(default\)    Spin-Up trace:", ":")
    cont_reg_sdown1 = grep_number(res_files['out'], r"l--> 1 Contour\(default\)    Spin-Down trace:", ":")

    assert abs(cont_reg_sup2 - 4.0815) <= 0.0005
    assert abs(cont_reg_sdown2 - 1.8776) <= 0.0005
    assert abs(cont_reg_sup1 - 3.1921) <= 0.0005
    assert abs(cont_reg_sdown1 - 3.2135) <= 0.0005


@pytest.mark.bulk
@pytest.mark.greensfunction
@pytest.mark.magnetism
@pytest.mark.serial
@pytest.mark.ldau
@pytest.mark.soc
def test_GreensFunction_HoAtom_SQA_theta(execute_fleur, grep_number, grep_exists):
    """Fleur Ho atom Green's function
    Simple test of the green's function calculation in FLEUR with one step:
    1. Generate starting density, run 1 Iteration and calculate Green's function
       for d-orbitals.  Make sure that occupation matrix matches the right result
       for theta=pi/2, phi=0
    """
    test_file_folder = './inputfiles/GreensFunction_HoAtom_SQA_theta/'

    res_files = execute_fleur(test_file_folder)
    should_files = ['out']
    res_file_names = list(res_files.keys())
    if 'cdn.hdf' in res_file_names:
        should_files.append('greensf.hdf')
    for file1 in should_files:
        assert file1 in res_file_names

    assert grep_exists(res_files['out'], "it=  2  is completed")

    assert grep_exists(res_files['out'], '0.9929  0.0000')
    assert grep_exists(res_files['out'], '0.9936  0.0000')
    assert grep_exists(res_files['out'], '0.1116  0.0000')
    assert grep_exists(res_files['out'], '0.3110  0.0000')
    assert grep_exists(res_files['out'], '0.1893  0.0000')
    assert grep_exists(res_files['out'], '0.2687  0.0000')
    spinup_trace = grep_number(res_files['out'], "Spin-Up trace:", ": ")
    spindown_trace = grep_number(res_files['out'], "Spin-Down trace:", ": ")

    assert abs(spinup_trace - 6.9528) <= 0.0005
    assert abs(spindown_trace - 1.9920) <= 0.0005


@pytest.mark.bulk
@pytest.mark.greensfunction
@pytest.mark.magnetism
@pytest.mark.serial
@pytest.mark.ldau
@pytest.mark.soc
def test_GreensFunction_HoAtom_SQA_phi(execute_fleur, grep_number, grep_exists):
    """Fleur Ho atom Green's function
    Simple test of the green's function calculation in FLEUR with one step:
    1. Generate starting density, run 1 Iteration and calculate Green's function
       for d-orbitals. Make sure that occupation matrix matches the right result
       for theta=phi=pi/2
    """
    test_file_folder = './inputfiles/GreensFunction_HoAtom_SQA_phi/'

    res_files = execute_fleur(test_file_folder)
    should_files = ['out']
    res_file_names = list(res_files.keys())
    if 'cdn.hdf' in res_file_names:
        should_files.append('greensf.hdf')
    for file1 in should_files:
        assert file1 in res_file_names

    assert grep_exists(res_files['out'], "it=  2  is completed")

    assert grep_exists(res_files['out'], '0.9929  0.0000')
    assert grep_exists(res_files['out'], '0.9936  0.0000')
    assert grep_exists(res_files['out'], '0.1116  0.0000')
    assert grep_exists(res_files['out'], '0.3110  0.0000')
    assert grep_exists(res_files['out'], '0.0000  0.1893')
    assert grep_exists(res_files['out'], '0.0000 -0.1893')
    assert grep_exists(res_files['out'], '0.0000  0.2687')
    assert grep_exists(res_files['out'], '0.0000 -0.2687')
    spinup_trace = grep_number(res_files['out'], "Spin-Up trace:", ": ")
    spindown_trace = grep_number(res_files['out'], "Spin-Down trace:", ": ")

    assert abs(spinup_trace - 6.9528) <= 0.0005
    assert abs(spindown_trace - 1.9920) <= 0.0005


@pytest.mark.bulk
@pytest.mark.greensfunction
@pytest.mark.magnetism
@pytest.mark.serial
@pytest.mark.ldau
@pytest.mark.soc
@pytest.mark.non_collinear
def test_GreensFunction_rotated_SQA_noco(execute_fleur, grep_number, grep_exists):
    """Fleur Ho atom Green's function
    Simple test of the green's function calculation in FLEUR with one step:
    1. Generate starting density, run 1 Iteration and calculate Green's function
       for d-orbitals. Make sure that occupation matrix matches the right result
       for theta=phi=pi/2 (same as second variation results)
    """
    test_file_folder = './inputfiles/GreensFunction_rotated_SQA_noco/'

    res_files = execute_fleur(test_file_folder)
    should_files = ['out']
    res_file_names = list(res_files.keys())
    if 'cdn.hdf' in res_file_names:
        should_files.append('greensf.hdf')
    for file1 in should_files:
        assert file1 in res_file_names

    assert grep_exists(res_files['out'], "it=  2  is completed")

    assert grep_exists(res_files['out'], '0.9929  0.0000')
    assert grep_exists(res_files['out'], '0.9936  0.0000')
    assert grep_exists(res_files['out'], '0.1116  0.0000')
    assert grep_exists(res_files['out'], '0.3110  0.0000')
    assert grep_exists(res_files['out'], '0.0000  0.1893')
    assert grep_exists(res_files['out'], '0.0000 -0.1893')
    assert grep_exists(res_files['out'], '0.0000  0.2687')
    assert grep_exists(res_files['out'], '0.0000 -0.2687')
    spinup_trace = grep_number(res_files['out'], "Spin-Up trace:", ": ")
    spindown_trace = grep_number(res_files['out'], "Spin-Down trace:", ": ")

    assert abs(spinup_trace - 6.9528) <= 0.0005
    assert abs(spindown_trace - 1.9920) <= 0.0005


@pytest.mark.serial
@pytest.mark.greensfunction
@pytest.mark.non_collinear
@pytest.mark.bulk
def test_GreensFunction_mperp_xdir(execute_fleur, grep_number, grep_exists):
    """Fleur Fe bcc Green's function Radial Noco spin offdiagonal

    Simple test of the green's function calculation with radial dependence in FLEUR with one step:
    1. Generate starting density, run 1 Iteration and calculate Green's function
       for d-orbitals. Ensure that the occupations from the Green's function are
       close to the MT-charges obtained (also spin offdiagonal components) for the
       second atom rotated to alpha=0 beta=pi/2
    """
    test_file_folder = './inputfiles/GreensFunction_mperp_xdir/'
    res_files = execute_fleur(test_file_folder, only_copy=['inp.xml', 'JUDFT_WARN_ONLY'])
    res_file_names = list(res_files.keys())
    should_files = ['out']
    if 'cdn.hdf' in res_file_names:
        should_files.append('greensf.hdf')
    for file1 in should_files:
        assert (file1 in res_file_names), f'{file1} missing'

    assert grep_exists(res_files['out'], "it=  1  is completed")
    spinup_trace_atom1 = grep_number(res_files['out'], r"l--> 2 Contour\(default\)    Spin-Up trace:", ":", res_index=-2)
    spindn_trace_atom1 = grep_number(res_files['out'], r"l--> 2 Contour\(default\)    Spin-Down trace:", ":", res_index=-2)
    spinoffdx_trace_atom1 = grep_number(res_files['out'], r"l--> 2 Contour\(default\)    Spin-Offd trace \(x\):", ":", res_index=-2)
    spinoffdy_trace_atom1 = grep_number(res_files['out'], r"l--> 2 Contour\(default\)    Spin-Offd trace \(y\):", ":", res_index=-2)


    spinup_trace_atom2 = grep_number(res_files['out'], r"l--> 2 Contour\(default\)    Spin-Up trace:", ":")
    spindn_trace_atom2 = grep_number(res_files['out'], r"l--> 2 Contour\(default\)    Spin-Down trace:", ":")
    spinoffdx_trace_atom2 = grep_number(res_files['out'], r"l--> 2 Contour\(default\)    Spin-Offd trace \(x\):", ":")
    spinoffdy_trace_atom2 = grep_number(res_files['out'], r"l--> 2 Contour\(default\)    Spin-Offd trace \(y\):", ":")

    assert abs(spinup_trace_atom1 - 4.5242) <= 0.0005
    assert abs(spindn_trace_atom1 - 2.3679) <= 0.0005
    assert abs(spinoffdx_trace_atom1 - 0.5994) <= 0.0005
    assert abs(spinoffdy_trace_atom1 - 0.0000) <= 0.0005

    assert abs(spinup_trace_atom2 - 4.5214) <= 0.0005
    assert abs(spindn_trace_atom2 - 2.38) <= 0.0005
    assert abs(spinoffdx_trace_atom2 + 0.6027) <= 0.0005
    assert abs(spinoffdy_trace_atom2 - 0.0000) <= 0.0005


@pytest.mark.serial
@pytest.mark.greensfunction
@pytest.mark.non_collinear
@pytest.mark.bulk
def test_GreensFunction_mperp_ydir(execute_fleur, grep_number, grep_exists):
    """Fleur Fe bcc Green's function Radial Noco spin offdiagonal

    Simple test of the green's function calculation with radial dependence in FLEUR with one step:
    1. Generate starting density, run 1 Iteration and calculate Green's function
       for d-orbitals. Ensure that the occupations from the Green's function are
       close to the MT-charges obtained (also spin offdiagonal components) for the
       second atom rotated to alpha=beta=pi/2
    """
    test_file_folder = './inputfiles/GreensFunction_mperp_ydir/'
    res_files = execute_fleur(test_file_folder, only_copy=['inp.xml', 'JUDFT_WARN_ONLY'])
    res_file_names = list(res_files.keys())
    should_files = ['out']
    if 'cdn.hdf' in res_file_names:
        should_files.append('greensf.hdf')
    for file1 in should_files:
        assert (file1 in res_file_names), f'{file1} missing'

    assert grep_exists(res_files['out'], "it=  1  is completed")
    spinup_trace_atom1 = grep_number(res_files['out'], r"l--> 2 Contour\(default\)    Spin-Up trace:", ":", res_index=-2)
    spindn_trace_atom1 = grep_number(res_files['out'], r"l--> 2 Contour\(default\)    Spin-Down trace:", ":", res_index=-2)
    spinoffdx_trace_atom1 = grep_number(res_files['out'], r"l--> 2 Contour\(default\)    Spin-Offd trace \(x\):", ":", res_index=-2)
    spinoffdy_trace_atom1 = grep_number(res_files['out'], r"l--> 2 Contour\(default\)    Spin-Offd trace \(y\):", ":", res_index=-2)


    spinup_trace_atom2 = grep_number(res_files['out'], r"l--> 2 Contour\(default\)    Spin-Up trace:", ":")
    spindn_trace_atom2 = grep_number(res_files['out'], r"l--> 2 Contour\(default\)    Spin-Down trace:", ":")
    spinoffdx_trace_atom2 = grep_number(res_files['out'], r"l--> 2 Contour\(default\)    Spin-Offd trace \(x\):", ":")
    spinoffdy_trace_atom2 = grep_number(res_files['out'], r"l--> 2 Contour\(default\)    Spin-Offd trace \(y\):", ":")

    assert abs(spinup_trace_atom1 - 4.5242) <= 0.0005
    assert abs(spindn_trace_atom1 - 2.3679) <= 0.0005
    assert abs(spinoffdx_trace_atom1 - 0.0000) <= 0.0005
    assert abs(spinoffdy_trace_atom1 - 0.5994) <= 0.0005

    assert abs(spinup_trace_atom2 - 4.5214) <= 0.0005
    assert abs(spindn_trace_atom2 - 2.38) <= 0.0005
    assert abs(spinoffdx_trace_atom2 + 0.6027) <= 0.0005
    assert abs(spinoffdy_trace_atom2 - 0.0000) <= 0.0005