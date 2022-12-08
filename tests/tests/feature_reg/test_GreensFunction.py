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
    assert abs(spindown_trace - 2.64494) <= 0.0005


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

    # Stage 3 with only a HDLO on the orbital

    res_files = execute_fleur(test_file_folder, only_copy=[['inp-3.xml', 'inp.xml'], 'JUDFT_WARN_ONLY'])
    res_file_names = list(res_files.keys())
    should_files = ['out']
    if 'cdn.hdf' in res_file_names:
        should_files.append('greensf.hdf')
    for file1 in should_files:
        assert (file1 in res_file_names), f'{file1} missing'

    # now test output (Occupations are different because of different LO setup and it is the second iteration but consistent with mt charges)

    assert grep_exists(res_files['out'], "it=  1  is completed")
    cont_reg_sup2 = grep_number(res_files['out'], r"l--> 2 Contour\(default\)    Spin-Up trace:", ":", res_index=-2)
    cont_reg_sdown2 = grep_number(res_files['out'], r"l--> 2 Contour\(default\)    Spin-Down trace:", ":", res_index=-2)
    cont_reg_sup1 = grep_number(res_files['out'], r"l--> 2 Contour\(default\)    Spin-Up trace:", ":")
    cont_reg_sdown1 = grep_number(res_files['out'], r"l--> 2 Contour\(default\)    Spin-Down trace:", ":")

    assert abs(cont_reg_sup2 - 4.1237) <= 0.0005
    assert abs(cont_reg_sdown2 - 1.9298) <= 0.0005
    assert abs(cont_reg_sup1 - 4.1237) <= 0.0005
    assert abs(cont_reg_sdown1 - 1.9298) <= 0.0005

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


    assert grep_exists(res_files['out'], r'0\.9918 [\s\-]0\.0000')
    assert grep_exists(res_files['out'], r'0\.9907 [\s\-]0\.0000')
    assert grep_exists(res_files['out'], r'0\.4191 [\s\-]0\.0000')
    assert grep_exists(res_files['out'], r'0\.6783 [\s\-]0\.0000')
    assert grep_exists(res_files['out'], r'0\.1895 [\s\-]0\.0000')
    assert grep_exists(res_files['out'], r'0\.3421 [\s\-]0\.0000')
    spinup_trace = grep_number(res_files['out'], "Spin-Up trace:", ": ")
    spindown_trace = grep_number(res_files['out'], "Spin-Down trace:", ": ")

    assert abs(spinup_trace - 6.9375) <= 0.0005
    assert abs(spindown_trace - 3.9556) <= 0.0005


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

    assert grep_exists(res_files['out'], r'0\.9918 [\s\-]0\.0000')
    assert grep_exists(res_files['out'], r'0\.9907 [\s\-]0\.0000')
    assert grep_exists(res_files['out'], r'0\.4191 [\s\-]0\.0000')
    assert grep_exists(res_files['out'], r'0\.6783 [\s\-]0\.0000')
    assert grep_exists(res_files['out'], r'[\s\-]0\.0000  0\.1895')
    assert grep_exists(res_files['out'], r'[\s\-]0\.0000 -0\.1895')
    assert grep_exists(res_files['out'], r'[\s\-]0\.0000  0\.3421')
    assert grep_exists(res_files['out'], r'[\s\-]0\.0000 -0\.3421')
    spinup_trace = grep_number(res_files['out'], "Spin-Up trace:", ": ")
    spindown_trace = grep_number(res_files['out'], "Spin-Down trace:", ": ")

    assert abs(spinup_trace - 6.9375) <= 0.0005
    assert abs(spindown_trace - 3.9556) <= 0.0005

@pytest.mark.disabled
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

    assert grep_exists(res_files['out'], r'0\.9918 [\s\-]0\.0000')
    assert grep_exists(res_files['out'], r'0\.9907 [\s\-]0\.0000')
    assert grep_exists(res_files['out'], r'0\.4191 [\s\-]0\.0000')
    assert grep_exists(res_files['out'], r'0\.6793 [\s\-]0\.0000')
    assert grep_exists(res_files['out'], r'[\s\-]0\.0000  0\.1885')
    assert grep_exists(res_files['out'], r'[\s\-]0\.0000 -0\.1885')
    assert grep_exists(res_files['out'], r'[\s\-]0\.0000  0\.3424')
    assert grep_exists(res_files['out'], r'[\s\-]0\.0000 -0\.3424')
    spinup_trace = grep_number(res_files['out'], "Spin-Up trace:", ": ")
    spindown_trace = grep_number(res_files['out'], "Spin-Down trace:", ": ")

    assert abs(spinup_trace - 6.9375) <= 0.0005
    assert abs(spindown_trace - 3.9609) <= 0.0005

@pytest.mark.serial
@pytest.mark.greensfunction
@pytest.mark.non_collinear
@pytest.mark.bulk
def test_GreensFunction_mperp_xdir(execute_fleur, grep_number, grep_exists):
    """Fleur Fe bcc Green's function Radial Noco spin offdiagonal

    Simple test of the green's function calculation for spin-offdiagonal components in FLEUR with one step:
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
    assert abs(spindn_trace_atom1 - 2.3726) <= 0.0005
    assert abs(spinoffdx_trace_atom1 - 0.601) <= 0.0005
    assert abs(spinoffdy_trace_atom1 - 0.0000) <= 0.0005

    assert abs(spinup_trace_atom2 - 4.5242) <= 0.0005
    assert abs(spindn_trace_atom2 - 2.3726) <= 0.0005
    assert abs(spinoffdx_trace_atom2 + 0.601) <= 0.0005
    assert abs(spinoffdy_trace_atom2 - 0.0000) <= 0.0005


@pytest.mark.disabled
@pytest.mark.serial
@pytest.mark.greensfunction
@pytest.mark.non_collinear
@pytest.mark.bulk
def test_GreensFunction_mperp_ydir(execute_fleur, grep_number, grep_exists):
    """Fleur Fe bcc Green's function Radial Noco spin offdiagonal

    Simple test of the green's function calculation for spin-offdiagonal components in FLEUR with one step:
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
    assert abs(spindn_trace_atom1 - 2.3726) <= 0.0005
    assert abs(spinoffdx_trace_atom1 - 0.0000) <= 0.0005
    assert abs(spinoffdy_trace_atom1 + 0.601) <= 0.0005

    assert abs(spinup_trace_atom2 - 4.5242) <= 0.0005
    assert abs(spindn_trace_atom2 - 2.3726) <= 0.0005
    assert abs(spinoffdx_trace_atom2 + 0.601) <= 0.0005
    assert abs(spinoffdy_trace_atom2 - 0.0000) <= 0.0005


@pytest.mark.serial
@pytest.mark.greensfunction
@pytest.mark.bulk
def test_GreensFunction_InterOrbital(execute_fleur, grep_number, grep_exists):
    """Fleur GdCu Green's function interorbital elements

    Simple test of the green's function calculation for l-offdiagonal components in FLEUR with one step:
    1. Generate starting density, run 1 Iteration and calculate Green's function
       for d-orbitals. Ensure that the obtained occupation matrices look as expected
    """
    test_file_folder = './inputfiles/GreensFunction_InterOrbital/'
    res_files = execute_fleur(test_file_folder, only_copy=['inp.xml', 'JUDFT_WARN_ONLY'])
    res_file_names = list(res_files.keys())
    should_files = ['out']
    if 'cdn.hdf' in res_file_names:
        should_files.append('greensf.hdf')
    for file1 in should_files:
        assert (file1 in res_file_names), f'{file1} missing'

    assert grep_exists(res_files['out'], "it=  1  is completed")
    assert grep_exists(res_files['out'], r"Green's Function Elements: 9\s")

    #Check some matrix elements of the occupation matrix (f-p)
    assert grep_exists(res_files['out'], r" 0\.0008 [\s\-]0\.0000")
    assert grep_exists(res_files['out'], r"-0\.0013 [\s\-]0\.0000")
    assert grep_exists(res_files['out'], r" 0\.0010 [\s\-]0\.0000")
    assert grep_exists(res_files['out'], r"-0\.0021 [\s\-]0\.0000")
    assert grep_exists(res_files['out'], r" 0\.0034 [\s\-]0\.0000")
    assert grep_exists(res_files['out'], r"-0\.0027 [\s\-]0\.0000")


@pytest.mark.bulk
@pytest.mark.greensfunction
@pytest.mark.magnetism
@pytest.mark.serial
def test_GreensFunction_IntersiteSingleShell(execute_fleur, grep_number, grep_exists):
    """Fleur Greens Function intersite single shell
    Simple test of the intersite green's function calculation in FLEUR with one step:
    1. Generate starting density, run 1 Iteration and calculate Green's function
       for d-orbitals for the nearest neighbours. Ensure that the obtained occupation matrices
       look as expected
    """
    test_file_folder = './inputfiles/GreensFunction_IntersiteSingleShell/'

    res_files = execute_fleur(test_file_folder)
    should_files = ['out']
    res_file_names = list(res_files.keys())
    if 'cdn.hdf' in res_file_names:
        should_files.append('greensf.hdf')
    for file1 in should_files:
        assert file1 in res_file_names

    assert grep_exists(res_files['out'], "it=  1  is completed")

    #Check for the right shell being selected
    assert grep_exists(res_files['out'], r"Green's Function Elements: 9\s")
    #These are entries in the table of generated GF elements
    assert grep_exists(res_files['out'], r"2 [|] 2/2  [|]    1/    1 [|]       1 [|]      T [|]         1 [|]        \-1\(\-1\)      [|]      F\(F\)  [|]  1.000  1.000  1.000")
    assert grep_exists(res_files['out'], r"3 [|] 2/2  [|]    1/    1 [|]       1 [|]      T [|]         1 [|]         2\( 2\)      [|]      F\(F\)  [|] \-1.000 [\-\s]0.000 [\-\s]0.000")

    #Check trace of representative element
    spinup_trace_element1 = grep_number(res_files['out'], r"l--> 2 Contour\(default\)    Spin-Up trace:", ":", res_index=-8)
    spindn_trace_element1 = grep_number(res_files['out'], r"l--> 2 Contour\(default\)    Spin-Down trace:", ":", res_index=-8)
    #Check trace of a symmetry equivalent element
    spinup_trace_element2 = grep_number(res_files['out'], r"l--> 2 Contour\(default\)    Spin-Up trace:", ":")
    spindn_trace_element2 = grep_number(res_files['out'], r"l--> 2 Contour\(default\)    Spin-Down trace:", ":")

    assert abs(spinup_trace_element1 - 0.1484) <= 0.0005
    assert abs(spindn_trace_element1 - -0.0285) <= 0.0005

    assert abs(spinup_trace_element2 - 0.1484) <= 0.0005
    assert abs(spindn_trace_element2 - -0.0285) <= 0.0005

    #Check some matrix elements of the occupation matrix
    assert grep_exists(res_files['out'], " 0.1884 -0.1884")
    assert grep_exists(res_files['out'], " 0.1884  0.1884")
    assert grep_exists(res_files['out'], "-0.1884  0.1884")
    assert grep_exists(res_files['out'], "-0.1884 -0.1884")

    assert grep_exists(res_files['out'], " 0.0717 -0.0717")
    assert grep_exists(res_files['out'], " 0.0717  0.0717")
    assert grep_exists(res_files['out'], "-0.0717  0.0717")
    assert grep_exists(res_files['out'], "-0.0717 -0.0717")

    assert grep_exists(res_files['out'], r" 0\.0015 [\-\s]0\.0000")
    assert grep_exists(res_files['out'], r" 0\.0485 [\-\s]0\.0000")
    assert grep_exists(res_files['out'], r"\-0\.1951 [\-\s]0\.0000")
    assert grep_exists(res_files['out'], r" 0\.1206 [\-\s]0\.0000")


@pytest.mark.parametrize('filename',['inp_nogamma.xml','inp_gamma.xml'])
@pytest.mark.bulk
@pytest.mark.greensfunction
@pytest.mark.magnetism
@pytest.mark.serial
def test_GreensFunction_GammaNoGamma(execute_fleur, grep_number, grep_exists, filename):
    """Fleur Greens Function intersite single shell
    Simple test of the intersite green's function calculation in FLEUR with one step:
    1. Generate starting density, run 1 Iteration and calculate Green's function
       for d-orbitals for the nearest neighbours. Ensure that the obtained occupation matrices
       look as expected
    """
    test_file_folder = './inputfiles/GreensFunction_IntersiteGammaNoGamma/'

    res_files = execute_fleur(test_file_folder,only_copy=[[filename, 'inp.xml']])
    should_files = ['out']
    res_file_names = list(res_files.keys())
    if 'cdn.hdf' in res_file_names:
        should_files.append('greensf.hdf')
    for file1 in should_files:
        assert file1 in res_file_names

    assert grep_exists(res_files['out'], "it=  1  is completed")

    #Check for the right shell being selected
    assert grep_exists(res_files['out'], r"Green's Function Elements: 18\s")
    #These are entries in the table of generated GF elements
    assert grep_exists(res_files['out'], r"3 [|] 2/2  [|]    1/    2 [|]       1 [|]      T [|]         1 [|]        \-1\(\-1\)      [|]      F\(F\)  [|]  0.500  0.500  0.500")
    assert grep_exists(res_files['out'], r"4 [|] 2/2  [|]    1/    2 [|]       1 [|]      T [|]         1 [|]        \-1\(\-1\)      [|]      F\(F\)  [|]  0.500  0.500 \-0.500")
    assert grep_exists(res_files['out'], r"11 [|] 2/2  [|]    2/    1 [|]       1 [|]      T [|]         2 [|]        \-1\(\-1\)      [|]      F\(F\)  [|]  0.500  0.500  0.500")
    assert grep_exists(res_files['out'], r"12 [|] 2/2  [|]    2/    1 [|]       1 [|]      T [|]         2 [|]        \-1\(\-1\)      [|]      F\(F\)  [|]  0.500  0.500 \-0.500")


    #Check trace of representative element
    spinup_trace_element1 = grep_number(res_files['out'], r"l--> 2 Contour\(default\)    Spin-Up trace:", ":", res_index=-16)
    spindn_trace_element1 = grep_number(res_files['out'], r"l--> 2 Contour\(default\)    Spin-Down trace:", ":", res_index=-16)
    #Check trace of a symmetry equivalent element
    spinup_trace_element2 = grep_number(res_files['out'], r"l--> 2 Contour\(default\)    Spin-Up trace:", ":")
    spindn_trace_element2 = grep_number(res_files['out'], r"l--> 2 Contour\(default\)    Spin-Down trace:", ":")

    assert abs(spinup_trace_element1 - 0.2422) <= 0.0005
    assert abs(spindn_trace_element1 - 0.0351) <= 0.0005

    assert abs(spinup_trace_element2 - 0.2422) <= 0.0005
    assert abs(spindn_trace_element2 - 0.0351) <= 0.0005

    #Check some matrix elements of the occupation matrix
    assert grep_exists(res_files['out'], " 0.0611  0.0611")
    assert grep_exists(res_files['out'], " 0.0611 -0.0611")
    assert grep_exists(res_files['out'], "-0.0611  0.0611")
    assert grep_exists(res_files['out'], "-0.0611 -0.0611")

    assert grep_exists(res_files['out'], " 0.0599  0.0599")
    assert grep_exists(res_files['out'], " 0.0599 -0.0599")
    assert grep_exists(res_files['out'], "-0.0599  0.0599")
    assert grep_exists(res_files['out'], "-0.0599 -0.0599")

    assert grep_exists(res_files['out'], r" 0\.0142 [\-\s]0\.0000")
    assert grep_exists(res_files['out'], r" 0\.0713 [\-\s]0\.0000")
    assert grep_exists(res_files['out'], r"\-0\.2260 [\-\s]0\.0000")
    assert grep_exists(res_files['out'], r" 0\.1624 [\-\s]0\.0000")




@pytest.mark.bulk
@pytest.mark.greensfunction
@pytest.mark.magnetism
@pytest.mark.serial
def test_GreensFunction_IntersiteMultipleShells(execute_fleur, grep_number, grep_exists):
    """Fleur Greens Function intersite multiple shells
    Simple test of the intersite green's function calculation in FLEUR with one step:
    1. Generate starting density, run 1 Iteration and calculate Green's function
       for d-orbitals for the 5 nearest neighbours. This forces the shell construction algorithm to
       do a bit more work. Ensure that the obtained occupation matrices
       look as expected
    """
    test_file_folder = './inputfiles/GreensFunction_IntersiteMultipleShells/'

    res_files = execute_fleur(test_file_folder)
    should_files = ['out']
    res_file_names = list(res_files.keys())
    if 'cdn.hdf' in res_file_names:
        should_files.append('greensf.hdf')
    for file1 in should_files:
        assert file1 in res_file_names

    assert grep_exists(res_files['out'], "it=  1  is completed")

    #Check for the right shell being selected
    assert grep_exists(res_files['out'], r"Green's Function Elements: 59\s")
    #These are entries in the table of generated GF elements
    assert grep_exists(res_files['out'], r"2 [|] 2/2  [|]    1/    1 [|]       1 [|]      T [|]         1 [|]        \-1\(\-1\)      [|]      F\(F\)  [|]  1.000  1.000  1.000")
    assert grep_exists(res_files['out'], r"3 [|] 2/2  [|]    1/    1 [|]       1 [|]      T [|]         1 [|]         2\( 2\)      [|]      F\(F\)  [|] \-1.000 [\-\s]0.000 [\-\s]0.000")
    assert grep_exists(res_files['out'], r"52 [|] 2/2  [|]    1/    1 [|]       1 [|]      T [|]         1 [|]        -1\(-1\)      [|]      F\(F\)  [|]  2.000  2.000  2.000")
    assert grep_exists(res_files['out'], r"59 [|] 2/2  [|]    1/    1 [|]       1 [|]      T [|]         1 [|]        52\(22\)      [|]      F\(F\)  [|] [\-\s]0.000 [\-\s]0.000  2.000")

    #Check first shell again
    #Check trace of representative element
    spinup_trace_element1 = grep_number(res_files['out'], r"l--> 2 Contour\(default\)    Spin-Up trace:", ":", res_index=-58)
    spindn_trace_element1 = grep_number(res_files['out'], r"l--> 2 Contour\(default\)    Spin-Down trace:", ":", res_index=-58)
    #Check trace of a symmetry equivalent element
    spinup_trace_element2 = grep_number(res_files['out'], r"l--> 2 Contour\(default\)    Spin-Up trace:", ":",res_index=-51)
    spindn_trace_element2 = grep_number(res_files['out'], r"l--> 2 Contour\(default\)    Spin-Down trace:", ":",res_index=-51)

    assert abs(spinup_trace_element1 - 0.1484) <= 0.0005
    assert abs(spindn_trace_element1 - -0.0285) <= 0.0005

    assert abs(spinup_trace_element2 - 0.1484) <= 0.0005
    assert abs(spindn_trace_element2 - -0.0285) <= 0.0005

    #Check some matrix elements of the occupation matrix (Next to the corner m=3 mp=2 and diagonal elements)
    assert grep_exists(res_files['out'], " 0.0717 -0.0717")
    assert grep_exists(res_files['out'], " 0.0717  0.0717")
    assert grep_exists(res_files['out'], "-0.0717  0.0717")
    assert grep_exists(res_files['out'], "-0.0717 -0.0717")

    assert grep_exists(res_files['out'], " 0.1884 -0.1884")
    assert grep_exists(res_files['out'], " 0.1884  0.1884")
    assert grep_exists(res_files['out'], "-0.1884  0.1884")
    assert grep_exists(res_files['out'], "-0.1884 -0.1884")

    assert grep_exists(res_files['out'], r" 0\.0015 [\-\s]0\.0000")
    assert grep_exists(res_files['out'], r" 0\.0485 [\-\s]0\.0000")
    assert grep_exists(res_files['out'], r"\-0\.1951 [\-\s]0\.0000")
    assert grep_exists(res_files['out'], r" 0\.1206 [\-\s]0\.0000")

    #Check outermost shell
    #Check trace of representative element
    spinup_trace_element1 = grep_number(res_files['out'], r"l--> 2 Contour\(default\)    Spin-Up trace:", ":", res_index=-8)
    spindn_trace_element1 = grep_number(res_files['out'], r"l--> 2 Contour\(default\)    Spin-Down trace:", ":", res_index=-8)
    #Check trace of a symmetry equivalent element
    spinup_trace_element2 = grep_number(res_files['out'], r"l--> 2 Contour\(default\)    Spin-Up trace:", ":")
    spindn_trace_element2 = grep_number(res_files['out'], r"l--> 2 Contour\(default\)    Spin-Down trace:", ":")

    assert abs(spinup_trace_element1 -  0.0242) <= 0.0005
    assert abs(spindn_trace_element1 - -0.1257) <= 0.0005

    assert abs(spinup_trace_element2 -  0.0242) <= 0.0005
    assert abs(spindn_trace_element2 - -0.1257) <= 0.0005

    #Check some matrix elements of the occupation matrix
    #Check some matrix elements of the occupation matrix (Next to the corner m=3 mp=2 and diagonal elements)
    assert grep_exists(res_files['out'], " 0.0296 -0.0296")
    assert grep_exists(res_files['out'], " 0.0296  0.0296")
    assert grep_exists(res_files['out'], "-0.0296  0.0296")
    assert grep_exists(res_files['out'], "-0.0296 -0.0296")

    assert grep_exists(res_files['out'], " 0.0147 -0.0147")
    assert grep_exists(res_files['out'], " 0.0147  0.0147")
    assert grep_exists(res_files['out'], "-0.0147  0.0147")
    assert grep_exists(res_files['out'], "-0.0147 -0.0147")

    assert grep_exists(res_files['out'], r"\-0\.0024 [\-\s]0\.0000")
    assert grep_exists(res_files['out'], r" 0\.0157 [\-\s]0\.0000")
    assert grep_exists(res_files['out'], r"\-0\.0371 [\-\s]0\.0000")
    assert grep_exists(res_files['out'], r"\-0\.0172 [\-\s]0\.0000")


@pytest.mark.bulk
@pytest.mark.greensfunction
@pytest.mark.magnetism
@pytest.mark.serial
@pytest.mark.outxml_parser_xfail
def test_GreensFunction_IntersiteShellConstruction(execute_fleur, grep_exists):
    """Fleur Greens Function intersite shell construction
    Simple test of the intersite green's function initialization
    1. Generate the needed Green's function to calculate the J_ij for 60 shells at a different
       distance. Only run till after initialization (-check flag) and check that
       the Green's function table in the out file looks as expected
    """
    test_file_folder = './inputfiles/GreensFunction_IntersiteShellConstruction/'

    res_files = execute_fleur(test_file_folder,cmdline_param=['-check'], only_copy=['inp.xml']) #Only run the initializations
    should_files = ['out']
    res_file_names = list(res_files.keys())
    for file1 in should_files:
        assert file1 in res_file_names

    #Check for the right shell being selected
    assert grep_exists(res_files['out'], r"Green's Function Elements: 2445\s")
    #These are entries in the table of generated GF elements
    assert grep_exists(res_files['out'], r"2 [|] 2/2  [|]    1/    1 [|]       1 [|]      T [|]         1 [|]        \-1\(\-1\)      [|]      F\(F\)  [|]  1.000  1.000  1.000")
    assert grep_exists(res_files['out'], r"3 [|] 2/2  [|]    1/    1 [|]       1 [|]      T [|]         1 [|]         2\( 2\)      [|]      F\(F\)  [|] \-1.000 [\-\s]0.000  0.000")
    assert grep_exists(res_files['out'], r"2110 [|] 2/2  [|]    1/    1 [|]       1 [|]      T [|]         1 [|]        \-1\(\-1\)      [|]      F\(F\)  [|]  5.000 [\-\s]0.000 \-4.000")
    #This element was missed by the previous completeness detection
    assert grep_exists(res_files['out'], r"2125 [|] 2/2  [|]    1/    1 [|]       1 [|]      T [|]         1 [|]      2110\(31\)      [|]      F\(F\)  [|] \-9.000 \-5.000 \-5.000")
    assert grep_exists(res_files['out'], r"2422 [|] 2/2  [|]    1/    1 [|]       1 [|]      T [|]         1 [|]        \-1\(\-1\)      [|]      F\(F\)  [|]  4.000  4.000 \-4.000")
    assert grep_exists(res_files['out'], r"2445 [|] 2/2  [|]    1/    1 [|]       1 [|]      T [|]         1 [|]      2422\(46\)      [|]      F\(F\)  [|] \-8.000 \-8.000 \-4.000")

    res_files = execute_fleur(test_file_folder,cmdline_param=['-check'], only_copy=[['inp_start_from.xml', 'inp.xml']]) #Only run the initializations
    should_files = ['out']
    res_file_names = list(res_files.keys())
    for file1 in should_files:
        assert file1 in res_file_names

    #Check for the right shell being selected
    assert grep_exists(res_files['out'], r"Green's Function Elements: 45\s")
    #These are entries in the table of generated GF elements
    assert grep_exists(res_files['out'], r"2 [|] 2/2  [|]    1/    1 [|]       1 [|]      T [|]         1 [|]        \-1\(\-1\)      [|]      F\(F\)  [|]  1.000 [\-\s]0.000 \-1.000")
    assert grep_exists(res_files['out'], r"6 [|] 2/2  [|]    1/    1 [|]       1 [|]      T [|]         1 [|]         2\( 8\)      [|]      F\(F\)  [|]  1.000  1.000  2.000")
    assert grep_exists(res_files['out'], r"38 [|] 2/2  [|]    1/    1 [|]       1 [|]      T [|]         1 [|]        \-1\(\-1\)      [|]      F\(F\)  [|]  2.000  2.000  2.000")
    assert grep_exists(res_files['out'], r"45 [|] 2/2  [|]    1/    1 [|]       1 [|]      T [|]         1 [|]        38\(22\)      [|]      F\(F\)  [|] [\-\s]0.000 [\-\s]0.000  2.000")


@pytest.mark.bulk
@pytest.mark.greensfunction
@pytest.mark.magnetism
@pytest.mark.serial
@pytest.mark.outxml_parser_xfail
def test_GreensFunction_IntersiteShellConstructionFilm(execute_fleur, grep_exists):
    """Fleur Greens Function intersite shell construction
    Simple test of the intersite green's function initialization
    1. Generate the needed Green's function to calculate the J_ij for 60 shells at a different
       distance. Only run till after initialization (-check flag) and check that
       the Green's function table in the out file looks as expected
    """
    test_file_folder = './inputfiles/GreensFunction_IntersiteShellConstructionFilm/'

    res_files = execute_fleur(test_file_folder,cmdline_param=['-check'], only_copy=['inp.xml']) #Only run the initializations
    should_files = ['out']
    res_file_names = list(res_files.keys())
    for file1 in should_files:
        assert file1 in res_file_names

    #Check for the right shell being selected
    assert grep_exists(res_files['out'], r"Green's Function Elements: 373\s")
    #These are entries in the table of generated GF elements
    assert grep_exists(res_files['out'], r"2 [|] 2/2  [|]    1/    1 [|]       1 [|]      T [|]         1 [|]        \-1\(\-1\)      [|]      F\(F\)  [|]  1.000 [\-\s]0.000 [\-\s]0.000")
    assert grep_exists(res_files['out'], r"3 [|] 2/2  [|]    1/    1 [|]       1 [|]      T [|]         1 [|]         2\(\ 2\)      [|]      F\(F\)  [|] -1.000 [\-\s]0.000 [\-\s]0.000")
    assert grep_exists(res_files['out'], r"366 [|] 2/2  [|]    1/    1 [|]       1 [|]      T [|]         1 [|]        \-1\(\-1\)      [|]      F\(F\)  [|]  9.000  6.000 [\-\s]0.000")
    assert grep_exists(res_files['out'], r"373 [|] 2/2  [|]    1/    1 [|]       1 [|]      T [|]         1 [|]       366\( 8\)      [|]      F\(F\)  [|]  6.000  9.000 [\-\s]0.000")
