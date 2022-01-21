import pytest

@pytest.mark.soc
@pytest.mark.edsolver
def test_Gd_Hubbard1(execute_fleur, grep_number, grep_exists):
    """Fleur Gd Hubbard 1 SOC

    Simple test the Hubbard 1 method with SOC:
    1. Generate starting density, run 2 Iterations with one Hubbard iteration in between
        for f-orbitals. Ensure that the density matrix is reasonable
    """
    test_file_folder = './inputfiles/Gd_Hubbard1/'
    res_files = execute_fleur(test_file_folder)
    res_file_names = list(res_files.keys())
    should_files = ['out', 'Hubbard1/hubbard1.cfg', 'Hubbard1/hloc.cfg', 'Hubbard1/eigval7part.dat']
    for file1 in should_files:
        assert file1 in res_file_names
    
    assert grep_exists(res_files['out'], "it=  1  is completed")
    assert grep_exists(res_files['out'], "Hubbard 1 it=  1  is completed")

    xisoc = grep_number(res_files['Hubbard1/hloc.cfg'], "xiSOC", "xiSOC")

    mumatch = grep_number(res_files['out'], "muMatch =", ":")
    nmmp_occ = grep_number(res_files['out'], "nmmp occupation distance:", ":")
    nmmp_el = grep_number(res_files['out'], "nmmp element distance:", ":")   
    tenergy = grep_number(res_files['out'], "total energy=", "=")
    dist = grep_number(res_files['out'], "distance of spin densities for it=    3", ":")

    assert abs(xisoc - 0.21978) <= 0.00001
    assert abs(mumatch - -3.1063151E-002) <= 0.0001
    assert abs(nmmp_occ - 6.99999999668645) <= 0.00001
    assert abs(nmmp_el - 0.997526430566068) <= 0.00001
    assert abs(tenergy - -22560.5841855308) <= 0.00001
    assert abs(dist - 25.034287) <= 0.0001

@pytest.mark.soc
@pytest.mark.edsolver
def test_Gd_Hubbard1_noSYM(execute_fleur, grep_number, grep_exists):
    """Fleur Gd Hubbard 1 no Symmetry

    Simple test of the Hubbard 1 method without Symmetries 
    (and an addional onsite Green's Function on the d-orbitals to ensure the mapping is working correctly):
    1. Generate starting density, run 2 Iterations with one Hubbard iteration in between
       for f-orbitals. Ensure that the density matrix is reasonable
    """
    test_file_folder = './inputfiles/Gd_Hubbard1_noSYM/'
    res_files = execute_fleur(test_file_folder)
    res_file_names = list(res_files.keys())
    should_files = ['out', 'Hubbard1/1_f/hubbard1.cfg', 'Hubbard1/1_f/hloc.cfg', 'Hubbard1/1_f/eigval7part.dat',
    'Hubbard1/2_f/hubbard1.cfg', 'Hubbard1/2_f/hloc.cfg', 'Hubbard1/2_f/eigval7part.dat']
    for file1 in should_files:
        assert file1 in res_file_names
    
    assert grep_exists(res_files['out'], "it=  1  is completed")
    assert grep_exists(res_files['out'], "Hubbard 1 it=  1  is completed")

    nmmp_occ = grep_number(res_files['out'], "nmmp occupation distance:", ":")
    nmmp_el = grep_number(res_files['out'], "nmmp element distance:", ":")
    tenergy = grep_number(res_files['out'], "total energy=", "=")
    dist = grep_number(res_files['out'], "distance of spin densities for it=    3", ":")

    assert abs(nmmp_occ - 7.00000) <= 0.0005
    assert abs(nmmp_el - 0.9999) <= 0.0005
    assert abs(tenergy - -22560.55679) <= 0.00005
    assert abs(dist - 25.0582) <= 0.0002