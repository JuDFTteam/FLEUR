
import pytest

@pytest.mark.bulk
@pytest.mark.magnetism
@pytest.mark.hdf
def DEACTIVATED_test_Fe_bcc_FlipcdnXLDA(execute_fleur, grep_number, grep_exists):
    """FeBCCFlipX: Flipcdn and noco in MT test
    Simple test of Fleur with two steps:
    1.Generate a rotated starting density and run 1 iteration
    2.Calculate magnetization and check it's value. 
    """
    test_file_folder = './inputfiles/Fe_bcc_FlipcdnXLDA/'
 
    # Stage 1
    files = ['inp.xml', 'kpts.xml', 'sym.xml', 'JUDFT_WARN_ONLY']
    res_files = execute_fleur(test_file_folder, only_copy=files)
    should_files = ['out']
    res_file_names = list(res_files.keys())
    for file1 in should_files:
        assert file1 in res_file_names
    
    assert grep_exists(res_files['out'], "flip")

    # Stage 2
    res_files = execute_fleur(test_file_folder, only_copy=[['inp2.xml', 'inp.xml', 'JUDFT_WARN_ONLY']])
    mx = grep_number(res_files['out'], "mx=", "mx=")
    assert abs(mx - 2.116) <= 0.001


@pytest.mark.bulk
@pytest.mark.magnetism
@pytest.mark.hdf
def DEACTIVATED_test_Fe_bcc_FlipcdnYGGA(execute_fleur, grep_number, grep_exists):
    """FeBCCFlipY: Flipcdn and noco in MT test
    Simple test of Fleur with two steps:
    1.Generate a rotated starting density and run 1 iteration
    2.Calculate magnetization and check it's value. 
    """
    test_file_folder = './inputfiles/Fe_bcc_FlipcdnYGGA/'

    # Stage 1
    files = ['inp.xml', 'kpts.xml', 'sym.xml', 'JUDFT_WARN_ONLY']
    res_files = execute_fleur(test_file_folder, only_copy=files)
    should_files = ['out']
    res_file_names = list(res_files.keys())
    for file1 in should_files:
        assert file1 in res_file_names

    assert grep_exists(res_files['out'], "flip")

    # Stage 2
    res_files = execute_fleur(test_file_folder, only_copy=[['inp2.xml', 'inp.xml'], 'JUDFT_WARN_ONLY'])
    my = grep_number(res_files['out'], "my=", "my=")
    assert abs(my + 2.180) <= 0.001

@pytest.mark.bulk
@pytest.mark.magnetism
@pytest.mark.hdf
def test_Fe_bcc_SF_LDA(execute_fleur, grep_number, grep_exists):
    """FeBCCSFLDA: Sourcefree magnetism and magnetization scaling

    Simple test of Fleur with one step:
    1.Generate a starting density in z-direction, calculate the potential and make it sourcefree. 
    Check for the correct transformed magnetic field. And the correct resulting magnetization.
    """
    test_file_folder = './inputfiles/Fe_bcc_SF_LDA/'

    res_files = execute_fleur(test_file_folder)
    should_files = ['out']
    res_file_names = list(res_files.keys())
    for file1 in should_files:
        assert file1 in res_file_names

    assert grep_exists(res_files['out'], "SF")
    magmom = grep_number(res_files['out'], "Magmom before SF", "mz=")
    bfield_before = grep_number(res_files['out'], "Bfield before SF", "Bz=")
    bfield_after = grep_number(res_files['out'], "Bfield after SF", "Bz=")
    mz = grep_number(res_files['out'], "local frame: mx=", "mz=")

    assert abs(magmom - 2.37929) <= 0.001
    assert abs(bfield_before - -1.17291) <= 0.001
    assert abs(bfield_after - -1.06968) <= 0.001
    assert abs(mz - 1.98655) <= 0.001

@pytest.mark.bulk
@pytest.mark.magnetism
@pytest.mark.soc
def test_Fe_bcc_orbital_polarization_correction(execute_fleur, grep_number, grep_exists):
    """
    Test of the orbital polarization correction

    Simple Fleur test with one step:
    1.Run 12 iterations with a orbital polarization correction added on the 3d states and
      check for the correct orbital moment
    """
    test_file_folder = './inputfiles/Fe_bcc_OPC/'

    res_files = execute_fleur(test_file_folder)
    should_files = ['out']
    res_file_names = list(res_files.keys())
    for file1 in should_files:
        assert file1 in res_file_names

    assert grep_exists(res_files['out'], "orb. magnetic moments in the spheres:")
    orbital_mom1 = grep_number(res_files['out'], "--> mm       1")
    orbital_mom2 = grep_number(res_files['out'], "--> mm       2")

    assert abs(orbital_mom1 - 0.07252) <= 0.001
    assert abs(orbital_mom2 - 0.07252) <= 0.001

