
import pytest

@pytest.mark.serial
def test_Fe_Tetra_noSYM(execute_fleur, grep_number, grep_exists):
    """Fleur Fe tetra no symmetries

    Simple test of the linear tetrahedron method with bloechl corrections:
    1. Run 15 iterations with bloechl corrections and a small k-mesh. Then
       test for fermi energy, total energy and distance
    """
    test_file_folder = './inputfiles/Fe_Tetra_noSYM/'
    res_files = execute_fleur(test_file_folder)
    res_file_names = list(res_files.keys())
    should_files = ['out']
    for file1 in should_files:
        assert (file1 in res_file_names), f'{file1} missing'

    assert grep_exists(res_files['out'], "it= 15  is completed")
    fermi = grep_number(res_files['out'], "fermi energy =", "=")
    tenergy = grep_number(res_files['out'], "total energy=", "=")
    dist = grep_number(res_files['out'], "distance of spin densities for it=   15", ":")

    assert abs(fermi - 0.36646) <= 0.005
    assert abs(tenergy - -1272.8002015793) <= 0.00005
    assert abs(dist - 0.074) <= 0.002
