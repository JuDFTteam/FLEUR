
import pytest

@pytest.mark.xml
@pytest.mark.ldau
def test_NiO_ldauXML(execute_fleur, grep_number, grep_exists):
    """NiO with LDA+U (XML codepath)


    Simple test of Fleur with one step:
    1.Generate a starting density and run 9 iterations and compare convergence, fermi-energy & total energy
    """
    test_file_folder = './inputfiles/NiO_ldauXML/'
    res_files = execute_fleur(test_file_folder)
    res_file_names = list(res_files.keys())
    should_files = ['out']
    for file1 in should_files:
        assert (file1 in res_file_names), f'{file1} missing'
    
    assert grep_exists(res_files['out'], "it=  9  is completed")
    fermi = grep_number(res_files['out'], "new fermi energy", ":")
    tenergy = grep_number(res_files['out'], "total energy=", "=")
    dist = grep_number(res_files['out'], "distance of charge densities for it=  9", ":")
    dist = grep_number(res_files['out'], "mm       1", "mm       1")

    assert abs(fermi - 0.262) <= 0.001
    assert abs(tenergy - -3191.938) <= 0.001
    assert abs(dist - 0.25) <= 0.01
    assert abs(dist - 1.72) <= 0.03