
import pytest

@pytest.mark.hybrid
@pytest.mark.xml
def test_SiHybrid8kpt_nosym(execute_fleur, grep_number, grep_exists):
    """Fleur Si Hybrid Gamma - XML

    Simple test of Fleur with XML input with one step:
    1.Generate a starting density and run 1 iteration and compare fermi-energy & total energy
    """
    test_file_folder = './inputfiles/SiHybrid8kpt_nosym/'
    # What is testrun_seq in the original test?
    res_files = execute_fleur(test_file_folder)
    res_file_names = list(res_files.keys())
    should_files = ['out']
    for file1 in should_files:
        assert (file1 in res_file_names), f'{file1} missing'
    
    assert grep_exists(res_files['out'], "it=  1  is completed")
    tenergy = grep_number(res_files['out'], "    total energy=", "=")
    hfenergy = grep_number(res_files['out'], "HF total energy=", "=")
    dist13 = grep_number(res_files['out'], "distance of charge densities for it  13", ":")
    dist14 = grep_number(res_files['out'], "distance of charge densities for it  14", ":")
    bandgap = grep_number(res_files['out'], "bandgap", ":")

    assert abs(tenergy - -569.8053210563) <= 0.0001
    assert abs(hfenergy - -580.0376501381) <= 0.0001
    assert abs(dist13 - 0.000002) <= 0.0001
    assert abs(dist14 - 0.512727) <= 0.0001
    assert abs(bandgap - 0.089933) <= 0.000001


@pytest.mark.hybrid
def test_SiHybrid8kpt_sym(execute_fleur, grep_number, grep_exists):
    """Fleur Si Hybrid Gamma - XML

    Simple test of Fleur with XML input with one step:
    1.Generate a starting density and run 1 iteration and compare fermi-energy & total energy
    """
    test_file_folder = './inputfiles/SiHybrid8kpt_sym/'
    res_files = execute_fleur(test_file_folder)
    res_file_names = list(res_files.keys())
    should_files = ['out']
    for file1 in should_files:
        assert (file1 in res_file_names), f'{file1} missing'
    
    assert grep_exists(res_files['out'], "it=  1  is completed")
    tenergy = grep_number(res_files['out'], "    total energy=", "=")
    hfenergy = grep_number(res_files['out'], "HF total energy=", "=")
    dist13 = grep_number(res_files['out'], "distance of charge densities for it  13", ":")
    dist14 = grep_number(res_files['out'], "distance of charge densities for it  14", ":")
    bandgap = grep_number(res_files['out'], "bandgap", ":")

    assert abs(tenergy - -569.8004567175) <= 0.0001
    assert abs(hfenergy - -580.0362470542) <= 0.0001
    assert abs(dist13 - 0.000002) <= 0.0001
    assert abs(dist14 - 0.508222) <= 0.0001
    assert abs(bandgap - 0.088478) <= 0.000001

@pytest.mark.hybrid
def test_SiHybridGammaNoInv(execute_fleur, grep_number, grep_exists):
    """Fleur Si Hybrid Gamma - XML

    Simple test of Fleur with XML input with one step:
    1.Generate a starting density and run 1 iteration and compare fermi-energy & total energy
    """
    test_file_folder = './inputfiles/SiHybridGammaNoInv/'
    res_files = execute_fleur(test_file_folder)
    res_file_names = list(res_files.keys())
    should_files = ['out']
    for file1 in should_files:
        assert (file1 in res_file_names), f'{file1} missing'
    
    assert grep_exists(res_files['out'], "it=  1  is completed")
    tenergy = grep_number(res_files['out'], "    total energy=", "=")
    hfenergy = grep_number(res_files['out'], "HF total energy=", "=")
    dist23 = grep_number(res_files['out'], "distance of charge densities for it  23", ":")
    dist24 = grep_number(res_files['out'], "distance of charge densities for it  24", ":")
    bandgap = grep_number(res_files['out'], "bandgap", ":")

    assert abs(tenergy - -579.5474574652) <= 0.0001
    assert abs(hfenergy - -579.5486158288) <= 0.0001
    assert abs(dist23 - 0.000007) <= 0.0001
    assert abs(dist24 - 0.050063) <= 0.0001
    assert abs(bandgap - 0.191654) <= 0.000001