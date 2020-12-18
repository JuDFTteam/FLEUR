

def test_CuBulkXML(execute_fleur, grep_exists, grep_number):
    """
    Simple test of Fleur with XML input with one step:
    1.Generate a starting density and run 1 iteration and compare fermi-energy & total energy
    """
    test_file_folder = './inputfiles/inpgen/Si_full_para/files/'

    cmd_params = []
    res_files = execute_fleur(cmd_params, test_file_folder)
    should_files = ['out.xml', 'out', 'usage.json']
    res_file_names = list(res_files.keys())
    
    # Test if all files are there
    for file1 in should_files:
        assert file1 in res_file_names
    
    assert grep_exists(res_files['out'], "it=  1  is completed")
    fermi = grep_number(res_files['out'], "new fermi energy", ":")
    tenergy = grep_number(res_files['out'], "total energy=", "=")
    dist = grep_number(res_files['out'], "distance of charge densitie", "1:")
    
    assert abs(fermi - 0.4233) <= 0.001
    assert abs(tenergy + 3305.016) <= 0.001
    assert abs(dist - 45.8259) <= 0.001


    assert False

'''
def test_CuBandXML():
    """
    """
    assert False

def test_CuBulkLibXC():
    """
    """
    assert False

def test_CuDOSXML():
    """
    """
    assert False
'''
