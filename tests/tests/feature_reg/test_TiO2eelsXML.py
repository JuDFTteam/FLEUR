import pytest

@pytest.mark.eels
@pytest.mark.xml
@pytest.mark.skip('Test needs to be updated. Outdated input file')
def test_TiO2eelsXML(execute_fleur, grep_number, grep_exists):
    """Fleur TiO2 EELS spectrum XML

    Simple test of Fleur with two steps:
    1.Generate a starting density and run a single iteration
    2.Generate EELS spectrum and compare some values
    """
    test_file_folder = './inputfiles/TiO2eelsXML/'

    # Stage 1
    res_files = execute_fleur(test_file_folder, only_copy=['inp.xml'])
    res_file_names = list(res_files.keys())
    should_files = ['out']
    for file1 in should_files:
        assert (file1 in res_file_names), f'{file1} missing'
    assert ('cdn.hdf' in res_file_names) or ('cdn1' in res_file_names)

    # Stage 2
    res_files = execute_fleur(test_file_folder, only_copy=[['inp2.xml', 'inp.xml']])
    res_file_names = list(res_files.keys())
    should_files = ['out', 'fort.37']
    for file1 in should_files:
        assert (file1 in res_file_names), f'{file1} missing'

    assert grep_exists(res_files['fort.37'], "435.09")
    assert grep_exists(res_files['fort.37'], "1.917")
    assert grep_exists(res_files['fort.37'], "438.09")
    assert grep_exists(res_files['fort.37'], "1.050")
    assert grep_exists(res_files['fort.37'], "2.236")
    assert grep_exists(res_files['fort.37'], "3.299")
    assert grep_exists(res_files['fort.37'], "    0    1 430.09")
    assert grep_exists(res_files['fort.37'], "    1    1 446.09")
