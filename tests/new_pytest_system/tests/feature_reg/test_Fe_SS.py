import pytest

@pytest.mark.bulk
@pytest.mark.magnetism
@pytest.mark.spinspiral
@pytest.mark.skip('Test is probably outdated (Was not used in long time)')
def test_fe_bulk_SS_conv(execute_fleur, grep_number, grep_exists):
    """Fe monolayer fcc (110): spin spiral energy

    Simple test of Fleur which converges a spin spiral in
    bulk Fe.
    This test was callibrated via MaX-R4 release because
    MaX-R4 was shown to give the same results as version 26d.
    Uses: spin spiral force theorem
    """
    test_file_folder = './inputfiles/Fe_bulk_SS_conv/'
    res_files = execute_fleur(test_file_folder)
    should_files = ['out']
    res_file_names = list(res_files.keys())
    for file1 in should_files:
        assert file1 in res_file_names
    
    tenergy = grep_number(res_files['out'], "total energy=", "=")
    mm = grep_number(res_files['out'], "mm       1", "mm       1")

    assert abs(tenergy - -1270.4698105523) <= 0.00001
    assert abs(mm - 3.44232) <= 0.00001


@pytest.mark.bulk
@pytest.mark.magnetism
@pytest.mark.spinspiral
@pytest.mark.skip('Test is probably outdated (Was not used in long time)')
def test_Fe_film_SS_conv(execute_fleur, grep_number, grep_exists):
    """Fe monolayer fcc (110): spin spiral energy
    
    Simple test of Fleur that converges a spin spiral.
    This test was callibrated via MaX-R4 release because
    MaX-R4 was shown to give the same results as version 26d.
    Uses: spin spiral force theorem
    """
    test_file_folder = './inputfiles/Fe_film_SS_conv/'
    res_files = execute_fleur(test_file_folder)
    should_files = ['out']
    res_file_names = list(res_files.keys())
    for file1 in should_files:
        assert file1 in res_file_names
    
    tenergy = grep_number(res_files['out'], "total energy=", "=")
    mm = grep_number(res_files['out'], "mm       1", "mm       1")

    assert abs(tenergy - -1270.5495629406) <= 0.00001
    assert abs(mm - 3.79655) <= 0.00001

@pytest.mark.bulk
@pytest.mark.magnetism
@pytest.mark.spinspiral
@pytest.mark.skip('Test is probably outdated (Was not used in long time)')
def test_Fe_film_SSFT(execute_fleur, grep_number, grep_exists, validate_out_xml_file):
    """Fe monolayer fcc (110): spin spiral energy
    
    Simple test of Fleur with two steps:
    1.Generate the reference charge density
    2.Calculate spin spiral energy via the force theorem
    This test was callibrated via MaX-R4 release because
    MaX-R4 was shown to give the same results as version 26d.
    Uses: spin spiral force theorem
    WARNING: the test is aimed to calculate energy difference
    between two spin spirals. If it starts to fail, please
    check if energy difference between q1 and q2 remains the
    same. If it is true, just adjust q1 and q2 values in test.run2.
    Otherwise there seems to be a bug.
    """
    test_file_folder = './inputfiles/Fe_film_SSFT/'

    # Stage 1
    res_files = execute_fleur(test_file_folder, only_copy=['inp.xml'])
    should_files = ['out']
    res_file_names = list(res_files.keys())
    for file1 in should_files:
        assert file1 in res_file_names
    
    tenergy = grep_number(res_files['out'], "total energy=", "=")
    assert abs(tenergy - -1270.42430970000) <= 0.00001

    # Stage 2
    files = [['inp2.xml', 'inp.xml'], 'JUDFT_WARN_ONLY']
    res_files = execute_fleur(test_file_folder, only_copy=files, rm_files='mixing_history')
    should_files = ['out.xml']
    res_file_names = list(res_files.keys())
    for file1 in should_files:
        assert file1 in res_file_names

    assert validate_out_xml_file(res_files['out.xml'])

    assert grep_exists(res_files['out.xml'], "Forcetheorem_SSDISP qvectors=")
    ev_q1 = grep_number(res_files['out.xml'], 'Entry q= "1', "ev-sum=.")
    ev_q2 = grep_number(res_files['out.xml'], 'Entry q= "2', "ev-sum=.")

    assert abs(ev_q1 - -2.2243925) <= 0.000001
    assert abs(ev_q2 - -2.2178465) <= 0.000001