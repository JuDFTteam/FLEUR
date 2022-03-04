
import pytest

@pytest.mark.serial
@pytest.mark.film
@pytest.mark.magnetism
@pytest.mark.spinspiral
def test_FePt_film_SSFT(execute_fleur, grep_number, grep_exists, validate_out_xml_file):
    """Fe monolayer fcc (110): spin spiral energy

    Simple test of Fleur with two steps:
    1.Generate the reference charge density
    2.Calculatie spin spiral energy via the force theorem
    This test was callibrated via MaX-R4 release becuase
    MaX-R4 was shown to give the same results as version 26d.
    Uses: spin spiral force theorem
    WARNING: the test is aimed to calculate	energy difference
    between	two spin spirals. If it	starts to fail,	please 
    check if energy	difference between q1 and q2 remains the
    same. If it is true, just adjust q1 and	q2 values in test.run2.
    Otherwise there	seems to be a bug.
    """
    test_file_folder = './inputfiles/FePt_film_SSFT/'

    # Stage 1
    res_files = execute_fleur(test_file_folder, only_copy=['inp.xml'])
    should_files = ['out']
    res_file_names = list(res_files.keys())
    for file1 in should_files:
        assert file1 in res_file_names
    
    tenergy = grep_number(res_files['out'], "total energy=", "=")
    assert abs(tenergy - -19707.2347264178) <= 0.0001

    # Stage 2
    files = [['inp2.xml', 'inp.xml'], 'JUDFT_WARN_ONLY']
    res_files = execute_fleur(test_file_folder, only_copy=files, rm_files=['mixing_history'])
    should_files = ['out.xml']
    res_file_names = list(res_files.keys())
    for file1 in should_files:
        assert file1 in res_file_names
    
    assert validate_out_xml_file(res_files['out.xml'])

    assert grep_exists(res_files['out.xml'], "Forcetheorem_SSDISP qvectors=")
    ev_q1 = grep_number(res_files['out.xml'], 'Entry q="1"', "ev-sum=")
    ev_q2 = grep_number(res_files['out.xml'], 'Entry q="2"', "ev-sum=")

    assert abs(ev_q1 - -5.1963027) <= 0.000001
    assert abs(ev_q2 - -5.1706272) <= 0.000001


@pytest.mark.serial
@pytest.mark.film
@pytest.mark.magnetism
@pytest.mark.spinspiral
@pytest.mark.lo
def test_FePt_film_SSFT_LO(execute_fleur, grep_number, grep_exists, validate_out_xml_file):
    """Fe monolayer fcc (110): spin spiral energy with LO

    Simple test of Fleur with two steps:
    1.Generate the reference charge density
    2.Calculatie spin spiral energy via the force theorem
    This test was callibrated via MaX-R4 release becuase
    MaX-R4 was shown to give the same results as version 26d.
    Uses: spin spiral force theorem
    WARNING: the test is aimed to calculate	energy difference
    between	two spin spirals. If it	starts to fail,	please 
    check if energy	difference between q1 and q2 remains the
    same. If it is true, just adjust q1 and	q2 values in test.run2.
    Otherwise there	seems to be a bug.
    """
    test_file_folder = './inputfiles/FePt_film_SSFT_LO/'

    # Stage 1
    res_files = execute_fleur(test_file_folder, only_copy=['inp.xml'])
    should_files = ['out']
    res_file_names = list(res_files.keys())
    for file1 in should_files:
        assert file1 in res_file_names
    
    tenergy = grep_number(res_files['out'], "total energy=", "=")
    assert abs(tenergy - -19706.9922262005) <= 0.0001

    # Stage 2
    files = [['inp2.xml', 'inp.xml'], 'JUDFT_WARN_ONLY']
    res_files = execute_fleur(test_file_folder, only_copy=files, rm_files=['mixing_history'])
    should_files = ['out.xml']
    res_file_names = list(res_files.keys())
    for file1 in should_files:
        assert file1 in res_file_names
    
    assert validate_out_xml_file(res_files['out.xml'])

    assert grep_exists(res_files['out.xml'], "Forcetheorem_SSDISP qvectors=")
    ev_q1 = grep_number(res_files['out.xml'], 'Entry q="1"', "ev-sum=")
    ev_q2 = grep_number(res_files['out.xml'], 'Entry q="2"', "ev-sum=")

    assert abs(ev_q1 - -37.3675172) <= 0.000001
    assert abs(ev_q2 - -37.3421763) <= 0.000001
