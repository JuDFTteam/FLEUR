
import pytest

@pytest.mark.serial
@pytest.mark.film
@pytest.mark.magnetism
@pytest.mark.spinspiral
def test_FePt_film_SSFT(default_fleur_test):
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
    assert default_fleur_test('FePt_film_SSFT/stage1')
    assert default_fleur_test('FePt_film_SSFT/stage2')


@pytest.mark.serial
@pytest.mark.film
@pytest.mark.magnetism
@pytest.mark.spinspiral
@pytest.mark.lo
def test_FePt_film_SSFT_LO(default_fleur_test):
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
    assert default_fleur_test('FePt_film_SSFT_LO/stage1')
    assert default_fleur_test('FePt_film_SSFT_LO/stage2')
