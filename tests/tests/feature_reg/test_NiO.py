import pytest


@pytest.mark.serial
@pytest.mark.ldau
@pytest.mark.slow
@pytest.mark.magnetism
@pytest.mark.skip("AMF broken")
def test_NiOldaUAMF(default_fleur_test):
    """NiO: test with LDA+U with Around mean field double counting and magnetism

    Simple test of Fleur with one steps:
    1.Generate a starting density and run 25 iterations with LDA+U with AMF double counting
    """
    assert default_fleur_test("NiOldaUAMF")
   
@pytest.mark.serial
@pytest.mark.ldau
@pytest.mark.slow
@pytest.mark.magnetism
def test_NiOldaUFLL(default_fleur_test):
    """NiO: test with LDA+U with fully localized double counting and magnetism

    Simple test of Fleur with one steps:
    1.Generate a starting density and run 12 iterations with LDA+U with FLL double counting
    """
    
    assert default_fleur_test("NiOldaUFLL")
    