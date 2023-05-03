
import pytest

@pytest.mark.serial
@pytest.mark.fleur
@pytest.mark.bulk
@pytest.mark.xml
@pytest.mark.dos
@pytest.mark.mcd
@pytest.mark.skip("MCD broken")
def test_CoMCDXML(default_fleur_test):
    """Fleur Co MCD XML
    Simple test of Fleur with two steps:
    1.Generate a starting density and run 2 iterations
    2.Calculate and verify MCD spectra
    """

    # Stage 1
    assert default_fleur_test("CoMCDXML/stage1")
    
    # Stage 2
    assert default_fleur_test("CoMCDXML/stage2",hdf_checks=["banddos.hdf"])

  

@pytest.mark.fleur
@pytest.mark.bulk
@pytest.mark.band
def test_CoUnfold(default_fleur_test):
    """Co band unfolding test
    Simple test of Fleur with one step:
    1.Generate a starting density and run 1 iteration with band unfolding and compare several quantities
    """
    assert default_fleur_test("CoUnfold") #No banddos.hdf as it seems unstable
