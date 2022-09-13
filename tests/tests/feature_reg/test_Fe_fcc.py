
import pytest

@pytest.mark.serial
@pytest.mark.bulk
@pytest.mark.xml
@pytest.mark.magnetism
def test_Fe_fccXML(default_fleur_test):
    """Fleur Fe fcc spin-spiral XML

    Simple test of Fleur with two steps:
    1.Generate a starting density
    2.Run 20 iterations and compare convergence, fermi-energy & total energy

    "SOC",0,"complex",1,"EVP",0
    """
    assert default_fleur_test('Fe_fccXML')


@pytest.mark.magnetism
@pytest.mark.soc
@pytest.mark.lo
@pytest.mark.hdf
def test_FeFFNLOsSOC(default_fleur_test):
    """Iron LO's and SOC test in FFN

    Simple test of Fleur with one step:
    1.Generate a starting density and run 1 iteration and check magnetization of Fe. 
    """
    assert default_fleur_test('FeFFNLOsSOC')