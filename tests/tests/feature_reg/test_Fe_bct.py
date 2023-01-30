
import pytest

@pytest.mark.bulk
@pytest.mark.magnetism
@pytest.mark.libxc
def test_Fe_bct_LibXC(default_fleur_test):
    """Fleur Fe bct non-collinear XML

    Simple test of Fleur with two steps:
    1.Generate a starting density
    2.Run 20 iterations and compare convergence, fermi-energy & total energy
    """
    assert default_fleur_test("Fe_bct_LibXC")

@pytest.mark.serial
@pytest.mark.bulk
@pytest.mark.magnetism
@pytest.mark.xml
def test_Fe_bct_LOXML(default_fleur_test):
    """Fleur Fe bct noco LOs XML

    Simple test of Fleur with
    1.Generate a starting density (+get some errors afterwards)
    2.Run 4 noco iterations and compare convergence, fermi-energy & total energy
    """
    assert default_fleur_test("Fe_bct_LOXML")


@pytest.mark.serial
@pytest.mark.bulk
@pytest.mark.magnetism
@pytest.mark.soc
def test_Fe_bct_SOCXML(default_fleur_test):
    """Fleur Fe bct noco SOC XML

    Simple test of Fleur with two steps:
    1.Generate a starting density
    2.Run 20 iterations and compare convergence, fermi-energy & total energy
    """
    assert default_fleur_test("Fe_bct_SOCXML")

@pytest.mark.serial
@pytest.mark.bulk
@pytest.mark.magnetism
@pytest.mark.xml
def test_Fe_bctXML(default_fleur_test):
    """Fleur Fe bct non-collinear XML

    Simple test of Fleur with two steps:
    1.Generate a starting density
    2.Run 20 iterations and compare convergence, fermi-energy & total energy
    """
    assert default_fleur_test("Fe_bctXML")
    