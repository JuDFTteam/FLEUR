
import pytest

@pytest.mark.serial
@pytest.mark.soc
@pytest.mark.film
@pytest.mark.xml
def test_Fe_1l_SOCXML(default_fleur_test):
    """Fleur Fe Monolayer SOC XML
    Simple test of Fleur with one steps:
    1.Generate a starting density and run 1 iteration and compare convergence, fermi-energy & total energy
    """

    assert default_fleur_test("Fe_1l_SOCXML")

@pytest.mark.serial
@pytest.mark.film
def test_Fe_1l_Tria(default_fleur_test):
    """Fleur Fe Monolayer Triangular method
    Simple test of Fleur with XML input with one step:
    1.Generate a starting density and run a single iteration and compare convergence, fermi-energy & total energy
    (with linear triangular method for fermi energy evaluation)
    """
    assert default_fleur_test("Fe_1l_Tria")
    

@pytest.mark.serial
@pytest.mark.film
@pytest.mark.xml
def test_Fe_1lXML(default_fleur_test):
    """Fleur Fe Monolayer
    Simple test of Fleur with XML input with one step:
    1.Generate a starting density and run a single iteration and compare convergence, fermi-energy & total energy
    """
    assert default_fleur_test("Fe_1lXML")
    
