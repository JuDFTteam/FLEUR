
import pytest

@pytest.mark.fleur
@pytest.mark.bulk
@pytest.mark.serial
@pytest.mark.xml
@pytest.mark.fast
def test_CuBulkXML(default_fleur_test):
    """
    Simple test of Fleur with XML input with one step:
    1.Generate a starting density and run 1 iteration and compare fermi-energy & total energy
    """
    assert default_fleur_test('CuBulkXML')

@pytest.mark.fleur
@pytest.mark.bulk
@pytest.mark.xml
@pytest.mark.band
@pytest.mark.fast
@pytest.mark.serial
def test_CuBandXML(default_fleur_test):
    """
    Simple test of Fleur band structure calculation with XML input with one step:
    1.Generate a starting density, run 1 iteration, and generate band structure. Ensure that the files are created.
    """
    assert default_fleur_test("CuBandXML")
    
@pytest.mark.fleur
@pytest.mark.bulk
@pytest.mark.libxc
@pytest.mark.fast
def test_CuBulkLibXC(default_fleur_test):
    """Simple test of Fleur with XML input with one step:
    1.Generate a starting density and run 1 iteration and compare fermi-energy & total energy

    If fleur is not complied and linked to libXC this test will fail
    """
    assert default_fleur_test('CuBulkLibXC')

@pytest.mark.disabled
@pytest.mark.serial
@pytest.mark.fleur
@pytest.mark.bulk
@pytest.mark.xml
@pytest.mark.dos
@pytest.mark.fast
def test_CuDOSXML(default_fleur_test):
    """Simple test of Fleur DOS calculation with XML input with one step:
    1.Generate a starting density, run 1 iteration, and generate DOS. Ensure that the files are created.

    """
    assert default_fleur_test("CuDOSXML",hdf_checks=["banddos.hdf"])
    
@pytest.mark.disabled
@pytest.mark.fleur
@pytest.mark.bulk
@pytest.mark.xml
@pytest.mark.dos
@pytest.mark.orbcomp
@pytest.mark.fast
@pytest.mark.serial
def test_CuOrb(default_fleur_test):
    """Simple test of Fleur DOS calculation with XML input with one step:
    1.Generate a starting density, run 1 iteration, and generate DOS. Ensure that the files are created.
    """
    # TODO these are the same checks as in the CuDOSXML test, is this correct?
    # we do not check orbital decom...
    assert default_fleur_test("CuOrb",hdf_checks=['banddos.hdf'])