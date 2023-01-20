import pytest

@pytest.mark.serial
@pytest.mark.magnetism
@pytest.mark.dos
@pytest.mark.jdos
def test_SmAtomjDOS(default_fleur_test):
    """Fleur Sm Atom jDOS

    Test the decomposition into total angular momentum states in three stages:
    1. Run 15 iterations in a non-spin-polarized system with SOC and test for the fermi energy &total energy
    2. Generate a spin-polarized density with no magnetic moment using swsp (no checks here except for crash)
    3. Generate the jDOS and check the height of the peak around the fermi level
    """
    
    assert default_fleur_test("SmAtomjDOS/stage1")
    
    assert default_fleur_test("SmAtomjDOS/stage2")
    
    assert default_fleur_test("SmAtomjDOS/stage3",hdf_checks=["banddos.hdf"])

