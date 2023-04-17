
import pytest

@pytest.mark.serial
def test_Fe_Tetra_noSYM(default_fleur_test):
    """Fleur Fe tetra no symmetries

    Simple test of the linear tetrahedron method with bloechl corrections:
    1. Run 15 iterations with bloechl corrections and a small k-mesh. Then
       test for fermi energy, total energy and distance
    """

    assert default_fleur_test("Fe_Tetra_noSYM")
    
