import pytest

@pytest.mark.serial
@pytest.mark.fleur
@pytest.mark.bulk
@pytest.mark.noco
def test_simpleNoco_1atom(default_fleur_test):
    """Fleur Noco: simple 1 atom 
    Simple test of the noco-feature of Fleur. We take sc-Fe and do with four steps:
    1. Magnetisation in z-direction
    2. Magnetisation in x-direction
    3. Magnetisation in z-direction
    4. Magnetisation in non-sym direction
    """

    assert default_fleur_test("noco_simple1/stage1")
    assert default_fleur_test("noco_simple1/stage2",clean=True)
    assert default_fleur_test("noco_simple1/stage3",clean=True)
    assert default_fleur_test("noco_simple1/stage4",clean=True)

@pytest.mark.serial
@pytest.mark.fleur
@pytest.mark.bulk
@pytest.mark.noco
def test_simpleNoco_1atomSOC(default_fleur_test):
    """Fleur Noco: simple 1 atom with SOC 
    Simple test of the noco-feature of Fleur. We take sc-Fe and do with four steps:
    1. Magnetisation in z-direction
    2. Magnetisation in x-direction
    3. Magnetisation in z-direction
    4. Magnetisation in non-sym direction
    """

    assert default_fleur_test("noco_simple1SOC/stage1")
    assert default_fleur_test("noco_simple1SOC/stage2",clean=True)
    assert default_fleur_test("noco_simple1SOC/stage3",clean=True)
    assert default_fleur_test("noco_simple1SOC/stage4",clean=True)

  