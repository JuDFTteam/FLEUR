
import pytest
"""
Boiler-plate code for testset
"""
from read_tests import read_tests
all_tests = read_tests("noco")

@pytest.mark.fleur
@pytest.mark.noco
@pytest.mark.parametrize(("dir","desc","cmdline","mpi_procs"), all_tests)
def test_noco(dir,desc,cmdline,mpi_procs,default_fleur_test):
    """
    """

    assert default_fleur_test(dir, cmdline_args=cmdline, mpi_procs=mpi_procs)

@pytest.mark.fleur
@pytest.mark.noco
@pytest.mark.soc
@pytest.mark.bulk
@pytest.mark.hdf
def test_Fe_sc_GGA_mtNocoPot(default_fleur_test):
    """
    GGA with l_mtNocoPot=T and strong SOC: the MT spin moment has to stay along z (#805)
    """
    from xml.etree import ElementTree
    res_files = default_fleur_test("noco/Fe_sc_GGA_mtNocoPot", cmdline_args=["-warn_only"], mpi_procs=2)
    moments = ElementTree.parse(res_files['out.xml']).findall(".//magneticMomentsInMTSpheres/globalMagMoment")
    mx, my, mz = [float(x) for x in moments[-1].attrib["vec"].split()]
    assert abs(mx) < 1e-4 and abs(my) < 1e-4
    assert mz > 1.0


@pytest.mark.fleur
@pytest.mark.noco
@pytest.mark.bulk
@pytest.mark.hdf
def test_BfieldTransverse(default_fleur_test):
    """
    B-field perpendicular to the moment with l_mtNocoPot=T: moment along z with field along x and
    moment along y with field along -x are related by a spin rotation and must give the same energy
    """
    from xml.etree import ElementTree
    energies = []
    for dir in ["noco/1atBfieldTransZ", "noco/1atBfieldTransY"]:
        res_files = default_fleur_test(dir, clean=True, mpi_procs=2)
        root = ElementTree.parse(res_files['out.xml'])
        energies.append(float(root.findall(".//totalEnergy")[-1].attrib["value"]))
        moments = root.findall(".//magneticMomentsInMTSpheres/localMagMoment")
        mx, my, mz = [float(x) for x in moments[-1].attrib["vec"].split()]
        assert mx > 0.1 or my > 0.1  # moment tilted towards the field
    assert abs(energies[0] - energies[1]) < 1e-7


@pytest.mark.fleur
@pytest.mark.noco
@pytest.mark.bulk
@pytest.mark.dos
@pytest.mark.hdf
def test_FeBccNocoDM(default_fleur_test, grep_number):
    """
    Band-resolved density matrix with all four spin blocks (unsymmetrized, l_mperp):
    the traces including the spin off-diagonal part must reproduce the l-resolved DOS weights
    """
    res_files = default_fleur_test("noco/FeBccNocoDM", mpi_procs=2)
    deviation = grep_number(res_files['out'], "from the l-resolved DOS weights:")
    assert deviation < 1e-10


@pytest.mark.fleur
@pytest.mark.noco
@pytest.mark.film
@pytest.mark.soc
@pytest.mark.spinspiral
@pytest.mark.forcetheorem
def test_FePt_film_DMI(default_fleur_test):
    """
    DMI force theorem for an Fe/Pt film with q along x. H_so is odd in q and in the rotation axis and
    vanishes for axes along x and z. For q>0 and axis +y (counter-clockwise spiral) it is negative, as
    for an explicit supercell spiral with the same sense of rotation (#663)
    """
    from xml.etree import ElementTree
    res_files = default_fleur_test("noco/FePt_film_DMI", mpi_procs=2)
    h = [float(e.attrib["H_so"]) for e in ElementTree.parse(res_files['out.xml']).findall(".//allAtoms")]
    # q=+x/4 then q=-x/4, each with the axes +y, -y, +x, +z
    assert len(h) == 8
    assert -5.0e-4 < h[0] < -3.5e-4
    assert h[1] == pytest.approx(-h[0], abs=1e-8)
    assert h[4] == pytest.approx(-h[0], abs=1e-8)
    assert h[5] == pytest.approx(h[0], abs=1e-8)
    assert all(abs(h[i]) < 1e-8 for i in (2, 3, 6, 7))
