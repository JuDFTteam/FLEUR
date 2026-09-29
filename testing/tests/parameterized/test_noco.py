
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
@pytest.mark.dos
def test_FeBccNocoDM(default_fleur_test, grep_number):
    """
    Band-resolved density matrix with all four spin blocks (unsymmetrized, l_mperp):
    the traces including the spin off-diagonal part must reproduce the l-resolved DOS weights
    """
    res_files = default_fleur_test("noco/FeBccNocoDM", mpi_procs=2)
    deviation = grep_number(res_files['out'], "from the l-resolved DOS weights:")
    assert deviation < 1e-10
