
import pytest
"""
Boiler-plate code for testset
"""
from read_tests import read_tests
all_tests = read_tests("forces")

@pytest.mark.fleur
@pytest.mark.forces
@pytest.mark.parametrize(("dir","desc","cmdline","mpi_procs"), all_tests)
def test_forces(dir,desc,cmdline,mpi_procs,default_fleur_test):
    """
    """

    assert default_fleur_test(dir, cmdline_args=cmdline, mpi_procs=mpi_procs)

@pytest.mark.fleur
@pytest.mark.forces
@pytest.mark.bulk
def test_H2ORelaxCG(default_fleur_test):
    """
    CG relaxation step starting from one history step in relax.xml (#309)
    """
    import math
    import re

    res_files = default_fleur_test("forces/H2ORelaxCG", files=["relax.xml"], mpi_procs=2)
    with open(res_files["relax.xml"]) as f:
        disp = [[float(x) for x in d.split()] for d in re.findall(r"<displace>(.*?)</displace>", f.read())]

    expected = [[0.0, -0.2552066403, 0.0], [0.0002905528, 0.1303774103, 0.0]]
    assert len(disp) == len(expected)
    for d, e in zip(disp, expected):
        assert all(math.isfinite(x) for x in d)
        assert d == pytest.approx(e, abs=1e-4)