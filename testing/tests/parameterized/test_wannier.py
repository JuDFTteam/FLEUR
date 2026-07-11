
import pytest
"""
Regression tests for the wannierlib feature (library-mode Wannier90 in FLEUR).
This round checks the total Wannier spread Omega (reproducibility, bit-for-bit at mpi=1)
in addition to the default out.xml comparison + schema validation.
Systems cover the distinct FLEUR paths: no-SOC, SOC (spinor), noco (jspins=2), AFM.
Operators / band interpolation: future round.
"""
from read_tests import read_tests
all_tests = read_tests("wannier")

# Reference total Wannier spread Omega (last "Sum of centres and spreads" line in 'out'),
# per test id. Deterministic at mpi=1 (I_MPI_CBWR=1). Values on 2x2x2 mesh, itmax=1.
EXPECTED_SPREAD = {
    "WannPt":       6.23779531,   # fcc Pt, no SOC (jspins=1)
    "WannPtSOC":    12.48703050,  # fcc Pt, SOC (jspins=1, spinor)
    "WannFeFM":     21.93403790,  # fcc Fe FM, noco (jspins=2), no SOC
    "WannFeAFM":    21.94606419,  # fcc Fe AFM, noco (jspins=2), no SOC
    "WannFeAFMSOC": 21.35218438,  # fcc Fe AFM, noco (jspins=2) + SOC
}
SPREAD_TOL = 1.0e-5


@pytest.mark.fleur
@pytest.mark.wannierlib
@pytest.mark.parametrize(("dir", "desc", "cmdline", "mpi_procs"), all_tests)
def test_wannier(dir, desc, cmdline, mpi_procs, default_fleur_test, grep_number):
    """Run the wannierlib test and, on top of the default out.xml checks,
    verify the final total spread Omega against the stored reference."""
    res = default_fleur_test(dir, cmdline_args=cmdline, mpi_procs=mpi_procs)

    test_id = dir.split("/")[-1]
    if test_id in EXPECTED_SPREAD:
        omega = grep_number(res["out"], "Sum of centres and spreads",
                            split=")", line_index=1, res_index=-1)
        assert abs(omega - EXPECTED_SPREAD[test_id]) < SPREAD_TOL, (
            f"Wannier spread {omega} deviates from reference "
            f"{EXPECTED_SPREAD[test_id]} (tol {SPREAD_TOL})")
