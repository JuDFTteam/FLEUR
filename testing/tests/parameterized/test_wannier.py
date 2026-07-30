
import pytest
"""
Regression tests for the wannierlib feature (library-mode Wannier90 in FLEUR).
On top of the default out.xml comparison + schema validation these check the
Wannier spread decomposition. Systems cover the distinct FLEUR paths: no-SOC,
SOC (spinor), noco (jspins=2), AFM. WannPtSOCOps additionally covers the coarse
t_matrixelement pass (<operators_r>). Band interpolation: future round.
"""
from read_tests import read_tests
all_tests = read_tests("wannier")

# Omega = Omega_I + Omega_D + Omega_OD.
#
# Omega_I is fixed by the disentanglement (the optimal subspace) and is invariant
# under the MLWF gauge; Omega_D + Omega_OD depend on which minimum the wannierise
# iteration falls into, and that basin is decided by last-bit rounding (it moves
# with the MPI rank count, the MKL code path and the node's vector width). Asserting
# the total to 1e-5 therefore asserts on the unstable part: the previously stored
# totals could not be reproduced by any commit in their own history -- they were
# recorded in a different environment, not invalidated by a code change.
#
# So: Omega_I is the regression criterion, the total is kept only as a loose
# sanity bound. Values measured at mpi=1 with MKL_CBWR=AVX2 / I_MPI_CBWR=1,
# 2x2x2 mesh, itmax=1.
EXPECTED_OMEGA_I = {
    "WannPt":        4.841101663,  # fcc Pt, no SOC (jspins=1)
    "WannPtSOC":     9.754734063,  # fcc Pt, SOC (jspins=1, spinor)
    "WannPtSOCOps":  9.754734063,  # same system + <operators_r>; identical to WannPtSOC
                                   # to the last digit -- the operator export is gauge-neutral
    "WannFeFM":     16.713566919,  # fcc Fe FM, noco (jspins=2), no SOC
    "WannFeAFM":    16.720861134,  # fcc Fe AFM, noco (jspins=2), no SOC
    "WannFeAFMSOC": 16.694344590,  # fcc Fe AFM, noco (jspins=2) + SOC
    "WannFeBccSOC":     5.297505637,  # bcc Fe FM, COLLINEAR (jspins=2, l_noco=F) + SOC
    "WannFeAFMColSOC": 12.795211616,  # fcc Fe AFM, COLLINEAR + SOC: two sublattices, so the
                                      # spin sums cancel exactly -- the strongest check here
    "WannFeAFMSOCOps": 16.694344590,  # fcc Fe AFM, noco + SOC, now with <operators_r>: the
                                      # only coverage of the spin operator on the noco branch
}
OMEGA_I_TOL = 1.0e-5

EXPECTED_OMEGA_TOTAL = {
    "WannPt":        6.233122221,
    "WannPtSOC":    12.446292158,
    "WannPtSOCOps": 12.446292158,
    "WannFeFM":     21.934037895,
    "WannFeAFM":    21.946064187,
    "WannFeAFMSOC": 21.362891487,
    "WannFeBccSOC":     7.161120494,
    "WannFeAFMColSOC": 17.193476140,
    "WannFeAFMSOCOps": 21.362891487,
}
# Loose on purpose: absorbs a basin change, still catches a gross regression.
OMEGA_TOTAL_RTOL = 0.01

# Real-space operator files written by <operators_r>, per test id. Their presence is
# asserted by the fixture; rspauli.1 is additionally checked for physical bounds below.
_OP_R_FILES = ["WF1_hr.dat", "rspauli.1", "anglmomrs.1", "rssocmat.1", "wig_vectors"]
OPERATOR_FILES = {
    "WannPtSOCOps": _OP_R_FILES,
    "WannFeBccSOC":     _OP_R_FILES,
    "WannFeAFMColSOC": _OP_R_FILES,
    "WannFeAFMSOCOps": _OP_R_FILES,
}

# <w_0n|sigma_a|w_0n> is a Pauli expectation value on a normalized Wannier function, so
# |.| <= 1 holds elementwise -- for any gauge, which makes it basin-independent. This is
# the guard for the coarse-pass spin index: handing calc_abc the wrong (zMat, jspin)
# pairing leaves the bound violated by an order of magnitude (observed sigma_z ~ 14-19)
# while the run still completes and every other check stays green.
PAULI_BOUND = 1.0 + 1.0e-6


# Sum rule: for a non-magnetic system (jspins=1) time reversal forces Tr[sigma_a] over the
# Wannier manifold to vanish. This is what the missing spinor component broke -- it left every
# WF reporting +|up|^2 ~ +0.5, so the total came out at +9 for 18 WFs while every individual
# value stayed comfortably inside the Pauli bound. Gauge- and basin-independent.
SPIN_SUM_TOL = 0.05
# Non-magnetic, or antiferromagnetic with a sublattice-symmetric manifold: every component
# of the sum vanishes. WannFeAFMColSOC is the sharper of the two -- its WFs reach |sigma| =
# 0.993, so the cancellation to zero is not a small number made out of small numbers.
# WannFeAFMSOCOps is deliberately absent: its 36 WFs come from disentangling 72 bands, and
# that manifold does not respect the sublattice symmetry, so no exact rule applies.
NONMAGNETIC = ("WannPtSOCOps", "WannFeAFMColSOC")

# Collinear magnet quantised along z (theta=phi=0): the longitudinal sum is the manifold's
# net moment and is free, but the two transverse ones must vanish. Same kind of exact,
# gauge-independent constraint as the non-magnetic rule, and it is the one that catches a
# mix-up between the spin channels in the jspins=2 branch of the coarse pass -- where the
# radial index is isp and the radial-integral slot is 2, neither of which WannPtSOCOps
# exercises. Components are 1=sigma_x, 2=sigma_y, 3=sigma_z.
COLLINEAR_Z = ("WannFeBccSOC",)


def _rspauli_r0_diagonal_sums(path):
    """Per-component sum of Re O_nn over the R=0 diagonal of a 'generic'-format O(R) file."""
    tot = {}
    with open(path) as fh:
        for line in fh:
            f = line.split()
            if len(f) < 8:
                continue
            if f[0] == f[1] == f[2] == "0" and f[3] == f[4]:
                tot[f[5]] = tot.get(f[5], 0.0) + float(f[6])
    return tot


# anglmomrs.1 holds L(R). Two things are checked, both reference-free.
#
# L(R=0) is a matrix of <w_0i|L|w_0j> and must be hermitian, whatever the gauge.
L_HERM_TOL = 1.0e-10
# The trace over the manifold is gauge-invariant (a unitary mixing of the WFs leaves it
# alone), so it is a basin-independent quantity like Omega_I, and symmetry fixes it:
# spin-orbit coupling ties L to S but does not by itself produce a net orbital moment.
# Breaking time reversal does. A collinear magnet along z may therefore carry L_z, but
# its transverse components must vanish; an antiferromagnet with a sublattice-symmetric
# manifold must give zero in all three.
L_SUM_TOL = 1.0e-4
L_TRANSVERSE_ZERO = ("WannFeBccSOC",)
L_SUM_ZERO = ("WannFeAFMColSOC",)
# Both rules need num_wann == num_bands. A disentangled manifold does not inherit time
# reversal, so the sum over it is not the physical moment: on w222 every k is its own
# time-reversal partner, the cancellation has to happen within each k, and that needs the
# selected subspace to be T-invariant. It is not -- with Kramers doublets, keeping one
# combination of a degenerate pair leaves Omega_I unchanged, so the minimum is degenerate
# exactly at the edge of the subspace and nothing steers the iteration to the symmetric
# solution. Measured on Pt at fixed num_wann = 20, changing only whether there is anything
# to select: 20 bands for 20 WFs gives 1e-8 in all three components (and in the spin sums),
# 36 bands for 20 WFs gives (+0.029, +0.037, +0.075). The individual |L_nn| are larger in
# the first case, so this is cancellation, not small numbers.
#
# Hence WannPtSOCOps (36 -> 18) and WannFeAFMSOCOps (72 -> 36) are excluded, and the two
# that are listed have num_wann == num_bands. The same caveat applies to the spin sums
# above: WannPtSOCOps only passes NONMAGNETIC on tolerance (0.008, 0.011, 0.001), not by
# cancelling.


def _anglmom_r0(path):
    """R=0 block of a 'generic'-format O(R) file as {comp: {(i,j): complex}}."""
    blocks = {}
    with open(path) as fh:
        for line in fh:
            f = line.split()
            if len(f) < 8 or not (f[0] == f[1] == f[2] == "0"):
                continue
            blocks.setdefault(f[5], {})[(int(f[3]), int(f[4]))] = complex(
                float(f[6]), float(f[7]))
    return blocks


def _anglmom_r0_traces(path):
    """Per-component trace of the R=0 block: sum_n <w_0n|L|w_0n>."""
    return {c: sum(v.real for (i, j), v in b.items() if i == j)
            for c, b in _anglmom_r0(path).items()}


def _anglmom_r0_hermiticity(path):
    """Largest |L_ij - conj(L_ji)| over the R=0 block, worst component."""
    worst = 0.0
    for b in _anglmom_r0(path).values():
        for (i, j), v in b.items():
            w = b.get((j, i))
            if w is not None:
                worst = max(worst, abs(v - w.conjugate()))
    return worst


def _rspauli_r0_diagonal_max(path):
    """Largest |Re O_nn| over the R=0 diagonal of a 'generic'-format O(R) file.
    Layout (see m_matrixelement_io): R1 R2 R3  i j comp  Re Im."""
    worst = 0.0
    with open(path) as fh:
        for line in fh:
            f = line.split()
            if len(f) < 8:
                continue
            if f[0] == f[1] == f[2] == "0" and f[3] == f[4]:
                worst = max(worst, abs(float(f[6])))
    return worst


@pytest.mark.fleur
@pytest.mark.wannierlib
@pytest.mark.parametrize(("dir", "desc", "cmdline", "mpi_procs"), all_tests)
def test_wannier(dir, desc, cmdline, mpi_procs, default_fleur_test, grep_number):
    """Run the wannierlib test and, on top of the default out.xml checks, verify the
    gauge-invariant spread Omega_I (tight) and the total Omega (loose). Tests that
    request <operators_r> also get their O(R) exports checked."""
    test_id = dir.split("/")[-1]
    res = default_fleur_test(dir, files=OPERATOR_FILES.get(test_id),
                             cmdline_args=cmdline, mpi_procs=mpi_procs)

    if test_id in EXPECTED_OMEGA_I:
        omega_i = grep_number(res["out"], "Omega I", split="=", res_index=-1)
        assert abs(omega_i - EXPECTED_OMEGA_I[test_id]) < OMEGA_I_TOL, (
            f"gauge-invariant spread Omega_I {omega_i} deviates from reference "
            f"{EXPECTED_OMEGA_I[test_id]} (tol {OMEGA_I_TOL})")

    if test_id in EXPECTED_OMEGA_TOTAL:
        ref = EXPECTED_OMEGA_TOTAL[test_id]
        omega = grep_number(res["out"], "Omega Total", split="=", res_index=-1)
        assert abs(omega - ref) < OMEGA_TOTAL_RTOL * ref, (
            f"total spread Omega {omega} deviates from reference {ref} by more than "
            f"{100 * OMEGA_TOTAL_RTOL}% -- more than a wannierise basin change")

    if test_id in OPERATOR_FILES:
        worst = _rspauli_r0_diagonal_max(res["rspauli.1"])
        assert worst < PAULI_BOUND, (
            f"rspauli.1: max |<w_0n|sigma|w_0n>| = {worst} exceeds the Pauli bound 1; "
            "the coarse operator pass is using an inconsistent (zMat, jspin) pairing")

    if test_id in NONMAGNETIC:
        sums = _rspauli_r0_diagonal_sums(res["rspauli.1"])
        for comp, total in sorted(sums.items()):
            assert abs(total) < SPIN_SUM_TOL, (
                f"rspauli.1: sum over the R=0 diagonal of component {comp} is {total}, "
                f"but a non-magnetic system must give 0 (tol {SPIN_SUM_TOL})")

    if test_id in COLLINEAR_Z:
        sums = _rspauli_r0_diagonal_sums(res["rspauli.1"])
        for comp in ("1", "2"):
            total = sums.get(comp, 0.0)
            assert abs(total) < SPIN_SUM_TOL, (
                f"rspauli.1: transverse spin sum (component {comp}) is {total}, but a "
                f"collinear magnet along z must give 0 (tol {SPIN_SUM_TOL})")

    if test_id in OPERATOR_FILES:
        worst = _anglmom_r0_hermiticity(res["anglmomrs.1"])
        assert worst < L_HERM_TOL, (
            f"anglmomrs.1: L(R=0) is off hermitian by {worst}; <w_0i|L|w_0j> and "
            f"<w_0j|L|w_0i>* must agree (tol {L_HERM_TOL})")

    if test_id in L_SUM_ZERO:
        traces = _anglmom_r0_traces(res["anglmomrs.1"])
        for comp, total in sorted(traces.items()):
            assert abs(total) < L_SUM_TOL, (
                f"anglmomrs.1: trace of component {comp} is {total}, but this manifold "
                f"carries no net orbital moment (tol {L_SUM_TOL})")

    if test_id in L_TRANSVERSE_ZERO:
        traces = _anglmom_r0_traces(res["anglmomrs.1"])
        for comp in ("1", "2"):
            total = traces.get(comp, 0.0)
            assert abs(total) < L_SUM_TOL, (
                f"anglmomrs.1: transverse orbital moment (component {comp}) is {total}, "
                f"but a collinear magnet along z must give 0 (tol {L_SUM_TOL})")
