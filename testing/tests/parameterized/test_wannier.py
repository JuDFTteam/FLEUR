
import os
import re
import subprocess
import pytest
"""
Regression tests for the wannierlib feature (library-mode Wannier90 in FLEUR).
On top of the default out.xml comparison + schema validation these check the
Wannier spread decomposition. Systems cover the distinct FLEUR paths: no-SOC,
SOC (spinor), noco (jspins=2), AFM, and collinear jspins=2 without SOC -- the
only path that wannierises each spin channel separately and writes the .2
operator files. WannPtSOCOps additionally covers the coarse t_matrixelement
pass (<operators_r>). WannFeAFMCol is the only case that reaches that collinear
path with more than one atom type, which is what distinguishes an index mix-up
between spin and atom type from a layout that happens to coincide in memory.
WannFeBccInterp covers the interpolation drivers, on three output domains.
"""
from read_tests import read_tests
from wannier_files import (anglmom_r0,
                            anglmom_r0_hermiticity,
                            anglmom_r0_traces,
                            as_tuple,
                            dat_rows,
                            last_n,
                            nonzero_entries,
                            outxml_eigenvalues,
                            rspauli_r0_diagonal_max,
                            rspauli_r0_diagonal_sums,
                            spin_sumrule_values)
all_tests = read_tests("wannier")

_REFERENCE_DIR = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "..",
                              "inputfiles")

#: How far a Wannier spread may sit from the one the reference records. The run reproduces
#: the reference exactly on the toolchain the reference was made with, so this is headroom
#: for another compiler rather than a measured spread -- and still four orders below the
#: 0.97 Ang^2 by which the stored report had drifted before it was refreshed.
SPREAD_TOL = 1.0e-4


def _reference_wannier_functions(dir):
    """How many Wannier functions the stored reference reports, over all of its
    <wannierlibReport> blocks -- a collinear case writes one per spin channel."""
    with open(os.path.join(_REFERENCE_DIR, dir, "out.xml")) as fh:
        return fh.read().count("<wannierFunction ")


# ---------------------------------------------------------------- optional feature
# Wannier90 is optional. A FLEUR configured without it carries no w90_library, and every
# case here would stop with "wannierlib requires Wannier90 module API" -- a feature that
# was not asked for, not a regression, so the suite skips rather than failing. That is
# what the CI builds, and what anyone who left -wannier out of configure.sh builds.
#
# The symbol is the only witness available: the build publishes no feature list, and the
# configure decides CPP_WANNLIB_API from a compile test it does not record. If nm cannot
# be run the tests go ahead, because a loud failure is better than a silent skip.
@pytest.fixture(autouse=True, scope="module")
def _skip_without_wannier90(fleur_binary):
    binary = fleur_binary[0]
    try:
        probe = subprocess.run(["nm", binary], capture_output=True, text=True, timeout=300)
    except (OSError, ValueError, subprocess.SubprocessError):
        return
    if probe.returncode == 0 and "w90_disentangle" not in probe.stdout:
        pytest.skip("FLEUR was built without the Wannier90 module API; "
                    "configure with -wannier TRUE to run these")

# Omega = Omega_I + Omega_D + Omega_OD.
#
# Omega_I is fixed by the disentanglement (the optimal subspace) and is invariant
# under the MLWF gauge; Omega_D + Omega_OD depend on which minimum the wannierise
# iteration falls into, and that basin is decided by last-bit rounding -- it moves with
# the MPI rank count, the MKL code path and the node's vector width. So the invariant
# part is what is asserted here and the total is not. Values measured at mpi=1 with
# MKL_CBWR=AVX2 / I_MPI_CBWR=1, 2x2x2 mesh, itmax=1.
EXPECTED_OMEGA_I = {
    "WannPt":        4.841073617,  # fcc Pt, no SOC (jspins=1)
    "WannPtSOC":     9.754340102,  # fcc Pt, SOC (jspins=1, spinor)
    "WannPtSOCOps":  9.754340102,  # same system + <operators_r>; identical to WannPtSOC
                                   # to the last digit -- the operator export is gauge-neutral
    "WannFeFM":     16.711628612,  # fcc Fe FM, noco (jspins=2), no SOC
    # Same system as WannFeFM with the moment rotated to y (alpha = beta = pi/2).
    # Omega_I is identical to the last digit, as it must be: without SOC the physics is
    # isotropic and turning the moment cannot change the optimal subspace.
    "WannFeFMy":    16.711628612,
    "WannFeAFM":    16.718923683,  # fcc Fe AFM, noco (jspins=2), no SOC
    # fcc Fe AFM, noco (jspins=2) + SOC. This is the only case that is both noco and SOC, so
    # it is the only one that exercises hsmt_soc_offdiag -- the SOC block between the two
    # spin channels, which needs the full spinor structure to exist at all. The value moved
    # from 16.691730205 when that routine was fixed; every other case here is unchanged.
    "WannFeAFMSOC": 16.719550182,
    "WannFeBccSOC":     5.297166213,  # bcc Fe FM, COLLINEAR (jspins=2, l_noco=F) + SOC
    "WannFeAFMColSOC": 12.793658453,  # fcc Fe AFM, COLLINEAR + SOC: two sublattices, so the
                                      # spin sums cancel exactly -- the strongest check here
    # bcc Fe FM, COLLINEAR without SOC (jspins=2, l_noco=F, l_soc=F): the two channels are
    # separate eigenproblems, so each is wannierised on its own and there is one Omega per
    # channel. Values in channel order. This is the only test of that combination.
    "WannFeBcc": (2.606816455, 2.686252738),
    # fcc Fe AFM, COLLINEAR without SOC: same combination as WannFeBcc but with two atom
    # types, so the coefficient arrays are indexed by both spin and type and a transpose
    # between them no longer coincides in memory. The two channels are the two sublattices
    # exchanged, which is why their values agree to six digits without being identical.
    "WannFeAFMCol": (6.392856068, 6.392887652),
    "WannFeAFMSOCOps": 16.719550182,  # fcc Fe AFM, noco + SOC, now with <operators_r>: the
                                      # only coverage of the spin operator on the noco branch.
                                      # Same system as WannFeAFMSOC, so it moved with it.
}
OMEGA_I_TOL = 1.0e-5


# Real-space operator files written by <operators_r>, per test id. Their presence is
# asserted by the fixture; their contents are checked below.
_OP_R_FILES = ["WF1_hr.dat", "rspauli.1", "anglmomrs.1", "rssocmat.1", "wig_vectors"]
# The collinear no-SOC path writes one Hamiltonian per spin channel, and no spin-orbit
# operator. Spin and orbital are single files: melem_rspauli_collinear and
# melem_anglmom_collinear each assemble one 2N matrix out of both channels once they are
# wannierised. L has the same 2N shape as the spin operator -- block-diagonal, one block
# per gauge, and the cross-spin block identically zero because L acts on the spatial part
# alone.
_OP_R_FILES_2CH = ["WF1_hr.dat", "WF2_hr.dat", "anglmomrs.1",
                   "rspauli.1", "wig_vectors"]
OPERATOR_FILES = {
    "WannPtSOCOps": _OP_R_FILES,
    "WannFeBccSOC":     _OP_R_FILES,
    "WannFeAFMColSOC": _OP_R_FILES,
    "WannFeAFMSOCOps": _OP_R_FILES,
    # Only <operator name="spin"/>: this case exists for the spin sum rule below, and
    # asking for the rest would cost time without adding coverage.
    "WannFeFMy": ["rspauli.1", "wig_vectors"],
    "WannFeBcc":        _OP_R_FILES_2CH,
    "WannFeAFMCol":     _OP_R_FILES_2CH,
}

# The operator files in the generic O(R) format, the ones whose values can be read without
# knowing how many index columns they carry.
GENERIC_OP_FILES = ("rspauli.1", "anglmomrs.1", "rssocmat.1")

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
# Without spin-orbit coupling H commutes with sigma.n, so every Bloch state is a pure
# spinor and |<s>| = 1 exactly, whatever the direction of the moment. That is physics, not a
# measured value, so this needs no reference number and never goes stale. FLEUR already
# prints the per-band expectation values; this only reads them.
#
# WannFeFMy is the one case with the moment off the xz-plane (alpha = pi/2). It is here
# because the azimuthal rotation of the spin operator carried the wrong sign and no test
# noticed: every other noco test has alpha = 0, where the broken term does not contribute.
# With that sign wrong the same run reported |<s>| = 0.0313.
PURE_SPINOR = ("WannFeFMy",)
PURE_SPINOR_TOL = 1.0e-3


NONMAGNETIC = ("WannPtSOCOps", "WannFeAFMColSOC", "WannFeAFMCol")

# Collinear magnet quantised along z (theta=phi=0): the longitudinal sum is the manifold's
# net moment and is free, but the two transverse ones must vanish. Same kind of exact,
# gauge-independent constraint as the non-magnetic rule, and it is the one that catches a
# mix-up between the spin channels in the jspins=2 branch of the coarse pass -- where the
# radial index is isp and the radial-integral slot is 2, neither of which WannPtSOCOps
# exercises. Components are 1=sigma_x, 2=sigma_y, 3=sigma_z.
COLLINEAR_Z = ("WannFeBccSOC", "WannFeBcc")

# Without spin-orbit coupling the 2N Pauli is assembled from two separately wannierised
# channels, so sigma_z is block-diagonal and orthonormality within each channel fixes its
# R=0 diagonal at exactly +/-1 -- nw of each. That is sharper than the Pauli bound: it pins
# the value instead of bounding it, and it is what a broken gauge rotation of the cross-spin
# overlap would fail. Only sigma_z: the transverse components live entirely in the
# off-diagonal blocks, so their diagonal is zero by construction.
# --- Wannier interpolation ----------------------------------------------------------
#
# The assertion that carries it is reference-free and basin-independent: Wannier
# interpolation is EXACT on the mesh it was built from. H_W(k) = V^dagger diag(eig) V is a
# unitary rotation of the input eigenvalues, so its spectrum IS the input spectrum, and the
# k -> R -> k round trip is the identity on the original mesh. Interpolating onto w222 --
# the wannierisation mesh itself -- must therefore reproduce the eigenvalues that the same
# run wrote into out.xml, whatever basin the wannierisation fell into.
#
# It holds only without disentanglement: with num_bands > num_wann the Hamiltonian is built
# from the projected eigval2 rather than from eig. WannFeBccInterp has 6 bands for 6 Wannier
# functions, which is why it is the case that carries this.
INTERP_EXACT = {"WannFeBccInterp": (5, 10)}          # test id -> (minBand, maxBand)
# Worst measured deviation 4.3e-9 over 8 k-points and both spin channels -- half of the last
# digit of the f14.8 output. The residual is print rounding, not arithmetic.
INTERP_EXACT_TOL = 1.0e-8

# One file per interpolation driver, for every declared output domain and every wannierised
# spin channel. Absence is the check: a driver that silently stops writing is otherwise
# invisible, since every value test below would pass on a file that is not there.
INTERP_FILES = {
    "WannFeBccInterp": [
        "bands_wann_%s%s_spin%d.dat" % (base, dom, ch)
        for base in ("interpol", "interpol_ev", "orbmom", "velocity", "berrycurv",
                     "eigenstates")
        for dom in ("", "_plane", "_grid")
        for ch in (1, 2)
    ] + [# The only A(R) written: the WYSV A^(W)(R) that feeds the Berry curvature and the
         # anomalous velocity. Presence only, no value is checked here.
         "WF1_r.dat", "berry_centre_check.dat",
         # B(R) = <0n|H r|Rm>. NO numerical reference exists: nothing here says the
         # values are right, only that the operator still runs and still writes
         # something. What would pin it down is that B is linear in the eigenvalues,
         # so shifting the spectrum by a constant must move B by that constant times
         # A(R) -- and A is anchored. Until that test exists, this is presence only.
         "WF1_bmn.dat", "WF2_bmn.dat",
         # C(R) and F(R) of the modern theory. Presence and non-zero only, for the
         # same reason as B: no reference exists. What this does buy is that the
         # code paths behind them run at all -- until now nothing executed them.
         "WF1_cmn.dat", "WF2_cmn.dat", "WF1_fmn.dat", "WF2_fmn.dat"],
}

# The velocity must be checked on the FINE path and nowhere else: every point of w222 is a
# high-symmetry point, where dE/dk vanishes by symmetry, so a file of zeros is the correct
# answer there. The _plane domain is the 240-point path-2 list.
VELOCITY_FINE = {"WannFeBccInterp": "bands_wann_velocity_plane_spin1.dat"}

SIGMA_Z_UNIT = ("WannFeBcc", "WannFeAFMCol")
SIGMA_Z_TOL = 1.0e-8


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
# WannFeBcc has no spin-orbit coupling at all, so nothing ties L to the lattice and all
# three traces vanish -- for a sharper reason than the antiferromagnet's cancellation: the
# Bloch states are real, and L is imaginary in a real basis, so the gauge-invariant trace is
# exactly zero even though the individual |<w_n|L|w_n>| are not (the Wannier gauge is
# complex). num_wann == num_bands == 6, so the caveat below does not apply.
L_SUM_ZERO = ("WannFeAFMColSOC", "WannFeBcc", "WannFeAFMCol")
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


@pytest.mark.fleur
@pytest.mark.wannierlib
@pytest.mark.parametrize(("dir", "desc", "cmdline", "mpi_procs"), all_tests)
def test_wannier(dir, desc, cmdline, mpi_procs, default_fleur_test, grep_number):
    """Run the wannierlib test and, on top of the default out.xml checks, verify the
    gauge-invariant spread Omega_I (tight) and the total Omega (bounded, not reproduced). Tests that
    request <operators_r> also get their O(R) exports checked."""
    test_id = dir.split("/")[-1]
    want_files = list(OPERATOR_FILES.get(test_id, ()))
    want_files += INTERP_FILES.get(test_id, [])
    spreads = [["wannierFunction", "spread", i, SPREAD_TOL, None]
               for i in range(_reference_wannier_functions(dir))]
    res = default_fleur_test(dir, files=want_files or None, checks=spreads,
                             cmdline_args=cmdline, mpi_procs=mpi_procs)

    omega_i_got = ()
    if test_id in EXPECTED_OMEGA_I:
        # A tuple means the run wannierises more than once -- one collinear spin channel
        # after the other -- so every value is checked, in the order they are written.
        refs = as_tuple(EXPECTED_OMEGA_I[test_id])
        omega_i_got = last_n(grep_number(res["out"], "Omega I", split="=", res_index=None),
                              len(refs))
        for ch, (ref, omega_i) in enumerate(zip(refs, omega_i_got), start=1):
            assert abs(omega_i - ref) < OMEGA_I_TOL, (
                f"gauge-invariant spread Omega_I {omega_i} of channel {ch} deviates from "
                f"reference {ref} (tol {OMEGA_I_TOL})")

    # Omega = Omega_I + Omega_OD with Omega_OD >= 0. An identity, so it holds whatever
    # basin the run found -- and it is what stops a run that computed nothing from passing
    # a comparison against spreads it produced itself.
    if omega_i_got:
        got = last_n(grep_number(res["out"], "Omega Total", split="=", res_index=None),
                      len(omega_i_got))
        for ch, (omega_i, omega) in enumerate(zip(omega_i_got, got), start=1):
            assert omega >= omega_i - OMEGA_I_TOL, (
                f"total spread Omega {omega} of channel {ch} is below its own invariant part "
                f"Omega_I {omega_i}, which is impossible: Omega = Omega_I + Omega_OD")

    if test_id in OPERATOR_FILES:
        # Every other rule here is an upper bound -- the Pauli bound, the vanishing spin
        # sums, the vanishing orbital traces -- and a file of zeros satisfies all of them,
        # so an operator that computes nothing passes every check made on its output.
        for name in [f for f in OPERATOR_FILES[test_id] if f in GENERIC_OP_FILES]:
            assert nonzero_entries(res[name]) > 0, (
                f"{name}: every entry is zero, so the operator wrote a correctly shaped "
                "file with nothing in it")

    if "rspauli.1" in OPERATOR_FILES.get(test_id, ()):
        worst = rspauli_r0_diagonal_max(res["rspauli.1"])
        assert worst < PAULI_BOUND, (
            f"rspauli.1: max |<w_0n|sigma|w_0n>| = {worst} exceeds the Pauli bound 1; "
            "the coarse operator pass is using an inconsistent (zMat, jspin) pairing")

    if test_id in NONMAGNETIC:
        sums = rspauli_r0_diagonal_sums(res["rspauli.1"])
        for comp, total in sorted(sums.items()):
            assert abs(total) < SPIN_SUM_TOL, (
                f"rspauli.1: sum over the R=0 diagonal of component {comp} is {total}, "
                f"but a non-magnetic system must give 0 (tol {SPIN_SUM_TOL})")

    if test_id in COLLINEAR_Z:
        sums = rspauli_r0_diagonal_sums(res["rspauli.1"])
        for comp in ("1", "2"):
            total = sums.get(comp, 0.0)
            assert abs(total) < SPIN_SUM_TOL, (
                f"rspauli.1: transverse spin sum (component {comp}) is {total}, but a "
                f"collinear magnet along z must give 0 (tol {SPIN_SUM_TOL})")

    if test_id in SIGMA_Z_UNIT:
        diag = [v for (i, j), v in anglmom_r0(res["rspauli.1"]).get("3", {}).items() if i == j]
        assert diag, "rspauli.1: no sigma_z entries on the R=0 diagonal"
        worst = max(abs(abs(v.real) - 1.0) for v in diag)
        assert worst < SIGMA_Z_TOL, (
            f"rspauli.1: a sigma_z diagonal entry is off +/-1 by {worst}, but each Wannier "
            f"function lies wholly in one spin channel (tol {SIGMA_Z_TOL})")

    for name in [f for f in OPERATOR_FILES.get(test_id, ()) if f.startswith("anglmomrs")]:
        worst = anglmom_r0_hermiticity(res[name])
        assert worst < L_HERM_TOL, (
            f"{name}: L(R=0) is off hermitian by {worst}; <w_0i|L|w_0j> and "
            f"<w_0j|L|w_0i>* must agree (tol {L_HERM_TOL})")

    if test_id in L_SUM_ZERO:
        for name in [f for f in OPERATOR_FILES[test_id] if f.startswith("anglmomrs")]:
            for comp, total in sorted(anglmom_r0_traces(res[name]).items()):
                assert abs(total) < L_SUM_TOL, (
                    f"{name}: trace of component {comp} is {total}, but this manifold "
                    f"carries no net orbital moment (tol {L_SUM_TOL})")

    if test_id in INTERP_EXACT:
        # The eigenvalues the same run wrote, as the reference for its own interpolation.
        lo, hi = INTERP_EXACT[test_id]
        ref = outxml_eigenvalues(res["out.xml"])
        assert ref, "out.xml carries no <eigenvaluesAt> blocks to check the interpolation against"
        for ch in (1, 2):
            rows = dat_rows(res["bands_wann_interpol_spin%d.dat" % ch])
            assert rows, "bands_wann_interpol_spin%d.dat has no data rows" % ch
            for ik, row in enumerate(rows, start=1):
                want = sorted(ref[(ch, ik)][lo - 1:hi])
                have = sorted(row[1:])          # drop kdist
                assert len(want) == len(have), (
                    f"channel {ch}, k={ik}: {len(have)} interpolated bands against "
                    f"{len(want)} in the window {lo}..{hi}")
                worst = max(abs(a - b) for a, b in zip(want, have))
                assert worst < INTERP_EXACT_TOL, (
                    f"channel {ch}, k={ik}: the interpolated bands differ from the input "
                    f"eigenvalues by {worst}, but Wannier interpolation is exact on the "
                    f"mesh it was built from (tol {INTERP_EXACT_TOL})")

    if test_id in VELOCITY_FINE:
        # Shape, and that it is not a file of zeros. dE/dk is only forced to vanish on the
        # high-symmetry mesh, so on the fine path an all-zero velocity means the driver
        # produced nothing -- which every other check here would accept.
        rows = dat_rows(res[VELOCITY_FINE[test_id]])
        assert rows, "the velocity file on the fine path has no data rows"
        nw = (len(rows[0]) - 1) // 4          # per band: E, vx, vy, vz
        assert nw >= 1 and len(rows[0]) == 1 + 4 * nw, (
            f"velocity row has {len(rows[0])} columns, not 1 + 4*num_wann")
        worst = max(abs(r[1 + 4 * b + 1 + a]) for r in rows for b in range(nw)
                    for a in range(3))
        assert worst > 0.0, (
            "every velocity component on the fine path is exactly zero, so the driver "
            "wrote a correctly shaped file with nothing in it")

    if test_id in L_TRANSVERSE_ZERO:
        traces = anglmom_r0_traces(res["anglmomrs.1"])
        for comp in ("1", "2"):
            total = traces.get(comp, 0.0)
            assert abs(total) < L_SUM_TOL, (
                f"anglmomrs.1: transverse orbital moment (component {comp}) is {total}, "
                f"but a collinear magnet along z must give 0 (tol {L_SUM_TOL})")

    # Pure spinors: see PURE_SPINOR above. Reads what FLEUR already printed.
    if test_id in PURE_SPINOR:
        vals = spin_sumrule_values(res["out"])
        assert vals, ("no spin sum-rule block in the output -- this test needs "
                      "<operators_r> with the spin operator")
        worst = min(vals)
        assert abs(worst - 1.0) < PURE_SPINOR_TOL, (
            f"without SOC every Bloch state is a pure spinor, so |<s>| must be 1; the "
            f"smallest of {len(vals)} values is {worst:.6f}. The azimuthal rotation of the "
            f"spin operator is the usual cause -- it only shows up when alpha != 0.")
