
import os
import re
import subprocess
import pytest
"""
Regression tests for the wannierlib feature (library-mode Wannier90 in FLEUR).
On top of the default out.xml comparison + schema validation these check the
Wannier spread decomposition: Omega_I exactly, Omega within a bound. Systems cover the distinct FLEUR paths: no-SOC,
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
                            onsite_weight,
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
    # The only case whose atoms are not all in groups of their own: one species, two
    # equivalent atoms. abcof expands the matching coefficients about the LOCAL frame of
    # each atom so that they can share one radial set, and a quantity that carries no atom
    # index cannot be contracted against them in that form. Every other case here has one
    # atom per group, where that frame is the global one and the distinction cannot show.
    # The two atoms are also each other's inversion partner, so this reaches the shortcut
    # that fills one of them from the other instead of expanding it.
    "WannSiDiamond": 3.622612240,
    # The only film, and the only case whose point group has a mirror (det = -1 with
    # R != -I). That is the one operation where det(R).R is neither +I nor -I, which is
    # what the axial rotation of L rests on; no bulk case here can reach it. Three atoms
    # in two groups, the two S equivalent, so it also takes the local-frame path, and the
    # vacuum is exercised for the first time. Mo d plus S p, 11 Wannier functions,
    # jspins=2 without SOC because the orbital operator has no producer otherwise.
    # Its on-site share is 97.8 %, and it does not move when the same cell is run with
    # the symmetry switched off: the rotations land where not having to rotate does.
    "WannMoS2":      16.089990243,
    # The only non-symmorphic group here: wurtzite AlN has a screw axis, tau != 0. The
    # rotation of a coefficient uses mrot and never the translation, and this is the only
    # case that can show that ignoring tau is right. Both species sit in groups of two,
    # so it also exercises firstAtom(n)+na-1 with n > 1. Sixteen valence electrons give
    # eight occupied bands, isolated by the gap, and N s+p over two atoms gives exactly
    # eight Wannier functions: no disentanglement, so the window needs no measuring.
    # Its on-site share is 98.0 %, and the same cell without symmetry agrees to 1e-9.
    "WannAlN":        5.003795238,
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
    # Same run as WannFeBcc; this case exists for the interpolation files beside it, and the
    # two report Omega_I to the last digit alike.
    "WannFeBccInterp": (2.606816455, 2.686252738),
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
# <export orbitrans="T"/> and juwip="T" do not copy the operator files: they decide where
# <operatorsRs> writes them and under which name. The map is repeated here on purpose --
# this is the only thing that can fail if the Fortran side stops moving them, because a
# file that stays behind is simply a file the fixture does not find.
#
# orbitrans opens a relative input/ and INQUIREs for generic names; JuWIP opens w90files/
# and reads ours unchanged. A second channel keeps its own name either way: WF2_hr.dat has
# no counterpart in a code that wants one 2N basis.
ORBITRANS_DIR, JUWIP_DIR = "input", "w90files"
ORBITRANS_NAME = {"rspauli.1": "spin", "anglmomrs.1": "orbital", "rssocmat.1": "spin_orbit",
                  "WF1_hr.dat": "hoppings", "wig_vectors": "rpts"}
ORBITRANS_CASES = ("WannFeBcc",)
JUWIP_CASES = ("WannPtSOCOps",)


def op(test_id, canonical):
    """Where this case leaves that operator file."""
    if test_id in ORBITRANS_CASES:
        return f"{ORBITRANS_DIR}/{ORBITRANS_NAME.get(canonical, canonical)}"
    if test_id in JUWIP_CASES:
        return f"{JUWIP_DIR}/{canonical}"
    return canonical


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
    # Hamiltonian per spin channel and the orbital operator. No spin was asked for, so
    # there is no rspauli.1: this case is here for the mirror and the vacuum, and the
    # spin sum rule is already covered by the cases above.
    "WannMoS2": ["WF1_hr.dat", "WF2_hr.dat", "anglmomrs.1", "wig_vectors"],
    "WannAlN":  ["WF1_hr.dat", "WF2_hr.dat", "anglmomrs.1", "wig_vectors"],
}

# Which case exercises each flag of <export>, and what it has to leave behind. One flag per
# case, put where the case already holds what that writer needs: the spin operator for
# blochOperators, the projection list for orbitrans, a full-zone mesh for wannierberri and
# for the symmetrizer. Before this, six of the seven writers in the block ran in no test at
# all -- they were reachable only by hand, and only the symmetrizer was ever checked.
_W90_FILES = ["WF1.amn", "WF1.mmn", "WF1.eig"]
EXPORT_FILES = {
    # The only case that asks for the Wannier90 set on its own, and the only one of these
    # with two spin channels, so the per-channel seedname is exercised at the same time.
    "WannFeAFMCol":  _W90_FILES + ["WF2.amn", "WF2.mmn", "WF2.eig"],   # wannier90
    # orbitrans writes a directory, and the fixture lists what is inside it by relative
    # path. All four are named: orbitrans reads all four, and a missing one stops it
    # long after the run that should have produced it is gone.
    "WannFeBcc": [f"{ORBITRANS_DIR}/lattice", f"{ORBITRANS_DIR}/atom_centers",
                  f"{ORBITRANS_DIR}/wann_centers", f"{ORBITRANS_DIR}/N_wann_atom"],
}

# Every artefact the block can write. A case that asked for none of them must come out
# WITHOUT them, and that is the half which makes the other half an assertion: the .amn and
# .mmn used to be written by any serial run whether or not the flag was set, so demanding
# them where the flag is on would have passed with the flag removed.
#> plotWF has no case of its own: on the default 40x40x40 grid it costs one wann_real per
#> grid point per band per k-point, which took WannAlN from 20 seconds to over nine
#> minutes -- too much for a suite that runs on every change, for a picture that feeds no
#> other code. What stays here is the half that is free: no case asks for it, so no case
#> may draw. The writer itself is checked by hand, with <export plotGrid="8 8 8"/>.
_PLOT_FILES = ("WF1_001.xsf", "WF1_001.vts", "WF2_001.xsf", "WF2_001.vts")
EXPORT_ARTEFACTS = sorted({name for names in EXPORT_FILES.values() for name in names}
                          | set(_PLOT_FILES)
                          | {f"{ORBITRANS_DIR}/{n}" for n in ORBITRANS_NAME.values()}
                          | {op(c, n) for c in JUWIP_CASES for n in _OP_R_FILES})

# The operator files in the generic O(R) format, the ones whose values can be read without
# knowing how many index columns they carry.
GENERIC_OP_FILES = ("rspauli.1", "anglmomrs.1", "rssocmat.1")

# An operator that is conjugated with respect to the gauge it is paired with is invisible to
# every other check here: conj(O) is Hermitian when O is, keeps its eigenvalues, obeys the
# same sum rules and respects the same bounds. What it does break is the decay in real
# space, because the Wannier gauge is built in one convention and the operator arrives in
# the other. Measured over the thirteen cases, the on-site share runs 73-100 % for the spin,
# 86-99 % for the orbital and 73-91 % for the spin-orbit operator; the two conjugation bugs
# this guards against left 13 % and 2 %. Fifty is the gap between the two populations, not a
# tolerance on any one of them.
ONSITE_MIN = 50.0

# Localization in R is gauge-dependent by construction: Omega_I is the invariant part of the
# spread and Omega_D + Omega_OD are not. Where the minimization leaves the gauge undetermined,
# the share at R = 0 reports which gauge the run landed in, not whether the operator decays.
#
# WannFeAFMCol is such a case. Its two spin channels reach the same Omega to nine decimals
# from gauges that leave 32-41 % at R = 0 and ones that leave 88 %, with the same L(k)
# behind both, so there the share says which gauge the run found and nothing more. Both of
# its operators move.
#
# The guard is kept rather than dropped, because the two conjugation bugs it was built for
# left 13 % and 2 %, well under the 32 % the wobble reaches here. A case that wobbles gets
# its own floor, below the wobble and above the bugs. The gauge-invariant rules further
# down -- hermiticity, the vanishing traces -- apply unchanged either way.
ONSITE_MIN_GAUGE_LABILE = {"WannFeAFMCol": 25.0}

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


#: How far the band representation matrices may sit from the identities they must satisfy.
#:
#: Measured on Si over 14 of its 16 bands and 48+48 operations: unitarity 5.3e-14, identity
#: 2.7e-14, little-group products 2.7e-14 over 3186 pairs. Six orders of room, and the errors
#: this bound exists for are nowhere near it -- reading d_band transposed, or swapping one
#: pair of bands by hand, both land at 1e-01 and worse.
DBAND_TOL = 1e-8

#: The same, for a run with spin-orbit coupling.
#:
#: Two orders looser, and the two orders are physics rather than slack: with spin-orbit on,
#: the representation gates sit at 3.1e-08 instead of the 1e-14 they reach without it, and
#: that residue is LINEAR in socscale and falls to the machine floor when it is zero. It is
#: first order in the coupling, not a broken convention -- the group, the band count and the
#: window were each ruled out by measurement. The errors this bound exists for are still far
#: above it: the wrong spin convention on a 90 degree rotation lands at 9e-01.
DBAND_TOL_SOC = 1e-6

#: Cases that write the symmetrizer export: (d_band expected, time reversal expected).
#:
#: <export symmetrizer="T"> hands a reader the k-point stars, the space group and -- when it
#: can -- d_band, the matrices that represent each operation on the bands. Nothing about it
#: reaches out.xml, and the numbers below are read from the out file, where the export writes
#: them: the venv the test runner builds carries pytest and nothing else, so a check that
#: opened the HDF file itself would skip in every environment, and a check that always skips
#: is not a check. The stars need no assertion here -- wl_build_kstar gates its own partition
#: and stops the run -- so what is left is the part a reader depends on.
#:
#: The two cases sit on opposite sides of both flags on purpose. WannAlN is wurtzite: no
#: inversion, so the eigenvectors are complex and d_band can be built, and it is the
#: non-symmorphic one (tau != 0), while jspins=2 rules time reversal out. WannSiDiamond is
#: centrosymmetric, so FLEUR solves it in real eigenvectors, and the glide of the diamond
#: structure needs a complex phase they cannot carry: d_band is refused and the stars go out
#: alone, which is the behaviour this case pins down -- an export must not kill a run over a
#: part it cannot write. Its jspins=1 is also the time-reversal side of that flag.
#: Which systems check the symmetrizer export, and what each one must come out with:
#: (band representation matrices, time reversal, a basis with local orbitals).
#:
#: The three are picked to leave no branch unvisited. Si is the refusal -- real
#: eigenvectors and a non-symmorphic operation -- and carries time reversal without
#: d_band. AlN has d_band with neither time reversal nor local orbitals, and is
#: non-symmorphic and improper, which is what pins the Wigner convention. Pt has both
#: of the things AlN has not, so the antiunitary half of the muffin-tin rule and the
#: local orbital override are only measured here.
SYM_EXPORT = {"WannSiDiamond": (False, True, False, DBAND_TOL),
              "WannAlN": (True, False, False, DBAND_TOL),
              "WannPt": (True, True, True, DBAND_TOL),
              "WannPtSOC": (True, True, True, DBAND_TOL_SOC),
              "WannFeFM": (True, False, True, DBAND_TOL),
              "WannFeAFM": (True, False, True, DBAND_TOL),
              "WannFeFMy": (True, False, True, DBAND_TOL),
              #> Spin-orbit on a magnetisation the group does not preserve: soc_sym
              #> rejects eight of the sixteen operations and mapatom leaves their atom
              #> map at zero, so the export refuses the matrices and names that reason.
              "WannFeAFMSOC": (False, False, False, DBAND_TOL_SOC),
              "WannFeAFMSOCOps": (False, False, False, DBAND_TOL_SOC)}


def _check_symmetrizer_export(out, grep_number, expect_dband, expect_tr, expect_lo,
                              tol_rep=DBAND_TOL):
    """What the symmetrizer export must have written, in the order that localises a failure."""
    txt = open(out, errors="ignore").read()
    assert "symmetry export: k-point stars" in txt, (
        "the export wrote no star report, so it did not run at all")

    has_tr = "with time reversal" in txt
    assert has_tr == expect_tr, (
        f"the export {'carries' if has_tr else 'carries no'} time reversal, the opposite of "
        f"what this system calls for. It belongs in exactly the non-magnetic cases, and the "
        f"format it is written for is the reason this export exists next to the .dmn")

    #> Not "symmetry export: d_band", which is a prefix of the line that says it was NOT
    #> written. The pair count is written once per table and nowhere else.
    refused = "d_band NOT written" in txt
    has_dband = "little-group pairs checked" in txt
    assert has_dband != refused, (
        "the out file both reports a d_band table and says it was not written; one of the "
        "two messages is being printed in the wrong branch")
    assert has_dband == expect_dband, (
        f"the export {'carries' if has_dband else 'carries no'} band representation "
        f"matrices, the opposite of what this system allows. Which runs can have them is "
        f"decided in wannierlib_main, which writes the reason to this file")
    if not expect_dband:
        assert "no atom map" in txt or "eigenvectors are real" in txt, (
            "d_band was refused without naming a reason, and the reason is the only thing "
            "that tells a refusal from a failure")
        return

    #> The muffin-tin rotation rule, against the coefficients rebuilt from the rotated
    #> eigenvector by the routine that owns that job. It is the one check here that is not
    #> an identity d_band satisfies by construction: it says the atom map, the
    #> non-symmorphic phase and the Wigner convention are right rather than consistent.
    #> The two halves are reported apart because the group table never pairs a unitary
    #> operation with an antiunitary one, so nothing else in this list sees the second.
    worst = grep_number(out, "rebuilt A and B, unitary ops", split="ops")
    assert abs(worst) < DBAND_TOL, (
        f"the muffin-tin rotation rule disagrees with the rebuilt A and B by {worst} "
        f"(tol {DBAND_TOL}); the rotation is being applied in the wrong convention")
    if expect_tr:
        worst = grep_number(out, "rebuilt A and B, antiunitary ops", split="ops")
        assert abs(worst) < DBAND_TOL, (
            f"the rule disagrees with the rebuild by {worst} on the antiunitary half "
            f"(tol {DBAND_TOL}); conjugating a state is not conjugating its coefficients, "
            f"and that half of the rule carries the m flip and the sign that go with it")

    #> And that the local orbital branch is not a no-op. Their rows of the eigenvector are
    #> zero in the rotated copy, so the rule has to supply the whole l that carries one, and
    #> this number is how far that moved the coefficients. Asserted from BELOW: a small one
    #> would mean the branch is running and changing nothing, which is the one failure a
    #> tolerance on the gates above cannot see, because they would all stay green.
    overrode = "how far the rebuild was where overridden" in txt
    assert overrode == expect_lo, (
        f"the export {'overrode' if overrode else 'did not override'} muffin-tin "
        f"coefficients for local orbitals, which is not what this basis has")
    if expect_lo:
        moved = grep_number(out, "where overridden", split="overridden")
        assert moved > 1e-6, (
            f"the local orbital override moved the coefficients by {moved}, which is "
            f"nothing: the rows it exists to replace never reached them")

    assert "NO leading band set is closed" not in txt, (
        "no leading set of bands is closed under the group, so the matrices are not a change "
        "of basis of anything and nothing downstream can symmetrise with them")
    worst = grep_number(out, "closed up to band", split="=")
    assert abs(worst) < tol_rep, (
        f"d_band is not unitary on the set it reports as closed: |D D^+ - 1| = {worst} "
        f"(tol {tol_rep})")

    err = grep_number(out, r"D\(identity\)", split="max")
    assert abs(err) < DBAND_TOL, (
        f"d_band of the identity operation is off the identity by {err} (tol {DBAND_TOL}); "
        f"the overlaps are not being taken against the states they belong to")

    pairs = grep_number(out, "little-group pairs checked", split="=")
    assert pairs > 0, (
        "no pair of operations in any little group, so the group claim measured nothing. "
        "Either the stars are wrong or every representative has a trivial little group")
    worst = grep_number(out, r"D\(i\.j\)", split="=")
    assert abs(worst) < tol_rep, (
        f"d_band does not represent the little group: D(i.j) and D(i)D(j) differ by {worst} "
        f"over {pairs:.0f} pairs (tol {tol_rep}). This is the check that sees the (m,n) "
        f"order and the non-symmorphic phase; unitarity above is blind to both")


@pytest.mark.fleur
@pytest.mark.wannierlib
@pytest.mark.parametrize(("dir", "desc", "cmdline", "mpi_procs"), all_tests)
def test_wannier(dir, desc, cmdline, mpi_procs, default_fleur_test, grep_number):
    """Run the wannierlib test and, on top of the default out.xml checks, verify the
    gauge-invariant spread Omega_I, tight, with Omega >= Omega_I beside it. Tests that
    request <operators_r> or <export> also get what those wrote checked."""
    test_id = dir.split("/")[-1]
    want_files = [op(test_id, f) for f in OPERATOR_FILES.get(test_id, ())]
    want_files += INTERP_FILES.get(test_id, [])
    want_files += EXPORT_FILES.get(test_id, [])
    if test_id in SYM_EXPORT:
        #> wannierberri="T" writes the Wannier90 set as well, and asking for it here is what
        #> checks that: WannierBerri wannierises from those three and cannot use the matrices
        #> without them. WannSiDiamond carries no wannier90="T" of its own, so if the
        #> inclusion ever stops holding this is where it shows.
        want_files += ["WF1_sym.hdf"] + _W90_FILES
    res = default_fleur_test(dir, files=want_files or None,
                             cmdline_args=cmdline, mpi_procs=mpi_procs)

    expected = set(EXPORT_FILES.get(test_id, ()))
    #> A consumer flag renames and moves what <operatorsRs> asked for, and only that. A
    #> run that asks for three operators gets three files: nothing stands in for the one
    #> it did not ask for. WannFeBcc asks for no spinOrbit, so input/spin_orbit is in the
    #> census above and in no expected set, which is what turns that into an assertion.
    expected.update(op(test_id, f) for f in OPERATOR_FILES.get(test_id, ()))
    if test_id in SYM_EXPORT:
        expected.update(_W90_FILES)
    #> How many channels a run writes is a property of the run -- a collinear jspins=2 case
    #> wannierises each channel on its own -- and not of the flag, so a case that exports the
    #> set at all may leave the second channel beside the first. Asserted by name where it is
    #> expected (WannFeAFMCol), tolerated here rather than listed per case.
    if expected & set(_W90_FILES):
        expected.update("WF2." + ext for ext in ("amn", "mmn", "eig"))
    for name in EXPORT_ARTEFACTS:
        if name not in expected:
            assert name not in res, (
                f"{name} was written although this case asks for no such export. The flag "
                f"in <export> is what decides, and a writer that runs anyway puts files no "
                f"one wanted beside the run -- the .mmn of a large mesh is hundreds of "
                f"megabytes of text, which is the whole reason these are opt-in")
    if test_id in SYM_EXPORT:
        _check_symmetrizer_export(res["out"], grep_number, *SYM_EXPORT[test_id])

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

    # Omega_total is NOT compared against the reference, at any tolerance. Which of the
    # near-degenerate MLWF minima the iteration reaches follows the order of the
    # floating-point reductions, and that order is a property of the BUILD. Measured on
    # WannFeAFMColSOC, same source and same input: the CI's gfortran build lands at
    # 17.1917 and an ifx 2026.1 build at 22.1193, 2.5 % and 25.4 % from the reference.
    # Both are clean minima -- no Wannier function above 3.47 Ang^2, against the 4 Ang^2
    # the house asks -- and Omega_I is identical to its last printed digit in both. So a
    # bound here cannot separate a bad run from a different compiler, and it reddens the
    # pipeline for whoever does not build the way the reference was made.
    #
    # What carries the weight instead: Omega_I above, tight, because it is gauge-invariant
    # and the subspace is what the physics fixes; and Omega >= Omega_I, which is an
    # identity and so holds in any basin while still catching a run that computed nothing.

    if test_id in OPERATOR_FILES:
        # Every other rule here is an upper bound -- the Pauli bound, the vanishing spin
        # sums, the vanishing orbital traces -- and a file of zeros satisfies all of them,
        # so an operator that computes nothing passes every check made on its output.
        onsite_floor = ONSITE_MIN_GAUGE_LABILE.get(test_id, ONSITE_MIN)
        for name in [f for f in OPERATOR_FILES[test_id] if f in GENERIC_OP_FILES]:
            assert onsite_weight(res[op(test_id, name)]) > onsite_floor, (
                f"{name}: only {onsite_weight(res[op(test_id, name)]):.2f} % of |O(R)|^2 sits at R = 0 "
                f"(expected above {onsite_floor}). The operator no longer decays in real "
                f"space, which is what a conjugation convention that disagrees with the "
                f"gauge looks like -- and nothing else in this file would see it.")
            assert nonzero_entries(res[op(test_id, name)]) > 0, (
                f"{name}: every entry is zero, so the operator wrote a correctly shaped "
                "file with nothing in it")

    if "rspauli.1" in OPERATOR_FILES.get(test_id, ()):
        worst = rspauli_r0_diagonal_max(res[op(test_id, "rspauli.1")])
        assert worst < PAULI_BOUND, (
            f"rspauli.1: max |<w_0n|sigma|w_0n>| = {worst} exceeds the Pauli bound 1; "
            "the coarse operator pass is using an inconsistent (zMat, jspin) pairing")

    if test_id in NONMAGNETIC:
        sums = rspauli_r0_diagonal_sums(res[op(test_id, "rspauli.1")])
        for comp, total in sorted(sums.items()):
            assert abs(total) < SPIN_SUM_TOL, (
                f"rspauli.1: sum over the R=0 diagonal of component {comp} is {total}, "
                f"but a non-magnetic system must give 0 (tol {SPIN_SUM_TOL})")

    if test_id in COLLINEAR_Z:
        sums = rspauli_r0_diagonal_sums(res[op(test_id, "rspauli.1")])
        for comp in ("1", "2"):
            total = sums.get(comp, 0.0)
            assert abs(total) < SPIN_SUM_TOL, (
                f"rspauli.1: transverse spin sum (component {comp}) is {total}, but a "
                f"collinear magnet along z must give 0 (tol {SPIN_SUM_TOL})")

    if test_id in SIGMA_Z_UNIT:
        diag = [v for (i, j), v in anglmom_r0(res[op(test_id, "rspauli.1")]).get("3", {}).items() if i == j]
        assert diag, "rspauli.1: no sigma_z entries on the R=0 diagonal"
        worst = max(abs(abs(v.real) - 1.0) for v in diag)
        assert worst < SIGMA_Z_TOL, (
            f"rspauli.1: a sigma_z diagonal entry is off +/-1 by {worst}, but each Wannier "
            f"function lies wholly in one spin channel (tol {SIGMA_Z_TOL})")

    for name in [f for f in OPERATOR_FILES.get(test_id, ()) if f.startswith("anglmomrs")]:
        worst = anglmom_r0_hermiticity(res[op(test_id, name)])
        assert worst < L_HERM_TOL, (
            f"{name}: L(R=0) is off hermitian by {worst}; <w_0i|L|w_0j> and "
            f"<w_0j|L|w_0i>* must agree (tol {L_HERM_TOL})")

    if test_id in L_SUM_ZERO:
        for name in [f for f in OPERATOR_FILES[test_id] if f.startswith("anglmomrs")]:
            for comp, total in sorted(anglmom_r0_traces(res[op(test_id, name)]).items()):
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
        traces = anglmom_r0_traces(res[op(test_id, "anglmomrs.1")])
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
