<!-- Copyright (c) 2026 Mike Towler. MIT licence; see LICENSE in this directory. -->
<!-- Copied from the tree wien2kqmc, the reference writer: "this tree" below means wien2kqmc.
     FLEUR writes this format with output/lapwfn (src/fleur/io/writeLapwfn.F90; see README.md). -->

# SPEC: lapwfn.data, a code-neutral file of (L)APW orbitals for CASINO

Format version 1, 20 October 2026. The file carries what
CASINO's LAPW orbital type needs and nothing a code could not supply from
its own eigenstates: a writer inside elk or exciting, using its own radial
functions, local orbitals and matching coefficients, produces the same
file. The reference writer is wien2kqmc (this tree), from WIEN2k 24.1.
Clarified (20 October 2026) after an independent check, the layout
unchanged: what version 1 carries (LAPW only), the kinds, the join
radius, the core normalisation, the Fermi energy and functional, and
what gen_blip says of pwfn.data. Amended 21 October 2026, the layout
unchanged and every version-1 file still valid: a degenerate level at
k = G/2 may be written as real combinations of its members, as the Elk
and exciting writers write it (section 1, "Real orbitals").

## 1. The representation

Hartree atomic units throughout: lengths in bohr, energies in hartree.
Vectors are Cartesian in the crystal's frame; there are no local frames.

A crystal of N atoms per primitive cell, each of a type a with a sphere of
radius R_a about its centre R_alpha. An orbital n at Bloch vector k is

    outside every sphere     psi_nk(r) = sum_G c_G exp(i (k + G).r)
    inside sphere alpha      psi_nk(r) = exp(i k.T) sum_{l<=L} sum_m sum_j
                                          c^{alpha,j}_{lm} u_{a,l,j}(|p|) Y_lm(p^),
                                          p = r - R_alpha - T,  |p| < R_a,

where T is the lattice vector that brings the sphere's image containing r
onto the sphere listed, G runs over the reciprocal lattice vectors of the
primitive cell, u_{a,l,j} (j = 1..n_{a,l}) are the radial functions of type
a and angular momentum l, tabulated on the type's mesh, and Y_lm are the
complex spherical harmonics with the Condon-Shortley phase, in the order
(0,0), (1,-1), (1,0), (1,1), (2,-2), ..., (L,L):

    Y_lm = (-1)^m sqrt((2l+1)/(4 pi) (l-m)!/(l+m)!) P_l^m(cos theta) e^{i m phi}  (m >= 0),
    Y_{l,-m} = (-1)^m Y_lm^*.

The two parts are scaled together: with both the orbital is normalised to
one over the primitive cell when the source code's eigenvector is. The
plane-wave part is NOT in this file: it is pwfn.data in the same directory,
in CASINO's own format, with the same k points and bands in the same order,
its coefficients being c_G above with no prefactor (CASINO's gen_blip turns
it into blips, which the orbital type reads for the interstitial). The
sphere expansion and the plane-wave sum match in value at the sphere
surface to the accuracy of the matching (LAPW: value and slope per channel;
the plane waves beyond L are the residual); the file carries no blend. A
reader that joins the two smoothly does so itself (DESIGN_wien2k.md).

VERSION 1 CARRIES LAPW ORBITALS ONLY: every channel's orbital is
continuous in value and slope at the sphere surface. A source whose
orbitals are APW+lo in any channel (continuous in value only, with a
kink at the surface) is not written in version 1; the writer refuses it.

pwfn.data under gen_blip. Because the plane waves are the orbital only
outside the spheres, gen_blip warns, for every file of this kind, that
"Your plane-wave orbitals are not normalized over the primitive cell"
and that they "are nonorthogonal", and the kinetic energies it prints
are those of the interstitial expansion alone. All three are expected
and harmless: the blips reproduce the plane-wave sum as written (alpha
0.99999991 to 1 for the reference writer's diamond and BN).

Core states are atomic: for each type a and core level c, one radial
function u_{a,c}(r) with its l, tabulated on an extended mesh that reaches
beyond R_a to where the function is negligible, so that the core orbitals
of atom alpha are u_{a,c}(|r - R_alpha|) Y_lm(r - R_alpha^), m = -l..l,
with no Bloch sum written (a reader sums images if it wants them; they are
below 1e-12 of the maximum by the end of the mesh).

The inner series. Every radial function, valence and core, is regular at
the nucleus: inside its join radius r_j its table IS the series

    u(r) = A r^l (1 - Z r/(l+1) + c2 r^2 + c3 r^3),   r <= r_j,

joined to the source code's function at r_j: the Kato condition fixes
the r^{l+1} term; for s and p channels A, c2 and c3 match the value, slope
and curvature at r_j (c2 is then close to the regular solution's
(Z^2/(l+1) + V_0 - E)/(2l+3)); for l >= 2, where the two terms are small
corrections to r^l, c3 = 0 and A, c2 match value and slope, the curvature
jump being below 2e-3 of itself. The file gives (r_j, A, c2, c3) so that a
reader may evaluate the series instead of interpolating the table there.
The join radius is the writer's choice, at most 0.05 bohr and inside the
first tenth of the mesh's radius; the reference writer's rule (the
docstring of settle_index in wien2kqmc/convert.py) puts the valence
functions' joins of the five reference decks at 1.5e-4 to 2.0e-2 bohr,
their p channels at 4.4e-3 to 2.0e-2 (each deck's MEASUREMENTS.txt), and
refuses a join beyond the bound. Readers need only (r_j, A, c2, c3). (WIEN2k's own tables start differently over
their innermost 150 or so points, HOW_WIEN2K_WORKS.md section 3.)

The core's tail. Beyond the sphere radius a core function is
u(R_a) exp(-kappa (r - R_a)), kappa = -u'(R_a)/u(R_a), continuous in value
and slope; the file gives kappa and tabulates the tail on the extended
mesh. A core function is normalised to one over all space, its tail
included: the integral of u(r)^2 r^2 dr from 0 to the end of its mesh is
1 (the reference writer to 5e-7; boron's 1s at a sphere radius of 1.35
bohr to 8e-6). Core states carry no coefficients, so unlike the valence
functions their normalisation is stated.

Real orbitals. At a k point equivalent to -k (k = G/2) the writer
multiplies each orbital by a constant phase that makes it real wherever
one exists (a non-degenerate level), and says so per band in the
coefficient block; a degenerate level that no phase makes real is written
either complex, its real and imaginary parts then each a Bloch function at
k, or as real combinations of its own members orthonormal in the basis's
inner product, each band's eigenvalue then the expectation value of its
orbital; either leaves the level's span, and so the density of a wholly
occupied level, as the source gives it, and a writer may leave a level
that the Fermi energy divides as the source gives it. The pwfn.data
coefficients carry the same phase or combination, written as complex
pairs always.

## 2. The layout

Plain text, one item per line after its comment line, in the style of
pwfn.data: section titles underlined with dashes, numbers free-format,
complex numbers as two reals (real part, imaginary part). Line one is

    LAPWFN.DATA: LAPW orbitals, format version 1

and a reader checks the version. Then:

    Title
    <one line>

    BASIC INFO
    ----------
    Generated by:
    <writer, source code and version, case, date>
    Method:
    DFT
    DFT functional:
    <name>
    Relativistic treatment:
    none
    Spin polarized:
    F
    Units:
    hartree, bohr
    Total energy (au per primitive cell)
    <E>
    Fermi energy (au)
    <E_F>                           (the converged SCF's, see below)
    Number of electrons per primitive cell (core included)
    <N_e>
    Number of valence electrons per primitive cell (in the WAVE FUNCTION section)
    <N_v>
    Interstitial part:
    pwfn.data: same k points and bands in the same order; coefficients multiply exp(i(k+G).r) with no prefactor

    GEOMETRY
    --------
    Number of atoms per primitive cell
    <N>
    Atomic number, type, position (au) of each atom
    <Z> <type> <x> <y> <z>          (N lines; the position is the expansion's centre)
    Primitive lattice vectors (au)
    <a1x> <a1y> <a1z>
    <a2x> <a2y> <a2z>
    <a3x> <a3y> <a3z>

    ATOM TYPES
    ----------
    Number of atom types
    <n_t>
    [for each type:]
    Type, label, atomic number, sphere radius (au), number of mesh points
    <t> <label> <Z> <R> <n_r>
    Mesh points (au)
    <n_r values, four per line>
    Maximum l of the sphere expansion
    <L>
    Number of radial functions for l = 0 to L
    <L+1 integers>
    [for each l = 0..L, for each function j = 1..n_l:]
    Radial function: l, j, kind, energy (au)
    <l> <j> <u|udot|lo> <E>
    Inner series: join radius, A, c2, c3
    <r_j> <A> <c2> <c3>
    Values u(r) on the mesh
    <n_r values, four per line>

    CORE STATES
    -----------
    Number of core states of each type
    <n_t integers>
    [for each type, for each core state:]
    Core state: type, label, n, l, energy (au), occupation, number of mesh points
    <t> <label> <n> <l> <E> <occ> <n_c>
    Inner series: join radius, A, c2, c3
    <r_j> <A> <c2> <c3>
    Tail beyond the sphere: kappa
    <kappa>
    Mesh points (au)
    <n_c values, four per line>
    Values u(r) on the mesh
    <n_c values, four per line>

    WAVE FUNCTION
    -------------
    Number of k-points
    <n_k>
    [for each k point:]
    k-point # ; # of bands (up spin/down spin) ; k-point coords (au)
    <ik> <n_up> <n_down> <kx> <ky> <kz>
    [for each band:]
    Band, spin, eigenvalue (au), real (T/F)
    <n> <s> <E> <T|F>
    [for each atom alpha = 1..N:]
    Sphere coefficients of atom <alpha>: l, m, then (re, im) for each radial function of that l
    <l> <m> <re_1> <im_1> ... <re_{n_l}> <im_{n_l}>     ((L+1)^2 lines, l and m ascending as above)

The DFT functional is that of the self-consistent calculation the
orbitals come from, as the source code names it (the reference writer
reads WIEN2k's case.in0: XC_PBE gives PBE, XC_LDA LDA). The total energy
and the Fermi energy are that calculation's, from its last, converged
cycle (WIEN2k: the last :ENE and :FER of case.scf), in hartree; a reader
may set occupations from the Fermi energy. The Fermi energy is the
converged calculation's, and a band-mode eigenvalue of the valence
maximum can lie above it by rounding (FLEUR: 1e-10 Ha); a reader that
sets occupations from the Fermi energy compares with a tolerance,
1e-6 Ha say.

Kinds of radial function: u, the solution at the channel's energy; udot,
its energy derivative; lo, a local orbital's solution at its own energy.
The kind and the energy are labels for the reader's information and
readers do not use them: the coefficients and the tabulated functions
carry the orbital whatever the functions are. A writer whose extra
functions are of another sort (local orbitals built from second energy
derivatives, say, for which "energy" has no meaning) writes them as lo
with the energy of the channel they belong to.
The radial functions' normalisation is the source code's and need not be
stated (WIEN2k: u normalised to one in the sphere, udot orthogonal to u and
unnormalised, lo normalised); the coefficients absorb it.

Spin polarisation (version 1 writes F only): with T, each k point lists
n_up and n_down bands, the spin index of a band being 1 or 2, and every
radial function and core state is given once per spin in the ATOM TYPES
and CORE STATES sections, the spin-1 set first; version 1 readers may
refuse T.

## 3. What a reader may assume

- Spheres do not overlap; a point is in at most one sphere.
- The mesh is logarithmic, r_i = r_1 (r_2/r_1)^(i-1), its last point the
  sphere radius for valence functions; the core mesh continues it beyond
  R_a with the same ratio.
- Inside the join radius the table equals the series to rounding.
- Each band's coefficients give the SAME function as the plane-wave sum
  would inside the sphere to the matching's accuracy, so a reader may use
  either near the surface; a blend over a shell inside the sphere is the
  reader's choice.
- k points are given in Cartesian atomic units; the primitive reciprocal
  vectors follow from the lattice vectors.
