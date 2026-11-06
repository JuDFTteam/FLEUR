<!-- Copyright (c) 2026 Mike Towler. MIT licence; see LICENSE in this directory. -->

# fleurqmc: notes for FLEUR's developers

Mike Towler, October 2026 (for FLEUR MaX-R8.1). Ported to FLEUR develop in November 2026; the
port is described where it changed something.

## What it is

A writer, inside FLEUR, of the all-electron LAPW orbitals of a converged ground state in a
code-neutral file, `lapwfn.data` (format version 1, `SPEC.md`), with the interstitial plane waves in
`pwfn.data`. The quantum Monte Carlo code CASINO reads the pair, and reads the same pair written by
WIEN2k and by Elk.

The file carries everything a reader needs to evaluate the orbitals anywhere in the crystal:
- the radial functions on each type's mesh
- the sphere coefficients in the global frame
- the plane-wave coefficients
- the core states, with an exponential tail beyond the sphere

It carries nothing of the code that wrote it.

## How it hooks in

One module, `src/fleur/io/writeLapwfn.F90` (`m_writeLapwfn`, public `write_lapwfn`). It is switched
on by the optional element `output/lapwfn` of `inp.xml`:
- attributes `totalEnergy`, `realise`, `title`
- read into `t_banddos` (`l_lapwfn`, `lapwfn_etot`, `lapwfn_realise`, `lapwfn_title`)
- schema type `LapwfnOutputType`

`fleur.F90` calls `write_lapwfn` on rank 0 after `eigen` and second variation, beside the SPEX
interface, and then ends the run with `juDFT_end`. The writer reads the eigenvectors with
`read_eig` and takes all bands of every k point of the selected list. No band mode is needed, so
the list may be of any kind.

The non-relativistic limit the file requires comes from `expertModes/@warp_factor="10000.0"`. On
MaX-R8.1 that attribute was never read and the contribution patched `constants.f90`; on develop it
works (#833). The writer tests `c_light(1.0)` and refuses a relativistic calculation.

On MaX-R8.1 the writer was called from `cdnval` per k point of a band-mode run started as
`fleur -lapwfn`, with a control file `lapwfn.in`. Both are gone: `lapwfn.in`'s keys are the
attributes above, and the functional's name now comes from `xcpot`.

The checks at points use the writer's own probe (FLEUR's coefficients and radial tables evaluated by
the writer at given points), which shares the writer's conventions. The conventions themselves (i^l,
the global frame, the plane-wave normalisation) are proved in two ways:
- by the continuity of the orbital across the sphere surface, where the plane-wave side is FLEUR's
  eigenvector with no sphere convention in it
- by the rotation test

## What it takes from FLEUR, and what it writes

- **Radial functions** come from `t_radfun` (u, udot and the local orbitals in the slots of
  `atoms%slot_of_lo`, divided by r), with the energy parameters, on each type's logarithmic mesh.
  The first points of each table are replaced by the regular power series that SPEC.md prescribes,
  joined at or below 0.05 bohr.
- **Sphere coefficients** come from `t_abc%calc_abc` with a copy of `sym` whose `ngopr` is 1 and
  `invsat` 0, so that they are in the global frame rather than each atom's local one. They are
  multiplied by i^l, which `hsmt_ab` leaves out. The optional cache of matching coefficients
  (`-abcoeff_store`) is switched off for this, because it holds local-frame coefficients.
- **Inversion partners** (`invsat` 2) are computed directly in the global frame. On develop the
  partner's local-orbital data point to the representative's eigenvector entries (`t_lapw`),
  whereas on MaX-R8.1 the writer had to form them from the partner's.
  - The writer checks the relation A'(l,m) = p (-1)^(l+m) conj A(l,-m), with p = exp(ik.T), on the
    LAPW and on the local-orbital coefficients, and stops if it fails.
  - Residual 0 in every case checked; `lapwfn_measurements.out`, `inversion` lines.
- **Core states** are computed by `cored`'s own recipe: `differ` on the spherical potential, with
  the confining extension beyond the sphere. Their occupation is written as an integer.
  - Beyond the sphere the function is u(R) exp(-kappa (r - R)) with kappa = -u'(R)/u(R),
    tabulated to where it is below 1e-12 of its maximum.
  - The whole function is normalised to one.
- **Eigenvectors** are divided by sqrt(omtil). The plane-wave set is the union over the k points of
  the list, ordered by |G|.
- **Bands at k = G/2** are made real and flagged:
  - A non-degenerate band is made real by one phase.
  - A degenerate level (eigenvalues within 1e-6 Ha) is replaced by real orthonormal combinations
    of its bands. They span the level exactly (its projector, and the occupied density, are
    unchanged) and carry their expectation energies as labels.
  - A level that does not close under conjugation, or that the Fermi energy divides, is left as
    FLEUR gives it, flagged complex, and the measurements say so.
  - The highest band computed is never written, since its level may be incomplete; the list ends at
    the last level boundary below it.
- **Energies:** the Fermi energy comes from the stored density (`readPrevEFermi`), the eigenvalues
  from `results%eig`, and the total energy from `output/lapwfn/@totalEnergy`.

Files written:
- `lapwfn.data`, `pwfn.data`
- `lapwfn_measurements.out`: the joins, the core states, the inversion relation, the phases, the
  levels
- `lapwfn_ylm.out`: `ylm4` at fifty directions
- `lapwfn_points.out`, when `lapwfn_points.in` is present: FLEUR's own orbital values at given
  points, from the same coefficients and tables, which is the reference the checks use

## What was checked

### On develop (the port)

Four of the contributor's cases were re-run end to end on develop: SCF with `warp_factor="10000.0"`,
then the writer. The setup differs from MaX-R8.1 in two ways:
- the built-in `l91` LDA instead of libxc `lda_x` + `lda_c_pw`
- trigonal carbon on a coarser SCF mesh

Checkers used: `checks/spec_check.py`, `check_probe.py` and `check_scf.py`. The others need the
external `wien2kqmc` reader.

| Case | Result |
|---|---|
| diamond | SCF -75.5944073791 Ha (MaX-R8.1 with libxc: -75.5944091505), the same 12 iterations and last distance; spec_check 0 failures; frame test identical to MaX-R8.1 (1.350e-3, 1.018e-3); band, Fermi and total energies against the SCF 7e-11, 5e-11, 0; core 1s 3e-10 Ha |
| diamond with a 3p local orbital | inversion relation residual 0 on the LAPW and on the local-orbital coefficients at every k point; local-orbital coefficients non-zero and equal on both atoms; spec_check 0 failures; frame test identical (1.387e-3) |
| silicon | 2p core from both kappa solutions merged (1.65e-7 Ha apart); core energies against out.xml below 1e-10 Ha; spec_check 0 failures; frame test identical (3.892e-3) |
| trigonal carbon (P-6) | right frame 6.4e-4, transposed frame 1.74 (MaX-R8.1: 6.3e-4, 1.74); spec_check 0 failures |

The regression test `testing/tests/parameterized/test_lapwfn.py` (inputs
`testing/inputfiles/lapwfn/DiamondLO`) runs a small version of the local-orbital case.

**The near-nucleus start of the radial solver changed on develop** (b156539a5), and the inner
series joins move with it:
- s functions: the table is now regular down to the first point (join at 1.0e-4 bohr, MaX-R8.1
  2.2e-3).
- l >= 1, with `warp_factor="10000.0"`: the first points are wrong (the new Frobenius start
  diverges for large c, issue #834). The writer's inner series replaces them (join at 4.6e-3 bohr
  for p), within the 0.05 bohr the format allows.

### On MaX-R8.1 (the contribution)

Cases: diamond (a = 6.7407 bohr, rmt 1.3, lmax 8, rgkmax 9) and five variants of it:
- a doubled radial mesh with lmax 10
- a 3p local orbital on the carbon
- the p energy parameter fixed at 0.40 Ha by an enpara file
- zincblende BN, wurtzite BN
- a trigonal carbon of space group P-6

and silicon (a 2p core with both kappa solutions). All were LDA (libxc `lda_x`, `lda_c_pw`),
non-relativistic, converged to a distance of 1e-7, at Gamma and the k points of a cubic (or
doubled) simulation cell.

| Check | Result |
|---|---|
| FLEUR's own orbitals at random points inside the spheres, per band (the probe) | below 6e-7 (standard mesh), 5e-8 (fine mesh), 7e-8 (silicon) |
| At points on the mesh radii | 5e-15 |
| In the interstitial | 6e-15 |
| FLEUR's own values at points just inside the sphere surface | 3e-10 |
| Every band at k = G/2 flagged real, and real at the points | imaginary part below 2e-12 of the maximum, every case |
| Degenerate levels at k = G/2 as real combinations (closure, real span, projector, labels, occupied density) | closure below 5e-13; projectors 6e-15, labels 6e-15 Ha, density 7e-16 |
| Value and slope across the surface, channel by channel | value 1e-12; slope 1.5e-6 by the surface checker's rule, 3e-8 by polynomial fits |
| Pointwise jump of the orbital across the surface (truncation at lmax) | 1.7e-3 of the band's maximum (diamond, lmax 8) to 1.4e-2 (zincblende BN, its highest unoccupied bands) |
| Harmonics, `ylm4` against the reader's | 5e-13 |
| Global frame: the probe rotated by the atom's operation against the plane waves, and by its transpose (trigonal carbon) | 6e-4 and 1.7 |
| Band energies at Gamma against the SCF's, Fermi energy, total energy | 1.5e-10 Ha, 1e-10 Ha, equal |
| Core energies against out.xml | 5e-10 Ha |
| Core norm over all space by an independent quadrature | 1 to 5e-7 |
| A general k point (twist 0.13, 0.27, 0.41) | sphere 1.7e-7, interstitial 2e-15, surface value 5e-13; every band complex |
| Diamond with p at 0.40 Ha against WIEN2k and Elk | occupied density at Gamma within 9e-4 of its maximum at 200 points against WIEN2k, 1.6e-3 against Elk |
| Two writer runs from the same stored density | byte-identical but for the time stamp (also on develop) |

Against WIEN2k's and Elk's files for the same diamond and BN:
- The occupied band energies agree to a few 1e-4 Ha after a constant shift (the potential zero).
- The density of the occupied states at Gamma differs by one to two per cent between 0.2 and
  0.5 bohr from the carbon nucleus.

The three codes' own orbitals are each reproduced by their files. The difference is FLEUR's default
p energy parameter (see README.md).

## Found in MaX-R8.1, and the state on develop

1. **`expertModes@warp_factor` had no effect** (the constants object was in neither the input
   traversal nor the broadcast). Fixed on develop (#833).
2. **The default make target built `fleuriste`**, which installs Python packages. It is gone from
   develop.
3. **Built-in functional names:** inpgen's `pz` is Perdew-Zunger 1981, and `l91` is the
   Perdew-Wang uniform-gas correlation. Unchanged.
4. **Every run rewrites `cdn.hdf` at its start** (`qfix`, `writeDensity` in `fleur.F90`). This
   still holds; see README.md.
5. **Band mode warns for a k-point list not of kind `path`.** It no longer matters, because the
   writer does not need band mode.
6. **`t_lapw%nkvec` is unallocated without local orbitals.** Develop always allocates it.
7. **The start of the scalar-relativistic radial solver** (develop only, b156539a5) diverges for
   large c at l >= 1 (#834; see above).
