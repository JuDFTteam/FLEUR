<!-- Copyright (c) 2026 Mike Towler. MIT licence; see LICENSE in this directory. -->

# FLEUR orbitals for CASINO (fleurqmc)

FLEUR writes the all-electron LAPW orbitals of a converged ground state in the code-neutral format
`lapwfn.data` (format version 1, `SPEC.md` here), with the interstitial plane waves in CASINO's
`pwfn.data`. CASINO's LAPW orbital type reads the pair. WIEN2k's converter (wien2kqmc) and Elk's
writer (elkqmc) write the same files, and a reader cannot tell which code they came from.
Contributed by Mike Towler (MIT licence, `LICENSE`).

What is here:

- `README.md` (this file): how to run the writer.
- `SPEC.md`: the file format.
- `NOTES_FOR_DEVELOPERS.md`: what the writer takes from FLEUR and how it was checked.
- `tools/fleurqmc_klist.py`: the k points of CASINO's simulation cell and twist as a FLEUR k-point list,
  and the `inp.xml` settings of a writer run.
- `checks/`: the contributor's checking scripts for written files (see `checks/README.md`).

The writer itself is `src/fleur/io/writeLapwfn.F90`, switched on by the `output/lapwfn` element of
`inp.xml`.

## Running the writer

1. **Converge the ground state non-relativistically**, with `minDistance` set. The file carries
   non-relativistic orbitals, so the speed of light is scaled up for the SCF and the writer run:

   ```xml
   <expertModes ... warp_factor="10000.0"/>
   ```

   Write the number as `10000.0`. FLEUR's expression parser does not accept `1.0e4`. The writer
   refuses a run with c below 1e5, and a stored density whose last distance is above `minDistance`.

2. **Prepare the writer run** in a copy of the converged directory:

   ```
   fleurqmc_klist.py M11 M12 M13 M21 M22 M23 M31 M32 M33 [t1 t2 t3] --inp inp.xml --kpts kpts.xml \
                     --numbands N --etot E [--title TEXT]
   ```

   The rows of M are the supercell vectors in the primitive basis (integers), t the twist, N the
   number of bands wanted (more than the occupied ones) and E the converged total energy
   (`totalEnergy` in the SCF's `out.xml`). The tool adds the k points as a list `qmc` to
   `kpts.xml`, selects it, sets `numbands` to N + 1, and adds

   ```xml
   <output ...>
      <lapwfn totalEnergy="E" realise="T" title="TEXT"/>
   ```

   The same can be done by hand. `totalEnergy` is required: a run that does not iterate cannot know
   it. `realise` (default `T`) replaces each degenerate level at k = G/2 by real combinations (see
   below); `F` leaves it as FLEUR gives it. `title` is written to both files.

   FLEUR is asked for N + 1 bands because the writer never writes the highest band it computed,
   whose degenerate level may be incomplete. The list ends at the last level boundary below it
   (eigenvalues within 1e-6 Ha make a level), so at most N bands are written at a k point. Fewer are
   written when band N + 1 is degenerate with those below it, and `lapwfn_measurements.out` says
   which bands were dropped. FLEUR raises a request below its own floor of bands to that floor.

3. **Optionally,** put `lapwfn_points.in` (a count, then Cartesian points in bohr) in the directory.
   FLEUR's own orbital values at those points are then written to `lapwfn_points.out`, the
   reference for a reader's test (`checks/make_points.py` makes such a file).

4. **Run `fleur`** (one process). After the first diagonalisation on the stored density it writes
   these files and stops:
   - `lapwfn.data`, `pwfn.data`
   - `lapwfn_measurements.out`: the joins of the inner series, the core states, the inversion
     relation, the levels and phases
   - `lapwfn_ylm.out`: FLEUR's harmonics at fifty directions

   Each file is written under its name with `.part` appended and renamed when complete. A refusal
   removes the `.part` files, so a directory that holds a good set of files keeps it through a run
   that fails.

   Every FLEUR run rewrites `cdn.hdf` at its start: the stored density is written back after the
   charge normalisation, re-rounded. Copy the SCF directory first if the stored density matters.
   Two writer runs from the same stored density give files that are byte-identical but for the
   time stamp. A run on a rewritten `cdn.hdf` differs in the last digit of some coefficients
   (1e-15 relative).

## Refusals

The writer refuses, with a message naming the reason:
- the calculation setup: spin polarisation (jspins 2), non-collinear and spin-spiral runs,
  spin-orbit coupling, film geometry
- the basis: an APW channel (an APW+lo local orbital, or l above `lmaxAPW`), a relativistic
  (Dirac) local orbital
- a relativistic calculation (c below 1e5, see `warp_factor`)
- the stored density: not converged (no SCF distance, or one above `minDistance`), or written with
  another spin count than the input's
- missing `totalEnergy`, more than one MPI process, a frozen core, `kcrel` other than 0,
  Hubbard-1 states, a fractional number of valence electrons
- at a k point: fewer bands left than occupied ones, or every written band at or below the Fermi
  energy within 1e-6 Ha. A metal's occupied bands are never cut, so `numbands` must leave a band
  above the Fermi energy at every k point.

All but the last two are made before any file is opened. The last two are made at each k point and
remove the `.part` files.

Guards on what FLEUR hands the writer refuse in the same way:
- inversion partners with complex eigenvectors, or partners whose coefficients break the inversion
  relation
- more radial functions at one l than the format allows
- core states: a fractional occupation, a solution that fails, vanishes or grows at the sphere, or
  two kappa solutions of one core state that differ by more than 1e-6 Ha
- a series fit that fails, and a G vector missing from the union

The Fermi energy the file carries is the converged SCF's, as the format asks. The eigenvalue of the
valence maximum in the writer run can lie above it by about 1e-10 Ha. A reader that sets occupations
from it should therefore compare with a tolerance, as the writer's own guard does.

## Real orbitals at k = G/2

At a k point with 2k a reciprocal lattice vector, every band written is real and flagged T:
- A non-degenerate band is made real by one phase.
- A degenerate level (consecutive eigenvalues within 1e-6 Ha) is replaced by real orthonormal
  combinations of its bands. They span the level exactly, so its projector is unchanged (to 1e-12
  in the cases checked) and so is the occupied density. They are labelled by their expectation
  energies, which lie within the level's spread, and they are what CASINO's real-orbital path needs
  for a level that is wholly occupied.

The level must close under the conjugation of a real function:
- plane waves: c(G) to conj c(-G-2k)
- sphere coefficients: c(l,m) to (-1)^m conj c(l,-m), in the global frame

It closes when FLEUR's eigenvectors are exact and the basis is complete under conjugation, as it
was at every level of the cases checked (residual 1e-14). Three kinds of level are written as FLEUR
gives them, flagged F:
- a level that does not close
- a level divided by the Fermi energy
- a level left complex by `realise="F"`

The `level` lines of `lapwfn_measurements.out` say which applies: the closure residual, the singular
values of the real span, the projector and label checks, and the reason. A level that FLEUR's own
vectors already make real band by band, as every Gamma level of the checked cases is, is kept as
FLEUR's.

## The p energy parameter

FLEUR's default energy parameters (inpgen's `energyParameters`, the band-centre search) put
diamond's p linearisation energy 0.57 Ha below the top of the valence band:
- The total energy and the eigenvalues feel this only to 0.1 mHa.
- The occupied p orbitals at Gamma, however, are about two per cent low inside the sphere, between
  0.2 and 0.5 bohr from the nucleus.

The file carries those orbitals. This is FLEUR's density, reproduced to 6e-7, not the writer's
doing.

For QMC, either fix the p parameter near the occupied p states with an `enpara` file, or add a p
local orbital. With p at 0.40 Ha:
- the density agrees with WIEN2k's to 6e-4 between 0.2 and 0.5 bohr from the nucleus (9e-4 of its
  maximum at 200 random points)
- the total energy moves by 0.08 mHa

The same holds for any element whose occupied states lie far from the band centre FLEUR chooses.
