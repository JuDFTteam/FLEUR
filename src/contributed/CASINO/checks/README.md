<!-- Copyright (c) 2026 Mike Towler. MIT licence; see LICENSE in src/contributed/CASINO. -->

# Checks of written files

The contributor's checking scripts for `lapwfn.data` / `pwfn.data`, as he used them to validate
the writer. A "deck" is a directory holding a writer run's files:
- `lapwfn.data`, `pwfn.data`, `lapwfn_measurements.out`
- `lapwfn_points.out` and `lapwfn_ylm.out`, for the checks that use them

Every check ends with a `RESULT: PASS` or `RESULT: FAIL ...` line (`spec_check.py`: a
`# n failures` line) and exits non-zero on a failure.

## Dependencies

- Python 3 with numpy, and scipy for some of the checks.
- All but `spec_check.py`, `check_probe.py` and `check_scf.py` need the reference reader
  `wien2kqmc`, which is not part of FLEUR. The environment variable `WIEN2KQMC` names its tree
  (default `~/PROJECTS/wien2kqmc`).

## Preparing a deck

1. Converge the SCF and prepare the writer run (`../README.md`).
2. For the point checks, make points from a first `lapwfn.data` with `make_points.py`. It writes per
   atom points on mesh radii, random sphere points and points just inside the surface, plus
   interstitial points.
3. Run the writer again with that `lapwfn_points.in` in the directory. Its `lapwfn.data` is
   byte-identical to the first but for the time stamp.
4. For `check_levels.py`, run the writer once more with `realise="F"` into another directory.

## The checks

- **`spec_check.py <deck>`** checks the file against SPEC version 1, clause by clause. It is an
  independent reader of the format's text. It also checks every band's norm and overlaps, and the
  realness relations of every band flagged real.
- **`check_probe.py <lapwfn_points.out> <lapwfn_points.in> NM NS NSF nat`** is the frame test, from
  FLEUR's own values alone. At the surface points, the local-frame value with the right rotation must
  match the plane-wave sum and the transposed one must not, where the two differ.
- **`check_scf.py <lapwfn.data> <SCF out.xml> <lapwfn_measurements.out>`** checks the file against
  the SCF:
  - the eigenvalues at Gamma
  - the Fermi and total energies
  - the core energies and leakage
  - the electron count
- **`check_ylm.py <lapwfn_ylm.out>`** checks FLEUR's `ylm4` against the reference reader's
  harmonics, l <= 10.
- **`check_points.py`** compares the reference evaluator (wien2kqmc) on the written files with
  FLEUR's own orbital values at the points:
  - on mesh radii, at random sphere points, at the surface and in the interstitial, per band by
    its own maximum
  - a band flagged real must evaluate real at the points
- **`check_core.py <deck>`** checks the core states:
  - the norm over all space (the series from the origin, the table, the tail)
  - the tail's norm three ways, and the tail against the exponential at every written point
  - the last point against the floor
- **`check_mesh.py <deck>`** checks the radial mesh's adequacy for a cubic-spline reader.
- **`check_levels.py <deck> <realise F deck> <lapwfn_points.in> NOCC`** compares the deck level by
  level with the `realise="F"` deck:
  - the labels and the level projectors
  - the closure under conjugation and the singular values of the real span
  - the real flags, and the occupied density at the points
- **`slope_check.py <deck>`** checks the slope at the sphere surface channel by channel, from
  polynomial fits in ln r through the last 4, 6, 8 and 12 mesh points.
- **`check_wien2k.py <deck A> <deck B> NOCC NB`** and **`density_profile.py <deck A> <deck B> NOCC`**
  compare with another code's file for the same crystal:
  - the electron counts
  - the Gamma band energies after a constant shift
  - the occupied density at random points and along radial lines

`../NOTES_FOR_DEVELOPERS.md` lists the results on MaX-R8.1 and on develop.
