# Co (ferromagnet)

Reference: Table I and Table III of F. Tran, G. Baudesson, J. Carrete, G. K. H. Madsen, P. Blaha, K. Schwarz, and D. J. Singh, "Shortcomings of meta-GGA functionals when describing magnetism", Phys. Rev. B 102, 024407 (2020); arXiv:2004.04543.

## Structure

- Structure type: hcp
- Space group: P6_3/mmc (#194)
- Lattice: hexagonal, a = 2.507 Å, c = 4.070 Å
- Atoms: Co: 2c (1/3,2/3,1/4), (2/3,1/3,3/4)

## Magnetic order

Simple collinear ferromagnet -- no special magnetic cell needed. Run a
standard `jspins=2` spin-polarized calculation, let the moment relax
self-consistently (the paper used the fixed-spin-moment, FSM, method to
scan the energy vs. moment curve and pick the minimum for the meta-GGA
functionals -- try plain SCF first, and only fall back to FSM if SCAN
converges to a spurious/metastable moment).

## Target value

- **SCAN (Tran et al. 2020, Table III): 1.79 muB/f.u.**
- Experiment: 1.52 / 1.55-1.62 / 1.58 (footnotes b-d)

Simple collinear FM, no supercell needed.

## Files

- `structure.inpgen` -- best-effort `inpgen` freeform input (title line +
  `&lattice` namelist + atom list with atomic number and fractional
  coordinates). Verify the `latsys` code and Wyckoff-position expansion
  against your inpgen version before trusting it.
- `inp.xml.notes` -- what to change in `inp.xml` after running `inpgen`.
