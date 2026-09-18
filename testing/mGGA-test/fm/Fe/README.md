# Fe (ferromagnet)

Reference: Table I and Table III of F. Tran, G. Baudesson, J. Carrete, G. K. H. Madsen, P. Blaha, K. Schwarz, and D. J. Singh, "Shortcomings of meta-GGA functionals when describing magnetism", Phys. Rev. B 102, 024407 (2020); arXiv:2004.04543.

## Structure

- Structure type: bcc (alpha-Fe)
- Space group: Im-3m (#229)
- Lattice: cubic, a = 2.867 Å (conventional cell; primitive bcc cell edge = a*sqrt(3)/2)
- Atoms: Fe: 2a (0,0,0) [primitive bcc cell: 1 Fe at (0,0,0)]

## Magnetic order

Simple collinear ferromagnet -- no special magnetic cell needed. Run a
standard `jspins=2` spin-polarized calculation, let the moment relax
self-consistently (the paper used the fixed-spin-moment, FSM, method to
scan the energy vs. moment curve and pick the minimum for the meta-GGA
functionals -- try plain SCF first, and only fall back to FSM if SCAN
converges to a spurious/metastable moment).

## Target value

- **SCAN (Tran et al. 2020, Table III): 2.63 muB/f.u.**
- Experiment: 1.98 / 2.05 / 2.08 (three literature values, footnotes a-c)

Simple collinear FM, no supercell needed.

## Files

- `structure.inpgen` -- best-effort `inpgen` freeform input (title line +
  `&lattice` namelist + atom list with atomic number and fractional
  coordinates). Verify the `latsys` code and Wyckoff-position expansion
  against your inpgen version before trusting it.
- `inp.xml.notes` -- what to change in `inp.xml` after running `inpgen`.
