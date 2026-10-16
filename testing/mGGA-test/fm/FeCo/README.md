# FeCo (ferromagnet)

Reference: Table I and Table III of F. Tran, G. Baudesson, J. Carrete, G. K. H. Madsen, P. Blaha, K. Schwarz, and D. J. Singh, "Shortcomings of meta-GGA functionals when describing magnetism", Phys. Rev. B 102, 024407 (2020); arXiv:2004.04543.

## Structure

- Structure type: B2 (CsCl-type)
- Space group: Pm-3m (#221)
- Lattice: cubic, a = 2.857 Å
- Atoms: Fe: 1a (0,0,0);  Co: 1b (1/2,1/2,1/2)

## Magnetic order

Simple collinear ferromagnet -- no special magnetic cell needed. Run a
standard `jspins=2` spin-polarized calculation, let the moment relax
self-consistently (the paper used the fixed-spin-moment, FSM, method to
scan the energy vs. moment curve and pick the minimum for the meta-GGA
functionals -- try plain SCF first, and only fall back to FSM if SCAN
converges to a spurious/metastable moment).

## Target value

- **SCAN (Tran et al. 2020, Table III): 4.86 muB/f.u.**
- Experiment: 4.54

Total moment 4.86 muB/f.u.; paper's Bader-partitioned atomic moments in parentheses are Fe: 2.99, Co: 1.87 muB -- useful if you want to compare atom-resolved MT moments too.

## Files

- `structure.inpgen` -- best-effort `inpgen` freeform input (title line +
  `&lattice` namelist + atom list with atomic number and fractional
  coordinates). Verify the `latsys` code and Wyckoff-position expansion
  against your inpgen version before trusting it.
- `inp.xml.notes` -- what to change in `inp.xml` after running `inpgen`.
