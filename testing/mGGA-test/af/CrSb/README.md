# CrSb (antiferromagnet)

Reference: Table I and Table VI of F. Tran, G. Baudesson, J. Carrete, G. K. H. Madsen, P. Blaha, K. Schwarz, and D. J. Singh, "Shortcomings of meta-GGA functionals when describing magnetism", Phys. Rev. B 102, 024407 (2020); arXiv:2004.04543.

## Structure

- Structure type: NiAs-type, A-type AFM along c
- Parent (chemical) space group: P6_3/mmc (#194)
- Magnetic space group: P-3m1 (#164)
- Lattice: hexagonal a = 4.122 Å, c = 5.464 Å

## Magnetic structure notes

CrSb is the NiAs structure type (P6_3/mmc, #194), which already
contains **two** Cr atoms per hexagonal cell at (0,0,0) and
(0,0,1/2) -- no supercell is needed for the A-type AFM order (ferromagnetic
(0001) planes stacking antiferromagnetically along c). Assigning
opposite spins to these two already-inequivalent Wyckoff sites lowers
the symmetry from P6_3/mmc (#194) to P-3m1 (#164),
matching Table I's "(P6_3/mmc (#194),P-3m1 (#164))" notation.


## Target value

- **SCAN (Tran et al. 2020, Table VI), Bader spin moment on Cr: 3.32 muB**
- Experiment (total mu_S + mu_L): 3.0 (footnote p)
- **Comparison quantity in FLEUR: use the muffin-tin projected spin
  magnetic moment on the Cr atom** (per your instruction), not a
  Bader-volume moment -- FLEUR reports this directly per atom type in
  `out.xml`. Expect a systematic offset from the 3.32 muB target
  because muffin-tin and Bader volumes differ; judge success by the
  trend (SCAN moment clearly larger than PBE/LDA) and by the
  experimental value as the ultimate check, rather than by exact
  agreement with the paper's Bader number.

## Files

- `structure.inpgen` -- best-effort `inpgen` freeform input. See the
  derivation above for how the magnetic cell/spin assignment was
  constructed, and verify before trusting it.
- `inp.xml.notes` -- what to change in `inp.xml` after running `inpgen`.
