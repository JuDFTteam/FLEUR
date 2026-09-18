# CrSb2 (antiferromagnet)

Reference: Table I and Table VI of F. Tran, G. Baudesson, J. Carrete, G. K. H. Madsen, P. Blaha, K. Schwarz, and D. J. Singh, "Shortcomings of meta-GGA functionals when describing magnetism", Phys. Rev. B 102, 024407 (2020); arXiv:2004.04543.

## Structure

- Structure type: marcasite, collinear AFM
- Parent (chemical) space group: Pnnm (#58)
- Magnetic space group: P2_1/c (#14)
- Lattice: orthorhombic a = 6.028 Å, b = 6.874 Å, c = 3.272 Å

## Magnetic structure notes

CrSb2 (marcasite structure, orthorhombic) has two Cr sites per
cell -- (0,0,0) and (1/2,1/2,1/2) in the standard Pnnm setting -- which
become magnetically inequivalent under AFM order, lowering the symmetry
from Pnnm (#58) to P2_1/c (#14). **The exact
Wyckoff multiplicity/origin choice for the Sb positions was
reconstructed from the single representative position given in Table I
(Cr(0,0,0), Sb(0.1835,0.3165,0.32)) and was not
cross-checked against a crystallographic database in this session** --
please verify against e.g. the ICSD/COD entry for CrSb2 before trusting
this test.


## Target value

- **SCAN (Tran et al. 2020, Table VI), Bader spin moment on Cr: 3.18 muB**
- Experiment (total mu_S + mu_L): 1.94 (footnote g)
- **Comparison quantity in FLEUR: use the muffin-tin projected spin
  magnetic moment on the Cr atom** (per your instruction), not a
  Bader-volume moment -- FLEUR reports this directly per atom type in
  `out.xml`. Expect a systematic offset from the 3.18 muB target
  because muffin-tin and Bader volumes differ; judge success by the
  trend (SCAN moment clearly larger than PBE/LDA) and by the
  experimental value as the ultimate check, rather than by exact
  agreement with the paper's Bader number.

## Files

- `structure.inpgen` -- best-effort `inpgen` freeform input. See the
  derivation above for how the magnetic cell/spin assignment was
  constructed, and verify before trusting it.
- `inp.xml.notes` -- what to change in `inp.xml` after running `inpgen`.
