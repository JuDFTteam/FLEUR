# CuO (antiferromagnet)

Reference: Table I and Table VI of F. Tran, G. Baudesson, J. Carrete, G. K. H. Madsen, P. Blaha, K. Schwarz, and D. J. Singh, "Shortcomings of meta-GGA functionals when describing magnetism", Phys. Rev. B 102, 024407 (2020); arXiv:2004.04543.

## Structure

- Structure type: tenorite (monoclinic), collinear AFM approximant
- Parent (chemical) space group: C2/c (#15)
- Magnetic space group: P2_1/c (#14)
- Lattice: monoclinic a = 4.684 Å, b = 3.423 Å, c = 5.129 Å, beta = 99.54 deg

## Magnetic structure notes -- read before use

CuO (tenorite) has a genuinely **incommensurate, non-collinear**
magnetic structure below its Neel temperature; the paper's
"(C2/c (#15),P2_1/c (#14))" notation refers to a
commonly used *collinear approximant* (breaking the C-centering of C2/c
down to the primitive P2_1/c cell, assigning opposite collinear spins to
the two C-centering-related Cu sites). This is the roughest
reconstruction in this test suite -- **verify the Cu/O site assignment
independently (e.g. against MAGNDATA or the original CuO neutron
diffraction literature) before trusting this test case**; treat a
moment mismatch here as more likely to be a structure-setup issue than
a SCAN-implementation bug.


## Target value

- **SCAN (Tran et al. 2020, Table VI), Bader spin moment on Cu: 0.57 muB**
- Experiment (total mu_S + mu_L): 0.65 (footnote j)
- **Comparison quantity in FLEUR: use the muffin-tin projected spin
  magnetic moment on the Cu atom** (per your instruction), not a
  Bader-volume moment -- FLEUR reports this directly per atom type in
  `out.xml`. Expect a systematic offset from the 0.57 muB target
  because muffin-tin and Bader volumes differ; judge success by the
  trend (SCAN moment clearly larger than PBE/LDA) and by the
  experimental value as the ultimate check, rather than by exact
  agreement with the paper's Bader number.

## Files

- `structure.inpgen` -- best-effort `inpgen` freeform input. See the
  derivation above for how the magnetic cell/spin assignment was
  constructed, and verify before trusting it.
- `inp.xml.notes` -- what to change in `inp.xml` after running `inpgen`.
