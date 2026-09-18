# Cr2O3 (antiferromagnet)

Reference: Table I and Table VI of F. Tran, G. Baudesson, J. Carrete, G. K. H. Madsen, P. Blaha, K. Schwarz, and D. J. Singh, "Shortcomings of meta-GGA functionals when describing magnetism", Phys. Rev. B 102, 024407 (2020); arXiv:2004.04543.

## Structure

- Structure type: corundum, AFM along c
- Parent (chemical) space group: R-3c (#167)
- Magnetic space group: R3 (#146)
- Lattice: hexagonal a = 4.953 Å, c = 13.588 Å

## Magnetic structure notes

Cr2O3 is corundum-derived: parent (nonmagnetic-symmetry) space
group R-3c (#167), magnetic space group
R3 (#146) once the AFM order along c is imposed (spins
alternate along the c-axis sequence of Cr sites at
z, z+1/2, -z, -z+1/2 within the rhombohedral cell). The atomic positions
below are taken directly from Table I of the paper:
Cr(0,0,0.3475), O(0.3058,0,1/4).

**The up/down spin assignment to the two crystallographically related
Cr sites (z vs z+1/2) follows the commonly reported G-type-like
AFM order for this family of corundum oxides but was not independently
verified against a magnetic-structure database in this session** --
please double check before treating a mismatch as a code bug rather
than a wrong spin assignment.


## Target value

- **SCAN (Tran et al. 2020, Table VI), Bader spin moment on Cr: 2.73 muB**
- Experiment (total mu_S + mu_L): 2.44 / 2.48 / 2.76 (footnotes k-m)
- **Comparison quantity in FLEUR: use the muffin-tin projected spin
  magnetic moment on the Cr atom** (per your instruction), not a
  Bader-volume moment -- FLEUR reports this directly per atom type in
  `out.xml`. Expect a systematic offset from the 2.73 muB target
  because muffin-tin and Bader volumes differ; judge success by the
  trend (SCAN moment clearly larger than PBE/LDA) and by the
  experimental value as the ultimate check, rather than by exact
  agreement with the paper's Bader number.

## Files

- `structure.inpgen` -- best-effort `inpgen` freeform input. See the
  derivation above for how the magnetic cell/spin assignment was
  constructed, and verify before trusting it.
- `inp.xml.notes` -- what to change in `inp.xml` after running `inpgen`.
