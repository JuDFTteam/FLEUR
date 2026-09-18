# CoO (antiferromagnet)

Reference: Table I and Table VI of F. Tran, G. Baudesson, J. Carrete, G. K. H. Madsen, P. Blaha, K. Schwarz, and D. J. Singh, "Shortcomings of meta-GGA functionals when describing magnetism", Phys. Rev. B 102, 024407 (2020); arXiv:2004.04543.

## Structure

- Structure type: rock salt, type-II AFM
- Parent (chemical) space group: Fm-3m (#225)
- Magnetic space group: R-3m (#166)
- Lattice: cubic parent a = 4.254 Å

## Derivation of the type-II AFM magnetic cell

CoO is rock salt (Fm-3m, #225) with the Co cation on the
fcc sublattice and O on the interpenetrating fcc sublattice offset by
(1/2,0,0). Table I of the paper labels the structure "(Fm-3m (#225),R-3m (#166))"
which -- with no separate magnetic lattice parameters listed -- means
the *chemical* lattice stays the ideal cubic rock salt lattice, and the
reduced symmetry (R-3m, #166) comes purely from the type-II
antiferromagnetic order (ferromagnetic (111) planes, alternating in
sign from one (111) plane to the next).

Doubling the primitive fcc cell along the 3-fold [111] axis while
keeping 3-fold symmetry gives the rhombohedral magnetic cell (Cartesian,
units of the cubic lattice constant a = 4.254 Angstrom):

```
a1 = a (0.5, 0.5, 1.0)
a2 = a (1.0, 0.5, 0.5)
a3 = a (0.5, 1.0, 0.5)
```

with basis

```
Co (spin up)    (0,    0,    0   )
Co (spin down)  (0.5,  0.5,  0.5 )
O                       (0.25, 0.25, 0.25)
O                       (0.75, 0.75, 0.75)
```

This cell has volume 0.5 a^3 (2 formula units, as required for two
magnetically inequivalent Co sites) and a rhombohedral angle
alpha = arccos(5/6) = 33.5573 deg -- the well-known ideal angle quoted
for the undistorted type-II AFM rock-salt cell (real materials show a
small additional magnetostrictive distortion away from this ideal
value, which the paper's Table I does not include, so it is not
reproduced here either).

**This construction was derived by hand in this session** (cross-checked
via the well-known 33.5573 deg angle, which it reproduces exactly) and
was not checked against the `metaGGA` branch or a magnetic-structure
database -- please verify the site assignment before trusting the test.


## Target value

- **SCAN (Tran et al. 2020, Table VI), Bader spin moment on Co: 2.6 muB**
- Experiment (total mu_S + mu_L): 3.35 / 3.8 / 3.98 (footnotes e-g)
- **Comparison quantity in FLEUR: use the muffin-tin projected spin
  magnetic moment on the Co atom** (per your instruction), not a
  Bader-volume moment -- FLEUR reports this directly per atom type in
  `out.xml`. Expect a systematic offset from the 2.6 muB target
  because muffin-tin and Bader volumes differ; judge success by the
  trend (SCAN moment clearly larger than PBE/LDA) and by the
  experimental value as the ultimate check, rather than by exact
  agreement with the paper's Bader number.

## Files

- `structure.inpgen` -- best-effort `inpgen` freeform input. See the
  derivation above for how the magnetic cell/spin assignment was
  constructed, and verify before trusting it.
- `inp.xml.notes` -- what to change in `inp.xml` after running `inpgen`.
