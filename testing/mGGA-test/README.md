# FLEUR metaGGA (SCAN) validation test suite

Source paper: F. Tran, G. Baudesson, J. Carrete, G. K. H. Madsen, P. Blaha, K. Schwarz, and D. J. Singh, "Shortcomings of meta-GGA functionals when describing magnetism", Phys. Rev. B 102, 024407 (2020); arXiv:2004.04543.

This directory contains input decks intended to check the `metaGGA`
branch of FLEUR (`iffgit.fz-juelich.de/fleur/fleur`, branch `metaGGA`)
against the **SCAN** meta-GGA results reported for itinerant
ferromagnets (Table III) and antiferromagnetic insulators/semimetals
(Table VI) of the paper above. Only the SCAN column is used, as
requested -- the other 13 functionals tabulated in the paper (LDA, PBE,
HLE16, mBJLDA, TPSS, revTPSS, MGGA_MS2, MVS, TM, HLE17, TASK, SCAN-L,
BR89) are not part of this suite.

## Layout

```
mGGA-test/
  README.md              <- this file
  reference_values.csv   <- machine-readable target values for all 16 cases
  fm/<Material>/          7 ferromagnetic test cases (Table III)
  af/<Material>/          9 antiferromagnetic test cases (Table VI)
```

Each `<Material>` directory contains:

- `README.md`      -- structure, magnetic order, and target values for that material
- `structure.inpgen` -- best-effort FLEUR input-generator (`inpgen`) freeform
  file (title line, `&lattice` namelist, atom list). **Treat this as a
  starting point, not a ready-to-run file** -- see caveat 1 below for the
  XC-functional part, and each AFM README for the magnetic-cell
  derivation.
- `inp.xml.notes`  -- what to add/change by hand once `inpgen` has produced
  `inp.xml` (XC functional switch, initial magnetic moments per atom
  group, spin-polarization, suggested k-mesh density).

## Ferromagnets (Table III, SCAN column)

| Material | Structure | SCAN mu_S (muB/f.u.) | Expt. |
|---|---|---|---|
| Fe | bcc (alpha-Fe) | 2.63 | 1.98 / 2.05 / 2.08 (three literature values, footnotes a-c) |
| Co | hcp | 1.79 | 1.52 / 1.55-1.62 / 1.58 (footnotes b-d) |
| Ni | fcc | 0.76 | 0.52 / 0.55 |
| FeCo | B2 (CsCl-type) | 4.86 | 4.54 |
| YFe2 | C15 Laves phase | 3.88 | 2.90 |
| ZrZn2 | C15 Laves phase | 1.08 | 0.178 |
| Ni3Al | L1_2 (Cu3Au-type) | 0.95 | 0.24 |

## Antiferromagnets (Table VI, SCAN column)

Values are the Bader-volume spin moment on the transition-metal atom, as
reported by the paper; see caveat 4 below for how to compare them to
FLEUR's muffin-tin moment.

| Material | Structure | TM atom | SCAN mu_S (muB) | Expt. (mu_S+mu_L) |
|---|---|---|---|---|
| MnO | rock salt, type-II AFM | Mn | 4.53 | 4.58 (footnote a) |
| FeO | rock salt, type-II AFM | Fe | 3.62 | 3.32 / 4.2 / 4.6 (footnotes b-d) |
| CoO | rock salt, type-II AFM | Co | 2.6 | 3.35 / 3.8 / 3.98 (footnotes e-g) |
| NiO | rock salt, type-II AFM | Ni | 1.6 | 1.9 / 2.2 (footnotes h,i) |
| Cr2O3 | corundum, AFM along c | Cr | 2.73 | 2.44 / 2.48 / 2.76 (footnotes k-m) |
| Fe2O3 | corundum (hematite), AFM along c | Fe | 4.01 | 4.17 / 4.22 (footnotes n,o) |
| CuO | tenorite (monoclinic), collinear AFM approximant | Cu | 0.57 | 0.65 (footnote j) |
| CrSb | NiAs-type, A-type AFM along c | Cr | 3.32 | 3.0 (footnote p) |
| CrSb2 | marcasite, collinear AFM | Cr | 3.18 | 1.94 (footnote g) |

## Important caveats -- please read before trusting these tests

This suite was generated without direct network access to
`iffgit.fz-juelich.de` (blocked by this sandbox's egress policy) and
without being able to download the paper's PDF/HTML as a single document
(only fetched in fragments through a web-fetch tool). Two consequences:

1. **XC-functional input.** The `inp.xml` files use the input the `metaGGA`
   branch expects:

   ```xml
   <xcFunctional name="LibXC" relativisticCorrections="F">
      <LibXCName exchange="mgga_x_scan" correlation="mgga_c_scan"/>
      <AuxGGA exchange="117" correlation="130"/>
   </xcFunctional>
   ```

   `<AuxGGA>` is required: it selects the GGA (here RPBE exchange with PBE
   correlation, the choice of Doumont et al., PRB 105, 195138 (2022)) used for
   the radial basis functions and the core states. fleur also requires
   `ctail="F"`, `frcor="F"` and `l_f="F"` for MetaGGA calculations. With inpgen,
   use `xctyp='LibXC: Exch: mgga_x_scan, Cor: mgga_c_scan, AuxExchID: 117, AuxCorID: 130'`;
   inpgen then sets `ctail="F"` itself.

2. **Structural data (Table I) and the SCAN/experimental reference values
   (Table III, Table VI) were extracted from the paper via automated PDF
   text extraction**, not by looking at typeset tables directly. The
   numbers were cross-checked between two independent fetches of the same
   tables and are internally consistent, but a manual spot check against
   the published PRB version (or the arXiv PDF) is recommended before
   using this as a "golden" reference suite.

3. **Antiferromagnetic magnetic cells.** The paper's Table I gives the
   *crystallographic* structure (columns "space group" list both the
   chemical space group and, in parentheses-like notation, the magnetic
   one, e.g. "MnO (225,166)"). It does not spell out the magnetic
   supercell vectors. For the rock-salt oxides (MnO, FeO, CoO, NiO) the
   type-II AFM rhombohedral cell used here was reconstructed by hand
   (see `af/<X>/README.md` for the derivation) and cross-checked via the
   well-known ideal rhombohedral angle of the type-II AFM cell,
   alpha = arccos(5/6) = 33.56 deg, which the derivation reproduces
   exactly. For Cr2O3, Fe2O3, CrSb, CrSb2 and CuO the atomic positions
   from Table I are used directly, but the assignment of which
   crystallographic transition-metal site carries "up" vs "down" spin is
   a best-effort reconstruction of the known magnetic structure (G-type /
   A-type order along the relevant axis) and is **not** independently
   verified against a magnetic-structure database (e.g. MAGNDATA) in this
   session. Please spot-check before relying on exact numbers, especially
   for CuO, whose real magnetic order is an incommensurate spin spiral;
   the collinear approximant used here is a common simplification in the
   DFT literature but is an approximation.

4. **AFM comparison quantity.** Per your instruction, the AFM tests
   compare against FLEUR's own **muffin-tin projected spin magnetic
   moment** (the standard `momentsAtom`/`magneticMomentsInMTSpheres`
   entry in `out.xml`), *not* the paper's Bader-volume moment used in
   Table VI. These two quantities differ by construction (different
   integration volume), so do not expect an exact numerical match to the
   Table VI number quoted in each README -- treat it as the right
   order of magnitude / right trend (SCAN moment noticeably larger than
   PBE), and use the "Expt." column as the ultimate sanity check instead
   if the MT-vs-Bader offset is a concern.

5. Only Fe and Ni have been run so far (October 2026, serial build, inputs as
   given here). The total spin moments per formula unit agree with the
   self-consistent WIEN2k values of Doumont et al., PRB 105, 195138 (2022),
   Table V, which are the better targets for the ferromagnets than the
   fixed-spin-moment values of the 2020 paper used in the table above:

   | | FLEUR PBE | WIEN2k PBE | FLEUR SCAN | WIEN2k SCAN |
   |---|---|---|---|---|
   | Fe | 2.23 | 2.22 | 2.60 | 2.60 |
   | Ni | 0.63 | 0.64 | 0.76 | 0.77 |

   The other cases are still input decks without verified results.
