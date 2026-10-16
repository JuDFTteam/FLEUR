# (a) Theory page — FLEUR-specific background

Optional. Write one only if FLEUR does something the reader cannot learn from a textbook or the
original paper: its own formalism, representation, approximation, convention or algorithm, and
the reader needs it to understand the input parameters or interpret the results.

## The test for every paragraph

> Would this paragraph be equally true in the manual of any other DFT code?

If yes, cut it and cite a review or the original paper instead. Keep only what is specific to
FLEUR and the FLAPW method:

- **Representation**: how the quantity is expanded in the muffin-tin spheres, interstitial and
  vacuum regions (radial functions u, u̇, local orbitals, lm expansion, star functions,
  2D vacuum stars).
- **FLEUR conventions**: units (Hartree, bohr), sign and phase conventions, the order of indices
  (e.g. the `n_mmp` convention ⟨l'm's|ρ|lms'⟩), the global and local spin and orbital frames,
  the reciprocal-lattice units of input vectors.
- **Approximations as implemented**: what is neglected, e.g. spherical MT potential only,
  first-order perturbation, no SOC in the generalized Bloch theorem, and when that is still
  valid.
- **Algorithmic choices** that affect results: symmetrization, degenerate-state averaging,
  energy windows, mixing, and the default choice FLEUR makes.
- **Link from symbols to input**: every symbol the user controls is tied to its XPath.
- **Limitations**: which combinations (film, noco, SOC, LOs, MPI) are unsupported or untested.

## Sources

The code (comments, variable names, the actual formula implemented) and the developer who
wrote it. Ask the developer when the physics intent is not clear from the code; do not
reconstruct it from general knowledge. DOIs must come from code comments, `CITATION.cff`,
existing pages or the developer. Never write a DOI from memory without saying so to the user.

## Format (match `spinSpirals.md`, `nocoMag.md`, `socSecVar.md`)

- File: `docs/documentation/Theory/<topic>.md` (camelCase name). Nav: under
  `User Guide → Theoretical background` in `mkdocs.yml`; add it to `latex/include.list` after the
  last `documentation/Theory/` entry.
- Model pages are in `docs/documentation/Theory/`. Links to other kinds go across directories:
  `../Calculations/<feature>.md`, `../Reference/<feature>Setup.md`.
- `##` page title, no front matter.
- Inline math `$...$`; display math `\begin{equation} ... \end{equation}` or `$$ ... $$`.
  Define every symbol right after the equation it first appears in.
- Input switches go in `!!! note ""` blocks (body indented 4 spaces) using XPath notation,
  followed by an ```` ```xml ```` snippet if that helps.
- End with `#### Further reading`, a short list of linked references.
- No inline HTML. Images: PNG in `Theory/` with a name unique across `docs/documentation/`
  (the PDF converter copies all images into one directory),
  `![cap=Caption.,width=0.8\textwidth](file.png)`.

## Skeleton

```markdown
## <Topic as a noun phrase>

<One paragraph: which quantity or method this is and in which FLEUR calculations it appears.
Point to the howto page for running it: [..](../Calculations/<feature>.md).>

### Representation in FLEUR

<Equations as implemented, each symbol defined; the FLAPW-specific parts.>

!!! note ""
    <Which XPath sets which symbol, e.g. $\beta_\mu$ is `atomGroups/atomGroup/nocoParams/@beta`.>

### Approximations and conventions

<What is approximated, frames and units, defaults FLEUR chooses and why.>

### Limitations

<Unsupported combinations, known accuracy limits.>

#### Further reading

* <Original method paper>
* <Paper describing the FLEUR implementation, if any>
```

Link it from the howto page ("background: [..](../Theory/<topic>.md)") and from the tutorial's
Introduction (`https://www.flapw.de/rel/documentation/Theory/<topic>/`).
