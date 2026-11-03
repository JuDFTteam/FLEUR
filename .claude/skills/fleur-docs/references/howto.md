# (b) Howto page — how to perform the calculation

Task-oriented: the reader wants to do X with FLEUR and needs the steps, the input to set, and
what output to expect. Theory only as far as needed to make a choice (one or two sentences,
then link the theory page).

## Content

1. **What you get and when to use it** — one paragraph, including prerequisites (converged
   density, `l_noco`, SOC, a special k-point set, a film setup, ...).
2. **Steps**, in the order the user performs them. Each step says what to do and why, and gives:
   - the input as XPath in a `!!! note ""` block, plus an ```` ```xml ```` snippet,
   - commands where needed: `inpgen -f <file>`, `inpgen -inp.xml -kpt grid=...`, `fleur`,
     `fleur_MPI`, `fleur -xmlXPath "<path>=<value>"`.
3. **Output** — every file, out.xml element and HDF dataset the user needs, with how to read it:
   `grep <tag> out.xml`, a short Python/h5py or masci-tools snippet, the meaning of columns.
4. **Example** — a complete minimal case, taken from `testing/inputfiles/<dir>/<test>/` when one
   exists (say which test). The inpgen input or inp.xml excerpt, the changes, and **expected
   numbers from a real run** (the test's reference values or a run you did; state system and
   parameters). Never write example numbers that did not come out of FLEUR.
5. **Practical hints** — convergence parameters that matter, typical pitfalls, cost, MPI notes.
6. Close with:

```markdown
!!! note "More info"
    All parameters are listed in the [<feature> reference](../Reference/<feature>Setup.md).
    An interactive example is the tutorial <Xn> "<title>" in the FLEUR tutorial container.
```

## Format (match `dos.md`, `ldaU.md`, `bandstructure.md`, `phonons.md`)

- File: `docs/documentation/Calculations/<feature>.md` (camelCase), for both nav sections:
  `User Guide → Basic calculations` for standard tasks, `Expert knowledge and troubleshooting`
  for specialised ones. A topic with many pages gets its own subdirectory
  (`Calculations/<topic>/`, like `Calculations/wannier/`). Add the page to `latex/include.list`
  next to its nav neighbours.
- Model pages are in `docs/documentation/Calculations/`. Links: same directory `dos.md`, other
  kinds `../Reference/<feature>Setup.md`, `../Theory/<topic>.md`.
- `##` page title, `###` for steps and sections. No front matter.
- File names, tags and switches in backticks; XPath notation as in `userGuideOverview.md`.
- `!!! note ""` for text that tells the user to change the input; `!!! warning "..."` for traps.
- Images: PNG in `Calculations/` next to the page, with a name unique across
  `docs/documentation/` (prefix it with the feature), `![cap=Caption.,width=0.8\textwidth](file.png)`.
  Plots must be made from real output.
- No inline HTML (do not copy the `<div>` style of `greensFunctionTutorial.md`).

## Skeleton

```markdown
## <Calculating X>

<What X is in one or two sentences, what the calculation gives, prerequisites.
Background: [<theory page title>](../Theory/<topic>.md).>

### Setting up the input

!!! note ""
    Set `<xpath>/@<attr>` to `T` ...

```xml
<element attr="..."/>
```

### Running the calculation

<Commands, number of runs, what happens in each.>

### Output

<Files / out.xml tags / HDF datasets and how to read them, e.g.:>

```bash
grep '<tag' out.xml
```

### Example: <system>

<Input taken from testing/inputfiles/..., expected values from a real run.>

### Hints

<Convergence, pitfalls, cost.>

!!! note "More info"
    ...
```
