---
name: fleur-docs
description: Write the user documentation for a FLEUR feature - up to four separate documents - and deliver them as merge requests. (a) FLEUR-specific theory, (b) a howto with steps, input and expected output, and (c) a concise input/output reference go to the www.flapw.de repo (MR into develop); (d) an interactive Jupyter tutorial goes to fleur_tutorial_container (MR into master). Use whenever a change is user-visible (new or changed inp.xml input, new output in out.xml, files or HDF5, new command-line option) and before such work is called finished; when asked to document a feature or write a tutorial; or to turn a draft in the fleur repo's docs/ into proper site pages.
---

# FLEUR user documentation

A feature is not finished until users can find out how to use it. This skill writes up to four
documents, each with one job, and puts them where FLEUR users look:

| Kind | Job | Where |
|---|---|---|
| a) Theory | FLEUR-specific background only: representation, conventions, approximations | www.flapw.de `docs/documentation/Theory/<topic>.md`, nav "Theoretical background" |
| b) Howto | steps of the calculation, required input, expected output, example | www.flapw.de `docs/documentation/Calculations/<feature>.md`, nav "Basic calculations" or "Expert knowledge and troubleshooting" |
| c) Reference | terse tables: every input attribute, every output file/tag/dataset | www.flapw.de `docs/documentation/Reference/<feature>Setup.md`, nav "Reference" |
| d) Tutorial | interactive notebook on a small system | fleur_tutorial_container `tutorial/<chapter>/<code>/<Name>.ipynb` |

The site layout under `docs/` on `develop`:

| Directory | Holds | Use for feature docs? |
|---|---|---|
| `documentation/Theory/` | theory pages + their PNGs, `theoryOverview.md` | yes, kind a |
| `documentation/Calculations/` | howto pages for both nav sections "Basic calculations" and "Expert knowledge" + PNGs, `calculationsOverview.md`; large topics get a subdirectory (`Calculations/wannier/`) | yes, kind b |
| `documentation/Reference/` | reference pages, `referenceOverview.md`, `references.md` (literature) | yes, kind c |
| `documentation/` (top level) | overview, quick start, installation, running, `expertKnowledge.md`, `schema.md` | no |
| `home/`, `examples/`, `unused/`, `index.md` | home pages, old tutorials, retired pages | no |

Detailed rules and a skeleton for each kind are in `references/theory.md`, `references/howto.md`,
`references/reference.md` and `references/tutorial.md`. Read the matching file before writing
that document.

Scripts (Python standard library only; if `python3` cannot parse XML because of a broken
pyexpat, use `/usr/bin/python3` or `python3.13`):

```bash
SK=<fleur>/.claude/skills/fleur-docs/scripts
python3 $SK/schema_diff.py diff <base> [<head>]   # input/output schema changes between git refs
python3 $SK/schema_diff.py show <xpath>           # one element as xml example + attribute table
python3 $SK/schema_diff.py find <name>            # all XPaths of an element/attribute name
python3 $SK/check_pages.py <www-worktree>         # nav, links, fences, math, xml/XPath vs schema
python3 $SK/new_notebook.py <nb> --title ... --register   # tutorial skeleton + Outline entry
python3 $SK/check_notebook.py <nb> [--execute ...]        # structure, registration, execution
```

## 1. Scope: what changed for the user, and which documents it needs

Find the user-visible changes of the work (base = where the feature branch left `develop`, or
the commit before the feature):

```bash
python3 $SK/schema_diff.py diff <base>                      # inp.xml and out.xml schema changes
git diff <base> --stat -- src/libraries/fleurinput/ src/fleur/io/ src/fleur/main/fleur_help.F90
git log --oneline <base>..HEAD
```

Also look for new files or HDF5 datasets written by the feature (`h5dcreate_f`, `io_write_*`,
`OPEN(` in its modules) and new tests in `testing/tests/` with inputs in `testing/inputfiles/`.

Decide per kind:

- **a) Theory**: only if FLEUR has its own formalism, representation, approximation or convention
  here that a reader needs to understand the input or the result. Not for a new output switch of
  a known quantity.
- **b) Howto**: whenever a user performs a calculation with the feature.
- **c) Reference**: whenever input or output changed. That is almost always.
- **d) Tutorial**: when a small system shows the feature in minutes in the tutorial container.

Tell the user the selection, one line of reason per kind, and the planned page names and nav
positions. Wait for confirmation before writing.

## 2. Inventory what exists

Look before writing, so you extend pages instead of duplicating them:

```bash
grep -rli '<feature words>\|<new xml names>' <www>/docs/documentation/ <www>/docs/unused/ <tutorials>/tutorial/
ls <fleur>/docs/                                   # drafts, e.g. wannierlib.md, densityMatrixDOS.md
```

- Prefer extending an existing page (e.g. a new DOS mode → `dos-modes.md`, Green's-function
  options → `greensFunctionSetup.md`). Keep that page's format. A retired page in `docs/unused/`
  (e.g. `hybrids.md`, `exchangeMatrix.md`) may be a starting point. Moving it back into
  `documentation/` is a change the user should agree to.
- Drafts in the fleur repo's `docs/` (and the feature's commit messages) are source material. They
  usually mix theory, howto and reference, so split their content across a–c. Do not delete
  or move the drafts without asking.

## 3. Gather facts from the code, not from memory

Each statement in the documents must trace back to one of these sources:

| Fact | Source |
|---|---|
| element/attribute names, types, schema defaults | `schema_diff.py show` on `FleurInputSchema.xsd` |
| real defaults and accepted values | `read_xml` in `src/libraries/fleurinput/types_*.F90` (or the reading module) |
| out.xml tags | `FleurOutputSchema.xsd` (`--schema output`) and the `writeXMLElement*` calls |
| files and HDF5 layout | the writer code: dataset paths, shapes, index order, units |
| example input | `testing/inputfiles/<dir>/<test>/` (name the test) |
| expected numbers | the test's reference values in `testing/tests/parameterized/test_*.py`, or a run you made |
| physics intent, references | code comments, the commit message, the developer. Ask the user when unclear |

Never invent attribute names, defaults, numbers, plots or DOIs. If the schema and the code
disagree on a default, tell the user. The code wins in the docs, and the mismatch is a
`fleur-bugfix` case.

## 4. Set up worktrees in the two repos

Leave the user's clones and their working trees alone (they may have local work or an unfinished
merge). Work in worktrees on a fresh branch:

```bash
WWW=${FLEUR_WWW_REPO:-<fleur>/../www.flapw.de}
TUT=${FLEUR_TUTORIAL_REPO:-<fleur>/../fleur_tutorial_container}
git -C $WWW fetch origin develop
git -C $WWW worktree add -b docs/<feature> $WWW.worktrees/docs-<feature> origin/develop
git -C $TUT fetch origin master
git -C $TUT worktree add -b docs/<feature> $TUT.worktrees/docs-<feature> origin/master
```

If a clone is missing, clone it next to the fleur checkout
(`glab repo clone fleur/www.flapw.de`, `glab repo clone fleur/fleur_tutorial_container`).
Only create the worktrees the selection needs.

Facts about the target repos:

- **www.flapw.de** (`fleur/www.flapw.de`, project id 1018). `develop` holds the complete site:
  `docs/`, `mkdocs.yml`, the theme (`web_scheme/`), the schema browser (`fleurinput-xml/`) and the
  PDF manual (`latex/`). CI builds it as is and publishes it at `https://www.flapw.de/develop/`.
  Release branches (`MaX-x.y`) are published at `https://www.flapw.de/<branch>/`, and `rel` points
  to the current release. Feature pages go to `develop` only.
- **PDF manual**: `latex/convert.pl` includes only the pages listed in `latex/include.list`
  (paths relative to `docs/`, in reading order). It copies all PNGs of
  `docs/documentation/*/` into **one flat directory**, so an image name must be unique across
  `Theory/`, `Calculations/` and `Reference/`.
- **fleur_tutorial_container** (`fleur/fleur_tutorial_container`, project id 1102). `master` is
  the maintained branch; CI publishes it as the `judft/future:latest` image.
- `www.flapw.de/tutorial_to_webpage.pl` converts each chapter's `Outline.ipynb` entries to site
  pages. That is why the Outline line format matters.

## 5. Write the documents

Follow the matching `references/<kind>.md`. Cross-link the documents:

- theory ↔ howto
- howto → reference ("More info" note) and → tutorial
- reference → howto
- tutorial "Learn more" → `https://www.flapw.de/rel/documentation/<Section>/<page>/`, e.g.
  `.../documentation/Calculations/dos/`

Links between site pages are relative to the page: within a directory just `ldaUSetup.md`, across
directories `../Reference/ldaUSetup.md`, `../Theory/nocoMag.md`, `../Calculations/dos.md`. Anchors
use the heading slug: `../Reference/greensFunctionSetup.md#energy-contour-setup`.

Rules for all pages:

- **Site pages:**
  - Put each page in the directory for its kind (table above). Add it to the `nav` in the
    worktree's `mkdocs.yml` under the matching section, and to `latex/include.list` at the matching
    position (after the neighbouring page of the same section) so it reaches the PDF manual.
  - Use XPath notation and `!!! note ""` blocks for input, as
    `documentation/userGuideOverview.md` describes.
  - No inline HTML.
  - Put PNGs in the same directory as the page that uses them, made from real output. Give them a
    name that is unique across all of `docs/documentation/`, e.g. prefixed with the feature name.
- **Tutorial:** create the skeleton with `new_notebook.py ... --register`, script all input
  changes, and keep the system small (see `references/tutorial.md`).
- **Readers:** write for FLEUR users, not developers. Leave out source file names, internal
  variable names and the history of the code.

## 6. Verify

```bash
python3 $SK/check_pages.py $WWW.worktrees/docs-<feature>          # all changed pages vs origin/develop
python3 $SK/check_notebook.py <nb> --www $WWW.worktrees/docs-<feature>
python3 $SK/check_notebook.py <nb> --execute --build-dir <fleur>/build.debug \
        --python <scratch>/nbenv/bin/python --workdir <scratch>
```

- Fix all errors from `check_pages.py`. Treat its schema warnings (unknown element, attribute or
  XPath) as errors unless you can show the name is valid.
- Run the notebook against a current build (see the `fleur-build-test` skill for which build
  directory). Put the values from the executed copy into the text, then commit with outputs
  cleared. If execution was skipped (exit 2), tell the user. Do not claim the notebook was tested.
- Optional but recommended for site pages: build the site the way CI does, directly in the
  worktree (the output goes to the scratchpad):

  ```bash
  python3.13 -m venv <scratch>/mkenv
  <scratch>/mkenv/bin/pip install -q mkdocs mkdocs-bootswatch python-markdown-math markdown-include
  cd $WWW.worktrees/docs-<feature>
  <scratch>/mkenv/bin/mkdocs build -d <scratch>/site 2>&1 | grep -E 'WARNING|ERROR' | grep '<page names>'
  ```

  The site has many old warnings (mostly in `examples/` and `unused/`), so filter for the new
  pages. Open `<scratch>/site/documentation/<Section>/<page>/index.html` to check the math and
  tables.

## 7. Deliver: commit, then push and open MRs only after the user says so

Commit in each worktree. Use an imperative subject that names the feature, and a body that lists
the pages and the fleur commit or MR they document. Keep it under 300 words, with no AI co-author
line.

```bash
git -C $WWW.worktrees/docs-<feature> add -A docs/documentation mkdocs.yml latex/include.list && git -C ... commit
git -C $TUT.worktrees/docs-<feature> add -A tutorial && git -C ... commit
```

Show the user the file list (`git show --stat`) and the selection of documents. **Pushing and
opening MRs is visible to all FLEUR developers, so ask first**, every time:

```bash
git -C $WWW.worktrees/docs-<feature> push -u origin docs/<feature>
glab mr create -R iffgit.fz-juelich.de/fleur/www.flapw.de --source-branch docs/<feature> \
    --target-branch develop --title "Document <feature>" --description "$(cat <body-file>)"
git -C $TUT.worktrees/docs-<feature> push -u origin docs/<feature>
glab mr create -R iffgit.fz-juelich.de/fleur/fleur_tutorial_container --source-branch docs/<feature> \
    --target-branch master --title "Tutorial <code>: <title>" --description "$(cat <body-file>)"
```

The www.flapw.de MR description lists:

- the new and changed pages, their nav positions and `latex/include.list` entries
- the fleur commit or MR being documented
- a link to the tutorial MR

The tutorial MR description notes that its "Learn more" links use the new site layout
(`documentation/<Section>/<page>/`) and resolve under `/rel/` once the release in `rel` is built
from the new structure. Until then they work under `https://www.flapw.de/develop/`.

Report the MR URLs. Keep the worktrees until the MRs are merged, since review comments are
addressed there. Remove them when the user says the MRs are merged:
`git -C $WWW worktree remove $WWW.worktrees/docs-<feature>`.

## Checklist

- [ ] User-visible changes listed (`schema_diff.py diff`, writers, CLI flags)
- [ ] Selection of a–d agreed with the user, with page names and nav positions
- [ ] Existing pages and drafts checked; existing pages extended where they fit
- [ ] Every name, default, unit and number traced to schema, parser, writer, test or a real run
- [ ] Worktrees on `docs/<feature>` from `origin/develop` (site) and `origin/master` (tutorial)
- [ ] Pages in `Theory/`, `Calculations/` or `Reference/`; in `mkdocs.yml` nav and `latex/include.list`; cross-linked
- [ ] Image names unique across `docs/documentation/`; `check_pages.py` without errors or schema warnings
- [ ] Notebook registered in `Outline.ipynb`; `check_notebook.py` clean; executed (or skip reported); outputs cleared
- [ ] Commits made; pushed and MRs opened only after the user agreed
- [ ] MR description lists pages, nav and `include.list` entries, and the documented fleur change
- [ ] User told: documents written, checks run (and skipped), MR URLs
