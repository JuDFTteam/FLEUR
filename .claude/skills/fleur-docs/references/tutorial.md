# (d) Tutorial — interactive Jupyter notebook

Shows the feature on a small system that a learner runs in the FLEUR tutorial container
(`judft/future` image: `fleur_MPI`, `inpgen`, masci-tools, bash kernel, `OMP_NUM_THREADS=2`, a
few cores). Aim for at most a few minutes per FLEUR cell. Leave it out when no small system
shows the feature, or when it only repeats an existing tutorial (extend that one instead).

## Placement in `fleur_tutorial_container` (branch `master`)

- Chapter: `1.BasicFLEUR` (codes F1…), `2.Magnetism` (MA1…), `3.AdvancedTopics` (A1…).
  Take the next free code: `ls tutorial/<chapter>/`.
- `tutorial/<chapter>/<code>/<CamelCaseName>.ipynb`
- Inputs in a subfolder named after the system: `<code>/<System>/inp_<System>.txt`
- Helper scripts and PNGs (expected results) in `<code>/`.
- Registration: one line `***<code>***: [<Title>](<code>/<Name>.ipynb)` in the chapter's
  `Outline.ipynb`. `www.flapw.de/tutorial_to_webpage.pl` parses exactly this pattern, so do not
  change it. `new_notebook.py --register` writes it.

## Create the skeleton

```bash
SK=<fleur>/.claude/skills/fleur-docs/scripts
python3 $SK/new_notebook.py tutorial/3.AdvancedTopics/A9/<Name>.ipynb --title "<Title>" \
    --system <System> --input inp_<System>.txt \
    --objective "<what the reader learns>" --objective "..." \
    --docs Calculations/<howtoPage> --docs Reference/<referencePage> --register
```

Use `--kernel sos` only if Python cells (plots with masci-tools/matplotlib) are needed. Then fill
in every `TODO`. `check_notebook.py` warns about each one left.

## Content rules

- **Sections**: `# Title`, `## Objectives` (bullets), `## Introduction`, the hands-on steps
  (`##`/`###` headings), `## Learn more` linking the howto, reference and theory pages as
  `https://www.flapw.de/rel/documentation/<Section>/<page>/`, where `<Section>` is `Calculations`,
  `Reference` or `Theory` (the directory of the page under `docs/documentation/`). These URLs
  follow the new site layout and work under `/develop/` now, and under `/rel/` once the current
  release is built from it. Older tutorials still link the flat `/rel/documentation/<page>/`;
  do not copy that form.
- **Cells**: plain shell commands (bash kernel, no `!` or `%%bash`), one action per code cell. A
  markdown cell before each code cell says what it does and why. The cell after it says what to
  look at in the output.
- **Admonitions** (MyST, converted for the website by `convert_admonition.pl`):
  ```` ```{attention} ```` for input changes, ```` ```{warning} ```` for runtime or reduced
  parameters, ```` ```{hint} ```` and ```` ```{note} ````. Do not nest fenced blocks inside an
  admonition unless the outer fence is longer (four backticks).
- **Input changes must be scripted**, so that the notebook runs unattended (and
  `check_notebook.py --execute` can verify it):
  - Attribute values: `fleur_MPI -xmlXPath "/fleurInput/<path>/@<attr>=<value>"`. Several
    settings are separated by `:`. This only affects that run, the element must already exist,
    and path and value are each limited to 100 characters.
  - Persistent changes or new elements: a visible `sed -i` (or heredoc) in the cell, followed by
    `grep` showing the result.
  - In an ```` ```{attention} ```` block, also say which XPaths change, so readers learn the input
    and could edit `inp.xml` by hand.
- **Small**: one or two atoms, a coarse k-mesh via `&kpt div1=...` in the inpgen input or
  `inpgen -inp.xml -kpt grid=n,n,n`, few iterations, or a shipped `cdn.hdf` to skip the SCF.
  Mention in a ```` ```{warning} ```` that parameters are reduced.
- **Results**: inspect with `grep`, `xmllint --xpath`, `masci-tools plot ...` or a helper script.
  State the value the reader should see: the number from your own execution, not an invented one.
- **Outputs cleared** in the committed notebook (`check_notebook.py --clear`).
- AiiDA-FLEUR parts are optional, in a separate `## AiiDA-FLEUR workflow` section. They cannot be
  executed by `check_notebook.py` here.

## Verify

```bash
python3 $SK/check_notebook.py <nb> --www <www-worktree>        # structure, links, registration
python3 $SK/check_notebook.py <nb> --execute --build-dir <fleur>/build.debug \
        --python <scratch>/nbenv/bin/python --workdir <scratch>
```

A working Jupyter with `bash_kernel` is needed. If the system `jupyter` is broken, make a scratch venv:

```bash
python3.13 -m venv <scratch>/nbenv
<scratch>/nbenv/bin/pip install -q nbconvert nbclient bash_kernel
<scratch>/nbenv/bin/python -m bash_kernel.install --sys-prefix
```

`--execute` runs a copy. `fleur_MPI`, `fleur` and `inpgen` are mapped to the build directory, the
stack limit is raised, and personal Jupyter config is ignored. The executed copy is kept for
inspection; take the expected values in the text from it. Exit code 2 means execution was skipped.
Say so to the user rather than claiming the notebook was tested.
