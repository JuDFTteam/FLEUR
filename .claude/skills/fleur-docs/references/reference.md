# (c) Reference page — concise input/output lookup

For readers who already know what they want and look up a name, default or file. Complete and
terse: one sentence of introduction, then examples and tables. No step-by-step text and no
theory.

## Input

Start from the schema, then correct it against the code:

```bash
SK=.claude/skills/fleur-docs/scripts
python3 $SK/schema_diff.py show <element-or-path> --schema input   # xml example + attribute table
python3 $SK/schema_diff.py find <name>                             # where else the name occurs
```

Then check every row against the parser (`read_xml` in
`src/libraries/fleurinput/types_*.F90`, or the module that reads the element):

- **Default**: the schema default and the parser's default must agree. If they differ, the
  code is what users get. Tell the user, since this is a bug for `fleur-bugfix`.
- **Units**: Hartree, bohr, reciprocal-lattice units, degrees or radians. Say which.
- **Allowed values** and their meaning, including aliases the parser accepts.
- **Interplay**: "only used if `@x="T"`", "ignored for films", "requires `l_noco`".
- The xml example comes from a real test input (`testing/inputfiles/...`), not from the
  generated skeleton with `...` placeholders.

Table format for new pages:

```markdown
|Attribute|Type|Default|Description|
|---|---|---|---|
|sigma|real|0.015|Gaussian broadening in Htr.|
```

When extending an existing page, keep its table format (most use `|Attribute|Description|` with
the default written in the description).

## Output

There is no out.xml reference on the site, so each reference page documents the output of its
feature:

- **Files** written to the working directory (name, when written, format/columns).
- **out.xml**: element XPath, attributes, units. Check them with
  `schema_diff.py show <element> --schema output` and in the writer code
  (`writeXMLElement*`, `openXMLElement*` calls of the feature's module).
- **HDF5** (`banddos.hdf`, `greensf.hdf`, ...): dataset path, shape and index order, dtype,
  meaning, units. Take them from the writer (`h5dcreate_f`, `io_write_*` calls).
- **Command-line options** of `fleur`/`inpgen` that belong to the feature (`fleur -h`,
  `src/fleur/main/fleur_help.F90`).

## Format (match `greensFunctionSetup.md`, `ldaUSetup.md`, `dos-modes.md`)

- File: `docs/documentation/Reference/<feature>Setup.md`. Nav: `User Guide → Reference`, before
  the entries "The XML-Schema for the inp.xml file" and "References". Add it to
  `latex/include.list` before `documentation/Reference/references.md`. If the feature belongs to
  an existing reference page (e.g. a new DOS mode → `Reference/dos-modes.md`), extend that page.
- Model pages are in `docs/documentation/Reference/`. Links: same directory `ldaUSetup.md`, howto
  `../Calculations/<feature>.md`.
- `###` page title, `####` per element or output kind. No front matter.
- One ```` ```xml ```` example per element, then its table.
- Link back to the howto page in the intro sentence.

## Skeleton

```markdown
### <Feature> setup

`<parent XPath>/<element>` switches on <feature>; see [<howto title>](../Calculations/<feature>.md) for usage.

```xml
<element attr="value">
   <child attr="value"/>
</element>
```

|Attribute|Type|Default|Description|
|---|---|---|---|
|...|...|...|...|

#### `<child>`

|Attribute|Type|Default|Description|
|---|---|---|---|

#### Output

|File|Content|
|---|---|
|`<file>`|...|

|out.xml element|Content|
|---|---|
|`/fleurOutput/scfLoop/iteration/<tag>`|...|

|HDF5 dataset (`<file>.hdf`)|Shape|Content|
|---|---|---|
|`/<group>/<dataset>`|(...)|...|
```
