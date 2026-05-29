# wannierlib in FLEUR

This page documents the new `wannierlib` input path in FLEUR and summarizes
how projection descriptors are interpreted.

The behavior described here is based on:
- `src/libraries/fleurinput/types_wannierlib.f90`
- `src/fleur/io/xml/FleurInputSchema.xsd`
- Wannier90 projection definitions:
  https://wannier90.readthedocs.io/en/latest/user_guide/wannier90/projections/#angular-functions

## 1. Where to put input in `inp.xml`

`wannierlib` is configured in two places:

1. Global workflow controls in:
   `/fleurInput/output/wannierlib`
2. Projection descriptors per species in:
   `/fleurInput/atomSpecies/species/wannierproj`

## 2. Global `wannierlib` block

### 2.1 Structure

```xml
<output>
  <wannierlib wannierize="T">
    <bands numBands="20" minBand="21" maxBand="40"/>
    <disentanglement
      disWinMin="-8.0" disWinMax="12.0"
      disFrozMin="-6.0" disFrozMax="4.0"
      numIter="200"
      mixRatio="0.5"
      convTol="1.0e-10"/>
  </wannierlib>
</output>
```

### 2.2 Attributes and meaning

#### `/fleurInput/output/wannierlib/@wannierize`
- Type: boolean (`F`/`T`)
- Schema default: `F`
- Effect in parser: if the node exists, internal flag is set true and then possibly
  overwritten by this attribute if present.

#### `/fleurInput/output/wannierlib/bands`
- `numBands` (positive integer, schema default `0`)
- `minBand` (positive integer, schema default `0`)
- `maxBand` (positive integer, schema default `0`)

Runtime consistency logic in `init_wannierlib`:
- If `minBand == 0`, it is auto-set to:
  - `atoms%nlotot + 1` (collinear), or
  - `2 * atoms%nlotot + 1` (noncollinear or SOC)
- If `maxBand == 0` and `numBands > 0`, then `maxBand = minBand + numBands - 1`
- If `numBands == 0` and both `minBand` and `maxBand` are set, then
  `numBands = maxBand - minBand + 1`
- If all three are set and inconsistent, FLEUR aborts with an error.

#### `/fleurInput/output/wannierlib/disentanglement`
All attributes are required by schema:
- `disWinMin`
- `disWinMax`
- `disFrozMin`
- `disFrozMax`
- `numIter`
- `mixRatio`
- `convTol`

These are passed to Wannier90 library-mode options as disentanglement controls.

## 3. Species projections: `<wannierproj .../>`

Each projection is attached to one species entry:

```xml
<atomSpecies>
  <species name="Fe-1" element="Fe" atomicNumber="26">
    ...
    <wannierproj l="2" m="0" spin="u" rwf="1" zona="1.0"/>
    <wannierproj l="1" m="3" spin="d" theta="0.0" phi="0.0"/>
  </species>
</atomSpecies>
```

### 3.1 XML attributes (schema + parser defaults)

Required:
- `l` (integer)

Optional:
- `m` (integer)
- `spin` (string)
- `theta`, `phi` (aliases for `beta`, `alpha` in parser)
- `alpha`, `beta`, `gamma`
- `rwf` (integer)
- `zona`, `regio`
- `j`, `mj`
- `weight`
- `shiftX`, `shiftY`, `shiftZ`

Internal defaults when optional attributes are missing (from parser):
- `m = 0` 
- `rwf = 0`
- `alpha = 0.0`, `beta = 0.0`, `gamma = 0.0`
- `zona = 0.0`
- `regio = 1.0`
- `j = -1.0`, `mj = 0.0`
- `weight = 1.0`
- `shift = (0.0, 0.0, 0.0)`
- `spin`: interpreted as
  - starts with `u`/`U` -> `+1`
  - starts with `d`/`D` -> `-1`
  - otherwise -> `0` (auto/unspecified)

Note on angle aliases:
- `alpha` can be given as `phi`
- `beta` can be given as `theta`

## 4. Projection expansion rules in FLEUR

After XML read-in, projections are expanded in `init_wannierlib`:

1. Species expansion:
- A projection declared on one species is replicated to all atoms of that species.

2. `m` expansion:
- If `m == 0`, FLEUR expands to all allowed `m` channels for the chosen `l`.
- If `m != 0`, only that specific channel is used.

3. Spin expansion:
- In noncollinear (`l_noco`) or SOC (`l_soc`) mode:
  - if spin is unspecified (`spin -> 0`), projection is duplicated to up and down.
- Otherwise spin is kept as specified.

4. Shift handling:
- The final per-atom shift is
  `proj_shift + atoms%pos(:, atom_index)`.
   !Attention! As the shift is on a 'per species' basis, this only makes sense if you do not have symmetry equivalent atoms with the same species.

The resulting expanded number of Wannier projections is stored as `num_wann`.

## 5. Allowed `l` values and number of channels

FLEUR uses the following `l` -> number of channels mapping for expansion (`m=0`):

- `l = 0` -> 1 channel
- `l = 1` -> 3 channels
- `l = 2` -> 5 channels
- `l = 3` -> 7 channels
- `l = -1` -> 2 channels  (sp)
- `l = -2` -> 3 channels  (sp2)
- `l = -3` -> 4 channels  (sp3)
- `l = -4` -> 5 channels  (sp3d)
- `l = -5` -> 6 channels  (sp3d2)

Any other value is rejected when expansion is needed.

## 6. Projection names from Wannier90 angular functions

Wannier90 defines real angular-function labels for `(l, mr)` that are useful when
choosing `l` and `m` in FLEUR.

### 6.1 `l >= 0` angular functions

- `l=0`: `m=1 -> s`
- `l=1`: `m=1 -> pz`, `m=2 -> px`, `m=3 -> py`
- `l=2`: `m=1 -> dz2`, `m=2 -> dxz`, `m=3 -> dyz`, `m=4 -> dx2-y2`, `m=5 -> dxy`
- `l=3`:
  - `m=1 -> fz3`
  - `m=2 -> fxz2`
  - `m=3 -> fyz2`
  - `m=4 -> fz(x2-y2)`
  - `m=5 -> fxyz`
  - `m=6 -> fx(x2-3y2)`
  - `m=7 -> fy(3x2-y2)`

### 6.2 Hybrid sets (`l < 0`)

- `l=-1` (`sp`): `m=1..2`
- `l=-2` (`sp2`): `m=1..3`
- `l=-3` (`sp3`): `m=1..4`
- `l=-4` (`sp3d`): `m=1..5`
- `l=-5` (`sp3d2`): `m=1..6`

This matches the internal FLEUR expansion counts used in `projection_m_count`.

For full mathematical expressions of each angular function and hybrid combination,
see Wannier90 documentation:
https://wannier90.readthedocs.io/en/latest/user_guide/wannier90/projections/#angular-functions

## 7. Complete minimal example

```xml
<fleurInput>
  ...
  <atomSpecies>
    <species name="Fe-1" element="Fe" atomicNumber="26">
      ...
      <!-- all d channels on all Fe-1 atoms, spin auto-expands in SOC/NOCO -->
      <wannierproj l="2" m="0" spin="" rwf="1" zona="1.2"/>
      <!-- one specific p channel with explicit spin -->
      <wannierproj l="1" m="3" spin="u" theta="0.0" phi="0.0"/>
    </species>
  </atomSpecies>

  <output>
    <wannierlib wannierize="T">
      <bands numBands="20" minBand="21" maxBand="40"/>
      <disentanglement
        disWinMin="-8.0" disWinMax="12.0"
        disFrozMin="-6.0" disFrozMax="4.0"
        numIter="200" mixRatio="0.5" convTol="1.0e-10"/>
    </wannierlib>
  </output>
  ...
</fleurInput>
```

## 8. Practical notes

- If you use `m=0`, verify your `l` is one of the supported values above.
- Keep `numBands`, `minBand`, `maxBand` consistent if all are set explicitly.
- currently in `disentanglement` all values are mandatory by schema for the `wannierlib` block.
