# wannierlib in FLEUR

`wannierlib` runs Wannier90 as a **library** inside FLEUR. There is no `seedname.win` and no
separate `wannier90.x` call: FLEUR builds the overlaps and projections, hands them to
Wannier90 through its module API, and post-processes the resulting gauge in the same
execution. Everything is driven from `inp.xml`.

Companion documents:
- `src/fleur/wannierlib/README.md` — how the pipeline is put together.
- `src/fleur/matrixelements/README.md` — how to add an operator.

## 0. Before you start

**Wannier90 3.x, built as a library.** Library mode needs the module API, which arrived in
Wannier90 3.x. Against 1.2 FLEUR still compiles, but the interface is simply absent and
nothing here works.

```bash
cd wannier90-3.1.0
cp config/make.inc.ifort make.inc      # then set F90 = mpiifort, COMMS = mpi
make default lib                       # -> libwannier.a

./configure.sh -libxc -wannier -hdf5 \
     -libdir <wannier90 dir> -includedir <hdf5 dir>/include -l wannlib
cd build.wannlib && make -j16
```

The classic route (`<wannier>`, `CPP_WANN`) and library mode (`<wannierlib>`,
`CPP_WANNLIB_API`) use different Wannier90 interfaces, so one binary serves one route. Keep
a build directory per route; `configure.sh -l <name>` writes `build.<name>`.

**The wannierisation mesh must cover the full zone**, gamma-centred and with symmetry off,
because the neighbour shell **b** is otherwise incomplete:

```bash
inpgen -inp.xml -kpt "wann8#gamma@grid=8,8,8" -noKsym
```

Then point `<kPointListSelection listName="wann8"/>` at that list. The SCF itself may keep
using the symmetry-reduced mesh.

**The number of k-points must divide the number of MPI ranks exactly.** If it does not,
FLEUR does *not* abort — it **hangs** on one k-point.

**Numbers in `inp.xml` are plain decimals.** FLEUR's expression parser reads digits, one
decimal point and a leading sign, and stops at anything else, so `1.0e-10` fails with
`Error in expression: Unknown character string found`. Write `0.0000000001`. (`10^-10` fails
too: the parser rejects an operator following an operator.)

**Working examples** live in `testing/inputfiles/wannier/` — fifteen cases covering
collinear, non-collinear, SOC, interpolation and the operator exports. Starting from one of
those and changing the structure is the shortest path in.

## 1. Where to put input in `inp.xml`

`wannierlib` is configured in two places:

1. Global workflow controls in `/fleurInput/output/wannierlib`
2. Projection descriptors per species in `/fleurInput/atomSpecies/species/wannierproj`

## 2. The `<wannierlib>` block

### 2.1 Structure

Only `<bands>` is required; every other sub-block is optional and independent.

```xml
<output>
  <wannierlib wannierize="T">
    <bands numBands="36" minBand="7" maxBand="42"/>

    <disentanglement disWinMin="0.097" disWinMax="2.670"
                     disFrozMin="0.097" disFrozMax="0.661"
                     numIter="3000" convTol="0.00001" mixRatio="0.5"/>

    <wannierization numIter="3000" convTol="0.00001"/>

    <operators_r>
      <operator name="hamiltonian"/>
      <operator name="spin"/>
    </operators_r>

    <interpolation useWsDistance="F">
      <domain listName="path-2"/>
      <operator name="hamiltonian"/>
    </interpolation>

    <export wannier90="T"/>
  </wannierlib>
</output>
```

### 2.2 `<wannierlib>` itself

| attribute | default | meaning |
|---|---|---|
| `wannierize` | `F` | run the wannierisation |
| `plotWF` | `F` | write the Wannier functions as XSF files, one per function, on a 40x40x40 grid over one unit cell. Costs a second pass over the k-points. Not available for spinor Wannier functions (noco or SOC). |

### 2.3 `<bands>` — which states enter

| attribute | default | |
|---|---|---|
| `numBands` | `0` | how many states are wannierised |
| `minBand` | `0` | first state |
| `maxBand` | `0` | last state |

Resolved in `init_wannierlib`:

- `minBand == 0` becomes `atoms%nlotot + 1` (collinear) or `2*atoms%nlotot + 1` (noco or SOC).
- `maxBand == 0` with `numBands > 0` gives `maxBand = minBand + numBands - 1`.
- `numBands == 0` with both bounds set gives `numBands = maxBand - minBand + 1`.
- All three set and inconsistent: FLEUR aborts.

`minBand` is worth measuring rather than guessing: it is the first band whose maximum over
the whole mesh enters the window, which is what excludes the semicore.

### 2.4 `<disentanglement>` — only when `numBands > num_wann`

| attribute | required | default | |
|---|---|---|---|
| `disFrozMax` | **yes** | | top of the frozen window |
| `numIter` | **yes** | | disentanglement iterations |
| `convTol` | **yes** | | convergence threshold |
| `mixRatio` | **yes** | | mixing |
| `disWinMin` | no | derived | bottom of the outer window |
| `disWinMax` | no | derived | top of the outer window |
| `disFrozMin` | no | `disWinMin` | bottom of the frozen window |
| `disFrozProj` | no | `F` | freeze by projectability instead of by energy |
| `disProjMin` | no | `0.01` | |
| `disProjMax` | no | `0.95` | |
| `spinBalanced` | no | `F` | |

All energies are in **Hartree**, like the rest of `inp.xml`.

> **State `disWinMin` and `disWinMax` explicitly for now.** They are optional in the schema
> and are meant to be derived from the band range on the wannierisation mesh, but that
> derivation currently cannot see the eigenvalues and the run stops with
> `the outer energy window cannot be derived from the bands`, printing what it did see into
> the `out` file. Giving them is the supported path.

Measure the outer window **on the wannierisation mesh**, not on the SCF one. The SCF mesh
need not contain Gamma, so its minimum comes out above the true one; states then fall
outside the window at some k and Wannier90 drops them without a word.

### 2.5 `<wannierization>` — the spread minimisation

| attribute | required | |
|---|---|---|
| `numIter` | **yes** | MLWF iterations |
| `convTol` | **yes** | convergence threshold |

`numIter="0"` skips the minimisation and keeps the projection gauge: the Wannier functions
are then the trial orbitals themselves, orthonormalised. That is a real choice, not a
degenerate case — the minimisation converges to maximally localised functions whatever it
starts from, so the initial projections only decide the gauge when it is switched off.

### 2.6 `<operators_r>` — O(R) for external post-processing

Writes the real-space operator matrices in the Wannier90/FLEUR standalone format. No band
interpolation: these are the files an external transport code reads.

```xml
<operators_r>
  <operator name="hamiltonian"/>   <!-- WF<n>_hr.dat   -->
  <operator name="position"/>      <!-- WF<n>_r.dat    -->
  <operator name="spin"/>          <!-- rspauli.1      -->
  <operator name="orbital"/>       <!-- anglmomrs.1    -->
  <operator name="spin_orbit"/>    <!-- rssocmat.1     -->
</operators_r>
```

Accepted names:

| name | what it is |
|---|---|
| `hamiltonian` | H(R) |
| `position` | A(R) = `<0n|r|Rm>` |
| `position_pw90` | the same in the postw90 convention (Eq. 44 of WYSV06 on the diagonal as well, hermitised). Both are needed at once: `position` carries the Wannier centres, while `berry.F90` refuses that form for the orbital magnetisation. Writes `WF<n>_rpw.dat`. |
| `bmn` | B(R) = `<0n|H(r-R)|Rm>` |
| `fmn` | F(R) = `<0n|r_a r_b|Rm>` — needs the pair overlap (uIu) |
| `cmn` | C(R) = `<0n|r_a H r_b|Rm>` — needs uHu |
| `spin` | S(R) |
| `orbital` | L(R) |
| `spin_orbit` | the spin-orbit operator |

`fmn` and `cmn` are refused on a film: they need the vacuum halves of the pair overlap and
of the momentum, which this layer does not reach. `spin` is refused on a film for the same
kind of reason — it sums muffin tins and interstitial, and a film has a third region
carrying spin density.

### 2.7 `<interpolation>` — bands and operators on any k-list

```xml
<interpolation useWsDistance="F">
  <domain listName="path-2"/>
  <domain listName="path-2" suffix="dense" npts="10"/>
  <operator name="hamiltonian"/>
</interpolation>
```

`<domain>` — one per set of output k-points, repeatable:

| attribute | required | |
|---|---|---|
| `listName` | **yes** | an existing kPointList, in `kpts.xml` or any included file |
| `suffix` | no | distinguishes this domain's `bands_wann_*.dat`; may be left off for at most one domain |
| `npts` | no | subdivides each segment of an ordered list into `npts` pieces (`npts<=1` = list as-is) |

Whether a list is a line, a plane or a mesh is a property of the list, not of FLEUR.

`useWsDistance` (default `F`) picks the minimum-distance replica of each R vector. It
matters when the Wannier centres sit far from the origin. **Wannier90's own default is
`T`**; the default here is deliberately the other one.

`<operator>` — what to interpolate, repeatable:

| name | output |
|---|---|
| `hamiltonian` | `bands_wann_interpol.dat`, `bands_wann_interpol_ev.dat` |
| `spin` | `bands_wann_spin.dat` |
| `orbital` | `bands_wann_orbmom.dat` |
| `spin_orbit` | `bands_wann_spin_orbit.dat` |
| `velocity` | `bands_wann_velocity.dat`, `bands_wann_berrycurv.dat` |
| `eigenstates` | `bands_wann_eigenstates.dat` |

`total="T"` (the default) asks for the site-summed projection; only `spin` and `orbital`
have anything to sum over. `comp` selects Cartesian components.

> The operator names are the **same** in `<operators_r>` and `<interpolation>` — in
> particular it is `spin_orbit` in both. The two blocks do not accept the same *set* of
> names, because not everything that has an O(R) has an interpolation driver and the other
> way round, but a name never means two different things.

### 2.8 `<export>` — hand out matrices the run already holds

One boolean per artefact, all `F` by default, all independent: they are not alternatives,
and asking for one does not turn the others off.

| attribute | writes |
|---|---|
| `wannier90` | `WF<n>.amn`, `.mmn`, `.eig` — what a standalone `wannier90.x` needs to repeat the same wannierisation. The `.mmn` of a 512-point mesh with 36 bands is a few hundred megabytes of text, which is why it is opt-in. |
| `wannierberri` | `WF<n>_basis.hdf` — the plane-wave index list and the interstitial eigenvector coefficients, for symmetry analysis. Needs a full-zone mesh and a build with HDF5. |
| `gauge` | `WF<n>_gauge.hdf` — `u_opt` and `u_mlwf`, the two factors Wannier90 returned. |
| `blochOperators` | `WF<n>_s0.dat` — the spin operator in the Bloch basis, before any gauge. |

What `<operators_r>` and `<interpolation>` write is deliberately not here: those compute
something first and writing it is the last step, while `<export>` only copies out.

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
      <!-- m="0" expands every m, so s+p+d is 1+3+5 = 9;
           with spinors these double on their own, 9 -> 18 -->
      <wannierproj l="0" m="0" spin=""/>
      <wannierproj l="1" m="0" spin=""/>
      <wannierproj l="2" m="0" spin=""/>
    </species>
  </atomSpecies>

  <output>
    <wannierlib wannierize="T">
      <bands numBands="36" minBand="9" maxBand="44"/>
      <disentanglement
        disWinMin="0.020" disWinMax="3.000"
        disFrozMin="0.020" disFrozMax="0.562"
        numIter="3000" convTol="0.00001" mixRatio="0.5"/>
      <wannierization numIter="3000" convTol="0.00001"/>
      <operators_r>
        <operator name="hamiltonian"/>
        <operator name="spin"/>
      </operators_r>
    </wannierlib>
  </output>
  ...
</fleurInput>
```

Energies in Hartree; `convTol` as a plain decimal, never `1.0e-5`.

## 8. Practical notes

- **Give `disWinMin` and `disWinMax`**, and measure them on the wannierisation mesh. See 2.4.
- **Plain decimals only.** `1.0e-10` is not a number to FLEUR's parser.
- **`n_k` must divide the MPI ranks exactly**, or the run hangs instead of failing.
- **The wannierisation mesh is full-zone and gamma-centred** (`inpgen ... -noKsym`).
- With `jspins=2` and no SOC the two spin channels wannierise **separately**, and the whole
  pipeline runs once per channel; with noco or SOC the states are spinors and there is one
  pass over a basis of twice the size.
- If you use `m=0`, check that `l` is one of the values in section 5.
- Keep `numBands`, `minBand` and `maxBand` consistent if you set all three.
- `Omega_I` is not a quality measure. It says how localised the subspace is, not how
  faithful: a narrower frozen window can give a better `Omega_I` and worse bands.
- The Fermi level of the Wannier Hamiltonian is **not** the Fermi level of the SCF — up to
  1 eV apart in Cu and Pt, with no warning. Shift the axis when reading the output; do not
  re-run.

## 9. Where to start

| you want | copy |
|---|---|
| collinear, no SOC | `testing/inputfiles/wannier/WannFeFM` |
| collinear, two channels, SOC as an operator | `WannFeAFMColSOC` |
| spinors with SOC | `WannPtSOC` |
| band interpolation | `WannFeBccInterp` |
| the operator exports | `WannPtSOCOps`, `WannFeAFMSOCOps` |
| non-collinear, several sites | `WannMn3IrNoco` |

Each directory holds a complete `inp.xml` with its `kpts.xml` and `sym.xml`.
