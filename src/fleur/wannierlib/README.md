# `wannierlib` — Wannier functions in library mode

FLEUR builds the overlaps and projections, calls Wannier90 as a **library** (no
`seedname.win`, no separate `wannier90.x` run), and post-processes the result in the same
execution. Everything is driven from `inp.xml`.

Companion file: `../matrixelements/README.md` — how to add an operator.

## The pipeline

| Stage | Where | Produces |
|---|---|---|
| 1. overlaps and projections | `wannierlib_mmnkb`, `wannierlib_amn` | `M_mn(k,b)`, `A_mn(k)` |
| 2. operator matrices on the coarse mesh | `../matrixelements/melem_coarse` | `O(k)` on the ab-initio states |
| 3. wannierisation | `wannierlib_w90_adapter` → Wannier90 | the gauge `u_opt`, `u_matrix` |
| 4. post-processing | `postproc/melem_run` | `O(R)` files and interpolated bands |

Stage 2 needs no gauge and runs before the wannierisation; stages 3 and 4 do. With
`jspins=2` and no SOC the two spin channels wannierise **separately**, so stages 1–4 run
once per channel.

## Turning it on

```xml
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
  <interpolation>
    <domain listName="path-2"/>
    <operator name="hamiltonian"/>
    <operator name="velocity"/>
  </interpolation>
</wannierlib>
```

`<disentanglement>` is optional and is left out when `numBands == num_wann`. The initial
guesses are `<wannierproj l=".." m=".." spin=".."/>` children of a `<species>`; their total
count over all atoms is `num_wann`.

The two operator blocks are independent and spell some names differently:
`<operators_r>` writes `O(R)` to disk and does no interpolation; `<interpolation>` projects
onto the interpolated bands. The accepted names of each live in
`fleurinput/types_melem_optable.f90`, and a name the tables do not have stops the run with
the accepted list in the message.

## Output domains

Inside `<interpolation>`, declare one `<domain>` per set of k-points you want the operators
evaluated on. There is no limit and no kinds: a domain is a **named `kPointList`** plus,
when there is more than one, a suffix that keeps their files apart.

```xml
<domain listName="path-2"/>                 <!-- bands_wann_*.dat        -->
<domain listName="kz0"    suffix="_plane"/> <!-- bands_wann_*_plane.dat  -->
<domain listName="w222"   suffix="_grid"/>  <!-- bands_wann_*_grid.dat   -->
```

Whether a list traces a line, covers a plane or fills the zone is a property of the list,
not something FLEUR is told: every domain is interpolated the same way. With a single
domain you write no suffix and the files keep their base names.

`listName` names the input, `suffix` names the output, and they are separate on purpose --
two domains may share a list (the same path subdivided differently by `@npts`), and the
file names are read by other tools, so you choose them. At most one domain may go without
a suffix, and a repeated one is rejected.

The lists themselves come from `inp.xml`: `inpgen` writes them into `kpts.xml`, which is
`<xi:include>`d, and any other included file works too -- they are found by name with the
standard `t_kpts%read_kpts_by_name`, the same routine DFPT uses. There is no way to point
at a k-point file of your own, and that is deliberate: all input lives in `inp.xml`.

Declaring operators to interpolate but no domain is an error, not a silent no-op.

## What comes out

Real space, from `<operators_r>`:

| File | Content |
|---|---|
| `WF<n>_hr.dat` | `H(R)` in eV, Wannier90 `hr` format |
| `WF<n>_r.dat` | `A(R) = <0n\|r\|Rm>` in Å, the non-centred WYSV form |
| `WF<n>_bmn.dat` | `B(R) = <0n\|H r\|Rm>` in eV·Å — **not numerically validated yet** |
| `rspauli.1` | spin, 3 components |
| `anglmomrs.<n>` | orbital moment, 3 components |
| `rssocmat.1` | spin-orbit, the 2×2 spinor blocks |
| `wig_vectors` | the Wigner-Seitz `R` mesh |
| `WF<n>.amn`, `WF<n>.mmn` | the projections and overlaps, **written on one rank only** |

Interpolated, from `<interpolation>`: `bands_wann_<what>[_domain][_spinN].dat`, one row per
k-point, `kdist` first and then the bands. `<n>` and `_spinN` are the collinear spin
channel; with spinors there is only `WF1` and no suffix.

## Running it

- The wannierisation needs the **full BZ, Γ-centred**: generate the mesh with
  `inpgen -noKsym`. An irreducible mesh silently gives the wrong neighbours.
- Reproducibility: the minimisation has several basins and rounding decides which one.
  Fix `MKL_CBWR=AVX2` and `I_MPI_CBWR=1` to get the same answer twice on the same rank
  count.
- Across **different** rank counts, expect `Omega_I` to agree and `Omega_total` not to:
  the reduction order changes, and a different basin is a legitimate outcome, not a bug.
  Never compare outputs byte-for-byte between different numbers of ranks.

## Adding something

**A new operator** — a contraction over the states: see `../matrixelements/README.md`.
Nothing in this directory changes.

**A new interpolated output** — the operator matrix already exists and you want another
quantity out of it:

1. a row in `WANNIERLIB_INTERP` (`fleurinput/types_melem_optable.f90`), naming which
   catalogue entry it needs and its output basename;
2. a `CASE` branch in `postproc/melem_run.F90`;
3. an `<xsd:enumeration>` in `FleurInputSchema.xsd`.

If the quantity is one number per band per k, `postproc/melem_interpolate_op.F90` already
does the whole pipeline — rotate to the Wannier gauge, transform to `R`, transform back on
the output domain, diagonalise and project. It only needs the component count and a file
name. A quantity that is not of that shape gets its own driver, as the velocity and the
currents do.

**A new on-disk format** — one `CASE` in `postproc/melem_io.F90`, which is the only file
that knows any layout.

## Layout

```
wannierlib/           overlaps, projections, the Wannier90 adapter, the driver
wannierlib/postproc/  everything that needs the gauge: transform to R,
                      interpolation, O(R) export, output domains, file formats
../matrixelements/    the operator matrices themselves; knows nothing of Wannier
```

The dependency runs one way: `postproc` uses `matrixelements`, never the reverse.
