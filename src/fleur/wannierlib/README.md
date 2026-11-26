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
| 4. post-processing | `postproc/wgauge_run` | `O(R)` files and interpolated bands |

Stage 2 needs no gauge and runs before the wannierisation; stages 3 and 4 do. With
`jspins=2` and no SOC the two spin channels wannierise **separately**, so stages 1–4 run
once per channel.

## Turning it on

```xml
<wannierlib wannierize="T">
  <bands numBands="36" minBand="7" maxBand="42"/>
  <disentanglement disWinMin="0.097" disWinMax="2.670"
                   disFrozMin="0.097" disFrozMax="0.661"
                   disNumIter="3000" disConvTol="0.00001" disMixRatio="0.5"/>
  <wannierization numIter="3000" convTol="0.00001"/>
  <operators>
    <hamiltonian RS="T" INT="T"/>
    <spin        RS="T"/>
    <velocity           INT="T"/>
  </operators>
  <interpolation>
    <domain listName="path-2"/>
  </interpolation>
  <export>
    <wannier90 on="T"/>
  </export>
</wannierlib>
```

`<disentanglement>` is optional and is left out when `numBands == num_wann`. The initial
guesses are `<wannierproj l=".." m=".."/>` children of a `<species>`; their total
count over all atoms is `num_wann`.

`<operators>` holds one element per operator and, on it, one switch per basis: `KS="T"`
hands over what the run holds in `k`, `RS="T"` writes `O(R)` to disk and does no
interpolation, and `INT="T"` projects onto the interpolated bands. The three are
independent -- any may be off, and asking for a quantity in one basis never obliges
another. They accept different *sets* of names, because not everything with an `O(R)` has
an interpolation driver and the other way round, and a switch set on a basis that does not
offer that name is refused when the input is read. But **a name never means two different
things** -- the spin-orbit operator is `spinOrbit` everywhere, and `bmn` is the same `B`
under `KS` and under `RS`, one Fourier sum apart. The accepted names of each live in
`fleurinput/types_melem_optable.f90`, and a name the tables do not have stops the run with
the accepted list in the message.

`KS="T"` carries two kinds of name and writes them to two files. `spin`, `orbital`,
`spinOrbit` and `exchangeTorque` are a matrix at a single `k`: they go out in the Bloch
basis, before any gauge exists, and that half costs no wannierisation. `bmn`, `cmn` and `fmn`
are contractions over *pairs* of neighbours, and summing over `b` needs a band index that
means the same state at `k+b1` and at `k+b2`, which the Bloch one does not -- the gauge that
enters is the one at the neighbours, not at `k`, so there is nothing to factor out and no
gauge-free `B(k)`, `C(k)` or `F(k)` exists. Those three therefore carry the Wannier gauge,
and the two kinds go to different files so that no consumer has to guess which it has.

`<export>` is a separate thing: it hands out matrices the run already holds, one element
per consumer (`<wannier90>`, `<irrep>`, `<wannierberri>`, `<plotWF>`, `<orbitrans>`,
`<juwip>`), each with an `on` switch, all off by
default and none of them exclusive. It computes nothing.

## Output domains

Inside `<interpolation>`, declare one `<domain>` per set of k-points you want the operators
evaluated on. There is no limit and no kinds: a domain is a **named `kPointList`** plus,
when there is more than one, a suffix that keeps their files apart.

```xml
<domain listName="path-2"/>                 <!-- group "default"         -->
<domain listName="kz0"    suffix="_plane"/> <!-- group "plane"           -->
<domain listName="w222"   suffix="_grid"/>  <!-- group "grid"            -->
```

Whether a list traces a line, covers a plane or fills the zone is a property of the list,
not something FLEUR is told: every domain is interpolated the same way. With a single
domain you write no suffix and the files keep their base names.

`listName` names the input, `suffix` names the output, and they are separate on purpose --
two domains may share a list (the same path subdivided differently by `@refinementPoints`), and the
file names are read by other tools, so you choose them. At most one domain may go without
a suffix, and a repeated one is rejected.

The lists themselves come from `inp.xml`: `inpgen` writes them into `kpts.xml`, which is
`<xi:include>`d, and any other included file works too -- they are found by name with the
standard `t_kpts%read_kpts_by_name`, the same routine DFPT uses. There is no way to point
at a k-point file of your own, and that is deliberate: all input lives in `inp.xml`.

Declaring operators to interpolate but no domain is an error, not a silent no-op.

## What comes out

Real space, from `RS="T"`:

| File | Content |
|---|---|
| `WF<n>_hr.dat` | `H(R)` in eV, Wannier90 `hr` format |
| `WF<n>_r.dat` | `A(R) = <0n\|r\|Rm>` in Å, the non-centred WYSV form |
| `WF<n>_bmn.dat` | `B(R) = <0n\|H r\|Rm>` in eV·Å; reproduced to 6e-08 by rebuilding it outside FLEUR from the Bloch-basis overlaps |
| `rspauli.1` | spin, 3 components |
| `anglmomrs.<n>` | orbital moment, 3 components |
| `rssocmat.1` | spin-orbit, the 2×2 spinor blocks |
| `wig_vectors` | the Wigner-Seitz `R` mesh |
| `WF<n>.amn`, `WF<n>.mmn`, `WF<n>.eig` | from `<export><wannier90 on="T"/></export>`; also written unasked on a serial run, as a path for inspection |

In `k`, from `KS="T"`:

| File | Content |
|---|---|
| `WF<n>_operators_ks.hdf` | `spin`, `orbital`, `spinOrbit`, `exchangeTorque` in the Bloch basis, no gauge, with the b-shells beside them. Carries the k-list every dataset is indexed by |
| `WF<n>_modern_k.hdf` | `B(k)`, `C(k)`, `F(k)` in the Wannier gauge, one Fourier sum from the three `.dat` files above and in the same units, so the step between them is a plain phase sum. Carries the k-list it was built on, so that sum needs no second file |

Each of the three HDF files above also carries, on its root group, what the Wannier90
checkpoint beside it cannot say about itself: the route, the band range, the energy
windows that were actually used, and `chkCRC32`, a CRC32 of the checkpoint's own bytes.
Wannier90's ten load-time checks answer "is this checkpoint shaped like this case"; these
answer "is it from this run". Checking the stamp costs a consumer one `zlib.crc32` call.
Strings are stored the way every FLEUR attribute is, as an array of one-character strings.

Interpolated, from `<interpolation>`: `WF<channel>_interp.hdf`, one group per domain.
Each group carries `kdist`, `energy` and a `kpts` subgroup, plus one dataset per operator
asked for. `<channel>` is the collinear spin channel; with spinors there is only `WF1`.

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
2. a `CASE` branch in `postproc/wgauge_run.F90`;
3. an `<xsd:enumeration>` in `FleurInputSchema.xsd`.

If the quantity is one number per band per k, `postproc/wgauge_interpolate_op.F90` already
does the whole pipeline — rotate to the Wannier gauge, transform to `R`, transform back on
the output domain, diagonalise and project. It only needs the component count and a file
name. A quantity that is not of that shape gets its own driver, as the velocity and the
currents do.

**A new on-disk format** — one `CASE` in `postproc/wgauge_io.F90`, which is the only file
that knows any layout.

## Layout

```
wannierlib/           the path to U: overlaps, projections, the Wannier90 adapter,
                      the driver that runs them in order
wannierlib/export/    what is done with U once it exists: the basis, the gauge, the
                      Bloch coefficients, the plots, the C and F tensors
wannierlib/postproc/  everything that needs the gauge: transform to R,
                      interpolation, O(R) export, output domains, file formats
../matrixelements/    the operator matrices themselves; knows nothing of Wannier
```

The dependency runs one way: `postproc` uses `matrixelements`, never the reverse. And of
everything in `wannierlib/`, only the driver names `export/` -- the rest is the path to U and
reads as that. Both are asserted in `testing/tests/structure/test_layering.py`.
