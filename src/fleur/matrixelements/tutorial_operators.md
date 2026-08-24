# Adding an operator to `matrixelements`

`README.md` next to this file is the checklist. This one holds what no single file can
hold: which of three routes you are on, the sequences that cross files, and the mistakes.

**The detail of each routine is in its own header, not here.** A contract restated away
from the code drifts, and this document has drifted before:

| For | Read the header of |
|---|---|
| the provider contract, and a kernel whose answer you know in advance | `types_matelements_template.F90` |
| the same shape with a real operator in it | `types_matelements_orbital.F90` |
| what the two tables are for and who owns them | `fleurinput/types_melem_optable.f90` |
| the Fourier core, and why it is operator-agnostic | `postproc/melem_ft.F90` |
| `H_W(k)`, and why every driver shares one | `postproc/melem_hamk.F90` |
| `A(R)`, `B(R)`, `F(R)`, `C(R)` and their references | `postproc/melem_coeff_{a,b,tensor}.F90` |
| every on-disk layout | `postproc/melem_io.F90` |
| the generic invariants, and why they warn instead of aborting | `melem_check.F90` |

---

## 1. First: which of three things are you adding?

This is the decision that saves the most work, and the one the tables in the code try to
make you take explicitly.

| What it is | Examples already in the tree | What you touch |
|---|---|---|
| **A real operator** — a contraction over the states, `O_mn(k) = <psi_m\|O\|psi_n>` | `spin`, `orbital`, `spin_orbit` | catalogue + provider + coarse slice + exposure (**route A**, §3) |
| **Only an exposure** — the coarse matrix already exists | `soc` (reuses `spin_orbit` under another name) | one table row + one branch (**route B**, §4) |
| **Not an operator at all** — it comes out of the neighbour overlaps or the eigenvalues | `hamiltonian`, `velocity`, `eigenstates`, `position`, `position_pw90`, `bmn`, `fmn`, `cmn` | a row with `operator=''` + a driver of its own (**route C**, §5) |

> **The rule in one line.** If it is not a contraction over the states, it does not belong
> in the catalogue — it would force the coarse pass to produce something nobody reads.

Route C is the majority: eight of the twelve exposed names go through it, and it is the
only route with no provider to copy.

---

## 2. The map: how a name travels

```
  inp.xml                types_wannierlib     t_melem_request      melem_coarse
  <operator name=".."/>  ---> reader ------->  validates ------->  builds O(k)
                                    ^               ^                   |
                                    |               |                   v
                          types_melem_optable ------+            postproc/melem_run
                          (both tables)                                 |
                                    |                        +----------+----------+
                                    +----------------------> |                     |
                                              melem_interpolate_op        melem_operators_r
                                              bands_wann_*.dat                 O(R)
```

Three consumers — the reader, the coarse pass and the dispatch — ask the same question of
the same two tables. The `operator` column is the pointer between them: it names the
catalogue entry needed to serve an exposed name, and is empty when none is.

Never copy an array size into prose: it is wrong on the next commit that adds a row, and
this document has been wrong that way twice. One line measures it:

```bash
grep -oE '(MELEM_OPERATORS|WANNIERLIB_INTERP|WANNIERLIB_OPR)\([0-9]+\)' \
     src/libraries/fleurinput/types_melem_optable.f90
```

`postproc` depends on `matrixelements` and the dependency never runs the other way;
`testing/tests/structure/test_layering.py` fails if it ever does.

---

## 3. Route A: a genuinely new operator

Copy `types_matelements_template.F90`, whose header carries the contract. What follows is
only what that header cannot see, because it spans more than one file.

### The two flags decide the shape

`init_mat` hands out one block or a 2x2 in spin space:
`nsp = MERGE(2, 1, spinorwavefcts .OR. spinoroperator)`. The three existing operators cover
three different combinations, which is why the comparison lives here and not in any of them:

| Operator | `spinoroperator` | `spinorwavefcts` | Result lands in |
|---|---|---|---|
| `spin` | `.TRUE.` | `.TRUE.` | `mat(2,2)` — the four blocks of sigma |
| `spin_orbit` | `.TRUE.` | `.FALSE.` | `mat(2,2)` |
| `orbital` | `.FALSE.` | `.FALSE.` | `comp(:,:,1:3)` — L is spin-diagonal |

`mat` carries spin structure; `comp(band, band, alpha)` carries Cartesian indices, which
`mat` cannot hold. **`comp` has no distributed counterpart**: an operator with components
runs k-parallel, never with eigenvector parallelism.

### The array you will spend all your time indexing

`types_abc.F90` names the axes and stops there, so the meanings are here:

```fortran
abc(s, ntyp)%cof(nu, lm, iOrd, iAtom)
!  nu    = band, within the selected window (NOT the global band index)
!  lm    = l(l+1) + m, COMPLEX spherical harmonics, m from -l to +l
!  iOrd  = 1 = u, 2 = udot, 3+ = local orbitals
!          and its bound is abc%n_r(l), which CHANGES with l
!  iAtom = the equivalent atom of that type
!  s     = spin component
```

Three things no declaration shows:

- **`SIZE(zmat)` is not `jspins`.** It is 1 when the record holds the whole spinor and 2
  when they are independent channels. A consumer that addresses a spin block by row offset
  needs the first case; with the second it reads one row past the end of the record, and
  what it finds there is plausible rather than obviously wrong.
- **`radfun%integral` is allocated `(.,.,.,jspins,jspins)`**, so a component has to be
  clamped to the sets that exist with `radial_slot(radfun, s)` — once, outside the loop.
- **The coefficients are in each atom's LOCAL frame.** A quantity that is not invariant
  under that rotation has to be rotated before the sites are summed (`t_abc%rotate`). This
  is what an antiferromagnet gets wrong if you forget.

And one that leaves a mark on disk: **assign, do not accumulate** when each element gets
exactly one contribution — accumulating onto a cleared matrix turns `-0.0` into `+0.0`.

### The coarse slice: three places in `melem_coarse.F90`

Declare `myo0(nb, nb, 3, nk_loc)` as `ALLOCATABLE` in `t_melem_coarse`; allocate it in
`%init` **behind `request%needs_op`**; fill it in `%calc` inside the loop over this rank's
k-points, calling the provider's `%init` **outside** that loop, because what an instance
binds to — a site, and a channel when there is more than one — does not depend on `k`.

> **Never allocate a `(1,1,1,1)` stub** for something nobody asked for: it is a *valid*
> array of the wrong shape, so indexing it wrongly returns whatever is there instead of
> failing (§7). Unallocated, the receiver is forced to declare it `ALLOCATABLE` and decide.

> Ask `request%needs_op('name')`, **not** the `l_spin` / `l_orbmom` / `l_socop` flags. Those
> summarise both exposure lists and cannot tell an interpolation request from a real-space
> one — exactly the distinction that matters, since the two are not producible in the same
> spin configurations.

Then expose it: §6.

---

## 4. Route B: only an exposure

The coarse matrix already exists and you want another quantity derived from it. No provider,
no new slice: one row pointing at the existing catalogue entry.

```fortran
t_melem_exposed('soc', 'spin_orbit', 1, 'bands_wann_soc', '', .FALSE.)
```

That row says three things at once: the name is accepted in `<interpolation>`, serving it
requires the `spin_orbit` matrix even though the user never writes that word, and the result
has one component.

Every column of `t_melem_exposed` is documented on the field itself, in
`types_melem_optable.f90`. The one worth repeating: **`out1`/`out2` are not decorative** —
the per-domain renaming reads them from the table, so leaving them empty makes every domain
write the same file and the second overwrite the first.

Then expose it: §6.

---

## 5. Route C: not an operator at all

`operator=''` — nothing in the catalogue serves it — and `ncomp=0`, which is how the table
says *this name has a driver of its own*. No provider, no coarse slice. The whole question
is where your driver gets its input.

| Source | What it gives you | Built by | Serves |
|---|---|---|---|
| the eigenvalues and the gauge | `H_W(k)` | `melem_hamk.F90` | `hamiltonian`, `eigenstates`, `velocity` |
| the neighbour overlaps and the b-shell weights | `A^(W)(R)` | `melem_coeff_a.F90` | `position`, `position_pw90`, `velocity`, `bmn` |
| the wavefunctions at **two** neighbours at once | `F(k)` / `C(k)` | `wannierlib_uiu` / `_uhu` | `fmn`, `cmn` |

The Fourier transforms are not yours to write: `melem_ft.F90` interpolates any `mat_k` on
the coarse mesh, and `melem_interp_util.F90` has the abscissa and the LAPACK query. Copy the
signature of any `melem_interpolate_*.F90`; the four are deliberately the same.

**The one rule that decides everything else** follows from the table above:

> A per-rank input (`mmn_loc`, `f0_loc`, a coarse slice) means your driver is **collective**
> and every rank must reach it. A complete gauge means **rank 0 only**, and the
> `IF (irank /= 0) RETURN` goes first, before anything touches `kfrac` — which is
> unallocated elsewhere, and reading it is undefined behaviour a two-rank suite can pass
> without showing. This layer has produced that bug three times.

**And one exception.** `fmn` and `cmn` need the gauge at two neighbouring k-points at once,
which the coarse pass cannot do (no gauge yet) and `postproc` cannot do (no wavefunctions).
Those two are built in `wannierlib_main.F90`, between the wannierisation and `melem_run`,
behind `IF (request%has_op_r(...))`. If your quantity needs the wavefunctions after the
gauge is known, you are in this case; if it does not, stay out of that file.

Then expose it: §6.

---

## 6. Exposing it

The same three steps for all three routes. What differs is only whether the branch calls
the generic driver (A and B) or yours (C).

### For interpolation — `<interpolation><operator name=".."/>`

A row in `WANNIERLIB_INTERP` and a branch in the `SELECT CASE` of `postproc/melem_run.F90`:

```fortran
CASE ('my_operator')
   CALL melem_interpolate_operator(manifold, cell, kpts, eig, u_matrix, u_opt, &
                                   coarse%myo0, gk_loc, 3, kfrac, outname(iop, 1), &
                                   irank, mpicm)
```

`kfrac` is the domain's k-set and `outname(iop, 1)` its basename with the suffix already on
it, both resolved by the loop around this `SELECT CASE` — neither is something a new
operator builds. For routes A and B that is all: the generic driver does
`O_W(k) = V^dagger O(k) V`, the transform to `R` and back, the diagonalisation of `H(k')`
and the projection, **for any operator**. A route-C name calls its own driver here instead.

### For the real-space export — `<operators_r><operator name=".."/>`

A row in `WANNIERLIB_OPR` and a branch in `postproc/melem_operators_r.F90`:

```fortran
CASE ('my_operator')
   CALL melem_op_rs_distributed(this, cell, kpts, vloc, coarse%myo0, gk_loc, 3, &
                                mpicm, irank, .FALSE., 'myfile.1')
```

`WANNIERLIB_OPR` has no output columns: these files carry a channel or a spin-block index,
so the writer builds the name — and `melem_write_realspace` does not append `.dat`. A new
on-disk layout is a `CASE` in `melem_io.F90`, the only file that knows any. Do not open a
file anywhere else.

### The schema

```xml
<xsd:enumeration value="my_operator"/>
```

In `FleurInputSchema.xsd`, inside `WannierlibOperatorType` or `WannierlibOperatorRType`.
Without it the input does not validate. **Add the row and the enumeration in the same
commit:** the two lists have to agree and nothing enforces it. If they drift, an input
passes the schema and then stops in `melem_request_init` — a clean failure, but the schema
is the first thing a new user reads and should not promise what does not exist.

---

## 7. Things that have actually cost us a run

| Symptom | Cause and fix |
|---|---|
| Error `#6401` at compile time | You named a dummy argument `mpi_comm` and it clashes with `use mpi`. Use `mpicm`. |
| Error `#6404` on a constant | A missing `USE m_constants, ONLY: ...` (e.g. `tpi_const`, `ImagUnit`). |
| Error `#6632` on a routine you only moved | The destination module does not `USE m_judft`, so `juDFT_error(..., calledby=...)` has keyword arguments with no explicit interface. Check what the **destination** imports. |
| **Works in serial; garbage or a hang with `np>1`** | An unallocated `ALLOCATABLE` reached a **collective** routine with a non-allocatable dummy. Undefined behaviour: a serial run cannot show it and a two-rank suite can pass with it present. |
| A layer `l = lmax` vanishes silently | `lm = lmax**2` instead of `lmax*(lmax + 2)`. |
| Unexpected precision | In FLEUR a bare `real` **is double precision**. |
| The operator comes out tiny instead of failing | A stub slice reached the export without ever being computed. |
| The link fails, unrelated to what you touched | The `fleurinput` CMakeLists lists files without the path prefix; this directory's lists them with it. |

Two about the build: **never compile on the login node**, and an `ifx` internal compiler
error has two causes that one test separates — if it survives at `-O0` it is a stale build
state, so reconfigure from scratch before you look at what you wrote.

---

## 8. Validating it

`melem_check` gives any provider three invariants for free — finite, non-zero, hermitian —
and its header says what each one catches and why it warns instead of aborting. Ask for
hermiticity only of the slices that have it: a cross-spin block does not.

> **Byte-identity at the same rank count; `Omega_I` across different rank counts.**

Between different rank counts byte-identity does not apply: the reduction tree changes and
the sums happen in another order. `MKL_CBWR` and `I_MPI_CBWR` fix the code path, not the
tree. What survives is `Omega_I`, which is gauge- and basin-invariant — exactly to the last
digit without disentanglement, approximately with it, while `Omega_total` always moves. A
change of basin is legitimate and is not a bug.

A free trick falls out: running at two rank counts **classifies** the outputs, because the
ones that survive identical are by construction the gauge-invariant ones.

A new operator is inert unless someone asks for it — no slice, no provider call, no file.
The only real failure mode is forgetting the gate, which is the stub in §3.

Which suite to run and where its inputs are: `README.md`.
