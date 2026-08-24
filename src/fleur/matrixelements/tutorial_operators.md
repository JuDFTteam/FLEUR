# Adding an operator to `matrixelements`

The long form. `README.md` next to this file is the checklist; this one carries the
reasoning, the contract of each routine you have to write, and a worked example.

---

## 1. First: which of three things are you adding?

This is the decision that saves the most work, and the one the tables in the code try to
make you take explicitly.

| What it is | Examples already in the tree | What you touch |
|---|---|---|
| **A real operator** — a contraction over the states, `O_mn(k) = <psi_m\|O\|psi_n>` | `spin`, `orbital`, `spin_orbit` | catalogue + provider + coarse slice + exposure (**route A**, §3) |
| **Only an exposure** — the coarse matrix already exists; you want another file or a derived quantity | `soc` (reuses `spin_orbit` under another name) | one table row + one branch (**route B**, §4) |
| **Not an operator at all** — it comes out of the neighbour overlaps or the eigenvalues | `hamiltonian`, `velocity`, `eigenstates`, `position`, `position_pw90`, `bmn`, `fmn`, `cmn` | a row with `operator=''` + a driver of its own (**route C**, §5) |

> **The rule in one line.** If it is not a contraction over the states, it does not belong
> in the catalogue. `hamiltonian` is not an operator — it is `V^dagger diag(E) V`.
> `position` is not either — it comes out of the `M_mn`. Putting them in the catalogue
> would force the coarse pass to produce something nobody reads.

Route C is the majority: eight of the twelve exposed names go through it. It is also the
route with no provider to copy, which is why §5 is the longest section here.

---

## 2. The map: two tables, and how a name travels

All the routing lives in one file, `fleurinput/types_melem_optable.f90`, and it holds
**two** tables with different owners.

```fortran
!> WHAT CAN BE COMPUTED. They know nothing about Wannier functions.
TYPE(t_melem_operator), PARAMETER, PUBLIC :: MELEM_OPERATORS(...) = [ &
   t_melem_operator('spin',       3), &
   t_melem_operator('orbital',    3), &
   t_melem_operator('spin_orbit', 4)]

!> HOW THE WANNIERISATION EXPOSES THINGS
!                   name          operator     ncomp  out1                out2  total
TYPE(t_melem_exposed), PARAMETER, PUBLIC :: WANNIERLIB_INTERP(...) = [ &
   t_melem_exposed('spin',        'spin',       3, 'bands_wann_spin',   '', .TRUE.), &
   ... ]
TYPE(t_melem_exposed), PARAMETER, PUBLIC :: WANNIERLIB_OPR(...) = [ ... ]
```

The array sizes are elided on purpose: a count copied into prose is wrong on the next
commit that adds a row, and this document has been wrong that way before. The file is the
authority, and one line reads it:

```bash
grep -oE '(MELEM_OPERATORS|WANNIERLIB_INTERP|WANNIERLIB_OPR)\([0-9]+\)' \
     src/libraries/fleurinput/types_melem_optable.f90
```

The `operator` column is the pointer from one table to the other: it says which catalogue
entry has to be built to serve that name, and it is empty when none does. That is what
turns *"which coarse matrix does this name need?"* from a fact to be remembered in four
places into one lookup, `melem_exposed_find`.

A name travels like this:

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

The three consumers of the tables — the reader, the coarse pass and the dispatch — ask the
same question and get the same answer. That is the whole point of centralising them.

Which directory holds what is in `README.md` under *Where things live*. The one invariant
worth repeating: `postproc` depends on `matrixelements`, and the dependency never runs the
other way. `testing/tests/structure/test_layering.py` fails if it ever does.

---

## 3. Route A: a genuinely new operator

### 3.1 The provider

Create `types_matelements_<name>.F90` extending `t_matelements`. The contract is two
procedures: `init`, which fixes what does not depend on `k`, and `calc_matrix_elements`,
which is what does.

**Start from `types_matelements_template.F90`.** It is a working skeleton whose kernel
computes the muffin-tin overlap — the identity restricted to the spheres — so its answer
is known in advance: the diagonal is the muffin-tin charge of the state, hence real,
positive and below one. You replace the angular factor, not the scaffolding.

### 3.2 The three variables that decide everything

Set in `init`:

```fortran
this%spinoroperator = .FALSE.   !> does the operator act in spin space?
this%spinorwavefcts = .FALSE.   !> are the wavefunctions spinors?
this%n_alpha        = 3         !> Cartesian components it produces
```

The first two decide whether `init_mat` hands out **one block or a 2x2** in spin space:
`nsp = MERGE(2, 1, spinorwavefcts .OR. spinoroperator)`. The three existing operators cover
three different combinations:

| Operator | `spinoroperator` | `spinorwavefcts` | Result lands in |
|---|---|---|---|
| `spin` | `.TRUE.` | `.TRUE.` | `mat(2,2)` — the four blocks of sigma |
| `spin_orbit` | `.TRUE.` | `.FALSE.` | `mat(2,2)` |
| `orbital` | `.FALSE.` | `.FALSE.` | `comp(:,:,1:3)` — L is spin-diagonal |

**The two stores are different and not interchangeable.** `mat` carries spin structure;
`comp(band, band, alpha)` carries Cartesian indices, which `mat` cannot hold. An operator
with two Cartesian indices uses `n_alpha` and `n_beta`, and the pairs are stored row-major.

> **`comp` has no distributed counterpart.** `init_mat` refuses Cartesian components
> together with an `mpi_subcomm`, and the factory checks it earlier so the error arrives
> before a sub-communicator has been chosen. An operator with components runs k-parallel,
> never with eigenvector parallelism.

### 3.3 The contract of `calc_matrix_elements`

Four things arrive, none of them optional:

```fortran
TYPE(t_mat),    INTENT(IN) :: zmat(:)    !> ONE matrix if it is a whole spinor,
                                         !> TWO if the records are independent channels
TYPE(t_abc),    INTENT(IN) :: abc(:,:)   !> (2, ntype) matching coefficients
TYPE(t_radfun), INTENT(IN) :: radfun(:)  !> (ntype) radial functions + MT integrals
TYPE(t_usdus),  INTENT(IN) :: usdus      !> values and derivatives at the MT boundary
```

The array you will spend all your time indexing:

```fortran
abc(s, ntyp)%cof(nu, lm, iOrd, iAtom)
!  nu    = band, within the selected window (NOT the global band index)
!  lm    = l(l+1) + m, COMPLEX spherical harmonics, m from -l to +l
!  iOrd  = 1 = u, 2 = udot, 3+ = local orbitals
!          and its bound is abc%n_r(l), which CHANGES with l
!  iAtom = the equivalent atom of that type
!  s     = spin component
```

Three things that are not visible from the declaration:

- **`SIZE(zmat)` is not `jspins`.** It is 1 when the record holds the whole spinor and 2
  when they are independent channels. A consumer that addresses a spin block by row offset
  needs the first case; with the second it reads one row past the end of the record, and
  what it finds there is plausible rather than obviously wrong.
- **`radfun%integral` is allocated `(.,.,.,jspins,jspins)`**, so a component has to be
  clamped to the sets that exist with `radial_slot(radfun, s)`. Resolve it **once, outside
  the loop**: inside it is the same number fetched a million times.
- **The coefficients are in each atom's LOCAL frame.** A quantity that is not invariant
  under that rotation has to be rotated before the sites are summed. `t_abc` has `rotate`
  and `rot_to_unrotated` for it. This is what an antiferromagnet gets wrong if you forget.

And one subtlety that leaves a mark in an exported file: **assign, do not accumulate**,
when each element receives exactly one contribution. Accumulating onto a cleared matrix
turns a `-0.0` result into `+0.0`, and the sign of zero is visible on disk.

```fortran
this%comp(i, j, 1) =  0.5 * (cp + cm)              ! L_x = (L+ + L-)/2
this%comp(i, j, 2) = -0.5 * ImagUnit * (cp - cm)   ! L_y = (L+ - L-)/(2i)
this%comp(i, j, 3) =  cz                           ! L_z
```

### 3.4 Registering it

```fortran
! types_melem_optable.f90
t_melem_operator('my_operator', 3)
```

```cmake
# src/fleur/matrixelements/CMakeLists.txt --- WITH the path prefix
src/fleur/matrixelements/types_matelements_myop.F90
```

> The `CMakeLists.txt` under `fleurinput` lists files **without** a path prefix, unlike this
> one. Both conventions live in the tree, and confusing them produces a link failure much
> later with no apparent connection to what you touched.

### 3.5 The coarse slice

Three places in `melem_coarse.F90`, always the same three:

```fortran
! (a) the slice, inside t_melem_coarse
COMPLEX, ALLOCATABLE :: myo0(:, :, :, :)   !< (nb, nb, 3, nk_loc)

! (b) allocate it in %init, ALWAYS behind the request
IF (request%needs_op('my_operator')) &
   ALLOCATE(this%myo0(nb, nb, 3, nkc_loc), source=cmplx(0.0, 0.0))

! (c) fill it in %calc, inside the loop over this rank's k-points
CALL myop%init(atoms, itype, iatom)
CALL matrix_element_factory(myop, eig_id, ikpt, input, atoms, sym, cell, noco, &
                            nococonv, enpara, lapw, vtot, fmpi, ev_list=ev_list, &
                            l_both_spinors=l_spinor_records, kpts=kpts)
this%myo0(:, :, 1:3, il) = myop%comp(:, :, 1:3)
```

> **Never allocate a `(1,1,1,1)` stub** for something nobody asked for. A stub is a *valid*
> array of the wrong shape: indexing it wrongly does not fail, it returns whatever is
> there. Unallocated there is nothing to index wrongly, and whoever receives it is forced
> to declare it `ALLOCATABLE` and decide what to do when it is absent.

> Ask `request%needs_op('name')`, **not** the `l_spin` / `l_orbmom` / `l_socop` flags. Those
> three summarise both exposure lists and cannot tell an interpolation request from a
> real-space one — which is exactly the distinction that matters, since the two are not
> producible in the same spin configurations.

The provider's `init` goes **outside** the k loop: what an instance binds to — a site, and
a channel when there is more than one — is the same at every k. The k dependence arrives
with the coefficients, and `init_mat` clears the result while reusing the allocation.

Then expose it: §6.

---

## 4. Route B: only an exposure

When the coarse matrix already exists and what you want is another quantity derived from
it, there is no provider and no new slice. One row pointing at the existing catalogue entry
is enough:

```fortran
t_melem_exposed('soc', 'spin_orbit', 1, 'bands_wann_soc', '', .FALSE.)
```

That row says three things at once: the name is accepted in `<interpolation>`, serving it
requires the `spin_orbit` matrix — even though the user never writes that word — and the
result has one component. The reader, the coarse pass and the dispatch all read it from
there.

| Column | Meaning |
|---|---|
| `name` | what the user writes in `inp.xml` |
| `operator` | the catalogue entry to build; `''` if none |
| `ncomp` | components the generic driver writes; `0` = it has a driver of its own |
| `out1`, `out2` | output basenames, without `.dat` and without the domain suffix |
| `honours_total` | whether it respects `total="T"`, i.e. has a site-summed projection |

`out1` and `out2` are not decorative: the per-domain renaming reads them from the table. An
operator that leaves them empty writes the same file every time and the second domain
overwrites the first.

Then expose it: §6.

---

## 5. Route C: not an operator at all

The row carries `operator=''` — nothing in the catalogue serves it — and `ncomp=0`, which
is how the table says *this name has a driver of its own*. Both columns are read: the
coarse pass builds no slice for you, and the generic driver will not be called.

So route C has no provider and no coarse slice. What it has instead is **a driver you
write**, and the whole question is where your driver gets its input.

### 5.1 The three sources, and which of them is yours

Every route-C quantity in the tree comes from one of three places. Finding yours here is
the equivalent of choosing `spinoroperator` in route A: it decides the shape of everything
that follows.

| Source | What it gives you | Built by | Serves |
|---|---|---|---|
| the eigenvalues and the gauge | `H_W(k)`, the Wannier-gauge Hamiltonian | `melem_hamk.F90` | `hamiltonian`, `eigenstates`, `velocity` |
| the neighbour overlaps `M_mn` and the b-shell weights | `A^(W)(R)`, the Berry connection | `melem_coeff_a.F90` | `position`, `position_pw90`, `velocity`, `bmn` |
| the wavefunctions at **two** neighbours at once | `F(k)` / `C(k)`, the geometric tensors | `wannierlib_uiu` / `wannierlib_uhu` | `fmn`, `cmn` |

The first two are free: they are built from what `melem_run` is already handed. The third
is not, and §5.5 is about the price.

### 5.2 What you never write yourself

The Fourier machinery is operator-agnostic and already there. `melem_ft.F90` is the core:

```fortran
melem_ft_interpolate(cell, mat_k, kpts, kfrac, mat_interp)        ! k -> R -> k', one call
melem_ft_to_real_reduce(cell, kpts, mat_loc, gk_loc, commw, ...)  ! COLLECTIVE, k -> R
melem_ft_rtok(mat_r, irvec, ndegen, nrpts, kfrac, mat_interp)     ! R -> k'
```

It is agnostic in the strict sense: feed it `H_W(k)` and you interpolate bands, feed it a
spin matrix and you interpolate spin. **A new route-C driver builds its own `mat_k` and
reuses this unchanged.** The Wigner-Seitz `R` set is cached inside, keyed on the mesh, so
asking twice costs nothing.

The two halves are also available separately (`melem_ft_to_real`), and the way back has a
velocity variant (`melem_ft_rtok_velocity`, which multiplies by `i R_cart` on the fly).

Two more pieces of scaffolding, in `melem_interp_util.F90`:

```fortran
melem_kpath(cell, kfrac, kdist)                      ! the abscissa every .dat starts with
melem_zheev_workspace('N'|'V', n, work, rwork, lwork)! the LAPACK query, done once
```

`melem_kpath` returns a cumulative distance. **Only a path gives that a physical meaning**;
a plane or a grid still gets a monotonically growing number, which is why those domains
also write the fractional coordinates.

### 5.3 An interpolation driver: rank 0 only

The four that exist share one signature, and it is worth copying literally:

```fortran
SUBROUTINE melem_interpolate_<what>(this, cell, kpts, eig, u_matrix, u_opt, kfrac, &
                                    out1[, out2], irank)
   REAL, ALLOCATABLE, INTENT(IN) :: kfrac(:, :)   !> (3, np) the domain's k-set
   ...
   IF (irank /= 0) RETURN          ! only the master holds the full U(k)
   np = SIZE(kfrac, 2)
   CALL melem_kpath(cell, kfrac, kdist)
   CALL melem_build_hamk(this, eig, u_matrix, u_opt, ham_k)
   CALL melem_ft_interpolate(cell, ham_k, kpts, kfrac, H_interp)
   ...                              ! your own step, then OPEN(TRIM(out1)//'.dat')
```

Three things about it that are not decoration:

- **The early return comes first, before anything touches `kfrac`.** `kfrac` is allocated
  on rank 0 only, because it is the k-set of a named `kPointList` that rank 0 read. Off
  rank 0 it is an unallocated `ALLOCATABLE` and reading it is undefined behaviour — which a
  two-rank test suite can pass without showing. This layer has produced that bug three
  times; the early return is the whole defence.
- **`out1`/`out2` arrive built**, suffix included: the driver appends `.dat` and nothing
  else. Never build an output name from a flag.
- **`melem_build_hamk` is shared on purpose.** The band driver needs `H_W(k)` for the
  bands, the operator drivers for the eigenvectors they project on, the velocity to
  differentiate it. A correction to `eigval2` has one place to reach.

### 5.4 A real-space driver: collective, and rank 0 writes

The `O(R)` side is the mirror image. The reduce is collective, so **every rank must reach
it**, and only the last step is rank 0:

```fortran
CALL melem_ft_to_real_reduce(cell, kpts, o_loc, gk_loc, mpicm, o_r, irvec, ndegen, nrpts)
IF (irank == 0) CALL melem_write_realspace(o_r, irvec, ndegen, nrpts, nw, ncomp, &
                                           'generic', 'myfile.dat', irank)
```

Whether your driver is collective or rank-0-only follows from §5.1 and nothing else: a
per-rank input (`mmn_loc`, `f0_loc`, a coarse slice) means collective, a complete gauge
means rank 0. Getting it backwards is the *"works in serial, garbage or a hang with
`np>1`"* line of §7, and it is the single most expensive mistake in this layer.

One consequence worth knowing before it surprises you: **`H(R)` is built on every rank**
even though only rank 0 writes it, because `melem_write_bmn` is collective and needs its
`R` set. Duplicated work on every rank was the cheap way out of handing an unallocated
array to a collective callee.

The layouts `melem_write_realspace` knows are `hr`, `r`, `bmn`, `soc`, `generic`, `cart2`
and `cart2e`. If yours is none of them, add a `CASE` there — **it is the only file that
knows any layout.** Do not open a file anywhere else. And note that the writer does not
append `.dat`: the caller passes the whole name, because some of these files are legacy
names with a channel index in them (`rspauli.1`, `anglmomrs.<n>`).

### 5.5 The one case that touches `wannierlib_main`

`fmn` and `cmn` need the gauge at two neighbouring k-points **at the same time**. That
cannot happen in the coarse pass, which runs before any gauge exists, and it cannot happen
in `postproc`, which can no longer reach the wavefunctions. So those two are built in
`wannierlib_main.F90`, between the wannierisation and `melem_run`, and only when asked for:

```fortran
IF (request%has_op_r('fmn')) &
   CALL wannierlib_uiu(manifold, bmesh, kpts, ..., u_matrix, u_opt, f0_loc)
```

`f0_loc` and `c0_loc` then travel into `melem_run` as `ALLOCATABLE` arguments, unallocated
unless they were requested — which is why `melem_run` declares them `ALLOCATABLE` rather
than assumed-shape.

**This is the exception to the rule that adding an operator does not touch
`wannierlib_main`.** If your quantity needs the wavefunctions after the gauge is known, you
are in this case and there is no way around it. If it does not, stay out of that file.

### 5.6 Worked example: `velocity`, which uses two sources at once

`v_alpha = dE_n/dk_alpha` is route C with both of the free sources in play:

1. `H_W(k)` from `melem_build_hamk` — the same one the band driver builds.
2. `A^(W)(R)` from `melem_build_berry_aw_r`, reduced from this rank's overlaps. Collective,
   so it runs on all ranks; the centre check that calibrates the sign runs on rank 0.
3. `v_W(k') = FT[i R_cart H_W]` through `melem_ft_rtok_velocity`, then `H(k')` diagonalised,
   then `<v>_n = [C^dagger v_W C]_nn`.

Two details of the dispatch are worth copying. `A^(W)(R)` is built **once and reused across
output domains** — the `IF (.NOT. ALLOCATED(aw_r))` in `melem_run` — because it does not
depend on where you evaluate. And the diagonal `<n|v|n>` needs no Berry-connection
correction, which is why this driver is exact where an off-diagonal one would not be.

Then expose it: §6.

---

## 6. Exposing it

The same three steps for all three routes. What differs is only whether the branch calls
the generic driver (A and B) or yours (C).

### For interpolation — `<interpolation><operator name=".."/>`

A row in `WANNIERLIB_INTERP` and a branch in the `SELECT CASE` of
`../wannierlib/postproc/melem_run.F90`:

```fortran
CASE ('my_operator')
   CALL melem_interpolate_operator(manifold, cell, kpts, eig, u_matrix, u_opt, &
                                   coarse%myo0, gk_loc, 3, kfrac, outname(iop, 1), &
                                   irank, mpicm)
```

`kfrac` is the current domain's k-set and `outname(iop, 1)` its output basename, both
resolved by the loop around this `SELECT CASE`: the name comes from the `out1` column of
your row with the domain suffix already on it, and the driver appends `.dat`. Neither is
something a new operator has to build.

For routes A and B that is all. `melem_interpolate_operator` is the generic driver and does
everything else **for any operator**: `O_W(k) = V^dagger O(k) V`, the distributed transform
to `R`, the transform back on the output domain, the diagonalisation of `H(k')` and the
projection. The only thing a new operator contributes is its component count and the name
in the table. A route-C name calls its own driver here instead — `CASE ('hamiltonian')` in
that same `SELECT CASE` is the shortest example of one.

### For the real-space export — `<operators_r><operator name=".."/>`

A row in `WANNIERLIB_OPR` and a branch in `postproc/melem_operators_r.F90`:

```fortran
CASE ('my_operator')
   CALL melem_op_rs_distributed(this, cell, kpts, vloc, coarse%myo0, gk_loc, 3, &
                                mpicm, irank, .FALSE., 'myfile.1')
```

`WANNIERLIB_OPR` has no output columns: these files carry a channel or a spin-block index
in the name, so the writer builds it.

### The schema

In `FleurInputSchema.xsd`, inside `WannierlibOperatorType` or `WannierlibOperatorRType`:

```xml
<xsd:enumeration value="my_operator"/>
```

Without it the input does not even validate. **Add the row and the enumeration in the same
commit.** The two lists have to agree and nothing enforces it; if they drift, an input
validates against the schema and then stops in `melem_request_init` with *"is not an
operator this layer can interpolate"*. The failure is clean and the message lists the
accepted names, so the cost is a confusing session rather than a wrong number — but the
schema is the first thing a new user reads, so it should not promise what does not exist.

---

## 7. Things that have actually cost us a run

| Symptom | Cause and fix |
|---|---|
| Error `#6401` at compile time | You named a dummy argument `mpi_comm` and it clashes with `use mpi`. Use `mpicm`, as the other routines do. |
| Error `#6404` on a constant | A missing `USE m_constants, ONLY: ...` (e.g. `tpi_const`, `ImagUnit`). |
| Error `#6632` on a routine you only moved | The destination module does not `USE m_judft`, so `juDFT_error(..., calledby=...)` has keyword arguments with no explicit interface. Check what the **destination** imports, not what the callers do. |
| **Works in serial; garbage or a hang with `np>1`** | You passed an unallocated `ALLOCATABLE` — allocated under `IF (irank == 0)` — to a **collective** routine with a non-allocatable dummy. That is undefined behaviour and a serial run can never show it; a two-rank suite can pass with it present. |
| A layer `l = lmax` vanishes silently | `lm = lmax**2` instead of `lmax*(lmax + 2)`. |
| Unexpected precision | In FLEUR a bare `real` **is double precision**. |
| The operator comes out tiny instead of failing | A stub slice reached the export without ever being computed. Small looks like numerical noise rather than like the absence of a calculation. |
| The link fails, unrelated to what you touched | The `fleurinput` CMakeLists lists files without the path prefix. |

Two about the build, kept here because nothing else in the tree records them: **never
compile on the login node**, and an internal compiler error from `ifx` has two causes that
one test separates — if it survives at `-O0` it is a stale build state and not the code, so
reconfigure from scratch before you go looking at what you wrote.

---

## 8. Validating it

### The generic invariants, for free

`melem_check_provider` runs three checks on any provider's result, at `k = 1`, and knows
nothing about which operator it was handed:

- **finite** — NaN or Inf. Reading eigenvector storage that was never written gives
  harmless zeros in serial but stale window memory under MPI-RMA.
- **non-zero** — a result that is identically zero, i.e. a slice never filled.
- **hermitian** — `|O - O^dagger| / max|O|`. This is the one that pays for the module:
  transpose the band indices, or put `CONJG` on the wrong factor, and the matrix stops
  being Hermitian while every other symptom is a plausible-looking number.

Wire your provider into the same place the three existing ones use and you get all three
without writing anything. Ask for hermiticity only of the slices that have it: a cross-spin
block does not, and `l_hermitian=.FALSE.` says so.

### The numerical criterion

> **Byte-identity at the same rank count; `Omega_I` across different rank counts.**

A change that does not intend to move numbers is validated by comparing byte for byte
against the previous output with the same number of processes. Between *different* rank
counts that does not apply: the reduction tree changes and the sums happen in another
order. `MKL_CBWR` and `I_MPI_CBWR` fix the code path, not the tree.

What must survive between `np=1` and `np=2` is `Omega_I`, which is gauge- and
basin-invariant. Without disentanglement it agrees to the last digit; with it, it moves,
and `Omega_total` always moves — a change of basin is legitimate and is not a bug. The
figures for particular cases are deliberately not quoted here: they are measurements of
one run each and they age faster than the criterion does.

A free trick falls out of this: running at two rank counts **classifies** the outputs,
because the ones that survive identical are by construction the gauge-invariant ones.

Which suite to run and where its inputs are: `README.md`, *The validation criterion*.

A new operator is **inert unless someone asks for it**: the slice is not allocated, the
provider is not called, nothing is written. The only real failure mode is forgetting the
gate, which is the stub in §3.5.

---

## 9. A full walkthrough: the orbital operator

Route A, end to end.

1. **Catalogue.** `t_melem_operator('orbital', 3)` — three Cartesian components.

2. **Provider.** `types_matelements_orbital.F90`. Both spinor flags `.FALSE.`, because L
   acts only on the spatial part: the cross-spin blocks vanish and the result is a single
   matrix. `n_alpha = 3`, because `L_x` and `L_y` both come out of `L_+` and `L_-` and
   therefore out of one pass.

3. **Bound to a site.** One instance covers *one* atom and *one* channel; whoever wants the
   total sums the sites. For L that sum is direct: it is spin-diagonal and its trace is
   invariant under the local-to-global rotation.

4. **Coarse slice.** `l0(nb,nb,3,nat,channel,nk_loc)` — the only slice with a channel index,
   because L is the quantity that exists in the same form for spinors and for two separate
   collinear channels. It is one array, not two.

5. **Exposure.** Two rows: `'orbital'` in `WANNIERLIB_INTERP` (with
   `honours_total = .TRUE.`) and in `WANNIERLIB_OPR`.

6. **Outputs.** `bands_wann_orbmom.dat` when interpolated; `anglmomrs.<channel>` in real
   space. Note the suffix carries the channel — that is why `WANNIERLIB_OPR` has no output
   column and the writer builds the name.
