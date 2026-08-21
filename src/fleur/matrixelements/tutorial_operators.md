# Adding an operator to `matrixelements`

The long form. `README.md` next to this file is the checklist; this one carries the
reasoning, the contract of each routine you have to write, and a worked example.

Files are named by symbol rather than by line number, because line numbers rot and
symbol names do not.

---

## 1. First: which of three things are you adding?

This is the decision that saves the most work, and the one the tables in the code try to
make you take explicitly.

| What it is | Examples already in the tree | What you touch |
|---|---|---|
| **A real operator** — a contraction over the states, `O_mn(k) = <psi_m\|O\|psi_n>` | `spin`, `orbital`, `spin_orbit` | catalogue + provider + coarse slice + exposure (**route A**) |
| **Only an exposure** — the coarse matrix already exists; you want another file or a derived quantity | `soc` (reuses `spin_orbit` under another name) | one table row + one branch (**route B**) |
| **Not an operator at all** — it comes out of the neighbour overlaps or the eigenvalues | `hamiltonian`, `position`, `velocity`, `bmn`, `eigenstates` | a row with `operator=''` + a driver of its own (**route C**) |

> **The rule in one line.** If it is not a contraction over the states, it does not belong
> in the catalogue. `hamiltonian` is not an operator — it is `V^dagger diag(E) V`.
> `position` is not either — it comes out of the `M_mn`. Putting them in the catalogue
> would force the coarse pass to produce something nobody reads.

---

## 2. The map: two tables, and how a name travels

All the routing lives in one file, `fleurinput/types_melem_optable.f90`, and it holds
**two** tables with different owners.

```fortran
!> WHAT CAN BE COMPUTED. Three entries, and they know nothing about Wannier functions.
TYPE(t_melem_operator), PARAMETER, PUBLIC :: MELEM_OPERATORS(3) = [ &
   t_melem_operator('spin',       3), &
   t_melem_operator('orbital',    3), &
   t_melem_operator('spin_orbit', 4)]

!> HOW THE WANNIERISATION EXPOSES THINGS
!                   name          operator     ncomp  out1                out2  total
TYPE(t_melem_exposed), PARAMETER, PUBLIC :: WANNIERLIB_INTERP(8) = [ &
   t_melem_exposed('spin',        'spin',       3, 'bands_wann_spin',   '', .TRUE.), &
   ... ]
TYPE(t_melem_exposed), PARAMETER, PUBLIC :: WANNIERLIB_OPR(6) = [ ... ]
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

### Where the files are

| Directory | Holds |
|---|---|
| `matrixelements/` | the providers, the factory, the neighbour overlaps, the coarse pass, the generic invariants |
| `../wannierlib/postproc/` | everything that needs the Wannier gauge: the rotation, the transform to `R`, the interpolation, the `O(R)` export, the file formats |
| `fleurinput/types_melem_optable.f90` | the two tables |
| `fleur/io/xml/FleurInputSchema.xsd` | which names the input accepts |

`postproc` depends on `matrixelements`; the dependency never runs the other way.

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

---

## 4. Exposing it

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
your `WANNIERLIB_INTERP` row with the domain suffix already on it, and the driver appends
`.dat`. Neither is something a new operator has to build.

That is all. `melem_interpolate_operator` is the generic driver and does everything else
**for any operator**: `O_W(k) = V^dagger O(k) V`, the distributed transform to `R`, the
transform back on the output domain, the diagonalisation of `H(k')` and the projection.
The only thing a new operator contributes is its component count and the name in the table.

### For the real-space export — `<operators_r><operator name=".."/>`

A row in `WANNIERLIB_OPR` and a branch in `postproc/melem_operators_r.F90`:

```fortran
CASE ('my_operator')
   CALL melem_op_rs_distributed(this, cell, kpts, vloc, coarse%myo0, gk_loc, 3, &
                                mpicm, irank, .FALSE., 'myfile.1')
```

If the on-disk layout is not one of the seven that exist (`hr`, `r`, `rti`, `bmn`, `soc`,
`generic`, `cart2`), add a `CASE` in `postproc/melem_io.F90` — **the only file that knows
any layout.** Do not open a file anywhere else.

### The schema

In `FleurInputSchema.xsd`, inside `WannierlibOperatorType` or `WannierlibOperatorRType`:

```xml
<xsd:enumeration value="my_operator"/>
```

Without it the input does not even validate.

---

## 5. Route B: only an exposure

When the coarse matrix already exists and what you want is another quantity derived from
it, there is no provider and no new slice. One row pointing at the existing catalogue entry
is enough:

```fortran
t_melem_exposed('soc', 'spin_orbit', 1, 'bands_wann_soc', '', .FALSE.)
```

That row says three things at once: the name is accepted in `<interpolation>`, serving it
requires the `spin_orbit` matrix — even though the user never writes that word — and the result
has three components. The reader, the coarse pass and the dispatch all read it from there.

| Column | Meaning |
|---|---|
| `name` | what the user writes in `inp.xml` |
| `operator` | the catalogue entry to build; `''` if none |
| `ncomp` | components the generic driver writes; `0` = it has a driver of its own |
| `out1`, `out2` | output basenames, without `.dat` and without the domain suffix |
| `honours_total` | whether it respects `total="T"`, i.e. has a site-summed projection |

`out1` and `out2` are not decorative: the per-domain renaming (`_plane`, `_grid`) reads them
from the table. An operator that leaves them empty writes the same file every time and the
second domain overwrites the first.

---

## 6. Things that have actually cost us a run

| Symptom | Cause and fix |
|---|---|
| Error `#6401` at compile time | You named a dummy argument `mpi_comm` and it clashes with `use mpi`. Use `mpicm`, as the other routines do. |
| Error `#6404` on a constant | A missing `USE m_constants, ONLY: ...` (e.g. `tpi_const`, `ImagUnit`). |
| **Works in serial; garbage or a hang with `np>1`** | You passed an unallocated `ALLOCATABLE` — allocated under `IF (irank == 0)` — to a **collective** routine with a non-allocatable dummy. That is undefined behaviour and a serial run can never show it. |
| A layer `l = lmax` vanishes silently | `lm = lmax**2` instead of `lmax*(lmax + 2)`. |
| Unexpected precision | In FLEUR a bare `real` **is double precision**. |
| The operator comes out tiny instead of failing | A stub slice reached the export without ever being computed. Small looks like numerical noise rather than like the absence of a calculation. |
| The link fails, unrelated to what you touched | The `fleurinput` CMakeLists lists files without the path prefix. |

Two more, about the build itself: **never compile on the login node**, and when the
compiler reports an internal error, suspect a stale build state before you look for the
cause in the code — configure from scratch first.

---

## 7. Validating it

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
without writing anything.

### The numerical criterion

> **Byte-identity at the same rank count; `Omega_I` across different rank counts.**

A change that does not intend to move numbers is validated by comparing byte for byte
against the previous output with the same number of processes. Between *different* rank
counts that does not apply: the reduction tree changes and the sums happen in another
order. `MKL_CBWR` and `I_MPI_CBWR` fix the code path, not the tree.

What must survive between `np=1` and `np=2` is `Omega_I`, which is gauge- and
basin-invariant. Measured over the current cases:

| Case | bands → WFs | `Omega_I` rel. | `Omega_total` rel. |
|---|---|---|---|
| `WannFeBcc` | 6 → 6 | 0 | 3.1e-2 |
| `WannFeBccInterp` | 6 → 6 | 0 | 3.1e-2 |
| `WannPtSOCOps` | 36 → 18 | 4.6e-6 | 2.3e-3 |
| `WannFeAFMSOCOps` | 72 → 36 | 2.9e-4 | 2.3e-3 |

Without disentanglement `Omega_I` agrees to the last digit. With it, it moves, and
`Omega_total` always moves: a change of basin is legitimate and is not a bug. A free trick
falls out of this — running at two rank counts **classifies** the outputs, because the ones
that survive identical are by construction the gauge-invariant ones.

- tests: `testing/tests/parameterized/test_wannier.py`
- inputs: `testing/inputfiles/wannier/`

### It cannot break what already works

A new operator is **inert unless someone asks for it**: the slice is not allocated, the
provider is not called, nothing is written. The only real failure mode is forgetting the
gate, which is the stub above.

---

## 8. A full walkthrough: the orbital operator

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

---

## Appendix: a name the schema accepts and the tables do not

The enumeration in the XSD and `WANNIERLIB_INTERP` are two lists that have to agree, and
nothing enforces it. If they drift, an input validates against the schema and then stops in
`melem_request_init` with *"is not an operator this layer can interpolate"*. The failure is
clean and the message lists the accepted names, so the cost is a confusing session rather
than a wrong number — but the schema is the first thing a new user reads, so it should not
promise what does not exist. Add the row and the enumeration in the same commit.
