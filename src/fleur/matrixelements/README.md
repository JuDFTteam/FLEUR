# `matrixelements` — adding an operator

This file is the checklist. `tutorial_operators.md`, next to it, carries what no single
file can: which of three routes you are on, the sequences that cross files, and the
mistakes. The contract of each routine is in that routine's own header.

## Where things live

This directory holds what computes matrix elements: the providers, the factory that
reads the states and their coefficients, the neighbour overlaps, and the coarse pass
that produces `O(k)` on the ab-initio eigenstates. Nothing here knows what a Wannier
function is.

Everything that needs the Wannier gauge `V = u_opt . u_matrix` — the rotation, the
Fourier transform to `R`, the band interpolation, the `O(R)` export and its file
formats — lives in `../wannierlib/postproc/`.

The seam between the two is two types, `t_melem_window` (which bands were selected)
and `t_melem_request` (which operators were asked for), both of which stay here.
`t_melem_manifold` extends the window with `num_wann` and the disentanglement edges
and lives in `postproc`, so nothing here knows what the bands were selected for.
`postproc` depends on this directory; this directory depends on nothing of it, and
`testing/tests/structure/test_layering.py` fails if that ever changes.

## First: which of the three are you adding?

- **A real operator** — a contraction over the states, `O_mn(k) = <psi_m|O|psi_n>`
  (`spin`, `orbital`, `spin_orbit`). Route A, the seven touchpoints below.
- **Only an exposure** — the coarse matrix already exists (`soc` reuses `spin_orbit`).
  Route B: one table row and one branch.
- **Not an operator at all** — built from the neighbour overlaps or the eigenvalues
  (`hamiltonian`, `velocity`, `eigenstates`, `position`, `position_pw90`, `bmn`, `fmn`,
  `cmn`). Route C: a row with `operator=''` and a driver of your own.

**The rule:** if it is not a contraction over the states, it does not go in the
catalogue. Route C is the majority — eight of the twelve exposed names — and the tutorial's §5 is
the one to read before writing a line of it.

## The seven touchpoints (route A)

| # | File | What to add |
|---|---|---|
| 1 | `types_matelements_<name>.F90` | the provider. **Start from `types_matelements_template.F90`** |
| 2 | `matrixelements/CMakeLists.txt` | one line, **with** the path prefix |
| 3 | `fleurinput/types_melem_optable.f90` | a row in `MELEM_OPERATORS` |
| 4 | `melem_coarse.F90` | three places: the slice, the `ALLOCATE` behind `request%needs_op`, the fill inside the k loop |
| 5 | `../wannierlib/postproc/melem_run.F90` and/or `melem_operators_r.F90` | a `CASE` branch |
| 6 | `types_melem_optable.f90` | a row in `WANNIERLIB_INTERP` and/or `WANNIERLIB_OPR` |
| 7 | `fleur/io/xml/FleurInputSchema.xsd` | one `<xsd:enumeration>` |

Route B is 5, 6 and 7. Route C is 5, 6 and 7 plus the driver itself — and, only if it
needs the wavefunctions after the gauge is known, one `IF` in `wannierlib_main.F90`.

New on-disk format? Add a `CASE` in `../wannierlib/postproc/melem_io.F90` — the only
file that knows the layout. Do not open a file anywhere else.

## Before you debug

The symptom-to-cause table is in `tutorial_operators.md` §7. The two that cost the most
runs, in short: a dummy argument named `mpi_comm` clashes with `use mpi` (`#6401`, use
`mpicm`), and handing an unallocated `ALLOCATABLE` to a **collective** routine is
undefined behaviour that a serial run — and a two-rank suite — can pass without showing.

## The validation criterion

**Byte-identity at the same rank count; `Omega_I` across different rank counts.**
Never byte-identity between different numbers of ranks: the reduction tree changes.
`MKL_CBWR`/`I_MPI_CBWR` fix the code path, not the tree.

- tests: `testing/tests/parameterized/test_wannier.py`
- inputs: `testing/inputfiles/wannier/`
- layering: `testing/tests/structure/test_layering.py`

## Why the template is compiled

`types_matelements_template.F90` is in `CMakeLists.txt` although nothing calls it.
A skeleton that is not compiled stops matching the interface it claims to
demonstrate and nothing says so — which is how `FleurInputSchema.xsd` once accepted
two operator names that neither exposure table had, until someone noticed by hand.
Compiled, a change to `t_matelements` breaks the build on the commit that makes it.

It costs one object file. The `CMakeLists.txt` entry is a commit of its own, so
dropping the pattern is a single revert — and if you do drop it, take the
"keep this file in the build" paragraph out of the template's header too.

## Reading order

1. `types_matelements_template.F90` — the skeleton, with the contract in the comments
2. `types_matelements_orbital.F90` — the same shape with a real operator in it
3. `types_matelements.F90` — the abstract type: `init_mat`, `add`, `full_matrix`
