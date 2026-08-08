# `matrixelements` — adding an operator

Full tutorial with the reasoning: `August/refactor/tutorial_operadores.pdf`.
This file is the checklist.

## First: which of the three are you adding?

| | Examples | Work |
|---|---|---|
| **A real operator** — a contraction over the states, `O_mn(k) = <psi_m|O|psi_n>` | `spin`, `orbital`, `spin_orbit` | catalogue + provider + coarse slice + exposure |
| **Only an exposure** — the coarse matrix already exists | `spinCurrent`, `orbitalCurrent` | one table row + one branch |
| **Not an operator at all** — built from neighbour overlaps or eigenvalues | `hamiltonian`, `position`, `velocity`, `bmn`, `eigenstates` | row with `operator=''` + a driver of its own |

**The rule:** if it is not a contraction over the states, it does not go in the
catalogue. `hamiltonian` is `V^dagger diag(E) V`; `position` comes from the `M_mn`.

## The seven touchpoints (a real operator)

| # | File | What to add |
|---|---|---|
| 1 | `types_matelements_<name>.F90` | the provider. **Start from `types_matelements_template.F90`** |
| 2 | `matrixelements/CMakeLists.txt` | one line, **with** the path prefix |
| 3 | `fleurinput/types_melem_optable.f90` | a row in `MELEM_OPERATORS` |
| 4 | `melem_coarse.F90` | three places: the slice, the `ALLOCATE` behind `request%needs_op`, the fill inside the k loop |
| 5 | `melem_run.F90` and/or `melem_operators_r.F90` | a `CASE` branch |
| 6 | `types_melem_optable.f90` | a row in `WANNIERLIB_INTERP` and/or `WANNIERLIB_OPR` |
| 7 | `fleur/io/xml/FleurInputSchema.xsd` | one `<xsd:enumeration>` |

New on-disk format? Add a `CASE` in `melem_io.F90` — the only file that knows the
layout. Do not open a file anywhere else.

## Things that have actually cost us a run

- `mpi_comm` as a dummy name clashes with `use mpi` (#6401). Use `mpicm`.
- **Passing an unallocated `ALLOCATABLE` to a collective routine with a
  non-allocatable dummy is undefined behaviour, and a serial run can never show it.**
  If the callee is collective, every rank must reach it with a valid actual.
- Never allocate a `(1,1,1,1)` stub for something nobody asked for. A stub is a
  *valid* array of the wrong shape: indexing it wrongly does not fail.
- `lm = lmax*(lmax + 2)`, not `lmax**2`.
- Bare `real` in FLEUR is double precision.
- `fleurinput/CMakeLists.txt` lists files **without** the path prefix; this directory's
  lists them **with** it.
- Ask `request%needs_op('name')`, not the `l_spin`/`l_orbmom`/`l_socop` flags — those
  summarise both lists and cannot tell interpolation from real-space export.

## Validating it

**Byte-identity at the same rank count; `Omega_I` across different rank counts.**
Never byte-identity between different numbers of ranks: the reduction tree changes.
`MKL_CBWR`/`I_MPI_CBWR` fix the code path, not the tree.

- tests: `testing/tests/parameterized/test_wannier.py`
- inputs: `testing/inputfiles/wannier/`

## Why the template is compiled

`types_matelements_template.F90` is in `CMakeLists.txt` although nothing calls it.
A skeleton that is not compiled stops matching the interface it claims to
demonstrate and nothing says so — which is how `FleurInputSchema.xsd` came to accept
`spinTorque` and `orbitalTorque` while neither exposure table has them. Compiled, a
change to `t_matelements` breaks the build on the commit that makes it.

It costs one object file. The `CMakeLists.txt` entry is a commit of its own, so
dropping the pattern is a single revert — and if you do drop it, take the
"keep this file in the build" paragraph out of the template's header too.

## Reading order

1. `types_matelements_template.F90` — the skeleton, with the contract in the comments
2. `types_matelements_orbital.F90` — the same shape with a real operator in it
3. `types_matelements.F90` — the abstract type: `init_mat`, `add`, `full_matrix`
