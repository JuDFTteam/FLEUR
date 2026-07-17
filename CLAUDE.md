# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## What this is

FLEUR is an all-electron DFT code implementing the FLAPW method, written in Fortran and built with CMake. The primary repository is at iffgit.fz-juelich.de/fleur/fleur; `develop` is the default branch and small changes go directly there (keep it compiling and passing tests). User documentation lives at https://www.flapw.de.

## Build

```bash
./configure.sh              # detects compiler, creates build/ and runs cmake
cd build && make -j
```

- `./configure.sh -l mylabel` configures into `build.mylabel`; `-d` adds a debug build (`build.debug`). Run `./configure.sh -h` for all options (compiler selection, building bundled external libraries, etc.). The many `build.*` directories at the repo root are independent build trees.
- Main targets: `fleur` (serial), `fleur_MPI` (built when MPI is available), and the input generator `inpgen` (plus `inpgen2`/`inpgen3` and the `fleuriste` TUI under `src/tools/`).
- CMake logic lives in `cmake/` — source lists in `cmake/Files_and_Targets.txt`, compiler detection in `cmake/CompilerConfig.txt`. Optional externals (ELPA, ScaLAPACK, libxc, HDF5, MAGMA, chASE, EDsolver, Wannier90, ...) are in `external/` and toggle preprocessor definitions and test availability.
- `.F90` files are CPP-preprocessed (`CPP_MPI`, `CPP_HDF`, ...); `.f90` files are not. `fleur` and `fleur_MPI` are built from the same sources with different definitions.

## Tests

The test suite is pytest-based and lives in `testing/` (head `conftest.py` there). Each build directory gets a `run_tests.sh` wrapper that points pytest at the right build and forwards all arguments:

```bash
cd build
./run_tests.sh                    # full suite
./run_tests.sh -k CuBulk          # single test / substring match
./run_tests.sh -m bulk            # by marker (bulk, film, soc, lo, ldau, noco, forces, hybrid, fast, slow, ...)
./run_tests.sh -x                 # stop at first failure
```

Equivalently: `pytest <repo>/testing --build_dir=<builddir>`. Set `juDFT_PYTHON` to choose a non-default Python.

- Regression tests are defined declaratively in markdown tables in `testing/tests/parameterized/tests.md`; a `+` in the first column enables a test. Inputs live in `testing/inputfiles/<testset>/<name>`. To add a regression test, add the input directory and a table row.
- Test runs execute in `build/Testing/work`; failed-test artifacts are preserved in `build/Testing/failed_test_results`.
- Tests requiring unavailable libraries (elpa, magma, wannier, ...) are excluded automatically based on the build configuration.

## Lint / style

- Linter: `fortitude check` (config in `fortitude.toml`: line length 192, autofix enabled).
- pre-commit hooks (installed by `configure.sh` if `pre-commit` is available) enforce, among others: no `STOP` statements, `implicit none`, and copyright headers.
- From CONTRIBUTING.md:
  - 3-space indentation.
  - One module per file, module named `m_<name>` matching the file name; modules start with `implicit none` and `private`.
  - Never use `stop` — use `judft_error`, `judft_warn`, `judft_end` from the juDFT library.
  - Do not read/write files ad hoc; file IO belongs in dedicated IO routines or the standard FLEUR output files.
  - Pass arrays assumed-shape (`x(:,:)`) or allocatable, not explicit-shape (`x(n,m)`).

## Architecture

Three source layers under `src/`:

- `src/libraries/juDFT` — infrastructure used everywhere: error handling (`judft_error`/`judft_warn`/`judft_end`), timing (`timestart`/`timestop`), HDF5 wrappers, command-line arguments, logging.
- `src/libraries/fleurinput` — everything representing the `inp.xml` input file: one derived type per XML section (`t_atoms`, `t_cell`, `t_input`, `t_noco`, `t_kpts`, `t_sym`, `t_wannier`, ...), the XML readers, and MPI broadcast of the input. These are aggregated into `t_fleurinput`, passed around as `fi` (so `fi%atoms`, `fi%input`, ...) and treated as read-only after setup.
- `src/fleur/` — the actual physics code, organized by subsystem directory.

Program flow: `src/fleur/main/fleur.F90` contains the SCF loop (`scfloop:` DO WHILE). Each iteration: `vgen` (potential generation, `src/fleur/vgen/`) → `eigen` (Hamiltonian/overlap setup and diagonalization, `src/fleur/eigen/` + `src/fleur/diagonalization/`) → Fermi energy (`src/fleur/fermi/`) → `cdngen` (new charge density, `src/fleur/cdn*/`) → `totale` (total energy) → `mix` (density mixing, `src/fleur/mix/`). Setup happens in `fleur_init.F90`; `fleur_job.F90` drives job scheduling around it.

State is carried in derived types, not globals. Run-time types live in `src/fleur/types/` — the most pervasive are `t_potden` (potentials and densities share one type), `t_lapw` (basis set), `t_mat`/`t_mpimat` (dense matrices, serial vs. ScaLAPACK-distributed), `t_nococonv` (non-collinear magnetism state), and `t_results`. When adding data that flows through the SCF loop, extend or follow the pattern of these types.

`src/fleur/diagonalization/` abstracts over eigensolvers (LAPACK, ScaLAPACK, ELPA, MAGMA, chASE, cuSOLVER) selectable at run time via `-diag <solver>`.

A FLEUR calculation is driven by `inp.xml` (generated from a simple structure file by `inpgen`), writes `out.xml`/`out`, and stores the density in `cdn.hdf` (or `cdn*` files without HDF5). The XML schema is versioned and generated at build time (`cmake/Generate_Schema.cmake`).
