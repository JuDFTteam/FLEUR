# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## What FLEUR Is

FLEUR is an all-electron DFT code implementing the full-potential linearized augmented plane-wave (FLAPW) method. It is a scientific HPC code written primarily in Fortran 90/95+, with Python tooling for input generation and testing. The primary executables are `fleur` (serial) and `fleur_MPI` (parallel), plus `inpgen`/`inpgen3` for input file generation.

## Building

FLEUR uses CMake via a `configure.sh` wrapper:

```bash
# Standard build (creates ./build/ directory)
./configure.sh

# Build with a label (creates ./build.<label>/)
./configure.sh <label>

# Debug build
./configure.sh -debug

# After configure, build from the build directory
cd build && make -j
```

The configure script auto-detects compilers (GFortran, Intel OneAPI, NVHPC) and optional dependencies (MPI, HDF5, ELSI, ELPA, libxc, ScaLAPACK). CMake configuration output goes to `configure.out`.

Useful configure flags not shown above:
- `-ninja` — use Ninja instead of GNU make
- `-make` — run `make` immediately after configure
- `-debug` — debug build
- `-gpu <spec>` — GPU support (e.g., `acc:cc80` for OpenACC on Compute Capability 8.0)
- `-elsi`, `-elpa`, `-chase`, `-magma` — optional diagonalization libraries (auto-downloaded if TRUE)
- `-libxc`, `-wannier` — enable libxc / Wannier interface
- `-libdir <path>`, `-includedir <path>`, `-link <flag>` — pass custom library paths/flags

Pre-configured machine setups for HPC clusters (JURECA, HAWK, SUPERMUC, CLAIX, etc.) live in `cmake/machines/`.

### Double-precision compilation

FLEUR source code uses plain, undecorated `real` declarations throughout (not explicit `real(kind=dp)` typing). Double precision is enforced by a compiler flag chosen per-compiler in `cmake/tests/test_precision.cmake`, applied per-target (not via global `CMAKE_Fortran_FLAGS`):

| Compiler | Flag |
|---|---|
| Intel / NAG | `-r8` |
| GNU (gfortran) | `-fdefault-real-8;-fdefault-double-8` |
| NVHPC / PGI | `-Mr8;-Mr8intrinsics` |
| IBM XL | `-qrealsize=8` |
| Cray | `-s real64` |

A normal FLEUR build always applies this flag — do not assume plain `real` means single precision when reading source.

## Running Tests

Tests use pytest and require a built FLEUR. The easiest way from the build directory:

```bash
./run_tests.sh                        # all tests
./run_tests.sh -k <substring>         # tests matching name substring
./run_tests.sh -m <marker>            # tests with a specific marker
./run_tests.sh -x                     # stop at first failure
./run_tests.sh testing/tests/feature_reg/test_CuBulk.py  # single file
```

Alternatively, from `testing/` with a non-default build dir:

```bash
cd testing
pytest --build_dir=../build.123
pytest tests/feature_reg/test_CuBulk.py --build_dir=../build
```

Common test markers: `bulk`, `film`, `collinear`, `non_collinear`, `soc`, `noco`, `forces`, `hybrid`, `dfpt`, `greensfunction`, `fast`, `slow`, `very_slow`, `spinspiral`, `lo`, `ldau`, `xml`, `mpi`, `serial`.

Additional pytest options:
- `--no-cleanup` — preserve all test work dirs (default: passing dirs deleted)
- `--runevery N` / `--testoffset M` — run every Nth test, skipping first M (useful for CI sampling)
- `--skipmarkers <list>` — comma-separated markers to skip

The build step generates `pytest_incl.py` in the build dir with marker exclusions based on which libraries were linked (e.g., libxc/ELPA tests auto-skipped if those weren't built). Failed test directories are preserved in `build/Testing/failed_test_results/`.

Set `juDFT_PYTHON` to override the Python interpreter used by `run_tests.sh`.

## Source Layout

```
src/
  fleur/           # Core DFT engine (Fortran)
    cdn/           # Charge density in interstitial region
    cdn_mt/        # Charge density in muffin-tin spheres
    core/          # Core electron calculations
    dfpt/          # Density functional perturbation theory (phonons)
    diagonalization/ # Eigenvalue solvers
    eigen/         # Hamiltonian/overlap matrix construction
    eigen_soc/     # Spin-orbit coupling
    fft/           # Fast Fourier transforms
    force/         # Hellmann-Feynman forces
    greensf/       # Green's functions
    hybrid/        # Hybrid functionals
    init/          # Initialization routines
    io/            # XML and HDF5 I/O
    main/          # Top-level SCF loop
    mix/           # SCF density mixing
    mpi/           # MPI parallelization
    types/         # Derived type definitions
    vgen/          # Potential generation
    wannier/       # Wannier function interface
    global/        # Global variables/parameters
    math/          # Mathematical utilities
  libraries/
    juDFT/         # Error handling, timing, MPI wrappers
    fleurinput/    # XML input parsing
  tools/
    inpgen2/       # Legacy input generator
    inpgen3/       # Current input generator (Fortran + Python)
    inpgen3/fleuriste/  # TUI/CLI parallelization helper (Python)
testing/
  tests/           # pytest test files (feature_reg/, inpgen/, masci_tools/, libxc/)
  inputfiles/      # Input files for tests
  helpers/         # Shared Python test utilities
  conftest.py      # Pytest fixtures and configuration
cmake/             # CMake modules and build configuration
external/          # Git submodules (libxc, HDF5, ELSI, ELPA, SCALAPACK, etc.)
```

## Coding Conventions (Fortran)

- **Indentation:** 3 spaces (no tabs)
- **Module naming:** prefix with `m_` (e.g., `m_sorad`, `m_types_fleur`)
- **File/module correspondence:** one module per file, names must match
- **Every module:** starts with `implicit none` and `private`
- **Error handling:** use `judft_error()`, `judft_warn()`, `judft_end()` — never `stop`
- **Array arguments:** use shape-assumed `real, intent(in) :: x(:,:)` or allocatable arrays; avoid explicit-size `real, intent(in) :: x(n,m)` which allows unsafe rank/size reinterpretation
- **No file I/O outside `io/`:** files are not substitutes for common blocks or status variables

Pre-commit hooks (`.pre-commit-config.yaml`) enforce: copyright header presence, `implicit none`, absence of `stop` statements, and validity of XML/YAML/TOML files and check for added large files. Install them with:

```bash
pre-commit install
```

For Fortran static analysis, `fortitude` is configured via `fortitude.toml` at the repo root (`fix = true`, line-length 192, ignores C003/C121/C131):

```bash
fortitude check src/        # check only
fortitude check --fix src/  # check and auto-fix
```

## Architecture Notes

The SCF loop lives in `src/fleur/main/`. Each iteration calls charge density generation (`cdn/`, `cdn_mt/`), potential generation (`vgen/`), Hamiltonian construction and diagonalization (`eigen/`, `diagonalization/`), and mixing (`mix/`). Parallelization is MPI-based via wrappers in `mpi/` and `libraries/juDFT/`.

Input/output uses XML (`inp.xml`, `out.xml`) parsed via `fleurinput/` and optionally HDF5 for eigenvalue storage. The `types/` directory contains the central derived types that are threaded through most subroutines.

DFPT (phonon calculations) lives in `src/fleur/dfpt/` and is a major subsystem with its own driver.

### Key derived types

Input-describing types live in `src/libraries/fleurinput/` and are aggregated in `t_fleurinput` (types_fleurinput.f90):
- `t_atoms` — atom types and sites, LDA+U parameters (types_atoms.F90)
- `t_cell` — Bravais matrix (amat), reciprocal matrix (bmat), cell volume (types_cell.f90)
- `t_sym` — symmetry operations (types_sym.f90)
- `t_kpts` — k-point set (types_kpts.F90)
- `t_noco` — non-collinear magnetism settings (types_noco.f90)
- `t_xcpot` — XC functional selection (types_xcpot.F90)

Calculation-state types live in `src/fleur/types/` and are collected in `m_types` (types.F90):
- `t_lapw` — LAPW basis (G-vectors, k-points) (types_lapw.F90)
- `t_potden` — potential and density arrays (interstitial + muffin-tin) (types_potden.F90)
- `t_mat` — Hamiltonian/overlap matrix (dense or MPI-distributed) (types_mat.F90)
- `t_abc` — LAPW matching coefficients A, B, C (types_abc.F90)
- `t_mpi` — MPI communicator and rank info (types_mpi.F90)

All input types extend `t_fleurinput_base` and expose an `mpi_bc` method for broadcasting over MPI.

### juDFT utilities

`src/libraries/juDFT/` provides more than error handling:
- **Timing:** `call timestart("label")` / `call timestop("label")` — produces a timing report; integrates with NVTX for NVIDIA profilers
- **HDF5 wrappers:** `hdf_tools*.F90` — high-level array read/write, attributes, stride I/O
- **NumPy interop:** `npy.F90` — read/write `.npy` files directly from Fortran
- **String utilities:** `string.f90` — `int2str` and other helpers
- **Argument parsing:** `args.F90`, `check_arguments.F90` — CLI flag extraction

### inpgen3

`src/tools/inpgen3/` operates in two modes:
1. Generate `inp.xml` from a simple input file: `inpgen -f <input>`
2. Add a k-point set to an existing `inp.xml`: `inpgen -kpt ...`

`FleurInpgen.py` provides a ctypes interface to the inpgen Fortran library; it searches standard paths and respects the `FLEUR_BUILDDIR` environment variable for library discovery. The `fleuriste/` subdirectory is a TUI/CLI tool for setting up parallelized FLEUR runs.

## Environment Variables

- `juDFT_MPI` — Custom MPI command template (use `{mpi_procs}` as placeholder)
- `juDFT_NPROCS` — Override number of MPI processes for tests
- `juDFT_PYTHON` — Override the Python interpreter used by `run_tests.sh`
- `juDFT_ARGS` — Extra command-line arguments passed to FLEUR
- `FLEUR_BUILDDIR` — Build directory for `FleurInpgen.py` ctypes library discovery

## Writing Tests

Tests in `testing/tests/` use fixtures from `testing/conftest.py`. Key helpers:

- `execute_fleur(...)` / `execute_inpgen(...)` — copy input files and run the binary
- `grep_exists(filepath, expression)` — assert a regex appears in an output file
- `grep_number(filepath, expression, ...)` — extract a float value from an output file
- `check_outxml(filepath, ...)` — compare values in `out.xml` against expected numbers
- `check_hdf(...)` — compare HDF5 output via `h5diff`

Test inputs live in `testing/inputfiles/` organised by feature (e.g. `noco/`, `forces/`, `dfpt/`, `hybrid/`). Multi-step tests use `stage1/`, `stage2/`, … subdirectories.

`out.xml` is validated against `src/fleur/io/xml/FleurOutputSchema.xsd`; `inp.xml` against `FleurInputSchema.xsd` in the same directory.

## Git Workflow

- `develop` — default branch; small changes and bugfixes go here directly
- `stable` — snapshots of develop considered stable (no direct commits)
- Feature branches — for larger developments; merge `develop` frequently
- CI pipeline on GitLab (`iffgit.fz-juelich.de/fleur/fleur`) must pass before merging; check `.gitlab-ci.yml` for pipeline details
- Bugs in the `release` branch require a dedicated bugfix branch + merge request into `release`, then merge the fix back into `develop`
