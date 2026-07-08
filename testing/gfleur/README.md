# gfleur smoke tests

Light-weight, self-contained regression tests for the Green's-function
embedding executables `gfleur` / `gfleur_MPI`.

They are **not** wired into the main fleur pytest fixtures (which are built
around the `fleur`/`inpgen` binaries). Each test locates the gfleur binary,
stages a minimal input set into a temporary directory, runs it, and checks
that the whole initialisation + driver chain completes cleanly.

## Running

```
# from a build that produced the gfleur target
GFLEUR_BINARY=/path/to/build/gfleur pytest testing/gfleur -v
```

If `GFLEUR_BINARY` is unset the tests look for `<repo>/build*/gfleur`. When no
binary (or no `FleurInputSchema.xsd`) is found the tests skip, so the file is
safe to collect in a plain fleur build.

The MPI executable can be checked the same way, e.g.

```
GFLEUR_BINARY="$(pwd)/build.mpi/gfleur_MPI" pytest testing/gfleur -v
```
(run it under `mpirun` yourself for a genuine multi-rank check; the fixture
runs the binary directly, which is a single rank).

## Fixtures

* `inputfiles/2layer_film/` — two identical 1-atom Fe film layers
  (`layer1_inp.xml`, `layer2_inp.xml`), a minimal `gf_inp` (tmat mode,
  `iter = 0`) and a 3-point energy contour `gf_en`.

## Scope

Only the ported infrastructure is exercised (per-layer `fleur_init`, the
k×layer×energy parallel setup, step functions, the 2D-matrix/T-matrix stores
and `gfleur_execute`). The computational kernels (T-matrices, transmission,
self-consistent charge) need converged per-layer potentials and physically
reviewed reference data that are not part of this fixture; those belong in a
follow-up once seed potentials are available.
