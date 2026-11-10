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

* `inputfiles/cu_supercell/` — a 4-atom Cu fcc supercell along (001), used by
  `test_gfprep.py`. The test converges it briefly with `fleur` and converts it with
  `gfprep` into four layers. It then checks the outputs and the accuracy of
  the conversion, and starts `gfleur` on the result. `gfprep`, `fleur` and
  `gfleur` are taken from the same build directory. Point to it with
  `GFPREP_BINARY=/path/to/build/gfprep`.

* `inputfiles/cu_film/` — a Cu(001) film of three atomic planes (RMT 1.6,
  Kmax 4.5 set by hand after inpgen), used by `test_gfprep_film.py`. The test
  converges it with `fleur`, converts it with `gfprep` into one or three
  layers with vacuum embedding planes on both sides (also with `zvac=1`), and
  compares the DOS peaks of `gfleur` at Gamma with the film eigenvalues of
  `fleur`. A surface run embeds two film planes between a bulk lead (CBS of
  `cu_supercell`) and the vacuum.

* `test_gfleur_transport.py` also starts from `cu_supercell` (ABAB stacking,
  two atomic planes per lattice period), for two MT radii: 2.35 (the
  embedding planes cut the MT spheres) and 1.6 (they do not). One period is
  cut as a single layer, and a CBS run (`&CBSMODE cbs=t`) writes the lead
  embedding potentials. The tests compare the propagating CBS states at Gamma
  with a `fleur` band run along (0,0,kz), check `spectral=f` against
  `spectral=t`, and require the transmission of one layer (`curr=-1`) and of
  a two-layer stack (`curr=3`, `curr=16`) to equal the number of open
  channels. At E=0.1 with RMT 2.35 the CBS tests are expected to fail
  (`xfail`) because of the cut spheres.

## Scope

Covered: the per-layer setup, gfprep for bulk and films (conversion
accuracy), T-matrices and Green functions with the gfprep potential, the
DOS with vacuum embedding, the CBS with the lead embedding potentials, and
the Landauer transmission. Not covered: self-consistent charge (`charge=t`),
the bias-dependent current, and MPI runs.

* `test_gfleur_noco.py` also covers spin-orbit coupling (first-variation SOC of
  the noco Hamiltonian, `l_soc="T"` with `l_noco="T"`): the Mn film with the
  moment along z (symmetry reduced to E and sigma_z, as FLEUR demands for SOC)
  reproduces the SOC splitting of a d-level pair at Gamma of the `fleur`
  eigenvalues, and sc Fe with SOC transmits every open channel. Collinear
  inputs with `l_soc="T"` (second variation) are rejected by gfprep and gfleur.

* `test_gfleur_charge.py` tests the charge mode (self-consistency): a stack of two Cu layers
  between bulk leads (fixed point from the bulk density, zero-mixing stationarity, Broyden
  convergence) and a semi-infinite Cu(001) surface (upper four planes of the 7-plane film in
  `inputfiles/cu_film7` between a shifted bulk lead and the frozen vacuum, converged with `imix=5`
  to the charges of the film). The module takes about twenty minutes.
