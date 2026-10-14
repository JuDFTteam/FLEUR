# gfleur

`gfleur` is FLEUR's Green's-function embedding executable for layered,
film-like systems.  It initializes one FLEUR calculation for each layer and
then evaluates the embedding Green function on a shared two-dimensional
k-point mesh and energy contour.

The executable is built only when HDF5 is available.  A serial build produces
`gfleur`; an MPI build produces `gfleur_MPI`.  The `gfprep` executable is built
alongside them and prepares a layered calculation from a converged bulk
supercell.

## Preparing a calculation

The usual workflow is:

1. Run a converged bulk FLEUR calculation for a supercell stacked along the
   z direction.
2. In a separate working directory, create `gfprep_inp` and run `gfprep`.
   See [`../prep/README.md`](../prep/README.md) for its input and output.
3. Adjust the generated `gf_inp`, provide an energy grid in `gf_en`, and run
   `gfleur` in that working directory.

`gfprep` writes the layer input files and an initial `gf_inp` template.  Its
generated HDF5 files, including `gf_pot.hdf` and `gf_cdn.hdf`, are needed for
modes that use the converted potential or perform charge self-consistency.

## Required input files

`gfleur` reads the following files from its current directory:

* `gf_inp` contains the Green's-function-specific namelists and the layer
  geometry.
* `gf_en` defines the complex energy grid.
* `<prefix><i>_inp.xml` is the normal FLEUR input for layer `i`.  With the
  default prefix, these are `layer1_inp.xml`, `layer2_inp.xml`, and so on.
  The input files must describe film geometries.
* `FleurInputSchema.xsd` and `FleurOutputSchema.xsd` must be available in the
  run directory, as for a normal FLEUR calculation.

The layer XML files are the source of the electronic-structure settings:
basis parameters, local orbitals, LDA+U, non-collinearity, spin-orbit
coupling, and the exchange-correlation potential are configured there, not in
`gf_inp`.  All layers must have the same in-plane lattice, two-dimensional
k-point set, number of spin channels, and non-collinear setting.

## Minimal `gf_inp`

The following is sufficient to initialize a two-layer calculation without
self-consistency:

```fortran
&LAYERS nlayers=2 /
&MODE charge=f, dos=f, tmat=t, gmat=f /
```

`&LAYERS` is mandatory.  Its `prefix` defaults to `layer`; optional arrays
`c1`, `c2`, `d`, `dt`, and `c` describe each layer's embedding geometry.
`gfprep` generates these values.  If supplied, `d` and `dt` should agree with
the layer XML `dVac` and `dTilda`, respectively.

`&MODE` is also mandatory.  It selects the calculation type.  The most common
switches are `tmat`, `gmat`, `charge`, `dos`, and `curr`.  Charge generation
and conductance calculations require `addemb=t` in the optional `&EMB`
namelist.  Self-consistency is controlled with the optional `&SC` namelist;
for example, `iter` selects the number of iterations.  Keep `iter=0` when
only checking that a newly prepared setup initializes correctly.

## Energy grid

`gf_en` may contain explicitly listed energies or generated intervals.  For
example, a three-point list is:

```text
list 3
-0.10 0.01
 0.00 0.01
 0.10 0.01
```

Each energy is specified by real and imaginary parts in the units expected by
FLEUR.  `inter` and `circle` records can instead generate linear and
semicircular contours; formats may be combined in one file.

## Embedding-potential I/O: `gf_io`

`gf_io` is an optional text file that routes two-dimensional matrices to and
from HDF5 files.  It is required when `&EMB addemb=t` is used to read
external embedding potentials.  The matrix data are stored in the referenced
`*.hdf` files; `gf_io` only identifies the file, boundary, and access mode.

For `addemb=t`, provide exactly one readable conventional embedding potential
(`mode='E'`) for each lead boundary:

* the left boundary: `region=1`, `side=1`;
* the right boundary: `region=<nlayers>`, `side=2`.

The first line of `gf_io` is the number of `&IO` records.  It is followed by
one `&FILE` record for every file ID used by an `&IO` record.  For a two-layer
calculation, a minimal read-only configuration is:

```fortran
2
&IO mode='E', region=1, side=1, rw='r', id=1 /
&IO mode='E', region=2, side=2, rw='r', id=2 /

&FILE name='embpot11', energy_shift=0.0, kpts_eps=1e-5, en_eps=(1e-5,1e-5) /
&FILE name='embpotN2', energy_shift=0.0, kpts_eps=1e-5, en_eps=(1e-5,1e-5) /
```

This configuration reads `embpot11.hdf` and `embpotN2.hdf`.  The principal
`&IO` fields are:

| field | meaning |
|---|---|
| `mode` | matrix category; use `E` for a conventional embedding potential |
| `region`, `side` | layer/region and its boundary (`1` left, `2` right) |
| `rw` | access mode: `r` read, `w` write, or `b` read/write |
| `id` | one-based index of the associated `&FILE` record |
| `spinrot`, `theta`, `phi` | optional spin rotation applied on read |
| `shift` | optional in-plane interface translation |

`&FILE name` is a base name; `gfleur` adds the `.hdf` suffix.  `kpts_eps` and
`en_eps` control the tolerances for matching the stored k-points and complex
energy grid.  HDF5 files contain the matrix data together with their energy
grid, in-plane k-points, reciprocal vectors, and symmetry metadata.

The special first line `default` causes `gfleur` to replace `gf_io` with a
two-boundary default that reads `embpot11.hdf` and `embpotN2.hdf`.  Do not use
it for an existing custom configuration because it overwrites the file.

### Practical requirements for T-matrix embedding potentials

`gfleur` obtains an interface embedding potential by propagating a known lead
boundary potential through the T matrices of a stack of layers.  To produce
HDF5 embedding-potential files from such a calculation:

1. Use a multi-layer setup (`nlayers > 1`) so that layer T matrices are
   constructed.
2. Enable T-matrix work with `tmat=t` in `&MODE`.
3. Configure readable `mode='E'` `gf_io` handles for the initial left and
   right lead boundary potentials.
4. Configure writable `mode='E'` handles for every layer/interface and side
   where the propagated potential should be saved.
5. Run a calculation mode that reaches the embedding-potential propagation
   stage.

The saved embedding potential can only be reused with compatible
embedding-plane geometry, in-plane lattice, two-dimensional reciprocal-vector
basis, k-point mesh, spin treatment, and energy grid.

## Running

Run the executable from the directory containing the files above:

```sh
/path/to/gfleur
```

For an MPI build, use the launcher appropriate to the local system:

```sh
mpirun -n 4 /path/to/gfleur_MPI
```

The setup writes its main log to `out`, reports `GF-SETUP done` when layer
initialization succeeds, and creates `gf_steps.hdf` for the analytical layer
step functions.  A setup-only run with `iter=0` is a useful first check before
running computational modes that require converged layer potentials.
