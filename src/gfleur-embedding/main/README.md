# gfleur

`gfleur` is FLEUR's Green's-function embedding executable for layered,
film-like systems.  It initializes one FLEUR calculation for each layer and
then evaluates the embedding Green function on a shared two-dimensional
k-point mesh and energy contour.

The executable is built only when HDF5 is available.  A serial build produces
`gfleur`; an MPI build produces `gfleur_MPI`.  The `gfprep` executable is built
alongside them and prepares a layered calculation from a converged bulk
supercell or film.

## Preparing a calculation

The usual workflow is:

1. Run a converged bulk FLEUR calculation for a supercell stacked along the
   z direction, or a film calculation.
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
basis parameters, local orbitals and the exchange-correlation potential are
configured there, not in `gf_inp`.  (LDA+U and spin-orbit coupling are
not supported by the Green's-function kernels, see below.)  All layers must have the same in-plane lattice, two-dimensional
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

## Basis of the layer Green function

The interstitial of a layer is described by plane waves of period `dt`
along z. `&BASIS napw` (default 15, one value per layer) sets their number;
thicker layers need more. With `&MODE spectral=t` (default) the Green
function is built from the eigenstates of the layer Hamiltonian. These plane
waves are an overcomplete set on the layer, so the overlap matrix becomes
ill-conditioned as `napw` grows. With `&MODE s_eps=<x>` (default 0, off),
eigenvectors of the overlap with eigenvalues below `x` times the largest one
are dropped before the eigenvalue problem is solved, e.g. `s_eps=1e-10`.

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

A boundary with vacuum embedding (`surface`, `vacleft`, see below) has no
lead and must not have a read handle.

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

The special first line `default` selects a two-boundary default that reads
`embpot11.hdf` (region 1, side 1) and `embpotN2.hdf` (region `<nlayers>`,
side 2).  In complex-band-structure mode (`&CBSMODE cbs=t`) the same default
writes these two files instead, so a CBS run of a bulk-like layer produces the
lead embedding potentials.  A side with vacuum embedding is left out of the
default.  The defaults are built in memory, `gf_io` itself is never modified.

Without a `gf_io` file, runs with `nlayers > 1` keep the layer T-matrices in
memory; `addemb=t` requires `gf_io` unless both sides have vacuum embedding.

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

## Lead embedding potentials from a complex band structure

The embedding potentials of the left and right lead come from a complex band
structure (CBS) calculation of one lattice period of the lead, run as a
separate `gfleur` calculation with `nlayers=1` (one period = one layer, cut
with `gfprep ... planes_per_layer=<planes of the period>`):

```fortran
&MODE tmat=t, gmat=t /
&CBSMODE cbs=t /
&SC iter=1 /
```

together with a `gf_io` whose first line is `default`, and a `gf_en` with the
same energies as the later calculation.  The run writes `embpot11.hdf` (left
lead) and `embpotN2.hdf` (right lead), the CBS (`CBS_eigenvalues.<spin>`,
`CBS_bloch.<spin>`, or `CBS.hdf` with `hdfio=t`) and, in a serial run with
`gf_cdn.hdf` present, the `Boundary` group needed for charge calculations.
Notes:

* Use the *primitive* period as lead.  The transfer matrix of a thick period
  loses the evanescent states numerically and gives inaccurate embedding
  potentials (a four-plane Cu period gave wrong stacks, the two-plane
  period is exact).
* Layers thinner than the muffin-tin diameter (one Cu plane, 3.4 a.u.) are
  less accurate than thicker ones.
* The lead and the stack must be prepared with the same bulk calculation,
  Gmax scale and `kgrid`, and use the same `gf_en`.

## Transport

Transport needs `addemb=t`, `iter=1` and a `gf_io` with two read handles for
the leads, left at `region=1, side=1`, right at `region=<nlayers>, side=2`:

```fortran
&MODE gmat=t, curr=3 /
&EMB addemb=t /
&SC iter=1 /
```

`curr=-1` (Landauer from the combined T matrix) is exact for `nlayers=1`; for
stacks it multiplies the raw layer T matrices, which loses accuracy, so use
`curr=3`, `4` or `16` (propagation of the embedding potential, `16` is the
cheapest) for `nlayers>1`.

`transcurrent` (`transcurrent.hdf` with `hdfio=t`) lists the transmission
for every k of the full 2D Brillouin zone, columns `spin, kx, ky, E, T...`.
Both edges of the grid (`k=-0.5` and `0.5`) are listed, so the conductance
`G(E)` in units of `e^2/h` is the average of `T` over the `k` points of one
energy with a weight of 1/2 for each coordinate that equals ±0.5 (a point
with both coordinates on the edge has weight 1/4).

For a perfect crystal, `T` at every k equals the number of right-moving Bloch
states of the CBS, which is a useful test.

## Vacuum embedding (films and surfaces)

A layer stack cut out of a FLEUR film by `gfprep` ends at vacuum embedding
planes. `&EMB surface=t` embeds the upper side of the last layer in the
vacuum, `vacleft=t` the lower side of layer 1. A side with vacuum embedding
needs no lead file: with the `default` `gf_io` only the remaining leads are
read, and without `gf_io` a free film (both sides vacuum) runs with
`addemb=t`. gfprep writes the `&EMB` line into the `gf_inp` template:

```fortran
&EMB addemb=T, surface=T, vacleft=T /
```

The embedding potential is diagonal in the 2D plane waves,
`Sigma_GG = 1/2 phi'/phi` at the plane, with `phi` the solution of the
planar-averaged potential beyond the plane that decays into the vacuum (or is
an outgoing wave above the vacuum level). The potential comes from
`gf_pot.hdf` (groups `vacemb_side1`, `vacemb_side2`, written by `gfprep`). The
warping of the potential beyond the plane is neglected; `zvac` in `gfprep`
moves the planes further out. With `simplevacuum=t` the potential beyond the
plane is the constant `vacuum_energy`.

The vacuum embedding is available in the `tmat`/`gmat`/`dos` modes, not in
`charge` (self-consistent surfaces) and `CBS` mode.

## Non-collinear magnetism

Layers with `l_noco="T"` (from a noco bulk or film via `gfprep`) are treated
in one spinor pass. The CBS, the layer and embedded DOS, the transmission
and the vacuum embedding work as for collinear layers. In the DOS files the
third spin is the off-diagonal (up-down) part; the MT DOS is given in the
local spin frame of each atom type. Not supported: spin spirals,
spin-orbit coupling with noco, and the `charge` (SCF) mode.

## Not supported

* Self-consistent surfaces (`charge=t` with vacuum embedding).
* Noco: spin spirals, spin-orbit coupling and the `charge` mode.
* Spin-orbit coupling, LDA+U and Hubbard-1 are written into the layer inputs
  by `gfprep` but are not used by the Green's-function kernels.

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
