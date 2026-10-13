# gfprep

`gfprep` turns a converged FLEUR bulk supercell, stacked along z, into the
layered setup used by `gfleur`. It is built with the gfleur targets (HDF5 is
required).

## Usage

Run it in an empty directory, next to the bulk calculation:

```
bulk/            inp.xml (+ sym.xml, kpts.xml), converged cdn.hdf
gf/gfprep_inp    &GFPREP bulkdir='../bulk/' /
cd gf && gfprep
```

The bulk directory is only read. The third lattice vector must be
perpendicular to the first two. Non-collinear bulk calculations are not
supported.

## Input: `&GFPREP` in `gfprep_inp` (optional)

| parameter | default | meaning |
|---|---|---|
| `bulkdir` | `'./'` | directory of the bulk calculation |
| `prefix` | `'layer'` | layer files are `<prefix><i>_inp.xml` |
| `plane_tol` | 0.2 | atoms closer than this in z (a.u.) form one atomic plane |
| `planes_per_layer` | 1 | atomic planes per layer (auto mode) |
| `first_plane`, `last_plane` | 1, one period | range of atomic planes (auto mode). Planes continue periodically, so ranges can span several periods |
| `zplanes` | – | explicit embedding planes (cartesian z, ascending). Overrides auto mode |
| `dt_addon` | 1.5 | `dt = d + dt_addon` |
| `vshift` | 0 | constant shift of the potential; 999 means `-V_int(G=0)` |
| `spinflip` | F | per layer: exchange the spins |
| `kgrid` | bulk grid | 2D k-point grid shared by all layers |
| `gmax_scale` | 1.5 | Gmax of the layers relative to the bulk |

In auto mode, the embedding planes lie midway between the atomic planes.

## Output

* `<prefix><i>_inp.xml` describes layer i as a film:
  * `dTilda = dt`, `dVac = d`;
  * the atoms whose MT spheres lie within `|z| <= d/2`;
  * the bulk symmetry operations that leave the layer invariant;
  * the same k-points in every layer.
* `gf_inp` is a template with the `&LAYERS` geometry (`c1`, `c2`, `d`, `dt`, `c`). An existing `gf_inp` is kept.
* `gf_setup.hdf` holds, per layer:
  * `stepinfo`: the MT spheres of the layer and the outer spheres cut by the embedding planes;
  * `embdesc`: the descriptors of the two embedding planes.
* `gf_pot.hdf` holds the total potential (`pot`) and the Coulomb potential (`coulpot`). `gf_cdn.hdf` holds the density.
* `gfprep_planes.txt` lists the atomic planes and the chosen layers.
* `gfprep.out` holds the log, including a self-check: bulk and layer functions are compared at interstitial points with `|z| < c/2` and inside the MT spheres.

## Method

* **Muffin-tins:** the lattice-harmonic expansion of the bulk atom is projected onto the lattice harmonics of the layer by angular quadrature. This handles the lower site symmetry of the layer and a change of the local frame.
* **Interstitial:** the bulk function is kept unchanged for `|z| < d/2`. In the buffer `d/2 < |z| < dt/2` it is blended smoothly into its periodic image, then expanded with the period `dt` and projected onto the layer stars. The blending avoids Gibbs oscillations from a hard cut.
