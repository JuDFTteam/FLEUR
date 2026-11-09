# gfprep

`gfprep` turns a converged FLEUR bulk supercell, stacked along z, or a film
into the layered setup used by `gfleur`. It is built with the gfleur targets
(HDF5 is required).

## Usage

Run it in an empty directory, next to the bulk calculation:

```
bulk/            inp.xml (+ sym.xml, kpts.xml), converged cdn.hdf
gf/gfprep_inp    &GFPREP bulkdir='../bulk/' /
cd gf && gfprep
```

The bulk directory is only read; it may also hold a film calculation. The
third lattice vector must be perpendicular to the first two.

Non-collinear bulks and films (`l_noco`, also with `l_mperp` and `l_mtNocoPot`)
are supported: the layer inp.xml gets the noco setup of the bulk and the
angles `alpha`, `beta` of the bulk type of each layer atom group. Not supported:
spin spirals, `l_relaxSQA`, constrained moments, `l_mtNocoPot` only for some
types, and `spinflip`.

## Input: `&GFPREP` in `gfprep_inp` (optional)

| parameter | default | meaning |
|---|---|---|
| `bulkdir` | `'./'` | directory of the bulk calculation |
| `prefix` | `'layer'` | layer files are `<prefix><i>_inp.xml` |
| `plane_tol` | 0.2 | atoms closer than this in z (a.u.) form one atomic plane |
| `planes_per_layer` | 1 | atomic planes per layer (auto mode) |
| `first_plane`, `last_plane` | 1, one period (film: the last plane) | range of atomic planes (auto mode). In bulk the planes continue periodically, so ranges can span several periods |
| `zplanes` | – | explicit embedding planes (cartesian z, ascending). Overrides auto mode |
| `dt_addon` | 1.5 | `dt = d + dt_addon` |
| `vshift` | 0 | constant shift of the potential; 999 means `-V_int(G=0)` |
| `spinflip` | F | per layer: exchange the spins |
| `kgrid` | bulk grid | 2D k-point grid shared by all layers |
| `gmax_scale` | 1.5 | Gmax of the layers relative to the bulk |
| `zvac` | 0 | film: the vacuum embedding planes lie at `±(z1+zvac)`, `z1` = dVac/2 of the film |

In auto mode, the embedding planes lie midway between the atomic planes.

For a film, `first_plane` and `last_plane` count the atomic planes of the
film. A layer range that includes the lowest (highest) atomic plane ends at
a vacuum embedding plane at `-(z1+zvac)` (`z1+zvac`). This plane must not cut
a MT sphere.

## Output

* `<prefix><i>_inp.xml` describes layer i as a film:
  * `dTilda = dt`, `dVac = d`;
  * the atoms whose MT spheres lie within `|z| <= d/2`;
  * the bulk symmetry operations that leave the layer invariant;
  * the same k-points in every layer.
* `gf_inp` is a template with the `&LAYERS` geometry (`c1`, `c2`, `d`, `dt`, `c`). For a film with vacuum embedding planes it contains `&EMB addemb=T, surface=..., vacleft=... /`. An existing `gf_inp` is kept.
* `gf_setup.hdf` holds, per layer:
  * `stepinfo`: the MT spheres of the layer and the outer spheres cut by the embedding planes;
  * `embdesc`: the descriptors of the two embedding planes.
* `gf_pot.hdf` holds the total potential (`pot`) and the Coulomb potential (`coulpot`). `gf_cdn.hdf` holds the density.
  For noco, `pot` also contains `pwod` (the interstitial off-diagonal part V21 or rho21, global frame) and, with `l_mperp`, `vrod` (the MT components 3:4, local frame of each atom).
* For a film, `gf_pot.hdf` also holds `vacemb_side1` and `vacemb_side2`: the planar-averaged potential `vz(z,spin)` (noco: `V11, V22, Re V21, Im V21`) beyond the vacuum embedding plane below layer 1 and above the last layer, from the plane to the end of the FLEUR vacuum grid with spacing `delz` (attribute).
* `gfprep_planes.txt` lists the atomic planes and the chosen layers.
* `gfprep.out` holds the log, including a self-check: bulk and layer functions are compared at interstitial points with `|z| < c/2` and inside the MT spheres. For noco there are additional lines `potential21` and `density21` for the off-diagonal parts.

## Method

* **Muffin-tins:** the lattice-harmonic expansion of the bulk atom is projected onto the lattice harmonics of the layer by angular quadrature. This handles the lower site symmetry of the layer and a change of the local frame.
* **Interstitial:** the bulk function is kept unchanged for `|z| < d/2`. In the buffer `d/2 < |z| < dt/2` it is blended smoothly into its periodic image, then expanded with the period `dt` and projected onto the layer stars. The blending avoids Gibbs oscillations from a hard cut.
* **Film:** the film function is the plane-wave expansion for `|z| < z1` and the vacuum expansion (planar average and warping) beyond. It is sampled in the layer window and projected by quadrature with the same blending. The plane-wave and vacuum expansions of FLEUR do not match exactly at `|z| = z1`; the self-check skips points within 0.1 a.u. of `z1`.
