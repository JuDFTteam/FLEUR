# Band-resolved local density matrix ("DensityMatrix" DOS)

With this DOS mode, FLEUR stores the full local density matrix of every eigenstate,
k-point and selected atom in `banddos.hdf`. Projections that otherwise need a dedicated
FLEUR type can be computed afterwards in a few lines of Python, from the stored matrices.
Examples are l-, m-, j- or j_eff-resolved weights, t2g/eg splitting, cubic or
crystal-field harmonics, arbitrary local frames, and orbital or spin moments per state.

---

## Definition

For eigenstate ν at k-point **k**, the muffin-tin part of the wave function of atom α is
expanded as

$$\psi_{\nu\mathbf{k}}^{s}(\mathbf{r}) = \sum_{lm}\sum_i c^{s}_{\nu,lm,i}\, u^{ls}_i(r)\, Y_{lm}(\hat{\mathbf{r}}),$$

where *i* runs over all radial functions of the channel: u, u̇ and all local orbitals.
For a block (l, l') and a spin pair (s, s'), the stored matrix is

$$n^{\nu\mathbf{k}}_{m m'}(l l', s s') = \sum_{ij} c^{s}_{\nu,l'm',i}\,\big(c^{s'}_{\nu,lm,j}\big)^{*}\,\langle u^{l's}_i | u^{ls'}_j\rangle
 = \langle l' m' s|\hat\rho_{\nu\mathbf{k}}|l m s'\rangle .$$

This is the convention of the LDA+U density matrix n_mmp in FLEUR, so that
$\sum_{\nu\mathbf{k}} w_{\nu\mathbf{k}}\, n^{\nu\mathbf{k}}(ll, ss)$ is the usual occupation matrix.
The basis consists of FLEUR's complex spherical harmonics $Y_{lm}$.

- **Radial overlaps:** diagonal blocks (l = l') use the radial overlaps of the density
  calculation. Blocks with l ≠ l' use the overlaps $\langle u^{l'}_i|u^{l}_j\rangle$ of
  the radial functions of both channels.
- **Blocks between different l:** only l < l' is stored. The partner block is
  $n(l'l, s's) = n(ll', ss')^{\dagger}$.
- **Spin blocks:** `11` for non-spin-polarized runs, `11 22` for collinear runs,
  and `11 22 21 12` for non-collinear runs (block `12` is the per-band conjugate
  transpose of `21`).
- **Charge and moments:** $n_{\rm charge}=\mathrm{Tr}(n_{11}+n_{22})$,
  $m_z=\mathrm{Tr}(n_{11}-n_{22})$, $m_x=2\,\mathrm{Re\,Tr}\,n_{21}$,
  $m_y=2\,\mathrm{Im\,Tr}\,n_{21}$. This is the convention of `denmat_to_mag` in FLEUR.

---

## Input

The mode is switched on by a `densityMatrix` element inside `bandDOS`. It needs
`dos="T"` or `band="T"` and `fleurInputVersion="0.39"`.

```xml
<output dos="T">
  <bandDOS minEnergy="-0.5" maxEnergy="0.5" sigma="0.015">
    <densityMatrix l="spd" lCross="pd" symmetrize="auto"/>
  </bandDOS>
</output>
```

| Attribute | Values | Meaning |
|---|---|---|
| `l` | subset of `spdf` (default `spd`) | l channels to store |
| `lCross` | `none` (default), `all`, or pairs such as `"sd pd"` | additional blocks between different l |
| `symmetrize` | `auto` (default), `T`, `F` | see below |

**Atom selection** follows the other DOS modes: the `banddos` attribute of the atom
positions (`all_atoms` in `bandDOS`).

**Local orbital frame:** the per-atom Euler angles `alpha`, `beta`, `gamma` or
`alignToSpin` rotate the matrices into that frame, exactly as for orbcomp and jDOS.
The angles used are written as attribute `eulerAngles`.

**Spin frame:** in non-collinear calculations, `globalspin="T"` in `bandDOS` rotates the
spin blocks from the local spin frame of each atom type into the global frame, as it
does for the Local DOS.

**Symmetrization:**
- `T`: each matrix is averaged over the site symmetry of the atom and over the equivalent
  atoms, which gives one element per atom type. The k-point weights of the irreducible
  wedge then represent the full Brillouin zone.
- `F`: every selected atom is stored separately. It is stored in the frame of its
  muffin-tin functions, which for equivalent atoms is rotated by the symmetry operation
  that maps the atom onto the representative one.
- `auto`: symmetrizes unless the k-points form a path (band structures).

In all modes, the matrices are averaged over degenerate states.

Memory and file size scale as 49 · (spin blocks) · neig · nkpt · (number of blocks) ·
(number of elements) complex numbers. The size is printed to `out`; keep an eye on it for
large k-point sets.

---

## Output in banddos.hdf

```
/DensityMatrix/DOS, /DensityMatrix/EV, /DensityMatrix/BS   scalar weights, see below
/DensityMatrix/matrix                   attrs: symmetrized, globalSpinFrame, nSpinBlocks
/DensityMatrix/matrix/eigenvalues       (spin, k, band) in Hartree, not shifted
/DensityMatrix/matrix/elem_<n>          attrs: atom, atomType, eulerAngles
/DensityMatrix/matrix/elem_<n>/pair_<l><l'>   e.g. pair_dd, pair_pd
```

**Dataset shape:** in h5py (C order), a `pair_<l><l'>` dataset has the shape
`(nkpt, neig, nSpinBlocks, 2l'+1, 2l+1, 2)`. Read in this order, the array is the ordinary
density matrix, with rows l' and columns l:
`pair[k, band, b, i, j, :] = (Re, Im) of <l', m'=i-l', s|rho|l, m=j-l, s'>`,
where (s, s') is spin block b.

**Scalar weights:** these are the traces of the diagonal blocks, named `DM:<atom><l>` per
spin for collinear runs. Non-collinear runs add `DM:<atom><l>mx`, `...my` and `...mz`.
They behave like any other DOS or band weight. For example, the masci-tools plotting
routines can use them.

**Built-in check:** after the calculation, FLEUR compares these traces with the
l-resolved weights `MT:<type><l>` of the Local DOS and writes the largest difference
to `out`:

```
 DensityMatrix DOS: max. deviation of the traces from the l-resolved DOS weights:   2.776E-16
```

Symmetrized data is compared per atom type. Unsymmetrized data is averaged over the
equivalent atoms and, in non-collinear runs, also compared for the spin off-diagonal
part.

---

## Python example

```python
import h5py, numpy as np

f = h5py.File("banddos.hdf")
m = f["DensityMatrix/matrix"]
a = np.asarray(m["elem_1/pair_dd"])
rho = a[..., 0] + 1j*a[..., 1]        # (k, band, spin block, m, m') = <m s|rho|m' s'>

# d charge of each state (non-collinear: spin blocks 11 and 22)
nd = np.trace(rho[:, :, :2], axis1=-2, axis2=-1).real.sum(axis=-1)

# real d orbitals (m = -2..2), the same combinations as orbcomp indices 5..9
h = np.sqrt(0.5)
T = np.array([[h, 0, 0, 0, -h],      # xy
              [0, h, 0, h, 0],       # yz
              [0, h, 0, -h, 0],      # zx
              [h, 0, 0, 0, h],       # x2-y2
              [0, 0, 1, 0, 0]])      # z2
w = np.einsum("am,kbmn,an->kba", T, rho[:, :, 0], np.conj(T)).real   # (k, band, orbital)
t2g = w[..., [0, 1, 2]].sum(-1)
eg = w[..., [3, 4]].sum(-1)
```

For the j-resolved weights of FLEUR's jDOS, combine the four spin blocks with the
Clebsch–Gordan coefficients of $l\otimes\tfrac12$, with a minus sign on the spin-down
component as in `calc_jDOS`. The regression tests `basic/CuDM` and `noco/FeBccNocoDM`
cover the definitions above. The jDOS reconstruction was checked against `jDOS:*`
to 1e-15.

---

## Limitations

- l ≤ 3.
- Symmetrization in non-collinear calculations uses the spin phases of the LDA+U code.
  Those are only meaningful if the magnetic structure is compatible with the symmetry
  operations. For non-collinear runs with symmetry, prefer `symmetrize="F"`.
- Non-collinear calculations without `l_mperp` compute both spins in a single pass when
  the mode is active, as jDOS does. The results are unchanged, but the `mtCharges` of both
  spins then appear in a single `out.xml` element.
