# j-resolved DOS with t2g-based j_eff decomposition

This document describes the extended j-resolved density-of-states (jDOS) feature
in FLEUR, which adds t2g-based effective angular momentum ($j_\text{eff}$) channels
for d-states and an optional per-atom reference-frame rotation.

---

## Physical background

For transition-metal oxides and similar systems with partially filled d-shells in an
octahedral crystal field, the d-states split into t2g (lower) and e_g manifolds.
Within the t2g triplet the effective orbital angular momentum is $l_\text{eff}=1$,
so spin-orbit coupling further splits the six t2g spinors into

| Channel | Degeneracy | $m_j$ values |
|---------|-----------|--------------|
| $j_\text{eff} = \tfrac{1}{2}$ | 2 | $-\tfrac{1}{2},\; +\tfrac{1}{2}$ |
| $j_\text{eff} = \tfrac{3}{2}$ | 4 | $-\tfrac{3}{2},\; -\tfrac{1}{2},\; +\tfrac{1}{2},\; +\tfrac{3}{2}$ |

The decomposition is carried out using Clebsch–Gordan coefficients for
$l_\text{eff}=1 \otimes s=\tfrac{1}{2}$, with the t2g real harmonics
($|xy\rangle,\,|yz\rangle,\,|xz\rangle$) mapped to $m_\text{eff}=+1,\,-1,\,0$,
respectively.

---

## Required input

### 1. Enable DOS and jDOS globally

In `inp.xml`, set both `dos` and `jDOS` flags:

```xml
<output dos="T">
  <bandDOS jDOS="T"
           minEnergy="-0.5" maxEnergy="0.5" sigma="0.02"
           numberPoints="1000" all_atoms="T" />
</output>
```

| Attribute | Required | Description |
|-----------|----------|-------------|
| `dos="T"` | yes | Activates DOS output |
| `jDOS="T"` | yes | Activates j-resolved decomposition |
| `minEnergy`, `maxEnergy` | yes | Energy window in Htr (relative to $E_\mathrm{F}$) |
| `sigma` | yes | Gaussian broadening in Htr |
| `numberPoints` | yes | Number of energy grid points |
| `all_atoms` | optional | Default `F`; set `T` to include every atom automatically |

### 2. Select atoms (when `all_atoms="F"`)

Mark individual atoms in the `atomGroup` position lists using the `banddos` attribute:

```xml
<atomGroup species="Fe-1">
  <relPos banddos="T">0.0 0.0 0.0</relPos>
</atomGroup>
```

### 3. Optional: rotate the local quantisation frame (per atom)

To project the spinor amplitudes onto a rotated reference frame before computing
the $j_\text{eff}$ decomposition, add Euler angles (in radians) to the same
position tag:

```xml
<atomGroup species="Fe-1">
  <relPos banddos="T" alpha="0.0" beta="0.7854" gamma="0.0">0.0 0.0 0.0</relPos>
</atomGroup>
```

| Attribute | Default | Description |
|-----------|---------|-------------|
| `alpha` | `0.0` | First Euler angle (rotation around global z before tilt) in radians |
| `beta`  | `0.0` | Polar tilt angle in radians |
| `gamma` | `0.0` | Azimuthal angle of the new local z-axis (see note below) in radians |

When all three angles are zero no rotation is applied and no extra cost is incurred.
The Wigner D-matrix rotation acts on all $l \leq 3$ channels (s through f), so
both the legacy j-channels and the new $j_\text{eff}$ channels are decomposed in
the rotated frame.

> **Tip:** The same `alpha`/`beta`/`gamma` attributes also control the reference
> frame for orbital-composition (`orbcomp`) DOS.  Setting them once affects both
> analyses simultaneously.

### 4. Optional: align the orbital frame to the local spin frame (per atom)

Instead of specifying Euler angles manually, you can request automatic alignment
to the noncollinear spin quantisation axis for each atom:

```xml
<atomGroup species="Fe-1">
  <relPos banddos="T" alignToSpin="T">0.0 0.0 0.0</relPos>
</atomGroup>
```

| Attribute | Default | Description |
|-----------|---------|-------------|
| `alignToSpin` | `F` | If `T`, uses `nococonv%alph/beta` (for the atom type) to define the orbital quantisation axis |

When `alignToSpin="T"` for an atom, jDOS uses the spin-frame mapping internally
(`beta = beta_spin`, `gamma = pi/2 - alpha_spin`, `alpha = 0`) and ignores any
manually provided `alpha`/`beta`/`gamma` values for that atom in the jDOS rotation.
Atoms with `alignToSpin="F"` keep the manual Euler-angle behavior.

### Euler-angle convention and relation to the noco spin frame

The `euler` routine in FLEUR uses a **ZXZ** convention:
$$R = R_z(\texttt{gamma}) \cdot R_x(-\texttt{beta}) \cdot R_z(\texttt{alpha})$$

The new local z-axis (quantisation axis) in global Cartesian coordinates is:
$$\hat{n}_\text{orb} = (\sin\gamma\sin\beta,\; \cos\gamma\sin\beta,\; \cos\beta)$$

The noco spin frame is defined by `nococonv%alph` (azimuthal, $\phi_s$) and
`nococonv%beta` (polar, $\theta_s$) per atom type, giving the local spin-up
direction:
$$\hat{n}_\text{spin} = (\sin\theta_s\cos\phi_s,\; \sin\theta_s\sin\phi_s,\; \cos\theta_s)$$

**For $\hat{n}_\text{orb} = \hat{n}_\text{spin}$ the mapping is:**
$$\texttt{beta} = \theta_s = \texttt{nococonv\%beta(itype)}$$
$$\texttt{gamma} = \tfrac{\pi}{2} - \phi_s = \tfrac{\pi}{2} - \texttt{nococonv\%alph(itype)}$$
$$\texttt{alpha} = \text{arbitrary (rotates within the local plane, does not move the axis)}$$

Setting `banddos%alpha = nococonv%alph` and `banddos%beta = nococonv%beta` as a
naive mapping does **not** align the frames; the azimuthal angle enters
`banddos%gamma` with a $\pi/2$ offset. Also note that `nococonv%alph/beta` are
indexed **per atom type** while `banddos%alpha/beta/gamma` are **per atom** — use
the same type-level angles for all atoms of the same type.

This convention mismatch is a pre-existing property of the `orbcomp` feature
(which uses the same `t_abc%rotate` function) and is not specific to jDOS.

---

## Output

FLEUR writes the DOS to `banddos.hdf` (post-processed from the `jDOS` weights).
Each selected atom contributes **15 weight channels** in the following order:

### Legacy j-resolved channels (7 per atom)

| Label | Channel |
|-------|---------|
| `jDOS:<N>s`     | s-state (not split by SOC) |
| `jDOS:<N>p1-2`  | $p_{1/2}$ ($j = \tfrac{1}{2}$) |
| `jDOS:<N>p3-2`  | $p_{3/2}$ ($j = \tfrac{3}{2}$) |
| `jDOS:<N>d3-2`  | $d_{3/2}$ ($j = \tfrac{3}{2}$) |
| `jDOS:<N>d5-2`  | $d_{5/2}$ ($j = \tfrac{5}{2}$) |
| `jDOS:<N>f5-2`  | $f_{5/2}$ ($j = \tfrac{5}{2}$) |
| `jDOS:<N>f7-2`  | $f_{7/2}$ ($j = \tfrac{7}{2}$) |

`<N>` is the global atom index.

### t2g $j_\text{eff}$ total channels (2 per atom)

| Label | Channel |
|-------|---------|
| `jDOS:<N>dj1/2` | Total $j_\text{eff} = \tfrac{1}{2}$ weight within t2g |
| `jDOS:<N>dj3/2` | Total $j_\text{eff} = \tfrac{3}{2}$ weight within t2g |

### t2g $j_\text{eff}$ $m_j$-resolved channels (6 per atom)

| Label | $j_\text{eff}$ | $m_j$ |
|-------|----------------|-------|
| `jDOS:<N>d1/2m-1` | $\tfrac{1}{2}$ | $-\tfrac{1}{2}$ |
| `jDOS:<N>d1/2m+1` | $\tfrac{1}{2}$ | $+\tfrac{1}{2}$ |
| `jDOS:<N>d3/2m-3` | $\tfrac{3}{2}$ | $-\tfrac{3}{2}$ |
| `jDOS:<N>d3/2m-1` | $\tfrac{3}{2}$ | $-\tfrac{1}{2}$ |
| `jDOS:<N>d3/2m+1` | $\tfrac{3}{2}$ | $+\tfrac{1}{2}$ |
| `jDOS:<N>d3/2m+3` | $\tfrac{3}{2}$ | $+\tfrac{3}{2}$ |

> **Note on normalisation:** All weights are expressed as a percentage of the
> muffin-tin sphere contribution (`qmtp`), consistent with the legacy j-channels.
> A value of 100 means the full MT charge is in that channel.

### Occupation numbers

The integrated occupations (sum over k-points and bands weighted by occupation)
are printed to stdout at the end of the calculation for all channels:
the 7 legacy `occ(l,jj,atom)` arrays, the 2 total $j_\text{eff}$ occupations
`occ_jeff_d(1:2, atom)`, and the 6 $m_j$-resolved occupations
`occ_jeff_d_mj(1:6, atom)`.

---

## Minimal working example

```xml
<output dos="T">
  <bandDOS jDOS="T" minEnergy="-0.8" maxEnergy="0.8" sigma="0.02"
           numberPoints="2000" all_atoms="F" />
</output>

<!-- In the atomGroups section: -->
<atomGroup species="Ir-1">
  <!-- Standard quantisation axis (z along [001]) -->
  <relPos banddos="T">0.0 0.0 0.0</relPos>
</atomGroup>

<atomGroup species="Ir-2">
  <!-- Rotate to quantise along [111]: beta = arccos(1/sqrt(3)) ≈ 0.9553 rad
       [111] has phi=pi/4, so gamma = pi/2 - pi/4 = pi/4 ≈ 0.7854 rad -->
  <relPos banddos="T" alpha="0.0" beta="0.9553" gamma="0.7854">0.5 0.5 0.5</relPos>
</atomGroup>
```

After running FLEUR, open `banddos.hdf` with the `masci-tools` Python library or
the `fleur_inpgen` post-processor to plot any of the 15 channels per atom.
