# CaIrO3 RIXS reciprocal-coordinate mapping

This note records the coordinate mapping required before finite-Q RIXS physics
runs. It is specific to the validated CaIrO3 production cell and must be rechecked if
the cell representation changes.

## Coordinate convention in FLEUR

The source establishes the following representation:

- XML `relPos` values are stored in `atoms%taual`.
- `t_atoms%init` constructs Cartesian positions as
  `atoms%pos = matmul(cell%amat,atoms%taual)`.
- `t_cell%init` constructs `cell%bmat = 2*pi*inverse(cell%amat)` in FLEUR's
  stored matrix convention.
- `t_abc%calc_abc` evaluates its electronic site phase as
  `2*pi*dot_product(k+G,atoms%taual)`.

Therefore `atoms%taual` is the dimensionless fractional direct-cell coordinate
that must enter

```text
exp(+i 2*pi dot_product(Q_rlu,atoms%taual(:,iatom))).
```

`atoms%pos` is Cartesian and must not be used in this RLU dot product.

## Validated CaIrO3 production basis

The validated production input uses a diagonal Bravais matrix with lengths

```text
10.01276779  10.59822655  14.42913723 bohr
 5.298         5.608        7.636       Angstrom.
```

These match the experimental/manuscript Pbnm lattice parameters `(a,b,c)` in
the same order. There is no axis permutation or supercell transformation, so
the reciprocal-coordinate mapping is the identity:

```text
(H,K,L)_experimental Pbnm = (Q1,Q2,Q3)_FLEUR RLU.
```

In particular, use full momenta `(1,0,1)` and `(1/2,-1/2,1)` unchanged. The
four Ir representatives in this input are

```text
(0,0,0), (0,0,1/2), (1/2,1/2,0), (1/2,1/2,1/2).
```

The full Q is retained in the photon/site phase. Only the electronic k-point
lookup uses Q reduced modulo reciprocal lattice vectors.

For `Q=(1,0,1)`, the four site phases for the representatives above are
`(+1,-1,-1,+1)`. Thus the reduced electronic transfer is zero while the
basis-site interference is still nontrivial.

## Production-run gate

Before using another CaIrO3 input, compare its `cell%amat` columns with the
experimental Pbnm `(a,b,c)` axes.  If they are permuted, primitive, rotated, or
form a supercell, derive and record the reciprocal-basis transformation before
setting `momentumTransfer`.  Do not infer the mapping from the calculation name
or space-group label alone.

## k-mesh commensurability gate

The electronic mapping needs `k_v + Q` to land on another point of the same
full k list.  For a Gamma-centred `n1 x n2 x n3` mesh this requires the reduced
momentum `q = Q mod G` to satisfy

```text
q_i * n_i = integer   for i = 1,2,3.
```

If it does not, `rixs_build_kq_map` aborts with "found no k+q image for
k-point"; it never silently substitutes a nearby point.

For the two CaIrO3 targets:

```text
Q = (1,0,1)        -> q = (0,0,0)      commensurate with any mesh
Q = (1/2,-1/2,1)   -> q = (1/2,1/2,0)  needs even n1 and n2
```

The older 100-point `nx=5 ny=5 nz=4` list therefore supports `Q = (1,0,1)`
but **cannot** run `Q = (1/2,-1/2,1)`: 1/2 is not a multiple of 1/5, so every
k point fails to map and the run aborts immediately. The dense production
`16x16x12` full-BZ mesh has even `nx` and `ny`, so it supports the required
`q=(1/2,1/2,0)` permutation exactly.

Note also that `Q = (1,0,1)` gives `q = 0`, so `k_n = k_v` for every k point and
the electronic mapping is vertical.  That case exercises the coherent four-site
full-Q structure factor but not the k-remapping machinery; a commensurate
non-zero q is needed to exercise both at once.

## XML values for the two production transfers

Use the full experimental coordinates, not the reduced electronic q:

```xml
momentumTransfer="1.0 0.0 1.0"
```

and

```xml
momentumTransfer="0.5 -0.5 1.0"
```

The complete RIXS blocks and output conventions are documented in
[the RIXS user instructions](fleur_rixs_user_instructions.md).
