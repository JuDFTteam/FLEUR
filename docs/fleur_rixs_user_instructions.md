# FLEUR independent-particle RIXS: user instructions

This note describes the XML-controlled independent-particle RIXS capability in
FLEUR. Two paths are available. Omitting `momentumTransfer` selects the legacy
same-k implementation. Supplying `momentumTransfer` opts into the finite-Q
first-variation spinor implementation. The two paths intentionally have
different absorber-site coherence and spin-support rules.

## Current scope

The legacy same-k path supports:

- scalar `jspins=1` calculations;
- first-variation noco spinor calculations with `jspins=2`;
- vertical electronic transitions at the same k point;
- an incoherent sum over absorber-site intensities.

The opt-in finite-Q path supports:

- first-variation noco spinors with `jspins=2`, with or without
  first-variation SOC;
- electronic `k -> k+q` mapping on a commensurate full k mesh;
- coherent interference of all matching absorber sites using the full photon
  momentum transfer;
- separate complex site-partial and coherently summed pair diagnostics.

Both paths support:

- edges `K`, `L2`, and `L3`;
- linear polarizations `x`, `y`, and `z`;
- absorber selection by `absorberZ`;
- optional rank-local contribution tables;
- optional valence/intermediate band windows;
- pure k-point and shared-k-point-subgroup MPI layouts;
- MPI reduction of the broadened spectrum;
- normalized Gaussian broadening of the energy-loss spectrum.

Current guards and limitations:

- no collinear spin-polarized scalar RIXS;
- no second-variation SOC RIXS; SOC is supported through the validated
  first-variation noco spinor path;
- full-k/no-star k meshes only, i.e. `nkpt == nkptf`;
- no explicit core hole;
- no multiplets;
- no coherent symmetry-star reconstruction yet;
- no combined `L23` edge; run `L2` and `L3` separately.

The finite-Q implementation is local electric-dipole E1-E1 RIXS. It does not
include intra-muffin-tin variation of the photon plane wave, E2 terms, an
explicit core hole, electron-hole attraction, or multiplet physics.

## Formula

### Legacy same-k path

For scalar bands, the legacy path uses the spin-degenerate S1 trace. The
final electron and valence-hole spin labels are orthogonal final-state labels
and are traced incoherently:

```text
I(Omega) = sum_{k,v,n,a} w_k f_vk (1 - f_nk)
           sum_{sigma_v,sigma_n}
           | sum_mj
               M_emit^a(v,k,mj,sigma_v; eps_out)
               M_abs^a(n,k,mj,sigma_n; eps_in)
             / (omegaIn - (eps_nk - eps_core^a) + i gammaCore)
           |^2
           g_etaLoss[Omega - (eps_nk - eps_vk)]
```

Here `v` is an occupied valence band, `n` is an intermediate unoccupied band,
`a` is an absorber atom, and `mj` labels the selected core sublevel. Occupation
filtering remains active: valence states require `f_vk > tol`, and intermediate
states require `1 - f_nk > tol`.

For legacy first-variation spinor bands, FLEUR instead stores the coherent complex
core-`mj` amplitude for each valence/intermediate band pair and local absorber
atom. The final spectrum remains an incoherent sum over distinct band-pair and
absorber-site final states.

### Finite-Q path

The XML value `momentumTransfer` is the full experimental transfer
`Q_rlu = (H,K,L)` in the reciprocal basis dual to the FLEUR direct-cell basis.
FLEUR forms the reduced transfer

```text
q = Q_rlu modulo reciprocal-lattice vectors
```

only for the electronic mapping

```text
k_n = (k_v + q) modulo reciprocal-lattice vectors.
```

It does not replace the full Q by q in the absorber-site phase. For a valence
band `v`, intermediate band `n`, valence k point `k_v`, and absorber site `s`,
the implemented first-variation spinor amplitude is

```text
A_vn,k(Q) = sum_s exp(+i 2*pi Q_rlu.tau_s)
                    sum_mj M_out(v,k_v,s,mj) M_in(n,k_n,s,mj)
                    / [omegaIn - (epsilon_n,k_n - epsilon_core,s)
                       + i gammaCore].
```

Here `tau_s` is the dimensionless fractional direct-cell coordinate stored in
`atoms%taual`. The matrix elements use

```text
M_in  = <psi_n,k_n | epsilon_in.r_s | c_s,mj>
M_out = <c_s,mj | (epsilon_out.r_s)^dagger | psi_v,k_v>.
```

The emission routine applies the required conjugation internally. Callers
must pass `epsilon_out` directly and must not conjugate it themselves.

The loss spectrum contains

```text
sum_{k,v,n} w_k f_v,k [1-f_n,k+q] |A_vn,k(Q)|^2
              g_etaLoss[loss - (epsilon_n,k+q - epsilon_v,k)].
```

Thus absorber sites interfere coherently within each band pair, while distinct
band-pair final states are summed incoherently. Orthogonal outgoing
polarizations are also summed only after squaring. For example, the
unanalysed signal for incoming `z` and outgoing `x,y` is `I_zx + I_zy`, never
the square of a sum of the `zx` and `zy` complex amplitudes.

## Minimal XML input

Add the RIXS block under `<output>`:

```xml
<output>
   <rixs l_rixs="T"
         absorberZ="26"
         edge="L3"
         omegaIn="30.0"
         gammaCore="0.20"
         lossMin="0.0"
         lossMax="10.0"
         nLoss="501"
         etaLoss="0.05"
         incomingPolarizations="x y z"
         outgoingPolarizations="x y z"
         outputPrefix="rixs"/>
</output>
```

All energies and broadenings in this block are in Hartree.

If no `<rixs>` section is present, RIXS is off.

If `<rixs>` is present but `momentumTransfer` is absent, FLEUR uses the legacy
same-k path. `momentumTransfer` is therefore an explicit opt-in switch for the
finite-Q path; a numerically zero value such as `momentumTransfer="0 0 0"`
still selects finite-Q behavior and coherent absorber-site amplitudes.

## Finite-Q Ir L3 examples

The following uses the actual XML attribute accepted by the schema for the
full experimental transfer `Q=(1,0,1)`:

```xml
<rixs l_rixs="T"
      absorberZ="77"
      edge="L3"
      omegaIn="412.1804"
      gammaCore="0.20"
      lossMin="0.0"
      lossMax="0.03"
      nLoss="601"
      etaLoss="0.0004"
      incomingPolarizations="z"
      outgoingPolarizations="x y"
      momentumTransfer="1.0 0.0 1.0"
      valenceBandMin="169"
      valenceBandMax="172"
      intermediateBandMin="173"
      intermediateBandMax="180"
      outputPrefix="cairo3_101"
      writeContributions="F"
      writeStateCharacter="F"/>
```

Electronically, `Q=(1,0,1)` reduces to `q=(0,0,0)`, so `k_n=k_v`. It is not
equivalent to the legacy path: the finite-Q amplitude retains the generally
nontrivial four-Ir factors `exp(+i 2*pi Q.tau_s)` and sums the site amplitudes
coherently.

For the CaIrO3 transfer `Q=(1/2,-1/2,1)`, change the full-Q attribute and use a
commensurate mesh:

```xml
<rixs l_rixs="T"
      absorberZ="77"
      edge="L3"
      omegaIn="412.1804"
      gammaCore="0.20"
      lossMin="0.0"
      lossMax="0.03"
      nLoss="601"
      etaLoss="0.0004"
      incomingPolarizations="z"
      outgoingPolarizations="x y"
      momentumTransfer="0.5 -0.5 1.0"
      valenceBandMin="169"
      valenceBandMax="172"
      intermediateBandMin="173"
      intermediateBandMax="180"
      outputPrefix="cairo3_half_half"
      writeContributions="F"
      writeStateCharacter="F"/>
```

For the documented CaIrO3 reciprocal basis this reduces electronically to
`q=(1/2,1/2,0)`. A Gamma-centred `5x5x4` mesh cannot represent translation by
one half in its first two directions: adding `1/2` places every point between
the available multiples of `1/5`. A `16x16x12` full-BZ mesh has even first and
second dimensions and maps exactly under this transfer. FLEUR requires an
exact one-to-one k-point permutation and aborts instead of substituting a
nearby point. See [CaIrO3 reciprocal-coordinate mapping](caIrO3_rixs_reciprocal_mapping.md)
for the cell-specific basis check.

## Recommended first physical example: MgO O K edge

MgO is a useful first physical example for the scalar RIXS branch: it has a
small cell, and the O K edge avoids the L-edge core spin-angular complications.
The validated MgO benchmark specifically uses scalar PBE without SOC or noco,
an explicit full-k/no-star list, 60 bands, and occupied bands `1 ... 8`.

The following block is a convenient starting point:

```xml
<rixs l_rixs="T"
      absorberZ="8"
      edge="K"
      omegaIn="18.94"
      gammaCore="0.20"
      lossMin="0.0"
      lossMax="2.0"
      nLoss="401"
      etaLoss="0.03"
      incomingPolarizations="x y z"
      outgoingPolarizations="x y z"
      outputPrefix="rixs_mgo_ok"
      writeContributions="T"/>
```

For the validated scalar MgO calculation, `omegaIn="18.94"` Ha was selected
from the calculated O K-edge XAS maximum. This incident energy and the band
indices below are benchmark-specific, not general O K-edge defaults. For a new
structure, functional, or k mesh, inspect its calculated XAS spectrum or
transition table and select the desired resonance again. The loss window and
number of loss points are analysis choices and should be converged for the
relevant valence-to-conduction excitation energies.

The validated calculation used occupied bands `1 ... 8` and unoccupied bands
`9 ... 60`. The corresponding optional RIXS band windows are:

```xml
<rixs l_rixs="T"
      absorberZ="8"
      edge="K"
      omegaIn="18.94"
      gammaCore="0.20"
      lossMin="0.0"
      lossMax="4.0"
      nLoss="801"
      etaLoss="0.03"
      incomingPolarizations="x"
      outgoingPolarizations="x"
      outputPrefix="rixs_mgo_ok_xx"
      writeContributions="T"
      valenceBandMin="1"
      valenceBandMax="8"
      intermediateBandMin="9"
      intermediateBandMax="60"/>
```

These windows remain subject to the normal occupation/vacancy filters. They
apply to the validated 60-band MgO fixture and must be adapted if the band
ordering, electron count, or requested number of bands changes. The strongest
O K-edge intermediate transitions in this benchmark are mainly in bands
`9 ... 13`, but bands through 60 were retained for the validated spectrum.

## Attributes

| Attribute | Required? | Meaning | Example |
|---|---:|---|---|
| `l_rixs` | no | Enables or disables RIXS. If the `<rixs>` element is present, the intended use is `T`. | `T` |
| `absorberZ` | yes when RIXS is enabled | Atomic number of absorber atoms. All matching atom types are summed. | `26` |
| `edge` | yes | Supported: `K`, `L2`, `L3`. | `L3` |
| `omegaIn` | yes | Incident photon energy in Hartree. | `30.0` |
| `gammaCore` | yes | Core-hole lifetime broadening in the intermediate denominator, in Hartree. | `0.20` |
| `lossMin`, `lossMax` | yes | Energy-loss output window in Hartree. | `0.0`, `10.0` |
| `nLoss` | yes | Number of points in the loss-energy grid. | `501` |
| `etaLoss` | yes | Standard deviation (sigma) of the normalized Gaussian loss broadening, in Hartree. | `0.05` |
| `incomingPolarizations` | no | Space-separated list of incoming linear polarizations. Default: `x`. | `x y z` |
| `outgoingPolarizations` | no | Space-separated list of outgoing linear polarizations. Default: `x`. | `x y z` |
| `outputPrefix` | no | Prefix for RIXS output files. Default: `rixs`. | `rixs` |
| `momentumTransfer` | no | Full experimental `Q_rlu=(H,K,L)` in the reciprocal basis dual to the FLEUR cell. Its presence selects finite-Q RIXS. | `0.5 -0.5 1.0` |
| `writeContributions` | no | Write rank-local valence/intermediate-band contribution tables. Default: `F`. | `T` |
| `writeStateCharacter` | no | Write rank-local HDF5 band-character sidecars for legacy first-variation spinor RIXS. Currently unsupported with `momentumTransfer`. Default: `F`. | `T` |
| `stateLigandZ` | required with `writeStateCharacter="T"` | Atomic number of candidate ligands used to construct the local structural frame. | `8` |
| `valenceBandMin`, `valenceBandMax` | no | Optional 1-based valence-band loop window. | `1`, `20` |
| `intermediateBandMin`, `intermediateBandMax` | no | Optional 1-based intermediate-band loop window. | `21`, `80` |

Supported edge aliases are the same as for XML XAS:

- `K`, `k`, `1s1/2`, `1S1/2`
- `L2`, `l2`, `2p1/2`, `2P1/2`
- `L3`, `l3`, `2p3/2`, `2P3/2`

`L23` intentionally aborts for now.

The normalized loss broadening is

```text
g_eta(x) = exp[-0.5*(x/etaLoss)^2] / [etaLoss*sqrt(2*pi)].
```

Consequently `etaLoss` is sigma, not the full width at half maximum. The
corresponding `FWHM` is approximately `2.35482*etaLoss`.

## Optional band windows

Band windows restrict the explicit valence (`v`) and intermediate (`n`) band
loops before occupation/vacancy filtering:

```xml
<rixs l_rixs="T"
      absorberZ="26"
      edge="L3"
      omegaIn="30.0"
      gammaCore="0.20"
      lossMin="0.0"
      lossMax="10.0"
      nLoss="501"
      etaLoss="0.05"
      incomingPolarizations="x"
      outgoingPolarizations="x"
      outputPrefix="rixs"
      writeContributions="T"
      valenceBandMin="1"
      valenceBandMax="20"
      intermediateBandMin="21"
      intermediateBandMax="80"/>
```

Partial windows are allowed:

- `valenceBandMax="20"` means valence bands `1 ... 20`;
- `valenceBandMin="5"` means valence bands `5 ... all`;
- `intermediateBandMin="21"` means intermediate bands `21 ... all`.

Bounds are 1-based FLEUR band indices. Non-positive bounds and `min > max`
abort. Upper bounds are clamped to the number of bands available at each k
point. If a lower bound is above the available band count at a k point, that
k point is skipped safely.

## Output files

Legacy spectrum files are written on rank 0 as

```text
<outputPrefix>_<edge>_<incoming>_<outgoing>_omega<omegaIn>_rixs.dat
```

For example:

```text
rixs_L3_x_x_omega30p000000_rixs.dat
rixs_L3_x_y_omega30p000000_rixs.dat
```

Each spectrum file contains:

```text
loss_energy_Ha loss_energy_eV intensity
```

Finite-Q spectrum names retain the full Q label:

```text
<outputPrefix>_<edge>_<incoming>_<outgoing>_Q_<H>_<K>_<L>_omega<omegaIn>_finiteq_rixs.dat
```

The three full-Q coordinates are formatted to six decimal places with `p` for
the decimal point and `m` for a minus sign. Separate files are emitted for
every requested incoming/outgoing polarization pair. FLEUR does not create a
combined unanalysed-polarization file; form sums such as `I_zx + I_zy` from
the separate intensity columns.

If `writeContributions="T"`, one rank-local contribution table is written per
requested polarization pair and MPI rank. For the legacy path the name is

```text
<outputPrefix>_<edge>_<incoming>_<outgoing>_omega<omegaIn>_contrib_rank0000.dat
```

The contribution table includes the k point, valence band, intermediate band,
absorber atom/type, occupations, denominator, loss energy, `amplitude_abs2`, and
`weighted_strength`.

For scalar S1 production, `amplitude_abs2` is already the spin-traced quantity:

```text
sum_{sigma_v,sigma_n} |A_{sigma_v sigma_n}|^2
```

`amplitude_real` and `amplitude_imag` are zero placeholders because no single
coherent complex amplitude represents this incoherent spin trace. The row
weight is:

```text
weighted_strength = k_weight * f_v * (1 - f_n) * amplitude_abs2
```

For first-variation spinor RIXS, `amplitude_real` and `amplitude_imag` store
the coherent core-`mj` amplitude for the band pair, and `amplitude_abs2` is its
squared modulus.

When contribution output is enabled, FLEUR also prints a
contribution-to-spectrum consistency check. It should report `PASS`.

For finite-Q, two rank-local tables are written for each polarization pair:

```text
<outputPrefix>_<edge>_<incoming>_<outgoing>_Q_<H>_<K>_<L>_omega<omegaIn>_finiteq_pair_rank0000.dat
<outputPrefix>_<edge>_<incoming>_<outgoing>_Q_<H>_<K>_<L>_omega<omegaIn>_finiteq_site_rank0000.dat
```

The pair table contains `k_v`, `k_n`, the integer reciprocal-lattice shift,
band energies and occupations, the coherently summed complex pair amplitude,
its squared modulus, and its weighted strength. Its shift obeys
`k_v + Q_full = k_n + reciprocal_shift` for the k representatives printed in
the row. The site table contains `tau_fractional`, the full-Q phase, the
intermediate denominator, the unphased local amplitude, and the phased complex
partial amplitude for every absorber site. Summing `phased_partial` as complex
numbers over sites reconstructs the pair amplitude. Summing sitewise squared
moduli does not reconstruct the finite-Q spectrum.

Both tables can be very large. They are intended for focused interference and
accounting audits, not routine dense production runs. The finite-Q pair table
also participates in the contribution-to-spectrum consistency check.

## State-character sidecars

With a first-variation spinor legacy calculation, `writeStateCharacter="T"`
and a positive `stateLigandZ` write one rank-local HDF5 sidecar per owner rank:

```text
<outputPrefix>_<edge>_state_character_rank0000.hdf
```

The sidecars annotate band states in the configured valence/intermediate
windows with local absorber-site d-spin, orbital, and `j_eff` character. These
are state descriptors, not an incoherent decomposition of RIXS intensity.
This option requires an HDF5-enabled build. It is deliberately guarded off in
the finite-Q driver because a sidecar convention for distinct valence and
intermediate k points has not yet been defined.

## MPI execution model

FLEUR has two MPI levels. `fmpi%irank/isize` identify ranks in the global
communicator, while `fmpi%n_rank/n_size` identify ranks in the subgroup that
collaborates on one k point. All ranks in a k-point subgroup share the same
`fmpi%k_list`; iterating over that list does not by itself assign unique RIXS
work to subgroup ranks.

RIXS currently evaluates each k point serially rather than distributing its
band-pair transitions within a subgroup. Consequently, only the subgroup root
(`fmpi%n_rank == 0`) evaluates RIXS matrix elements, contribution rows, and the
local spectrum for that subgroup. Other subgroup ranks retain zero local RIXS
arrays and still participate in the collective global reductions. Global rank
zero (`fmpi%irank == 0`) remains responsible for the final spectrum files and
summary output.

This gives correct results for both:

- pure k-point parallelism, for example `-pe_per_kpt 1`;
- layouts where multiple MPI ranks share each k point.

The latter layout currently provides correctness but no intra-k-point RIXS
speedup, because transition-pair parallelization has not been implemented.

## MPI validation and degenerate manifolds

The tracked validation package is in
`testing/validation/rixs_mpi`. It exercises 1-, 2-, and 4-rank pure k-point
layouts and the shared-k-point subgroup-root path. The validation requires
MPI-invariant spectra, complete transition coverage, contribution-spectrum
reconstruction, and correct row arithmetic.

Individual band-labelled contribution rows are not physical invariants when
an eigensolver chooses different bases inside degenerate valence or
intermediate subspaces. A unitary rotation can redistribute
`amplitude_abs2` and `weighted_strength` among those rows while leaving the
spectrum unchanged. The validator therefore:

1. Matches rows by the discrete identity `ikpt`, valence band, intermediate
   band, absorber atom, and absorber type.
2. Compares `loss_energy_Ha` separately within `1e-12` Ha.
3. Constructs valence and intermediate manifolds independently from
   `eps_v_Ha` and `eps_n_Ha`, using a fixed `1e-10` Ha degeneracy threshold.
4. Requires identical manifold membership and a complete Cartesian product of
   valence and intermediate states in every compared layout.
5. Compares manifold-summed `amplitude_abs2` and `weighted_strength` with the
   unchanged spectrum comparison tolerance.

Raw per-band differences and their worst transition identities remain in the
report as diagnostics. Manifolds are never chosen from loss energy alone.

The same point applies to finite-Q pair and site rows: compare degenerate
manifold sums and the final spectrum across MPI layouts, not arbitrary
eigenvector labels inside a degenerate subspace. The physical spectrum,
k-point coverage, and coherent site reconstruction remain the invariants.

## Setup summary

For the legacy path, the setup summary includes:

```text
Approximation            : direct same-k independent-particle
Scalar spin treatment    : spin-degenerate S1 trace
Symmetry treatment       : full-k only, no star reconstruction
Valence band window      : all
Intermediate band window : all
```

The finite-Q path additionally prints the full `Q_rlu`, the reduced electronic
`q`, and a reminder that absorber-site amplitudes are summed coherently before
squaring. Treat those two momentum lines as a production-run check: the first
must match the experimental transfer expressed in the FLEUR reciprocal basis,
and the second must match the intended electronic k-point displacement.

For active band windows the summary prints the requested range, for example:

```text
Valence band window      : 1 ... 20
Intermediate band window : 21 ... 80
```

## Ir L3-edge applications

Ir L3-edge materials, including CaIrO3-type applications, can use the validated
first-variation noco spinor path. The finite-Q path adds the coherent basis-site
structure factor and commensurate k-point displacement needed for momentum-
resolved calculations. Second-variation SOC and many-body RIXS physics remain
outside the implementation. Legacy state-character sidecars can provide local
`j_eff` diagnostics, but they are not available in the finite-Q branch and are
not themselves RIXS matrix-element weights.

## Practical validation checklist

For a new setup, check:

1. The setup summary matches the XML input.
2. The guard condition `nkpt == nkptf` is satisfied.
3. The requested edge/core state is found.
4. For cubic nonmagnetic scalar test systems, diagonal spectra such as `xx`,
   `yy`, and `zz` should obey the expected symmetry.
5. All spectrum intensities are non-negative up to numerical roundoff.
6. `omegaIn` was chosen using the calculated XAS spectrum or transition output,
   not only from an isolated nominal core-edge estimate.
7. Vary `gammaCore` and check that the spectra and integrated intensities
   respond smoothly.
8. Check convergence with respect to k mesh, band windows, `etaLoss`, and the
   loss-energy window.
9. If contribution tables are enabled, the contribution-spectrum check reports
   `PASS`.
10. If band windows are used, the `band_v` and `band_n` columns in contribution
   tables lie inside the requested windows, after any per-k clamping.
11. For finite-Q, the experimental HKL axes have been mapped explicitly to the
    FLEUR reciprocal basis and the printed full Q and reduced q are correct.
12. For finite-Q, the full mesh is closed under `k -> k+q`; half-grid transfers
    require even mesh dimensions along the shifted directions.
13. For a focused finite-Q audit with `writeContributions="T"`, the complex sum
    of all site partials reconstructs each pair amplitude and the pair table
    reconstructs the final Gaussian spectrum.
14. Sum unanalysed outgoing channels only after squaring.

## Interpretation notes

This is an independent-particle implementation. Absolute
core-edge energies generally require alignment to a reference calculation or to
experiment.

The current intensity is in arbitrary units. It is not normalized per unit-cell
volume, per absorber atom, or for a specific experimental scattering geometry.
Use relative trends only after checking convergence with respect to k mesh, band
window, `gammaCore`, `etaLoss`, and the loss-energy window.

Contribution tables are rank-local and can become very large, especially for
dense k meshes, broad band windows, and many incoming/outgoing polarization
pairs. Enable `writeContributions` only for runs where this diagnostic detail is
needed.

For finite-Q, `momentumTransfer` supplies a lattice-scale basis-site phase but
the dipole operator remains local within each muffin tin. Results should not be
interpreted as including the radial variation of `exp(i q.r)`, quadrupole
transitions, a core-hole potential, excitonic binding, or atomic multiplets.
