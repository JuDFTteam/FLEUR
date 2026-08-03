# Independent-particle RIXS prototype

This note describes the current XML-controlled, guarded scalar and
first-variation spinor RIXS prototype in FLEUR. It is intended for early
independent-particle calculations and diagnostics, not as a complete
production RIXS theory.

## Current scope

The implementation evaluates a same-k, vertical Kramers-Heisenberg-like
independent-particle RIXS spectrum for selected absorber atoms, edge, and
linear incoming/outgoing polarizations.

Currently supported:

- scalar `jspins=1` calculations;
- first-variation noco spinor calculations with `jspins=2`;
- edges `K`, `L2`, and `L3`;
- linear polarizations `x`, `y`, and `z`;
- absorber selection by `absorberZ`;
- optional rank-local contribution tables;
- optional valence/intermediate band windows;
- pure k-point and shared-k-point-subgroup MPI layouts;
- MPI reduction of the broadened spectrum.

Current guards and limitations:

- no collinear spin-polarized scalar RIXS;
- no second-variation SOC RIXS; SOC is supported through the validated
  first-variation noco spinor path;
- full-k/no-star k meshes only, i.e. `nkpt == nkptf`;
- no momentum transfer yet;
- no explicit core hole;
- no multiplets;
- no coherent symmetry-star reconstruction yet;
- no pseudospin or `j_eff` analysis yet;
- no combined `L23` edge; run `L2` and `L3` separately.

## Formula

For scalar bands, the production path uses the spin-degenerate S1 trace. The
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

For first-variation spinor bands, FLEUR instead stores the coherent complex
core-`mj` amplitude for each valence/intermediate band pair and local absorber
atom. The final spectrum remains an incoherent sum over distinct band-pair and
absorber-site final states.

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
| `etaLoss` | yes | Gaussian broadening of the final loss spectrum, in Hartree. | `0.05` |
| `incomingPolarizations` | no | Space-separated list of incoming linear polarizations. Default: `x`. | `x y z` |
| `outgoingPolarizations` | no | Space-separated list of outgoing linear polarizations. Default: `x`. | `x y z` |
| `outputPrefix` | no | Prefix for RIXS output files. Default: `rixs`. | `rixs` |
| `writeContributions` | no | Write rank-local valence/intermediate-band contribution tables. Default: `F`. | `T` |
| `valenceBandMin`, `valenceBandMax` | no | Optional 1-based valence-band loop window. | `1`, `20` |
| `intermediateBandMin`, `intermediateBandMax` | no | Optional 1-based intermediate-band loop window. | `21`, `80` |

Supported edge aliases are the same as for XML XAS:

- `K`, `k`, `1s1/2`, `1S1/2`
- `L2`, `l2`, `2p1/2`, `2P1/2`
- `L3`, `l3`, `2p3/2`, `2P3/2`

`L23` intentionally aborts for now.

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

Spectrum files are written on rank 0:

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

If `writeContributions="T"`, one rank-local contribution table is written per
requested polarization pair and MPI rank:

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

## Setup summary

When RIXS is enabled, FLEUR prints a setup summary including:

```text
Approximation            : direct same-k independent-particle
Scalar spin treatment    : spin-degenerate S1 trace
Symmetry treatment       : full-k only, no star reconstruction
Valence band window      : all
Intermediate band window : all
```

For active band windows the summary prints the requested range, for example:

```text
Valence band window      : 1 ... 20
Intermediate band window : 21 ... 80
```

## Ir L3-edge applications

Ir L3-edge materials, including CaIrO3-type applications, can use the validated
first-variation noco spinor path. Second-variation SOC, pseudospin, and `j_eff`
analysis remain outside the current implementation, so interpretation still
requires care beyond the independent-particle spectrum.

## Practical validation checklist

For a new scalar setup, check:

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

## Interpretation notes

This is an independent-particle diagnostic/prototype implementation. Absolute
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
