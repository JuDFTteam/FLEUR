# Scalar independent-particle RIXS prototype

This note describes the current XML-controlled, guarded RIXS prototype in FLEUR.
It is intended for early independent-particle calculations and diagnostics, not
as a complete production RIXS theory.

## Current scope

The implementation evaluates a same-k, vertical Kramers-Heisenberg-like
independent-particle RIXS spectrum for selected absorber atoms, edge, and
linear incoming/outgoing polarizations.

Currently supported:

- scalar `jspins=1` calculations only;
- edges `K`, `L2`, and `L3`;
- linear polarizations `x`, `y`, and `z`;
- absorber selection by `absorberZ`;
- optional rank-local contribution tables;
- optional valence/intermediate band windows;
- MPI reduction of the broadened spectrum.

Current guards and limitations:

- no SOC or noco RIXS yet;
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

When contribution output is enabled, FLEUR also prints a
contribution-to-spectrum consistency check. It should report `PASS`.

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

## Practical validation checklist

For a new scalar setup, check:

1. The setup summary matches the XML input.
2. The guard condition `nkpt == nkptf` is satisfied.
3. The requested edge/core state is found.
4. For cubic nonmagnetic scalar test systems, diagonal spectra such as `xx`,
   `yy`, and `zz` should obey the expected symmetry.
5. If contribution tables are enabled, the contribution-spectrum check reports
   `PASS`.
6. If band windows are used, the `band_v` and `band_n` columns in contribution
   tables lie inside the requested windows, after any per-k clamping.

## Interpretation notes

This is an independent-particle diagnostic/prototype implementation. Absolute
energies usually require alignment to a reference or experiment. Intensities and
polarization trends are most useful after checking convergence with respect to
k mesh, band window, `gammaCore`, and `etaLoss`.
