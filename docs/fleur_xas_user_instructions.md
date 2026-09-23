# FLEUR XAS input: user instructions

This note describes the first XML-controlled XAS interface implemented in the FLEUR prototype.

## Current scope

Implemented and validated:

- Independent-particle electric-dipole XAS.
- Edges: `K`, `L2`, `L3`.
- Linear polarizations: `x`, `y`, `z`.
- Absorber selection by chemical atomic number `absorberZ`.
- Scalar/collinear calculations.
- First-variation noncollinear SOC calculations, including unitary spatial-star symmetry reconstruction.
- MPI reduction path for the tested four-rank case.

Not implemented yet:

- Combined `L23` edge output. Run `L2` and `L3` separately for now.
- Second-variation SOC XAS path.
- Time-reversal star members for SOC/noncollinear XAS.
- Circular polarization as an XML production option.
- Core-hole, multiplet, and BSE-like many-body effects.

## Minimal XML input

Add the XAS block under `<output>`:

```xml
<output>
   <xas l_xas="T"
        absorberZ="26"
        edge="L3"
        eta="0.030"
        nEnergy="401"
        polarizations="x y z"
        outputPrefix="xas"/>
</output>
```

If no `<xas>` section is present, XAS is off and the calculation should behave as before.

## Optional explicit energy window

By default, the XAS driver determines the spectrum window automatically from the available transition energies. To override this, provide `eMin` and `eMax` in Hartree:

```xml
<xas l_xas="T"
     absorberZ="26"
     edge="L3"
     eta="0.030"
     nEnergy="401"
     eMin="-0.500"
     eMax="1.500"
     polarizations="x y z"
     outputPrefix="xas"/>
```

The output prints both the input request and the resolved energy grid actually used. For automatic windows, the output reports the final numerical `eMin/eMax` in Hartree and eV.

## Attributes

| Attribute | Required? | Meaning | Example |
|---|---:|---|---|
| `l_xas` | no | Enables or disables XAS. If the `<xas>` element is present, the intended use is `T`. | `T` |
| `absorberZ` | yes when XAS is enabled | Atomic number of absorber atoms. All matching atom types are summed. | `26` |
| `edge` | no | Absorption edge. Supported: `K`, `L2`, `L3`. Default: `L3`. | `L3` |
| `eta` | no | Gaussian broadening in Hartree. Default: `0.030`. | `0.030` |
| `nEnergy` | no | Number of points in the output spectrum. Default: `401`. | `401` |
| `eMin`, `eMax` | no | Explicit energy window in Hartree. If omitted, automatic transition range is used. | `-0.5`, `1.5` |
| `polarizations` | no | Space-separated list of linear polarizations. Supported now: `x y z`. | `x y z` |
| `outputPrefix` | no | Prefix for spectrum output files. Default: `xas`. | `xas` |

## Supported edge aliases

- `K`, `k`, `1s1/2`, `1S1/2`
- `L2`, `l2`, `2p1/2`, `2P1/2`
- `L3`, `l3`, `2p3/2`, `2P3/2`

`L23` currently aborts intentionally with a clear message. Use separate `L2` and `L3` runs.

## Output files

The output filename includes the prefix, edge, polarization, and broadening, for example:

```text
xas_L3_x_eta0p030.dat
xas_L3_y_eta0p030.dat
xas_L3_z_eta0p030.dat
```

Each file contains the broadened XAS spectrum for the selected edge and polarization.

## Output summary

When XAS is enabled, FLEUR prints a setup summary similar to:

```text
 ---------- XAS setup -------------------------------
 XAS enabled              : T
 Absorber Z               : 26
 Edge                     : L3
 Broadening eta           :   0.030000 Ha
 Number of energy points  : 401
 Energy window input      : automatic
 Polarizations            : x y z
 Output prefix            : xas
 ----------------------------------------------------
```

After the grid has been resolved, it prints the actual energy grid:

```text
 ---------- XAS resolved energy grid ----------------
 Energy window used       :    25.107212 ...    34.474637 Ha  (  683.202057 ...   938.102655 eV)
 Energy step              :     0.023419 Ha  (    0.637251 eV)
 Number of energy points  : 401
 Window source            : automatic transition range
 ----------------------------------------------------
```

For explicit `eMin/eMax`, the resolved window should match the input values.

## Practical validation checklist

Before trusting a new setup, check these lines in the output:

1. XAS setup summary matches your XML input.
2. Resolved energy window is reasonable.
3. The correct core state is found for the requested edge.
4. For high-symmetry nonmagnetic scalar tests, `Sx`, `Sy`, and `Sz` should be equal.
5. For symmetry-reduced runs, the star-weight summary should satisfy:

```text
selected = 1.0000000000E+00
expanded = 1.0000000000E+00
diff     = 0.0000E+00
```

## Notes for noncollinear SOC

The validated SOC path is the first-variation noncollinear SOC path, where both `l_noco` and `l_soc` are true and two local spinor components are available. The implementation keeps protective guards for unsupported SOC variants.

For inpgen-generated noncollinear SOC inputs, remember that the inpgen `&soc theta phi /` line is used for symmetry generation. In noco runtime input, the local spin direction is stored in atom-group `nocoParams` through `alpha` and `beta`. At present, the inpgen workflow may require care so that non-noco runtime SOC `theta/phi` values do not trigger the noco warning.

## Current limitations

- The calculation is independent-particle XAS: it does not include a self-consistent core hole, multiplet effects, electron-hole attraction, or lifetime physics beyond the chosen broadening.
- Absolute edge positions usually need alignment to experiment or a reference calculation.
- The output is best interpreted first in terms of relative spectral shapes and polarization dependence.
