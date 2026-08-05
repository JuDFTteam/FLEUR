# Legacy macOS reference provenance

These checksums describe the historical SQAy XAS spectra generated on
2026-07-08. The spectra are not canonical acceptance references for this MPI
test.

- source commit: `b3f89728e76db77fdb0bb8a4017c1ae59b56212d`, dirty tree;
- FLEUR compilation date: 2026-06-25 09:10:21;
- execution: 2026-07-08 23:13:37--23:13:38 CEST;
- platform: macOS/Darwin ARM64;
- command: `/Users/mlezaic/FLEUR/fleur/build/fleur -warn_only`;
- linked numerical library: OpenBLAS;
- XML library: macOS SDK libxml2;
- OpenMP threads: 8;
- exact compiler version: not recorded;
- original `cdn.hdf` SHA256: `05feb810fcce3d04e85e1f6e0a35ca4820a6db447d9006732edcd6ebb1c0d535`.

The automatic XAS energy window is derived from numerical transition extrema.
Tiny compiler- and platform-dependent changes in an endpoint can produce a
correspondingly tiny affine change in the energy grid. The July spectra differ
from the August JURECA calculation by up to `6.693e-12` Ha on that grid, while
integrated intensities agree to about `4e-17`. They are retained for provenance
and diagnostic comparison only, not as a bytewise or tolerance-gated reference.
