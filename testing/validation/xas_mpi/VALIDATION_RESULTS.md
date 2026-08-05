# XAS MPI ownership validation results

Status: **PASS**

The production fix is commit
`b3818dad42aa56fa3a835a9cf8e972e23ea2dd87` (`Fix duplicated XAS work in
shared k-point MPI groups`).

The approved clean JURECA validation used executable SHA256
`df4af8292c337614a5755910c7e59634fe638c7cf1ff076127dd501d577e65c4`
and identical isolated copies of the tracked SQAy fixture and restart. Source
status was clean in every job. The canonical x/y/z spectra are copied from the
`np1` run and recorded in `reference/canonical/SHA256SUMS`; complete build,
module, linkage, input, restart, launcher, job, and result provenance is stored
in `provenance/clean_jureca_manifest.json`.

## Accepted layouts

| Layout | MPI ranks | `-pe_per_kpt` | Transition rows per polarization/rank |
|---|---:|---:|---|
| `np1` | 1 | 1 | 2810 |
| `np2_pure` | 2 | 1 | disjoint pure-k-point subsets; 2810 total |
| `np4_pure` | 4 | 1 | disjoint pure-k-point subsets; 2810 total |
| `np4_shared` | 4 | 2 | 2144, 0, 666, 0; 2810 total |

Every layout completed successfully. For each x/y/z polarization the combined
transition count is 2810, with no duplicate discrete transition keys, complete
k-point coverage 1--27, and complete parent-star coverage 1--6. The collaborator
ranks in `np4_shared` retain valid header-only transition files.

Parent and star weights equal one in every layout. Additive strengths,
underflow counts, transition-key sets, and weighted transition strengths are
MPI-layout invariant. Transition weighted-strength arithmetic errors are zero.

## Numerical comparison

Across all layouts relative to `np1`:

- maximum energy-grid difference: `2.061e-13 Ha`;
- maximum absolute spectrum difference: `2.290e-15`;
- maximum relative spectrum difference: `1.466e-13`;
- maximum weighted-strength-sum difference: `3.166e-17`;
- transition arithmetic errors: zero.

The complete command

```bash
python3 scripts/compare_xas_mpi.py | tee xas_mpi_comparison_report.txt
```

finished with `OVERALL STATUS: PASS` without changing any tolerance.

The historical July 2026 macOS reference remains provenance only. Its tiny
automatic-window grid drift is compiler/platform dependent and is not used as
the canonical acceptance reference.
