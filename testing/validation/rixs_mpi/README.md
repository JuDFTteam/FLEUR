# First-variation spinor RIXS MPI validation

This package validates the transition-resolved spinor RIXS path with the same
one-atom SQAy Fe L3 calculation at 1, 2, and 4 MPI ranks. It compares production
spectra and combined rank-local contribution tables in both FLEUR MPI layouts:

- pure k-point parallelism, with one PE per k-point subgroup;
- shared-k-point execution, with all ranks collaborating in one subgroup.

Spinor RIXS is validated in serial and across 1-, 2-, and 4-rank pure k-point
MPI layouts, with additional successful validation of the shared-k-point
subgroup-root path. Spectra, transition coverage, reconstruction, and
degenerate-manifold-summed contributions are MPI invariant.

## Test definition

All three runs use:

- first-variation noco+SOC;
- a full 3x3x3 k-point list (`nkpt == nkptf == 27`);
- `absorberZ="26"` and `edge="L3"`;
- `omegaIn="25.2776419634"` Ha;
- `gammaCore="0.20"` Ha;
- a 0–4 Ha loss grid with 401 points and `etaLoss="0.03"` Ha;
- incoming x and outgoing x polarization;
- `writeContributions="T"`;
- one OpenMP thread per MPI rank.

The three directories contain independent copies of the same density restart,
k-point list, and symmetry file. Their `inp.xml` files differ only in
`outputPrefix`. The old fixture-local schemas were intentionally not copied
because they predate the RIXS XML extension; use the schema distributed with
the current FLEUR build.

The common `cdn.hdf` restart is identical, but FLEUR regenerates eigenvectors
during each one-shot job. Individual band-labelled contributions are not
invariant under rotations inside degenerate valence or intermediate subspaces.
The comparison therefore retains raw row comparisons as diagnostics and tests
the physical `amplitude_abs2` and `weighted_strength` invariants after summing
complete valence-manifold × intermediate-manifold Cartesian products.

## Directory layout

```text
rixs_mpi/
├── VALIDATION_RESULTS.md
├── SQAy_np1/
├── SQAy_np2/
├── SQAy_np4/
├── scripts/
│   ├── run_common.sh
│   ├── run_np1.slurm
│   ├── run_np2.slurm
│   ├── run_np4.slurm
│   ├── submit_all.sh
│   ├── compare_spinor_rixs_mpi.py
│   └── diagnose_spinor_rixs_manifolds.py
└── reference/
    ├── SQAy_stage2_serial_spectrum.dat
    └── seed_sha256.txt
```

`reference/seed_sha256.txt` records the common restart/input checksums and the
checksums of the three run-specific inputs.

## JURECA conventions found in the repository

The repository provides these usable conventions:

- FLEUR's Slurm tooling defaults to `srun fleur_MPI`.
- Generated Slurm jobs use explicit nodes, tasks, CPUs per task, output, and
  error directives.
- `cmake/machines/JURECA_DC_GCC.sh` describes a GCC/ParaStationMPI environment.
- `cmake/machines/JURECA_INTEL.sh` describes an Intel MPI toolchain.
- `cmake/machines/JURECA_GPU.sh` describes an NVHPC/OpenMPI environment.

Those machine files contain toolchain examples, but they do not establish the
current account, partition, or executable path for this validation. This
package therefore does not load modules or hard-code site allocation details.
Load the module environment matching the MPI-enabled FLEUR executable that you
intend to test.

No repository script provided a current JURECA account or partition.

## Required variables

Set these variables in the shell used for submission:

```bash
export FLEUR_BIN=/absolute/path/to/the/JURECA/MPI/fleur_MPI
export SLURM_ACCOUNT=your_jureca_allocation
export SLURM_PARTITION=your_jureca_partition
export SLURM_TIME=HH:MM:SS
```

`FLEUR_BIN` must be the MPI-enabled executable built from a source state that
contains the RIXS subgroup-root ownership fix. The scripts deliberately do not
guess its location.

The established PBE/noco SQAy fixture raises FLEUR's pre-existing
"Noco should be used only with LDA" warning. The run scripts pass `-warn_only`
only for that known warning. Other nonzero exits still terminate the job.

## MPI layouts and submission

FLEUR uses two MPI levels. `fmpi%irank/isize` describe the global communicator,
whereas `fmpi%n_rank/n_size` describe the ranks collaborating on a k point.
Every rank in such a subgroup sees the same `fmpi%k_list`. RIXS currently does
serial work per k point, so only `fmpi%n_rank == 0` evaluates the RIXS
transitions; all ranks still enter the global spectrum reductions.

The launch scripts select pure k-point parallelism by default:

```bash
export RIXS_MPI_LAYOUT=pure-kpoint
scripts/submit_all.sh
```

This passes `-pe_per_kpt 1`. To exercise the subgroup-root path explicitly, use
a clean copy of the package and run:

```bash
export RIXS_MPI_LAYOUT=shared-kpoint
scripts/submit_all.sh
```

For the shared layout, the scripts pass `-pe_per_kpt N` to the N-rank job. The
np2 and np4 jobs therefore use one k-point subgroup containing all ranks. The
scripts refuse to overwrite an earlier layout's generated output, so validate
the two layouts in separate clean copies or remove generated files first.

From this package directory, submit all three jobs with:

```bash
scripts/submit_all.sh
```

The exact individual commands are:

```bash
sbatch --account="$SLURM_ACCOUNT" --partition="$SLURM_PARTITION" \
       --time="$SLURM_TIME" --export=ALL scripts/run_np1.slurm

sbatch --account="$SLURM_ACCOUNT" --partition="$SLURM_PARTITION" \
       --time="$SLURM_TIME" --export=ALL scripts/run_np2.slurm

sbatch --account="$SLURM_ACCOUNT" --partition="$SLURM_PARTITION" \
       --time="$SLURM_TIME" --export=ALL scripts/run_np4.slurm
```

Each job uses one node. The scripts request 1, 2, or 4 tasks respectively and
one CPU per task. They set `OMP_NUM_THREADS=1` and the corresponding common
BLAS thread variables to one, and translate `RIXS_MPI_LAYOUT` into the explicit
`-pe_per_kpt` value described above.

Each script:

- derives the package and run directories from its own location;
- checks `FLEUR_BIN` and the allocated task count;
- records SLURM variables, loaded modules, launcher information, and dynamic
  libraries in `job_environment_JOBID.log`;
- captures FLEUR stdout and stderr separately;
- refuses to overwrite existing RIXS spectra or contribution files;
- checks that one spectrum, the expected number of rank-local tables, and a
  passing contribution-spectrum status were produced.

## Expected output

The run-specific prefixes are:

- `rixs_mpi_np1`
- `rixs_mpi_np2`
- `rixs_mpi_np4`

Each directory should contain one spectrum such as:

```text
rixs_mpi_npN_L3_x_x_omega25p277642_rixs.dat
```

and exactly N contribution files:

```text
rixs_mpi_npN_L3_x_x_omega25p277642_contrib_rank0000.dat
...
```

The run directory also receives:

```text
fleur_stdout.log
fleur_stderr.log
job_environment_JOBID.log
```

Slurm's `slurm-rixs-sqay-npN-JOBID.out` and `.err` files are written in the
directory from which `sbatch` is invoked.

## Comparison

After all three jobs finish, run from the package directory:

```bash
python3 scripts/compare_spinor_rixs_mpi.py --root . \
  | tee mpi_comparison_report.txt
```

The script uses only the Python standard library. It exits with:

- 0 when all required checks pass;
- 1 for a validation failure;
- 2 for missing or malformed inputs/outputs.

It checks:

- production grids, spectra, trapezoidal integrals, rectangular sums, extrema,
  negative values, and non-finite values;
- the expected number of rank-local files and rank IDs;
- row counts per rank and combined;
- complex-amplitude identity and weighted-strength arithmetic;
- non-negative `amplitude_abs2` and `weighted_strength`;
- per-rank and combined ikpt sets;
- missing k points and duplicated full transition keys;
- exact cross-run discrete transition coverage and a separate tolerant
  `loss_energy_Ha` comparison;
- raw band-labelled amplitude/strength differences as diagnostics, including
  their worst transition identities;
- independently constructed valence and intermediate manifolds using
  `eps_v_Ha` and `eps_n_Ha` with a fixed `1e-10` Ha degeneracy threshold;
- identical cross-layout manifold membership, complete Cartesian products,
  and manifold-summed `amplitude_abs2` and `weighted_strength`;
- the spinor setup summary, contribution-spectrum result, abort text, MPI
  errors, and NaN/Inf text in logs;
- agreement with the preserved Stage 2 serial spectrum.

Rows are matched exactly by the discrete identity:

```text
ikpt, band_v, band_n, absorber_atom, absorber_type
```

`loss_energy_Ha` is compared separately with an absolute tolerance of `1e-12`
Ha and is never used to construct a manifold. Repeated `ikpt` values are
expected. Duplication is checked with the same discrete transition identity.

Valence and intermediate manifolds are constructed independently for every
`ikpt, absorber_atom, absorber_type` context. States belong to the same
manifold when their corresponding `eps_v_Ha` or `eps_n_Ha` eigenvalues span no
more than `1e-10` Ha. Acceptance requires identical membership across layouts
and every expected valence-manifold x intermediate-manifold row.

The separate diagnostic can be used to inspect neighboring bands and alternative
degeneracy thresholds without changing acceptance:

```bash
python3 scripts/diagnose_spinor_rixs_manifolds.py --root .
```

## Acceptance criteria

The package requires:

- successful completion for all MPI sizes;
- production spectra satisfying
  `max_abs_diff < 1e-12 OR max_rel_diff < 1e-10`;
- `PASS` from every contribution-spectrum reconstruction;
- identical combined row counts;
- identical discrete transition coverage, no missing k points, and no
  duplicated transition keys;
- `loss_energy_Ha` agreement within `1e-12` Ha;
- amplitude identity and weighted arithmetic satisfying
  `max_abs_error < 1e-18 OR max_rel_error < 1e-12`;
- identical degenerate-manifold membership and complete Cartesian products;
- manifold-summed `amplitude_abs2` and `weighted_strength` satisfying the
  unchanged comparison tolerance of `max_abs_diff < 1e-12 OR
  max_rel_diff < 1e-10`;
- no spectrum values below `-1e-14`;
- no non-finite values.

MPI reduction order can produce ordinary floating-point differences, and
regenerated eigenvectors can rotate inside degenerate manifolds. The script
reports raw band-labelled differences but does not treat those gauge-dependent
rows as acceptance failures when complete manifold sums remain invariant.

## Three different integrated quantities

The report deliberately distinguishes:

1. The trapezoidal integral of the plotted finite loss grid.
2. FLEUR's rectangular finite-grid diagnostic:
   `sum(intensity) * constant_grid_spacing`.
3. `sum(weighted_strength)` over transition rows.

The first two differ by endpoint treatment. The third is the total weight of
normalized Gaussian transitions before truncation to the plotted 0–4 Ha
window, so it can differ because Gaussian tails and some transition centers
lie outside that finite window.

## Limitations

- This package validates one-node 1/2/4-rank pure k-point execution and the
  shared-k-point subgroup-root path.
- It does not test symmetry-star RIXS, second-variation SOC, coherent
  multi-atom amplitudes, GPUs, or multiple OpenMP threads per rank.
- The preserved serial spectrum was generated with a different local
  toolchain; tolerance-level rather than byte-level agreement is required on
  JURECA.
- No scalar MPI package is included because the new spinor contribution path
  is the priority and scalar MPI behavior was previously validated.
