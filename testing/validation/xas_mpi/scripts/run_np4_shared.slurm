#!/bin/bash
#SBATCH --job-name=xas-sqay-np4-shared
#SBATCH --nodes=1
#SBATCH --ntasks=4
#SBATCH --cpus-per-task=1
#SBATCH --output=slurm-%x-%j.out
#SBATCH --error=slurm-%x-%j.err

set -euo pipefail
: "${XAS_MPI_PACKAGE_ROOT:=${SLURM_SUBMIT_DIR:?Submit from the package root or use submit_all.sh.}}"
exec bash "${XAS_MPI_PACKAGE_ROOT}/scripts/run_common.sh" np4_shared 4 2
