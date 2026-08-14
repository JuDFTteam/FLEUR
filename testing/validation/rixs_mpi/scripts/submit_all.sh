#!/usr/bin/env bash

set -euo pipefail

: "${FLEUR_BIN:?Set FLEUR_BIN to the absolute path of the MPI-enabled FLEUR executable.}"
: "${SLURM_ACCOUNT:?Set SLURM_ACCOUNT to your JURECA allocation.}"
: "${SLURM_PARTITION:?Set SLURM_PARTITION to the intended JURECA partition.}"
: "${SLURM_TIME:?Set SLURM_TIME to the requested wall time, for example HH:MM:SS.}"

if [[ ! -x "${FLEUR_BIN}" ]]; then
    echo "FLEUR_BIN is not executable: ${FLEUR_BIN}" >&2
    exit 2
fi

SCRIPT_DIR=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)
PACKAGE_ROOT=$(cd -- "${SCRIPT_DIR}/.." && pwd)
cd "${PACKAGE_ROOT}"

export FLEUR_BIN
export RIXS_MPI_PACKAGE_ROOT="${PACKAGE_ROOT}"

for ranks in 1 2 4; do
    echo "Submitting SQAy ${ranks}-rank validation"
    sbatch \
        --account="${SLURM_ACCOUNT}" \
        --partition="${SLURM_PARTITION}" \
        --time="${SLURM_TIME}" \
        --chdir="${PACKAGE_ROOT}" \
        --export=ALL \
        "scripts/run_np${ranks}.slurm"
done
