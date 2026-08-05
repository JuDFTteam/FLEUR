#!/usr/bin/env bash

set -euo pipefail

: "${FLEUR_BIN:?Set FLEUR_BIN to the clean MPI-enabled FLEUR executable.}"
: "${FLEUR_SOURCE_DIR:?Set FLEUR_SOURCE_DIR to the clean Phase-1 checkout.}"
: "${SLURM_ACCOUNT:?Set SLURM_ACCOUNT to your JURECA allocation.}"
: "${SLURM_PARTITION:?Set SLURM_PARTITION to the intended JURECA partition.}"
: "${SLURM_TIME:?Set SLURM_TIME, for example 00:10:00.}"

SCRIPT_DIR=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)
PACKAGE_ROOT=$(cd -- "${SCRIPT_DIR}/.." && pwd)
EXPECTED_COMMIT=$(tr -d '[:space:]' < "${PACKAGE_ROOT}/provenance/PHASE1_COMMIT")
if [[ ! -d "${PACKAGE_ROOT}/runs" ]]; then
    echo "Run scripts/prepare_runs.sh before submission." >&2
    exit 2
fi
if [[ ! -x "${FLEUR_BIN}" ]]; then
    echo "FLEUR_BIN is not executable: ${FLEUR_BIN}" >&2
    exit 2
fi

source_commit=$(git -C "${FLEUR_SOURCE_DIR}" rev-parse HEAD)
source_status=$(git -C "${FLEUR_SOURCE_DIR}" status --short)
if [[ "${source_commit}" != "${EXPECTED_COMMIT}" ]]; then
    echo "Expected source commit ${EXPECTED_COMMIT}, found ${source_commit}." >&2
    exit 2
fi
if [[ -n "${source_status}" ]]; then
    echo "The FLEUR source checkout is not clean:" >&2
    printf '%s\n' "${source_status}" >&2
    exit 2
fi

fleur_version=$(cd -- "$(dirname -- "${FLEUR_BIN}")" && "${FLEUR_BIN}" -version 2>&1)
if ! grep -q "${EXPECTED_COMMIT}" <<<"${fleur_version}"; then
    echo "FLEUR version output does not contain expected commit ${EXPECTED_COMMIT}." >&2
    exit 2
fi

(cd "${PACKAGE_ROOT}/fixtures" && sha256sum --check SHA256SUMS)
for label in np1 np2_pure np4_pure np4_shared; do
    run_dir="${PACKAGE_ROOT}/runs/${label}"
    if [[ ! -d "${run_dir}" || ! -f "${run_dir}/SEED_SHA256" ]]; then
        echo "Prepared run directory or seed manifest is missing: ${run_dir}" >&2
        exit 2
    fi
    if ! diff -u "${PACKAGE_ROOT}/fixtures/SHA256SUMS" "${run_dir}/SEED_SHA256"; then
        echo "Prepared seed manifest differs for ${label}." >&2
        exit 2
    fi
    (cd "${run_dir}" && sha256sum --check SEED_SHA256)
    unexpected=$(find "${run_dir}" -maxdepth 1 -type f \
        ! -name cdn.hdf ! -name inp.xml ! -name kpts.xml ! -name sym.xml ! -name SEED_SHA256 -print)
    if [[ -n "${unexpected}" ]]; then
        echo "Refusing to submit ${label}; the run directory contains old or unexpected files:" >&2
        printf '%s\n' "${unexpected}" >&2
        exit 2
    fi
done

echo "Submission preflight passed for all four clean run directories."

export XAS_MPI_PACKAGE_ROOT="${PACKAGE_ROOT}"
cd "${PACKAGE_ROOT}"

for script in run_np1.slurm run_np2_pure.slurm run_np4_pure.slurm run_np4_shared.slurm; do
    echo "Submitting ${script} from ${PACKAGE_ROOT}"
    sbatch --account="${SLURM_ACCOUNT}" --partition="${SLURM_PARTITION}" \
        --time="${SLURM_TIME}" --chdir="${PACKAGE_ROOT}" \
        --export=ALL,XAS_MPI_PACKAGE_ROOT="${PACKAGE_ROOT}",FLEUR_BIN="${FLEUR_BIN}",FLEUR_SOURCE_DIR="${FLEUR_SOURCE_DIR}" \
        "scripts/${script}"
done
