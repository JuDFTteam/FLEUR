#!/usr/bin/env bash

set -euo pipefail

if [[ $# -ne 3 ]]; then
    echo "Usage: $0 RUN_DIRECTORY EXPECTED_MPI_RANKS OUTPUT_PREFIX" >&2
    exit 2
fi

RUN_DIRECTORY=$1
EXPECTED_MPI_RANKS=$2
OUTPUT_PREFIX=$3

RIXS_MPI_LAYOUT=${RIXS_MPI_LAYOUT:-pure-kpoint}
case "${RIXS_MPI_LAYOUT}" in
    pure-kpoint)
        PES_PER_KPT=1
        ;;
    shared-kpoint)
        PES_PER_KPT=${EXPECTED_MPI_RANKS}
        ;;
    *)
        echo "RIXS_MPI_LAYOUT must be pure-kpoint or shared-kpoint, got: ${RIXS_MPI_LAYOUT}" >&2
        exit 2
        ;;
esac

SCRIPT_DIR=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)
PACKAGE_ROOT=$(cd -- "${SCRIPT_DIR}/.." && pwd)
if [[ "${RUN_DIRECTORY}" = /* ]]; then
    RUN_DIR=${RUN_DIRECTORY}
else
    RUN_DIR="${PACKAGE_ROOT}/${RUN_DIRECTORY}"
fi

: "${FLEUR_BIN:?Set FLEUR_BIN to the JURECA MPI-enabled FLEUR executable.}"

if [[ ! -x "${FLEUR_BIN}" ]]; then
    echo "FLEUR_BIN is not executable: ${FLEUR_BIN}" >&2
    exit 2
fi
if [[ ! -d "${RUN_DIR}" ]]; then
    echo "Validation run directory not found: ${RUN_DIR}" >&2
    exit 2
fi
if [[ ${SLURM_NTASKS:-${EXPECTED_MPI_RANKS}} -ne ${EXPECTED_MPI_RANKS} ]]; then
    echo "Expected ${EXPECTED_MPI_RANKS} SLURM tasks, got ${SLURM_NTASKS:-unset}." >&2
    exit 2
fi

export OMP_NUM_THREADS=1
export MKL_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export BLIS_NUM_THREADS=1

JOB_TAG=${SLURM_JOB_ID:-manual}
ENV_LOG="${RUN_DIR}/job_environment_${JOB_TAG}.log"
STDOUT_LOG="${RUN_DIR}/fleur_stdout.log"
STDERR_LOG="${RUN_DIR}/fleur_stderr.log"

exec > >(tee -a "${ENV_LOG}") 2>&1

on_exit() {
    rc=$?
    trap - EXIT
    echo "Finished at $(date --iso-8601=seconds 2>/dev/null || date)"
    echo "Exit status: ${rc}"
    exit "${rc}"
}
trap on_exit EXIT

shopt -s nullglob
old_spectra=("${RUN_DIR}/${OUTPUT_PREFIX}"*_rixs.dat)
old_contributions=("${RUN_DIR}/${OUTPUT_PREFIX}"*_contrib_rank*.dat)
if (( ${#old_spectra[@]} > 0 || ${#old_contributions[@]} > 0 )); then
    echo "Refusing to overwrite existing RIXS validation output in ${RUN_DIR}."
    echo "Restore a clean validation directory before resubmitting."
    exit 2
fi

echo "RIXS MPI validation run: ${RUN_DIRECTORY}"
echo "Started at $(date --iso-8601=seconds 2>/dev/null || date)"
echo "Host: $(hostname)"
echo "Working directory: ${RUN_DIR}"
echo "FLEUR executable: ${FLEUR_BIN}"
echo "Expected MPI ranks: ${EXPECTED_MPI_RANKS}"
echo "RIXS MPI layout: ${RIXS_MPI_LAYOUT}"
echo "PEs per k-point subgroup: ${PES_PER_KPT}"
echo "OMP_NUM_THREADS=${OMP_NUM_THREADS}"
echo "SLURM variables:"
env | sort | grep '^SLURM_' || true
echo "Loaded modules:"
if command -v module >/dev/null 2>&1; then
    module -t list 2>&1 || true
else
    echo "Environment-modules command is unavailable."
fi
echo "Launcher:"
command -v srun
srun --version | head -n 1 || true
echo "FLEUR dynamic libraries:"
if command -v ldd >/dev/null 2>&1; then
    ldd "${FLEUR_BIN}" || true
fi

cd "${RUN_DIR}"

echo "Launching: srun --ntasks=${EXPECTED_MPI_RANKS} --cpus-per-task=1 ${FLEUR_BIN} -warn_only -pe_per_kpt ${PES_PER_KPT}"
srun --ntasks="${EXPECTED_MPI_RANKS}" --cpus-per-task=1 \
    "${FLEUR_BIN}" -warn_only -pe_per_kpt "${PES_PER_KPT}" >"${STDOUT_LOG}" 2>"${STDERR_LOG}"

spectra=("${RUN_DIR}/${OUTPUT_PREFIX}"*_rixs.dat)
contributions=("${RUN_DIR}/${OUTPUT_PREFIX}"*_contrib_rank*.dat)
if (( ${#spectra[@]} != 1 )); then
    echo "Expected one production spectrum, found ${#spectra[@]}."
    exit 3
fi
if (( ${#contributions[@]} != EXPECTED_MPI_RANKS )); then
    echo "Expected ${EXPECTED_MPI_RANKS} rank-local contribution files, found ${#contributions[@]}."
    exit 3
fi
if ! grep -q "status       = PASS" "${STDOUT_LOG}"; then
    echo "The FLEUR log does not contain a passing contribution-spectrum check."
    exit 3
fi

echo "Spectrum: ${spectra[0]}"
echo "Rank-local contribution files: ${#contributions[@]}"
echo "FLEUR stdout: ${STDOUT_LOG}"
echo "FLEUR stderr: ${STDERR_LOG}"
echo "Run completed successfully."
