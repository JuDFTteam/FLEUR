#!/usr/bin/env bash

set -euo pipefail

if [[ $# -ne 3 ]]; then
    echo "Usage: $0 RUN_LABEL MPI_RANKS PE_PER_KPT" >&2
    exit 2
fi

RUN_LABEL=$1
EXPECTED_MPI_RANKS=$2
EXPECTED_PE_PER_KPT=$3

: "${XAS_MPI_PACKAGE_ROOT:?XAS_MPI_PACKAGE_ROOT was not exported by submit_all.sh.}"
: "${FLEUR_SOURCE_DIR:?Set FLEUR_SOURCE_DIR to the clean Phase-1 checkout.}"
: "${FLEUR_BIN:?Set FLEUR_BIN to the clean MPI-enabled FLEUR executable.}"

PACKAGE_ROOT=$(cd -- "${XAS_MPI_PACKAGE_ROOT}" && pwd)
RUN_DIR="${PACKAGE_ROOT}/runs/${RUN_LABEL}"
EXPECTED_COMMIT=$(tr -d '[:space:]' < "${PACKAGE_ROOT}/provenance/PHASE1_COMMIT")
OUTPUT_PREFIX=xas_mpi

if [[ ! -x "${FLEUR_BIN}" ]]; then
    echo "FLEUR_BIN is not executable: ${FLEUR_BIN}" >&2
    exit 2
fi
if [[ ! -d "${RUN_DIR}" || ! -f "${RUN_DIR}/inp.xml" || ! -f "${RUN_DIR}/cdn.hdf" ]]; then
    echo "Prepared run directory is missing: ${RUN_DIR}" >&2
    exit 2
fi
if [[ ${SLURM_NTASKS:-${EXPECTED_MPI_RANKS}} -ne ${EXPECTED_MPI_RANKS} ]]; then
    echo "Expected ${EXPECTED_MPI_RANKS} tasks, got ${SLURM_NTASKS:-unset}." >&2
    exit 2
fi

source_commit=$(git -C "${FLEUR_SOURCE_DIR}" rev-parse HEAD)
source_status=$(git -C "${FLEUR_SOURCE_DIR}" status --short)
if [[ "${source_commit}" != "${EXPECTED_COMMIT}" ]]; then
    echo "Expected clean source commit ${EXPECTED_COMMIT}, found ${source_commit}." >&2
    exit 2
fi
if [[ -n "${source_status}" ]]; then
    echo "The FLEUR source checkout is not clean:" >&2
    printf '%s\n' "${source_status}" >&2
    exit 2
fi

expected_seed=$(awk '$2 == "cdn.hdf" {print $1}' "${PACKAGE_ROOT}/fixtures/SHA256SUMS")
actual_seed=$(sha256sum "${RUN_DIR}/cdn.hdf" | awk '{print $1}')
if [[ "${actual_seed}" != "${expected_seed}" ]]; then
    echo "Restart checksum mismatch before run: expected ${expected_seed}, found ${actual_seed}." >&2
    exit 2
fi

shopt -s nullglob
stale_outputs=()
for candidate in RUN_COMPLETE fleur_stdout.log fleur_stderr.log result_sha256.txt out out.xml; do
    if [[ -e "${RUN_DIR}/${candidate}" ]]; then
        stale_outputs+=("${RUN_DIR}/${candidate}")
    fi
done
generated_outputs=(
    "${RUN_DIR}"/job_environment_*.log
    "${RUN_DIR}/${OUTPUT_PREFIX}"_*_eta*.dat
    "${RUN_DIR}/${OUTPUT_PREFIX}"_*_transitions_rank*.dat
)
stale_outputs+=("${generated_outputs[@]}")
if (( ${#stale_outputs[@]} > 0 )); then
    echo "Refusing to mix a rerun with existing output in ${RUN_DIR}:" >&2
    printf '  %s\n' "${stale_outputs[@]}" >&2
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
FLEUR_MPI_ARGS=(-pe_per_kpt "${EXPECTED_PE_PER_KPT}")

exec > >(tee "${ENV_LOG}") 2>&1

on_exit() {
    rc=$?
    trap - EXIT
    echo "Finished at $(date --iso-8601=seconds 2>/dev/null || date)"
    echo "Exit status: ${rc}"
    exit "${rc}"
}
trap on_exit EXIT

echo "XAS MPI ownership validation"
echo "Started at $(date --iso-8601=seconds 2>/dev/null || date)"
echo "Run label: ${RUN_LABEL}"
echo "Host: $(hostname -f 2>/dev/null || hostname)"
echo "Working directory: ${RUN_DIR}"
echo "FLEUR source: ${FLEUR_SOURCE_DIR}"
echo "Source commit: ${source_commit}"
echo "Source status --short: CLEAN"
echo "FLEUR executable: ${FLEUR_BIN}"
echo "FLEUR executable SHA256: $(sha256sum "${FLEUR_BIN}" | awk '{print $1}')"
echo "FLEUR executable stat:"
stat "${FLEUR_BIN}" || true
echo "Input and original restart checksums before execution:"
(cd "${RUN_DIR}" && sha256sum inp.xml cdn.hdf kpts.xml sym.xml)
echo "Expected MPI ranks: ${EXPECTED_MPI_RANKS}"
echo "Requested PEs per k-point: ${EXPECTED_PE_PER_KPT}"
echo "OMP_NUM_THREADS=${OMP_NUM_THREADS}"
echo "SLURM job variables:"
env | sort | grep '^SLURM_' || true
echo "Loaded modules:"
module -t list 2>&1 || true
echo "Compiler versions:"
mpif90 --version 2>&1 | head -n 4 || true
mpicc --version 2>&1 | head -n 4 || true
if command -v ifx >/dev/null 2>&1; then ifx --version 2>&1 | head -n 4; fi
echo "MPI wrapper configuration:"
mpif90 -show 2>&1 || true
echo "MPI and launcher versions:"
if command -v mpirun >/dev/null 2>&1; then mpirun --version 2>&1 | head -n 6; fi
srun --version 2>&1 | head -n 2 || true
echo "MKL environment:"
echo "MKLROOT=${MKLROOT:-unset}"
module show imkl 2>&1 || true
echo "HDF5 configuration:"
if command -v h5pfc >/dev/null 2>&1; then h5pfc -showconfig 2>&1; elif command -v h5fc >/dev/null 2>&1; then h5fc -showconfig 2>&1; fi
echo "CMake version:"
cmake --version 2>&1 | head -n 2 || true
echo "libxml2 version:"
if command -v xml2-config >/dev/null 2>&1; then xml2-config --version; else module show libxml2 2>&1 || true; fi
echo "FLEUR version output:"
fleur_version=$(cd -- "$(dirname -- "${FLEUR_BIN}")" && "${FLEUR_BIN}" -version 2>&1)
printf '%s\n' "${fleur_version}"
if ! grep -q "${EXPECTED_COMMIT}" <<<"${fleur_version}"; then
    echo "FLEUR version output does not contain expected commit ${EXPECTED_COMMIT}." >&2
    exit 2
fi
echo "Relevant dynamic libraries:"
ldd "${FLEUR_BIN}" | grep -iE 'mpi|blacs|scalapack|hdf5|xml' || true

cd "${RUN_DIR}"
launcher=(srun --ntasks="${EXPECTED_MPI_RANKS}" --cpus-per-task=1 "${FLEUR_BIN}" -warn_only "${FLEUR_MPI_ARGS[@]}")
printf 'Exact launcher command:'
printf ' %q' "${launcher[@]}"
printf '\n'
"${launcher[@]}" >"${STDOUT_LOG}" 2>"${STDERR_LOG}"

spectra=("${RUN_DIR}/${OUTPUT_PREFIX}"_*_eta*.dat)
transitions=("${RUN_DIR}/${OUTPUT_PREFIX}"_*_transitions_rank*.dat)
if (( ${#spectra[@]} != 3 )); then
    echo "Expected three polarization spectra, found ${#spectra[@]}." >&2
    exit 3
fi
if (( ${#transitions[@]} != 3 * EXPECTED_MPI_RANKS )); then
    echo "Expected $((3 * EXPECTED_MPI_RANKS)) transition files, found ${#transitions[@]}." >&2
    exit 3
fi
if ! grep -q "XAS wrote spectrum to" "${STDOUT_LOG}"; then
    echo "XAS completion output was not found." >&2
    exit 3
fi

echo "Restart checksum after execution: $(sha256sum cdn.hdf | awk '{print $1}')"
sha256sum "${OUTPUT_PREFIX}"_*_eta*.dat "${OUTPUT_PREFIX}"_*_transitions_rank*.dat \
    fleur_stdout.log fleur_stderr.log out out.xml > result_sha256.txt
touch RUN_COMPLETE
echo "Spectra: ${#spectra[@]}"
echo "Rank-local transition files: ${#transitions[@]}"
echo "Run completed successfully."
