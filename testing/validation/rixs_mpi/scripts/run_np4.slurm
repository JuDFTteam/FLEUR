#!/bin/bash
#SBATCH --job-name=rixs-sqay-np4
#SBATCH --nodes=1
#SBATCH --ntasks=4
#SBATCH --cpus-per-task=1
#SBATCH --output=slurm-%x-%j.out
#SBATCH --error=slurm-%x-%j.err

set -euo pipefail

PACKAGE_ROOT=${RIXS_MPI_PACKAGE_ROOT:-${SLURM_SUBMIT_DIR:-$PWD}}
COMMON_SCRIPT="${PACKAGE_ROOT}/scripts/run_common.sh"
RUN_DIR="${PACKAGE_ROOT}/SQAy_np4"
if [[ ! -f "${COMMON_SCRIPT}" || ! -f "${RUN_DIR}/inp.xml" ]]; then
    echo "Cannot find RIXS MPI validation helper or input under ${PACKAGE_ROOT}" >&2
    exit 2
fi
cd "${RUN_DIR}"
exec bash "${COMMON_SCRIPT}" "${RUN_DIR}" 4 rixs_mpi_np4
