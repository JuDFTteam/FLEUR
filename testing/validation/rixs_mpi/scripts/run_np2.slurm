#!/bin/bash
#SBATCH --job-name=rixs-sqay-np2
#SBATCH --nodes=1
#SBATCH --ntasks=2
#SBATCH --cpus-per-task=1
#SBATCH --output=slurm-%x-%j.out
#SBATCH --error=slurm-%x-%j.err

set -euo pipefail

PACKAGE_ROOT=${RIXS_MPI_PACKAGE_ROOT:-${SLURM_SUBMIT_DIR:-$PWD}}
COMMON_SCRIPT="${PACKAGE_ROOT}/scripts/run_common.sh"
RUN_DIR="${PACKAGE_ROOT}/SQAy_np2"
if [[ ! -f "${COMMON_SCRIPT}" || ! -f "${RUN_DIR}/inp.xml" ]]; then
    echo "Cannot find RIXS MPI validation helper or input under ${PACKAGE_ROOT}" >&2
    exit 2
fi
cd "${RUN_DIR}"
exec bash "${COMMON_SCRIPT}" "${RUN_DIR}" 2 rixs_mpi_np2
