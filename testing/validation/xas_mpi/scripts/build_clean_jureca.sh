#!/usr/bin/env bash

set -euo pipefail

: "${FLEUR_SOURCE_DIR:?Set FLEUR_SOURCE_DIR to the clean FLEUR checkout.}"
EXPECTED_COMMIT=b3818dad42aa56fa3a835a9cf8e972e23ea2dd87
BUILD_LABEL=xas_mpi_b3818dad
BUILD_DIR="${FLEUR_SOURCE_DIR}/build.${BUILD_LABEL}"

source_commit=$(git -C "${FLEUR_SOURCE_DIR}" rev-parse HEAD)
source_status=$(git -C "${FLEUR_SOURCE_DIR}" status --short)
if [[ "${source_commit}" != "${EXPECTED_COMMIT}" ]]; then
    echo "Expected commit ${EXPECTED_COMMIT}, found ${source_commit}." >&2
    exit 2
fi
if [[ -n "${source_status}" ]]; then
    echo "The FLEUR source checkout is not clean:" >&2
    printf '%s\n' "${source_status}" >&2
    exit 2
fi
if [[ -e "${BUILD_DIR}" ]]; then
    echo "Refusing to reuse existing build directory: ${BUILD_DIR}" >&2
    exit 2
fi

export FC=mpif90
export CC=mpicc
export CXX=mpicxx

cd "${FLEUR_SOURCE_DIR}"
./configure.sh -c mpi -l "${BUILD_LABEL}" -mpi true -scalapack true \
    -link '-lmkl_scalapack_lp64;-lmkl_blacs_intelmpi_lp64'
cmake --build "${BUILD_DIR}" --target fleur_MPI -j 8

FLEUR_BIN="${BUILD_DIR}/fleur_MPI"
"${FLEUR_BIN}" -version
ldd "${FLEUR_BIN}" | grep -iE 'mpi|blacs|scalapack|hdf5|xml'
sha256sum "${FLEUR_BIN}"

echo "Build complete. Export this path before submission:"
echo "export FLEUR_BIN=${FLEUR_BIN}"
