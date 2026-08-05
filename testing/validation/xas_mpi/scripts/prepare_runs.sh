#!/usr/bin/env bash

set -euo pipefail

SCRIPT_DIR=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)
PACKAGE_ROOT=$(cd -- "${SCRIPT_DIR}/.." && pwd)
FIXTURE_DIR="${PACKAGE_ROOT}/fixtures"
RUNS_DIR="${PACKAGE_ROOT}/runs"

if [[ -e "${RUNS_DIR}" ]]; then
    echo "Refusing to overwrite existing run tree: ${RUNS_DIR}" >&2
    echo "Archive or remove that complete tree before preparing a fresh validation." >&2
    exit 2
fi

(cd "${FIXTURE_DIR}" && sha256sum --check SHA256SUMS)

for label in np1 np2_pure np4_pure np4_shared; do
    run_dir="${RUNS_DIR}/${label}"
    mkdir -p "${run_dir}"
    cp "${FIXTURE_DIR}/cdn.hdf" "${run_dir}/cdn.hdf"
    cp "${FIXTURE_DIR}/inp.xml" "${run_dir}/inp.xml"
    cp "${FIXTURE_DIR}/kpts.xml" "${run_dir}/kpts.xml"
    cp "${FIXTURE_DIR}/sym.xml" "${run_dir}/sym.xml"
    (cd "${run_dir}" && sha256sum cdn.hdf inp.xml kpts.xml sym.xml > SEED_SHA256)
    diff -u "${FIXTURE_DIR}/SHA256SUMS" "${run_dir}/SEED_SHA256"
done

echo "Prepared four isolated run directories under ${RUNS_DIR}."
echo "Each run has an independent copy of the original restart."
