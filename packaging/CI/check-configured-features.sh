#!/usr/bin/env bash

set -euo pipefail

build_dir=$1
require_elpa=${2:-false}
configure_output="$build_dir/configure.out"

if [ ! -r "$configure_output" ]; then
    echo "Missing configure summary: $configure_output" >&2
    exit 1
fi

require_feature() {
    local feature=$1

    if ! grep -Eq "$feature.*TRUE" "$configure_output"; then
        echo "Required CI feature was not enabled: $feature" >&2
        exit 1
    fi
}

require_feature "LibXC Library"
require_feature "Wannier90 module API"

if [ "$require_elpa" = true ]; then
    require_feature "ELPA \\(distributed\\) Library"
fi
