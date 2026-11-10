#!/usr/bin/env bash

# Build the external libraries FLEUR can compile in-place (libxc, Wannier90 and,
# with MPI, ELPA) from a fresh clone of FLEUR and install them into a prefix:
#   <prefix>/lib, <prefix>/include          serial libraries and modules
#   <prefix>/mpi/lib, <prefix>/mpi/include  MPI libraries and modules (self-contained)
# Use them with: configure.sh -libdir <dir>/lib -includedir <dir>/include
#
# Usage: build-fleur-dependencies.sh <compiler> <prefix>
# The source is taken from $FLEUR_REPOSITORY (branch $FLEUR_BRANCH, default develop);
# $FLEUR_BUILD_JOBS limits the parallel compile jobs (default nproc).

set -euo pipefail

compiler=$1
prefix=$2
repository=${FLEUR_REPOSITORY:-https://iffgit.fz-juelich.de/fleur/fleur.git}
branch=${FLEUR_BRANCH:-develop}
jobs=${FLEUR_BUILD_JOBS:-$(nproc)}

work_dir=$(mktemp -d)
trap 'rm -rf "$work_dir"' EXIT
source_dir=$work_dir/fleur

copy_required_file() {
    local build_dir=$1
    local filename=$2
    local destination=$3
    local target_name=${4:-$filename}
    local source_file

    source_file=$(find "$build_dir" -type f -name "$filename" -print -quit)
    if [ -z "$source_file" ]; then
        echo "Required dependency artifact $filename was not built in $build_dir" >&2
        exit 1
    fi

    install -Dm644 "$source_file" "$destination/$target_name"
}

copy_modules() {
    local modules_dir=$1
    local destination=$2

    if ! find "$modules_dir" -maxdepth 1 -type f -name '*.mod' -print -quit | grep -q .; then
        echo "No Fortran modules were built in $modules_dir" >&2
        exit 1
    fi

    install -d "$destination"
    find "$modules_dir" -maxdepth 1 -type f -name '*.mod' -exec cp {} "$destination" \;
}

configure_dependencies() {
    local label=$1
    shift

    cd "$source_dir"
    if [ "$compiler" = "mpi" ] && [ "$label" = "ci-deps-mpi" ]; then
        FLEUR_LIBRARIES="-lscalapack-openmpi" ./configure.sh -l "$label" -c "$compiler" -libxc TRUE -wannier TRUE -cmake_opts -DBUILD_SHARED_LIBS=OFF "$@" >&2
    else
        ./configure.sh -l "$label" -c "$compiler" -libxc TRUE -wannier TRUE -cmake_opts -DBUILD_SHARED_LIBS=OFF "$@" >&2
    fi
    printf '%s\n' "$source_dir/build.$label"
}

if [ "$compiler" = "intel" ] && [ -r /opt/intel/oneapi/setvars.sh ]; then
    # The oneAPI base image requires its compiler environment during docker build.
    # shellcheck disable=SC1091
    set +u
    source /opt/intel/oneapi/setvars.sh --force
    set -u
fi

git clone --depth 1 --branch "$branch" "$repository" "$source_dir"
echo "Building FLEUR dependencies from $repository $branch $(git -C "$source_dir" rev-parse HEAD)"

install -d "$prefix/lib" "$prefix/include"

serial_build=$(configure_dependencies ci-deps -mpi FALSE -hdf5 FALSE)
cmake --build "$serial_build" --parallel "$jobs" --target xcf03 Wannier90_lib
copy_required_file "$serial_build" libxc.a "$prefix/lib"
copy_required_file "$serial_build" libxcf03.a "$prefix/lib"
copy_required_file "$serial_build" libwannier90.a "$prefix/lib"
copy_modules "$serial_build/modules/external" "$prefix/include"
copy_modules "$serial_build/external/wannier90/fortran_mods" "$prefix/include"

if [ "$compiler" = "mpi" ]; then
    install -d "$prefix/mpi/lib" "$prefix/mpi/include"
    mpi_build=$(configure_dependencies ci-deps-mpi -hdf5 FALSE -elpa TRUE)
    cmake --build "$mpi_build" --parallel "$jobs" --target xcf03 Wannier90_lib ELPA
    copy_required_file "$mpi_build" libxc.a "$prefix/mpi/lib"
    copy_required_file "$mpi_build" libxcf03.a "$prefix/mpi/lib"
    copy_required_file "$mpi_build" libwannier90_mpi.a "$prefix/mpi/lib"
    # Installed as libelpa.a, the name configure checks for by default
    copy_required_file "$mpi_build" libelpa_openmp.a "$prefix/mpi/lib" libelpa.a
    copy_modules "$mpi_build/modules/external" "$prefix/mpi/include"
    copy_modules "$mpi_build/external/wannier90/fortran_mods" "$prefix/mpi/include"
    copy_modules "$mpi_build/elpa/src/ELPA/modules" "$prefix/mpi/include"
fi
