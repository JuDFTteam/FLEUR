#!/usr/bin/env bash

set -euo pipefail

# Create the docker images for the CI and optionally push them to the registry.
# Usage: create-docker-images-for-CI.sh [-push] [image ...]
# The dependency stage of each image clones FLEUR itself; the current develop
# revision is passed as build argument so a new commit invalidates the cache.
script_dir=$(cd "$(dirname "$0")" && pwd)
registry=iffregistry.fz-juelich.de/fleur/fleur
fleur_repository=https://iffgit.fz-juelich.de/fleur/fleur.git
fleur_branch=develop

push=false
if [ "${1:-}" = "-push" ]; then
   push=true
   shift
fi
images=("$@")
if [ ${#images[@]} -eq 0 ]; then
   images=(gfortran oneAPI LLVM)
fi

fleur_revision=$(git ls-remote "$fleur_repository" "refs/heads/$fleur_branch" | cut -f1)

for image in "${images[@]}"
do
   podman build -f "$script_dir/$image/Dockerfile" \
      --build-arg FLEUR_BRANCH="$fleur_branch" --build-arg FLEUR_REVISION="$fleur_revision" \
      -t "$registry:$image" "$script_dir"
done

if $push; then
   for image in "${images[@]}"
   do
      podman push "$registry:$image"
   done
fi
