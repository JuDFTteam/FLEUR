#!/usr/bin/env bash

# This is the consistent JURECA environment used by the successful August 2026
# validation. Review the versions against `module spider` before use if JURECA
# has moved to a newer stage.
module purge
module load Stages/2026
module load Intel/2025.2.0
module load imkl/2025.2.0
module load ParaStationMPI/5.13.0-1
module load libxml2/2.14.3
module load HDF5/1.14.6
module load CMake/4.0.3
module load git/2.50.1

module -t list
