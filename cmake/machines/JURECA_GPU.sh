#Configuration for JURECA@FZJ; use NVHPC toolchain
#-------------------------------------------------------------------------------
# Copyright (c) 2023 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
# This file is part of FLEUR and available as free software under the conditions
# of the MIT license as expressed in the LICENSE file in more detail.
#--------------------------------------------------------------------------------
ml NVHPC OpenMPI CMake HDF5
echo "PGI toolchain used"
export FC=mpif90
export CC=mpicc
export CXX=mpic++ 
export FLEUR_LIBRARIES="${NVHPC}/Linux_x86_64/23.1/comm_libs/openmpi4/openmpi-4.0.5/lib/libscalapack.a"
