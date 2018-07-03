echo "Using special config for IFF-cluster"
echo "Wannier is switched off-manually"

#Set the compilers to mpiifort, mpiicc, mpiicpc
export FC=mpiifort
export CC=mpiicc
export CXX=mpiicpc

FLEUR_LIBRARIES="-lxml2;-lscalapack_lp64;-lmkl_blacs_intelmpi_lp64;-mt_mpi"
CLI_USE_WANNIER=FALSE
CLI_USE_CHASE=FALSE

#Set environment variables usefull for external dependencies, e.g. ELPA
export CFLAGS=-mkl
export CMAKE_Fortran_FLAGS="$CMAKE_Fortran_FLAGS -mkl"
export FCFLAGS=-mkl
export LIBS="-mkl -lmkl_scalapack_lp64 -lmkl_blacs_intelmpi_lp64"
export SCALAPACK_LDFLAGS="-lmkl_scalapack_lp64 -lmkl_blacs_intelmpi_lp64"

