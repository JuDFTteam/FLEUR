echo "Using special config for IFF-cluster"

#Set the compilers to mpiifort, mpiicc, mpiicpc
export FC=mpiifort
export CC=mpiicc
export CXX=mpiicpc

FLEUR_LIBRARIES="-lxml2;-lmkl_scalapack_lp64;-lmkl_blacs_intelmpi_lp64;-mt_mpi"
echo "Wannier is switched off-manually"
CLI_USE_WANNIER=FALSE
echo "ChASE is switched off-manually"
CLI_USE_CHASE=FALSE
echo "Broken linker can be used"
CLI_WARN_ONLY=1

#Set environment variables usefull for external dependencies, e.g. ELPA
export CFLAGS=-mkl
export CMAKE_Fortran_FLAGS="$CMAKE_Fortran_FLAGS -mkl"
export FCFLAGS=-mkl
export LIBS="-mkl -lmkl_scalapack_lp64 -lmkl_blacs_intelmpi_lp64"
export SCALAPACK_LDFLAGS="-lmkl_scalapack_lp64 -lmkl_blacs_intelmpi_lp64"

