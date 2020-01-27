#Set the compilers to mpiifort, mpiicc, mpiicpc
export FC=${FC:=mpiifort}
export CC=${CC:=mpiicc}
export CXX=${CXX:=mpiicpc}

#Set environment variables usefull for external dependencies, e.g. ELPA
export CMAKE_Fortran_FLAGS="$CMAKE_Fortran_FLAGS -mkl"
export FCFLAGS=-mkl
export LIBS="-mkl -lmkl_scalapack_lp64 -lmkl_blacs_intelmpi_lp64"
export SCALAPACK_LDFLAGS="-lmkl_scalapack_lp64 -lmkl_blacs_intelmpi_lp64"

