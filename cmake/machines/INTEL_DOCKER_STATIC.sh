#Set the compilers to mpiifort, mpiicc, mpiicpc
export FC=${FC:=mpiifort}
export CC=${CC:=gcc}
export CXX=${CXX:=mpicxx}

#Special flags for docker image
export FLEUR_LIBRARIES="/usr/local/intel/compilers_and_libraries/linux/mkl/lib/intel64/libmkl_scalapack_lp64.a -Wl,--start-group /usr/local/intel/compilers_and_libraries/linux/mkl/lib/intel64/libmkl_intel_lp64.a /usr/local/intel/compilers_and_libraries/linux/mkl/lib/intel64/libmkl_intel_thread.a /usr/local/intel/compilers_and_libraries/linux/mkl/lib/intel64/libmkl_core.a /usr/local/intel/compilers_and_libraries/linux/mkl/lib/intel64/libmkl_blacs_intelmpi_lp64.a /usr/local/intel/compilers_and_libraries/linux/lib/intel64/libiomp5.a -Wl,--end-group -Wl,-Bstatic -liomp5  -Wl,-Bdynamic  -lpthread -lm -ldl -Wl,-Bstatic"

#export CMAKE_Fortran_FLAGS="-axSSE2 -qopenmp-link static -static_mpi -static-intel"



