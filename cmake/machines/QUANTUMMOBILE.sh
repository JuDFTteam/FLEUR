#Set the compilers to mpiifort, mpiicc, mpiicpc
export FC=${FC:=mpif90.openmpi}
export CC=${CC:=mpiicc}
export CXX=${CXX:=mpiicpc}

#Set environment variables usefull for external dependencies, e.g. ELPA
export CFLAGS=-mkl
export CMAKE_Fortran_FLAGS="$CMAKE_Fortran_FLAGS -I/usr/include/hdf5/serial"
export FCFLAGS=-mkl
export LIBS="-lhdf5_serial_fortran -lhdf5_serial -lblacsF77init-openmpi -lblacs-openmpi -lscalapack-openmpi"
export SCALAPACK_LDFLAGS="-lblacsF77init-openmpi -lblacs-openmpi -lscalapack-openmpi"

