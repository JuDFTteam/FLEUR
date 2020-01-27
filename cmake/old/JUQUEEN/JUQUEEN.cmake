#Set the compiler names
set(CMAKE_Fortran_COMPILER mpixlf2008_r)
set(CMAKE_C_COMPILER mpixlc)
#Add include pathes
set(FLEUR_Fortran_FLAGS "-I$ENV{ELPA_INCLUDE} -I$ENV{HDF5_DIR}/include")
#Add linker stuff
set(FLEUR_LIBRARIES "-L$ENV{SCALAPACK_ROOT}/lib;-lelpa;-lscalapack;-L/bgsys/local/lapack/3.3.0_g/lib;-llapack;-L/bgsys/local/lib;-qessl;-lesslsmpbg;-L$ENV{XML2LIB};-lxml2;-L$ENV{HDF5_DIR}/lib;-lhdf5_fortran;-lhdf5;-L/bgsys/local/zlib/lib/;-lz;-L/bgsys/local/szip/lib/;-lsz")
set(FLEUR_DEFINITIONS  "CPP_AIX" )
set(FLEUR_MPI_DEFINITIONS  "CPP_AIX" )

