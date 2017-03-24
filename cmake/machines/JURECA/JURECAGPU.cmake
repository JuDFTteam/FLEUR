#Set the compiler names
set(CMAKE_Fortran_COMPILER mpif90)
set(CMAKE_C_COMPILER mpicc)
#Add include pathes
#set(CMAKE_C_FLAGS " -I$ENV{XML2_ROOT}/include")
set(FLEUR_LIBRARIES "-lstdc++;-L$ENV{MKLROOT}/lib/intel64;-lmkl_scalapack_lp64;-lmkl_intel_lp64;-lmkl_pgi_thread;-lmkl_core;-lmkl_blacs_intelmpi_lp64"")
if (ENV{MAGMA_ROOT})
set(FLEUR_Fortran_FLAGS " -I$ENV{MAGMA_ROOT}/include")
set(FLEUR_LIBRARIES "${FLEUR_LIBARIES};-L$ENV{MAGMA_ROOT}/lib;-lmagma")
endif()
set(FLEUR_DEFINITIONS ${FLEUR_DEFINITIONS} "CPP_AIX")
set(FLEUR_MPI_DEFINITIONS ${FLEUR_MPI_DEFINITIONS} "CPP_AIX")
