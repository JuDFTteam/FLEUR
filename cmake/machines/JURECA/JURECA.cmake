#Set the compiler names
set(CMAKE_Fortran_COMPILER mpif90)
set(CMAKE_C_COMPILER mpicc)
#Add include pathes
#Add linker stuff
set(FLEUR_Fortran_FLAGS " -mkl")
#set(FLEUR_Fortran_FLAGS "-I$ENV{ELPA_MODULES} -I$ENV{EBROOTHDF5}/include -mkl")
#set(FLEUR_LIBRARIES "-L$ENV{ELPA_LIB};-lelpa_openmp;-lmkl_scalapack_lp64;-lmkl_blacs_intelmpi_lp64;-L$ENV{EBROOTHDF5}/lib;-lhdf5;-lhdf5_fortran")
set(FLEUR_LIBRARIES "-lmkl_scalapack_lp64;-lmkl_blacs_intelmpi_lp64")
