#Set the compiler names
set(CMAKE_Fortran_COMPILER mpif90)
#set(CMAKE_C_COMPILER mpiicc)
#Add include pathes
set(FLEUR_Fortran_FLAGS "")
#Add linker stuff
set(FLEUR_LIBRARIES "-lmkl_scalapack_lp64;-lmkl_blacs_intelmpi_lp64")