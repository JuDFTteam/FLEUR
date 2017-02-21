#Set the compiler names
set(CMAKE_Fortran_COMPILER mpiifort)
set(CMAKE_C_COMPILER mpiicc)
#Add include pathes
#set(FLEUR_Fortran_FLAGS "")
#Add linker stuff
set(FLEUR_LIBRARIES ${FLEUR_LIBRARIES} "-lmkl_scalapack_lp64;-lmkl_blacs_intelmpi_lp64")
