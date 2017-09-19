#Set the compiler names
set(CMAKE_Fortran_COMPILER mpif90)
#set(CMAKE_C_COMPILER mpiicc)
#Add include pathes
#set(FLEUR_Fortran_FLAGS "")
#Add linker stuff
set(FLEUR_LIBRARIES ${FLEUR_LIBRARIES} "/usr/lib/x86_64-linux-gnu/libxml2.so;/usr/lib/x86_64-linux-gnu/libscalapack-openmpi.so;/usr/lib/x86_64-linux-gnu/libblacsF77init-openmpi.so;/usr/lib/x86_64-linux-gnu/libblacs-openmpi.so;/usr/lib/liblapack.so;/usr/lib/libblas.so")
