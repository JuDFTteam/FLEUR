#Set the compiler names
set(CMAKE_Fortran_COMPILER mpif90)
#set(CMAKE_C_COMPILER mpiicc)
#Add include pathes
#set(FLEUR_Fortran_FLAGS "")
#Add linker stuff
set(FLEUR_LIBRARIES ${FLEUR_LIBRARIES} "-L/usr/lib;-L/usr/lib/x86_64-linux-gnu;-lxml2;-lscalapack-openmpi;-lblacsF77init-openmpi;-lblacs-openmpi;-llapack;-lblas")
