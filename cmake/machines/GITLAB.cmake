#Set the compiler names
set(CMAKE_Fortran_COMPILER mpif90)
#set(CMAKE_C_COMPILER mpiicc)
#Add include pathes
set(FLEUR_Fortran_FLAGS "-I/opt/include")
#Add linker stuff
set(FLEUR_LIBRARIES ${FLEUR_LIBRARIES} "-L/opt/lib;-lxcf03;-lxc;-L/usr/lib;-L/usr/lib/x86_64-linux-gnu;-lxml2;-lscalapack-openmpi;-lblacsF77init-openmpi;-lblacs-openmpi;-llapack;-lblas")
