##    This file can help you fixing compilation problems
##    You should modify it and put it into the working-directory 
##    in which you call the configuration scipt 
##    Note: the full file is a comment currently, lines starting with 
##    single # should be modified. Remove the # afterwards.

##    Set the compiler names, often cmake is good at finding the C-compiler
##    but not the fortran compiler you what to use
#set(CMAKE_Fortran_COMPILER mpiifort)
#set(CMAKE_C_COMPILER mpiicc)

##    Set options for the FORTRAN compiler:
##    You at least will need something like -r8 to promote real variables to double
##    precision. You can check cmake/compilerflags to see what is used for known compilers
##    Add also include pathes here. This might be nescessary for libraries with a F90 interface
##    such as ELPA,HDF5,...
#set(FLEUR_Fortran_FLAGS "-r8 -Isomepath")
##    Add linker stuff. Here you should add the -L and -l options needed for the linker to
##    find libraries. Please mind the format!
#set(FLEUR_LIBRARIES ${FLEUR_LIBRARIES} "-L$ENV{HOME/somepath;-lmkl_scalapack_lp64;-lmkl_blacs_intelmpi_lp64")
