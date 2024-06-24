#This file contains specific compiler flags for
#individual files and compilers
#E.G. it is used to switch of optimization for some files
#The compiler flags are added at the end and hence can be used
#to overwrite previous settings


if (CMAKE_Fortran_COMPILER_ID MATCHES "GNU")
   set_source_files_properties(${CMAKE_SOURCE_DIR}/src/fleur/vgen/vgen_coulomb.F90 PROPERTIES COMPILE_FLAGS -O0)
endif()

if (CMAKE_Fortran_COMPILER_ID MATCHES "Intel")
   #set_source_files_properties(io/eig66_mpi.F90 PROPERTIES COMPILE_FLAGS -O0)
   #set_source_files_properties(cdn/pwden.F90 PROPERTIES COMPILE_FLAGS -O0)
   #set_source_files_properties(eigen/apws.F90 PROPERTIES COMPILE_FLAGS -O0)
   set_source_files_properties(${CMAKE_SOURCE_DIR}/src/libraries/juDFT/time.F90 PROPERTIES COMPILE_FLAGS -O0)
   set_source_files_properties(${CMAKE_SOURCE_DIR}/src/fleur/io/eig66_mpi.F90 PROPERTIES COMPILE_FLAGS -g0)
   set_source_files_properties(${CMAKE_SOURCE_DIR}/src/fleur/vgen/mt_tofrom_grid.F90 PROPERTIES COMPILE_FLAGS -O1)
   set_source_files_properties(${CMAKE_SOURCE_DIR}/src/fleur/vgen/psqpw.F90 PROPERTIES COMPILE_FLAGS -O1)
   set_source_files_properties(${CMAKE_SOURCE_DIR}/src/fleur/cdn/cdnovlp.F90 PROPERTIES COMPILE_FLAGS -O1)
   
   if (CMAKE_Fortran_COMPILER_VERSION VERSION_LESS "14.1.0.0")
      set_source_files_properties(${CMAKE_SOURCE_DIR}/src/fleur/vgen/vmtxcg.F90 PROPERTIES COMPILE_FLAGS -no-openmp)
   endif()
endif()

if (CMAKE_Fortran_COMPILER_ID MATCHES "PGI")
set_source_files_properties(${CMAKE_SOURCE_DIR}/src/fleur/vgen/mkgylm.f90 PROPERTIES COMPILE_FLAGS "-O0 -Mvect=nosimd")
set_source_files_properties(${CMAKE_SOURCE_DIR}/src/fleur/eigen/hsmt_nonsph.F90 PROPERTIES COMPILE_FLAGS "-O1 -Mvect=nosimd -nomp")
endif()    
