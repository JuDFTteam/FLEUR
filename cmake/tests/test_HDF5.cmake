#first try if hdf already works
try_compile(FLEUR_USE_HDF5 ${CMAKE_BINARY_DIR} ${CMAKE_SOURCE_DIR}/cmake/tests/test_HDF5.f90
	    LINK_LIBRARIES ${FLEUR_LIBRARIES}
            )
#now try to find the library by adding the -l stuff to the FLEUR_LIBRARIES
foreach(ADD_String "-lhdf5_fortran;-lhdf5" 
                   "-lhdf5_fortran;-lhdf5_f90cstub;-lhdf5"
		   "-lhdf5_fortran;-lhdf5;-ldl" 
                   "-lhdf5_fortran;-lhdf5_f90cstub;-lhdf5;-ldl"
                   "-lhdf5_fortran;-lhdf5;-lz" 
                   "-lhdf5_fortran;-lhdf5_f90cstub;-lhdf5;-lz"
		   "-lhdf5_fortran;-lhdf5;-ldl;-lz" 
                   "-lhdf5_fortran;-lhdf5_f90cstub;-lhdf5;-ldl;-lz")
   if (NOT FLEUR_USE_HDF5)
     set(TEST_LIBRARIES "${FLEUR_LIBRARIES};${ADD_String}")
     try_compile(FLEUR_USE_HDF5 ${CMAKE_BINARY_DIR} ${CMAKE_SOURCE_DIR}/cmake/tests/test_HDF5.f90
   	    LINK_LIBRARIES ${TEST_LIBRARIES}
            )
     if (FLEUR_USE_HDF5)
          set(FLEUR_LIBRARIES ${TEST_LIBRARIES})
     endif()       	    
   endif()
endforeach()

#check if HDF is parallel
if ( FLEUR_USE_HDF5)
   try_compile(FLEUR_USE_HDF5MPI ${CMAKE_BINARY_DIR} ${CMAKE_SOURCE_DIR}/cmake/tests/test_HDF5MPI.f90
            LINK_LIBRARIES ${FLEUR_LIBRARIES}
            )
endif()

message("HDF5 Library found:${FLEUR_USE_HDF5}")

if (DEFINED CLI_FLEUR_USE_HDF5)
   if (CLI_FLEUR_USE_HDF5)
       if (NOT FLEUR_USE_HDF5)
           message(WARNING "You asked for HDF5 but cmake couldn't find it. We will try to download and compile HDF5 along with FLEUR")
           if(NOT EXISTS "${PROJECT_SOURCE_DIR}/external/hdf5-git/src" )
    	     find_package(Git REQUIRED)
    	     execute_process(COMMAND ${GIT_EXECUTABLE} submodule init external/hdf5-git WORKING_DIRECTORY ${PROJECT_SOURCE_DIR} RESULT_VARIABLE _res_init OUTPUT_QUIET ERROR_QUIET)
    	     execute_process(COMMAND ${GIT_EXECUTABLE} submodule update  WORKING_DIRECTORY ${PROJECT_SOURCE_DIR} RESULT_VARIABLE _res_update OUTPUT_QUIET ERROR_QUIET)
    	     if( ${_res_init} GREATER 0 OR ${_res_update} GREATER 0 )
               message(FATAL_ERROR "HDF5 source could not be downloaded.\n"
                            "We tried: 'git submodule init external/hdf5-git && git submodule update' and resulted in error" )
             endif()
	   endif()
	   set(HDF5_EXTERNALLY_CONFIGURED 1)
	   set(HDF5_EXPORTED_TARGETS "hdf5_fortran-static")
	   set(HDF5_BUILD_FORTRAN ON CACHE BOOL "Build FORTRAN support")
	   set(HDF5_BUILD_CPP_LIB OFF CACHE BOOL "Build HDF5 C++ Library")
	   set(HDF5_BUILD_HL_LIB OFF CACHE BOOL "Build HIGH Level HDF5 Library")
	   set(BUILD_SHARED_LIBS OFF CACHE BOOL "Build Shared Libraries")
	   set(HDF5_BUILD_TOOLS OFF CACHE BOOL "Build HDF5 Tools")
	   set(HDF5_BUILD_EXAMPLES OFF CACHE BOOL "Build HDF5 Library Examples")
	   set(BUILD_TESTING OFF CACHE BOOL "Build HDF5 Unit Testing")
	   if (FLEUR_USE_MPI)
	   set(HDF5_ENABLE_PARALLEL ON CACHE BOOL "Enable parallel build (requires MPI)")
	   else()
		set(HDF5_ENABLE_PARALLEL OFF CACHE BOOL "Enable parallel build (requires MPI)")
	   endif()
	   set(CMAKE_Fortran_MODULE_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/modules/hdf5")
	   add_subdirectory (external/hdf5-git EXCLUDE_FROM_ALL)
	   set(FLEUR_USE_HDF5 TRUE)
           set(FLEUR_USE_HDF5MPI FLEUR_USE_MPI)
           set(FLEUR_COMPILE_HDF true)
           include_directories("${CMAKE_CURRENT_BINARY_DIR}/modules/hdf5/static")
       endif()
   else()
       if (FLEUR_USE_HDF5)
           message("HDF5 library found, but you explicitely asked not to use it")
	   set(FLEUR_USE_HDF5 FALSE)
       endif()
   endif()	   
endif()      

if (FLEUR_USE_HDF5)
   message("Parallel HDF5 Library found:${FLEUR_USE_HDF5MPI}")
   set(FLEUR_DEFINITIONS ${FLEUR_DEFINITIONS} "CPP_HDF") 
   set(FLEUR_MPI_DEFINITIONS ${FLEUR_MPI_DEFINITIONS} "CPP_HDF")
   if (FLEUR_USE_HDF5MPI)
   if (FLEUR_USE_MPI)
      set(FLEUR_MPI_DEFINITIONS ${FLEUR_MPI_DEFINITIONS} "CPP_HDFMPI")
   endif()
   endif()
endif()
