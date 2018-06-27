#first try if hdf already works
try_compile(FLEUR_USE_HDF5 ${CMAKE_BINARY_DIR} ${CMAKE_SOURCE_DIR}/cmake/tests/test_HDF5.f90
	    LINK_LIBRARIES ${FLEUR_LIBRARIES}
            )
#now try to find the library by adding the -l stuff to the FLEUR_LIBRARIES
if (NOT FLEUR_USE_HDF5)
     set(TEST_LIBRARIES "${FLEUR_LIBRARIES};-lhdf5_fortran;-lhdf5")
     try_compile(FLEUR_USE_HDF5 ${CMAKE_BINARY_DIR} ${CMAKE_SOURCE_DIR}/cmake/tests/test_HDF5.f90
   	    LINK_LIBRARIES ${TEST_LIBRARIES}
            )
     if (FLEUR_USE_HDF5)
          set(FLEUR_LIBRARIES ${TEST_LIBRARIES})
     endif()       	    
endif()

#now try to find the library by adding the -l stuff to the FLEUR_LIBRARIES
if (NOT FLEUR_USE_HDF5)
     set(TEST_LIBRARIES "${FLEUR_LIBRARIES};-lhdf5_fortran;-lhdf5_f90cstub;-lhdf5")
     try_compile(FLEUR_USE_HDF5 ${CMAKE_BINARY_DIR} ${CMAKE_SOURCE_DIR}/cmake/tests/test_HDF5.f90
   	    LINK_LIBRARIES ${TEST_LIBRARIES}
            )
     if (FLEUR_USE_HDF5)
          set(FLEUR_LIBRARIES ${TEST_LIBRARIES})
     endif()       	    
endif()

#now try the find_package feature
if (NOT FLEUR_USE_HDF5)
      find_package(HDF5)
      if (NOT HDF5_LIBRARIES MATCHES "NOTFOUND")
          set(TEST_LIBRARIES ${HDF5_Fortran_LIBRARIES} ${FLEUR_LIBRARIES})
	  set(STORE_FLAGS ${CMAKE_Fortran_FLAGS})
          set(CMAKE_Fortran_FLAGS "-I${HDF5_INCLUDE_LIBRARIES}" ${CMAKE_Fortran_FLAGS})

try_compile(FLEUR_USE_HDF5 ${CMAKE_BINARY_DIR} ${CMAKE_SOURCE_DIR}/cmake/tests/test_HDF5.f90
            LINK_LIBRARIES ${TEST_LIBRARIES}
            )
	    if (${FLEUR_USE_HDF5})
	       set(FLEUR_LIBRARIES ${HDF5_Fortran_LIBRARIES} ${FLEUR_LIBRARIES})
	       set(FLEUR_MPI_LIBRARIES ${HDF5_Fortran_LIBRARIES} ${FLEUR_MPI_LIBRARIES})
	    else()
               set(CMAKE_Fortran_FLAGS ${STORE_FLAGS})
	    endif()   
      endif()	
endif()       

#check if HDF is parallel
if ( FLEUR_USE_HDF5)
   try_compile(FLEUR_USE_HDF5MPI ${CMAKE_BINARY_DIR} ${CMAKE_SOURCE_DIR}/cmake/tests/test_HDF5MPI.f90
            LINK_LIBRARIES ${FLEUR_LIBRARIES}
            )
endif()

message("HDF5 Library found:${FLEUR_USE_HDF5}")

if (DEFINED CLI_FLEUR_USE_HDF5})
   if (CLI_FLEUR_USE_HDF5})
       if (NOT FLEUR_USE_HDF5)
           message(FATAL_ERROR "You asked for HDF5 but cmake couldn't find it. Please set HDF5_ROOT and or give additional compiler/linker flags")
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
