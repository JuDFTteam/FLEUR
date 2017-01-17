#first try if hdf already works
try_compile(FLEUR_USE_HDF5 ${CMAKE_BINARY_DIR} ${CMAKE_SOURCE_DIR}/cmake/tests/test_HDF5.f90
	    LINK_LIBRARIES ${FLEUR_LIBRARIES}
            )
#now try to find the library using HDF5_ROOT variable
if (NOT FLEUR_USE_HDF5)
   if (DEFINED ENV{HDF5_ROOT})
     find_path(HDF5_INCLUDE hdf5.mod PATHS $ENV{HDF5_ROOT} PATH_SUFFIXES include NO_DEFAULT_PATH)
     find_path(HDF5_LIB libhdf5_fortran.a PATHS $ENV{HDF5_ROOT} PATH_SUFFIXES lib NO_DEFAULT_PATH)
     if (HDF5_LIB) 
        set(TEST_LIBRARIES "-L${HDF5_LIB};-lhdf5_fortran;-lhdf5;${FLEUR_LIBRARIES}")
     endif()
     if (HDF5_INCLUDE)
        set(STORE_FLAGS ${CMAKE_Fortran_FLAGS})
        set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -I${HDF5_INCLUDE}")
     endif()
	try_compile(FLEUR_USE_HDF5 ${CMAKE_BINARY_DIR} ${CMAKE_SOURCE_DIR}/cmake/tests/test_HDF5.f90
   	    LINK_LIBRARIES ${TEST_LIBRARIES}
            )
       if (NOT FLEUR_USE_HDF5)
	       set(CMAKE_Fortran_FLAGS ${STORE_FLAGS})
       else()
               set(FLEUR_LIBRARIES "-L${HDF5_LIB};-lhdf5_fortran;-lhdf5;${FLEUR_LIBRARIES}")
               set(FLEUR_MPI_LIBRARIES "-L${HDF5_LIB};-lhdf5_fortran;-lhdf5;${FLEUR_MPI_LIBRARIES}")
       endif()       	    
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
if (FLEUR_USE_HDF5)
   set(FLEUR_DEFINITIONS ${FLEUR_DEFINITIONS} "CPP_HDF") 
   set(FLEUR_MPI_DEFINITIONS ${FLEUR_MPI_DEFINITIONS} "CPP_HDF")
   if (FLEUR_USE_HDF5MPI)
   if (FLEUR_USE_MPI)
      set(FLEUR_MPI_DEFINITIONS ${FLEUR_MPI_DEFINITIONS} "CPP_HDFMPI")
   endif()
   endif()
endif()
