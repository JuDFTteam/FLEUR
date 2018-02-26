#first try if Wannier90 already works
try_compile(FLEUR_USE_WANN ${CMAKE_BINARY_DIR} ${CMAKE_SOURCE_DIR}/cmake/tests/test_Wannier90.f90
	    LINK_LIBRARIES ${FLEUR_LIBRARIES}
            )

message("Wannier90 1.2 Library found:${FLEUR_USE_WANN}")

if (DEFINED ENV{FLEUR_USE_WANNIER})
   if (ENV{FLEUR_USE_WANNIER})
       if (NOT FLEUR_USE_WANN)
           message(FATAL_ERROR "You asked for Wannier90 but cmake couldn't find it. Please check your Fortran compiler settings")
       endif()
   else()
       if (FLEUR_USE_WANN)
           message("Wannier library found, but you explicitely asked not to use it")
	   set(FLEUR_USE_WANN FALSE)
       endif()
   endif()	   
endif()  

if (FLEUR_USE_WANN)
   set(FLEUR_DEFINITIONS ${FLEUR_DEFINITIONS} "CPP_WANN") 
   set(FLEUR_MPI_DEFINITIONS ${FLEUR_MPI_DEFINITIONS} "CPP_WANN")
endif()
