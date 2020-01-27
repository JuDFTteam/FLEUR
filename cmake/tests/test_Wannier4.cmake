#first try if Wannier90 modification 4 already works
try_compile(FLEUR_USE_WANN4 ${CMAKE_BINARY_DIR} ${CMAKE_SOURCE_DIR}/cmake/tests/test_Wannier4.f90
	    LINK_LIBRARIES ${FLEUR_LIBRARIES}
            )

message("Wannier90-4 1.2 Library found:${FLEUR_USE_WANN4}")
if (DEFINED CLI_FLEUR_USE_WANNIER)
   if (${CLI_FLEUR_USE_WANNIER})
       if (NOT FLEUR_USE_WANN4)
           message(FATAL_ERROR "You asked for Wannier90 but cmake couldn't find it. Please check your Fortran compiler settings")
       endif()
   else()
       if (FLEUR_USE_WANN4)
           message("Wannier library found, but you explicitely asked not to use it")
	   set(FLEUR_USE_WANN4 FALSE)
       endif()
   endif()	   
endif()  
if (FLEUR_USE_WANN4)
   set(FLEUR_DEFINITIONS ${FLEUR_DEFINITIONS} "CPP_WANN4") 
   set(FLEUR_MPI_DEFINITIONS ${FLEUR_MPI_DEFINITIONS} "CPP_WANN4")
endif()
