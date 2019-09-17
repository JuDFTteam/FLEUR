#Check if we can compile with MPI
try_compile(FLEUR_USE_MPI ${CMAKE_BINARY_DIR} ${CMAKE_SOURCE_DIR}/cmake/tests/test_MPI.f90
LINK_LIBRARIES ${FLEUR_LIBRARIES}
            )

message("MPI Library found:${FLEUR_USE_MPI}")

if (DEFINED CLI_FLEUR_USE_MPI)
   if (${CLI_FLEUR_USE_MPI})
       if (NOT FLEUR_USE_MPI)
           message(FATAL_ERROR "You asked for MPI but cmake couldn't find it. Please check your Fortran compiler settings")
       endif()
   else()
       if (FLEUR_USE_MPI)
           message("MPI library found, but you explicitely asked not to use it")
	   set(FLEUR_USE_MPI FALSE)
       endif()
   endif()	   
endif()


If (FLEUR_USE_MPI)
   set(FLEUR_MPI_DEFINITIONS ${FLEUR_MPI_DEFINITIONS} "CPP_MPI")
   set(FLEUR_USE_SERIAL FALSE)
else()
   set(FLEUR_USE_SERIAL TRUE)
endif()

if (DEFINED CLI_FLEUR_USE_SERIAL)
  set(FLEUR_USE_SERIAL CLI_FLEUR_USE_SERIAL)
endif()
