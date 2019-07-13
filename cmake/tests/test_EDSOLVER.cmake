#first try if EdSolver already works
try_compile(FLEUR_USE_EDSOLVER ${CMAKE_BINARY_DIR} ${CMAKE_SOURCE_DIR}/cmake/tests/test_EDSOLVER.f90
	    LINK_LIBRARIES ${FLEUR_LIBRARIES}) 

message("EDSolver Library found:${FLEUR_USE_EDSOLVER}")

if (DEFINED CLI_FLEUR_USE_EDSOLVER)
   if (${CLI_FLEUR_USE_EDSOLVER})
      if (NOT FLEUR_USE_EDSOLVER)
         message(FATAL ERROR "You explicitly asked for EDSolver but cmake couldn't find it")
      endif()
   else()
      if (FLEUR_USE_EDSOLVER)
         message("EDSolver library found, but you explicitely asked not to use it")
         set(FLEUR_USE_EDSOLVER FALSE)
      endif()
   endif()
endif()

if (FLEUR_USE_EDSOLVER)
   set(FLEUR_DEFINITIONS ${FLEUR_DEFINITIONS} "CPP_EDSOLVER") 
   set(FLEUR_MPI_DEFINITIONS ${FLEUR_MPI_DEFINITIONS} "CPP_EDSOLVER")
endif()
