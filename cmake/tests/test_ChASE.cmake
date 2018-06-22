#First check if we can compile with ChASE
try_compile(FLEUR_USE_CHASE ${CMAKE_BINARY_DIR} ${CMAKE_SOURCE_DIR}/cmake/tests/test_ChASE.f90 LINK_LIBRARIES ${FLEUR_LIBRARIES})

#if (NOT FLEUR_USE_CHASE)
#   find_package(chase)
#   find_package(chase-fleur)
#
#   set(TEST_LIBRARIES ${FLEUR_LIBRARIES} ${CHASE_LIBRARIES} ${CHASE-FLEUR_LIBRARIES})
#   try_compile(FLEUR_USE_CHASE ${CMAKE_BINARY_DIR} ${CMAKE_SOURCE_DIR}/cmake/tests/test_ChASE.f90 LINK_LIBRARIES ${TEST_LIBRARIES})
#
#   if (FLEUR_USE_CHASE)
#      set(FLEUR_MPI_LIBRARIES ${FLEUR_MPI_LIBRARIES} ${CHASE_LIBRARIES} ${CHASE-FLEUR_LIBRARIES})
#   endif()
#endif()

message("ChASE Library found:${FLEUR_USE_CHASE}")
if (DEFINED CLI_FLEUR_USE_CHASE)
   if (${CLI_FLEUR_USE_CHASE})
      if (NOT FLEUR_USE_CHASE)
         set(FLEUR_USE_CHASE TRUE)
	 message("Test for Chase failed, but you specified to use it anyway...")
      endif()
   else()
      if (FLEUR_USE_CHASE)
         set(FLEUR_USE_CHASE FALSE)
         message("Test for Chase succeeded, but you specified not to use it.")
      endif()
   endif()
endif()

if (FLEUR_USE_CHASE)
#   set(FLEUR_DEFINITIONS ${FLEUR_DEFINITIONS} "CPP_CHASE") 
   set(FLEUR_MPI_DEFINITIONS ${FLEUR_MPI_DEFINITIONS} "CPP_CHASE")
endif()

