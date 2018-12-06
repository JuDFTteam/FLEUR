#First check if we can compile with ELPA
try_compile(FLEUR_USE_ELPA_ONENODE ${CMAKE_BINARY_DIR} ${CMAKE_SOURCE_DIR}/cmake/tests/test_ELPA.f90
LINK_LIBRARIES ${FLEUR_LIBRARIES})

if (NOT FLEUR_USE_ELPA_ONENODE)
   if (DEFINED CLI_ELPA_OPENMP)
      if (FLEUR_USES_GPU)
      set(TEST_LIBRARIES "-lelpa_onenode;${FLEUR_LIBRARIES}")
      else()
      set(TEST_LIBRARIES "-lelpa_onenode_openmp;${FLEUR_LIBRARIES}")
      endif()
   endif()
   try_compile(FLEUR_USE_ELPA_ONENODE ${CMAKE_BINARY_DIR} ${CMAKE_SOURCE_DIR}/cmake/tests/test_ELPA.f90
LINK_LIBRARIES ${TEST_LIBRARIES})
   if (FLEUR_USE_ELPA_ONENODE)
      set(FLEUR_LIBRARIES "${TEST_LIBRARIES}")
   endif()
endif()

message("ELPA (one node) Library found:${FLEUR_USE_ELPA_ONENODE}")

#Now check for version of elpa
if (FLEUR_USE_ELPA_ONENODE)
    set(FLEUR_USE_ELPA_ONENODE false)
try_compile(FLEUR_USE_ELPA_ONENODE_20180525 ${CMAKE_BINARY_DIR} ${CMAKE_SOURCE_DIR}/cmake/tests/test_ELPA_20180525.f90
LINK_LIBRARIES ${FLEUR_LIBRARIES})
    message("Version check for ELPA:")
    message("20180525  ELPA: ${FLEUR_USE_ELPA_ONENODE_20180525}")
   if (FLEUR_USE_ELPA_ONENODE_20180525)
       set(FLEUR_USE_ELPA_ONENODE TRUE)
       set(FLEUR_DEFINITIONS ${FLEUR_DEFINITIONS} "CPP_ELPA_ONENODE")
   endif()
endif() 
