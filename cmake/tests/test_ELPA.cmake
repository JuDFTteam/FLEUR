#First check if we can compile with ELPA
try_compile(FLEUR_USE_ELPA ${CMAKE_BINARY_DIR} ${CMAKE_SOURCE_DIR}/cmake/tests/test_ELPA.f90
LINK_LIBRARIES ${FLEUR_LIBRARIES})

if (NOT FLEUR_USE_ELPA)
   if (DEFINED ENV{ELPA_MODULES})
      set(STORE_FLAGS ${CMAKE_Fortran_FLAGS})
      set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -I$ENV{ELPA_MODULES}")
   endif()
   if (DEFINED ENV{ELPA_LIB})
      set(TEST_LIBRARIES "-L$ENV{ELPA_LIB};-lelpa_openmp;-lstdc++;${FLEUR_LIBRARIES}")
   endif()
   try_compile(FLEUR_USE_ELPA ${CMAKE_BINARY_DIR} ${CMAKE_SOURCE_DIR}/cmake/tests/test_ELPA.f90
LINK_LIBRARIES ${TEST_LIBRARIES})
   if (FLEUR_USE_ELPA)
      set(FLEUR_LIBRARIES "${TEST_LIBRARIES}")
   else()
      set(CMAKE_Fortran_FLAGS ${STORE_FLAGS})
   endif()
endif()


message("ELPA Library found:${FLEUR_USE_ELPA}")

#Now check for version of elpa
if (FLEUR_USE_ELPA)
    set(FLEUR_USE_ELPA false)
    try_compile(FLEUR_USE_ELPA_OLD ${CMAKE_BINARY_DIR} ${CMAKE_SOURCE_DIR}/cmake/tests/test_ELPA_OLD.f90
LINK_LIBRARIES ${FLEUR_LIBRARIES})
try_compile(FLEUR_USE_ELPA_NEW ${CMAKE_BINARY_DIR} ${CMAKE_SOURCE_DIR}/cmake/tests/test_ELPA_NEW.f90
LINK_LIBRARIES ${FLEUR_LIBRARIES})
try_compile(FLEUR_USE_ELPA_201605003 ${CMAKE_BINARY_DIR} ${CMAKE_SOURCE_DIR}/cmake/tests/test_ELPA_201605003.f90
LINK_LIBRARIES ${FLEUR_LIBRARIES})
    message("Version check for ELPA:")
    message("OLD ELPA      : ${FLEUR_USE_ELPA_OLD}")
    message("NEW ELPA      : ${FLEUR_USE_ELPA_NEW}")
    message("201605003 ELPA: ${FLEUR_USE_ELPA_201605003}")
#Set preprocessor switches
   if (FLEUR_USE_ELPA_OLD)
       set(FLEUR_USE_ELPA true)
       set(FLEUR_MPI_DEFINITIONS ${FLEUR_MPI_DEFINITIONS} "CPP_ELPA" "CPP_ELPA2")
   endif()
   if (FLEUR_USE_ELPA_NEW)
       set(FLEUR_USE_ELPA true)
       set(FLEUR_MPI_DEFINITIONS ${FLEUR_MPI_DEFINITIONS} "CPP_ELPA" "CPP_ELPA2" "CPP_ELPA_NEW")
   endif()
   if (FLEUR_USE_ELPA_201605003)
       set(FLEUR_USE_ELPA true)
       set(FLEUR_MPI_DEFINITIONS ${FLEUR_MPI_DEFINITIONS} "CPP_ELPA" "CPP_ELPA2" "CPP_ELPA_201605003")
   endif()
endif() 
