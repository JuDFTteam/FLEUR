#First check if we can compile with ELPA
try_compile(FLEUR_USE_ELPA ${CMAKE_BINARY_DIR} ${CMAKE_SOURCE_DIR}/cmake/tests/test_ELPA.f90
LINK_LIBRARIES ${FLEUR_LIBRARIES}
            )

message("ELPA Library found:${FLEUR_USE_ELPA}")

if (FLEUR_USE_ELPA)
   set(FLEUR_MPI_DEFINITIONS ${FLEUR_MPI_DEFINITIONS} "CPP_ELPA" "CPP_ELPA_201605003")
endif() 
