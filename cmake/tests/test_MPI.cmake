#Check if we can compile with MPI
try_compile(FLEUR_USE_MPI ${CMAKE_BINARY_DIR} ${CMAKE_SOURCE_DIR}/cmake/tests/test_MPI.f90
LINK_LIBRARIES ${FLEUR_LIBRARIES}
            )

message("MPI Library found:${FLEUR_USE_MPI}")

if (FLEUR_USE_MPI)
   set(FLEUR_MPI_DEFINITIONS ${FLEUR_MPI_DEFINITIONS} "CPP_MPI")
endif()