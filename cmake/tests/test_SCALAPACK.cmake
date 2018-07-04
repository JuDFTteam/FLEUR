#First check if we can compile with LAPACK
try_compile(FLEUR_USE_SCALAPACK ${CMAKE_BINARY_DIR} ${CMAKE_SOURCE_DIR}/cmake/tests/test_SCALAPACK.f90
LINK_LIBRARIES ${FLEUR_LIBRARIES}
            )

#Try typical mkl string
if (NOT FLEUR_USE_SCALAPACK)
     set(TEST_LIBRARIES "${FLEUR_LIBRARIES};-lmkl_scalapack_lp64;-lmkl_blacs_intelmpi_lp64")
     try_compile(FLEUR_USE_SCALAPACK ${CMAKE_BINARY_DIR} ${CMAKE_SOURCE_DIR}/cmake/tests/test_SCALAPACK.f90
           LINK_LIBRARIES ${TEST_LIBRARIES}
            )
     if (FLEUR_USE_SCALAPACK)
          set(FLEUR_LIBRARIES ${TEST_LIBRARIES})
     endif()
endif()

message("SCALAPACK Library found:${FLEUR_USE_SCALAPACK}")
if (FLEUR_USE_SCALAPACK)
   set(FLEUR_MPI_DEFINITIONS ${FLEUR_MPI_DEFINITIONS} "CPP_SCALAPACK")
endif()
