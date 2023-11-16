#First check if we can compile with ScaLAPACK
try_compile(FLEUR_USE_SCALAPACK ${CMAKE_BINARY_DIR} ${CMAKE_SOURCE_DIR}/cmake/tests/test_SCALAPACK.f90
LINK_LIBRARIES ${FLEUR_LIBRARIES}
            )

#Try typical mkl string
foreach(test_string "-lmkl_scalapack_lp64;-lmkl_blacs_intelmpi_lp64" "-lscalapack_openmpi" "-lscalapack-openmpi" "-lscalapack" "-Mscalapack")
if (NOT FLEUR_USE_SCALAPACK)
     message("Test for SCALAPACK with:${test_string}")
     set(TEST_LIBRARIES "${test_string};${FLEUR_LIBRARIES}")
     try_compile(FLEUR_USE_SCALAPACK ${CMAKE_BINARY_DIR} ${CMAKE_SOURCE_DIR}/cmake/tests/test_SCALAPACK.f90
           LINK_LIBRARIES ${TEST_LIBRARIES} OUTPUT_VARIABLE compile_output
            )
    if ("$ENV{VERBOSE}")
        message("Scalapack compile test: ${FLEUR_USE_SCALAPACK}\nLINK_LIBRARIES ${TEST_LIBRARIES}\n${compile_output}")
     endif()

     if (FLEUR_USE_SCALAPACK)
          set(FLEUR_LIBRARIES ${TEST_LIBRARIES})
     endif()
endif()
endforeach()



message("SCALAPACK Library found:${FLEUR_USE_SCALAPACK}")
if (FLEUR_USE_SCALAPACK)
   set(FLEUR_MPI_DEFINITIONS ${FLEUR_MPI_DEFINITIONS} "CPP_SCALAPACK")
endif()
