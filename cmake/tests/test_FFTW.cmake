#include_directories("/usr/include")
#First check if we can compile with FFT from FFTW
try_compile(FLEUR_USE_FFTW ${CMAKE_BINARY_DIR} ${CMAKE_SOURCE_DIR}/cmake/tests/test_FFTW.f90
LINK_LIBRARIES ${FLEUR_LIBRARIES}  OUTPUT_VARIABLE compile_output
            )

foreach (teststring "-lfftw3" "-lfftw3;-ldl")
    if (NOT FLEUR_USE_FFTW)
      set(TEST_LIBRARIES "${FLEUR_LIBRARIES};${teststring}")
      try_compile(FLEUR_USE_FFTW ${CMAKE_BINARY_DIR} ${CMAKE_SOURCE_DIR}/cmake/tests/test_FFTW.f90
      LINK_LIBRARIES ${TEST_LIBRARIES}  OUTPUT_VARIABLE compile_output
      )
      if ("$ENV{VERBOSE}")
      message("FFTW compile test: ${FLEUR_USE_FFTW}\nLINK_LIBRARIES ${TEST_LIBRARIES}\n${compile_output}")
      endif()
      if (FLEUR_USE_FFTW)
         set(FLEUR_LIBRARIES ${TEST_LIBRARIES})
      endif()
    endif()
  endforeach()

message("FFT from FFTW found:${FLEUR_USE_FFTW}")
if (FLEUR_USE_FFTW)
   set(FLEUR_MPI_DEFINITIONS ${FLEUR_MPI_DEFINITIONS} "CPP_FFTW")
   set(FLEUR_DEFINITIONS ${FLEUR_DEFINITIONS} "CPP_FFTW")
endif()
