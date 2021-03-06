#First check if we can compile with FFT from MKL
try_compile(FLEUR_USE_FFTMKL ${CMAKE_BINARY_DIR} ${CMAKE_SOURCE_DIR}/cmake/tests/test_FFTMKL.f90
LINK_LIBRARIES ${FLEUR_LIBRARIES}
            )
message("FFT from MKL found:${FLEUR_USE_FFTMKL}")
if (CLI_PATCH_INTEL)
   if (FLEUR_USE_FFTMKL)
      message("MKL with FFT found but not used as it is broken with our INTEL_PATCH")
   endif()
   set(FLEUR_USE_FFTMKL false)
endif()

if (FLEUR_USE_FFTMKL)
   set(FLEUR_MPI_DEFINITIONS ${FLEUR_MPI_DEFINITIONS} "CPP_FFT_MKL")
   set(FLEUR_DEFINITIONS ${FLEUR_DEFINITIONS} "CPP_FFT_MKL")
endif()
