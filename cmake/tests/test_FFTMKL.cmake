#First check if we can compile and link with FFT from MKL (interface module available)
try_compile(FLEUR_USE_FFTMKL ${CMAKE_BINARY_DIR} ${CMAKE_SOURCE_DIR}/cmake/tests/test_FFTMKL.f90
LINK_LIBRARIES ${FLEUR_LIBRARIES})

#otherwise compile the interface module from the MKL source. Finding the source is not
#enough, MKL also has to be linked. The sources are compiled in the given order
#(try_compile does not scan module dependencies), the interface has to be first.
if (NOT FLEUR_USE_FFTMKL)
   find_file(MKL_FFT_INTERFACE_FILE mkl_dfti.f90 HINTS "$ENV{MKLROOT}" PATH_SUFFIXES "/include")
   if (MKL_FFT_INTERFACE_FILE)
      try_compile(FLEUR_USE_FFTMKL ${CMAKE_BINARY_DIR}
                  SOURCES ${MKL_FFT_INTERFACE_FILE} ${CMAKE_SOURCE_DIR}/cmake/tests/test_FFTMKL.f90
                  LINK_LIBRARIES ${FLEUR_LIBRARIES})
   endif()
   message("Try to compile MKL-FFT interface at:${MKL_FFT_INTERFACE_FILE}")
endif()


message("FFT from MKL found:${FLEUR_USE_FFTMKL}")
if (CLI_PATCH_INTEL)
   if (FLEUR_USE_FFTMKL)
      message("MKL with FFT found but not used as it is broken with our INTEL_PATCH")
   endif()
   set(FLEUR_USE_FFTMKL false)
endif()

if (FLEUR_USE_FFTMKL)
   set(FLEUR_DEFINITIONS ${FLEUR_DEFINITIONS} "CPP_FFT_MKL")
endif()
