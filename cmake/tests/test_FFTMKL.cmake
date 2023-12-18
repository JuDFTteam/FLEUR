#First check if we can compile with FFT from MKL
try_compile(FLEUR_USE_FFTMKL ${CMAKE_BINARY_DIR} ${CMAKE_SOURCE_DIR}/cmake/tests/test_FFTMKL.f90
LINK_LIBRARIES ${FLEUR_LIBRARIES})

#compile interface module
if (NOT FLEUR_USE_FFTMKL)
	find_file(MKL_FFT_INTERFACE_FILE mkl_dfti.f90 HINTS "$ENV{MKLROOT}" PATH_SUFFIXES "/include")
   if (MKL_FFT_INTERFACE_FILE)
      #add target for interface
      #add_library(MKL_FFT_INTERFACE OBJECT MKL_FFT_INTERFACE_FILE)
      set(FLEUR_USE_FFTMKL True)
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
   set(FLEUR_MPI_DEFINITIONS ${FLEUR_MPI_DEFINITIONS} "CPP_FFT_MKL")
   set(FLEUR_DEFINITIONS ${FLEUR_DEFINITIONS} "CPP_FFT_MKL")
endif()
