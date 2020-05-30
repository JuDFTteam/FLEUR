#include_directories("/usr/include")
#First check if we can compile with FFT from FFTW
try_compile(FLEUR_USE_FFTW ${CMAKE_BINARY_DIR} ${CMAKE_SOURCE_DIR}/cmake/tests/test_FFTW.f90
LINK_LIBRARIES ${FLEUR_LIBRARIES} 
            )

message("FFT from FFTW found:${FLEUR_USE_FFTW}")
if (FLEUR_USE_FFTW)
   set(FLEUR_MPI_DEFINITIONS ${FLEUR_MPI_DEFINITIONS} "CPP_FFTW")
   set(FLEUR_DEFINITIONS ${FLEUR_DEFINITIONS} "CPP_FFTW")
endif()
