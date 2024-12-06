set(conffile "${CMAKE_BINARY_DIR}/conf.elpa.sh")
file(WRITE ${conffile} "export FC=${CMAKE_Fortran_COMPILER}\n")
file(APPEND ${conffile} "export CC=${CMAKE_C_COMPILER}\n")
if (CMAKE_CXX_COMPILER)
    file(APPEND ${conffile} "export CXX=${CMAKE_CXX_COMPILER}\n")
endif()    
set(LDFLAGS "")
if (FLEUR_LIBRARIES)
    set(LDFLAGS ${FLEUR_LIBRARIES})
    list(JOIN LDFLAGS " " LDFLAGS)
endif()
file(APPEND ${conffile} "export FCFLAGS")
file(APPEND ${conffile} ="${CMAKE_Fortran_FLAGS} ${CMAKE_Fortran_FLAGS_RELEASE} ${LDFLAGS}" "\n")
file(APPEND ${conffile} "export CFLAGS")
file(APPEND ${conffile} ="${CMAKE_C_FLAGS} ${CMAKE_C_FLAGS_RELEASE}" "\n")

set(elpa_flags " --enable-shared=no --enable-c-tests=no --enable-cpp-tests=no --enable-single-precision  --enable-runtime-threading-support-checks --enable-allow-thread-limiting --without-threading-support-check-during-build")
if (FLEUR_USE_OPENMP)
    set(elpa_flags "${elpa_flags} --enable-openmp")
endif()
if (FLEUR_USE_GPU)
    set(elpa_flags "${elpa_flags} --with-nvidia-gpu-support-only")
endif()
file(APPEND ${conffile} "echo 'Building ELPA in: ' $PWD \n")

file(APPEND ${conffile} "./configure ${elpa_flags}")

include(ExternalProject)
ExternalProject_Add(ELPA
PREFIX elpa/
BUILD_IN_SOURCE true
#GIT_REPOSITORY https://gitlab.mpcdf.mpg.de/elpa/elpa.git
#GIT_TAG release_2024_05_001
URL https://elpa.mpcdf.mpg.de/software/tarball-archive/Releases/2024.05.001/elpa-2024.05.001.tar.gz
    
CONFIGURE_COMMAND sh ${conffile}
BUILD_COMMAND make
INSTALL_COMMAND ""
)
#Now make ELPA known to FLEUR
if (FLEUR_USE_OPENMP)
    set(FLEUR_LIBRARIES "-L${CMAKE_BINARY_DIR}/elpa/src/ELPA/.libs;-lelpa_openmp;${FLEUR_LIBRARIES}")
else()
    set(FLEUR_LIBRARIES "-L${CMAKE_BINARY_DIR}/elpa/src/ELPA/.libs;-lelpa;${FLEUR_LIBRARIES}")
endif()    
include_directories("${CMAKE_CURRENT_BINARY_DIR}/elpa/src/ELPA/modules")
set(FLEUR_USE_ELPA TRUE)
set(FLEUR_MPI_DEFINITIONS ${FLEUR_MPI_DEFINITIONS} "CPP_ELPA")
set(FLEUR_USE_INTERNAL_ELPA true)
