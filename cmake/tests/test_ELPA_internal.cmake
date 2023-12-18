set(conffile "${CMAKE_BINARY_DIR}/conf.elpa.sh")
file(WRITE ${conffile} "./autogen.sh\n")
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

set(elpa_flags " --enable-shared=no --enable-c-tests=no --enable-cpp-tests=no")
if (FLEUR_USE_OPENMP)
    set(elpa_flags "${elpa_flags} --enable-openmp")
endif()
if (FLEUR_USE_GPU)
    set(elpa_flags "${elpa_flags} --enable-nvidia-gpu")
endif()

file(APPEND ${conffile} "./configure ${elpa_flags}  --prefix=${CMAKE_BINARY_DIR}/elpa")

include(ExternalProject)
ExternalProject_Add(ELPA
SOURCE_DIR elpa/
BINARY_DIR elpa/
GIT_REPOSITORY https://gitlab.mpcdf.mpg.de/elpa/elpa.git
GIT_TAG new_release_2023.05.001.rc1
CONFIGURE_COMMAND sh ${conffile}
BUILD_COMMAND make
)
#Now make ELPA known to FLEUR
if (FLEUR_USE_OPENMP)
    set(FLEUR_LIBRARIES "-L${CMAKE_BINARY_DIR}/elpa/lib;-lelpa_openmp;${FLEUR_LIBRARIES}")
else()
    set(FLEUR_LIBRARIES "-L${CMAKE_BINARY_DIR}/elpa/lib;-lelpa;${FLEUR_LIBRARIES}")
endif()    
include_directories("${CMAKE_CURRENT_BINARY_DIR}/elpa/include/elpa-2023.05.001.rc1/modules")
set(FLEUR_USE_ELPA TRUE)
set(FLEUR_MPI_DEFINITIONS ${FLEUR_MPI_DEFINITIONS} "CPP_ELPA" "CPP_ELPA2" "CPP_ELPA_201705003")
set(FLEUR_USE_INTERNAL_ELPA true)
