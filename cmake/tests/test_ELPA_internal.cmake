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
if (FLEUR_COMPILE_SCALAPACK)
    set(LDFLAGS "${LDFLAGS} -L${CMAKE_BINARY_DIR}/external/SCALAPACK-git/lib -lscalapack")
endif()   

file(APPEND ${conffile} "export FCFLAGS")
file(APPEND ${conffile} ="${CMAKE_Fortran_FLAGS} ${CMAKE_Fortran_FLAGS_RELEASE} ${LDFLAGS}" "\n")
file(APPEND ${conffile} "export CFLAGS")
file(APPEND ${conffile} ="${CMAKE_C_FLAGS} ${CMAKE_C_FLAGS_RELEASE} ${LDFLAGS}" "\n")

set(elpa_flags "--disable-c-tests --disable-cpp-tests --enable-shared=no --enable-c-tests=no --enable-cpp-tests=no --enable-single-precision --disable-avx512-kernels ")
set(elpa_flags "${elpa_flags} --enable-runtime-threading-support-checks --enable-allow-thread-limiting --without-threading-support-check-during-build")

if (DEFINED ENV{ELPA_CONF})
    set(elpa_flags "${elpa_flags} $ENV{ELPA_CONF}")
endif()
if (FLEUR_USE_OPENMP)
    set(elpa_flags "${elpa_flags} --enable-openmp")
endif()
if (FLEUR_USE_GPU)
    #check the cc string
     #Check if a CC is given
     if ("${FLEUR_CC_MODE}" GREATER_EQUAL 80)
        set(elpa_flags "${elpa_flags} --disable-sse-assembly --disable-avx --disable-avx2 --disable-sse --enable-nvidia-gpu-kernels --enable-nvidia-sm80-gpu-kernels --with-NVIDIA-sm_80-support-only --enable-gpu-streams=nvidia --with-NVIDIA-GPU-compute-capability=sm_80")
    else()
        set(elpa_flags "${elpa_flags} --disable-sse-assembly --disable-avx --disable-avx2 --disable-sse --enable-nvidia-gpu-kernels --enable-gpu-streams=nvidia --with-NVIDIA-GPU-compute-capability=sm_${FLEUR_CC_MODE}")
    endif()        
endif()
file(APPEND ${conffile} "echo 'Building ELPA in: ' $PWD \n")

file(APPEND ${conffile} "./configure ${elpa_flags} \n")

include(ExternalProject)
if (FLEUR_COMPILE_SCALAPACK)
    ExternalProject_Add(ELPA
    PREFIX elpa/
    BUILD_IN_SOURCE true
    #GIT_REPOSITORY https://gitlab.mpcdf.mpg.de/elpa/elpa.git
    #GIT_TAG release_2024_05_001
    URL https://elpa.mpcdf.mpg.de/software/tarball-archive/Releases/2024.05.001/elpa-2024.05.001.tar.gz
    PATCH_COMMAND patch -p1 < ${CMAKE_SOURCE_DIR}/cmake/patch/elpa-make-transform-public.patch
    CONFIGURE_COMMAND sh ${conffile}
    BUILD_COMMAND make
    INSTALL_COMMAND ""
    DEPENDS scalapack
    )
else()
    ExternalProject_Add(ELPA
    PREFIX elpa/
    BUILD_IN_SOURCE true
    #GIT_REPOSITORY https://gitlab.mpcdf.mpg.de/elpa/elpa.git
    #GIT_TAG release_2024_05_001
    URL https://elpa.mpcdf.mpg.de/software/tarball-archive/Releases/2024.05.001/elpa-2024.05.001.tar.gz
    PATCH_COMMAND patch -p1 < ${CMAKE_SOURCE_DIR}/cmake/patch/elpa-make-transform-public.patch
    CONFIGURE_COMMAND sh ${conffile}
    BUILD_COMMAND make
    INSTALL_COMMAND ""
    )
endif()
#Now create the elpa target
    add_library(elpa STATIC IMPORTED GLOBAL)
    if (FLEUR_USE_OPENMP)
        set_property(TARGET elpa PROPERTY IMPORTED_LOCATION "${CMAKE_BINARY_DIR}/elpa/src/ELPA/.libs/libelpa_openmp.a")
    else()
        set_property(TARGET elpa PROPERTY IMPORTED_LOCATION "${CMAKE_BINARY_DIR}/elpa/src/ELPA/.libs/libelpa.a")
    endif()  
    

#Now make ELPA known to FLEUR
#if (FLEUR_USE_OPENMP)
#    set(FLEUR_LIBRARIES "-L${CMAKE_BINARY_DIR}/elpa/src/ELPA/.libs;-lelpa_openmp;${FLEUR_LIBRARIES}")
#else()
#    set(FLEUR_LIBRARIES "-L${CMAKE_BINARY_DIR}/elpa/src/ELPA/.libs;-lelpa;${FLEUR_LIBRARIES}")
#endif()    
include_directories("${CMAKE_CURRENT_BINARY_DIR}/elpa/src/ELPA/modules")
set(FLEUR_USE_ELPA TRUE)
set(FLEUR_MPI_DEFINITIONS ${FLEUR_MPI_DEFINITIONS} "CPP_ELPA" "CPP_ELPA_PATCH" "CPP_ELPA_SP")
set(FLEUR_USE_INTERNAL_ELPA true)
