#Check if we can compile with LIBXC
try_compile(FLEUR_USE_LIBXC ${CMAKE_BINARY_DIR} ${CMAKE_SOURCE_DIR}/cmake/tests/test_LibXC.f90
LINK_LIBRARIES ${FLEUR_LIBRARIES}
            )

foreach (teststring "-lxcf90;-lxc")
if (NOT FLEUR_USE_LIBXC)
    set(TEST_LIBRARIES "${FLEUR_LIBRARIES};${teststring}")
    try_compile(FLEUR_USE_LIBXC ${CMAKE_BINARY_DIR} ${CMAKE_SOURCE_DIR}/cmake/tests/test_LibXC.f90
    LINK_LIBRARIES ${TEST_LIBRARIES} OUTPUT_VARIABLE compile_output
    )
    if ("$ENV{VERBOSE}")
       message("LIBXC compile test: ${FLEUR_USE_LIBXC}\nLINK_LIBRARIES ${TEST_LIBRARIES}\n${compile_output}")
    endif()
    if (FLEUR_USE_LIBXC)
       set(FLEUR_LIBRARIES ${TEST_LIBRARIES})
    endif()
endif()
endforeach()

if (DEFINED CLI_FLEUR_USE_LIBXC)
    if (CLI_FLEUR_USE_LIBXC)
       if (NOT FLEUR_USE_LIBXC)
           message("You asked for LibXC support but cmake couldn't find it. We will try to download and compile along with FLEUR")
	   if(NOT EXISTS "${PROJECT_SOURCE_DIR}/external/libxc-git/src" )
    	     find_package(Git REQUIRED)
    	     execute_process(COMMAND ${GIT_EXECUTABLE} submodule init external/libxc-git WORKING_DIRECTORY ${PROJECT_SOURCE_DIR} RESULT_VARIABLE _res_init OUTPUT_QUIET ERROR_QUIET)
    	     execute_process(COMMAND ${GIT_EXECUTABLE} submodule update  WORKING_DIRECTORY ${PROJECT_SOURCE_DIR} RESULT_VARIABLE _res_update OUTPUT_QUIET ERROR_QUIET)
    	     if (_res_init GREATER 0 OR _res_update GREATER 0)
               message(FATAL_ERROR "LibXC source could not be downloaded.\n"
                            "We tried: 'git submodule init external/libxc-git && git submodule update' and resulted in error" )
             endif()
	   endif()
	   message("libxc was patched")


           set(ENABLE_FORTRAN ON CACHE BOOL "Build Fortran 90 interface")
           set(ENABLE_FORTRAN03 ON CACHE BOOL "Build Fortran 2003 interface")
	   set(CMAKE_Fortran_MODULE_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/modules/external")
	   add_subdirectory (external/libxc-git EXCLUDE_FROM_ALL)
	   include_directories("${CMAKE_CURRENT_BINARY_DIR}/modules/external")
	   set(FLEUR_USE_LIBXC TRUE)
	   set(FLEUR_LINK_LIBRARIES "${FLEUR_LINK_LIBRARIES};xcf90")
           set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -I${CMAKE_CURRENT_BINARY_DIR}/external/libxc-git")
       endif()
    else()
        if (FLEUR_USE_LIBXC)
           message("LibXC found but you explicitely asked not to use it")
	   set(FLEUR_USE_LIBXC FALSE)
        endif()
    endif()
endif()


message("Libxc Library found:${FLEUR_USE_LIBXC}")

if (FLEUR_USE_LIBXC)
   set(FLEUR_MPI_DEFINITIONS ${FLEUR_MPI_DEFINITIONS} "CPP_LIBXC")
   set(FLEUR_DEFINITIONS ${FLEUR_DEFINITIONS} "CPP_LIBXC")
endif()
