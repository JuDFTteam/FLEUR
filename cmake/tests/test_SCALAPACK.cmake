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

#By default we always compile SCALAPACK if not found
if (NOT DEFINED CLI_FLEUR_USE_SCALAPACK)
  set(CLI_FLEUR_USE_SCALAPACK FLEUR_USE_GITVERSION)
endif() 


#Now download SCALAPACK and compile it if REQUIRED
if (DEFINED CLI_FLEUR_USE_SCALAPACK)
   if (CLI_FLEUR_USE_SCALAPACK)
       if (NOT FLEUR_USE_SCALAPACK)
           if (NOT EXISTS "${PROJECT_SOURCE_DIR}/.git")
                message(FATAL_ERROR "You asked for SCALAPACK to be used, but it cannot be found.\n"
                "Please either provide correct include and link directories for SCALAPACK manually, or use a git-version of FLEUR to download SCALAPACK automatically")
           endif()     
           message(WARNING "You asked for SCALAPACK but cmake couldn't find it. We will try to download and compile SCALAPACK along with FLEUR")
           if(NOT EXISTS "${PROJECT_SOURCE_DIR}/external/SCALAPACK-git/src" )
    	     find_package(Git REQUIRED)
    	     execute_process(COMMAND ${GIT_EXECUTABLE} submodule init external/SCALAPACK-git WORKING_DIRECTORY ${PROJECT_SOURCE_DIR} RESULT_VARIABLE _res_init OUTPUT_QUIET ERROR_QUIET)
    	     execute_process(COMMAND ${GIT_EXECUTABLE} submodule update  WORKING_DIRECTORY ${PROJECT_SOURCE_DIR} RESULT_VARIABLE _res_update OUTPUT_QUIET ERROR_QUIET)
    	     if (_res_init GREATER 0 OR _res_update GREATER 0)
               message(FATAL_ERROR "SCALAPACK source could not be downloaded.\n"
                            "We tried: 'git submodule init external/SCALAPACK-git && git submodule update' and resulted in error" )
             endif()
            endif()
            # 
            add_subdirectory(external/SCALAPACK-git EXCLUDE_FROM_ALL)
            set(FLEUR_COMPILE_SCALAPACK TRUE)
            set(FLEUR_USE_SCALAPACK TRUE)
        endif()
    endif()

endif()        




message("SCALAPACK Library found:${FLEUR_USE_SCALAPACK}")
if (FLEUR_USE_SCALAPACK)
   set(FLEUR_MPI_DEFINITIONS ${FLEUR_MPI_DEFINITIONS} "CPP_SCALAPACK")
endif()
