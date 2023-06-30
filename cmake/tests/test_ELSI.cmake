#First check if we can compile with ELSI
try_compile(FLEUR_USE_ELSI ${CMAKE_BINARY_DIR} ${CMAKE_SOURCE_DIR}/cmake/tests/test_ELSI.f90 LINK_LIBRARIES ${FLEUR_LIBRARIES} OUTPUT_VARIABLE compile_output)
if ("$ENV{VERBOSE}")
    message("ELSI compile test: ${FLEUR_USE_ELSI}\nLINK_LIBRARIES ${TEST_LIBRARIES}\n${compile_output}")
endif()

#Now download ELSI and compile it if REQUIRED
if (DEFINED CLI_FLEUR_USE_ELSI)
   if (CLI_FLEUR_USE_ELSI)
       if (NOT FLEUR_USE_ELSI)
           if (NOT EXISTS "${PROJECT_SOURCE_DIR}/.git")
                message(FATAL_ERROR "You asked for ELSI to be used, but it cannot be found.\n"
                "Please either provide correct include and link directories for ELSI manually, or use a git-version of FLEUR to download ELSI automatically")
           endif()     
           message(WARNING "You asked for ELSI but cmake couldn't find it. We will try to download and compile ELSI along with FLEUR")
           if(NOT EXISTS "${PROJECT_SOURCE_DIR}/external/elsi-git/src" )
    	     find_package(Git REQUIRED)
    	     execute_process(COMMAND ${GIT_EXECUTABLE} submodule init external/elsi-git WORKING_DIRECTORY ${PROJECT_SOURCE_DIR} RESULT_VARIABLE _res_init OUTPUT_QUIET ERROR_QUIET)
    	     execute_process(COMMAND ${GIT_EXECUTABLE} submodule update  WORKING_DIRECTORY ${PROJECT_SOURCE_DIR} RESULT_VARIABLE _res_update OUTPUT_QUIET ERROR_QUIET)
    	     if (_res_init GREATER 0 OR _res_update GREATER 0)
               message(FATAL_ERROR "ELSI source could not be downloaded.\n"
                            "We tried: 'git submodule init external/elsi-git && git submodule update' and resulted in error" )
             endif()
            endif()
            #
            set(ENABLE_CHASE ON CACHE BOOL "Build Chase support") 
            if (FLEUR_USE_GPU)
               set(USE_GPU_CUDA ON CACHE BOOL "Build GPU support") 
               STRING(REGEX MATCH ".*:cc(..)" CC_MODE "${CLI_FLEUR_USE_GPU}")
               message("ELSI CC:${CC_MODE} ${CMAKE_MATCH_1}")
               if (CC_MODE)
                  set(CMAKE_CUDA_ARCHITECTURES "${CMAKE_MATCH_1}")
                endif()
               
            endif()   
            add_subdirectory(external/elsi-git EXCLUDE_FROM_ALL)
            set(FLEUR_COMPILE_ELSI TRUE)
            set(FLEUR_USE_ELSI TRUE)
        endif()
    endif()

endif()        

message("ELSI Library found:${FLEUR_USE_ELSI}")

if (FLEUR_USE_ELSI)
   set(FLEUR_DEFINITIONS ${FLEUR_DEFINITIONS} "CPP_ELSI")
endif()
