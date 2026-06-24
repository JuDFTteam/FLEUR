#First check if we can compile with ELSI
try_compile(FLEUR_USE_ELSI ${CMAKE_BINARY_DIR} ${CMAKE_CURRENT_LIST_DIR}/test_ELSI.f90 LINK_LIBRARIES ${FLEUR_LIBRARIES} OUTPUT_VARIABLE compile_output)
if ("$ENV{VERBOSE}")
    message("ELSI compile test: ${FLEUR_USE_ELSI}\nLINK_LIBRARIES ${TEST_LIBRARIES}\n${compile_output}")
endif()

#By default we always compile ELSI if not found
#if (NOT DEFINED CLI_FLEUR_USE_ELSI)
#  set(CLI_FLEUR_USE_ELSI FLEUR_USE_GITVERSION)
#endif()


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

            include(ExternalProject)
            set(ELSI_PREFIX_DIR "${CMAKE_BINARY_DIR}/elsi")
            set(ELSI_BINARY_DIR "${ELSI_PREFIX_DIR}/build")
            set(ELSI_CMAKE_ARGS
                -DCMAKE_BUILD_TYPE=${CMAKE_BUILD_TYPE}
                -DCMAKE_Fortran_COMPILER=${CMAKE_Fortran_COMPILER}
                -DCMAKE_C_COMPILER=${CMAKE_C_COMPILER}
                -DCMAKE_CXX_COMPILER=${CMAKE_CXX_COMPILER}
                -DCMAKE_Fortran_FLAGS=${CMAKE_Fortran_FLAGS}
                -DCMAKE_C_FLAGS=${CMAKE_C_FLAGS}
                -DCMAKE_CXX_FLAGS=${CMAKE_CXX_FLAGS}
                -DENABLE_CHASE=ON
            )
            if (FLEUR_USE_GPU)
               list(APPEND ELSI_CMAKE_ARGS -DUSE_GPU_CUDA=ON)
               if (CC_MODE)
                  list(APPEND ELSI_CMAKE_ARGS -DCMAKE_CUDA_ARCHITECTURES=${CMAKE_MATCH_1})
               endif()
            endif()
            if (FLEUR_COMPILE_SCALAPACK)
               list(APPEND ELSI_CMAKE_ARGS
                    -DSCALAPACK_DIR=${CMAKE_BINARY_DIR}/external/SCALAPACK-git
                    -DSCALAPACK_LIBRARIES=${CMAKE_BINARY_DIR}/external/SCALAPACK-git/lib/libscalapack.a)
               ExternalProject_Add(ELSI
                   PREFIX elsi/
                   SOURCE_DIR ${PROJECT_SOURCE_DIR}/external/elsi-git
                   BINARY_DIR ${ELSI_BINARY_DIR}
                   CMAKE_ARGS ${ELSI_CMAKE_ARGS}
                   BUILD_COMMAND ${CMAKE_COMMAND} --build <BINARY_DIR> --target elsi
                   INSTALL_COMMAND ""
                   BUILD_BYPRODUCTS ${ELSI_BINARY_DIR}/lib/libelsi.a
                   DEPENDS scalapack
               )
            else()
               ExternalProject_Add(ELSI
                   PREFIX elsi/
                   SOURCE_DIR ${PROJECT_SOURCE_DIR}/external/elsi-git
                   BINARY_DIR ${ELSI_BINARY_DIR}
                   CMAKE_ARGS ${ELSI_CMAKE_ARGS}
                   BUILD_COMMAND ${CMAKE_COMMAND} --build <BINARY_DIR> --target elsi
                   INSTALL_COMMAND ""
                   BUILD_BYPRODUCTS ${ELSI_BINARY_DIR}/lib/libelsi.a
               )
            endif()

            add_library(elsi STATIC IMPORTED GLOBAL)
            add_dependencies(elsi ELSI)
            set_property(TARGET elsi PROPERTY IMPORTED_LOCATION "${ELSI_BINARY_DIR}/lib/libelsi.a")
            include_directories("${ELSI_BINARY_DIR}/include")
            include_directories("${ELSI_BINARY_DIR}/modules")
            include_directories("${ELSI_BINARY_DIR}/mod")
            set(FLEUR_ELSI_IMPORTED TRUE)
            
            set(FLEUR_COMPILE_ELSI TRUE)
            set(FLEUR_USE_ELSI TRUE)
        endif()
    endif()

endif()        

message("ELSI Library found:${FLEUR_USE_ELSI}")

if (FLEUR_USE_ELSI)
   set(FLEUR_DEFINITIONS ${FLEUR_DEFINITIONS} "CPP_ELSI")
endif()
