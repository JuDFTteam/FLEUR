#Test if the ARPACK library is present (prerequisite for EDsolver)

set(FLEUR_COMPILE_EDSOLVER FALSE)

try_compile(FLEUR_USE_ARPACK ${CMAKE_BINARY_DIR} ${CMAKE_SOURCE_DIR}/cmake/tests/test_ARPACK.f
            LINK_LIBRARIES ${FLEUR_LIBRARIES})

set(ARPACK_ADDED FALSE)
if(FLEUR_USE_ARPACK)
  set(FLEUR_ARPACK_LIBRARIES "")
else()
  #Try to find the library by adding linker options
  foreach(ADD_STRING "-larpack_ifort"
                     "-larpack_gfortran"
                     "-larpack")
     if (NOT FLEUR_USE_ARPACK)
        set(TEST_LIBRARIES "${FLEUR_LIBRARIES};${ADD_STRING}")
        try_compile(FLEUR_USE_ARPACK ${CMAKE_BINARY_DIR} ${CMAKE_SOURCE_DIR}/cmake/tests/test_ARPACK.f
                    LINK_LIBRARIES ${TEST_LIBRARIES})
       if (FLEUR_USE_ARPACK)
            set(FLEUR_ARPACK_LIBRARIES ${ADD_STRING})
            set(ARPACK_ADDED TRUE)
       endif()
     endif()
  endforeach()
endif()
message("ARPACK Library found:${FLEUR_USE_ARPACK}")

#Test if the EDsolver library is already present
try_compile(FLEUR_USE_EDSOLVER ${CMAKE_BINARY_DIR} ${CMAKE_SOURCE_DIR}/cmake/tests/test_EDsolver.f90
            LINK_LIBRARIES "${FLEUR_LIBRARIES};${FLEUR_ARPACK_LIBRARIES}")
if (NOT FLEUR_USE_EDSOLVER)
   if(FLEUR_USE_ARPACK)
      #try adding -lEDsolver to the linker options (before ARPACK)
      set(TEST_LIBRARIES "${FLEUR_LIBRARIES};-lEDsolver;${FLEUR_ARPACK_LIBRARIES}")
      try_compile(FLEUR_USE_EDSOLVER ${CMAKE_BINARY_DIR} ${CMAKE_SOURCE_DIR}/cmake/tests/test_EDsolver.f90
                  LINK_LIBRARIES ${TEST_LIBRARIES})
      if (FLEUR_USE_EDSOLVER)
         #If we are here the EDsolver library is compiled somewhere
         set(FLEUR_LIBRARIES ${TEST_LIBRARIES})
      endif()
   endif()
endif()
message("EDsolver Library found:${FLEUR_USE_EDSOLVER}")

if (DEFINED CLI_FLEUR_USE_EDSOLVER)
   if (CLI_FLEUR_USE_EDSOLVER)
      if (NOT FLEUR_USE_EDSOLVER)
         if (NOT FLEUR_USE_ARPACK)
            message(WARNING "You asked for the EDsolver library but cmake couldn't find the ARPACK library, \
                             which is a prerequisite for the EDsolver. \
                             Please check your configuration or install the ARPACK library.")
         else()
            message(WARNING "You asked for the EDsolver library but cmake couldn't find it. \
                             We will try to download and compile the EDsolver library along with FLEUR")
            if(NOT EXISTS "${PROJECT_SOURCE_DIR}/external/edsolver-library/src")
               find_package(Git REQUIRED)
               execute_process(COMMAND ${GIT_EXECUTABLE} submodule init external/edsolver-library
                               WORKING_DIRECTORY ${PROJECT_SOURCE_DIR} RESULT_VARIABLE _res_init)
               execute_process(COMMAND ${GIT_EXECUTABLE} submodule update
                               WORKING_DIRECTORY ${PROJECT_SOURCE_DIR} RESULT_VARIABLE _res_update)
               if (_res_init GREATER 0 OR _res_update GREATER 0)
                  message(FATAL_ERROR "EDsolver source could not be downloaded. \
                                       We tried: 'git submodule init external/edsolver-library \
                                       && git submodule update' and this resulted in error" )
               endif()
               if(NOT EXISTS "${PROJECT_SOURCE_DIR}/external/edsolver-library/src")
                  #If someone has no access to the repository but tries to to git submodule init/update
                  #It will complete with no error but nothing will happen
                  message(FATAL_ERROR "It seems that you asked for the EDsolver library to be pulled from git. \
                                       This is a private repository. \
                                       If you already have access, please configure your git to log you in automatically. \
                                       If not contact he.janssen@fz-juelich.de")
               endif()
            endif()
         endif()
      endif()
      if (FLEUR_USE_ARPACK)
         set(EDSOLVER_BUILD_FLEUR ON CACHE BOOL "Build EDsolver with juDFT")
         set(EDSOLVER_BUILD_UTILS ON CACHE BOOL "Build auxiliary programs for the EDsolver")
         if (CMAKE_BUILD_TYPE MATCHES Debug)
            set(EDSOLVER_DEBUG_BUILD ON CACHE BOOL "Build a debug version")
         endif()
         add_subdirectory (external/edsolver-library)
         set(FLEUR_USE_EDSOLVER TRUE)
         set(FLEUR_COMPILE_EDSOLVER TRUE)
         include_directories("${CMAKE_CURRENT_BINARY_DIR}/external/edsolver-library/")
         include_directories("${CMAKE_CURRENT_BINARY_DIR}/external/edsolver-library/modules")
      endif()
   else()
      if (FLEUR_USE_EDSOLVER)
         message("EDsolver library found, but you explicitly asked not to use it")
         set(FLEUR_USE_EDSOLVER FALSE)
      endif()
   endif()
endif()

if (FLEUR_USE_EDSOLVER)
   if (FLEUR_COMPILE_EDSOLVER)
      if (ARPACK_ADDED)
        set(FLEUR_LIBRARIES "${FLEUR_LIBRARIES};${FLEUR_ARPACK_LIBRARIES}")
      endif()
   endif()
   set(FLEUR_DEFINITIONS ${FLEUR_DEFINITIONS} "CPP_EDSOLVER")
   set(FLEUR_MPI_DEFINITIONS ${FLEUR_MPI_DEFINITIONS} "CPP_EDSOLVER")
endif()
