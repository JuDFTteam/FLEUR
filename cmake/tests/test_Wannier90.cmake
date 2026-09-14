#first try if Wannier90 already works
set(_old_try_compile_target_type ${CMAKE_TRY_COMPILE_TARGET_TYPE})
set(CMAKE_TRY_COMPILE_TARGET_TYPE EXECUTABLE)

try_compile(FLEUR_USE_WANN ${CMAKE_BINARY_DIR} ${CMAKE_SOURCE_DIR}/cmake/tests/test_Wannier90.f90
            LINK_LIBRARIES ${FLEUR_LIBRARIES}
            )
set(FLEUR_WANNIER90_LIBRARIES ${FLEUR_LIBRARIES}) 
foreach(ADD_String "-lwannier" )
   if (NOT FLEUR_USE_WANN)
     set(TEST_LIBRARIES "${ADD_String};${FLEUR_LIBRARIES}")

     try_compile(FLEUR_USE_WANN ${CMAKE_BINARY_DIR} ${CMAKE_SOURCE_DIR}/cmake/tests/test_Wannier90.f90
            LINK_LIBRARIES ${TEST_LIBRARIES} OUTPUT_VARIABLE compile_output
            )
     if ("$ENV{VERBOSE}")
        message("Wannier90 compile test: ${FLEUR_USE_WANN}\nLINK_LIBRARIES ${TEST_LIBRARIES}\n${compile_output}")
     endif()

     if (FLEUR_USE_WANN)
          set(FLEUR_WANNIER90_LIBRARIES ${TEST_LIBRARIES})
     endif()
   endif()
endforeach()

set(CMAKE_TRY_COMPILE_TARGET_TYPE ${_old_try_compile_target_type})
unset(_old_try_compile_target_type)

message("Wannier90 1.2 Library found:${FLEUR_USE_WANN}")

set(_WANNIER90_SUBPROJECT_OK FALSE)
if (EXISTS "${PROJECT_SOURCE_DIR}/external/wannier90/CMakeLists.txt")
   set(_WANNIER90_SUBPROJECT_OK TRUE)
endif()

set(FLEUR_USE_WANNLIB_API FALSE)

if (DEFINED CLI_FLEUR_USE_WANNIER)
   if (CLI_FLEUR_USE_WANNIER)
       if (NOT FLEUR_USE_WANN)
           if (NOT _WANNIER90_SUBPROJECT_OK)
              if (NOT EXISTS "${PROJECT_SOURCE_DIR}/.git")
                 message(FATAL_ERROR "You asked for Wannier90 to be used, but it cannot be found.\n"
                                     "Please either provide correct include and link directories for Wannier90 manually, or use a git-version of FLEUR to download Wannier90 automatically")
              endif()
              find_package(Git REQUIRED)
              execute_process(COMMAND ${GIT_EXECUTABLE} submodule init external/wannier90
                              WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}
                              RESULT_VARIABLE _res_init
                              OUTPUT_QUIET
                              ERROR_QUIET)
              execute_process(COMMAND ${GIT_EXECUTABLE} submodule update --init external/wannier90
                              WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}
                              RESULT_VARIABLE _res_update
                              OUTPUT_QUIET
                              ERROR_QUIET)
              if (_res_init GREATER 0 OR _res_update GREATER 0)
                 message(FATAL_ERROR "Wannier90 source could not be downloaded.\n"
                                     "We tried: 'git submodule init external/wannier90 && git submodule update --init external/wannier90' and it failed")
              endif()
              if (EXISTS "${PROJECT_SOURCE_DIR}/external/wannier90/CMakeLists.txt")
                 set(_WANNIER90_SUBPROJECT_OK TRUE)
              endif()
           endif()

           if (_WANNIER90_SUBPROJECT_OK)
              message(WARNING "You asked for Wannier90 but cmake couldn't find it. We will try to compile Wannier90 along with FLEUR")
              set(WANNIER90_MPI ${FLEUR_USE_MPI} CACHE BOOL "Build internal Wannier90 with MPI support" FORCE)
              set(WANNIER90_MPIH FALSE CACHE BOOL "Build internal Wannier90 using mpi.h" FORCE)
              set(WANNIER90_SHARED_LIBS FALSE CACHE BOOL "Build internal Wannier90 as shared library" FORCE)
              set(WANNIER90_INSTALL FALSE CACHE BOOL "Install internal Wannier90" FORCE)
              set(WANNIER90_TEST FALSE CACHE BOOL "Build internal Wannier90 tests" FORCE)
              if (NOT TARGET Wannier90_lib)
                 add_subdirectory(external/wannier90)
              endif()
              set(FLEUR_WANNIER90_LIBRARIES Wannier90_lib ${FLEUR_LIBRARIES})
              # Internal Wannier90 provides the module API (w90_library).
              # Keep legacy API disabled to avoid unresolved wannier_setup/run symbols.
              set(FLEUR_USE_WANN FALSE)
              set(FLEUR_USE_WANNLIB_API TRUE)
           endif()
       endif()
   else()
       if (FLEUR_USE_WANN)
           message("Wannier library found, but you explicitely asked not to use it")
	   set(FLEUR_USE_WANN FALSE)
       endif()
   endif()
endif()

unset(_WANNIER90_SUBPROJECT_OK)

if (FLEUR_USE_WANN)
   if (NOT FLEUR_USE_WANNLIB_API)
      set(_old_try_compile_target_type ${CMAKE_TRY_COMPILE_TARGET_TYPE})
      set(CMAKE_TRY_COMPILE_TARGET_TYPE EXECUTABLE)

      try_compile(FLEUR_USE_WANNLIB_API ${CMAKE_BINARY_DIR} ${CMAKE_SOURCE_DIR}/cmake/tests/test_Wannier90_module.f90
                  LINK_LIBRARIES ${FLEUR_WANNIER90_LIBRARIES}
                  OUTPUT_VARIABLE compile_output_wannlib_api)

      set(CMAKE_TRY_COMPILE_TARGET_TYPE ${_old_try_compile_target_type})
      unset(_old_try_compile_target_type)

      if ("$ENV{VERBOSE}")
         message("Wannier90 module API compile test: ${FLEUR_USE_WANNLIB_API}\nLINK_LIBRARIES ${FLEUR_WANNIER90_LIBRARIES}\n${compile_output_wannlib_api}")
      endif()
   endif()
endif()

message("Wannier90 module API found:${FLEUR_USE_WANNLIB_API}")

if (FLEUR_USE_WANN OR FLEUR_USE_WANNLIB_API)
   set(FLEUR_LIBRARIES ${FLEUR_WANNIER90_LIBRARIES})
endif()

if (FLEUR_USE_WANN)
   set(FLEUR_DEFINITIONS ${FLEUR_DEFINITIONS} "CPP_WANN")
   set(FLEUR_MPI_DEFINITIONS ${FLEUR_MPI_DEFINITIONS} "CPP_WANN")
endif()

if (FLEUR_USE_WANNLIB_API)
   set(FLEUR_DEFINITIONS ${FLEUR_DEFINITIONS} "CPP_WANNLIB_API")
   set(FLEUR_MPI_DEFINITIONS ${FLEUR_MPI_DEFINITIONS} "CPP_WANNLIB_API")
endif()
