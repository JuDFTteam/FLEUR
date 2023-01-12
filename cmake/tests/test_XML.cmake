#First check if we can compile with XML2
try_compile(FLEUR_USE_XML ${CMAKE_BINARY_DIR} ${CMAKE_SOURCE_DIR}/cmake/tests/test_XML.f90
LINK_LIBRARIES ${FLEUR_LIBRARIES}
            )

if (NOT FLEUR_USE_XML)
      find_package(LibXml2)
      list(TRANSFORM LIBXML2_INCLUDE_DIRS PREPEND -I)
      string (REPLACE ";" " " LIBXML2_INCLUDE_DIRS_STR "${LIBXML2_INCLUDE_DIRS}")
      set(CMAKE_C_FLAGS "${LIBXML2_INCLUDE_DIRS_STR} ${CMAKE_C_FLAGS}")
      if (LIBXML2_LIBRARIES)
          set(TEST_LIBRARIES ${FLEUR_LIBRARIES} ${LIBXML2_LIBRARIES})
      endif()
try_compile(FLEUR_USE_XML ${CMAKE_BINARY_DIR} ${CMAKE_SOURCE_DIR}/cmake/tests/test_XML.f90
	    LINK_LIBRARIES ${TEST_LIBRARIES} OUTPUT_VARIABLE compile_output 
            )
	    if	("$ENV{VERBOSE}")
            	message("XML F90 compile test: ${FLEUR_USE_XML}\nLINK_LIBRARIES ${TEST_LIBRARIES}\n${compile_output}")
     	    endif()
       if (FLEUR_USE_XML)
              set(FLEUR_LIBRARIES ${TEST_LIBRARIES} )
       endif()
endif()

#Try to simply add -lxml2
if (NOT FLEUR_USE_XML)
      set(TEST_LIBRARIES ${FLEUR_LIBRARIES} -lxml2)

try_compile(FLEUR_USE_XML ${CMAKE_BINARY_DIR} ${CMAKE_SOURCE_DIR}/cmake/tests/test_XML.f90
            LINK_LIBRARIES ${TEST_LIBRARIES} OUTPUT_VARIABLE compile_output 
            )
	    if ("$ENV{VERBOSE}")
            	message("XML F90 compile test: ${FLEUR_USE_XML}\nLINK_LIBRARIES ${TEST_LIBRARIES}\n${compile_output}")
     	    endif()
       if (FLEUR_USE_XML)
              set(FLEUR_LIBRARIES -lxml2 ${FLEUR_LIBRARIES})
       endif()
endif()


message("XML Library found for linking:${FLEUR_USE_XML}")

if (FLEUR_USE_XML)
   try_compile(FLEUR_USE_XML ${CMAKE_BINARY_DIR} ${CMAKE_SOURCE_DIR}/cmake/tests/test_XML.c
   CMAKE_FLAGS "-DCMAKE_C_LINK_EXECUTABLE='echo no linking'" LINK_LIBRARIES "-lxml2" OUTPUT_VARIABLE compile_output )
   if ("$ENV{VERBOSE}")
       message("XML C compile test: ${FLEUR_USE_XML}\n${compile_output}")
   endif()
   if (NOT FLEUR_USE_XML)
      find_package(LibXml2)
      list(TRANSFORM LIBXML2_INCLUDE_DIRS PREPEND -I)
      string (REPLACE ";" " " LIBXML2_INCLUDE_DIRS_STR "${LIBXML2_INCLUDE_DIRS}")
      set(CMAKE_C_FLAGS "${LIBXML2_INCLUDE_DIRS_STR} ${CMAKE_C_FLAGS}")
      try_compile(FLEUR_USE_XML ${CMAKE_BINARY_DIR} ${CMAKE_SOURCE_DIR}/cmake/tests/test_XML.c OUTPUT_VARIABLE compile_output
     )
     if ("$ENV{VERBOSE}")
       message("XML C compile test: ${FLEUR_USE_XML}\n INCLUDE DIRECTORIES: ${CMAKE_C_FLAGS}\n${compile_output}")
     endif()
   endif()
endif()

message("XML Library found for C:${FLEUR_USE_XML}")


if (NOT FLEUR_USE_XML)
   if (DEFINED CLI_FLEUR_COMPILE_LIBXML2)
       if (NOT EXISTS "${PROJECT_SOURCE_DIR}/.git")
            message(FATAL_ERROR "You asked for libXML2 to be compiled. This requires the use of the git-version of FLEUR\n")
       endif()     
       message(WARNING "You asked for libxml2 to be compiled. We will try to download and compile libxml2 along with FLEUR")
       if(NOT EXISTS "${PROJECT_SOURCE_DIR}/external/libxml2-git/include" )
    	    find_package(Git REQUIRED)
    	    execute_process(COMMAND ${GIT_EXECUTABLE} submodule init external/libxml2-git WORKING_DIRECTORY ${PROJECT_SOURCE_DIR} RESULT_VARIABLE _res_init OUTPUT_QUIET ERROR_QUIET)
    	    execute_process(COMMAND ${GIT_EXECUTABLE} submodule update  WORKING_DIRECTORY ${PROJECT_SOURCE_DIR} RESULT_VARIABLE _res_update OUTPUT_QUIET ERROR_QUIET)
    	    if (_res_init GREATER 0 OR _res_update GREATER 0)
               message(FATAL_ERROR "libxml2 source could not be downloaded.\n"
                            "We tried: 'git submodule init external/libxml2-git && git submodule update' and resulted in error" )
             endif()
	endif()
	set(CMAKE_Fortran_MODULE_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/modules/external")
	add_subdirectory (external/libxml2-git EXCLUDE_FROM_ALL)
	set(FLEUR_USE_XML TRUE)
       set(JUDFT_COMPILE_LIBXML2 ON CACHE BOOL "compile") 
       set(BUILD_SHARED_LIBS OFF CACHE BOOL "Build shared libraries")
       set(LIBXML2_WITH_AUTOMATA ON)
       set(LIBXML2_WITH_DEBUG OFF CACHE BOOL "Add the debugging module")
       set(LIBXML2_WITH_HTML OFF CACHE BOOL "Add the HTML support")
       set(LIBXML2_WITH_HTTP OFF CACHE BOOL "Add the HTTP support")
       set(LIBXML2_WITH_LZMA OFF CACHE BOOL "Use liblzma")
       set(LIBXML2_WITH_MODULES OFF CACHE BOOL "Add the dynamic modules support")
       set(LIBXML2_WITH_OUTPUT OFF CACHE BOOL "Add the serialization support")
       set(LIBXML2_WITH_PROGRAMS OFF CACHE BOOL "Build programs")
       set(LIBXML2_WITH_PYTHON OFF CACHE BOOL "Build Python bindings")
       set(LIBXML2_WITH_SAX1 OFF CACHE BOOL "Add the older SAX1 interface")
       set(LIBXML2_WITH_TESTS OFF CACHE BOOL  "Build tests")
       set(LIBXML2_WITH_ZLIB OFF CACHE BOOL "Use libz")


       include_directories("${CMAKE_CURRENT_BINARY_DIR}/modules/external")
       include_directories("${CMAKE_CURRENT_BINARY_DIR}/modules/external/static")
    endif()
endif()


