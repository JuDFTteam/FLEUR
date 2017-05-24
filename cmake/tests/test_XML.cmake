#First check if we can compile with XML2
try_compile(FLEUR_USE_XML ${CMAKE_BINARY_DIR} ${CMAKE_SOURCE_DIR}/cmake/tests/test_XML.f90
LINK_LIBRARIES ${FLEUR_LIBRARIES}
            )

if (NOT FLEUR_USE_XML)
      find_package(LibXml2)
      set(CMAKE_C_FLAGS "-I${LIBXML2_INCLUDE_DIR}")
      set(TEST_LIBRARIES ${FLEUR_LIBRARIES} ${LIBXML2_LIBRARIES})
 
try_compile(FLEUR_USE_XML ${CMAKE_BINARY_DIR} ${CMAKE_SOURCE_DIR}/cmake/tests/test_XML.f90
	    LINK_LIBRARIES ${TEST_LIBRARIES}
            )
       if (FLEUR_USE_XML)
              set(FLEUR_LIBRARIES ${LIBXML2_LIBRARIES} ${FLEUR_LIBRARIES})
	      set(FLEUR_MPI_LIBRARIES ${LIBXML2_LIBRARIES} ${FLEUR_MPI_LIBRARIES})
       endif()
endif()       

message("XML Library found for linking:${FLEUR_USE_XML}")

if (FLEUR_USE_XML)
   try_compile(FLEUR_USE_XML ${CMAKE_BINARY_DIR} ${CMAKE_SOURCE_DIR}/cmake/tests/test_XML.c
   LINK_LIBRARIES "-lxml2")
   if (NOT FLEUR_USE_XML)
      find_package(LibXml2)
      set(CMAKE_C_FLAGS "-I${LIBXML2_INCLUDE_DIR}")
      try_compile(FLEUR_USE_XML ${CMAKE_BINARY_DIR} ${CMAKE_SOURCE_DIR}/cmake/tests/test_XML.c
      LINK_LIBRARIES ${LIBXML2_LIBRARIES})
   endif()
endif()

message("XML Library found for C:${FLEUR_USE_XML}")


if (FLEUR_USE_XML)
   set(FLEUR_DEFINITIONS ${FLEUR_DEFINITIONS} "CPP_XML") 
   set(FLEUR_MPI_DEFINITIONS ${FLEUR_MPI_DEFINITIONS} "CPP_XML")
endif()
