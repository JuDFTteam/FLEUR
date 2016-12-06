#First check if we can compile with LAPACK
try_compile(FLEUR_USE_LAPACK ${CMAKE_BINARY_DIR} ${CMAKE_SOURCE_DIR}/cmake/tests/test_LAPACK.f90
LINK_LIBRARIES ${FLEUR_LIBRARIES}
            )

if (NOT FLEUR_USE_LAPACK)
      find_package(LAPACK)
      set(TEST_LIBRARIES ${FLEUR_LIBRARIES} ${LAPACK_LIBRARIES})
 
try_compile(FLEUR_USE_LAPACK ${CMAKE_BINARY_DIR} ${CMAKE_SOURCE_DIR}/cmake/tests/test_LAPACK.f90
   	    LINK_LIBRARIES ${TEST_LIBRARIES}
            )
       if (FLEUR_USE_LAPACK)
              set(FLEUR_LIBRARIES ${FLEUR_LIBRARIES} ${LAPACK_LIBRARIES})
              set(FLEUR_MPI_LIBRARIES ${FLEUR_MPI_LIBRARIES} ${LAPACK_LIBRARIES})
       endif()
endif()       

message("LAPACK Library found:${FLEUR_USE_LAPACK}")