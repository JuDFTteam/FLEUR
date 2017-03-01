#First check if we can compile with ELPA
try_compile(FLEUR_USE_MAGMA ${CMAKE_BINARY_DIR} ${CMAKE_SOURCE_DIR}/cmake/tests/test_MAGMA.f90
LINK_LIBRARIES ${FLEUR_LIBRARIES}
            )

message("MAGMA Library found:${FLEUR_USE_MAGMA}")

if (FLEUR_USE_MAGMA)
   set(FLEUR_DEFINITIONS ${FLEUR_DEFINITIONS} "CPP_MAGMA")
endif() 
