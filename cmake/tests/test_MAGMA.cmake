#First check if we can compile with MAGMA
if ($ENV{FLEUR_USE_MAGMA})
   message("Set FLEUR_USE_MAGMA to environment, skipping test")
   set(FLEUR_USE_MAGMA $ENV{FLEUR_USE_MAGMA})
else()
try_compile(FLEUR_USE_MAGMA ${CMAKE_BINARY_DIR} ${CMAKE_SOURCE_DIR}/cmake/tests/test_MAGMA.f90
LINK_LIBRARIES ${FLEUR_LIBRARIES}
            )
endif()
message("MAGMA Library found:${FLEUR_USE_MAGMA}")

if (FLEUR_USE_MAGMA)
   set(FLEUR_DEFINITIONS ${FLEUR_DEFINITIONS} "CPP_MAGMA")
endif() 
