#first try if Wannier90 already works
try_compile(FLEUR_USE_WANN ${CMAKE_BINARY_DIR} ${CMAKE_SOURCE_DIR}/cmake/tests/test_Wannier90.f90
	    LINK_LIBRARIES ${FLEUR_LIBRARIES}
            )

message("Wannier90 1.2 Library found:${FLEUR_USE_WANN}")

if (FLEUR_USE_WANN)
   set(FLEUR_DEFINITIONS ${FLEUR_DEFINITIONS} "CPP_WANN") 
   set(FLEUR_MPI_DEFINITIONS ${FLEUR_MPI_DEFINITIONS} "CPP_WANN")
endif()
