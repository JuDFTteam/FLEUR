#first try if Wannier90 modification 4 already works
try_compile(FLEUR_USE_WANN4 ${CMAKE_BINARY_DIR} ${CMAKE_SOURCE_DIR}/cmake/tests/test_Wannier4.f90
	    LINK_LIBRARIES ${FLEUR_LIBRARIES}
            )

message("Wannier90-4 1.2 Library found:${FLEUR_USE_WANN4}")

if (FLEUR_USE_WANN4)
   set(FLEUR_DEFINITIONS ${FLEUR_DEFINITIONS} "CPP_WANN4") 
   set(FLEUR_MPI_DEFINITIONS ${FLEUR_MPI_DEFINITIONS} "CPP_WANN4")
endif()
