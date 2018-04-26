#Check if we can compile with LIBXC
try_compile(FLEUR_USE_LIBXC ${CMAKE_BINARY_DIR} ${CMAKE_SOURCE_DIR}/cmake/tests/test_LibXC.f90
LINK_LIBRARIES ${FLEUR_LIBRARIES}
            )

message("Libxc Library found:${FLEUR_USE_LIBXC}")

if (FLEUR_USE_LIBXC)
   set(FLEUR_MPI_DEFINITIONS ${FLEUR_MPI_DEFINITIONS} "CPP_LIBXC")
   set(FLEUR_DEFINITIONS ${FLEUR_DEFINITIONS} "CPP_LIBXC")
endif()
