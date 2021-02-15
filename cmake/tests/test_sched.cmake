#Test if sched.h can be used to test OpenMP parallelism
try_compile(FLEUR_USE_SCHED ${CMAKE_BINARY_DIR} ${CMAKE_SOURCE_DIR}/cmake/tests/test_sched.c
LINK_LIBRARIES ${FLEUR_LIBRARIES} )

if(CMAKE_Fortran_COMPILER_ID MATCHES "PGI")
	set(FLEUR_USE_SCHED FALSE)
endif()

if (FLEUR_USE_SCHED)
   message("sched.h used")
   set(FLEUR_MPI_DEFINITIONS ${FLEUR_MPI_DEFINITIONS} "CPP_SCHED")
   set(FLEUR_DEFINITIONS ${FLEUR_DEFINITIONS} "CPP_SCHED")
else()
   message("sched.h NOT used")
endif()
