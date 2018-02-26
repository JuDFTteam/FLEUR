#Check if we can compile with GPU
if ($ENV{FLEUR_USE_GPU})
   #No check is done
   set(FLEUR_USE_GPU TRUE)
   set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -ta=tesla:cuda8.0,cc60 -Mcuda:kepler+  -Minfo=accel -acc ")
   set(FLEUR_MPI_DEFINITIONS ${FLEUR_MPI_DEFINITIONS} "CPP_GPU")
   set(FLEUR_DEFINITIONS ${FLEUR_DEFINITIONS} "CPP_GPU")
else()
   set(FLEUR_USE_GPU FALSE)
endif()