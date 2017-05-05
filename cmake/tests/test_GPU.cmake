#Check if we can compile with GPU
if ($ENV{FLEUR_USE_GPU})
   #No check is done
   set(FLEUR_USE_GPU TRUE)
   set(FLEUR_MPI_DEFINITIONS ${FLEUR_MPI_DEFINITIONS} "CPP_GPU")
   set(FLEUR_DEFINITIONS ${FLEUR_DEFINITIONS} "CPP_GPU")
else()
   set(FLEUR_USE_GPU FALSE)
endif()