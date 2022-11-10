#Check if we can compile with GPU
if (CLI_FLEUR_USE_GPU)
   #No check is done
   set(FLEUR_USE_GPU TRUE)
   message("GPU:${CLI_FLEUR_USE_GPU}")
   if(${CLI_FLEUR_USE_GPU} MATCHES "acc-gcc")
      set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -fopenacc ")
   elseif(${CLI_FLEUR_USE_GPU} MATCHES "acc")
      set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -acc -Mcuda -Mcudalib=cublas,cufft,cusolver -Minfo=accel -lnvToolsExt")
   elseif(${CLI_Fortran_FLAGS} MATCHES "omp")
      #We try to use OMP offloading
      set(FLEUR_MPI_DEFINITIONS ${FLEUR_MPI_DEFINITIONS} "CPP_OMP_GPU='$omp'")
      set(FLEUR_DEFINITIONS ${FLEUR_DEFINITIONS} "CPP_OMP_GPU='$omp'")
  else()
      message(ERROR,"Choose a GPU mode")
   endif()
   #Check if a CC is given
   STRING(REGEX MATCH ".*:(cc..)" CC_MODE "${CLI_FLEUR_USE_GPU}")
   message("CC:${CC_MODE} ${CMAKE_MATCH_1}")
   if (CC_MODE)
       set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -ta=tesla:${CMAKE_MATCH_1}")
   endif()
#Now check for cusolverDN library
#   set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -Mcuda -ta=tesla,cuda9.1 ")
#   try_compile(FLEUR_USE_CUSOLVER ${CMAKE_BINARY_DIR} ${CMAKE_SOURCE_DIR}/cmake/tests/test_cusolver.c
#	    LINK_LIBRARIES "-lcusolver"
#            )
#   if (FLEUR_USE_CUSOLVER)
#     set(FLEUR_LIBRARIES "${FLEUR_LIBRARIES};-lcusolver")
#     set(FLEUR_MPI_DEFINITIONS ${FLEUR_MPI_DEFINITIONS} "CPP_CUSOLVER")
#     set(FLEUR_DEFINITIONS ${FLEUR_DEFINITIONS} "CPP_CUSOLVER")
#   endif()
else()
   set(FLEUR_USE_GPU FALSE)
   #if we do not use GPU-code we should use OpenMP-on the CPU instead
   if (FLEUR_USE_OPENMP)
      set(FLEUR_MPI_DEFINITIONS ${FLEUR_MPI_DEFINITIONS} "CPP_OMP_CPU='$omp'")
      set(FLEUR_DEFINITIONS ${FLEUR_DEFINITIONS} "CPP_OMP_CPU='$omp'")
   endif()   
endif()
