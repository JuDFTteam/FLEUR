set(fleur_F90 ${fleur_F90}
   src/fleur/mpi/mpi_reduce_tool.F90
   src/fleur/mpi/mpi_bc_xcpot.F90
   src/fleur/mpi/setupMPI.F90
   src/fleur/mpi/omp_checker.F90
   src/fleur/mpi/mpi_bc_potden.F90
)

if (FLEUR_USE_MPI)
   set(fleur_F90 ${fleur_F90}
   src/fleur/mpi/mpi_bc_coreDen.F90
   src/fleur/mpi/mpi_bc_st.F90
   src/fleur/mpi/mpi_bc_pot.F90
   src/fleur/mpi/mpi_col_den.F90
   src/fleur/mpi/mpi_reduce_potden.F90
   src/fleur/mpi/mpi_make_groups.F90
   src/fleur/mpi/mpi_dist_forcetheorem.F90
)
endif()

if (FLEUR_USE_SCHED)
set(fleur_F90 ${fleur_F90}
   src/fleur/mpi/wrapper_sched_getcpu.c)
endif()
