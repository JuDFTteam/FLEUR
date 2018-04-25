set(fleur_F90 ${fleur_F90}
   mpi/setupMPI.F90)
if (${FLEUR_USE_MPI})
   set(fleur_F90 ${fleur_F90}
   mpi/mingeselle.F90
   mpi/mpi_bc_all.F90
   mpi/mpi_bc_potden.F90
   mpi/mpi_bc_coreDen.F90
   mpi/mpi_bc_st.F90
   mpi/mpi_bc_pot.F90
   mpi/mpi_col_den.F90
   mpi/mpi_make_groups.F90
   mpi/mpi_dist_forcetheorem.F90
   mpi/mpi_bc_xcpot.F90
)
endif()
