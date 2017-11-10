set(fleur_F77 ${fleur_F77}
)
set(fleur_F90 ${fleur_F90}
diagonalization/eigen_diag.F90
diagonalization/lapack_diag.F90
diagonalization/magma.F90 )

if (${FLEUR_USE_MPI})
   set(fleur_F90 ${fleur_F90}
      diagonalization/ssubredist1.F90
      diagonalization/ssubredist2.F90)
   if(${FLEUR_USE_ELPA}) 
       set(fleur_F90 ${fleur_F90} diagonalization/elpa.F90)
   endif()
   if(${FLEUR_USE_SCALAPACK}) 
       set(fleur_F90 ${fleur_F90} diagonalization/chani.F90)
   endif()
   if(${FLEUR_USE_ELEMENTAL}) 
       set(fleur_F90 ${fleur_F90} diagonalization/elemental.F90)
   endif()
endif()
