set(fleur_F77 ${fleur_F77}
)
set(fleur_F90 ${fleur_F90}
diagonalization/eigen_diag.F90
diagonalization/franza.F90
diagonalization/geneigbrobl.F90
diagonalization/locrectify.F90
diagonalization/mod_loccoeff.f90
diagonalization/zsymsecloc.F90
diagonalization/packed_to_full.f90
diagonalization/magma.F90 )

if (${Fleur_uses_MPI})
   set(fleur_F90 ${fleur_F90}
      diagonalization/ssubredist1.F90
      diagonalization/ssubredist2.F90)
   if(${Fleur_uses_ELPA}) 
       set(fleur_F90 ${fleur_F90} diagonalization/elpa.F90)
   endif()
   if(${Fleur_uses_SCALAPACK}) 
       set(fleur_F90 ${fleur_F90} diagonalization/chani.F90)
   endif()
   if(${Fleur_uses_ELEMENTAL}) 
       set(fleur_F90 ${fleur_F90} diagonalization/elemental.F90)
   endif()
endif()
