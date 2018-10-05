set(fleur_F77 ${fleur_F77}
)
set(fleur_F90 ${fleur_F90}
diagonalization/eigen_diag.F90
diagonalization/lapack_diag.F90
diagonalization/magma.F90
diagonalization/scalapack.F90
diagonalization/chase_diag.F90
diagonalization/symmetrize_matrix.f90
diagonalization/elemental.F90)
if (FLEUR_USE_ELPA_20180525)
 set(fleur_F90 ${fleur_F90}
  diagonalization/elpa_20180525.F90
)
else()
 set(fleur_F90 ${fleur_F90}
  diagonalization/elpa.F90
)
endif()
