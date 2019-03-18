set(fleur_F77 ${fleur_F77}
)
set(fleur_F90 ${fleur_F90}
diagonalization/eigen_diag.F90
diagonalization/lapack_diag.F90
diagonalization/magma.F90
diagonalization/scalapack.F90
diagonalization/chase_diag.F90
diagonalization/symmetrize_matrix.f90
diagonalization/cusolver_diag.F90
diagonalization/elemental.F90
diagonalization/available_solvers.F90
diagonalization/writeout.F90
diagonalization/elpa_20180525_onenode.F90)
if (FLEUR_USE_ELPA_20180525)
 set(fleur_F90 ${fleur_F90}
  diagonalization/elpa_20180525.F90
)
else()
 set(fleur_F90 ${fleur_F90}
  diagonalization/elpa.F90
)
endif()

set(diag_test_files ${diag_test_files}
diagonalization/eigen_diag.F90
diagonalization/lapack_diag.F90
diagonalization/magma.F90
diagonalization/scalapack.F90
diagonalization/chase_diag.F90
diagonalization/symmetrize_matrix.f90
diagonalization/cusolver_diag.F90
diagonalization/elemental.F90
diagonalization/available_solvers.F90
diagonalization/writeout.F90
diagonalization/elpa_20180525_onenode.F90)


add_executable(diag_test ${diag_test_files})
target_compile_definitions(diag_test PUBLIC ${FLEUR_DEFINITIONS})
target_link_libraries(diag_test ${FLEUR_LIBRARIES})
set_target_properties(diag_test PROPERTIES Fortran_MODULE_DIRECTORY diag_test_modules COMPILE_OPTIONS -Idiag_test_modules)

