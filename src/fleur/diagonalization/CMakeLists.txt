set(fleur_F77 ${fleur_F77}
)
set(fleur_F90 ${fleur_F90}
src/fleur/diagonalization/eigen_diag.F90
src/fleur/diagonalization/lapack_diag.F90
src/fleur/diagonalization/lapack_singlePrec_diag.F90
src/fleur/diagonalization/dummy_diag.F90
src/fleur/diagonalization/magma.F90
src/fleur/diagonalization/scalapack.F90
src/fleur/diagonalization/chase_diag.F90
src/fleur/diagonalization/symmetrize_matrix.f90
src/fleur/diagonalization/cusolver_diag.F90
src/fleur/diagonalization/elemental.F90
src/fleur/diagonalization/available_solvers.F90
src/fleur/diagonalization/writeout.F90
src/fleur/diagonalization/elpa_20180525_onenode.F90
src/fleur/diagonalization/elsi.F90)
if (FLEUR_USE_ELPA_20180525)
 set(fleur_F90 ${fleur_F90} src/fleur/diagonalization/elpa_20180525.F90)
 set(diag_test_files src/fleur/diagonalization/elpa_20180525.F90)
else()
 set(fleur_F90 ${fleur_F90} src/fleur/diagonalization/elpa.F90)
 set(diag_test_files src/fleur/diagonalization/elpa.F90)
endif()

if (FLEUR_USE_HDF5)	
set(diag_test_files ${diag_test_files}
src/fleur/diagonalization/diag_test.F90
src/fleur/diagonalization/eigen_diag.F90
src/fleur/diagonalization/lapack_diag.F90
src/fleur/diagonalization/magma.F90
src/fleur/diagonalization/scalapack.F90
#src/fleur/diagonalization/chase_diag.F90
#src/fleur/diagonalization/symmetrize_matrix.f90
src/fleur/diagonalization/cusolver_diag.F90
src/fleur/diagonalization/elemental.F90
src/fleur/diagonalization/available_solvers.F90
src/fleur/diagonalization/writeout.F90
src/fleur/diagonalization/elpa_20180525_onenode.F90
src/fleur/types/types_mat.F90
src/fleur/types/types_mpimat.F90
src/fleur/types/types_setup.F90
src/fleur/io/io_matrix.F90
src/fleur/io/iomatrix_hdf.F90
)

add_executable(diag_test EXCLUDE_FROM_ALL ${diag_test_files}  )
target_link_libraries(diag_test ${FLEUR_LIBRARIES})
target_link_libraries(diag_test juDFT)
if (FLEUR_USE_MPI)
target_compile_definitions(diag_test PUBLIC ${FLEUR_MPI_DEFINITIONS})
else()
target_compile_definitions(diag_test PUBLIC ${FLEUR_DEFINITIONS})
endif()
set_target_properties(diag_test PROPERTIES Fortran_MODULE_DIRECTORY modules/diag_test COMPILE_OPTIONS -Imodules/diag_test RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/tools")

endif()
