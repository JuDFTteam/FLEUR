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
 set(fleur_F90 ${fleur_F90} diagonalization/elpa_20180525.F90)
 set(diag_test_files diagonalization/elpa_20180525.F90)
else()
 set(fleur_F90 ${fleur_F90} diagonalization/elpa.F90)
 set(diag_test_files diagonalization/elpa.F90)
endif()

if (FLEUR_USE_HDF5)	
set(diag_test_files ${diag_test_files}
diagonalization/diag_test.F90
diagonalization/eigen_diag.F90
diagonalization/lapack_diag.F90
diagonalization/magma.F90
diagonalization/scalapack.F90
#diagonalization/chase_diag.F90
#diagonalization/symmetrize_matrix.f90
diagonalization/cusolver_diag.F90
diagonalization/elemental.F90
diagonalization/available_solvers.F90
diagonalization/writeout.F90
diagonalization/elpa_20180525_onenode.F90
types/types_mat.F90
types/types_mpimat.F90
types/types_gpumat.F90
types/types_setup.F90
io/io_matrix.F90
io/iomatrix_hdf.F90
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
