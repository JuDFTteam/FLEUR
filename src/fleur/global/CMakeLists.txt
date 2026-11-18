set(fleur_F77 ${fleur_F77}
src/fleur/global/fitchk.f
src/fleur/global/radsra.f
src/fleur/global/radsrd.F
src/fleur/global/radsrdn.f
src/fleur/global/soc_sym.f
src/fleur/global/vacudz.f
src/fleur/global/vacuz.f
src/fleur/global/abcrot.f
src/fleur/global/differ.f
)
set(fleur_F90 ${fleur_F90}
src/fleur/global/savewigner.f90
#src/fleur/global/differ.f90
#src/fleur/global/matrix_copy.F90
src/fleur/global/checkdop.F90
src/fleur/global/checkdopall.f90
src/fleur/global/sort.f90
src/fleur/global/genMTBasis.f90
src/fleur/global/chkmt.f90
src/fleur/global/convn.f90
src/fleur/global/phasy1.f90
src/fleur/global/qfix.f90
src/fleur/global/radflo.F90
src/fleur/global/relLO_dirac.f90
src/fleur/global/starf.f90
src/fleur/global/utility.F90
src/fleur/global/find_enpara.f90
src/fleur/global/rotMMPmat.f90
src/fleur/global/symMMPmat.f90
src/fleur/global/spin_sym.f90
)
if (FLEUR_USE_GPU)
  set(fleur_F90 ${fleur_F90}

  )
endif()

#Regression test for the spin rotation in rotMMPmat_dwgn (read-after-write in the in-place loop).
#Not built by default; testing/tests/invariants/test_rotmmpmat_invariants.py builds and runs it.
set(rotmmpmat_spin_test_files
${CMAKE_SOURCE_DIR}/src/fleur/global/rotMMPmat.f90
${CMAKE_SOURCE_DIR}/src/fleur/types/types_nococonv.F90
${CMAKE_SOURCE_DIR}/src/fleur/types/types_potden.F90
${CMAKE_SOURCE_DIR}/src/fleur/types/types_sphhar.f90
${CMAKE_SOURCE_DIR}/src/fleur/types/types_stars.f90
${CMAKE_SOURCE_DIR}/src/fleur/init/boxdim.f90
${CMAKE_SOURCE_DIR}/src/fleur/init/spgrot.f90
${CMAKE_SOURCE_DIR}/src/fleur/init/compile_descr.F90
${CMAKE_SOURCE_DIR}/src/fleur/global/sort.f90
${CMAKE_SOURCE_DIR}/src/fleur/global/utility.F90
${CMAKE_SOURCE_DIR}/src/fleur/math/intgr.F90
${CMAKE_SOURCE_DIR}/src/fleur/math/pol_angle.f90
${CMAKE_SOURCE_DIR}/src/fleur/io/xmlOutput.F90
${CMAKE_SOURCE_DIR}/testing/small-testcodes/rotmmpmat_spin_test.F90
)
add_executable(rotmmpmat_spin_test EXCLUDE_FROM_ALL ${rotmmpmat_spin_test_files})
add_dependencies(rotmmpmat_spin_test buildinfo)
target_link_libraries(rotmmpmat_spin_test fleurinput juDFT ${FLEUR_LIBRARIES} ${FLEUR_LINK_LIBRARIES})
target_compile_definitions(rotmmpmat_spin_test PUBLIC ${FLEUR_DEFINITIONS} ${FLEUR_PARALLEL_DEFINITIONS})
set_target_properties(rotmmpmat_spin_test PROPERTIES
   Fortran_MODULE_DIRECTORY modules/rotmmpmat_spin_test
   COMPILE_OPTIONS -Imodules/rotmmpmat_spin_test
   RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/tools")
target_compile_options(rotmmpmat_spin_test PUBLIC $<$<COMPILE_LANGUAGE:Fortran>:${FLEUR_PRECISION_OPTION}>)
target_compile_options(rotmmpmat_spin_test PRIVATE $<$<COMPILE_LANGUAGE:Fortran>:${FLEUR_COMPILE_OPTIONS}>)
