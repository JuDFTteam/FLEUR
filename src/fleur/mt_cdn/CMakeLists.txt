set(fleur_F77 ${fleur_F77}
)
set(fleur_F90 ${fleur_F90}
src/fleur/mt_cdn/alignSpinAxisMagn.F90
src/fleur/mt_cdn/cdncore.F90
src/fleur/mt_cdn/denMultipoleExp.f90
src/fleur/mt_cdn/efg.F90
src/fleur/mt_cdn/l-like.F90
src/fleur/mt_cdn/magDiMom.f90
src/fleur/mt_cdn/MagMoments.F90
src/fleur/mt_cdn/magMultipoles.F90
src/fleur/mt_cdn/magnMomFromDen.f90
src/fleur/mt_cdn/nIJ_mat.f90
src/fleur/mt_cdn/n_mat.f90
src/fleur/mt_cdn/radfun.f90
src/fleur/mt_cdn/radfun_integrals.f90
src/fleur/mt_cdn/resMoms.f90
src/fleur/mt_cdn/types_denmatrix.F90
src/fleur/mt_cdn/types_radfun.f90
src/fleur/mt_cdn/types_moessbauerParams.F90
)

set(radfun_integrals_test_files
${CMAKE_SOURCE_DIR}/src/fleur/math/intgr.F90
${CMAKE_SOURCE_DIR}/src/fleur/mt_cdn/radfun_integrals.f90
${CMAKE_SOURCE_DIR}/testing/small-testcodes/radfun_integrals_test.f90
)

add_executable(radfun_integrals_test EXCLUDE_FROM_ALL ${radfun_integrals_test_files})
target_link_libraries(radfun_integrals_test juDFT ${FLEUR_LIBRARIES} ${FLEUR_LINK_LIBRARIES})
set_target_properties(radfun_integrals_test PROPERTIES
   Fortran_MODULE_DIRECTORY modules/radfun_integrals_test
   COMPILE_OPTIONS -Imodules/radfun_integrals_test
   RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/tools")
target_compile_options(radfun_integrals_test PUBLIC $<$<COMPILE_LANGUAGE:Fortran>:${FLEUR_PRECISION_OPTION}>)
target_compile_options(radfun_integrals_test PRIVATE $<$<COMPILE_LANGUAGE:Fortran>:${FLEUR_COMPILE_OPTIONS}>)

set(n_mat_radial_overlap_test_files
${CMAKE_SOURCE_DIR}/src/fleur/math/intgr.F90
${CMAKE_SOURCE_DIR}/src/fleur/mt_cdn/radfun_integrals.f90
${CMAKE_SOURCE_DIR}/testing/small-testcodes/n_mat_test_types.f90
${CMAKE_SOURCE_DIR}/src/fleur/mt_cdn/n_mat.f90
${CMAKE_SOURCE_DIR}/testing/small-testcodes/n_mat_radial_overlap_test.f90
)

add_executable(n_mat_radial_overlap_test EXCLUDE_FROM_ALL ${n_mat_radial_overlap_test_files})
target_link_libraries(n_mat_radial_overlap_test juDFT ${FLEUR_LIBRARIES} ${FLEUR_LINK_LIBRARIES})
set_target_properties(n_mat_radial_overlap_test PROPERTIES
   Fortran_MODULE_DIRECTORY modules/n_mat_radial_overlap_test
   COMPILE_OPTIONS -Imodules/n_mat_radial_overlap_test
   RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/tools")
target_compile_options(n_mat_radial_overlap_test PUBLIC $<$<COMPILE_LANGUAGE:Fortran>:${FLEUR_PRECISION_OPTION}>)
target_compile_options(n_mat_radial_overlap_test PRIVATE $<$<COMPILE_LANGUAGE:Fortran>:${FLEUR_COMPILE_OPTIONS}>)
