set(fleur_F90 ${fleur_F90}
${CMAKE_SOURCE_DIR}/src/fleur/propcalc/outcdn.f90
${CMAKE_SOURCE_DIR}/src/fleur/propcalc/plot.F90
${CMAKE_SOURCE_DIR}/src/fleur/propcalc/mcd_init.f90
${CMAKE_SOURCE_DIR}/src/fleur/propcalc/local_d_spin_projector_core.F90
${CMAKE_SOURCE_DIR}/src/fleur/propcalc/local_d_spin_projector.F90
${CMAKE_SOURCE_DIR}/src/fleur/propcalc/local_d_spin_analysis.F90
${CMAKE_SOURCE_DIR}/src/fleur/propcalc/local_structural_frame.F90
${CMAKE_SOURCE_DIR}/src/fleur/propcalc/local_d_spin_frame_transform.F90
${CMAKE_SOURCE_DIR}/src/fleur/propcalc/local_coordination_core.F90
${CMAKE_SOURCE_DIR}/src/fleur/propcalc/local_coordination.F90
)

set(local_d_spin_projector_test_files
${CMAKE_SOURCE_DIR}/src/fleur/math/intgr.F90
${CMAKE_SOURCE_DIR}/src/fleur/propcalc/local_d_spin_projector_core.F90
${CMAKE_SOURCE_DIR}/testing/small-testcodes/local_d_spin_projector_test.F90
)

add_executable(local_d_spin_projector_test EXCLUDE_FROM_ALL ${local_d_spin_projector_test_files})
target_link_libraries(local_d_spin_projector_test juDFT ${FLEUR_LIBRARIES} ${FLEUR_LINK_LIBRARIES})
set_target_properties(local_d_spin_projector_test PROPERTIES
   Fortran_MODULE_DIRECTORY modules/local_d_spin_projector_test
   COMPILE_OPTIONS -Imodules/local_d_spin_projector_test
   RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/tools")
target_compile_options(local_d_spin_projector_test PUBLIC $<$<COMPILE_LANGUAGE:Fortran>:${FLEUR_PRECISION_OPTION}>)
target_compile_options(local_d_spin_projector_test PRIVATE $<$<COMPILE_LANGUAGE:Fortran>:${FLEUR_COMPILE_OPTIONS}>)

set(local_d_spin_analysis_test_files
${CMAKE_SOURCE_DIR}/src/fleur/math/intgr.F90
${CMAKE_SOURCE_DIR}/src/fleur/propcalc/local_d_spin_projector_core.F90
${CMAKE_SOURCE_DIR}/src/fleur/propcalc/local_d_spin_analysis.F90
${CMAKE_SOURCE_DIR}/testing/small-testcodes/local_d_spin_analysis_test.F90
)

add_executable(local_d_spin_analysis_test EXCLUDE_FROM_ALL ${local_d_spin_analysis_test_files})
target_link_libraries(local_d_spin_analysis_test juDFT ${FLEUR_LIBRARIES} ${FLEUR_LINK_LIBRARIES})
set_target_properties(local_d_spin_analysis_test PROPERTIES
   Fortran_MODULE_DIRECTORY modules/local_d_spin_analysis_test
   COMPILE_OPTIONS -Imodules/local_d_spin_analysis_test
   RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/tools")
target_compile_options(local_d_spin_analysis_test PUBLIC $<$<COMPILE_LANGUAGE:Fortran>:${FLEUR_PRECISION_OPTION}>)
target_compile_options(local_d_spin_analysis_test PRIVATE $<$<COMPILE_LANGUAGE:Fortran>:${FLEUR_COMPILE_OPTIONS}>)

set(local_d_spin_frame_test_files
${CMAKE_SOURCE_DIR}/src/fleur/math/intgr.F90
${CMAKE_SOURCE_DIR}/src/fleur/propcalc/local_d_spin_projector_core.F90
${CMAKE_SOURCE_DIR}/src/fleur/propcalc/local_d_spin_analysis.F90
${CMAKE_SOURCE_DIR}/src/fleur/propcalc/local_structural_frame.F90
${CMAKE_SOURCE_DIR}/src/fleur/propcalc/local_d_spin_frame_transform.F90
${CMAKE_SOURCE_DIR}/testing/small-testcodes/local_d_spin_frame_test.F90
)

add_executable(local_d_spin_frame_test EXCLUDE_FROM_ALL ${local_d_spin_frame_test_files})
target_link_libraries(local_d_spin_frame_test juDFT ${FLEUR_LIBRARIES} ${FLEUR_LINK_LIBRARIES})
set_target_properties(local_d_spin_frame_test PROPERTIES
   Fortran_MODULE_DIRECTORY modules/local_d_spin_frame_test
   COMPILE_OPTIONS -Imodules/local_d_spin_frame_test
   RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/tools")
target_compile_options(local_d_spin_frame_test PUBLIC $<$<COMPILE_LANGUAGE:Fortran>:${FLEUR_PRECISION_OPTION}>)
target_compile_options(local_d_spin_frame_test PRIVATE $<$<COMPILE_LANGUAGE:Fortran>:${FLEUR_COMPILE_OPTIONS}>)

set(local_coordination_test_files
${CMAKE_SOURCE_DIR}/src/fleur/propcalc/local_coordination_core.F90
${CMAKE_SOURCE_DIR}/src/fleur/propcalc/local_coordination.F90
${CMAKE_SOURCE_DIR}/src/fleur/propcalc/local_structural_frame.F90
${CMAKE_SOURCE_DIR}/testing/small-testcodes/local_coordination_test.F90
)

add_executable(local_coordination_test EXCLUDE_FROM_ALL ${local_coordination_test_files})
target_link_libraries(local_coordination_test fleurinput juDFT ${FLEUR_LIBRARIES} ${FLEUR_LINK_LIBRARIES})
set_target_properties(local_coordination_test PROPERTIES
   Fortran_MODULE_DIRECTORY modules/local_coordination_test
   COMPILE_OPTIONS -Imodules/local_coordination_test
   RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/tools")
target_compile_options(local_coordination_test PUBLIC $<$<COMPILE_LANGUAGE:Fortran>:${FLEUR_PRECISION_OPTION}>)
target_compile_options(local_coordination_test PRIVATE $<$<COMPILE_LANGUAGE:Fortran>:${FLEUR_COMPILE_OPTIONS}>)

set(spin_orbital_multipoles_test_files
${CMAKE_SOURCE_DIR}/src/fleur/propcalc/spin_orbital_multipoles.F90
${CMAKE_SOURCE_DIR}/src/fleur/propcalc/local_d_spin_multipoles.F90
${CMAKE_SOURCE_DIR}/testing/small-testcodes/spin_orbital_multipoles_test.F90
)

add_executable(spin_orbital_multipoles_test EXCLUDE_FROM_ALL ${spin_orbital_multipoles_test_files})
target_link_libraries(spin_orbital_multipoles_test juDFT ${FLEUR_LIBRARIES} ${FLEUR_LINK_LIBRARIES})
set_target_properties(spin_orbital_multipoles_test PROPERTIES
   Fortran_MODULE_DIRECTORY modules/spin_orbital_multipoles_test
   COMPILE_OPTIONS -Imodules/spin_orbital_multipoles_test
   RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/tools")
target_compile_options(spin_orbital_multipoles_test PUBLIC $<$<COMPILE_LANGUAGE:Fortran>:${FLEUR_PRECISION_OPTION}>)
target_compile_options(spin_orbital_multipoles_test PRIVATE $<$<COMPILE_LANGUAGE:Fortran>:${FLEUR_COMPILE_OPTIONS}>)

set(fleur_F77 ${fleur_F77}
${CMAKE_SOURCE_DIR}/src/fleur/propcalc/dr2fdr.f 
${CMAKE_SOURCE_DIR}/src/fleur/propcalc/nabla.f
)
