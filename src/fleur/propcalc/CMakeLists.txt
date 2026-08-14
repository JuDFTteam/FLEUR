set(fleur_F90 ${fleur_F90}
${CMAKE_SOURCE_DIR}/src/fleur/propcalc/outcdn.f90
${CMAKE_SOURCE_DIR}/src/fleur/propcalc/plot.F90
${CMAKE_SOURCE_DIR}/src/fleur/propcalc/mcd_init.f90
${CMAKE_SOURCE_DIR}/src/fleur/propcalc/local_d_spin_projector_core.F90
${CMAKE_SOURCE_DIR}/src/fleur/propcalc/local_d_spin_projector.F90
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

set(fleur_F77 ${fleur_F77}
${CMAKE_SOURCE_DIR}/src/fleur/propcalc/dr2fdr.f 
${CMAKE_SOURCE_DIR}/src/fleur/propcalc/nabla.f
)
