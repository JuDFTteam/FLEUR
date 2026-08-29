set(fleur_F90 ${fleur_F90}
${CMAKE_SOURCE_DIR}/src/fleur/propcalc/xas/xas_angular.F90
${CMAKE_SOURCE_DIR}/src/fleur/propcalc/xas/xas_amplitudes.F90
${CMAKE_SOURCE_DIR}/src/fleur/propcalc/xas/xas_core.F90
${CMAKE_SOURCE_DIR}/src/fleur/propcalc/xas/xas_transitions.F90
${CMAKE_SOURCE_DIR}/src/fleur/propcalc/xas/rixs_spectrum.F90
${CMAKE_SOURCE_DIR}/src/fleur/propcalc/xas/rixs_io.F90
${CMAKE_SOURCE_DIR}/src/fleur/propcalc/xas/rixs_state_character_core.F90
${CMAKE_SOURCE_DIR}/src/fleur/propcalc/xas/rixs_state_character.F90
${CMAKE_SOURCE_DIR}/src/fleur/propcalc/xas/rixs_driver.F90
${CMAKE_SOURCE_DIR}/src/fleur/propcalc/xas/xas_driver.F90
${CMAKE_SOURCE_DIR}/src/fleur/propcalc/xas/xas_io.F90
${CMAKE_SOURCE_DIR}/src/fleur/propcalc/xas/xas_matrixelements.F90
${CMAKE_SOURCE_DIR}/src/fleur/propcalc/xas/xas_radial.F90
${CMAKE_SOURCE_DIR}/src/fleur/propcalc/xas/xas_spectrum.F90
${CMAKE_SOURCE_DIR}/src/fleur/propcalc/xas/xas_symmetry.F90
)

if(FLEUR_USE_HDF5)
   set(rixs_state_character_test_files
      ${CMAKE_SOURCE_DIR}/src/fleur/propcalc/xas/rixs_state_character_core.F90
      ${CMAKE_SOURCE_DIR}/testing/small-testcodes/rixs_state_character_test.F90)
   add_executable(rixs_state_character_test EXCLUDE_FROM_ALL ${rixs_state_character_test_files})
   target_link_libraries(rixs_state_character_test juDFT ${FLEUR_LIBRARIES} ${FLEUR_LINK_LIBRARIES})
   target_compile_definitions(rixs_state_character_test PUBLIC ${FLEUR_DEFINITIONS})
   set_target_properties(rixs_state_character_test PROPERTIES
      Fortran_MODULE_DIRECTORY modules/rixs_state_character_test
      COMPILE_OPTIONS -Imodules/rixs_state_character_test
      RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/tools")
   target_compile_options(rixs_state_character_test PUBLIC $<$<COMPILE_LANGUAGE:Fortran>:${FLEUR_PRECISION_OPTION}>)
   target_compile_options(rixs_state_character_test PRIVATE $<$<COMPILE_LANGUAGE:Fortran>:${FLEUR_COMPILE_OPTIONS}>)
endif()
