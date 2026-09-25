set(fleur_F77 ${fleur_F77}
src/fleur/math/angle.f
src/fleur/math/clebsch.f
src/fleur/math/dsphbs.f
src/fleur/math/gaussp.f
src/fleur/math/inwint.f
src/fleur/math/outint.f
src/fleur/math/points.f
src/fleur/math/qsf.f
src/fleur/math/sphbes.f
src/fleur/math/sphpts.f
src/fleur/math/rcerf.f
src/fleur/math/diflgr.f
)
set(fleur_F90 ${fleur_F90}
src/fleur/math/inv3.f90
src/fleur/math/qranf.f90
src/fleur/math/differentiate.f90
src/fleur/math/SphBessel.f90
src/fleur/math/DoubleFactorial.f90
src/fleur/math/LattHarmsSphHarmsConv.f90
src/fleur/math/ExpSave.f90
src/fleur/math/intgr.F90
src/fleur/math/trapz.f90
src/fleur/math/lagrange3.f90
src/fleur/math/ylm4.F90
src/fleur/math/gradYlm.f90
src/fleur/math/pol_angle.f90
src/fleur/math/BfieldtoVmat.f90
src/fleur/math/matmul_dgemm.F90
)

#Regression test for the Clebsch-Gordan selection rule, see issue #804.
#Built on demand; testing/tests/invariants/test_math_invariants.py runs it.
set(clebsch_test_files
${CMAKE_SOURCE_DIR}/src/fleur/math/clebsch.f
${CMAKE_SOURCE_DIR}/testing/small-testcodes/clebsch_test.f90
)
add_executable(clebsch_test EXCLUDE_FROM_ALL ${clebsch_test_files})
set_target_properties(clebsch_test PROPERTIES
   Fortran_MODULE_DIRECTORY modules/clebsch_test
   RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/tools")
target_compile_options(clebsch_test PUBLIC $<$<COMPILE_LANGUAGE:Fortran>:${FLEUR_PRECISION_OPTION}>)
target_compile_options(clebsch_test PRIVATE -Imodules/clebsch_test)
