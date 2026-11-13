#Check that the Fortran compiler can build a trivial program with the given flags.
#project() tests the compiler without CMAKE_Fortran_FLAGS (see CMakeLists.txt), so
#broken flags (e.g. FFLAGS=-lbla) would otherwise only show up as failing library tests.
try_compile(FLEUR_COMPILER_OK ${CMAKE_BINARY_DIR} ${CMAKE_SOURCE_DIR}/cmake/tests/test_compiler.f90
            OUTPUT_VARIABLE compile_output)

if (FLEUR_COMPILER_OK)
   message("Fortran compiler works with the given flags")
else()
   message(FATAL_ERROR "\nThe Fortran compiler fails to compile and link a trivial program.\n"
                       "Compiler: ${CMAKE_Fortran_COMPILER}\n"
                       "Flags:    ${CMAKE_Fortran_FLAGS}\n"
                       "Check the flags set via FFLAGS, the -flags option of configure.sh or config.cmake.\n"
                       "Compiler output:\n${compile_output}")
endif()
