#First check if we can compile with LAPACK
try_compile(FLEUR_USE_LAPACK ${CMAKE_BINARY_DIR} ${CMAKE_SOURCE_DIR}/cmake/tests/test_LAPACK.f90
LINK_LIBRARIES ${FLEUR_LIBRARIES}
            )

if (NOT FLEUR_USE_LAPACK)
   find_package(LAPACK)
   set(LAPACK_LIBRARIES_dedup ${LAPACK_LIBRARIES})
   list(REMOVE_DUPLICATES LAPACK_LIBRARIES_dedup)
   if ("$ENV{VERBOSE}")
            	message("LAPACK_LIBRARIES: ${LAPACK_LIBRARIES}")
    endif()
   foreach(TEST_LIBRARIES "${LAPACK_LIBRARIES_dedup}" "${LAPACK_LIBRARIES}" "-lopenblas")
      if (NOT FLEUR_USE_LAPACK)
         try_compile(FLEUR_USE_LAPACK ${CMAKE_BINARY_DIR} ${CMAKE_SOURCE_DIR}/cmake/tests/test_LAPACK.f90
                     LINK_LIBRARIES ${FLEUR_LIBRARIES} ${TEST_LIBRARIES} OUTPUT_VARIABLE compile_output 
		     
         )
	 if ("$ENV{VERBOSE}")
            	message("LAPACK compile test: ${FLEUR_USE_LAPACK}\nLINK_LIBRARIES ${TEST_LIBRARIES}\n${compile_output}")
     	    endif()
         if (FLEUR_USE_LAPACK)
            set(FLEUR_LIBRARIES ${FLEUR_LIBRARIES} ${TEST_LIBRARIES})
         endif()
      endif()
   endforeach()
endif()

message("LAPACK Library found:${FLEUR_USE_LAPACK}")

# Verify the complex-return BLAS calling convention (zdotc/zdotu/cdotc/cdotu).
# Apple's Accelerate framework returns these via the f2c/g77 convention, which
# mismatches gfortran's native convention and yields silently wrong results.
# We can only check this by actually running a tiny program, so skip it when
# cross-compiling.
if (FLEUR_USE_LAPACK AND NOT CMAKE_CROSSCOMPILING)
   try_run(ZDOTC_RUN_RESULT ZDOTC_COMPILE_RESULT ${CMAKE_BINARY_DIR}
           ${CMAKE_SOURCE_DIR}/cmake/tests/test_zdotc.f90
           LINK_LIBRARIES ${FLEUR_LIBRARIES})
   if (ZDOTC_COMPILE_RESULT AND NOT ZDOTC_RUN_RESULT EQUAL 0)
      message(WARNING
         "The linked BLAS returns WRONG results for complex-valued functions "
         "(zdotc/zdotu/cdotc/cdotu): it uses the f2c/g77 calling convention, "
         "which is incompatible with this Fortran compiler. This is the "
         "well-known Apple Accelerate framework issue on macOS.\n"
         "FLEUR's own sources avoid these functions, but any other library you "
         "link against them may misbehave. The recommended fix is to link "
         "OpenBLAS instead, e.g.:\n"
         "  brew install openblas\n"
         "  -DLAPACK_LIBRARIES=<openblas-prefix>/lib/libopenblas.dylib\n"
         "or set -DBLA_VENDOR=OpenBLAS before configuring.")
   endif()
endif()
