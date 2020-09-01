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
