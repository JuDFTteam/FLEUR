#Try if the linker gives correct error messages
try_compile(THIS_SHOULD_FAIL ${CMAKE_BINARY_DIR} ${CMAKE_SOURCE_DIR}/cmake/tests/test_linker.f90
	    LINK_LIBRARIES ${FLEUR_LIBRARIES}
            )

if ( NOT THIS_SHOULD_FAIL )
   message("Linker seems OK")
else()
  message(WARNING "Your linker is broken and reports no error even if dependencies are not available.\n" "You might want to complain with your system admin.")
  if ( CLI_WARN_ONLY )
     message(WARNING "You choose to ignore the warning, your final linking may fail.")
  else()
     message(FATAL_ERROR "\n\nYou can use the -warn_only option to ignore this error." "If you do so you might fail at final linking. :-)\n")		      
  endif()
endif()