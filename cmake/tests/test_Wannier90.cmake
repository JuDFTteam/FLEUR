#first try if Wannier90 already works
try_compile(FLEUR_USE_WANN ${CMAKE_BINARY_DIR} ${CMAKE_SOURCE_DIR}/cmake/tests/test_Wannier90.f90
	    LINK_LIBRARIES ${FLEUR_LIBRARIES}
            )
            
            
message("Wannier90 1.2 Library found:${FLEUR_USE_WANN}")


foreach(ADD_String "-lwannier;-lmkl_intel_lp64;-lmkl_sequential;-lmkl_core" )



   if (NOT FLEUR_USE_WANN)
     set(TEST_LIBRARIES "${FLEUR_LIBRARIES};${ADD_String}")

     message("compilation test:${CMAKE_BINARY_DIR} ${CMAKE_SOURCE_DIR}/cmake/tests/test_Wannier90.f90
            LINK_LIBRARIES ${TEST_LIBRARIES}")


     try_compile(FLEUR_USE_WANN ${CMAKE_BINARY_DIR} ${CMAKE_SOURCE_DIR}/cmake/tests/test_Wannier90.f90
            LINK_LIBRARIES ${TEST_LIBRARIES} OUTPUT_VARIABLE TELL_ME
            )

     if(DEFINED ENV{VERBOSE})
        message("TELL_ME=${TELL_ME}")
     endif()


     if (FLEUR_USE_WANN)
          set(FLEUR_WANNIER90_LIBRARIES ${TEST_LIBRARIES})
          set(  FLEUR_LIBRARIES "${FLEUR_LIBRARIES};${ADD_String}"  )
     endif()
   endif()
endforeach()            
            
            

message("Wannier90 1.2 Library found:${FLEUR_USE_WANN}")

if (DEFINED CLI_FLEUR_USE_WANNIER)
   if (${CLI_FLEUR_USE_WANNIER})
       if (NOT FLEUR_USE_WANN)
           message("You asked for Wannier90 but cmake couldn't find it. We will download the gitlab sources and use it during compilation")
           set(FLEUR_USE_WANN TRUE)
           if(NOT EXISTS "${PROJECT_SOURCE_DIR}/external/wannier90/src" )
             find_package(Git REQUIRED)
    	     execute_process(COMMAND ${GIT_EXECUTABLE} submodule init external/wannier90 WORKING_DIRECTORY ${PROJECT_SOURCE_DIR} RESULT_VARIABLE _res_init OUTPUT_QUIET ERROR_QUIET)
    	     execute_process(COMMAND ${GIT_EXECUTABLE} submodule update  WORKING_DIRECTORY ${PROJECT_SOURCE_DIR} RESULT_VARIABLE _res_update OUTPUT_QUIET ERROR_QUIET)
    	     if( ${_res_init} GREATER 0 OR ${_res_update} GREATER 0 )
               message(FATAL_ERROR "Wannier90 source could not be downloaded.\n"
                        "We tried: 'git submodule init external/wannier90 && git submodule update' and resulted in error" )
             endif()			
          endif()
       endif()
   else()
       if (FLEUR_USE_WANN)
           message("Wannier library found, but you explicitely asked not to use it")
	   set(FLEUR_USE_WANN FALSE)
       endif()
   endif()	   
endif()  

if (FLEUR_USE_WANN)
   set(FLEUR_DEFINITIONS ${FLEUR_DEFINITIONS} "CPP_WANN") 
   set(FLEUR_MPI_DEFINITIONS ${FLEUR_MPI_DEFINITIONS} "CPP_WANN")
endif()
