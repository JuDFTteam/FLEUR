# -*- Mode: cmake -*-

include(CheckLanguage)
check_language(CUDA)

add_library( "chase_fleur" "./fleur_chase.cpp" )

if(CMAKE_CUDA_COMPILER)
	message("Use distributed multi-GPU version of ChASE")
	target_compile_definitions("chase_fleur" PRIVATE CHASE_USE_GPU)
        target_link_libraries("chase_fleur" PUBLIC
                chase_mpi chase_cuda
                )
else()
  	target_link_libraries("chase_fleur" PUBLIC
    		chase_mpi
    		)
endif()


install( TARGETS chase_fleur
  EXPORT chase-fleur
  LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
  INCLUDES DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}
  ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
  )

#install(EXPORT chase-fleur-target
#  NAMESPACE ChASE::
#  FILE chase-fleur-config.cmake
#  EXPORT_LINK_INTERFACE_LIBRARIES
#    DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/share)
