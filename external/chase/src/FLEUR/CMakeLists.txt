# -*- Mode: cmake -*-

add_library( "chase_fleur" "./fleur_chase.cpp" )
target_link_libraries("chase_fleur" PUBLIC
  chase_mpi
  )

#install( TARGETS chase_fleur
#  EXPORT chase-fleur-target
#  LIBRARY DESTINATION lib
#  INCLUDES DESTINATION include
#  ARCHIVE DESTINATION lib
#  )

#install(EXPORT chase-fleur-target
#  NAMESPACE ChASE::
#  FILE chase-fleur-config.cmake
#  EXPORT_LINK_INTERFACE_LIBRARIES
#    DESTINATION share)
