# -*- Mode: cmake -*-
cmake_minimum_required( VERSION 3.8 )

project( CHASE-DRIVERS LANGUAGES C CXX )

add_subdirectory("./ChASE" "chase_interface")

add_definitions( "-DCHASE_OUTPUT" )

add_subdirectory("./FLEUR")
#add_subdirectory("./ChASE/ChASE-MPI")
