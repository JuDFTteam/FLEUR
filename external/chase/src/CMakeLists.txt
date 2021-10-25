# -*- Mode: cmake -*-
cmake_minimum_required( VERSION 3.8 )

project( CHASE-DRIVERS LANGUAGES C CXX )

set(CMAKE_INSTALL_PREFIX ${CMAKE_CURRENT_BINARY_DIR}/install)

add_subdirectory("./ChASE" "chase_interface")

add_definitions( "-DCHASE_OUTPUT" )

add_subdirectory("./FLEUR")
