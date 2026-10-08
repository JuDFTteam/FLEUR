cmake_minimum_required(VERSION 3.5)
#Some stuff might be defined in config file
if (EXISTS "${CMAKE_BINARY_DIR}/config.cmake")
   include("${CMAKE_BINARY_DIR}/config.cmake")
endif()

if(${CMAKE_VERSION} VERSION_LESS "3.10.0")
    message("Your cmake Version is rather old. Please consider to update cmake")
    message("More modern cmake versions might be needed to make sure your compiler settings are detected correctly")
endif()

# sometimes cmake clears CMAKE_Fortran_FLAGS during project()
set(CMAKE_Fortran_FLAGS_backup ${CMAKE_Fortran_FLAGS})
set(CMAKE_Fortran_FLAGS "")
project(FLEUR LANGUAGES C CXX Fortran)
string(STRIP "${CMAKE_Fortran_FLAGS} ${CMAKE_Fortran_FLAGS_backup}" CMAKE_Fortran_FLAGS)

# CMake >= 4.0 rejects cmake_minimum_required(VERSION <3.5) in bundled
# submodules (libxc, BLACS test project of SCALAPACK, ...). The environment
# variable also covers sub-configures started via execute_process.
if (NOT DEFINED CMAKE_POLICY_VERSION_MINIMUM)
   set(CMAKE_POLICY_VERSION_MINIMUM 3.5)
endif()
if (NOT DEFINED ENV{CMAKE_POLICY_VERSION_MINIMUM})
   set(ENV{CMAKE_POLICY_VERSION_MINIMUM} ${CMAKE_POLICY_VERSION_MINIMUM})
endif()

include("cmake/CompilerConfig.txt")

include("cmake/rpath.cmake")

include("cmake/Generate_Schema.cmake")

include("cmake/Files_and_Targets.txt")

include("cmake/filespecific.cmake")
include("cmake/pytest_setup.cmake")

include("cmake/ReportConfig.txt")

#The target is renamed on output ('fleur' -> 'fleur_MPI' for an MPI build);
#install(TARGETS) honours RUNTIME_OUTPUT_NAME.
install(TARGETS fleur   RUNTIME DESTINATION bin)
