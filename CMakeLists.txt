cmake_minimum_required(VERSION 3.5)
#Some stuff might be defined in config file
if (EXISTS "${CMAKE_BINARY_DIR}/config.cmake")
   include("${CMAKE_BINARY_DIR}/config.cmake")
endif()

if(${CMAKE_VERSION} VERSION_LESS "3.10.0")
    message("Your cmake Version is rather old. Please consider to update cmake")
    message("More modern cmake versions might be needed to make sure your compiler settings are detected correctly")
endif()

# project() sets the cache entry CMAKE_Fortran_FLAGS from $ENV{FFLAGS}, but a normal variable
# of the same name hides it. Flags set before (config.cmake, -D) and the ones from the
# environment are combined, each only once.
set(CMAKE_Fortran_FLAGS_backup "${CMAKE_Fortran_FLAGS}")
set(CMAKE_Fortran_FLAGS "")
project(FLEUR LANGUAGES C CXX Fortran)
get_property(FLEUR_FFLAGS_INIT CACHE CMAKE_Fortran_FLAGS PROPERTY VALUE)
string(STRIP "${FLEUR_FFLAGS_INIT}" FLEUR_FFLAGS_INIT)
set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS_backup}")
if (NOT "${FLEUR_FFLAGS_INIT}" STREQUAL "")
   string(FIND "${CMAKE_Fortran_FLAGS_backup}" "${FLEUR_FFLAGS_INIT}" FLEUR_FFLAGS_POS)
   if (FLEUR_FFLAGS_POS EQUAL -1)
      set(CMAKE_Fortran_FLAGS "${FLEUR_FFLAGS_INIT} ${CMAKE_Fortran_FLAGS_backup}")
   endif()
endif()
string(STRIP "${CMAKE_Fortran_FLAGS}" CMAKE_Fortran_FLAGS)

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

# Also build inpgen
add_subdirectory("src/tools/inpgen2")
# Ensure buildinfo.h is generated before inpgen2 is compiled.
# The 'buildinfo' target (defined in cmake/compileenv.txt) generates the header
# independently, so we no longer need to depend on the full fleur/fleur_MPI builds.
add_dependencies("inpgen2" buildinfo)

include("cmake/ReportConfig.txt")

#Both targets are renamed on output ('fleur' -> 'fleur_MPI' for an MPI build,
#'inpgen2' -> 'inpgen'); install(TARGETS) honours RUNTIME_OUTPUT_NAME.
install(TARGETS fleur   RUNTIME DESTINATION bin)
install(TARGETS inpgen2 RUNTIME DESTINATION bin)
