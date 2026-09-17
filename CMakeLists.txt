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

include("cmake/CompilerConfig.txt")

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
