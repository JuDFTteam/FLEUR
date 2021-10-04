cmake_minimum_required(VERSION 3.12)
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
include("cmake/ReportConfig.txt")

#install(TARGETS fleur inpgen DESTINATION bin)
install(PROGRAMS ${CMAKE_BINARY_DIR}/fleur DESTINATION bin OPTIONAL)
install(PROGRAMS ${CMAKE_BINARY_DIR}/fleur_MPI DESTINATION bin OPTIONAL)
install(PROGRAMS ${CMAKE_BINARY_DIR}/inpgen DESTINATION bin)