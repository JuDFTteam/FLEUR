cmake_minimum_required(VERSION 3.0)
#Some stuff might be defined in config file
if (EXISTS "${CMAKE_BINARY_DIR}/config.cmake")
   include("${CMAKE_BINARY_DIR}/config.cmake")
endif()

project(FLEUR LANGUAGES C Fortran)

include("cmake/CompilerConfig.txt")

include("cmake/ReportConfig.txt")
include("cmake/Generate_Schema.cmake")
include("cmake/Files_and_Targets.txt")

include("cmake/filespecific.cmake")

#install(TARGETS fleur inpgen DESTINATION bin)
install(PROGRAMS ${CMAKE_BINARY_DIR}/fleur 
        CONFIGURATIONS Debug
        DESTINATION bin
        RENAME fleur_dbg)
