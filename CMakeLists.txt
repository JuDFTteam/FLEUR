cmake_minimum_required(VERSION 2.8)
project(FLEUR LANGUAGES NONE)

include("cmake/Architectures.txt")
include("cmake/ReportConfig.txt")

#Here the targets and the files are defined
include("cmake/Files_and_Targets.txt")

#install(TARGETS fleur inpgen DESTINATION bin)
install(PROGRAMS ${CMAKE_BINARY_DIR}/fleur 
        CONFIGURATIONS Debug
        DESTINATION bin
        RENAME fleur_dbg)
