cmake_minimum_required(VERSION 3.0)
#Some stuff might be defined in config file
if (EXISTS "${CMAKE_BINARY_DIR}/config.cmake")
   include("${CMAKE_BINARY_DIR}/config.cmake")
endif()

set(tmp ${CMAKE_Fortran_FLAGS})
project(FLEUR LANGUAGES C Fortran)
set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} ${tmp}")

include("cmake/CompilerConfig.txt")

include("cmake/Generate_Schema.cmake")

include("cmake/Files_and_Targets.txt")

include("cmake/filespecific.cmake")

include("cmake/ReportConfig.txt")




#install(TARGETS fleur inpgen DESTINATION bin)
install(PROGRAMS ${CMAKE_BINARY_DIR}/fleur 
        CONFIGURATIONS Debug
        DESTINATION bin
        RENAME fleur_dbg)
