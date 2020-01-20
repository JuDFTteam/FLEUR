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

add_subdirectory("inpgen2")


#install(TARGETS fleur inpgen DESTINATION bin)
install(PROGRAMS ${CMAKE_BINARY_DIR}/fleur
        CONFIGURATIONS Debug
        DESTINATION bin
        RENAME fleur_dbg)

find_program(XXD_PROG xxd)

if (XXD_PROG)
    ADD_CUSTOM_COMMAND(
    OUTPUT ${CMAKE_SOURCE_DIR}/inpgen2/default_econfig.h
    COMMAND ${XXD_PROG} -i default.econfig default_econfig.h
    DEPENDS ${CMAKE_SOURCE_DIR}/inpgen2/default.econfig
    WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}/inpgen2/
    COMMENT "Putting current default.econfig into default_econfig.h")
else()
    ADD_CUSTOM_COMMAND(
    OUTPUT ${CMAKE_SOURCE_DIR}/inpgen2/default_econfig.h
    COMMAND cp  ${CMAKE_SOURCE_DIR}/inpgen2/default_econfig.h.backup ${CMAKE_SOURCE_DIR}/inpgen2/default_econfig.h
    COMMENT "No xxd found using backup")
    message("No xxd command found! Using backup of default_econfig.h")
endif()
