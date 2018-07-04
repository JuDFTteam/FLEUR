cmake_minimum_required(VERSION 3.0)
#Some stuff might be defined in config file
if (EXISTS "${CMAKE_BINARY_DIR}/config.cmake")
   include("${CMAKE_BINARY_DIR}/config.cmake")
endif()

project(FLEUR LANGUAGES C Fortran)



#some variables might be set in the environment
set(FLEUR_LIBRARIES ${FLEUR_LIBRARIES} $ENV{FLEUR_LIBRARIES})
set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} $ENV{FLEUR_Fortran_FLAGS} $ENV{CMAKE_Fortran_FLAGS}")
message(STATUS "FLEUR_Fortran_FLAGS: ${FLEUR_Fortran_FLAGS}")

if (DEFINED ENV{FLEUR_USE_SERIAL})
	set(FLEUR_USE_SERIAL ENV{FLEUR_USE_SERIAL})
else()
	set(FLEUR_USE_SERIAL TRUE)
endif()


include("cmake/CompilerConfig.txt")

include("cmake/ReportConfig.txt")

include("cmake/Files_and_Targets.txt")

include("cmake/filespecific.cmake")

#install(TARGETS fleur inpgen DESTINATION bin)
install(PROGRAMS ${CMAKE_BINARY_DIR}/fleur 
        CONFIGURATIONS Debug
        DESTINATION bin
        RENAME fleur_dbg)
