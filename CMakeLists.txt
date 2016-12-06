cmake_minimum_required(VERSION 2.8)
project(FLEUR LANGUAGES C Fortran)

#some variables might be set in the environment
set(FLEUR_LIBRARIES ${FLEUR_LIBRARIES} $ENV{FLEUR_LIBRARIES})
set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} $ENV{CMAKE_Fortran_FLAGS}")

if (DEFINED ENV{FLEUR_NO_SERIAL})
	set(FLEUR_USE_SERIAL false)
else()
	set(FLEUR_USE_SERIAL true)
endif()


include("cmake/CompilerConfig.txt")

include("cmake/ReportConfig.txt")

include("cmake/Files_and_Targets.txt")

#install(TARGETS fleur inpgen DESTINATION bin)
install(PROGRAMS ${CMAKE_BINARY_DIR}/fleur 
        CONFIGURATIONS Debug
        DESTINATION bin
        RENAME fleur_dbg)
