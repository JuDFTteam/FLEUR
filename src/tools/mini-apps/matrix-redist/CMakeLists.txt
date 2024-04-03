cmake_minimum_required(VERSION 3.12)
project(FLEUR-mini1 LANGUAGES C Fortran)

set(FLEUR ${CMAKE_SOURCE_DIR}/../../..)

#We need juDFT
set(JUDFT_USE_MPI TRUE CACHE BOOL "Compile with MPI, will also work in serial")
add_subdirectory(${FLEUR}/libraries/juDFT juDFT)
include_directories("${CMAKE_CURRENT_BINARY_DIR}/modules/juDFT")

add_executable(fleur_mini1 main.f90 constants.f90 
            ${FLEUR}/fleur/types/types_mat.F90 ${FLEUR}/fleur/types/types_mpimat.F90 ${FLEUR}/fleur/hybrid/divide_most_evenly.f90)
target_compile_definitions(fleur_mini1 PUBLIC "CPP_MPI;CPP_SCALAPACK")
target_link_libraries(fleur_mini1 juDFT)
