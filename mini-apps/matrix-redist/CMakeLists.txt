cmake_minimum_required(VERSION 3.12)
project(FLEUR-mini1 LANGUAGES C Fortran)

set(FLEUR ../../)

#We need juDFT
set(JUDFT_USE_MPI TRUE CACHE BOOL "Compile with MPI, will also work in serial")
add_subdirectory(${FLEUR}/juDFT juDFT)
include_directories("${CMAKE_CURRENT_BINARY_DIR}/juDFT/modules/juDFT")

add_executable(fleur_mini1 main.f90 ${FLEUR}/types/types_mat.F90 ${FLEUR}/types/types_mpimat.F90)
target_link_libraries(fleur_mini1 juDFT)
