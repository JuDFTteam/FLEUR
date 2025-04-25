cmake_minimum_required(VERSION 3.12)
project(FLEUR-mini2 LANGUAGES C Fortran)

set(FLEUR ${CMAKE_SOURCE_DIR}/../../..)
set(CMAKE_TESTS ${CMAKE_SOURCE_DIR}/../../../../cmake/tests/)

include(${CMAKE_TESTS}test_precision.cmake)
include(${CMAKE_TESTS}test_MPI.cmake)
include(${CMAKE_TESTS}test_SCALAPACK.cmake)
include(${CMAKE_TESTS}test_OpenMP.cmake)

include(${CMAKE_TESTS}test_ELPA.cmake)
include(${CMAKE_TESTS}test_ELSI.cmake)


#We need juDFT
add_subdirectory(${FLEUR}/libraries/juDFT juDFT)
include_directories("${CMAKE_CURRENT_BINARY_DIR}/modules/juDFT")

#We compile stuff in "diagnalization"
set(diag_files
${FLEUR}/fleur/diagonalization/eigen_diag.F90
${FLEUR}/fleur/diagonalization/lapack_diag.F90
${FLEUR}/fleur/diagonalization/magma.F90
${FLEUR}/fleur/diagonalization/scalapack.F90
${FLEUR}/fleur/diagonalization/cusolver_diag.F90
${FLEUR}/fleur/diagonalization/available_solvers.F90 
#${FLEUR}/fleur/diagonalization/elpa_20180525.F90
${FLEUR}/fleur/diagonalization/elsi.F90
${FLEUR}/fleur/diagonalization/magma.F90
${FLEUR}/fleur/diagonalization/writeout.F90
${FLEUR}/fleur/io/io_matrix.F90
${FLEUR}/fleur/diagonalization/lapack_singlePrec_diag.F90
${FLEUR}/fleur/diagonalization/dummy_diag.F90
${FLEUR}/fleur/diagonalization/chase_diag.F90
${FLEUR}/fleur/diagonalization/elemental.F90
${FLEUR}/fleur/diagonalization/elpa_20180525_onenode.F90
)

if (FLEUR_USE_ELPA_20180525)
 set(diag_files ${diag_files} ${FLEUR}/fleur/diagonalization/elpa_20180525.F90)
else()
  set(diag_files ${diag_files} ${FLEUR}/fleur/diagonalization/elpa.F90)
endif()

add_executable(fleur_mini2 main.F90 constants.f90 ${diag_files}
            ${FLEUR}/fleur/types/types_mat.F90 ${FLEUR}/fleur/types/types_mpimat.F90 ${FLEUR}/fleur/hybrid/divide_most_evenly.f90)
target_link_libraries(fleur_mini2 juDFT)
target_compile_options(fleur_mini2 PUBLIC $<$<COMPILE_LANGUAGE:Fortran>:${FLEUR_PRECISION_OPTION}>)
target_compile_options(fleur_mini2 PUBLIC $<$<COMPILE_LANGUAGE:Fortran>:${FLEUR_COMPILE_OPTIONS}>)
target_compile_definitions(fleur_mini2 PUBLIC ${FLEUR_MPI_DEFINITIONS})
target_link_libraries(fleur_mini2 ${FLEUR_LINK_LIBRARIES})
target_link_libraries(fleur_mini2 ${FLEUR_LIBRARIES})
