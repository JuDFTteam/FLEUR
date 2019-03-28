if(NOT EXISTS "${PROJECT_SOURCE_DIR}/juDFT/CMakeLists.txt" )
    find_package(Git REQUIRED)
    execute_process(COMMAND ${GIT_EXECUTABLE} submodule init juDFT WORKING_DIRECTORY ${PROJECT_SOURCE_DIR} RESULT_VARIABLE _res_init OUTPUT_QUIET ERROR_QUIET)
    execute_process(COMMAND ${GIT_EXECUTABLE} submodule update  WORKING_DIRECTORY ${PROJECT_SOURCE_DIR} RESULT_VARIABLE _res_update OUTPUT_QUIET ERROR_QUIET)
    if( ${_res_init} GREATER 0 OR ${_res_update} GREATER 0 )
       message(FATAL_ERROR "HDF5 source could not be downloaded.\n"
                     "We tried: 'git submodule init external/libxc-git && git submodule update' and resulted in error" )
    endif()
endif()
set(JUDFT_USE_MPI ${FLEUR_USE_MPI} CACHE BOOL "Compile with MPI, will also work in serial")
set(JUDFT_USE_HDF5 ${FLEUR_USE_HDF5} CACHE BOOL "Compile with HDF5")
if (DEFINED FLEUR_USE_HDF5MPI)
set(JUDFT_USE_HDF5MPI FLEUR_USE_HDF5MPI CACHE BOOL "Is the HDF5 version able to do parallel IO" )
endif()
#In addition you might want to set
set(JUDFT_COMPILEOPTS ${FLEUR_PRECISION_OPTION})

add_subdirectory(juDFT)

include_directories("${CMAKE_CURRENT_BINARY_DIR}/juDFT/modules/juDFT")
