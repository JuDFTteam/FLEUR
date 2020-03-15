#Compile the kplib-libray as external dependency
if (DEFINED CLI_FLEUR_USE_KPLIB)
   #first get the source code if required
   if(NOT EXISTS "${PROJECT_SOURCE_DIR}/external/kplib/kplib/src" )
    find_package(Git REQUIRED)
    execute_process(COMMAND ${GIT_EXECUTABLE} submodule init external/kplib/kplib WORKING_DIRECTORY ${PROJECT_SOURCE_DIR} RESULT_VARIABLE _res_init OUTPUT_QUIET ERROR_QUIET)
    execute_process(COMMAND ${GIT_EXECUTABLE} submodule update external/kplib/kplib WORKING_DIRECTORY ${PROJECT_SOURCE_DIR} RESULT_VARIABLE _res_update OUTPUT_QUIET ERROR_QUIET)
    if( ${_res_init} GREATER 0 OR ${_res_update} GREATER 0 )
        message(FATAL_ERROR "kplib source could not be downloaded.\n"
            "We tried: 'git submodule init external/kplib/kplibt && git submodule update' and resulted in error" )
    endif()
   endif()
   #Now add library
   add_subdirectory("${PROJECT_SOURCE_DIR}/external/kplib/kplib")
   set(INPGEN_USE_kplib ON CACHE BOOL "Use kplib library")
   message("Using the kplib interface")
endif()