#this file sets some preprocessor variables that are used in
#init/compile_descr.F90 to determine the program version and
#some compilation environment description

cmake_host_system_information(RESULT compile_host QUERY HOSTNAME)
set(compile_user $ENV{USER})
string(TIMESTAMP compile_time)
set(git_hash unknown)
set(git_describe unknown)
set(git_branch unknown)
if (EXISTS "${CMAKE_SOURCE_DIR}/.git")
   execute_process(COMMAND git describe --tags --match "MaX*" --dirty
                   WORKING_DIRECTORY ${CMAKE_SOURCE_DIR} OUTPUT_VARIABLE git_describe ERROR_QUIET)
   execute_process(COMMAND git rev-parse --abbrev-ref HEAD
                   WORKING_DIRECTORY ${CMAKE_SOURCE_DIR} OUTPUT_VARIABLE git_branch ERROR_QUIET)
   execute_process(COMMAND git rev-parse HEAD
                   WORKING_DIRECTORY ${CMAKE_SOURCE_DIR} OUTPUT_VARIABLE git_hash ERROR_QUIET)
elseif (EXISTS "${CMAKE_SOURCE_DIR}/version")
   file(READ ${CMAKE_SOURCE_DIR}/version git_describe)
endif()

#normalize the strings, fix for problems in git commands above
if (git_hash)
   string(STRIP ${git_hash} git_hash)
else()
   set(git_hash unknown)
endif()
if (git_describe)
   string(STRIP ${git_describe} git_describe)
else()
   set(git_describe unknown)
endif()
if (git_branch)
   string(STRIP ${git_branch} git_branch)
else()
   set(git_branch unknown)
endif()

file(WRITE "${BI_FILE}" "compile_date=\"${compile_time}\"\ncompile_user=\"${compile_user}\"\ncompile_host=\"${compile_host}\"\ngitdesc=\"${git_describe}\"\ngitbranch=\"${git_branch}\"\ngithash=\"${git_hash}\"")
