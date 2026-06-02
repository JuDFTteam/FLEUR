set(FLEUR_USE_GITVERSION FALSE)
set(FLEUR_GIT_FAILURE_REASON "")
find_package(Git)
if (Git_FOUND)
    execute_process(
        COMMAND ${GIT_EXECUTABLE} submodule update
        WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}
        RESULT_VARIABLE _res_update
        OUTPUT_VARIABLE _git_update_stdout
        ERROR_VARIABLE _git_update_stderr
    )
    if ( _res_update EQUAL 0)
        set(FLEUR_USE_GITVERSION TRUE)
    else()
        string(STRIP "${_git_update_stderr}" _git_update_stderr)
        string(STRIP "${_git_update_stdout}" _git_update_stdout)
        if (_git_update_stderr)
            set(FLEUR_GIT_FAILURE_REASON "git submodule update failed: ${_git_update_stderr}")
        elseif (_git_update_stdout)
            set(FLEUR_GIT_FAILURE_REASON "git submodule update failed: ${_git_update_stdout}")
        else()
            set(FLEUR_GIT_FAILURE_REASON "git submodule update failed with exit code ${_res_update}")
        endif()
        set(FLEUR_WARN_MESSAGE "${FLEUR_WARN_MESSAGE}\n${FLEUR_GIT_FAILURE_REASON}")
    endif()
else()
    set(FLEUR_GIT_FAILURE_REASON "Git executable not found")
    set(FLEUR_WARN_MESSAGE "${FLEUR_WARN_MESSAGE}\n${FLEUR_GIT_FAILURE_REASON}")
endif()    