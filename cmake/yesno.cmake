#macro to ask user for confirmation
#Provide question as argument
#only used if FLEUR_USE_INTERACTIVE is TRUE, otherwise yesno_answer variable is always 'no'
macro(yesno)

if (FLEUR_USE_INTERACTIVE)
find_program(BASH NAMES bash bash.exe
PATHS "c:/msys64/usr/bin" "$ENV{PROGRAMFILES}/Git/bin"
)
# message("#BASH ${BASH}")
message("---------- User confirmation required ----------\n${ARGV}\n---------- Please answer yes/no [no]  ----------")
if (BASH)
execute_process(COMMAND "${BASH}" "-c"
[=[
read var
if [[ ( "$var" =~ 'y' ) || ( "$var" =~ 'Y' ) ]]; then
    echo "yes"
else
    echo "no"
    exit 1
fi
]=]
                WORKING_DIRECTORY ${CMAKE_CURRENT_LIST_DIR}
                OUTPUT_VARIABLE yesno_answer
                #OUTPUT_FILE /dev/stdout
                RESULT_VARIABLE ret_code
                OUTPUT_STRIP_TRAILING_WHITESPACE
                ERROR_STRIP_TRAILING_WHITESPACE
                )
else()
message("User IO not supported, assuming 'no'")
set(yesno_answer "no")
endif()
else()
set(yesno_answer "no")
endif()

endmacro()

