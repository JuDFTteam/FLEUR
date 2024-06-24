#this file sets up the FLEURist script

set(Python3_FIND_VIRTUALENV FIRST)
#This option prevents cmake finding a newer global python version
#if a virtualenv with a older version is explicitely activated
set(Python3_FIND_STRATEGY LOCATION)
find_package(Python3)
message("Python3 found:${Python3_FOUND}")
message("Python3 path:${Python3_EXECUTABLE}")
message("The python executable used for FLEURist"
        "can be overwritten with the juDFT_PYTHON environment variable")

if( Python3_FOUND )
    set(FLEUR_PYTHON ${Python3_EXECUTABLE})
else()
    set(FLEUR_PYTHON "python3")
endif()

#write build script
file(GENERATE OUTPUT ${CMAKE_BINARY_DIR}/FLEURist CONTENT
"#!/usr/bin/env bash
PYTHON_EXECUTABLE=\"${FLEUR_PYTHON}\"
if [[ ! -z \"\${juDFT_PYTHON}\" ]]; then
  PYTHON_EXECUTABLE=\${juDFT_PYTHON}
fi
$PYTHON_EXECUTABLE \"${CMAKE_SOURCE_DIR}/src/tools/fleurist/fleurist.py\" \"$@\"  
")

add_custom_target(fleurist ALL
                  COMMAND chmod +x FLEURist
                  WORKING_DIRECTORY ${CMAKE_BINARY_DIR}
                  COMMENT "Making FLEURist script executable")


