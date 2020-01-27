# add a target to generate API documentation with Doxygen
if(FLEUR_USE_DOXYGEN)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/docs/Doxyfile.in ${CMAKE_CURRENT_BINARY_DIR}/docs/Doxyfile @ONLY)
add_custom_target(doc
${DOXYGEN_EXECUTABLE} ${CMAKE_CURRENT_BINARY_DIR}/docs/Doxyfile
WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/docs
COMMENT "Generating API documentation with Doxygen" VERBATIM
)
endif()
