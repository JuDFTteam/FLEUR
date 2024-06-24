file(REMOVE ${CMAKE_SOURCE_DIR}/src/fleur/io/xml/inputSchema.h)
find_program(XXD_PROG xxd)
if (XXD_PROG)
  ADD_CUSTOM_COMMAND(
        OUTPUT ${CMAKE_BINARY_DIR}/include/inputSchema.h
        COMMAND ${XXD_PROG} -i FleurInputSchema.xsd ${CMAKE_BINARY_DIR}/include/inputSchema.h
        DEPENDS ${CMAKE_SOURCE_DIR}/src/fleur/io/xml/FleurInputSchema.xsd
	WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}/src/fleur/io/xml/
        COMMENT "Putting current inputSchema into inputSchema.h")
  ADD_CUSTOM_COMMAND(
        OUTPUT ${CMAKE_BINARY_DIR}/include/outputSchema.h
        COMMAND ${XXD_PROG} -i FleurOutputSchema.xsd ${CMAKE_BINARY_DIR}/include/outputSchema.h
        DEPENDS ${CMAKE_SOURCE_DIR}/src/fleur/io/xml/FleurOutputSchema.xsd
  WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}/src/fleur/io/xml/
        COMMENT "Putting current outputSchema into outputSchema.h")
else()
  ADD_CUSTOM_COMMAND(
        OUTPUT ${CMAKE_BINARY_DIR}/include/inputSchema.h
        COMMAND cp ${CMAKE_SOURCE_DIR}/src/fleur/io/xml/inputSchema.h.backup ${CMAKE_BINARY_DIR}/include/inputSchema.h
        COMMENT "No xxd found using backup")
  message("No xxd command found! Using backup of inputSchema.h")
  ADD_CUSTOM_COMMAND(
        OUTPUT ${CMAKE_BINARY_DIR}/include/outputSchema.h
        COMMAND cp ${CMAKE_SOURCE_DIR}/src/fleur/io/xml/outputSchema.h.backup ${CMAKE_BINARY_DIR}/include/outputSchema.h
        COMMENT "No xxd found using backup")
  message("No xxd command found! Using backup of outputSchema.h")
endif()     
