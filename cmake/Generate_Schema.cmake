find_program(XXD_PROG xxd)
if (XXD_PROG)
  ADD_CUSTOM_COMMAND(
        OUTPUT ${CMAKE_BINARY_DIR}/include/inputSchema.h
        COMMAND ${XXD_PROG} -i FleurInputSchema.xsd ${CMAKE_BINARY_DIR}/include/inputSchema.h
        DEPENDS ${CMAKE_SOURCE_DIR}/io/xml/FleurInputSchema.xsd
	WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}/io/xml/
        COMMENT "Putting current Schema into inputSchema.h")
else()
  ADD_CUSTOM_COMMAND(
        OUTPUT ${CMAKE_BINARY_DIR}/include/inputSchema.h
        COMMAND cp ${CMAKE_SOURCE_DIR}/io/xml/inputSchema.h.backup ${CMAKE_BINARY_DIR}/include/inputSchema.h
        COMMENT "No xxd found using backup")
  message("No xxd command found! Using backup of inputSchema.h")
endif()     
