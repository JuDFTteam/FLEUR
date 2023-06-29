#First check if we can compile with ELSI

try_compile(FLEUR_USE_ELSI ${CMAKE_BINARY_DIR} ${CMAKE_SOURCE_DIR}/cmake/tests/test_ELSI.f90 LINK_LIBRARIES ${FLEUR_LIBRARIES} OUTPUT_VARIABLE compile_output)

if ("$ENV{VERBOSE}")
    message("ELSI compile test: ${FLEUR_USE_ELSI}\nLINK_LIBRARIES ${TEST_LIBRARIES}\n${compile_output}")
endif()


message("ELSI Library found:${FLEUR_USE_ELSI}")

if (FLEUR_USE_ELSI)
   set(FLEUR_DEFINITIONS ${FLEUR_DEFINITIONS} "CPP_ELSI")
endif()
