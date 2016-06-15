enable_testing()

add_test(CuBulk ${CMAKE_CURRENT_SOURCE_DIR}/tests/test.pl CuBulk "${CMAKE_BINARY_DIR}/fleur")
#add_test(NAME test2 WORKING_DIRECTORY tests COMMAND $ENV{HOME}/bin/testcmake a1 a2 )
add_test(test $ENV{HOME}/bin/testcmake arg1 arg2 )
