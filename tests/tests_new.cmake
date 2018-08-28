#Loop over testcases subdirectories and add test.cmake 

file(GLOB dirs "tests/testcases/*")
foreach(dir ${dirs})
    include(${dir}/test.cmake)
endforeach()