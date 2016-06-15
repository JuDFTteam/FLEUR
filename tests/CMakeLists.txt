enable_testing()

# all the following lines should have the following structure:
#add_test(TESTNAME ${CMAKE_CURRENT_SOURCE_DIR}/tests/test.pl TESTDIR "${CMAKE_BINARY_DIR}/BINDIR MPIEXEC)
#where:
#    TESTNAME is the name of the test (your choice)
#    TESTDIR is the subdirectory of tests/tests in the source directory in which the test is defined
#    BINDIR is the name of the executable to test (i.e. fleur or fleur_INVS or fleur_SOC or ...)
#    MPIEXEC (optional) is the additional stuff you have to add to start an MPI programm


add_test(CuBulk ${CMAKE_CURRENT_SOURCE_DIR}/tests/test.pl CuBulk "${CMAKE_BINARY_DIR}/fleur")
