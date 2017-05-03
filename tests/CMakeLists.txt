enable_testing()

set(Testdirs CuBulk CuBulkXML Fe_1l Fe_1lXML Fe-Atom CuBand CuBandXML CuDOS CuDOSXML Fe_bct_LO Fe_bct_LOXML Fe_bct Fe_bctXML NiO_ldau NiO_ldauXML PTO PTOXML Fe_fcc Fe_fccXML Bi2Te3 Bi2Te3XML Fe_1l_SOC Fe_1l_SOCXML PTO-SOC PTO-SOCXML Fe_bct_SOC Fe_bct_SOCXML) 
#The serial tests
if (${FLEUR_USE_SERIAL})
   #tests of fleur
   foreach(test ${Testdirs})
    add_test("FLEUR:${test}" ${CMAKE_CURRENT_SOURCE_DIR}/tests/test.pl ${test} "${CMAKE_BINARY_DIR}/fleur")
   endforeach(test)
endif()

#The parallel tests
if (${FLEUR_USE_MPI})
   #tests of fleur
   foreach(test ${Testdirs})
    add_test("FLEUR_MPI:${test}" ${CMAKE_CURRENT_SOURCE_DIR}/tests/test.pl ${test} "${CMAKE_BINARY_DIR}/fleur_MPI" "${MPIEXEC} ${MPI_NUMPROC_FLAGS 2")
   endforeach(test)
endif()
