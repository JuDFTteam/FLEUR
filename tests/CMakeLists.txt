enable_testing()

set(Testdirs CuBulk CuBulkXML Fe_1l Fe-Atom CuBand CuBandXML CuDOS CuDOSXML Fe_bct_LO Fe_bct NiO_ldau PTO Fe_fcc)
set(Testdirs_INVS CuBulk CuBulkXML Fe_1l Fe-Atom CuBand CuBandXML CuDOS CuDOSXML)
set(Testdirs_SOC Bi2Te3 Fe_1l_SOC PTO-SOC Fe_bct_SOC)

#The serial tests
if (${Fleur_uses_serial})
   #tests of fleur
   foreach(test ${Testdirs})
    add_test(${test} ${CMAKE_CURRENT_SOURCE_DIR}/tests/test.pl ${test} "${CMAKE_BINARY_DIR}/fleur")
   endforeach(test)
   #test of fleur_INVS
   foreach(test ${Testdirs_INVS})
    add_test("I:${test}" ${CMAKE_CURRENT_SOURCE_DIR}/tests/test.pl ${test} "${CMAKE_BINARY_DIR}/fleur_INVS")
   endforeach(test)
   #test of fleur_SOC
   foreach(test ${Testdirs_SOC})
    add_test("S:${test}" ${CMAKE_CURRENT_SOURCE_DIR}/tests/test.pl ${test} "${CMAKE_BINARY_DIR}/fleur_SOC")
   endforeach(test)
endif()

#The parallel tests
if (${Fleur_uses_MPI})
   #tests of fleur
   foreach(test ${Testdirs})
    add_test("MPI-${test}" ${CMAKE_CURRENT_SOURCE_DIR}/tests/test.pl ${test} "${CMAKE_BINARY_DIR}/fleur_MPI")
   endforeach(test)
   #test of fleur_INVS
   foreach(test ${Testdirs_INVS})
    add_test("MPI-I:${test}" ${CMAKE_CURRENT_SOURCE_DIR}/tests/test.pl ${test} "${CMAKE_BINARY_DIR}/fleur_INVS_MPI")
   endforeach(test)
   #test of fleur_SOC
   foreach(test ${Testdirs_SOC})
    add_test("MPI-S:${test}" ${CMAKE_CURRENT_SOURCE_DIR}/tests/test.pl ${test} "${CMAKE_BINARY_DIR}/fleur_SOC_MPI")
   endforeach(test)
endif()



