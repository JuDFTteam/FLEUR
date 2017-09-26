enable_testing()

set(SerialParallelTests CuBulk CuBulkXML Fe_1l Fe_1lXML Fe-Atom CuBand CuBandXML CuDOS CuDOSXML
Fe_bct Fe_bctXML PTO PTOXML Fe_1l_SOC Fe_1l_SOCXML PTO-SOC PTO-SOCXML Fe_bct_SOC Fe_bct_SOCXML
GaAsMultiUForceXML SiFilmSlicePlotXML)

set(SerialOnlyTests Fe_bct_LO Fe_bct_LOXML Fe_fcc Fe_fccXML TiO2eels TiO2eelsXML)

set(Testdirs ${SerialParallelTests} ${SerialOnlyTests})
set(ParTestdirs ${SerialParallelTests})

#Check if all tests (including those running for a long time) should be executed
if (all_tests)
    set(Testdirs ${Testdirs} Bi2Te3 Bi2Te3XML NiO_ldau NiO_ldauXML)
    set(ParTestdirs ${Testdirs} Bi2Te3 Bi2Te3XML NiO_ldau NiO_ldauXML)
endif()

#Add Wannier tests if fleur is compiled with Wannier support
if (${FLEUR_USE_WANN})
    set(Testdirs ${Testdirs} Cwann)
    set(ParTestdirs ${ParTestdirs} Cwann)
endif()

#The serial tests
if (${FLEUR_USE_SERIAL})
   foreach(test ${Testdirs})
    add_test("FLEUR:${test}" ${CMAKE_CURRENT_SOURCE_DIR}/tests/test.pl ${test} "${CMAKE_BINARY_DIR}/fleur")
   endforeach(test)
endif()

#The parallel tests
if (${FLEUR_USE_MPI})
   if (MPIEXEC)
      set(mpi_exec "${MPIEXEC} ${MPI_NUMPROC_FLAGS} 2")
   else()
      set(mpi_exec "mpirun -n 2")
   endif()
   foreach(test ${ParTestdirs})
    add_test("FLEUR_MPI:${test}" ${CMAKE_CURRENT_SOURCE_DIR}/tests/test.pl
${test} "${CMAKE_BINARY_DIR}/fleur_MPI" "${mpi_exec}")
   endforeach(test)
endif()
