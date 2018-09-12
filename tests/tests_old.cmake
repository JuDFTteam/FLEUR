enable_testing()

set(SerialParallelTests CuBulk CuBulkXML SiLOXML Fe_1l Fe_1lXML Fe-Atom CuBand CuBandXML CuDOS CuDOSXML
Fe_bct Fe_bctXML PTO PTOXML Fe_1l_SOC Fe_1l_SOCXML PTO-SOC PTO-SOCXML Fe_bct_SOC Fe_bct_SOCXML Fe_fccXML
GaAsMultiUForceXML SiFilmPlotXML SiFilmSlicePlotXML CoMCDXML Fe_Kerker Fe_bct_LOXML)

set(SerialOnlyTests Fe_bct_LO Fe_fcc)# TiO2eels TiO2eelsXML)
set(InpgenTests Si_plain Si_plain_explicit Si_full_para)# Si_kpt Si_kden Si_round_trip) 

set(Testdirs ${SerialParallelTests} ${SerialOnlyTests})
set(ParTestdirs ${SerialParallelTests})
set(InpTestdirs ${InpgenTests})

#Check if all tests (including those running for a long time) should be executed
if (all_tests)
    set(Testdirs ${Testdirs} Bi2Te3 Bi2Te3XML NiO_ldauXML)
    set(ParTestdirs ${Testdirs} Bi2Te3 Bi2Te3XML NiO_ldauXML)
endif()

#The inpgen tests
#if (${INPGEN})
foreach(test ${InpTestdirs})
 add_test("INPGEN:${test}" ${CMAKE_CURRENT_SOURCE_DIR}/tests/test.pl "inpgen/${test}" "${CMAKE_BINARY_DIR}/inpgen")
endforeach(test)
#endif()

#Add Wannier tests if fleur is compiled with Wannier support
if (${FLEUR_USE_WANN})
    set(Testdirs ${Testdirs} Cwann CwannXML)
    set(ParTestdirs ${ParTestdirs} Cwann CwannXML)
endif()

#Tests for LibXC
if (${FLEUR_USE_LIBXC})
    set(Testdirs ${Testdirs} CuBulkLibXC Fe_bct_LibXC)
    set(ParTestdirs ${ParTestdirs} CuBulkLibXC Fe_bct_LibXC)
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

