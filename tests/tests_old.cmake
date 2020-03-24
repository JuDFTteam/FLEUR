enable_testing()

set(SerialParallelTests  CuBulkXML SiLOXML  Fe_1lXML
   CuBandXML  CuDOSXML  Fe_bctXML  PTOXML Fe_1l_SOCXML
   PTO-SOCXML  Fe_bct_SOCXML Fe_fccXML GaAsMultiUForceXML
   CoMCDXML Fe_Kerker Fe_bct_LOXML
   FePt_film_SSFT FePt_film_SSFT_LO
   GaAsHybridPBE0 FeHybridPBE0 MnHybridNoinv
   Fe_bcc_GreensFunction GreensFunction_MultiContour 
   Fe_bcc_FlipcdnXLDA Fe_bcc_FlipcdnYGGA FeFFNLOsSOC
    Fe_bcc_SF_LDA SmAtomjDOS) 

#Currently disabled Tests (Hybrid+Greenfct)
# SiHybridGammaNoInv SiHybrid8kpt_sym  SiHybrid8kpt_nosym
# Fe_1l_GreensFunction
#  CoUnfold

#Tests to renable soon
#RelaxMTFeature
#SiFilmPlotXML PlotDenandPot SiFilmSlicePlotXML

set(SerialOnlyTests  )
set(InpgenTests Si_plain Si_plain_explicit Si_full_para)# Si_kpt Si_kden Si_round_trip)


#Removed GW tests
#if (${FLEUR_USE_HDF5})
#   set(SerialOnlyTests ${SerialOnlyTests} gw1Interface gw2Interface)
#endif()

set(Testdirs ${SerialParallelTests} ${SerialOnlyTests})

set(ParTestdirs ${SerialParallelTests})
set(InpTestdirs ${InpgenTests})

#Check if all tests (including those running for a long time) should be executed
if (all_tests)
    set(Testdirs ${Testdirs} Bi2Te3 Bi2Te3XML NiO_ldauXML)
    set(ParTestdirs ${ParTestdirs} Bi2Te3 Bi2Te3XML NiO_ldauXML)
endif()

#The inpgen tests
#if (${INPGEN})
foreach(test ${InpTestdirs})
 add_test("INPGEN:${test}" ${CMAKE_CURRENT_SOURCE_DIR}/tests/test.pl "inpgen/${test}" "${CMAKE_BINARY_DIR}/inpgen2/inpgen2")
endforeach(test)
#endif()

#Add Wannier tests if fleur is compiled with Wannier support
if (${FLEUR_USE_WANN})
    set(Testdirs ${Testdirs} Cwann CwannXML)
    set(ParTestdirs ${ParTestdirs} Cwann CwannXML)
endif()

#Tests for LibXC
if (${FLEUR_USE_LIBXC})
   set(Testdirs ${Testdirs} CuBulkLibXC Fe_bct_LibXC Diamond_SCAN)
   set(ParTestdirs ${ParTestdirs} CuBulkLibXC Fe_bct_LibXC Diamond_SCAN)
endif()

#Tests for EDsolver
if (${FLEUR_USE_EDSOLVER})
   set(Testdirs ${Testdirs} Gd_Hubbard1 Gd_Hubbard1_SOC Gd_Hubbard1_noSYM)
   set(ParTestdirs ${ParTestdirs} Gd_Hubbard1 Gd_Hubbard1_SOC Gd_Hubbard1_noSYM)
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
    add_test("FLEUR_MPI:${test}" ${CMAKE_CURRENT_SOURCE_DIR}/tests/test.pl ${test} "${CMAKE_BINARY_DIR}/fleur_MPI" "${mpi_exec}")
   endforeach(test)
   set(mpi_exec "sequential")
   foreach(test ${SerialOnlyTests})
    add_test("FLEUR_MPI:${test}" ${CMAKE_CURRENT_SOURCE_DIR}/tests/test.pl ${test} "${CMAKE_BINARY_DIR}/fleur_MPI" "${mpi_exec}")
   endforeach(test)
endif()

#Hybrid tests
add_test("KClHybridPBE0" ${CMAKE_CURRENT_SOURCE_DIR}/tests/tests/KClHybridPBE0/test.py --bindir ${CMAKE_BINARY_DIR} --testdir ${CMAKE_BINARY_DIR}/Testing)
