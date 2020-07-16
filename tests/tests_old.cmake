enable_testing()

set(SerialParallelTests CuBulkXML SiLOXML  Fe_1lXML
    Fe_bctXML  PTOXML Fe_1l_SOCXML
   PTO-SOCXML  Fe_bct_SOCXML Fe_fccXML GaAsMultiUForceXML H2ORelaxBFGS
     Fe_Kerker Fe_bct_LOXML SiFilmPlotXML SiFilmSlicePlotXML
   FePt_film_SSFT FePt_film_SSFT_LO
   Fe_bcc_GreensFunction GreensFunction_MultiContour Fe_1l_GreensFunction
   Fe_1l_Tria GreensFunctionRadial GreensFunctionRadial_LO )

#DOS related tests, need to be enabled again!
set(SerialParallelTests ${SerialParallelTests}
  CuBandXML  CuDOSXML CoMCDXML SmAtomjDOS )

#Currently disabled Tests (Hybrid)
# SiHybridGammaNoInv SiHybrid8kpt_sym  SiHybrid8kpt_nosym
#  CoUnfold
#Other disabled tests
#FePt_film_SSFT FePt_film_SSFT_LO

set(SerialOnlyTests  )

set(InpgenTests Si_plain Si_plain_explicit Si_full_para)# Si_kpt Si_kden Si_round_trip)

set(HybridTests
   KClHybridPBE0
   GaAsHybridPBE0
   FeHybridPBE0
   MnHybridNoinv
)

set(FFNTests
   Fe_bcc_FlipcdnXLDA Fe_bcc_FlipcdnYGGA FeFFNLOsSOC
   PlotDenandPot PlotOnlyMT Noncollinear_downward_compatible
   RelaxMTFeature Fe_bcc_SF_LDA
)


if (FLEUR_USE_HDF5)
    set(SerialParallelTests ${SerialParallelTests} ${FFNTests})
endif()

#Check if all tests (including those running for a long time) should be executed
if (all_tests)
   set(SerialParallelTests ${SerialParallelTests} Bi2Te3 Bi2Te3XML NiO_ldauXML ${FFNTests})
endif()

#Add Wannier tests if fleur is compiled with Wannier support
if (FLEUR_USE_WANN)
   set(SerialOnlyTests ${SerialOnlyTests} CwannXML)
endif()

#Tests for LibXC
if (FLEUR_USE_LIBXC)
   set(SerialParallelTests ${SerialParallelTests} CuBulkLibXC Fe_bct_LibXC Diamond_SCAN)
endif()

#Tests for EDsolver
if (FLEUR_USE_EDSOLVER)
   set(SerialParallelTests ${SerialParallelTests} Gd_Hubbard1 Gd_Hubbard1_noSYM)
endif()


#The inpgen tests
#if (INPGEN)
foreach(test ${InpgenTests})
   add_test("INPGEN:${test}" ${CMAKE_CURRENT_SOURCE_DIR}/tests/test.pl "inpgen/${test}" "${CMAKE_BINARY_DIR}/inpgen2/inpgen2")
endforeach(test)
#endif()

#The serial tests
if (FLEUR_USE_SERIAL)
   foreach(test ${SerialParallelTests} ${SerialOnlyTests})
    add_test("FLEUR:${test}" ${CMAKE_CURRENT_SOURCE_DIR}/tests/test.pl ${test} "${CMAKE_BINARY_DIR}/fleur")
   endforeach(test)
endif()

#The parallel tests
if (FLEUR_USE_MPI)
   if (MPIEXEC)
      set(mpi_exec "${MPIEXEC} ${MPI_NUMPROC_FLAGS} 2")
   else()
      set(mpi_exec "mpirun -n 2")
   endif()
   foreach(test ${SerialParallelTests})
    add_test("FLEUR_MPI:${test}" ${CMAKE_CURRENT_SOURCE_DIR}/tests/test.pl ${test} "${CMAKE_BINARY_DIR}/fleur_MPI" "${mpi_exec}")
   endforeach(test)
   set(mpi_exec "sequential")
   foreach(test ${SerialOnlyTests})
    add_test("FLEUR_MPI:${test}" ${CMAKE_CURRENT_SOURCE_DIR}/tests/test.pl ${test} "${CMAKE_BINARY_DIR}/fleur_MPI" "${mpi_exec}")
   endforeach(test)
endif()

#Hybrid tests
foreach(test ${HybridTests})
   add_test("${test}" ${CMAKE_CURRENT_SOURCE_DIR}/tests/tests/${test}/test.py --bindir ${CMAKE_BINARY_DIR} --testdir ${CMAKE_BINARY_DIR}/Testing/${test})
endforeach()
