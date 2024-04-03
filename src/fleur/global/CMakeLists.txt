set(fleur_F77 ${fleur_F77}
src/fleur/global/fitchk.f
src/fleur/global/radsra.f
src/fleur/global/radsrd.F
src/fleur/global/radsrdn.f
src/fleur/global/soc_sym.f
src/fleur/global/vacudz.f
src/fleur/global/vacuz.f
src/fleur/global/abcrot.f
src/fleur/global/differ.f
)
set(fleur_F90 ${fleur_F90}
src/fleur/global/savewigner.f90
#src/fleur/global/differ.f90
#src/fleur/global/matrix_copy.F90
src/fleur/global/checkdop.F90
src/fleur/global/checkdopall.f90
src/fleur/global/sort.f90
src/fleur/global/genMTBasis.f90
src/fleur/global/chkmt.f90
src/fleur/global/convn.f90
src/fleur/global/phasy1.f90
src/fleur/global/qfix.f90
src/fleur/global/radflo.F90
src/fleur/global/starf.f90
src/fleur/global/utility.F90
src/fleur/global/find_enpara.f90
src/fleur/global/rotMMPmat.f90
src/fleur/global/symMMPmat.f90
)
if (FLEUR_USE_GPU)
  set(fleur_F90 ${fleur_F90}

  )
endif()
