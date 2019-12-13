set(fleur_F77 ${fleur_F77}
global/fitchk.f
global/nmat_rot.f
global/od_cylbes.f
global/od_phasy.f
global/radsra.f
global/radsrd.F
global/radsrdn.f
global/soc_sym.f
global/starf.f
global/triang.f
global/vacudz.f
global/vacuz.f
global/abcrot.f
global/differ.f
)
set(fleur_F90 ${fleur_F90}
global/savewigner.f90
#global/differ.f90
#global/matrix_copy.F90
global/checkdop.F90
global/checkdopall.f90
global/sort.f90
global/genMTBasis.f90
global/chkmt.f90
global/convn.f90
global/phasy1.f90
global/qfix.f90
global/radflo.F90
global/utility.F90
global/find_enpara.f90
)
if (FLEUR_USE_GPU)
  set(fleur_F90 ${fleur_F90}
  global/nvtx.F90
  )
endif()
