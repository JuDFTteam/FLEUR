set(fleur_F77 ${fleur_F77}
)
set(fleur_F90 ${fleur_F90}
eigen/apws.F90
eigen/eigen.F90
eigen/hlomat.F90
eigen/hs_int.F90
eigen/hsmt_fjgj.F90
eigen/hsmt_ab.f90
eigen/hsmt_soc_init.F90
eigen/hsmt_sph.F90
eigen/hsmt_nonsph.F90
eigen/hsmt_spinor.F90
eigen/hsmt_distspin.F90
eigen/hsmt_lo.F90
eigen/mt_setup.F90
eigen/eigen_hssetup.F90
eigen/hsmt.F90
eigen/hsvac.F90
eigen/lodpot.f90
eigen/od_hsvac.F90
eigen/od_vacfun.f90
eigen/orthoglo.F90
eigen/rad_ovlp.f90
eigen/setabc1lo.f90
eigen/setabc1locdn.f90
eigen/setabc1locdn1.f90
eigen/slomat.F90
eigen/tlmplm.F90
eigen/tlmplm_cholesky.F90
eigen/tlmplm_store.F90
eigen/tlo.f90
eigen/vacfun.f90
eigen/vec_for_lo.f90
)
if (FLEUR_USE_GPU)
   set(fleur_F90 ${fleur_F90} eigen/hsmt_nonsph_GPU.F90)
endif()
