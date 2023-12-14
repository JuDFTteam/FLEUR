set(fleur_F77 ${fleur_F77}
)
set(fleur_F90 ${fleur_F90}
eigen/hsmt_mtNocoPot_offdiag.F90
eigen/eigen.F90
eigen/hslomat.F90
eigen/hs_int.F90
eigen/hsmt_fjgj.F90
eigen/hsmt_ab.F90
eigen/hsmt_sph.F90
eigen/hsmt_nonsph.F90
eigen/hsmt_spinor.F90
eigen/hsmt_distspin.F90
eigen/hsmt_offdiag.F90
eigen/hsmt_soc_offdiag.F90
eigen/hsmt_lo.F90
eigen/local_hamiltonian.F90
eigen/eigen_hssetup.F90
eigen/hsmt.F90
eigen/hsvac.F90
eigen/orthoglo.F90
eigen/rad_ovlp.f90
eigen/setabc1lo.f90
#eigen/slomat.F90
eigen/tlmplm.F90
#eigen/tlmplm_store.F90
eigen/tlo.f90
eigen/vacfun.F90
eigen/eigen_redist_matrix.F90
eigen/hs_int_direct.F90
)
#if (FLEUR_USE_GPU)
#   set(fleur_F90 ${fleur_F90}
#eigen/hsmt_nonsph_GPU.F90)
#endif()
