set(fleur_F77 ${fleur_F77}
)
set(fleur_F90 ${fleur_F90}
src/fleur/eigen/hsmt_mtNocoPot_offdiag.F90
src/fleur/eigen/eigen.F90
src/fleur/eigen/hlomat.F90
src/fleur/eigen/hs_int.F90
src/fleur/eigen/hsmt_fjgj.F90
src/fleur/eigen/hsmt_ab.F90
src/fleur/eigen/hsmt_sph.F90
src/fleur/eigen/hsmt_nonsph.F90
src/fleur/eigen/hsmt_spinor.F90
src/fleur/eigen/hsmt_distspin.F90
src/fleur/eigen/hsmt_offdiag.F90
src/fleur/eigen/hsmt_soc_offdiag.F90
src/fleur/eigen/hsmt_lo.F90
src/fleur/eigen/local_hamiltonian.F90
src/fleur/eigen/eigen_hssetup.F90
src/fleur/eigen/hsmt.F90
src/fleur/eigen/hsvac.F90
src/fleur/eigen/orthoglo.F90
src/fleur/eigen/rad_ovlp.f90
src/fleur/eigen/setabc1lo.f90
src/fleur/eigen/slomat.F90
src/fleur/eigen/tlmplm.F90
#src/fleur/eigen/tlmplm_store.F90
src/fleur/eigen/tlo.f90
src/fleur/eigen/vacfun.F90
src/fleur/eigen/eigen_redist_matrix.F90
src/fleur/eigen/hs_int_direct.F90
)
#if (FLEUR_USE_GPU)
#   set(fleur_F90 ${fleur_F90}
#src/fleur/eigen/hsmt_nonsph_GPU.F90)
#endif()
