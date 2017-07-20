set(fleur_F90 ${fleur_F90}
hybrid/abcrot.F90          hybrid/exchange_core.F90   hybrid/mixedbasis.F90      hybrid/symmetrizeh.F90
hybrid/checkolap.F90       hybrid/exchange_val_hf.F90 hybrid/read_core.F90       hybrid/wavefproducts.F90
hybrid/coulombmatrix.F90   hybrid/gen_wavf.F90        hybrid/spmvec.F90          hybrid/hybrid.F90
hybrid/HF_init.F90   hybrid/hsfock.F90          hybrid/subvxc.F90          hybrid/add_Vnonlocal.F90
hybrid/hf_setup.F90  hybrid/kp_perturbation.F90 hybrid/symm_hf.F90 hybrid/trafo.F90
hybrid/exponential_integral.f90
hybrid/hsefunctional.F90
hybrid/olap.F90
)

set(fleur_F77 ${fleur_F77}
hybrid/wrapper.F
)
