set(fleur_F77 ${fleur_F77}
math/angle.f
math/cfft.F
math/clebsch.f
math/cylbes.f
math/cylpts.f
math/d_wigner.F
math/dcylbs.f
math/dsphbs.f
math/gaussp.f
math/grule.f
math/intgr.F
math/inv3.f
math/inwint.f
math/matmul.f
math/outint.f
math/points.f
math/qranf.f
math/qsf.f
math/rfft.F
math/sphbes.f
math/sphpts.f
math/util.F
math/ylm4.f
)
set(fleur_F90 ${fleur_F90}
math/differentiate.f90
math/fft2d.F90
math/fft3d.f90
math/fft_interface.F90
)
if (FLEUR_USE_FFTMKL)
   set(fleur_F90 ${fleur_F90} math/mkl_dfti.f90)
endif()
