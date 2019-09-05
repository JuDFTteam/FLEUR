set(fleur_F77 ${fleur_F77}
math/angle.f
math/clebsch.f
math/cylbes.f
math/cylpts.f
math/d_wigner.F
math/dcylbs.f
math/dsphbs.f
math/gaussp.f
math/grule.f
math/inv3.f
math/inwint.f
math/matmul.f
math/outint.f
math/points.f
math/qsf.f
math/sphbes.f
math/sphpts.f
math/util.F
math/difcub.f
)
set(fleur_F90 ${fleur_F90}
math/cfft.F90
math/rfft.f90
math/qranf.f90
math/differentiate.f90
math/fft2d.F90
math/fft3d.f90
math/fft_interface.F90
math/SphBessel.f90
math/DoubleFactorial.f90
math/ExpSave.f90
math/intgr.F90
math/ylm4.F90
)
if (FLEUR_USE_FFTMKL)
   set(fleur_F90 ${fleur_F90} math/mkl_dfti.f90)
endif()
