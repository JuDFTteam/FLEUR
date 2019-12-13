set(fleur_F77 ${fleur_F77}

)
set(fleur_F90 ${fleur_F90}
fft/cfft.F90
fft/rfft.f90
fft/fft2d.F90
fft/fft3d.f90
fft/fft_interface.F90
fft/selectFFT.F90
)
if (FLEUR_USE_FFTMKL)
   set(fleur_F90 ${fleur_F90} fft/mkl_dfti.f90)
endif()
