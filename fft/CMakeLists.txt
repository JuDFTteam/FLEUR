set(fleur_F77 ${fleur_F77}

)
set(fleur_F90 ${fleur_F90}
fft/cfft.F90
fft/rfft.f90
fft/fft2d.F90
fft/fft3d.f90
fft/fft_interface.F90
fft/selecFFT.F90
fft/fftw3.F90
)
if (MKL_FFT_INTERFACE_FILE)
   set(fleur_F90 ${fleur_F90} ${MKL_FFT_INTERFACE_FILE})
endif()
