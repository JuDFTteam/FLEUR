set(fleur_F77 ${fleur_F77}

)
set(fleur_F90 ${fleur_F90}
src/fleur/fft/cfft.F90
src/fleur/fft/rfft.f90
src/fleur/fft/fft2d.F90
src/fleur/fft/fft3d.f90
src/fleur/fft/fft_interface.F90
src/fleur/fft/selecFFT.F90
src/fleur/fft/fftw3.F90
)
if (MKL_FFT_INTERFACE_FILE)
   set(fleur_F90 ${fleur_F90} ${MKL_FFT_INTERFACE_FILE})
endif()
