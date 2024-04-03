module FFTW3
#ifdef CPP_FFTW
   use, intrinsic :: iso_c_binding
   include 'fftw3.f03'
#endif
end module FFTW3