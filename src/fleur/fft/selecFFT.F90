MODULE m_selecFFT

   USE m_juDFT

   IMPLICIT NONE

   INTEGER, PARAMETER :: defaultFFT_const = 0
   INTEGER, PARAMETER :: mklFFT_const     = 1
   INTEGER, PARAMETER :: spFFT_const      = 2
   INTEGER, PARAMETER :: fftw_const       = 3
   integer, parameter :: cufft_const      = 4

#ifdef CPP_FFT_MKL
   LOGICAL, PARAMETER :: mklFFT_available = .TRUE.
#else
   LOGICAL, PARAMETER :: mklFFT_available = .FALSE.
#endif

#ifdef CPP_FFTW
   LOGICAL, PARAMETER :: fftw_available = .TRUE.
#else
   LOGICAL, PARAMETER :: fftw_available = .FALSE.
#endif

#ifdef CPP_SPFFT
   LOGICAL, PARAMETER :: spFFT_available = .TRUE.
#else
   LOGICAL, PARAMETER :: spFFT_available = .FALSE.
#endif

#ifdef _OPENACC
   LOGICAL, PARAMETER :: cuFFT_available = .TRUE.
#else
   LOGICAL, PARAMETER :: cuFFT_available = .FALSE.
#endif


CONTAINS

   FUNCTION selecFFT(l_sparse, l_gpu)

      INTEGER             :: selecFFT
      LOGICAL, INTENT(IN) :: l_sparse
      logical,optional, intent(in) :: l_gpu

      INTEGER :: fftRoutine

      fftRoutine = defaultFFT_const
#ifdef CPP_FFTW
      fftRoutine = fftw_const
#endif

#ifdef CPP_FFT_MKL
      fftRoutine = mklFFT_const
#endif

#ifdef CPP_SPFFT
      IF(l_sparse) fftRoutine = spFFT_const
#endif

#ifdef _OPENACC
      if(present(l_gpu)) then
         if(l_gpu) fftRoutine = cufft_const
      endif
#endif

      IF (TRIM(juDFT_string_for_argument("-fft"))=="fftw") THEN
#ifdef CPP_FFTW
         fftRoutine = fftw_const
#else
         CALL juDFT_error("Selected fftw is not available!", calledby="selecFFT")
#endif

      END IF
      IF (TRIM(juDFT_string_for_argument("-fft"))=="mkl") THEN
#ifdef CPP_FFT_MKL
         fftRoutine = mklFFT_const
#else
         CALL juDFT_error("Selected FFT (mkl) is not available!", calledby="selecFFT")
#endif
      END IF

      IF (TRIM(juDFT_string_for_argument("-fft"))=="spfft") THEN
#ifdef CPP_SPFFT
         IF(l_sparse) fftRoutine = spFFT_const
#else
         CALL juDFT_error("Selected FFT (spfft) is not available!", calledby="selecFFT")
#endif
      END IF

      IF (TRIM(juDFT_string_for_argument("-fft"))=="cufft") THEN
#ifdef _OPENACC
         if(present(l_gpu)) then
            IF(l_gpu) fftRoutine = cuFFT_const
         endif
#else
         CALL juDFT_error("Selected FFT (spfft) is not available!", calledby="selecFFT")
#endif
      END IF

      IF (TRIM(juDFT_string_for_argument("-fft"))=="inbuilt") THEN
         IF(l_sparse) fftRoutine = defaultFFT_const
      END IF

      selecFFT = fftRoutine

   END FUNCTION

END MODULE m_selecFFT
