MODULE m_selectFFT

   USE m_juDFT
   USE m_types_fft
   use m_types_fft_mkl
   use m_types_fft_singleton
   use m_types_fft_spfft
   use m_types_fft_fftw
   !use m_types_fft_cufft

   IMPLICIT NONE
#if CPP_FFT_TYPE==t_fft_singleton
   INTEGER,PARAMETER:: defaultlist(1)=[fft_singleton]
#elif CPP_FFT_TYPE==t_fft_mkl
   INTEGER,PARAMETER:: defaultlist(1)=[fft_mkl]
#elif CPP_FFT_TYPE==t_fft_fftw
   INTEGER,PARAMETER:: defaultlist(1)=[fft_fftw]
#elif CPP_FFT_TYPE==t_fft_spfft
   INTEGER,PARAMETER:: defaultlist(1)=[fft_spfft]
#else
   INTEGER,PARAMETER:: defaultlist(3)=[fft_mkl,fft_fftw,fft_singleton]
#endif



   TYPE t_fftbox
     class(t_fft),allocatable :: fft
   end type

   PUBLIC selectFFT
   CONTAINS

   FUNCTION selectFFT(wishlist)
     INTEGER,optional            :: wishlist(:)
     class(t_fft),ALLOCATABLE    :: selectFFT

     INTEGER,ALLOCATABLE :: list(:)
     INTEGER :: n

      if (present(wishlist)) THEN
        list=wishlist
      else
        list=defaultlist
      endif

      DO n=1,size(list)
        select case(list(n))
        case (fft_mkl)
          if (fft_mkl>0) then
            allocate(t_fft_mkl::selectFFT)
            RETURN
          endif
        case (fft_singleton)
          if (fft_singleton>0) then
            allocate(t_fft_singleton::selectFFT)
            RETURN
          endif
        case (fft_fftw)
          if (fft_fftw>0) then
            allocate(t_fft_fftw::selectFFT)
            RETURN
          endif
        case (fft_spfft)
          if (fft_spfft>0) then
            allocate(t_fft_spfft::selectFFT)
            RETURN
          endif
        end select
      enddo
      call  judft_error("No FFT routine found")

   END FUNCTION

END MODULE m_selectFFT
