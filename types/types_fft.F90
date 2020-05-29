module m_types_fft
#ifdef CPP_FFT_MKL
   USE mkl_dfti
#endif
   USE m_selecFFT
   use m_judft

   type t_fft 
      logical :: initialized = .False.
      integer :: backend = -1
      integer :: length(3) = [-1,-1,-1]
      logical :: forw
      ! cfft storage
      real, allocatable :: afft(:), bfft(:)
#ifdef CPP_FFT_MKL
      ! mkl
      type(dfti_descriptor), pointer :: dfti_handle
#endif
   contains 
      procedure :: init => t_fft_init
      procedure :: exec => t_fft_exec
      procedure :: free => t_fft_free
   end type t_fft
contains
   subroutine t_fft_init(fft, length, forw, indices)
      implicit none       
      class(t_fft)                  :: fft
      integer, intent(in)           :: length(3) !length of data in each direction
      logical, intent(in)           :: forw          !.true. for the forward transformation, .false. for the backward one
      INTEGER, OPTIONAL, INTENT(IN) :: indices(:)    !array of indices of relevant/nonzero elements in the FFT mesh

      integer :: size_dat, ok

      fft%initialized = .True.
      fft%backend = defaultFFT_const
      fft%backend = selecFFT(PRESENT(indices))
      fft%length  = length
      fft%forw    = forw

      select case(fft%backend)
      case(spFFT_const)
         call juDFT_error("not yet implemented")
      case(mklFFT_const)
#ifdef CPP_FFT_MKL
         ok = DftiCreateDescriptor(fft%dfti_handle, dfti_double, dfti_complex, 3, length)
         if (ok /= 0) call juDFT_error("cant create descriptor", calledby="fft_interface")
         ok = DftiCommitDescriptor(fft%dfti_handle)
         if (ok /= 0) call juDFT_error("can't commit descriptor", calledby="fft_interface")
#endif
      case default 
         size_dat = product(length)
         allocate (fft%afft(size_dat), fft%bfft(size_dat), stat=ok)
         if (ok /= 0) call juDFT_error("can't alloc afft & bfft", calledby="fft_interface")
      end select
   end subroutine

   subroutine t_fft_exec(fft, dat)
      USE m_cfft
      implicit none 
      class(t_fft), intent(inout) :: fft
      complex, intent(inout)      :: dat(:) 
      integer :: isn, size_dat, ok

      size_dat = product(fft%length)

      select case(fft%backend)
      case(spFFT_const)
         call juDFT_error("not yet implemented")
      case(mklFFT_const)
#ifdef CPP_FFT_MKL
         if (fft%forw) then
            ok = DftiComputeForward(fft%dfti_handle, dat)
         else
            ok = DftiComputeBackward(fft%dfti_handle, dat)
         end if
#endif
      case default
         fft%afft = real(dat)
         fft%bfft = aimag(dat)
   
         isn = merge(-1, 1, fft%forw)
         CALL cfft(fft%afft, fft%bfft, size_dat, fft%length(1),fft% length(1), isn)
         CALL cfft(fft%afft, fft%bfft, size_dat, fft%length(2), fft%length(1)*fft%length(2), isn)
         CALL cfft(fft%afft, fft%bfft, size_dat, fft%length(3), size_dat, isn)
         dat = cmplx(fft%afft, fft%bfft)
      end select
   end subroutine t_fft_exec

   subroutine t_fft_free(fft)
      implicit none 
      class(t_fft) :: fft 
      integer :: ok

      fft%initialized = .False.
      fft%backend    = -1
      fft%length     = [-1,-1,-1]

      if(allocated(fft%afft)) deallocate(fft%afft)
      if(allocated(fft%bfft)) deallocate(fft%bfft)

      
      select case(fft%backend)
      case(spFFT_const)
         call juDFT_error("not yet implemented")
      case(mklFFT_const)
#ifdef CPP_FFT_MKL
         ok = DftiFreeDescriptor(fft%dfti_handle)
#endif
      case default
         
      end select
   end subroutine t_fft_free
end module m_types_fft
