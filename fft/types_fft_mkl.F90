module m_types_fft_mkl
   use m_types_fft
#ifndef CPP_FFT_MKL
   integer,parameter:: fft_mkl=-2
   type,extends(t_fft):: t_fft_mkl
   end type
#else
   USE mkl_dfti
   USE iso_c_binding
   integer,parameter:: fft_mkl=2

   type,extends(t_fft):: t_fft_mkl
      type(dfti_descriptor), pointer :: dfti_handle
   contains
      procedure :: init => t_fft_init
      procedure :: exec => t_fft_exec
      procedure :: free => t_fft_free
   end type t_fft
contains
   subroutine t_fft_init(fft, length, forward,mode, indices)
      implicit none
      class(t_fft_mkl)                  :: fft
      integer, intent(in)           :: length(:) !length of data in each direction
      logical, intent(in)           :: forward          !.true. for the forward transformation, .false. for the backward one
      INTEGER, OPTIONAL, INTENT(IN) :: mode,indices(:)    !array of indices of relevant/nonzero elements in the FFT mesh

      integer ::  ok

      fft%initialized = .True.
      fft%length  = product(length)
      fft%forw    = forward


      ok = DftiCreateDescriptor(fft%dfti_handle, dfti_double, dfti_complex, size(length), length)

      if (ok /= 0) call juDFT_error("cant create descriptor", calledby="types_fft_mkl")
      ok = DftiCommitDescriptor(fft%dfti_handle)
      if (ok /= 0) call juDFT_error("can't commit descriptor", calledby="types_fft_mkl")
   end subroutine

   subroutine t_fft_exec_c(fft, dat)
      USE m_cfft
      implicit none
      class(t_fft_mkl), intent(inout) :: fft
      complex, intent(inout)      :: dat(:)
      integer :: ok

      if (size(dat)/=fft%length) call judft_error("Wrong dimension in FFT")

      if (fft%forw) then
        ok = DftiComputeForward(fft%dfti_handle, dat)
      else
        ok = DftiComputeBackward(fft%dfti_handle, dat)
      end if
      if (ok /= 0 ) call judft_error("MKL FFT failed")
   end subroutine t_fft_exec_c

   subroutine t_fft_free(fft)
     implicit none
     integer      :: ok
     class(t_fft_mkl) :: fft

     ok = DftiFreeDescriptor(fft%dfti_handle)

      fft%initialized = .False.
    end subroutine t_fft_free
#endif
end module m_types_fft_mkl
