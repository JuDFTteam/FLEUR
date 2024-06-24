MODULE m_fft_interface
   USE m_juDFT

   IMPLICIT NONE

CONTAINS

   subroutine fft_interface(dimen, length, dat, forw, indices)
      use m_types_fft
#ifdef _OPENACC
      use openacc 
#endif      
      implicit none
      ! provides interfaces to fft subroutines

      integer, intent(in)           :: dimen         !dimension of fft transformation
      integer, intent(in)           :: length(dimen) !length of data in each direction
      complex, intent(inout)        :: dat(:)        !data to be transformed, size(dat) should be product(length)
      logical, intent(in)           :: forw          !.true. for the forward transformation, .false. for the backward one
      INTEGER, OPTIONAL, INTENT(IN) :: indices(:)    !array of indices of relevant/nonzero elements in the FFT mesh

      integer :: size_dat
      logical :: l_gpu
      type(t_fft) :: fft
      call timestart("FFT_interface")
      size_dat = product(length)
      if (size(dat) .ne. size_dat) call juDFT_error('array bounds are inconsistent', calledby='fft_interface')
      if (dimen .ne. 3) call juDFT_error('sorry, not implemented yet for this value of dimen', calledby='fft_interface')

#ifdef _OPENACC
      l_gpu=acc_is_present(dat)
#else
      l_gpu=.false.
#endif      

      call fft%init(length, forw, indices,l_gpu=l_gpu)
      call fft%exec(dat)
      call fft%free()
        
      call timestop("FFT_interface")
   end subroutine fft_interface

 

END MODULE m_fft_interface
