MODULE m_fft_interface
   USE m_juDFT

   IMPLICIT NONE

CONTAINS

   subroutine fft_interface(dimen, length, dat, forw, indices)
      use m_types_fft
      implicit none
      ! provides interfaces to fft subroutines

      integer, intent(in)           :: dimen         !dimension of fft transformation
      integer, intent(in)           :: length(dimen) !length of data in each direction
      complex, intent(inout)        :: dat(:)        !data to be transformed, size(dat) should be product(length)
      logical, intent(in)           :: forw          !.true. for the forward transformation, .false. for the backward one
      INTEGER, OPTIONAL, INTENT(IN) :: indices(:)    !array of indices of relevant/nonzero elements in the FFT mesh

      integer :: size_dat
      type(t_fft) :: fft

      size_dat = product(length)
      if (size(dat) .ne. size_dat) call juDFT_error('array bounds are inconsistent', calledby='fft_interface')
      if (dimen .ne. 3) call juDFT_error('sorry, not implemented yet for this value of dimen', calledby='fft_interface')

      call fft%init(length, forw, indices)
      call fft%exec(dat)
      call fft%free()

   end subroutine fft_interface

 

END MODULE m_fft_interface
