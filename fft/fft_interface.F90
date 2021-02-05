MODULE m_fft_interface

   USE m_types_fft
   USE m_selectFFT
   IMPLICIT NONE

CONTAINS

   subroutine fft_interface( length, dat, forw, indices)
      use m_types_fft
      implicit none
      ! provides interfaces to fft subroutines

      integer, intent(in)           :: length(:) !length of data in each direction
      complex, intent(inout)        :: dat(:)        !data to be transformed, size(dat) should be product(length)
      logical, intent(in)           :: forw          !.true. for the forward transformation, .false. for the backward one
      INTEGER, OPTIONAL, INTENT(IN) :: indices(:)    !array of indices of relevant/nonzero elements in the FFT mesh

      class(t_fft),ALLOCATABLE :: fft

      fft=selectFFT()

      call fft%init(length, forw,fft_c2c, indices)
      call fft%exec(dat)
      call fft%free()


   end subroutine fft_interface

END MODULE m_fft_interface
