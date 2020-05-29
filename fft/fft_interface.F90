MODULE m_fft_interface

   USE m_cfft
   USE m_juDFT
   USE m_selecFFT
#ifdef CPP_FFT_MKL
   USE mkl_dfti
#endif
#ifdef CPP_SPFFT
   USE iso_c_binding
   USE spfft
#endif

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

   contains
      subroutine cfft_wrapper(length, dat, forw)
         use m_juDFT
         implicit none
         integer, intent(in)           :: length(3)
         complex, intent(inout)        :: dat(:)
         logical, intent(in)           :: forw

         real, allocatable :: afft(:), bfft(:)
         integer :: isn, size_dat, ok

         size_dat = product(length)
         allocate (afft(size_dat), bfft(size_dat), stat=ok)
         if (ok /= 0) call juDFT_error("can't alloc afft & bfft", calledby="fft_interface")

         afft = real(dat)
         bfft = aimag(dat)

         isn = merge(-1, 1, forw)
         CALL cfft(afft, bfft, size_dat, length(1), length(1), isn)
         CALL cfft(afft, bfft, size_dat, length(2), length(1)*length(2), isn)
         CALL cfft(afft, bfft, size_dat, length(3), size_dat, isn)

         dat = cmplx(afft, bfft)
      end subroutine
   end subroutine fft_interface

   function g2fft(fft_size, g_in) result(g_idx)
      implicit none 
      integer, intent(in)        :: g_in(3), fft_size(3)
      integer                    :: g_idx
  
      integer ::  shifted_g(3)
  
      ! the fft-grid starts at g=0, not -g_max
      ! therefore all negative g's need to be shifted
      
      shifted_g = merge(g_in+fft_size,g_in, g_in<0)
  
      ! map it to 1d
      g_idx = shifted_g(1) &
            + shifted_g(2) * fft_size(1) &
            + shifted_g(3) * fft_size(1) * fft_size(2)
    end function g2fft

END MODULE m_fft_interface
