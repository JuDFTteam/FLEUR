module m_types_fft_singleton
   use m_types_fft
   private
   integer,public,parameter:: fft_singleton=1
   type,public,extends(t_fft):: t_fft_singleton

   contains
      procedure :: init => t_fft_init
      procedure :: exec => t_fft_exec_c
      procedure :: free => t_fft_free
   end type t_fft_singleton
contains
   subroutine t_fft_init(fft, length, forward,mode, indices)
      implicit none
      class(t_fft_singleton)        :: fft
      integer, intent(in)           :: length(:) !length of data in each direction
      logical, intent(in)           :: forward          !.true. for the forward transformation, .false. for the backward one
      integer,intent(in),optional   :: mode
      INTEGER, OPTIONAL, INTENT(IN) :: indices(:)    !array of indices of relevant/nonzero elements in the FFT mesh

      integer :: size_dat

      call fft%t_fft%init(length, forward,mode, indices)

      end subroutine

   subroutine t_fft_exec_c(fft, dat)
      USE m_cfft
      implicit none
      class(t_fft_singleton), intent(inout) :: fft
      complex, intent(inout)      :: dat(:)
      integer      :: isn, size_dat
      ! cfft storage
      real, allocatable :: afft(:), bfft(:)
      call fft%t_fft%exec(dat)
      size_dat=product(fft%length)
      afft = real(dat)
      bfft = aimag(dat)
      isn = merge(-1, 1, fft%forw)
      CALL cfft(afft, bfft, size_dat, fft%length(1),fft%length(1), isn)
      if (size(fft%length)>1) THEN
        CALL cfft(afft, bfft, size_dat, fft%length(2), fft%length(1)*fft%length(2), isn)
        if (size(fft%length)>2) THEN
          CALL cfft(afft, bfft, size_dat, fft%length(3), product(fft%length), isn)
        endif
      endif
      dat = cmplx(afft, bfft)
   end subroutine t_fft_exec_c

   subroutine t_fft_free(fft)
      implicit none
      class(t_fft_singleton) :: fft

      fft%initialized = .False.
   end subroutine t_fft_free
end module m_types_fft_singleton
