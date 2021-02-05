module m_types_fft
   use m_judft
   private
   integer,public,parameter :: fft_c2c=1 !Different fft modes

   type,public:: t_fft
     logical :: initialized=.false.
     integer :: mode
     integer,allocatable:: length(:)
     logical :: forw
   contains
      procedure :: init => t_fft_init
      procedure :: exec => t_fft_exec_c  !more modes could be added, like c2r...
      procedure :: free => t_fft_free
   end type t_fft
contains
   subroutine t_fft_init(fft, length, forward,mode, indices)
      implicit none
      class(t_fft)                  :: fft
      integer, intent(in)           :: length(:) !length of data in each direction
      logical, intent(in)           :: forward          !.true. for the forward transformation, .false. for the backward one
      integer,intent(in),optional   :: mode
      INTEGER, OPTIONAL, INTENT(IN) :: indices(:)    !array of indices of relevant/nonzero elements in the FFT mesh

      fft%initialized=.true.
      if (present(mode)) THEN
        fft%mode=mode
      else
        fft%mode=fft_c2c
      endif
      fft%length=length
      fft%forw=forward

   end subroutine

   subroutine t_fft_exec_c(fft, dat)
      implicit none
      class(t_fft), intent(inout) :: fft
      complex, intent(inout)      :: dat(:)

      if (.not.fft%initialized) call judft_error("DFT not initialized")
      if (size(dat).ne.product(fft%length)) call judft_error("Wrong size of data in FFT")
      if (fft%mode.ne.fft_c2c) call judft_error("Wrong data for choosen FFT mode")
   end subroutine t_fft_exec_c

   subroutine t_fft_free(fft)
      implicit none
      class(t_fft) :: fft
      fft%initialized=.false.
   end subroutine t_fft_free
end module m_types_fft
