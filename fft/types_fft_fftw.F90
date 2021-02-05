module m_types_fft_fftw
  use m_types_fft
#ifndef CPP_FFTW
  INTEGER,PUBLIC,PARAMETER :: fft_fftw=-4
  type,PUBLIC,extends(t_fft):: t_fft_fftw
  end type
#else
  USE iso_c_binding
  use fftw3
  use m_juDFT
  PRIVATE
  INTEGER,PUBLIC,PARAMETER :: fft_fftw=4
  type,public,extends(t_fft):: t_fft_fftw
      ! cfft storage

      type(c_ptr) :: plan, ptr_in, ptr_out
      complex(C_DOUBLE_COMPLEX), pointer :: in(:), out(:)
   contains
      procedure :: init => t_fft_init
      procedure :: exec => t_fft_exec
      procedure :: free => t_fft_free
   end type t_fft_fftw
contains
   subroutine t_fft_init(fft, length, forward,mode, indices)
      implicit none
      class(t_fft_fftw)                  :: fft
      integer, intent(in)           :: length(:) !length of data in each direction
      logical, intent(in)           :: forward          !.true. for the forward transformation, .false. for the backward one
      INTEGER, OPTIONAL, INTENT(IN) :: mode,indices(:)    !array of indices of relevant/nonzero elements in the FFT mesh

      call fft%t_fft%init(length, forward,mode, indices)
      !$OMP critical
      fft%ptr_in = fftw_alloc_complex(int(product(length), C_SIZE_T))
      call c_f_pointer(fft%ptr_in, fft%in, [product(length)])
      fft%ptr_out = fftw_alloc_complex(int(product(length), C_SIZE_T))
      call c_f_pointer(fft%ptr_out, fft%out, [product(length)])
      select case(size(length))
      case(1)
        fft%plan = fftw_plan_dft_1d(fft%length(1),fft%in, fft%out, fft%forw,FFTW_MEASURE)
      case(2)
        fft%plan = fftw_plan_dft_2d(fft%length(2), fft%length(1),fft%in, fft%out, fft%forw,FFTW_MEASURE)
      case(3)
        fft%plan = fftw_plan_dft_3d(fft%length(3), fft%length(2), fft%length(1),&
        fft%in, fft%out, fft%forw,FFTW_MEASURE)
      case default
        call judft_error("fftw interface only for 1-3 dimensions")
      end select
      !$OMP end critical
   end subroutine

   subroutine t_fft_exec(fft, dat)
      USE m_cfft
      implicit none
      class(t_fft_fftw), intent(inout) :: fft
      complex, intent(inout)      :: dat(:)

      call fft%t_fft%exec(dat)
      fft%in = dat
      call fftw_execute_dft(fft%plan, fft%in, fft%out)
      dat = fft%out
   end subroutine t_fft_exec

   subroutine t_fft_free(fft)
      implicit none
      integer      :: ok
      class(t_fft_fftw) :: fft

      !$OMP critical
      call fftw_destroy_plan(fft%plan)
      call fftw_free(fft%ptr_in)
      call fftw_free(fft%ptr_out)
      !$OMP end critical
      fft%plan  = c_null_ptr

      fft%initialized = .False.
   end subroutine t_fft_free
#endif
end module m_types_fft_fftw
