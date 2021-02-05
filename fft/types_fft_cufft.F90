module m_types_fft_cufft
  use m_types_fft
  PRIVATE
#ifndef CPP_CuFFT
  INTEGER,PUBLIC,PARAMETER :: fft_cufft=-6
  type,PUBLIC,extends(t_fft):: t_fft_cufft
  end type
#else
  use OPENACC
  use cufft
  INTEGER,PUBLIC,PARAMETER :: fft_cufft=6
  type,PUBLIC,extends(t_fft):: t_fft_cufft
    integer :: plan
    complex,ALLOCATABLE:: out(:)
   contains
      procedure :: init => t_fft_init
      procedure :: exec => t_fft_exec
      procedure :: free => t_fft_free
   end type t_fft_cufft
contains
   subroutine t_fft_init(fft, length, forward,mode, indices)
      implicit none
      class(t_fft_cufft)                  :: fft
      integer, intent(in)           :: length(:) !length of data in each direction
      logical, intent(in)           :: forward          !.true. for the forward transformation, .false. for the backward one
      INTEGER, OPTIONAL, INTENT(IN) :: mode,indices(:)    !array of indices of relevant/nonzero elements in the FFT mesh

      integer:: ierr
      call fft%t_fft%init(length, forward,mode, indices)

      select case(size(length))
      case(1)
        ierr=cufftPlan1D(fft%plan,length(1),CUFFT_Z2Z)
      case(2)
        ierr=cufftPlan2D(fft%plan,length(1),length(2),CUFFT_Z2Z)
      case(3)
        ierr=cufftPlan3D(fft%plan,length(1),length(2),length(3),CUFFT_Z2Z)
      end select

      allocate(fft%out(product(length)))
      !$acc enter data create(fft,fft%out)
      if (ierr/=0) call judft_error("CuFFT plan failed")
    end subroutine

   subroutine t_fft_exec(fft, dat)
      USE m_cfft
      implicit none
      class(t_fft_cufft), intent(inout) :: fft
      complex, intent(inout)      :: dat(:)
      integer                     :: ierr

      call fft%t_fft%exec(dat)
      if (.not.acc_is_present(dat)) call judft_error("FFT-data must be present on GPU for cuFFT")

      direction=merge(CUFFT_FORWARD,CUFFT_INVERSE,fft%forw)
      !$acc host_data use_device(dat,fft%out)
      ierr=cufftExecZ2Z(fft%plan,dat,fft%out,direction)
      !$acc end host_data
      !$acc kernels present(dat,fft,fft%out)
      dat=fft%out
      !$acc end kernels
      if (ierr/=0) call judft_error("CuFFT failed")
   end subroutine t_fft_exec

   subroutine t_fft_free(fft)
      implicit none
      integer      :: ierr
      class(t_fft_cufft) :: fft

      !$acc exit data delete(fft%out,fft)
      deallocate(fft%out)

      ierr=cufftDestroy(fft%plan)
      if (ierr/=0) call judft_error("cuFFTDestroy failed")
      fft%initialized = .False.
   end subroutine t_fft_free
#endif
end module m_types_fft_cufft
