module m_types_fft
#ifdef CPP_FFT_MKL
   USE mkl_dfti
#endif
   USE m_selecFFT
   use m_judft
   USE iso_c_binding
#ifdef CPP_SPFFT
   USE spfft
#endif
#ifdef CPP_FFTW 
   use fftw3
#endif
   !$ use omp_lib

#ifdef CPP_FFT_MKL
   type ptr_container
      type(dfti_descriptor), pointer :: dfti_handle
   endtype ptr_container
#endif

   type t_fft 
      logical :: initialized = .False.
      integer :: backend = -1
      integer :: batch_size = 1
      integer :: length(3) = [-1,-1,-1]
      logical :: forw
      ! cfft storage
      real, allocatable :: afft(:), bfft(:)
#ifdef CPP_FFT_MKL
      ! mkl
      type(ptr_container), allocatable :: container(:)
#endif
#ifdef CPP_SPFFT
      !SpFFT
      integer, allocatable :: indices(:)
      type(c_ptr)          :: transform = c_null_ptr, realSpacePtr = c_null_ptr
      integer              :: xyPlanesize
      COMPLEX(C_DOUBLE_COMPLEX), ALLOCATABLE :: recSpaceFunction(:)
      COMPLEX(C_DOUBLE_COMPLEX), POINTER     :: externalRealSpaceMesh(:, :, :)
#endif
#ifdef CPP_FFTW
      type(c_ptr), allocatable :: plan(:)
      type(c_ptr)              :: ptr_in, ptr_out
      complex(C_DOUBLE_COMPLEX), pointer :: in(:,:), out(:,:)       
#endif
   contains 
      procedure :: init => t_fft_init
      procedure :: exec => t_fft_exec_single
      procedure :: exec_batch => t_fft_exec_batched
      procedure :: free => t_fft_free
   end type t_fft
contains
   subroutine t_fft_init(fft, length, forw, indices, batch_size)
      implicit none       
      class(t_fft)                  :: fft
      integer, intent(in)           :: length(3) !length of data in each direction
      logical, intent(in)           :: forw          !.true. for the forward transformation, .false. for the backward one
      INTEGER, OPTIONAL, INTENT(IN) :: indices(:)    !array of indices of relevant/nonzero elements in the FFT mesh
      integer, optional, intent(in) :: batch_size

      INTEGER, PARAMETER :: numOMPThreads = 1
      integer :: size_dat, ok, fftMeshIndex, maxNumLocalZColumns
      integer :: temp, x, y, z, xCoord, yCoord, zCoord, i
      INTEGER, ALLOCATABLE :: sparseCoords(:)
      LOGICAL, ALLOCATABLE :: nonzeroArea(:, :)
      type(c_ptr)          :: grid = c_null_ptr
      logical :: in_openmp = .false.
      integer :: max_threads = 1, thread_id = 0
      integer :: n_plans

      !$ thread_id   = omp_get_thread_num()
      !$ max_threads = omp_get_max_threads()
      !$ in_openmp   = omp_in_parallel()

      if(present(batch_size)) then
         fft%batch_size = batch_size
      else
         fft%batch_size = 1
      endif

      if(in_openmp .and. fft%batch_size > 1) then
         call juDFT_error("Can't create FFTs inside OpenMP and use batch at the same time")
      endif

      fft%initialized = .True.
      fft%backend = defaultFFT_const
      fft%backend = selecFFT(PRESENT(indices))
      fft%length  = length
      fft%forw    = forw

      select case(fft%backend)
#ifdef CPP_FFTW
      case(FFTW_const)
         n_plans = min(max_threads, batch_size)
         allocate(fft%plan(n_plans))

         fft%ptr_in = fftw_alloc_complex(int(n_plans * product(length), C_SIZE_T))
         call c_f_pointer(fft%ptr_in, fft%in, [product(length), n_plans])

         fft%ptr_out = fftw_alloc_complex(int(n_plans * product(length), C_SIZE_T))
         call c_f_pointer(fft%ptr_out, fft%out, [product(length), n_plans])

         do i = 1,n_plans
            !$omp critical
            if(fft%forw) then
               fft%plan(i) = fftw_plan_dft_3d(fft%length(3), fft%length(2), fft%length(1),&
                                             fft%in(:,i), fft%out(:,i), FFTW_FORWARD,FFTW_MEASURE) 
            else
               fft%plan(i) = fftw_plan_dft_3d(fft%length(3), fft%length(2), fft%length(1),&
                                             fft%in(:,i), fft%out(:,i), FFTW_BACKWARD,FFTW_MEASURE) 
            endif
            !$omp end critical
         enddo
#endif
      case(mklFFT_const)
#ifdef CPP_FFT_MKL
         n_plans = min(max_threads, batch_size)
         allocate(fft%container(n_plans))
         do i = 1,n_plans
            ok = DftiCreateDescriptor(fft%container(i)%dfti_handle, dfti_double, dfti_complex, 3, length)
            if (ok /= 0) call juDFT_error("cant create descriptor", calledby="fft_interface")
            ok = DftiCommitDescriptor(fft%container(i)%dfti_handle)
            if (ok /= 0) call juDFT_error("can't commit descriptor", calledby="fft_interface")
         enddo
#endif

#ifdef CPP_SPFFT
      case(spFFT_const)
         fft%indices = indices
         ALLOCATE(sparseCoords(3*SIZE(fft%indices)))
         if(.not. allocated(fft%recSpaceFunction)) ALLOCATE(fft%recSpaceFunction(SIZE(fft%indices)))
         ALLOCATE(nonzeroArea(0:length(1) - 1, 0:length(2) - 1))
         nonzeroArea(:, :) = .FALSE.
         fft%xyPlaneSize = fft%length(1)*fft%length(2)
         DO i = 1, SIZE(fft%indices)
            zCoord = fft%indices(i)/fft%xyPlaneSize
            temp = MOD(fft%indices(i), fft%xyPlaneSize)
            yCoord = temp/length(1)
            xCoord = MOD(temp, length(1))

            sparseCoords(3*(i - 1) + 3) = zCoord
            sparseCoords(3*(i - 1) + 2) = yCoord
            sparseCoords(3*(i - 1) + 1) = xCoord

            nonzeroArea(xCoord, yCoord) = .TRUE.
         END DO

         maxNumLocalZColumns = COUNT(nonzeroArea)
         IF (fft%forw) THEN
            ok = spfft_grid_create(grid, length(1), length(2), length(3), &
                                          maxNumLocalZColumns, SPFFT_PU_HOST, numOMPThreads); 
            IF (ok /= SPFFT_SUCCESS) CALL juDFT_error("Error in creating spFFT grid! (1)")

            ok = spfft_transform_create(fft%transform, grid, SPFFT_PU_HOST, SPFFT_TRANS_C2C, &
                                               length(1), length(2), length(3), length(3), &
                                               size(fft%recSpaceFunction), SPFFT_INDEX_TRIPLETS, sparseCoords)
            IF (ok /= SPFFT_SUCCESS) CALL juDFT_error("Error in creating spFFT transform! (1)")

            ok = spfft_grid_destroy(grid)
            IF (ok /= SPFFT_SUCCESS) CALL juDFT_error("Error in destroying spFFT grid! (1)")

            ok = spfft_transform_get_space_domain(fft%transform, SPFFT_PU_HOST, fft%realSpacePtr)
            IF (ok /= SPFFT_SUCCESS) CALL juDFT_error("Error in obtaining spFFT space domain! (1)")

            CALL C_F_POINTER(fft%realSpacePtr, fft%externalRealSpaceMesh, [length(1), length(2), length(3)])
         ELSE
            ok = spfft_grid_create(grid, length(1), length(2), length(3), &
                                          maxNumLocalZColumns, SPFFT_PU_HOST, numOMPThreads); 
            IF (ok /= SPFFT_SUCCESS) CALL juDFT_error("Error in creating spFFT grid! (2)")

            ok = spfft_transform_create(fft%transform, grid, SPFFT_PU_HOST, SPFFT_TRANS_C2C, &
                                               length(1), length(2), length(3), length(3), &
                                               size(fft%recSpaceFunction), SPFFT_INDEX_TRIPLETS, sparseCoords)
            IF (ok /= SPFFT_SUCCESS) CALL juDFT_error("Error in creating spFFT transform! (2)")

            ok = spfft_grid_destroy(grid)
            IF (ok /= SPFFT_SUCCESS) CALL juDFT_error("Error in destroying spFFT grid! (2)")

            ok = spfft_transform_get_space_domain(fft%transform, SPFFT_PU_HOST, fft%realSpacePtr)
            IF (ok /= SPFFT_SUCCESS) CALL juDFT_error("Error in obtaining spFFT space domain! (2)")
         END IF
#endif
      case default 
         size_dat = product(length)
         allocate (fft%afft(size_dat), fft%bfft(size_dat), stat=ok)
         if (ok /= 0) call juDFT_error("can't alloc afft & bfft", calledby="fft_interface")
      end select
   end subroutine


   subroutine t_fft_exec_batched(fft, dat)
      USE m_cfft
      implicit none 
      class(t_fft), intent(inout) :: fft
      complex, intent(inout)      :: dat(:,:) 
      integer      :: isn, size_dat 
      INTEGER      ::  i, x, y, z, fftMeshIndex, ok, me

      size_dat = product(fft%length)

      select case(fft%backend)
#ifdef CPP_FFTW
      case(fftw_const)
         me = 1
         !$omp parallel do default(none) private(me, i) shared(fft, dat)
         do i = 1,size(dat,2)
            !$ me = omp_get_thread_num() + 1
            fft%in(:,me) = dat(:,i)
            call fftw_execute_dft(fft%plan(me), fft%in(:,me), fft%out(:,me))
            dat(:,i) = fft%out(:,me)
         enddo
         !$omp end parallel do
#endif
      case(mklFFT_const)
#ifdef CPP_FFT_MKL
         me = 1
         !$omp parallel do default(none) private(me, i, ok) shared(fft, dat)
         do i = 1,size(dat,2)
            !$ me = omp_get_thread_num() + 1
            if (fft%forw) then
               ok = DftiComputeForward(fft%container(me)%dfti_handle, dat(:,i))
            else
               ok = DftiComputeBackward(fft%container(me)%dfti_handle, dat(:,i))
            end if
            if(ok /= 0) juDFT_error("problem executing dft")
         enddo
         !$omp end parallel do
#endif
#ifdef CPP_SPFFT
      case(spFFT_const)
         IF (fft%forw) THEN
            DO z = 1, SIZE(fft%externalRealSpaceMesh, 3)
               DO y = 1, SIZE(fft%externalRealSpaceMesh, 2)
                  DO x = 1, SIZE(fft%externalRealSpaceMesh, 1)
                     fftMeshIndex = (x - 1) + (y - 1)*fft%length(1) + (z - 1)*fft%xyPlaneSize + 1
                     fft%externalRealSpaceMesh(x, y, z) = dat(fftMeshIndex)
                  END DO
               END DO
            END DO
            ok = spfft_transform_forward(fft%transform, SPFFT_PU_HOST, fft%recSpaceFunction, SPFFT_NO_SCALING)!SPFFT_FULL_SCALING)
            IF (ok /= SPFFT_SUCCESS) THEN
               CALL juDFT_error("Error in spFFT forward fft%transform! (1)", calledby="fft_interface")
            END IF
            dat(:) = CMPLX(0.0, 0.0)
            DO i = 1, SIZE(fft%indices)
               dat(fft%indices(i) + 1) = fft%recSpaceFunction(i)
            END DO

         ELSE
            DO i = 1, SIZE(fft%indices)
               fft%recSpaceFunction(i) = dat(fft%indices(i) + 1)
            END DO
            ok = spfft_transform_backward(fft%transform, fft%recSpaceFunction, SPFFT_PU_HOST)
            IF (ok /= SPFFT_SUCCESS) THEN
               CALL juDFT_error("Error in spFFT backward fft%transform! (2)", calledby="fft_interface")
            END IF

            CALL C_F_POINTER(fft%realSpacePtr, fft%externalRealSpaceMesh, [fft%length(1), fft%length(2), fft%length(3)])

            DO z = 1, SIZE(fft%externalRealSpaceMesh, 3)
               DO y = 1, SIZE(fft%externalRealSpaceMesh, 2)
                  DO x = 1, SIZE(fft%externalRealSpaceMesh, 1)
                     fftMeshIndex = (x - 1) + (y - 1)*fft%length(1) + (z - 1)*fft%xyPlaneSize + 1
                     dat(fftMeshIndex) = fft%externalRealSpaceMesh(x, y, z)
                  END DO
               END DO
            END DO
         END IF
#endif
      case default
         do i = 1,size(dat,2)
            fft%afft = real(dat(:,i))
            fft%bfft = aimag(dat(:,i))
      
            isn = merge(-1, 1, fft%forw)
            CALL cfft(fft%afft, fft%bfft, size_dat, fft%length(1),fft% length(1), isn)
            CALL cfft(fft%afft, fft%bfft, size_dat, fft%length(2), fft%length(1)*fft%length(2), isn)
            CALL cfft(fft%afft, fft%bfft, size_dat, fft%length(3), size_dat, isn)
            dat(:,i) = cmplx(fft%afft, fft%bfft)
         enddo
      end select
   end subroutine t_fft_exec_batched

   subroutine t_fft_exec_single(fft, dat)
      implicit none
      class(t_fft), intent(inout)   :: fft
      complex, intent(inout),target :: dat(:) 
      integer :: i
      type(c_ptr)      :: ptr
      complex, pointer :: tmp_2d(:,:) 

      ! if an array is 1D just pretend it's 2d
      ptr = c_loc(dat)
      call c_f_pointer(ptr, tmp_2d, [size(dat), 1])

      call t_fft_exec_batched(fft, tmp_2d)
   end subroutine t_fft_exec_single

   subroutine t_fft_free(fft)
      implicit none 
      integer      :: ok
      class(t_fft) :: fft 
      logical :: in_openmp = .false.
      integer :: i
      !$ in_openmp   = omp_in_parallel()

      if(allocated(fft%afft)) deallocate(fft%afft)
      if(allocated(fft%bfft)) deallocate(fft%bfft)
      select case(fft%backend)
#ifdef CPP_FFTW
      case(FFTW_const)
         call fftw_free(fft%ptr_in)
         call fftw_free(fft%ptr_out)

         do i=1,size(fft%plan)
            !$omp critical
            call fftw_destroy_plan(fft%plan(i))
            !$omp end critical
            fft%plan(i)  = c_null_ptr
         enddo     
         deallocate(fft%plan)
#endif
      case(mklFFT_const)
#ifdef CPP_FFT_MKL
         do i=1,size(fft%container)
            ok = DftiFreeDescriptor(fft%container(i)%dfti_handle)
         enddo
         deallocate(fft%container)
#endif
#ifdef CPP_SPFFT
      case(spFFT_const)
         ok = spfft_transform_destroy(fft%transform)
         IF (ok /= SPFFT_SUCCESS) CALL juDFT_error("Error in destroying spFFT fft%transform! (1)")
         fft%transform    = c_null_ptr
         fft%realSpacePtr = c_null_ptr
#endif
      case default
         
      end select

      fft%initialized = .False.
      fft%backend    = -1
      fft%length     = [-1,-1,-1]
      fft%batch_size = -1
   end subroutine t_fft_free
end module m_types_fft
