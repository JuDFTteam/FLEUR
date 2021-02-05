module m_types_fft_spfft
  use m_types_fft
#ifndef CPP_FFT_SPFFT
    integer,PARAMETER ::fft_spfft=-3
    type,extends(t_fft):: t_fft_spfft
    end type
#else
   use m_judft
   USE iso_c_binding
   USE spfft
   integer,PARAMETER :: fft_spfft=3
   type,extends(t_fft):: t_fft_spfft
      !SpFFT
      integer, allocatable :: indices(:)
      type(c_ptr)          :: transform = c_null_ptr, realSpacePtr = c_null_ptr
      integer              :: xyPlanesize
      COMPLEX(C_DOUBLE_COMPLEX), ALLOCATABLE :: recSpaceFunction(:)
      COMPLEX(C_DOUBLE_COMPLEX), POINTER     :: externalRealSpaceMesh(:, :, :)
   contains
      procedure :: init => t_fft_init
      procedure :: exec => t_fft_exec_c
      procedure :: free => t_fft_free
   end type t_fft
contains
   subroutine t_fft_init(fft, length, forward,mode, indices)
      implicit none
      class(t_fft_spfft)                  :: fft
      integer, intent(in)           :: length(:) !length of data in each direction
      logical, intent(in)           :: forward          !.true. for the forward transformation, .false. for the backward one
      INTEGER, OPTIONAL, INTENT(IN) :: mode,indices(:)    !array of indices of relevant/nonzero elements in the FFT mesh

      INTEGER, PARAMETER :: numOMPThreads = 1
      integer :: size_dat, ok, fftMeshIndex, maxNumLocalZColumns
      integer :: temp, x, y, z, xCoord, yCoord, zCoord, i
      INTEGER, ALLOCATABLE :: sparseCoords(:)
      LOGICAL, ALLOCATABLE :: nonzeroArea(:, :)
      type(c_ptr)          :: grid = c_null_ptr

      call fft%t_fft%init(fft, length, forward,mode, indices)

      if (size(length).ne.3) call judft_error("Only 3d FFT in SpFFT interface")

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
    end subroutine

    subroutine t_fft_exec_c(fft, dat)
      USE m_cfft
      implicit none
      class(t_fft_spfft), intent(inout) :: fft
      complex, intent(inout)      :: dat(:)
      integer      :: isn, size_dat
      INTEGER      ::  i, x, y, z, fftMeshIndex, ok

      call fft%t_fft%execute(dat)
      size_dat = product(fft%length)

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
   end subroutine t_fft_exec_c

   subroutine t_fft_free(fft)
      implicit none
      integer      :: ok
      class(t_fft_spfft) :: fft

      call fft%t_fft%free()
      ok = spfft_transform_destroy(fft%transform)
      IF (ok /= SPFFT_SUCCESS) CALL juDFT_error("Error in destroying spFFT fft%transform! (1)")
      fft%transform    = c_null_ptr
      fft%realSpacePtr = c_null_ptr
    end subroutine t_fft_free
#endif
end module m_types_fft_spfft
