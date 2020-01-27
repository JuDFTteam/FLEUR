MODULE m_fft_interface

   USE m_cfft
   USE m_juDFT
   USE m_selectFFT
#ifdef CPP_FFT_MKL        
   USE mkl_dfti
#endif
#ifdef CPP_SPFFT
   USE iso_c_binding
   USE spfft
#endif

IMPLICIT NONE

CONTAINS

    subroutine fft_interface(dimen,length,dat,forw,indices)
    ! provides interfaces to fft subroutines

    integer, intent(in)           :: dimen         !dimension of fft transformation
    integer, intent(in)           :: length(dimen) !length of data in each direction
    complex, intent(inout)        :: dat(:)        !data to be transformed, size(dat) should be sum(length)
    logical, intent(in)           :: forw          !.true. for the forward transformation, .false. for the backward one
    INTEGER, OPTIONAL, INTENT(IN) :: indices(:)    !array of indices of relevant/nonzero elements in the FFT mesh

#ifdef CPP_FFT_MKL
    type(dfti_descriptor),pointer :: dfti_handle
    integer :: dfti_status
#endif
#ifdef CPP_SPFFT
    INTEGER,                   ALLOCATABLE :: sparseCoords(:)
    COMPLEX(C_DOUBLE_COMPLEX), ALLOCATABLE :: recSpaceFunction(:)
    type(c_ptr)                            :: grid = c_null_ptr
    type(c_ptr)                            :: transform = c_null_ptr
    type(c_ptr)                            :: realSpacePtr
    INTEGER, PARAMETER                     :: numOMPThreads = -1 ! -1 gives you the default number of OMP threads
    LOGICAL, ALLOCATABLE                   :: nonzeroArea(:,:)
    INTEGER                                :: xCoord, yCoord, zCoord, maxNumLocalZColumns
    INTEGER                                :: errorCode, x, y, z, fftMeshIndex
    COMPLEX(C_DOUBLE_COMPLEX), POINTER     :: externalRealSpaceMesh(:,:,:)
#endif

    ! default variables
    real,allocatable :: afft(:),bfft(:)
    integer :: isn

    integer :: size_dat,i
    INTEGER :: fftRoutine, xyPlaneSize, temp
    LOGICAL :: l_sparse

       size_dat = 1
       do i = 1,dimen
          size_dat = size_dat * length(i)
       enddo
       if (size(dat) .ne. size_dat) call juDFT_error('array bounds are inconsistent',calledby ='fft_interface')
       if (dimen .ne. 3 )  call juDFT_error('sorry, not implemented yet for this value of dimen',calledby ='fft_interface')

       l_sparse = PRESENT(indices)

       fftRoutine = defaultFFT_const
       fftRoutine = selectFFT(l_sparse)

       IF(fftRoutine.EQ.spFFT_const) THEN
#ifdef CPP_SPFFT
          ALLOCATE (sparseCoords(3*SIZE(indices)))
          ALLOCATE (recSpaceFunction(SIZE(indices)))
          ALLOCATE (nonzeroArea(0:length(1)-1,0:length(2)-1))
          nonzeroArea(:,:) = .FALSE.
          xyPlaneSize = length(1)*length(2)
          DO i = 1, SIZE(indices)
             zCoord = indices(i) / xyPlaneSize
             temp = MOD(indices(i),xyPlaneSize)
             yCoord = temp/length(1)
             xCoord = MOD(temp,length(1))

             sparseCoords(3*(i-1)+3) = zCoord
             sparseCoords(3*(i-1)+2) = yCoord
             sparseCoords(3*(i-1)+1) = xCoord

             nonzeroArea(xCoord,yCoord) = .TRUE.
          END DO

          maxNumLocalZColumns = COUNT(nonzeroArea)

          IF(forw) THEN
             errorCode = spfft_grid_create(grid, length(1), length(2), length(3), &
                                           maxNumLocalZColumns, SPFFT_PU_HOST, numOMPThreads);
             IF (errorCode.NE.SPFFT_SUCCESS) THEN
                CALL juDFT_error("Error in creating spFFT grid! (1)", calledby="fft_interface")
             END IF

             errorCode = spfft_transform_create(transform, grid, SPFFT_PU_HOST, SPFFT_TRANS_C2C, &
                                                length(1), length(2), length(3), length(3), &
                                                size(recSpaceFunction), SPFFT_INDEX_TRIPLETS, sparseCoords)
             IF (errorCode.NE.SPFFT_SUCCESS) THEN
                CALL juDFT_error("Error in creating spFFT transform! (1)", calledby="fft_interface")
             END IF

             errorCode = spfft_grid_destroy(grid)
             IF (errorCode.NE.SPFFT_SUCCESS) THEN
                CALL juDFT_error("Error in destroying spFFT grid! (1)", calledby="fft_interface")
             END IF

             errorCode = spfft_transform_get_space_domain(transform, SPFFT_PU_HOST, realSpacePtr)
             IF (errorCode.NE.SPFFT_SUCCESS) THEN
                CALL juDFT_error("Error in obtaining spFFT space domain! (1)", calledby="fft_interface")
             END IF

             CALL C_F_POINTER(realSpacePtr, externalRealSpaceMesh, [length(1),length(2),length(3)])

             DO z = 1, SIZE(externalRealSpaceMesh, 3)
                DO y = 1, SIZE(externalRealSpaceMesh, 2)
                   DO x = 1, SIZE(externalRealSpaceMesh, 1)
                      fftMeshIndex = (x-1) + (y-1)*length(1) + (z-1)*xyPlaneSize + 1
                      externalRealSpaceMesh(x,y,z) = dat(fftMeshIndex)
                   END DO
                END DO
             END DO

             errorCode = spfft_transform_forward(transform, SPFFT_PU_HOST, recSpaceFunction, SPFFT_NO_SCALING)!SPFFT_FULL_SCALING)
             IF (errorCode.NE.SPFFT_SUCCESS) THEN
                CALL juDFT_error("Error in spFFT forward transform! (1)", calledby="fft_interface")
             END IF

             errorCode = spfft_transform_destroy(transform)
             IF (errorCode.NE.SPFFT_SUCCESS) THEN
                CALL juDFT_error("Error in destroying spFFT transform! (1)", calledby="fft_interface")
             END IF

             dat(:) = CMPLX(0.0,0.0)
             DO i = 1, SIZE(indices)
                dat(indices(i)+1) = recSpaceFunction(i)
             END DO


          ELSE
             DO i = 1, SIZE(indices)
                recSpaceFunction(i) = dat(indices(i)+1)
             END DO

             errorCode = spfft_grid_create(grid, length(1), length(2), length(3), &
                                           maxNumLocalZColumns, SPFFT_PU_HOST, numOMPThreads);
             IF (errorCode.NE.SPFFT_SUCCESS) THEN
                CALL juDFT_error("Error in creating spFFT grid! (2)", calledby="fft_interface")
             END IF

             errorCode = spfft_transform_create(transform, grid, SPFFT_PU_HOST, SPFFT_TRANS_C2C, &
                                                length(1), length(2), length(3), length(3), &
                                                size(recSpaceFunction), SPFFT_INDEX_TRIPLETS, sparseCoords)
             IF (errorCode.NE.SPFFT_SUCCESS) THEN
                CALL juDFT_error("Error in creating spFFT transform! (2)", calledby="fft_interface")
             END IF

             errorCode = spfft_grid_destroy(grid)
             IF (errorCode.NE.SPFFT_SUCCESS) THEN
                CALL juDFT_error("Error in destroying spFFT grid! (2)", calledby="fft_interface")
             END IF

             errorCode = spfft_transform_get_space_domain(transform, SPFFT_PU_HOST, realSpacePtr)
             IF (errorCode.NE.SPFFT_SUCCESS) THEN
                CALL juDFT_error("Error in obtaining spFFT space domain! (2)", calledby="fft_interface")
             END IF

             errorCode = spfft_transform_backward(transform, recSpaceFunction, SPFFT_PU_HOST)
             IF (errorCode.NE.SPFFT_SUCCESS) THEN
                CALL juDFT_error("Error in spFFT backward transform! (2)", calledby="fft_interface")
             END IF

             CALL C_F_POINTER(realSpacePtr, externalRealSpaceMesh, [length(1),length(2),length(3)])

             DO z = 1, SIZE(externalRealSpaceMesh, 3)
                DO y = 1, SIZE(externalRealSpaceMesh, 2)
                   DO x = 1, SIZE(externalRealSpaceMesh, 1)
                      fftMeshIndex = (x-1) + (y-1)*length(1) + (z-1)*xyPlaneSize + 1
                      dat(fftMeshIndex) = externalRealSpaceMesh(x,y,z)
                   END DO
                END DO
             END DO

             errorCode = spfft_transform_destroy(transform)
             IF (errorCode.NE.SPFFT_SUCCESS) THEN
                CALL juDFT_error("Error in destroying spFFT transform! (2)", calledby="fft_interface")
             END IF

          END IF

#else
          CALL juDFT_error("Invalid state(1) in fft_interface", calledby="fft_interface")
#endif
       ELSE IF(fftRoutine.EQ.mklFFT_const) THEN
#ifdef CPP_FFT_MKL
          !using MKL library
          dfti_status = DftiCreateDescriptor(dfti_handle,dfti_double,dfti_complex,3,length)
          dfti_status = DftiCommitDescriptor(dfti_handle)
          if (forw) then
              dfti_status = DftiComputeForward(dfti_handle,dat)
          else
              dfti_status = DftiComputeBackward(dfti_handle,dat)
          end if
          dfti_status = DftiFreeDescriptor(dfti_handle)
#else
          CALL juDFT_error("Invalid state(2) in fft_interface", calledby="fft_interface")
#endif
       ELSE

          allocate(afft(size_dat),bfft(size_dat))
          afft = real(dat)
          bfft = aimag(dat)
          if (forw) then 
             isn = -1
          else
             isn = 1
          end if
          CALL cfft(afft,bfft,size_dat,length(1),length(1),isn)
          CALL cfft(afft,bfft,size_dat,length(2),length(1)*length(2),isn)
          CALL cfft(afft,bfft,size_dat,length(3),size_dat,isn)

          dat = cmplx(afft,bfft)

       END IF

    end subroutine fft_interface

END MODULE m_fft_interface
