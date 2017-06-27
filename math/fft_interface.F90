MODULE m_fft_interface

 use m_cfft
 use m_juDFT
#ifdef CPP_FFT_MKL        
 use mkl_dfti
#endif

IMPLICIT NONE
CONTAINS

    subroutine fft_interface(dimen,length,dat,forw)
    ! provides interfaces to ftt subroutines

    integer, intent(in)    :: dimen         !dimension of fft transformation
    integer, intent(in)    :: length(dimen) !length of data in each direction
    complex, intent(inout) :: dat(:)        !data to be transformed, size(dat) should be sum(length)
    logical, intent(in)    :: forw          !.true. for the forward transformation, .false. for the backward one

#ifdef CPP_FFT_MKL        
    type(dfti_descriptor),pointer :: dfti_handle
    integer :: dfti_status
#else
    real,allocatable :: afft(:),bfft(:)
    integer :: isn
#endif
    integer :: size_dat,i
    
       size_dat = 1
       do i = 1,dimen
          size_dat = size_dat * length(i)
       enddo
       if (size(dat) .ne. size_dat) call juDFT_error('array bounds are inconsistent',calledby ='fft_interface')
       if (dimen .ne. 3 )  call juDFT_error('sorry, not implemented yet for this value of dimen',calledby ='fft_interface')
 
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
#endif

    end subroutine fft_interface

END MODULE m_fft_interface
