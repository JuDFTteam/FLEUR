PROGRAM test 
    use cusolverdn
    INTEGER,PARAMETER ::N=100
    INTEGER                 :: istat,ne,ne_found,lwork_d,devinfo
    real                    :: eig(n),h(N,N)
    type(cusolverDnHandle)  :: handle 

    ne=10
    istat = cusolverDnCreate(handle)
    if (istat /= CUSOLVER_STATUS_SUCCESS) write(*,*) 'handle creation failed'
    !$ACC DATA copyin(H,eig)
    !$ACC HOST_DATA USE_DEVICE(H,eig)
    istat = cusolverDnDsygvdx_bufferSize(handle, CUSOLVER_EIG_TYPE_1, CUSOLVER_EIG_MODE_VECTOR, CUSOLVER_EIG_RANGE_I, CUBLAS_FILL_MODE_UPPER, N, H, N, &
        H, N, 0.0, 0.0, 1, ne, ne_found, eig, lwork_d)
    !$ACC END HOST_DATA
    !$ACC END DATA
END