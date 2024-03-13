!--------------------------------------------------------------------------------
! Copyright (c) 2016 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------
MODULE m_cusolver_diag
  USE m_types_mat
  USE m_types_mpimat
  USE m_judft
#ifdef CPP_CUSOLVER
  use cusolverDn  
#endif  
  IMPLICIT NONE
  PRIVATE
#ifdef CPP_CUSOLVER
  type(cusolverDnHandle)  :: handle        
#endif  
  PUBLIC cusolver_diag

CONTAINS
  SUBROUTINE cusolver_diag(hmat,smat,ne,eig,zmat)
    !!Simple driver to solve Generalized Eigenvalue Problem using CuSolverDN
    IMPLICIT NONE
    CLASS(t_mat),INTENT(INOUT) :: hmat,smat 
    INTEGER,INTENT(INOUT)      :: ne
    CLASS(t_mat),ALLOCATABLE,INTENT(OUT)    :: zmat
    REAL,INTENT(OUT)           :: eig(:)

#ifdef CPP_CUSOLVER
    INTEGER                 :: istat,ne_found,lwork_d,devinfo(1)
    real,allocatable        :: work_d(:),eig_tmp(:)
    complex,allocatable     :: work_c(:)

    logical :: firstcall=.true.
    if (firstcall) istat=cusolverDnCreate(handle)
    if (istat /= CUSOLVER_STATUS_SUCCESS) call judft_error('handle creation failed')

    ALLOCATE(t_mat::zmat)
    ALLOCATE(eig_tmp(hmat%matsize1))
    CALL zmat%alloc(hmat%l_real,hmat%matsize1,ne)
    !!$acc Data copyin(hmat,smat)
    IF (hmat%l_real) THEN
      associate(h=>hmat%data_r,s=>smat%data_r)
        !$ACC DATA copyin(s)COPY(h)COPYOUT(eig_tmp)
        !$ACC HOST_DATA USE_DEVICE(s,h,eig_tmp)
        istat = cusolverDnDsygvdx_bufferSize(handle, CUSOLVER_EIG_TYPE_1, CUSOLVER_EIG_MODE_VECTOR, CUSOLVER_EIG_RANGE_I, CUBLAS_FILL_MODE_UPPER, hmat%matsize1, h, hmat%matsize1, &
            s, smat%matsize1, 0.0, 0.0, 1, ne, ne_found, eig_tmp, lwork_d)
        !$acc end host_data
        if (istat /= CUSOLVER_STATUS_SUCCESS) call judft_error('cusolverDnZhegvdx_buffersize failed')
        allocate(work_d(lwork_d))
        !$ACC DATA create(work_d,devinfo)
        !$ACC HOST_DATA USE_DEVICE(s,h,eig_tmp,work_d,devinfo)
        istat = cusolverDnDsygvdx(handle, CUSOLVER_EIG_TYPE_1, CUSOLVER_EIG_MODE_VECTOR, CUSOLVER_EIG_RANGE_I, CUBLAS_FILL_MODE_UPPER, hmat%matsize1, h, hmat%matsize1, &
        s, smat%matsize1, 0.0, 0.0, 1, ne, ne_found, eig_tmp, work_d,lwork_d,devinfo(1))
        !$ACC END HOST_DATA
        !$ACC END DATA
        !$ACC END DATA
        if (istat /= CUSOLVER_STATUS_SUCCESS) call judft_error('cusolverDnZhegvdx failed')
        ne=ne_found
        CALL zmat%alloc(hmat%l_real,hmat%matsize1,ne_found)
        zmat%data_r=h(:,:ne_found)
        eig=eig_tmp(:ne)
      end associate
    ELSE
      associate(h=>hmat%data_c,s=>smat%data_c)
        !$ACC DATA copyin(s) COPY(h) COPYOUT(eig_tmp)
        !$ACC HOST_DATA USE_DEVICE(s,h,eig_tmp)
        istat = cusolverDnZhegvdx_bufferSize(handle, CUSOLVER_EIG_TYPE_1, CUSOLVER_EIG_MODE_VECTOR, CUSOLVER_EIG_RANGE_I, CUBLAS_FILL_MODE_UPPER, hmat%matsize1, h, hmat%matsize1, &
          s, smat%matsize1, 0.0, 0.0, 1, ne, ne_found, eig_tmp, lwork_d)
        !$acc end host_data
        if (istat /= CUSOLVER_STATUS_SUCCESS) write(*,*) 'cusolverDnZhegvdx_buffersize failed'
        allocate(work_c(lwork_d))
        !$ACC DATA create(work_c,devinfo)
        !$ACC HOST_DATA USE_DEVICE(s,h,eig_tmp,work_c,devinfo)
        istat = cusolverDnZhegvdx(handle, CUSOLVER_EIG_TYPE_1, CUSOLVER_EIG_MODE_VECTOR, CUSOLVER_EIG_RANGE_I, CUBLAS_FILL_MODE_UPPER, hmat%matsize1, h, hmat%matsize1, &
          s, smat%matsize1, 0.0, 0.0, 1, ne, ne_found, eig_tmp, work_c,lwork_d,devinfo(1))
        !$ACC END HOST_DATA
        !$acc update self(devinfo)
        if (istat /= CUSOLVER_STATUS_SUCCESS) THEN 
             write(*,*) devinfo
             call judft_error('cusolverDnZhegvdx failed')
        endif     
        !$ACC END DATA
        !$ACC END DATA
        ne  =ne_found
        CALL zmat%alloc(hmat%l_real,hmat%matsize1,ne_found)
        zmat%data_c=h(:,:ne_found)
        eig=eig_tmp(:ne)
    
      end associate  
    END IF
#endif
       
  END SUBROUTINE cusolver_diag

    
END MODULE m_cusolver_diag
