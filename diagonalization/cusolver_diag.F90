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
 PUBLIC cusolver_diag

CONTAINS
  SUBROUTINE cusolver_diag(hmat,smat,ne,eig,zmat)
    !Simple driver to solve Generalized Eigenvalue Problem using CuSolverDN
    IMPLICIT NONE
    CLASS(t_mat),INTENT(INOUT) :: hmat,smat
    INTEGER,INTENT(INOUT)      :: ne
    CLASS(t_mat),ALLOCATABLE,INTENT(OUT)    :: zmat
    REAL,INTENT(OUT)           :: eig(:)

#ifdef CPP_CUSOLVER
    INTEGER                 :: istat,ne_found,lwork_d,devinfo
    real,allocatable        :: work_d(:)
    type(cusolverDnHandle)  :: handle        

    istat = cusolverDnCreate(handle)
    if (istat /= CUSOLVER_STATUS_SUCCESS) call judft_error('handle creation failed')

    ALLOCATE(t_mat::zmat)
    CALL zmat%alloc(hmat%l_real,hmat%matsize1,ne)
    IF (hmat%l_real) THEN
      !$ACC DATA copyin(smat%data_r)COPY(hmat%data_r)COPYOUT(eig)
      !$ACC HOST_DATA USE_DEVICE(smat%data_r,hmat%data_r,eig)
      istat = cusolverDnDsygvdx_bufferSize(handle, CUSOLVER_EIG_TYPE_1, CUSOLVER_EIG_MODE_VECTOR, CUSOLVER_EIG_RANGE_I, CUBLAS_FILL_MODE_UPPER, hmat%matsize1, hmat%data_r, hmat%matsize1, &
      smat%data_r, smat%matsize1, 0.0, 0.0, 1, ne, ne_found, eig, lwork_d)
      if (istat /= CUSOLVER_STATUS_SUCCESS) call judft_error('cusolverDnZhegvdx_buffersize failed')
      allocate(work_d(lwork_d))
      !$ACC DATA create(work_d)
      !$ACC HOST_DATA USE_DEVICE(work_d)
      istat = cusolverDnDsygvdx(handle, CUSOLVER_EIG_TYPE_1, CUSOLVER_EIG_MODE_VECTOR, CUSOLVER_EIG_RANGE_I, CUBLAS_FILL_MODE_UPPER, hmat%matsize1, hmat%data_r, hmat%matsize1, &
      smat%data_r, smat%matsize1, 0.0, 0.0, 1, ne, ne_found, eig, work_d,lwork_d,devinfo)
      !$ACC END HOST_DATA
      !$ACC END HOST_DATA
      !$ACC END DATA
      !$ACC END DATA
      if (istat /= CUSOLVER_STATUS_SUCCESS) call judft_error('cusolverDnZhegvdx failed')
      ne=ne_found
      CALL zmat%alloc(hmat%l_real,hmat%matsize1,ne_found)
      zmat%data_c=hmat%data_c(:,:ne_found)
    
    ELSE
      !$ACC DATA copyin(smat%data_r)COPY(hmat%data_r)COPYOUT(eig)
      !$ACC HOST_DATA USE_DEVICE(smat%data_r,hmat%data_r,eig)
      istat = cusolverDnZhegvdx_bufferSize(handle, CUSOLVER_EIG_TYPE_1, CUSOLVER_EIG_MODE_VECTOR, CUSOLVER_EIG_RANGE_I, CUBLAS_FILL_MODE_UPPER, hmat%matsize1, hmat%data_c, hmat%matsize1, &
      smat%data_c, smat%matsize1, 0.0, 0.0, 1, ne, ne_found, eig, lwork_d)
      if (istat /= CUSOLVER_STATUS_SUCCESS) write(*,*) 'cusolverDnZhegvdx_buffersize failed'
      allocate(work_d(lwork_d))
      !$ACC DATA create(work_d)
      !$ACC HOST_DATA USE_DEVICE(work_d)
      istat = cusolverDnZhegvdx(handle, CUSOLVER_EIG_TYPE_1, CUSOLVER_EIG_MODE_VECTOR, CUSOLVER_EIG_RANGE_I, CUBLAS_FILL_MODE_UPPER, hmat%matsize1, hmat%data_c, hmat%matsize1, &
      smat%data_c, smat%matsize1, 0.0, 0.0, 1, ne, ne_found, eig, work_d,lwork_d,devinfo)
      !$ACC END HOST_DATA
      !$ACC END HOST_DATA
      !$ACC END DATA
      !$ACC END DATA
      if (istat /= CUSOLVER_STATUS_SUCCESS) call judft_error('cusolverDnZhegvdx failed')
      ne=ne_found
      CALL zmat%alloc(hmat%l_real,hmat%matsize1,ne_found)
      zmat%data_c=hmat%data_c(:,:ne_found)
    END IF
#endif
       
  END SUBROUTINE cusolver_diag

    
END MODULE m_cusolver_diag
