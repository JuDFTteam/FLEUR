!--------------------------------------------------------------------------------
! Copyright (c) 2016 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------

MODULE m_magma
  use m_juDFT
  INTEGER,SAVE :: Magma_NumGPU=1
  !**********************************************************
  !     Solve the generalized eigenvalue problem
  !     using the MAGMA library for multiple GPUs
  !**********************************************************
CONTAINS
  SUBROUTINE magma_diag(hmat,smat,ne,eig,zmat)
#ifdef CPP_MAGMA
    use magma
#endif
    use m_types
    IMPLICIT NONE

    ! ... Arguments ...
    TYPE(t_mat),INTENT(INOUT)  :: hmat,smat
    INTEGER,INTENT(INOUT)      :: ne
    CLASS(t_mat),ALLOCATABLE,INTENT(OUT)    :: zmat
    REAL,INTENT(OUT)           :: eig(:)

#ifdef CPP_MAGMA

    ! ... Local Variables ..
    INTEGER :: lwork,liwork,lrwork,error,mout(1),i
    REAL    :: eigTemp(hmat%matsize1)
    LOGICAL :: initialized=.false.


    REAL,    ALLOCATABLE :: rwork(:)
    INTEGER, ALLOCATABLE :: iwork(:)
    COMPLEX, ALLOCATABLE :: work(:)

#ifdef _OPENACC
    Magma_NumGPU=acc_get_num_devices(acc_device_nvidia)
#endif
    IF (.NOT.initialized) THEN
       initialized=.TRUE.
       CALL magmaf_init()
    ENDIF

    IF (hmat%l_real) THEN
       ALLOCATE(rwork(1),iwork(1))
       CALL magmaf_dsygvdx_m(Magma_numGPU,1,'v','i','U',hmat%matsize1,hmat%data_r,SIZE(hmat%data_r,1),smat%data_r,&
                           SIZE(smat%data_r,1),0.0,0.0,1,ne,mout,eigTemp,rwork,-1,iwork,-1,error)
       IF (error/=0) THEN
          WRITE(*,*) 'magmaf_dsygvdx error code: ', error
          CALL juDFT_error("Failed to query workspaces (1)",calledby="magma.F90")
       END IF
       lrwork=rwork(1)
       liwork=iwork(1)
       DEALLOCATE(rwork,iwork)
       ALLOCATE(rwork(lrwork),iwork(liwork))
       CALL magmaf_dsygvdx_m(Magma_numGPU,1,'v','i','U',hmat%matsize1,hmat%data_r,SIZE(hmat%data_r,1),smat%data_r,&
                           SIZE(smat%data_r,1),0.0,0.0,1,ne,mout,eigTemp,rwork,lrwork,iwork,liwork,error)
       IF (error/=0) THEN
          WRITE(*,*) 'magmaf_dsygvdx error code: ', error
          CALL juDFT_error("Magma failed to diagonalize Hamiltonian (1)",calledby="magma.F90")
       END IF
    ELSE
       !Query the workspace size
       ALLOCATE(work(1),rwork(1),iwork(1))
       !CALL magmaf_zhegvdx_2stage_m(NGPU_CONST,&
       CALL magmaf_zhegvdx_m(Magma_numGPU,1,'v','i','U',hmat%matsize1,hmat%data_c,SIZE(hmat%data_c,1),smat%data_c,&
                           SIZE(smat%data_c,1),0.0,0.0,1,ne,mout,eigTemp,work,-1,rwork,-1,iwork,-1,error)
       IF (error/=0) THEN
          WRITE(*,*) 'magmaf_zhegvdx error code: ', error
          CALL juDFT_error("Failed to query workspaces (2)",calledby="magma.F90")
       END IF
       lwork=work(1)
       lrwork=rwork(1)
       liwork=iwork(1)
       DEALLOCATE(work,rwork,iwork)
       ALLOCATE(work(lwork),rwork(lrwork),iwork(liwork))
       !Now the diagonalization
       !CALL magmaf_zhegvdx_2stage_m(NGPU_CONST,&
       CALL magmaf_zhegvdx_m(Magma_numGPU,1,'v','i','U',hmat%matsize1,hmat%data_c,SIZE(hmat%data_c,1),smat%data_c,&
                           SIZE(smat%data_c,1),0.0,0.0,1,ne,mout,eigTemp,work,lwork,rwork,lrwork,iwork,liwork,error)
       IF (error/=0) THEN
          WRITE(*,*) 'magmaf_zhegvdx error code: ', error
          CALL juDFT_error("Magma failed to diagonalize Hamiltonian (2)",calledby="magma.F90")
       END IF
    ENDIF

    ALLOCATE(t_mat::zmat)
    CALL zmat%alloc(hmat%l_real,hmat%matsize1,ne)
    DO i = 1, ne
       eig(i) = eigTemp(i)
       IF (hmat%l_real) THEN
          zmat%data_r(:,i)=hmat%data_r(:,i)
       ELSE
          zmat%data_c(:,i)=hmat%data_c(:,i)
       ENDIF
    END DO
#endif
  END SUBROUTINE magma_diag
END MODULE m_magma
