!--------------------------------------------------------------------------------
! Copyright (c) 2016 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------
MODULE m_lapack_diag
  USE m_types
  USE m_judft
IMPLICIT NONE
  CONTAINS
  SUBROUTINE lapack_diag(hmat,smat,ne,eig,zmat)
    !Simple driver to solve Generalized Eigenvalue Problem using LAPACK routine
    IMPLICIT NONE
    TYPE(t_mat),INTENT(INOUT)  :: hmat,smat
    INTEGER,INTENT(INOUT)      :: ne
    TYPE(t_mat),ALLOCATABLE,INTENT(OUT)    :: zmat
    REAL,INTENT(OUT)           :: eig(:)

    INTEGER            :: lwork,info,m
    INTEGER,ALLOCATABLE:: ifail(:),iwork(:)
    COMPLEX,ALLOCATABLE:: work(:)
    REAL,ALLOCATABLE   :: rwork(:)
    REAL               :: dumrwork(1),abstol
    COMPLEX            :: dumwork(1)
    REAL,external      :: dlamch

    ALLOCATE(t_mat::zmat)
    CALL zmat%alloc(hmat%l_real,hmat%matsize1,ne)
    abstol=2*dlamch('S')
    IF (hmat%l_real) THEN
       ALLOCATE(iwork(5*hmat%matsize1),ifail(hmat%matsize1))
       CALL dsygvx(1,'V','I','U', hmat%matsize1,hmat%data_r,SIZE(hmat%data_r,1),smat%data_r,SIZE(smat%data_r,1),&
            0.0,0.0,1,ne,abstol,m,eig,zmat%data_r,SIZE(zmat%data_r,1),dumrwork,-1, iwork, ifail, info)
       lwork=dumrwork(1)
       ALLOCATE(rwork(lwork))
       CALL dsygvx(1,'V','I','U', hmat%matsize1,hmat%data_r,SIZE(hmat%data_r,1),smat%data_r,SIZE(smat%data_r,1),&
            0.0,0.0,1,ne,abstol,m,eig,zmat%data_r,SIZE(zmat%data_r,1),rwork, lwork, iwork, ifail, info)
    ELSE
       ALLOCATE(rwork(7*hmat%matsize1),iwork(5*hmat%matsize1),ifail(hmat%matsize1))
       !Do a workspace query
       CALL zhegvx(1,'V','I','U',hmat%matsize1,hmat%data_c,SIZE(hmat%data_c,1),smat%data_c,SIZE(smat%data_c,1),&
            0.0,0.0,1,ne,abstol,m,eig,zmat%data_c,SIZE(zmat%data_c,1),dumwork,-1,rwork,iwork,ifail,info)
       lwork=dumwork(1)
       ALLOCATE(work(lwork))
       !Perform diagonalization
       CALL zhegvx(1,'V','I','U',hmat%matsize1,hmat%data_c,SIZE(hmat%data_c,1),smat%data_c,SIZE(smat%data_c,1),&
            0.0,0.0,1,ne,abstol,m,eig,zmat%data_c,SIZE(zmat%data_c,1),work,lwork,rwork,iwork,ifail,info)
    ENDIF
    IF (info.NE.0) CALL judft_error("Diagonalization via LAPACK failed")
    IF (m.NE.ne) CALL judft_error("Diagonalization via LAPACK failed failed without explicit errorcode.")
  END SUBROUTINE lapack_diag
END MODULE m_lapack_diag
