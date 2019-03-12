MODULE m_exact_diag
   !************************************************************************
   !*This modules contains a driver for the eigenvalue problem Ax=lambdax
   !*We don't use the standard routines to avoid setting up a big unity matrix
   !************************************************************************
   USE m_juDFT
   USE m_types
   IMPLICIT NONE
   !
   !TODO: add other diagonalization routines
   !
   CONTAINS

   SUBROUTINE exact_diag(mpi,hmat,ne,eig,ev)

      USE m_types 

      TYPE(t_mpi),                  INTENT(IN)     :: mpi 
      CLASS(t_mat),                 INTENT(INOUT)  :: hmat
      CLASS(t_mat),  ALLOCATABLE,   INTENT(OUT)    :: ev 
      INTEGER,                      INTENT(IN)     :: ne 
      REAL,                         INTENT(OUT)    :: eig(:)
      




   END SUBROUTINE exact_diag

   SUBROUTINE lapack_ex_diag(hmat,ne,eig,zmat)
      !Simple driver for solving eigenvalue problem with Overlap matrix being Unity

      IMPLICIT NONE
      TYPE(t_mat),               INTENT(INOUT)  :: hmat
      INTEGER,                   INTENT(INOUT)  :: ne
      CLASS(t_mat),ALLOCATABLE,  INTENT(OUT)    :: zmat
      REAL,                      INTENT(OUT)    :: eig(:)

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
      IF(hmat%l_real) THEN
         ALLOCATE(iwork(5*hmat%matsize1),ifail(hmat%matsize1))
         CALL dsyevx("V","I","U",hmat%matsize1,hmat%data_r,SIZE(hmat%data_r,1),0.0,0.0,1,ne,&
                     abstol,m,eig,zmat%data_r,SIZE(zmat%data_r,1),dumrwork,-1, iwork, ifail, info) 
         lwork=dumrwork(1)
         ALLOCATE(rwork(lwork))
         IF (info.NE.0) CALL judft_error("Diagonalization via LAPACK failed (Workspace)",no=info)
         CALL dsyevx("V","I","U",hmat%matsize1,hmat%data_r,SIZE(hmat%data_r,1),0.0,0.0,1,ne,&
                     abstol,m,eig,zmat%data_r,SIZE(zmat%data_r,1),rwork,lwork, iwork, ifail, info) 
      ELSE
         ALLOCATE(rwork(7*hmat%matsize1),iwork(5*hmat%matsize1),ifail(hmat%matsize1))
         !Do a workspace query
         CALL zheevx('V','I','U',hmat%matsize1,hmat%data_c,SIZE(hmat%data_c,1),&
            0.0,0.0,1,ne,abstol,m,eig,zmat%data_c,SIZE(zmat%data_c,1),dumwork,-1,rwork,iwork,ifail,info)
         lwork=dumwork(1)
         ALLOCATE(work(lwork))
         IF (info.NE.0) CALL judft_error("Diagonalization via LAPACK failed (Workspace)",no=info)
         !Perform diagonalization
         CALL zheevx('V','I','U',hmat%matsize1,hmat%data_c,SIZE(hmat%data_c,1),&
            0.0,0.0,1,ne,abstol,m,eig,zmat%data_c,SIZE(zmat%data_c,1),work,lwork,rwork,iwork,ifail,info)
      ENDIF
      IF (info.NE.0) CALL judft_error("Diagonalization via LAPACK failed(zhegvx/dsygvx)",no=info)
      IF (m.NE.ne) CALL judft_error("Diagonalization via LAPACK failed failed without explicit errorcode.")

   END SUBROUTINE lapack_ex_diag
END MODULE m_exact_diag