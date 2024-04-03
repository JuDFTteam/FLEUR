! Copyright (c) 2018 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions
! of the MIT license as expressed in the LICENSE file in more detail.
!
! @authors: Miriam Hinzen, Gregor Michalicek
! Added MPI implementation, DW 2018
!--------------------------------------------------------------------------------
MODULE m_lapack_singlePrec_diag
  USE m_judft
  USE m_constants
  IMPLICIT NONE
  PRIVATE
  integer,parameter:: sp=selected_real_kind(6)
  
  PUBLIC lapack_singlePrec_diag

CONTAINS

  SUBROUTINE lapack_singlePrec_diag(hmat,smat,ne,eig,zmat)

    USE m_types_mat
    USE m_judft

    !Simple driver to solve Generalized Eigenvalue Problem using the ChASE library
    IMPLICIT NONE

    TYPE(t_mat),               INTENT(INOUT) :: hmat,smat
    INTEGER,                   INTENT(INOUT) :: ne
    CLASS(t_mat), ALLOCATABLE, INTENT(OUT)   :: zmat
    REAL,                      INTENT(OUT)   :: eig(:)

    INTEGER            :: nev,lwork,liwork
    INTEGER            :: info

    
    
    ALLOCATE(t_mat::zmat)
    CALL zmat%alloc(hmat%l_real,hmat%matsize1,ne)

    
    IF (hmat%l_real) THEN

       ! --> start with Cholesky factorization of b ( so that b = l * l^t)
       ! --> b is overwritten by l
       CALL dpotrf('U',smat%matsize1,smat%data_r,SIZE(smat%data_r,1),info)
       IF (info.NE.0) THEN
          WRITE (*,*) 'Error in dpotrf: info =',info
          CALL juDFT_error("Diagonalization failed",calledby="lapack_singlePrec_diag")
       ENDIF

       ! --> now reduce a * z = eig * b * z to the standard form a' * z' = eig * z'
       ! --> where a' = (l)^-1 * a * (l^t)^-1 and z' = l^t * z
       CALL dsygst(1,'U',smat%matsize1,hmat%data_r,SIZE(hmat%data_r,1),smat%data_r,SIZE(smat%data_r,1),info)
       IF (info.NE.0) THEN
          WRITE (oUnit,*) 'Error in dsygst: info =',info
          CALL juDFT_error("Diagonalization failed",calledby="lapack_singlePrec_diag")
       ENDIF

       ! --> solve a' * z' = eig * z' for eigenvalues eig between lb und ub
       BLOCK
        REAL(kind=sp),allocatable:: h(:,:),z(:,:),eigval(:),work(:)
        integer,allocatable      :: iwork(:),ifail(:)
        Allocate(h(size(hmat%data_r,1),size(hmat%data_r,2)))
        Allocate(eigval(size(hmat%data_r,1)),ifail(size(hmat%data_r,1)))
        Allocate(z(size(hmat%data_r,1),ne))
        h=hmat%data_r

        allocate(work(1),iwork(1))
        call ssyevx('V','I','U',size(h,1),h,size(h,1),0.0,0.0,1,ne,0.0,nev,eigval,z,size(z,1),work,-1,iwork,liwork,ifail,info)
        lwork=work(1)
        liwork=iwork(1)
        deallocate(work,iwork)
        allocate(work(lwork),iwork(liwork))

        call ssyevx('V','I','U',size(h,1),h,size(h,1),0.0,0.0,1,ne,0.0,nev,eigval,z,size(z,1),work,lwork,iwork,liwork,ifail,info)
        
        eig=eigval(:ne)
        zmat%data_r=z(:,:ne)
        deallocate(h,z,eigval,work,iwork)
       END BLOCK
       ! --> recover the generalized eigenvectors z by solving z' = l^t * z
       CALL dtrtrs('U','N','N',hmat%matsize1,nev,smat%data_r,smat%matsize1,zMat%data_r,zmat%matsize1,info)
       IF (info.NE.0) THEN
          WRITE (oUnit,*) 'Error in dtrtrs: info =',info
          CALL juDFT_error("Diagonalization failed",calledby="lapack_singlePrec_diag")
       ENDIF
       

    ELSE

       ! --> start with Cholesky factorization of b ( so that b = l * l^t)
       ! --> b is overwritten by l
       CALL zpotrf('U',smat%matsize1,smat%data_c,SIZE(smat%data_c,1),info)
       IF (info.NE.0) THEN
          WRITE (*,*) 'Error in zpotrf: info =',info
          CALL juDFT_error("Diagonalization failed",calledby="chase_diag")
       ENDIF

       ! --> now reduce a * z = eig * b * z to the standard form a' * z' = eig * z'
       ! --> where a' = (l)^-1 * a * (l^t)^-1 and z' = l^t * z
       CALL zhegst(1,'U',smat%matsize1,hmat%data_c,SIZE(hmat%data_c,1),smat%data_c,SIZE(smat%data_c,1),info)
       IF (info.NE.0) THEN
          WRITE (oUnit,*) 'Error in zhegst: info =',info
          CALL juDFT_error("Diagonalization failed",calledby="chase_diag")
       ENDIF

       ! --> solve a' * z' = eig * z' for eigenvalues eig between lb und ub
       BLOCK
        COMPLEX(kind=sp),allocatable:: h(:,:),z(:,:),work(:)
        REAL(kind=sp),allocatable:: eigval(:),rwork(:)
        integer,allocatable      :: iwork(:),ifail(:)
        Allocate(h(size(hmat%data_c,1),size(hmat%data_c,2)))
        Allocate(eigval(size(hmat%data_c,1)),ifail(size(hmat%data_c,1)))
        Allocate(z(size(hmat%data_c,1),ne),rwork(7*size(hmat%data_c,1)))
        h=hmat%data_c

        allocate(work(1),iwork(5*size(hmat%data_c,1)))
        call cheevx('V','I','U',size(h,1),h,size(h,1),0.0,0.0,1,ne,0.0,nev,eigval,z,size(z,1),work,-1,rwork,iwork,ifail,info)
        lwork=work(1)
        deallocate(work)
        allocate(work(lwork))

        call cheevx('V','I','U',size(h,1),h,size(h,1),0.0,0.0,1,ne,0.0,nev,eigval,z,size(z,1),work,lwork,rwork,iwork,ifail,info)
        eig=eigval(:ne)
        zmat%data_c=z(:,:ne)
        deallocate(h,z,eigval,work,rwork,iwork)
        END BLOCK   

       ! --> recover the generalized eigenvectors z by solving z' = l^t * z
       CALL ztrtrs('U','N','N',hmat%matsize1,nev,smat%data_c,smat%matsize1,zMat%data_c,zmat%matsize1,info)
       IF (info.NE.0) THEN
          WRITE (oUnit,*) 'Error in ztrtrs: info =',info
          CALL juDFT_error("Diagonalization failed",calledby="chase_diag")
       ENDIF


    ENDIF
    IF (info.NE.0) CALL judft_error("Diagonalization via ChASE failed", calledby = 'chase_diag')
  END SUBROUTINE lapack_singlePREC_diag

  END MODULE m_lapack_singlePrec_diag
