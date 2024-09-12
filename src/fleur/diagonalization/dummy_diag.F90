! Copyright (c) 2018 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions
! of the MIT license as expressed in the LICENSE file in more detail.
!
! @authors: Miriam Hinzen, Gregor Michalicek
! Added MPI implementation, DW 2018
!--------------------------------------------------------------------------------
MODULE m_dummy_diag
  USE m_judft
  USE m_constants
  IMPLICIT NONE
  PRIVATE
  
  PUBLIC dummy_diag

CONTAINS

  SUBROUTINE dummy_diag(hmat,smat,ne,eig,zmat)
   !Dummy diver: does not solve actual eigenvalue problem but simply returns a set of orthogonal vectors. 
   !Could be useful for performance testing workloads in which we do not want to look at the diagonalization.
   ! A Cholesky decomp is still done to be able to do a back transform so that the resulting vector are orthonormal
   ! with respect to overlapp matrix.
    
    USE m_types_mat
    USE m_judft

    IMPLICIT NONE

    TYPE(t_mat),               INTENT(INOUT) :: hmat,smat
    INTEGER,                   INTENT(INOUT) :: ne
    CLASS(t_mat), ALLOCATABLE, INTENT(OUT)   :: zmat
    REAL,                      INTENT(OUT)   :: eig(:)

    INTEGER            :: nev,lwork,liwork,n
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


       ! --> solve a' * z' = eig * z' for eigenvalues eig between lb und ub
       zmat%data_r=0.0
       DO n=1,ne
         eig(ne)=-0.1+ne*1E-5
         zmat%data_r(ne,ne)=1.0
       enddo    
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

       
       ! --> solve a' * z' = eig * z' for eigenvalues eig between lb und ub
       zmat%data_c=0.0
       DO n=1,ne
         eig(ne)=-0.1+ne*1E-5
         zmat%data_c(ne,ne)=1.0
       enddo    
       
       ! --> recover the generalized eigenvectors z by solving z' = l^t * z
       CALL ztrtrs('U','N','N',hmat%matsize1,nev,smat%data_c,smat%matsize1,zMat%data_c,zmat%matsize1,info)
       IF (info.NE.0) THEN
          WRITE (oUnit,*) 'Error in ztrtrs: info =',info
          CALL juDFT_error("Diagonalization failed",calledby="chase_diag")
       ENDIF


    ENDIF
  END SUBROUTINE dummy_diag

  END MODULE m_dummy_diag
