!--------------------------------------------------------------------------------
! Copyright (c) 2016 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------

MODULE m_aline
  USE m_juDFT
CONTAINS
  SUBROUTINE aline(eig_id, nk,atoms,sym,&
       cell,input, jsp,el,usdus,lapw,tlmplm, noco, oneD,eig,ne,zMat,hmat,smat)
    !************************************************************************
    !*                                                                      *
    !*     eigensystem-solver for moderatly-well converged potentials       *
    !*     a*z=e*b*z is transformed to h*z'=e*s*z' , whereby                *
    !*     h=C^T*a*C, s=C^T*b*C and z'=C^(-1)*z, when C is z of the last    *
    !*     iteration (lapw%nv*ne-array)                                          *
    !*     For ne<<lapw%nv the matrixsize is significantly reduced               *
    !*     aline uses ESSL-calls (use LAPACK's reduc3, tred3, bisect,       *
    !*     tinvit, trback and rebk3  if no ESSL available):                 *
    !*     SSPMV:  matrix-vector multiplication for symmetric matrices      *
    !*             in packed storage.                                       *
    !*     SSYGV:  eigensystem-solver for symmetric, real h and positive    *
    !*             definite, real, symmetric s using Cholesky-factorisation *
    !*             tridiagonalisation and a QL-algorithm.                   *
    !*     For all eigenvalues are needed, DSYGV should perform better      *
    !*     then seclr4 (hope so)                                            *
    !*                                                     Gustav           * *                                                                      *
    !************************************************************************

#include"cpp_double.h"
    USE m_abcof
    USE m_hssrwu
    USE m_eig66_io
    USE m_types
    IMPLICIT NONE
    
    TYPE(t_oneD),INTENT(IN)        :: oneD
    TYPE(t_input),INTENT(IN)       :: input
    TYPE(t_noco),INTENT(IN)        :: noco
    TYPE(t_sym),INTENT(IN)         :: sym
    TYPE(t_cell),INTENT(IN)        :: cell
    TYPE(t_atoms),INTENT(IN)       :: atoms
    TYPE(t_usdus),INTENT(IN)       :: usdus
    TYPE(t_tlmplm),INTENT(IN)      :: tlmplm
    TYPE(t_lapw),INTENT(INOUT)     :: lapw
    TYPE(t_mat),INTENT(INOUT)      :: zMat

    !     ..
    !     .. Scalar Arguments ..
    INTEGER, INTENT (IN) :: eig_id
    INTEGER, INTENT (IN) :: jsp,nk
    INTEGER, INTENT (OUT):: ne
    !     ..
    !     .. Array Arguments ..
    REAL,    INTENT (IN)  :: el(0:atoms%lmaxd,atoms%ntype,input%jspins)
    REAL,    INTENT (OUT) :: eig(input%neig)
    TYPE(t_mat),INTENT(IN):: hmat,smat

    !     ..
    !     .. Local Scalars ..
    INTEGER lhelp
    INTEGER i,info,j 
    !     ..
    !     .. Local Arrays ..
    COMPLEX, ALLOCATABLE :: acof(:,:,:),bcof(:,:,:),ccof(:,:,:,:)

    REAL,    ALLOCATABLE :: help_r(:),h_r(:,:),s_r(:,:) 
    REAL     CPP_BLAS_sdot
    EXTERNAL CPP_BLAS_sdot,CPP_BLAS_sspmv

    COMPLEX,   PARAMETER :: one_c=(1.0,0.0), zro_c=(0.0,0.0)
    COMPLEX, ALLOCATABLE :: help_c(:),h_c(:,:),s_c(:,:) 
    COMPLEX  CPP_BLAS_cdotc
    EXTERNAL CPP_BLAS_cdotc,CPP_BLAS_chpmv
    REAL,    ALLOCATABLE :: rwork(:)

    LOGICAL:: l_real
    l_real=zMat%l_real


    lhelp= MAX(lapw%nmat,(input%neig+2)*input%neig)
    CALL read_eig(eig_id,nk,jsp,neig=ne, eig=eig,zmat=zmat)
    IF (l_real) THEN
       ALLOCATE ( h_r(input%neig,input%neig),s_r(input%neig,input%neig) )
       h_r = 0.0 ; s_r=0.0
       ALLOCATE ( help_r(lhelp) )
    ELSE
       !     in outeig z is complex conjugated to make it usable for abcof. Here we 
       !                       first have to undo this  complex conjugation for the 
       ! multiplication with a and b matrices.

       zmat%data_c=conjg(zmat%data_c)
       ALLOCATE ( h_c(input%neig,input%neig),s_c(input%neig,input%neig) )
       h_c = 0.0 ; s_c=0.0
       ALLOCATE ( help_r(lhelp) )
    ENDIF
    !
    DO i = 1,ne
       IF (l_real) THEN
          help_r=MATMUL(hmat%data_r,zmat%data_r(:,i))
       ELSE
          help_c=MATMUL(hmat%data_c,zmat%data_c(:,i))
       ENDIF
       DO j = i,ne
          IF (l_real) THEN
             h_r(j,i)=dot_PRODUCT(zmat%data_r(:,j),help_r)
          ELSE
             h_c(j,i)=dot_PRODUCT(zmat%data_c(:,j),help_c)
          ENDIF
       END DO
    END DO

    DO i = 1,ne
       IF (l_real) THEN
          help_r=MATMUL(smat%data_r,zmat%data_r(:,i))
       ELSE
          help_c=MATMUL(smat%data_c,zmat%data_c(:,i))
       ENDIF
       DO j = i,ne
          IF (l_real) THEN
             s_r(j,i) = dot_product(zmat%data_r(:,j),help_r)
          ELSE
             s_c(j,i) =dot_PRODUCT(zmat%data_c(:,j),help_c)
          ENDIF
       END DO
    END DO

    ALLOCATE ( acof(input%neig,0:atoms%lmaxd*(atoms%lmaxd+2),atoms%nat),bcof(input%neig,0:atoms%lmaxd*(atoms%lmaxd+2),atoms%nat) )
    ALLOCATE ( ccof(-atoms%llod:atoms%llod,input%neig,atoms%nlod,atoms%nat) ) 

    !     conjugate again for use with abcof; finally use cdotc to revert again
    IF (.NOT.l_real) zMat%data_c = CONJG(zMat%data_c)
    if (noco%l_soc)  CALL juDFT_error("no SOC & reduced diagonalization",calledby="aline")

    CALL abcof(input,atoms,sym,cell,lapw,ne,&
         usdus,noco,1,oneD,acof,bcof,ccof,zMat)  ! ispin = 1&


    !
    CALL timestart("aline: hssr_wu")
    IF (l_real) THEN
       CALL hssr_wu(atoms,sym, jsp,el,ne,usdus,lapw,input,&
            tlmplm, acof,bcof,ccof, h_r,s_r)
    ELSE
       CALL hssr_wu(atoms,sym, jsp,el,ne,usdus,lapw,input,&
            tlmplm, acof,bcof,ccof, h_c=h_c,s_c=s_c)
    ENDIF

    DEALLOCATE ( ccof, bcof, acof )
    CALL timestop("aline: hssr_wu")
    CALL timestart("aline: seclr4")

    !
    IF (l_real) THEN
       !---> LAPACK call
       CALL CPP_LAPACK_ssygv(1,'V','L',ne,h_r,input%neig,s_r,input%neig,eig,help_r,lhelp,info)
    ELSE
       ALLOCATE ( rwork(MAX(1,3*ne-2)) )
       CALL CPP_LAPACK_chegv(1,'V','L',ne,h_c,input%neig,s_c,input%neig,eig,help_c,lhelp,rwork,info)
       DEALLOCATE ( rwork )
    ENDIF
    IF (info /= 0) THEN
       WRITE (6,FMT=8000) info
       IF (i < 0) THEN
          WRITE(6,'(a7,i3,a22)') 'element',info,' has an illegal value'
       ELSEIF (i > ne) THEN
          WRITE(6,'(a2,i3,a22)') 's:',info-ne,' not positive definite'
       ELSE
          WRITE(6,'(a8,i3,a15)') 'argument',info,' not  converged'
       ENDIF
       CALL juDFT_error("Diagonalisation failed",calledby ='aline')
    ENDIF
8000 FORMAT (' AFTER CPP_LAPACK_ssygv: info=',i4)
    CALL timestop("aline: seclr4")

    DO i = 1,lapw%nmat
       IF (l_real) THEN
          help_r(:ne)=zMat%data_r(i,:ne)
       ELSE
          help_c(:ne)=zMat%data_c(i,:ne)
       END IF
       DO j = 1,ne
          IF (l_real) THEN
             !--->       for LAPACK call
             zMat%data_r(i,j) = CPP_BLAS_sdot(ne,help_r,1,h_r(1,j),1)
          ELSE
             zMat%data_c(i,j) = CPP_BLAS_cdotc(ne,help_c,1,h_c(1,j),1)
          ENDIF
       END DO
    END DO

  END SUBROUTINE aline
END MODULE m_aline
