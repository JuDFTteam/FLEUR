!--------------------------------------------------------------------------------
! Copyright (c) 2016 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------

MODULE m_aline
  USE m_juDFT
CONTAINS
  SUBROUTINE aline(eig_id, nk,atoms,DIMENSION,sym,&
       cell,input, jsp,el,usdus, a,b,lapw,tlmplm, noco, oneD, bkpt,z,eig,ne)
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
    TYPE(t_dimension),INTENT(IN)   :: DIMENSION
    TYPE(t_oneD),INTENT(IN)        :: oneD
    TYPE(t_input),INTENT(IN)       :: input
    TYPE(t_noco),INTENT(IN)        :: noco
    TYPE(t_sym),INTENT(IN)         :: sym
    TYPE(t_cell),INTENT(IN)        :: cell
    TYPE(t_atoms),INTENT(IN)       :: atoms
    TYPE(t_usdus),INTENT(IN)       :: usdus
    TYPE(t_tlmplm),INTENT(IN)      :: tlmplm
    TYPE(t_lapw),INTENT(INOUT)     :: lapw
    !     ..
    !     .. Scalar Arguments ..
    INTEGER, INTENT (IN) :: eig_id
    INTEGER, INTENT (IN) :: jsp,nk
    INTEGER, INTENT (OUT):: ne
    !     ..
    !     .. Array Arguments ..
    REAL,    INTENT (IN)  :: el(0:atoms%lmaxd,atoms%ntypd,DIMENSION%jspd)
    REAL,    INTENT (OUT) :: eig(DIMENSION%neigd),bkpt(3)
#ifdef CPP_INVERSION
    REAL,    INTENT (IN)  :: a(:),b(:)!(matsize)
    REAL,    INTENT (OUT) :: z(DIMENSION%nbasfcn,DIMENSION%neigd)
#else
    COMPLEX, INTENT (IN)  :: a(:),b(:)!(matsize)
    COMPLEX, INTENT (OUT) :: z(DIMENSION%nbasfcn,DIMENSION%neigd)
#endif
    !     ..
    !     .. Local Scalars ..
    INTEGER lhelp
    INTEGER i,info,j 
    !     ..
    !     .. Local Arrays ..
    INTEGER kveclo(atoms%nlotot)
    COMPLEX, ALLOCATABLE :: acof(:,:,:),bcof(:,:,:),ccof(:,:,:,:)
#ifdef CPP_INVERSION
    REAL,      PARAMETER :: one=1.0, zro=0.0
    REAL,    ALLOCATABLE :: help(:),h(:,:),s(:,:) 
    REAL     CPP_BLAS_sdot
    EXTERNAL CPP_BLAS_sdot,CPP_BLAS_sspmv
#else
    COMPLEX,   PARAMETER :: one=(1.0,0.0), zro=(0.0,0.0)
    COMPLEX, ALLOCATABLE :: help(:),h(:,:),s(:,:) 
    COMPLEX  CPP_BLAS_cdotc
    EXTERNAL CPP_BLAS_cdotc,CPP_BLAS_chpmv
    REAL,    ALLOCATABLE :: rwork(:)
#endif
#if ( defined(CPP_INVERSION) && !defined(CPP_SOC) )
     REAL     zhlp(DIMENSION%nbasfcn,DIMENSION%neigd)
#else
    COMPLEX  zhlp(DIMENSION%nbasfcn,DIMENSION%neigd)
#endif

    !     ..
    !     .. External Subroutines ..
    EXTERNAL CPP_LAPACK_ssygv
    !     ..
    !

    CALL read_eig(eig_id,nk,jsp,bk=bkpt,neig=ne,nv=lapw%nv(jsp),nmat=lapw%nmat,&
         &      eig=eig,kveclo=kveclo,z=z)
#ifndef CPP_INVERSION
    !     in outeig z is complex conjugated to make it usable for abcof. Here we 
    !                       first have to undo this  complex conjugation for the 
    z = CONJG(z)    ! multiplication with a and b matrices.
#endif

    !
    ALLOCATE ( h(DIMENSION%neigd,DIMENSION%neigd),s(DIMENSION%neigd,DIMENSION%neigd) )
    h = zro ; s=zro
    IF (lapw%nmat.GT.(DIMENSION%neigd+2)*DIMENSION%neigd) THEN
       ALLOCATE ( help(lapw%nmat) )
    ELSE
       ALLOCATE ( help((DIMENSION%neigd+2)*DIMENSION%neigd) )
    ENDIF
    lhelp= (DIMENSION%neigd+2)*DIMENSION%neigd

    !
    DO i = 1,ne
#ifdef CPP_INVERSION
       CALL CPP_BLAS_sspmv('U',lapw%nmat,one,a,z(1,i),1,zro,help,1)
#else
       CALL CPP_BLAS_chpmv('U',lapw%nmat,one,a,z(1,i),1,zro,help,1)
#endif
       DO j = i,ne
#ifdef CPP_INVERSION
          h(j,i) = CPP_BLAS_sdot(lapw%nmat,z(1,j),1,help,1)
#else
          h(j,i) = CPP_BLAS_cdotc(lapw%nmat,z(1,j),1,help,1)
#endif
       END DO
    END DO

    DO i = 1,ne
#ifdef CPP_INVERSION
       CALL CPP_BLAS_sspmv('U',lapw%nmat,one,b,z(1,i),1,zro,help,1)
#else
       CALL CPP_BLAS_chpmv('U',lapw%nmat,one,b,z(1,i),1,zro,help,1)
#endif
       DO j = i,ne
#ifdef CPP_INVERSION
          s(j,i) = CPP_BLAS_sdot(lapw%nmat,z(1,j),1,help,1)
#else
          s(j,i) = CPP_BLAS_cdotc(lapw%nmat,z(1,j),1,help,1)
#endif
       END DO
    END DO

    ALLOCATE ( acof(DIMENSION%neigd,0:DIMENSION%lmd,atoms%natd),bcof(DIMENSION%neigd,0:DIMENSION%lmd,atoms%natd) )
    ALLOCATE ( ccof(-atoms%llod:atoms%llod,DIMENSION%neigd,atoms%nlod,atoms%natd) ) 
#ifndef CPP_INVERSION
    !     conjugate again for use with abcof; finally use cdotc to revert again
    z = CONJG(z)
#endif
#ifdef CPP_SOC
    CALL juDFT_error("no SOC & reduced diagonalization",calledby="aline")
#else
    CALL abcof(input,atoms,dimension%neigd,sym,cell, bkpt,lapw,ne,z,&
         usdus,noco,1,kveclo,oneD,acof,bcof,ccof)  ! ispin = 1&
         
#endif
    !
    CALL timestart("aline: hssr_wu")
    CALL hssr_wu(atoms,DIMENSION,sym, jsp,el,ne,usdus,lapw,&
         tlmplm, acof,bcof,ccof, h,s)

    DEALLOCATE ( ccof, bcof, acof )
    CALL timestop("aline: hssr_wu")
    CALL timestart("aline: seclr4")

    !
#ifdef CPP_INVERSION
    !---> LAPACK call
    CALL CPP_LAPACK_ssygv(1,'V','L',ne,h,DIMENSION%neigd,s,DIMENSION%neigd,eig,help,lhelp,info)
#else
    ALLOCATE ( rwork(MAX(1,3*ne-2)) )
    CALL CPP_LAPACK_chegv(1,'V','L',ne,h,DIMENSION%neigd,s,DIMENSION%neigd,eig,help,lhelp,rwork,info)
    DEALLOCATE ( rwork )
#endif
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
       DO j = 1,ne
          help(j) = z(i,j)
       END DO
       DO j = 1,ne
#ifdef CPP_INVERSION
          !--->       for LAPACK call
          z(i,j) = CPP_BLAS_sdot(ne,help,1,h(1,j),1)
#else
          z(i,j) = CPP_BLAS_cdotc(ne,help,1,h(1,j),1)
#endif
       END DO
    END DO
    DEALLOCATE ( help,h,s )

  END SUBROUTINE aline
END MODULE m_aline
