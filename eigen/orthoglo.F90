!--------------------------------------------------------------------------------
! Copyright (c) 2016 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------

MODULE m_orthoglo
  USE m_juDFT
  !*********************************************************************
  ! Each G-vector corresponds to a vector of C-coeff. These vectors must
  ! be linearly independent. This is checked by this soubroutine for an
  ! atom that doesn't have an inversion partner.
  ! Philipp Kurz 99/04
  !*********************************************************************
CONTAINS
  SUBROUTINE orthoglo(l_real,atoms, nkvec,lo,l,linindq,l_lo2, cwork, linind)
    !
    !*************** ABBREVIATIONS ***************************************
    ! cwork   : contains the vectors of C-coeff.
    ! l_lo2   : changes this routine to old 'orthoglo2': same as orthoglo, 
    !           but for a pair of atoms that can be mapped onto eachother 
    !           by inversion.
    ! CF Replaced (unstable) Gram-Schmidt by diagonalization.
    !*********************************************************************
    !
#include"cpp_double.h"
    !
    USE m_types_setup
    IMPLICIT NONE
    TYPE(t_atoms),INTENT(IN)   :: atoms
    !     ..
    !     .. Scalar Arguments ..
    INTEGER, INTENT (IN) :: l,lo,nkvec
    REAL,    INTENT (IN) :: linindq
    LOGICAL, INTENT (IN) :: l_lo2,l_real
    LOGICAL, INTENT (OUT) :: linind
    !     ..
    !     .. Array Arguments ..
    COMPLEX,INTENT (INOUT):: cwork(-2*atoms%llod:2*atoms%llod+1,2*(2*atoms%llod+1) ,atoms%nlod)
    !     ..
    !     .. Local Scalars ..
    INTEGER dim,low,i,j
    !     ..
    !     .. Local Arrays ..
    REAL eig(nkvec),rwork(3*nkvec)
    REAL olap_r(nkvec,nkvec)
    EXTERNAL CPP_LAPACK_ssyev
    COMPLEX olap_c(nkvec,nkvec),work(2*nkvec)
    EXTERNAL CPP_LAPACK_cheev

    IF (l_lo2) THEN
       dim = 2* (2*l+1)
       low = -2*l
    ELSE
       dim = 2*l+1
       low = -l
    ENDIF

    DO i = 1,nkvec
       DO j = 1,nkvec
          IF (l_real) THEN
             olap_r(i,j) = DOT_PRODUCT(cwork(low:low+dim-1,i,lo), cwork(low:low+dim-1,j,lo))
          ELSE
             olap_c(i,j) = DOT_PRODUCT(cwork(low:low+dim-1,i,lo), cwork(low:low+dim-1,j,lo))
          ENDIF
       ENDDO
    ENDDO
    IF (l_real) THEN
       CALL CPP_LAPACK_ssyev('N','U',nkvec,olap_r,nkvec,eig, rwork,3*nkvec,i)
       IF(i/=0)  CALL juDFT_error("(S,D)SYEV failed.","orthoglo")
    ELSE
       CALL CPP_LAPACK_cheev('N','U',nkvec,olap_c,nkvec,eig, work,2*nkvec,rwork,i)
       IF(i/=0)  CALL juDFT_error("(C,Z)HEEV failed.","orthoglo")
    ENDIF
    IF(eig(1).LT.linindq) THEN
       linind = .FALSE.
    ELSE
       linind = .TRUE.
    ENDIF
    RETURN

  END SUBROUTINE orthoglo
END MODULE m_orthoglo
