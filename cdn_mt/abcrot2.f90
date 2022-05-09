!--------------------------------------------------------------------------------
! Copyright (c) 2016 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------

MODULE m_abcrot2
  PRIVATE
  PUBLIC :: abcrot2
CONTAINS
  SUBROUTINE abcrot2(itype,na,atoms,banddos,eigVecCoeffs,jsp,acof,bcof,ccof)
    USE m_dwigner
    USE m_types
    IMPLICIT NONE
    INTEGER,INTENT(IN)                 :: itype,na
    TYPE(t_atoms),INTENT(IN)           :: atoms
    TYPE(t_banddos),INTENT(IN)         :: banddos
    TYPE(t_eigVecCoeffs),INTENT(IN)   :: eigVecCoeffs
    COMPLEX, ALLOCATABLE,INTENT(INOUT) :: acof(:,:)
    COMPLEX, ALLOCATABLE,INTENT(INOUT) :: bcof(:,:)
    COMPLEX, ALLOCATABLE,INTENT(INOUT) :: ccof(:,:,:)
    !     ..
    !     .. Scalar Arguments ..
    INTEGER, INTENT (IN) :: jsp
    !     ..
    !     .. Local Scalars ..
    INTEGER ilo,i,l
    REAL amx(3,3,1),imx(3,3)
    COMPLEX  d_wgn(-atoms%lmaxd:atoms%lmaxd,-atoms%lmaxd:atoms%lmaxd,1:atoms%lmaxd,1)


    CALL euler(banddos%alpha(na),banddos%beta(na),banddos%gamma(na), amx)

    imx(:,:) = 0. ; imx(1,1) = 1. ; imx(2,2) = 1. ; imx(3,3) = 1.

    CALL d_wigner(1,amx,imx,atoms%lmaxd, d_wgn)

    DO l = 1, atoms%lmax(itype)
      DO i = 1, size(acof,1)
        acof(i,l**2:l*(l+2)) = MATMUL(CONJG(d_wgn(-l:l,-l:l,l,1)),&
        eigVecCoeffs%acof2(i,l**2:l*(l+2),0,na,jsp))
        bcof(i,l**2:l*(l+2)) = MATMUL(CONJG(d_wgn(-l:l,-l:l,l,1)),&
        eigVecCoeffs%acof2(i,l**2:l*(l+2),1,na,jsp))
      ENDDO
    ENDDO
    DO ilo = 1, atoms%nlo(itype)
      l = atoms%llo(ilo,itype)
      IF (l.GT.0) THEN
        DO i = 1 ,size(acof,1)
          ccof(-l:l,i,ilo) = MATMUL(CONJG(d_wgn(-l:l,-l:l,l,1)),&
          eigVecCoeffs%ccof(-l:l,i,ilo,na,jsp))
        ENDDO
      ENDIF
    ENDDO
  END SUBROUTINE abcrot2

  !********************************************************************
  !********************************************************************
  SUBROUTINE euler(alpha,beta,gamma,amx)
    IMPLICIT NONE

    REAL,    INTENT (IN)  :: alpha,beta,gamma
    REAL,    INTENT (OUT) :: amx(3,3,1)

    REAL  alph,bet,gamm
    REAL bmx(3,3),cmx(3,3),dmx(3,3),hmx(3,3)
    INTEGER nwf,i,j,ii

    !..define the D,C,B-matrices
    amx(:,:,:)=0.

    alph = alpha ; bet = beta ; gamm = gamma

    dmx(1,1) = COS(alph) ; dmx(1,2) = SIN(alph) ; dmx(1,3) = 0.
    dmx(2,1) =-SIN(alph) ; dmx(2,2) = COS(alph) ; dmx(2,3) = 0.
    dmx(3,1) = 0.        ; dmx(3,2) = 0.        ; dmx(3,3) = 1.

    cmx(1,1) = 1.  ; cmx(1,2) = 0.        ; cmx(1,3) = 0.
    cmx(2,1) = 0.  ; cmx(2,2) = COS(bet)  ; cmx(2,3) = SIN(bet)
    cmx(3,1) = 0.  ; cmx(3,2) =-SIN(bet)  ; cmx(3,3) = COS(bet)

    bmx(1,1) = COS(gamm) ; bmx(1,2) = SIN(gamm) ; bmx(1,3) = 0.
    bmx(2,1) =-SIN(gamm) ; bmx(2,2) = COS(gamm) ; bmx(2,3) = 0.
    bmx(3,1) = 0.        ; bmx(3,2) = 0.        ; bmx(3,3) = 1.

    hmx(:,:) = 0.
    DO i = 1,3
       DO j = 1,3
          DO ii = 1,3
             hmx(i,j) = hmx(i,j) + cmx(i,ii)*dmx(ii,j)
          ENDDO
       ENDDO
    ENDDO

    DO i = 1,3
       DO j = 1,3
          DO ii = 1,3
             amx(i,j,1) = amx(i,j,1) + bmx(i,ii)*hmx(ii,j)
          ENDDO
       ENDDO
    ENDDO

  END SUBROUTINE euler

END MODULE m_abcrot2
