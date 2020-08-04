!--------------------------------------------------------------------------------
! Copyright (c) 2016 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------

MODULE m_lhglpts
  !     **********************************************************
  !     calculates lattice harmonics on the gauss-legendre angular
  !     mesh - r.pentcheva Feb'96
  !     **********************************************************
CONTAINS
  SUBROUTINE lhglpts(&
       &                   sphhar,atoms,&
       &                   rx,nsp,&
       &                   sym,&
       &                   ylh,ylhmh)
    !
    USE m_ylm
    USE m_types_sym
    USE m_types_sphhar
    USE m_types_atoms

    IMPLICIT NONE

    TYPE(t_sym),INTENT(IN)         :: sym
    TYPE(t_sphhar),INTENT(IN)      :: sphhar
    TYPE(t_atoms),INTENT(IN)       :: atoms
    !     ..
    !     .. Scalar Arguments ..
    INTEGER, INTENT (IN) :: nsp
    !     ..
    !     .. Array Arguments ..
    REAL,    INTENT (IN) :: rx(:,:) !(3,dimension%nspd)
    REAL,    INTENT (OUT):: ylh(:,0:,:) !(dimension%nspd,0:sphhar%nlhd,sphhar%ntypsd)
    COMPLEX, OPTIONAL, INTENT (OUT):: ylhmh(:,0:,:)
    !     ..
    !     .. Local Scalars ..
    REAL s
    INTEGER k,lh,mem,nd,ll1,lm
    !     ..
    !     .. Local Arrays ..
    COMPLEX ylm( (atoms%lmaxd+1)**2 )
    !     ..
    IF(PRESENT(ylhmh)) ylhmh = 0.0
    DO  nd = 1,sym%nsymt
       DO  k = 1,nsp

          CALL ylm4(&
               &                atoms%lmaxd,rx(:,k),&
               &                ylm)

          DO lh = 0,sphhar%nlh(nd)
             s = 0
             ll1 = sphhar%llh(lh,nd) * ( sphhar%llh(lh,nd) + 1 ) + 1
             DO mem = 1,sphhar%nmem(lh,nd)
                lm = ll1 + sphhar%mlh(mem,lh,nd)
                s = s + REAL( sphhar%clnu(mem,lh,nd) * ylm(lm) )
                IF(PRESENT(ylhmh)) ylhmh(k,lm-1,nd) = ylm(lm)
             ENDDO
             ylh(k,lh,nd) = s
          ENDDO

       ENDDO
    ENDDO
    RETURN
  END SUBROUTINE lhglpts
END MODULE m_lhglpts
