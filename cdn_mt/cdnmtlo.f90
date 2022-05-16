!--------------------------------------------------------------------------------
! Copyright (c) 2016 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------

MODULE m_cdnmtlo
  !***********************************************************************
  ! Add the local orbital contributions to the charge density. The
  ! corresponding summation of the pure apw contribuions is done in
  ! cdnval.
  ! Philipp Kurz 99/04
  !***********************************************************************
CONTAINS
  SUBROUTINE cdnmtlo(itype,ilSpinPr,ilSpin,input,atoms,sphhar,sym, uloulopn,dulon,uulon,&
       ello,vr, denCoeffs, f,g, rho,moments,qmtllo, rhoIm)

    USE m_constants, ONLY : c_light,sfp_const
    USE m_types
    USE m_radsra
    USE m_radsrdn

    IMPLICIT NONE

    TYPE(t_input),INTENT(IN)    :: input
    TYPE(t_sphhar),INTENT(IN)   :: sphhar
    TYPE(t_atoms),INTENT(IN)    :: atoms
    TYPE(t_sym),INTENT(IN)      :: sym
    TYPE (t_denCoeffs),        INTENT(IN) :: denCoeffs

    !     ..
    !     .. Scalar Arguments ..
    INTEGER,    INTENT (IN) :: itype, ilSpinPr, ilSpin
    !     ..
    !     .. Array Arguments ..
    REAL,    INTENT (IN) :: dulon(:),uulon(:),vr(:)
    REAL,    INTENT (IN) :: uloulopn(:,:),ello(:)
    REAL,    INTENT (IN) :: f(:,:,0:),g(:,:,0:)
    REAL,    INTENT (INOUT) :: qmtllo(0:)
    REAL,    INTENT (INOUT) :: rho(:,0:)
    REAL, OPTIONAL, INTENT(INOUT) :: rhoIm(:,0:)
    TYPE(t_moments), INTENT(INOUT) :: moments

    INTEGER, PARAMETER :: lcf=3
    !     ..
    !     .. Local Scalars ..
    REAL dsdum,usdum ,c_1,c_2
    INTEGER j,l,lh,lo,lop,lp,nodedum,llp
    REAL dus,ddn,c
    COMPLEX :: ctemp
    !     ..
    !     .. Local Arrays ..
    REAL,    ALLOCATABLE :: flo(:,:,:),glo(:,:)
    REAL filo(atoms%jmtd,2)
    !     ..
    c = c_light(1.0)
    c_1 = 1.0 / atoms%neq(itype)
    c_2 = 1.0 /(atoms%neq(itype)*sfp_const)
    !
    DO lo = 1,atoms%nlo(itype)
       l = atoms%llo(lo,itype)
       ctemp = (denCoeffs%mt_ulo_coeff(lo,itype,0,ilSpinPr,ilSpin)+denCoeffs%mt_lou_coeff(lo,itype,0,ilSpinPr,ilSpin)) &
           & * uulon(lo) &
           & + (denCoeffs%mt_ulo_coeff(lo,itype,1,ilSpinPr,ilSpin)+denCoeffs%mt_lou_coeff(lo,itype,1,ilSpinPr,ilSpin)) &
           & * dulon(lo)
       qmtllo(l) = qmtllo(l) + REAL(ctemp) * c_1
       DO lop = 1,atoms%nlo(itype)
          IF (atoms%llo(lop,itype).EQ.l) THEN
             ctemp = denCoeffs%mt_lolo_coeff(lop,lo,itype,ilSpinPr,ilSpin) * uloulopn(lop,lo)
             qmtllo(l) = qmtllo(l) + REAL(ctemp) * c_1
          END IF
       END DO
    END DO
    ALLOCATE ( flo(atoms%jmtd,2,atoms%nlod),glo(atoms%jmtd,2) )

    !---> calculate the local ortital radial functions

    DO lo = 1,atoms%nlo(itype)
       l = atoms%llo(lo,itype)
       CALL radsra(ello(lo),l,vr,atoms%rmsh(1,itype),atoms%dx(itype),atoms%jri(itype),c,&
            usdum,dus,nodedum,flo(:,1,lo),flo(:,2,lo))
       !+apw+lo
       IF (atoms%l_dulo(lo,itype).or.atoms%ulo_der(lo,itype).ge.1) THEN
          !--->    calculate orthogonal energy derivative at e
          j = atoms%ulo_der(lo,itype)
          IF(atoms%l_dulo(lo,itype)) j = 1
          CALL radsrdn(ello(lo),l,vr,atoms%rmsh(1,itype),atoms%dx(itype),atoms%jri(itype),c,&
               usdum,dsdum,ddn,nodedum,glo,filo,flo(:,:,lo),dus,j) ! filo is a dummy array&
          DO j=1,atoms%jri(itype)
             flo(j,1,lo) = glo(j,1)
             flo(j,2,lo) = glo(j,2)
          ENDDO
          ddn = sqrt(ddn)
          IF(atoms%l_dulo(lo,itype)) ddn=1.0
          flo(:,:,lo) = flo(:,:,lo)/ddn ! Normalize ulo (flo) if APW+lo is not used
       ENDIF
       !-apw+lo
    END DO

    !---> add the contribution of the local orbitals and flapw - lo cross-
    !---> terms to the spherical chargedensity inside the muffin tins.

    DO lo = 1,atoms%nlo(itype)
       l = atoms%llo(lo,itype)
       llp = (l* (l+1))/2 + l
       DO j = 1,atoms%jri(itype)
          ctemp = (denCoeffs%mt_ulo_coeff(lo,itype,0,ilSpinPr,ilSpin)+denCoeffs%mt_lou_coeff(lo,itype,0,ilSpinPr,ilSpin)) &
                * (f(j,1,l)*flo(j,1,lo)+f(j,2,l)*flo(j,2,lo)) &
                + (denCoeffs%mt_ulo_coeff(lo,itype,1,ilSpinPr,ilSpin)+denCoeffs%mt_lou_coeff(lo,itype,1,ilSpinPr,ilSpin)) &
                * (g(j,1,l)*flo(j,1,lo)+g(j,2,l)*flo(j,2,lo))
          rho(j,0) = rho(j,0) + c_2 * REAL(ctemp)
          IF (PRESENT(rhoIm)) rhoIm(j,0) = rhoIm(j,0) + c_2 * AIMAG(ctemp)
          IF (l.LE.input%lResMax.AND.ilSpinPr.EQ.ilSpin) THEN
             moments%rhoLRes(j,0,llp,itype,ilSpin) = moments%rhoLRes(j,0,llp,itype,ilSpin) + c_2 * REAL(ctemp)
          END IF
       END DO
       DO lop = 1,atoms%nlo(itype)
          IF (atoms%llo(lop,itype).EQ.l) THEN
             DO j = 1,atoms%jri(itype)
                ctemp = c_2 * denCoeffs%mt_lolo_coeff(lop,lo,itype,ilSpinPr,ilSpin) &
                     * (flo(j,1,lop)*flo(j,1,lo)+flo(j,2,lop)*flo(j,2,lo))
                rho(j,0) = rho(j,0) + REAL(ctemp)
                IF (PRESENT(rhoIm)) rhoIm(j,0) = rhoIm(j,0) + AIMAG(ctemp)
                IF (l.LE.input%lResMax.AND.ilSpinPr.EQ.ilSpin) THEN
                   moments%rhoLRes(j,0,llp,itype,ilSpin) = moments%rhoLRes(j,0,llp,itype,ilSpin) + REAL(ctemp)
                END IF
             END DO
          END IF
       END DO
    END DO

    !---> add the contribution of the local orbitals and flapw - lo cross-
    !---> terms to the non-spherical chargedensity inside the muffin tins.

    DO lh = 1,sphhar%nlh(sym%ntypsy(sum(atoms%neq(:itype-1))+1))
       DO lp = 0,atoms%lmax(itype)
          DO lo = 1,atoms%nlo(itype)
             l = atoms%llo(lo,itype)
             IF(atoms%l_outputCFpot(itype).AND.atoms%l_outputCFremove4f(itype)&
                .AND.(l.EQ.lcf.AND.lp.EQ.lcf)) CYCLE !Exclude non-spherical contributions for CF
             llp = (MAX(l,lp)* (MAX(l,lp)+1))/2 + MIN(l,lp)
             DO j = 1,atoms%jri(itype)
                ctemp = c_1 &
                      * ((denCoeffs%nmt_ulo_coeff(lp,lo,lh,itype,0,ilSpinPr,ilSpin)+denCoeffs%nmt_lou_coeff(lp,lo,lh,itype,0,ilSpinPr,ilSpin)) &
                      * (f(j,1,lp)*flo(j,1,lo)+f(j,2,lp)*flo(j,2,lo)) &
                      +  (denCoeffs%nmt_ulo_coeff(lp,lo,lh,itype,1,ilSpinPr,ilSpin)+denCoeffs%nmt_lou_coeff(lp,lo,lh,itype,1,ilSpinPr,ilSpin)) &
                      * (g(j,1,lp)*flo(j,1,lo)+g(j,2,lp)*flo(j,2,lo)))
                rho(j,lh) = rho(j,lh) + REAL(ctemp)
                IF (PRESENT(rhoIm)) rhoIm(j,0) = rhoIm(j,0) + AIMAG(ctemp)
                IF ((l.LE.input%lResMax).AND.(lp.LE.input%lResMax).AND.ilSpinPr.EQ.ilSpin) THEN
                   moments%rhoLRes(j,lh,llp,itype,ilSpin) = moments%rhoLRes(j,lh,llp,itype,ilSpin) + REAL(ctemp)
                END IF
             END DO
          END DO
       END DO
       DO lo = 1,atoms%nlo(itype)
          l = atoms%llo(lo,itype)
          DO lop = 1,atoms%nlo(itype)
             lp = atoms%llo(lop,itype)
             llp = (MAX(l,lp)* (MAX(l,lp)+1))/2 + MIN(l,lp)
             DO j = 1,atoms%jri(itype)
                ctemp = c_1 * denCoeffs%nmt_lolo_coeff(lop,lo,lh,itype,ilSpinPr,ilSpin) * (flo(j,1,lop)*flo(j,1,lo)+flo(j,2,lop)*flo(j,2,lo))
                rho(j,lh) = rho(j,lh) + REAL(ctemp)
                IF (PRESENT(rhoIm)) rhoIm(j,0) = rhoIm(j,0) + AIMAG(ctemp)
                IF ((l.LE.input%lResMax).AND.(lp.LE.input%lResMax).AND.ilSpinPr.EQ.ilSpin) THEN
                   moments%rhoLRes(j,lh,llp,itype,ilSpin) = moments%rhoLRes(j,lh,llp,itype,ilSpin) + REAL(ctemp)
                END IF
             END DO
          END DO
       END DO
    END DO
    DEALLOCATE (flo,glo)

  END SUBROUTINE cdnmtlo
END MODULE m_cdnmtlo
