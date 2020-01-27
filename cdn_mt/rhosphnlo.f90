!--------------------------------------------------------------------------------
! Copyright (c) 2016 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------

MODULE m_rhosphnlo
  !***********************************************************************
  ! Add the local orbital contributions to the charge density. The 
  ! corresponding summation of the pure apw contribuions is done in
  ! cdnval.
  ! Philipp Kurz 99/04
  !***********************************************************************
CONTAINS
  SUBROUTINE rhosphnlo(itype,atoms,sphhar, uloulopn,dulon,uulon,&
       ello,vr, aclo,bclo,cclo,acnmt,bcnmt,ccnmt,f,g, rho,qmtllo)

    USE m_constants, ONLY : c_light,sfp_const
    USE m_radsra
    USE m_radsrdn
    USE m_types
    IMPLICIT NONE
    TYPE(t_sphhar),INTENT(IN)   :: sphhar
    TYPE(t_atoms),INTENT(IN)    :: atoms
    !     ..
    !     .. Scalar Arguments ..
    INTEGER,    INTENT (IN) :: itype 
    !     ..
    !     .. Array Arguments ..
    REAL,    INTENT (IN) :: aclo(atoms%nlod),bclo(atoms%nlod),cclo(atoms%nlod,atoms%nlod)
    REAL,    INTENT (IN) :: acnmt(0:atoms%lmaxd,atoms%nlod,sphhar%nlhd)
    REAL,    INTENT (IN) :: bcnmt(0:atoms%lmaxd,atoms%nlod,sphhar%nlhd)
    REAL,    INTENT (IN) :: ccnmt(atoms%nlod,atoms%nlod,sphhar%nlhd)
    REAL,    INTENT (IN) :: dulon(atoms%nlod),uulon(atoms%nlod),vr(atoms%jmtd)
    REAL,    INTENT (IN) :: uloulopn(atoms%nlod,atoms%nlod),ello(atoms%nlod)
    REAL,    INTENT (IN) :: f(atoms%jmtd,2,0:atoms%lmaxd),g(atoms%jmtd,2,0:atoms%lmaxd)
    REAL,    INTENT (INOUT) :: qmtllo(0:atoms%lmaxd)
    REAL,    INTENT (INOUT) :: rho(atoms%jmtd,0:sphhar%nlhd)
    !     ..
    !     .. Local Scalars ..
    REAL dsdum,usdum ,c_1,c_2
    INTEGER j,l,lh,lo,lop,lp,nodedum
    REAL dus,ddn,c
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
       qmtllo(l) = qmtllo(l) + (aclo(lo)*uulon(lo) +bclo(lo)*dulon(lo)) * c_1
       DO lop = 1,atoms%nlo(itype)
          IF (atoms%llo(lop,itype).EQ.l) THEN
             qmtllo(l) = qmtllo(l) + (cclo(lop,lo) *uloulopn(lop,lo)) * c_1
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
       DO j = 1,atoms%jri(itype)
          rho(j,0) = rho(j,0) + c_2 *&
               (aclo(lo) * ( f(j,1,l)*flo(j,1,lo) +f(j,2,l)*flo(j,2,lo) ) +&
               bclo(lo) * ( g(j,1,l)*flo(j,1,lo) +g(j,2,l)*flo(j,2,lo) ) )
       END DO
       DO lop = 1,atoms%nlo(itype)
          IF (atoms%llo(lop,itype).EQ.l) THEN
             DO j = 1,atoms%jri(itype)
                rho(j,0) = rho(j,0) + c_2 * cclo(lop,lo) *&
                     ( flo(j,1,lop)*flo(j,1,lo) +flo(j,2,lop)*flo(j,2,lo) )
             END DO
          END IF
       END DO
    END DO

    !---> add the contribution of the local orbitals and flapw - lo cross-
    !---> terms to the non-spherical chargedensity inside the muffin tins.

    DO lh = 1,sphhar%nlh(atoms%ntypsy(atoms%nat))
       DO lp = 0,atoms%lmax(itype)
          DO lo = 1,atoms%nlo(itype)
             DO j = 1,atoms%jri(itype)
                rho(j,lh) = rho(j,lh) + c_1 * (&
                     acnmt(lp,lo,lh) * (f(j,1,lp)*flo(j,1,lo) +f(j,2,lp)*flo(j,2,lo) ) +&
                     bcnmt(lp,lo,lh) * (g(j,1,lp)*flo(j,1,lo) +g(j,2,lp)*flo(j,2,lo) ) )
             END DO
          END DO
       END DO
       DO lo = 1,atoms%nlo(itype)
          DO lop = 1,atoms%nlo(itype)
             DO j = 1,atoms%jri(itype)
                rho(j,lh) = rho(j,lh) + c_1 * ccnmt(lop,lo,lh) *&
                     ( flo(j,1,lop)*flo(j,1,lo) +flo(j,2,lop)*flo(j,2,lo) )
             END DO
          END DO
       END DO
    END DO
    DEALLOCATE (flo,glo)

  END SUBROUTINE rhosphnlo
END MODULE m_rhosphnlo
