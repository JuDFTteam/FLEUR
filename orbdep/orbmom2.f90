MODULE m_orbmom2
  !     ***************************************************************
  !     perform the sum over m (for each l) and calculate the
  !     spherical contribution to orbital moment.                
  !     ***************************************************************
  !
CONTAINS
  SUBROUTINE orbmom2(atoms,itype,ispin,ddn,orb,uulon,dulon,uloulopn,clmom)

    !      USE m_types, ONLY : t_orb,t_orbl,t_orblo
    USE m_types
    USE m_constants
    IMPLICIT NONE

    TYPE(t_atoms),INTENT(IN)   :: atoms
    !     ..
    !     .. Scalar Arguments ..
    INTEGER, INTENT (IN) :: itype, ispin
    !     ..
    !     .. Array Arguments ..
    REAL,    INTENT (IN) :: ddn(0:atoms%lmaxd),uulon(atoms%nlod),dulon(atoms%nlod)
    REAL,    INTENT (IN) :: uloulopn(atoms%nlod,atoms%nlod)
    TYPE (t_orb),  INTENT (IN) :: orb
    REAL,    INTENT (OUT) :: clmom(3)
    !     ..
    !     .. Local Scalars ..
    INTEGER l , ilo, ilop,m
    REAL qmtt, qmttx, qmtty, sumlm
    COMPLEX orbp, orbm
    !     ..
    !     .. Local Arrays ..
    REAL qmtl(0:atoms%lmaxd),qmtlx(0:atoms%lmaxd),qmtly(0:atoms%lmaxd)

    qmtt = 0.
    qmttx = 0.
    qmtty = 0.
    DO l = 0,atoms%lmax(itype)
       !--->    lm-decomposed density for each atom type
       qmtl(l) = 0.
       qmtlx(l) = 0.
       qmtly(l) = 0.
       DO m = -l,l
          ! lz
          sumlm = m * (orb%uu(l,m,itype,ispin) + orb%dd(l,m,itype,ispin) * ddn(l) ) 
          ! lx,ly
          orbp = SQRT(REAL((l-m)*(l+m+1))) * ( orb%uup(l,m,itype,ispin) + orb%ddp(l,m,itype,ispin) * ddn(l) ) 

          orbm = SQRT(REAL((l+m)*(l-m+1))) * ( orb%uum(l,m,itype,ispin) + orb%ddm(l,m,itype,ispin) * ddn(l) )
          !+gu
          IF (m.EQ.l)  orbp = CMPLX(0.0,0.0)
          IF (m.EQ.-l) orbm = CMPLX(0.0,0.0)
          !+gu
          qmtl(l)  = qmtl(l)  + sumlm
          qmtlx(l) = qmtlx(l) + 0.5*( REAL(orbp)+ REAL(orbm))
          qmtly(l) = qmtly(l) + 0.5*(AIMAG(orbp)-AIMAG(orbm))
          ! 
       ENDDO
    ENDDO
    !
    ! --> LO contribution
    DO ilo = 1, atoms%nlo(itype)
       l = atoms%llo(ilo,itype)
       DO m = -l,l
          sumlm = m * (orb%uulo(ilo,m,itype,ispin) * uulon(ilo) + orb%dulo(ilo,m,itype,ispin) * dulon(ilo) )

          orbp = SQRT(REAL((l-m)*(l+m+1))) * ( orb%uulop(ilo,m,itype,ispin) * uulon(ilo) +&
               orb%dulop(ilo,m,itype,ispin) * dulon(ilo) )

          orbm = SQRT(REAL((l+m)*(l-m+1))) * ( orb%uulom(ilo,m,itype,ispin) * uulon(ilo) +&
               orb%dulom(ilo,m,itype,ispin) * dulon(ilo) )

          IF (m.EQ.l)  orbp = CMPLX(0.0,0.0)
          IF (m.EQ.-l) orbm = CMPLX(0.0,0.0)

          qmtl(l)  = qmtl(l)  + sumlm
          qmtlx(l) = qmtlx(l) + 0.5*( REAL(orbp)+ REAL(orbm))
          qmtly(l) = qmtly(l) + 0.5*(AIMAG(orbp)-AIMAG(orbm))
       ENDDO
       DO ilop = 1, atoms%nlo(itype)
          IF (atoms%llo(ilop,itype).EQ.l) THEN
             DO m = -l,l
                sumlm = m * orb%z(ilo,ilop,m,itype,ispin) * uloulopn(ilo,ilop)
                orbp = SQRT(REAL((l-m)*(l+m+1))) * orb%p(ilo,ilop,m,itype,ispin) * uloulopn(ilo,ilop)
                orbm = SQRT(REAL((l+m)*(l-m+1))) * orb%m(ilo,ilop,m,itype,ispin) * uloulopn(ilo,ilop)
                IF (m.EQ.l)  orbp = CMPLX(0.0,0.0)
                IF (m.EQ.-l) orbm = CMPLX(0.0,0.0)

                qmtl(l)  = qmtl(l)  + sumlm
                qmtlx(l) = qmtlx(l) + 0.5*( REAL(orbp)+ REAL(orbm))
                qmtly(l) = qmtly(l) + 0.5*(AIMAG(orbp)-AIMAG(orbm))
             ENDDO
          ENDIF
       ENDDO
    ENDDO
    !
    ! --> sum up & print
    DO l = 0,atoms%lmax(itype)
       qmtl(l)  = qmtl(l)  / atoms%neq(itype)
       qmtlx(l) = qmtlx(l) / atoms%neq(itype)
       qmtly(l) = qmtly(l) / atoms%neq(itype)
       qmtt =  qmtt  + qmtl(l)
       qmttx = qmttx + qmtlx(l)
       qmtty = qmtty + qmtly(l)
    ENDDO
    clmom(1) = qmttx
    clmom(2) = qmtty
    clmom(3) = qmtt

! The following output was commented out, because the subroutine is now  used in parallel.
! Jan. 2019   U.Alekseeva
!
!    WRITE (oUnit,FMT=8100) itype, (qmtl(l),l=0,3), qmtt
!    WRITE (oUnit,FMT=8100) itype, (qmtlx(l),l=0,3),qmttx
!    WRITE (oUnit,FMT=8100) itype, (qmtly(l),l=0,3),qmtty
!8100 FORMAT (' -->',i2,2x,4f9.5,2x,f9.5)

  END SUBROUTINE orbmom2
END MODULE m_orbmom2
