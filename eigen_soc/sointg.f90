MODULE m_sointg
  !*********************************************************************
  !     compute radial spin-orbit integrant
  !*********************************************************************
CONTAINS
  SUBROUTINE sointg(ntyp,e,vr,v0,atoms,input, vso)
    !
    USE m_differentiate,ONLY:diff3
    USE m_constants,ONLY: c_light
    USE m_types
    IMPLICIT NONE

    TYPE(t_input),INTENT(IN)   :: input
    TYPE(t_atoms),INTENT(IN)   :: atoms
    !     ..
    !     .. Scalar Arguments ..
    INTEGER,INTENT(IN)   :: ntyp
    REAL,    INTENT (IN) :: e
    !     ..
    !     .. Array Arguments ..
    REAL, INTENT (IN)  ::  v0(atoms%jri(ntyp))
    REAL, INTENT (IN)  ::  vr(atoms%jmtd,input%jspins)
    REAL, INTENT (OUT) :: vso(atoms%jmtd,2)
    !     ..
    !     .. Local Scalars ..
    REAL dr,dxv,r,cin2
    INTEGER i,jspin
    !     ..
    !     .. Local Arrays ..
    REAL dv(atoms%jri(ntyp)),xmrel(atoms%jri(ntyp))
    !     ..
    !     ..
    cin2 = 1.0/c_light(1.0)**2
    dr = EXP(atoms%dx(ntyp))
    r = atoms%rmsh(1,ntyp)
    !
    !---> relativistic mass:
    DO i = 1,atoms%jri(ntyp)
       xmrel(i) = 1. + (e-v0(i)/r)/2.*cin2
       r = r*dr
    END DO
    !---> potential derivative (on logarithmic mesh) : v0 := r*v

    CALL diff3(vr(:,1),atoms%dx(ntyp),dv)
    !
    r = atoms%rmsh(1,ntyp)
    DO i = 1,atoms%jri(ntyp)
       dxv = (dv(i) - vr(i,1))/(r**3)
       vso(i,2) = cin2*dxv/(4.*xmrel(i)**2)
       r = r*dr
    END DO

    IF (input%jspins.EQ.2) THEN
       CALL diff3(vr(:,input%jspins),atoms%dx(ntyp), dv)
       !
       r = atoms%rmsh(1,ntyp)
       DO i = 1,atoms%jri(ntyp)
          dxv = (dv(i) - vr(i,input%jspins))/(r**3)
          vso(i,1) = cin2*dxv/(4.*xmrel(i)**2)
          r = r*dr
       END DO
    ELSE 
       DO i = 1,atoms%jri(ntyp)
          vso(i,1) =  vso(i,2)
       END DO
    ENDIF

  END SUBROUTINE sointg
END MODULE m_sointg
