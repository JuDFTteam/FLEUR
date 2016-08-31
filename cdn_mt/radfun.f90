MODULE m_radfun
  USE m_juDFT
CONTAINS
  SUBROUTINE radfun(l,itype,jsp,e,vr,atoms, f,g,usdus,nodeu,noded,wronk)
    !*********************************************************************
    !     generates the scalar relativistic wavefunctions (function: f;
    !     energy derivative: g) at an energy e for angular momentum l.
    !     the values on the sphere boundaries are also returned.
    !             m. weinert   jan. 1987
    !     the solutions r*u(r) are on a log. mesh.
    !
    !      us ... u(R) if R is the muffin tin Radius
    !     dus ... u'(R)   (radial derivative)
    !             .               .
    !     uds ... u(R)   duds ... u'(R)  (energy derivative)
    !              . .
    !     ddn ... <u|u>  norm of u-dot
    !
    !*********************************************************************

    USE m_constants, ONLY : c_light
    USE m_radsra
    USE m_radsrd
    USE m_types
    IMPLICIT NONE
    TYPE(t_usdus),INTENT(INOUT):: usdus
    TYPE(t_atoms),INTENT(IN)   :: atoms
    !     ..
    !     .. Scalar Arguments ..
    INTEGER, INTENT (IN) :: l,itype,jsp
    INTEGER, INTENT (OUT):: noded,nodeu
    REAL,    INTENT (IN) :: e
    REAL,    INTENT (OUT):: wronk
    !     ..
    !     .. Array Arguments ..
    REAL,    INTENT (IN) :: vr(atoms%jmtd)
    REAL,    INTENT (OUT):: f(atoms%jmtd,2),g(atoms%jmtd,2)
    !     ..
    !     ..
    REAL :: c,us,dus
    IF (atoms%jri(itype)>atoms%jmtd)  CALL juDFT_error("atoms%jmtd too small","radfun")
    !
    c = c_light(1.0)
    !
    !--->    calculate normalized function at e
    CALL radsra(e,l,vr,atoms%rmsh(1,itype),atoms%dx(itype),atoms%jri(itype),&
         c, usdus%us(l,itype,jsp),usdus%dus(l,itype,jsp),nodeu,f(:,1),f(:,2))
    !
    !--->    calculate orthogonal energy derivative at e
    CALL radsrd(e,l,vr,atoms%rmsh(1,itype),atoms%dx(itype),atoms%jri(itype),&
         c, usdus%uds(l,itype,jsp),usdus%duds(l,itype,jsp),usdus%ddn(l,itype,jsp),&
         noded,g(:,1),g(:,2), f(:,1),f(:,2),usdus%dus(l,itype,jsp))
    !
    !--->    calculate wronskian
    wronk = usdus%uds(l,itype,jsp)*usdus%dus(l,itype,jsp) - usdus%duds(l,itype,jsp)*usdus%us(l,itype,jsp)

  END SUBROUTINE radfun
END MODULE m_radfun
