!--------------------------------------------------------------------------------
! Copyright (c) 2026 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------

MODULE m_relLO_dirac
  !********************************************************************************
  ! Subroutine for relativistic local orbital's (relLO) radial function: the large/small 
  ! components (P,Q) of the Dirac j=l-1/2 bound state d0, solved on a mesh extended 
  ! beyond the muffin-tin radius, returned normalized on the muffin-tin mesh 
  ! together with the large component's value and slope at R_mt.
  !********************************************************************************
  USE m_juDFT
  IMPLICIT NONE
  PRIVATE
  PUBLIC :: relLO_dirac_radial
CONTAINS
  SUBROUTINE relLO_dirac_radial(atoms,itype,l,nqn,e,vr, p,q,us,dus)
    USE m_constants,      ONLY : c_light
    USE m_differ
    USE m_intgr,          ONLY : intgr0
    USE m_differentiate,  ONLY : difcub
    USE m_types_atoms
    TYPE(t_atoms),INTENT(IN)    :: atoms
    INTEGER,      INTENT(IN)    :: itype,l,nqn
    REAL,         INTENT(INOUT) :: e                  ! energy guess in, converged Dirac eigenvalue out
    REAL,         INTENT(IN)    :: vr(atoms%jmtd)     ! spherical potential r*V on the muffin-tin mesh
    REAL,         INTENT(OUT)   :: p(atoms%jmtd),q(atoms%jmtd) ! normalized large/small comp. (1:jri filled)
    REAL,         INTENT(OUT)   :: us,dus             ! large-component value and slope at R_mt

    REAL, ALLOCATABLE :: vrd(:),fint(:),f_rel(:,:)
    REAL    :: d,rn,rmt,t1,t2,rr,fn,fl,fj,ddn,dp_dr
    INTEGER :: j,jri,msh,ierr

    jri = atoms%jri(itype)

    ! extended mesh: differ's bound-state search integrates a bit beyond the muffin-tin
    ! radius, even though only the 1:jri part is kept.
    d   = EXP(atoms%dx(itype))
    rmt = atoms%rmsh(1,itype) * d**(jri-1)
    rn  = rmt
    msh = jri
    DO WHILE (rn < rmt + 20.0)
       msh = msh + 1
       rn  = rn*d
    END DO
    rn = atoms%rmsh(1,itype) * d**(msh-1)
    ALLOCATE(vrd(msh),f_rel(msh,2))

    ! extend the potential beyond the muffin-tin radius (linear with slope t1 / a.u.)
    vrd(1:jri) = vr(1:jri)
    t1 = 0.125
    t2 = vrd(jri)/rmt - rmt*t1
    rr = rmt
    DO j = jri+1, msh
       rr     = d*rr
       vrd(j) = rr*( t2 + rr*t1 )
    END DO

    fn = REAL(nqn) ; fl = REAL(l) ; fj = fl - 0.5
    CALL differ(fn,fl,fj,c_light(1.0),atoms%zatom(itype),atoms%dx(itype),atoms%rmsh(1,itype),&
         rn,d,msh,vrd, e, f_rel(:,1),f_rel(:,2),ierr)
    IF (ierr/=0) CALL juDFT_error("Dirac solve for relLO did not converge",calledby="relLO_dirac_radial")
    p(1:jri) = f_rel(1:jri,1)
    q(1:jri) = f_rel(1:jri,2)
    DEALLOCATE(vrd,f_rel)

    ALLOCATE(fint(jri))
    fint = p(1:jri)**2 + q(1:jri)**2
    CALL intgr0(fint,atoms%rmsh(1,itype),atoms%dx(itype),jri,ddn)
    ddn = SQRT(ddn)
    p(1:jri) = p(1:jri)/ddn
    q(1:jri) = q(1:jri)/ddn
    DEALLOCATE(fint)

    ! value & slope of the large component R(r)=P(r)/r at R_mt, obtained numerically since
    ! the Dirac d0 does not satisfy radsra's closed-form us/dus relation. P stores r*R(r), so
    ! dR/dr = (dP/dr)/r - P/r**2 = (dP/dr - us)/r .
    rmt   = atoms%rmsh(jri,itype)
    us    = p(jri)/rmt
    dp_dr = difcub(atoms%rmsh(jri-3,itype),p(jri-3),rmt)
    dus   = (dp_dr - us)/rmt
  END SUBROUTINE relLO_dirac_radial
END MODULE m_relLO_dirac
