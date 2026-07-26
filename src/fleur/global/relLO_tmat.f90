!--------------------------------------------------------------------------------
! Copyright (c) 2026 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------

MODULE m_relLO_tmat
  !********************************************************************************
  ! SPHERICAL scalar-relativistic Hamiltonian couplings of a confined relativistic 
  ! local orbital (relLO) to the APW radial functions u and u-dot:
  !
  !   t_uulo = <u   | H_sph | d0>   (radial function d0 = the relLO's Dirac solution)
  !   t_dulo = <udot| H_sph | d0>
  !
  ! tlo.f90 uses them for the relativistic local orbital in the (first-variation)
  ! muffin-tin Hamiltonian: t_uulo goes into tuulo, t_dulo into tdulo.
  !
  ! Derivation: for a function pair (a, d0) with a in {u,udot} an SRA eigenfunction
  ! and d0 NOT an SRA eigenfunction, Hermiticity of H_sph gives the directed element
  !   <a|H_sph|d0> = <H_sph a|d0> + W(a,d0) = el*<a|d0> + W(a,d0)
  ! with the kinetic-energy Wronskian surface term at the muffin-tin radius
  !   W(a,d0) = -0.5 * R_mt^2 * ( a(R)*d0'(R) - a'(R)*d0(R) ).
  ! FLEUR builds a Hermitian Hamiltonian by symmetrizing the (a,LO) and (LO,a)
  ! blocks, which keeps HALF of this boundary term on each side, so the stored
  ! element carries 0.5*W. (For an ORDINARY LO the Wronskian identity
  ! W(a,LO)=(ello-el)*<a|LO> collapses this to the usual 0.5*(el+ello)*<a|LO>
  ! shortcut; for the Dirac relLO that identity fails and W must be used explicitly.)
  !********************************************************************************
  IMPLICIT NONE
  PRIVATE
  PUBLIC :: relLO_sph_tmat
CONTAINS
  PURE SUBROUTINE relLO_sph_tmat(atoms,usdus,enpara,ntyp,l,lo,jsp, t_uulo,t_dulo)
    USE m_types
    TYPE(t_atoms),  INTENT(IN)  :: atoms
    TYPE(t_usdus),  INTENT(IN)  :: usdus
    TYPE(t_enpara), INTENT(IN)  :: enpara
    INTEGER,        INTENT(IN)  :: ntyp,l,lo,jsp
    REAL,           INTENT(OUT) :: t_uulo,t_dulo

    REAL :: rmt2,wu,wudot,el

    rmt2  = atoms%rmt(ntyp)**2
    el    = enpara%el0(l,ntyp,jsp)
    wu    = -0.5*rmt2*( usdus%us(l,ntyp,jsp)*usdus%dulos(lo,ntyp,jsp) &
                      - usdus%dus(l,ntyp,jsp)*usdus%ulos(lo,ntyp,jsp) )
    wudot = -0.5*rmt2*( usdus%uds(l,ntyp,jsp)*usdus%dulos(lo,ntyp,jsp) &
                      - usdus%duds(l,ntyp,jsp)*usdus%ulos(lo,ntyp,jsp) )
    t_uulo = el*usdus%uulon(lo,ntyp,jsp) + 0.5*wu
    t_dulo = el*usdus%dulon(lo,ntyp,jsp) + usdus%uulon(lo,ntyp,jsp) + 0.5*wudot
  END SUBROUTINE relLO_sph_tmat
END MODULE m_relLO_tmat
