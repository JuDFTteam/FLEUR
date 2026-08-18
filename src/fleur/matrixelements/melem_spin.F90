!--------------------------------------------------------------------------------
! Copyright (c) 2026 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------
!>  What acts on the four spin-block overlaps once they exist: the Pauli components and
!>  the sum rule. The blocks themselves come from t_matelements_spin, which is the provider;
!>  this module is not one. On-site (b = 0, same k):
!>
!>      S0_alpha,mn(k) = < psi_mk | sigma_alpha | psi_nk > ,   alpha = x,y,z
!>
!>  built IN MEMORY from the two spinor components of the library-mode
!>  wavefunctions (no updown.mmn0 on disk). It is the Bloch-basis input O^(0)
!>  of the operator-interpolation pipeline: the driver rotates it to the Wannier
!>  gauge (V^dagger S0 V) and hands each Cartesian component to the generic core
!>  m_melem_ft.
!>
!>  The four spin-block overlaps  o_ab(m,n) = <phi^a_m|phi^b_n>  (a,b = global
!>  spin up=1/dn=2) are assembled by the spin operator itself, which contracts
!>  abc%cof with radfun%integral site by site. What is left here is what acts on
!>  those blocks once they exist: the Pauli components and the sum rule.
MODULE m_melem_spin
  USE m_juDFT
  USE m_constants, ONLY : ImagUnit, oUnit
  USE m_types_atoms
  USE m_types_abc
  USE m_types_radfun
  USE m_types_spinor_layout, ONLY: radial_slot
  USE m_types_nococonv
  USE m_types_stars
  USE m_types_lapw
  USE m_types_mat
  IMPLICIT NONE
  PRIVATE
  PUBLIC :: melem_pauli_from_blocks, melem_spin_sumrule
CONTAINS

  !> Assemble the three Pauli matrices at one k from the four global spin-block
  !> overlaps (interstitial + MT already summed into o_ab):
  !>   S_z = o_uu - o_dd ;  S_x = o_ud + o_du ;  S_y = -i (o_ud - o_du)
  SUBROUTINE melem_pauli_from_blocks(o_uu, o_dd, o_ud, o_du, s0)
    COMPLEX, INTENT(IN)  :: o_uu(:, :), o_dd(:, :), o_ud(:, :), o_du(:, :)  ! (nb,nb)
    COMPLEX, INTENT(OUT) :: s0(:, :, :)                                     ! (nb,nb,3) 1=x 2=y 3=z

    s0(:, :, 1) = o_ud + o_du                    ! sigma_x
    s0(:, :, 2) = -ImagUnit * (o_ud - o_du)      ! sigma_y
    s0(:, :, 3) = o_uu - o_dd                    ! sigma_z
  END SUBROUTINE melem_pauli_from_blocks


  !> Sum-rule / sanity check on the Bloch-basis spin matrices at one k.
  !> With interstitial + MT summed, spin-trace orthonormality gives
  !>   o_uu(m,m) + o_dd(m,m) = 1  (norm), and the per-band spin |<sigma>_m| <= 1.
  !> For a magnet with the moment along z: <sigma_z>_m ~ +/-1, <sigma_xy>_m ~ 0.
  SUBROUTINE melem_spin_sumrule(s0, o_uu, o_dd, ik, tol)
    COMPLEX, INTENT(IN) :: s0(:, :, :)         ! (nb,nb,3)
    COMPLEX, INTENT(IN) :: o_uu(:, :), o_dd(:, :)
    INTEGER, INTENT(IN) :: ik
    REAL,    INTENT(IN) :: tol

    INTEGER :: nb, m, nbad
    REAL    :: nrm, sx, sy, sz, smag

    nb = SIZE(s0, 1); nbad = 0
    WRITE(oUnit, '(a,i0)') 'wannierlib spin sum-rule check, k = ', ik
    WRITE(oUnit, '(a)')    '  band     norm    <sx>     <sy>     <sz>    |<s>|'
    DO m = 1, nb
      nrm = REAL(o_uu(m, m) + o_dd(m, m))
      sx  = REAL(s0(m, m, 1)); sy = REAL(s0(m, m, 2)); sz = REAL(s0(m, m, 3))
      smag = SQRT(sx*sx + sy*sy + sz*sz)
      IF (ABS(nrm - 1.0) > tol .OR. smag > 1.0 + tol) nbad = nbad + 1
      WRITE(oUnit, '(i6,5f9.4)') m, nrm, sx, sy, sz, smag
    END DO
    IF (nbad == 0) THEN
      WRITE(oUnit, '(a)') '  sum-rule OK: norms=1 and |<s>|<=1 for all bands'
    ELSE
      WRITE(oUnit, '(a,i0,a)') '  WARNING: ', nbad, ' bands violate norm=1 or |<s>|<=1'
    END IF
  END SUBROUTINE melem_spin_sumrule

END MODULE m_melem_spin
