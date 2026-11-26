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
!>  m_wgauge_ft.
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
  PUBLIC :: melem_spin_channel_check
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

  !> The within-channel spin block of a collinear run, against what it has to be.
  !>
  !> Over one spin channel sigma_z reduces to that channel's own overlap, so the block is
  !> the identity -- but it IS the integral, taken over the spheres, the interstitial and
  !> the vacuum like every other matrix element here, not a number written down. That is
  !> what makes this a check and not a tautology: a region left out of the sum, a bad
  !> normalisation, or a film whose vacuum never reached the operator all move the
  !> diagonal off one, and none of them would move a constant.
  !>
  !> The magnitude is tested and not the sign, because the sign is the channel label and
  !> is assigned rather than measured: channel 1 is up.
  !>
  !> Every k is judged; only l_report prints. The two are separate because they reach the
  !> reader by different routes. An abort leaves any rank, but a WRITE to oUnit only
  !> reaches `out` from rank 0, and the k are handed out in strides, so a check that
  !> printed at several k would show one of them and hide the rest -- which is how this
  !> routine started life, reporting three k and displaying one.
  SUBROUTINE melem_spin_channel_check(sz, ik, tol, l_report)
    COMPLEX, INTENT(IN) :: sz(:, :, :)        ! (nb,nb,2) the two channels' blocks
    INTEGER, INTENT(IN) :: ik
    REAL,    INTENT(IN) :: tol
    LOGICAL, INTENT(IN) :: l_report

    INTEGER :: nb, ch, m, n
    REAL    :: dev, offmax
    CHARACTER(LEN=200) :: msg

    nb = SIZE(sz, 1)
    IF (l_report) WRITE(oUnit, '(a,i0)') 'wannierlib collinear spin-channel check, k = ', ik
    DO ch = 1, SIZE(sz, 3)
      dev = 0.0; offmax = 0.0
      DO m = 1, nb
        dev = MAX(dev, ABS(ABS(sz(m, m, ch)) - 1.0))
        DO n = 1, nb
          IF (n /= m) offmax = MAX(offmax, ABS(sz(n, m, ch)))
        END DO
      END DO
      IF (l_report) WRITE(oUnit, '(a,i0,a,f10.6,a,f10.6)') &
        '  channel ', ch, ': max | |diag|-1 | = ', dev, '   largest off-diagonal = ', offmax
      IF (dev > tol .OR. offmax > tol) THEN
        WRITE(msg, '(a,i0,a,i0,a,es10.3,a,es10.3)') &
          'the spin block of collinear channel ', ch, ' at k = ', ik, &
          ' is not the unit matrix: max | |diag|-1 | = ', dev, ', largest off-diagonal = ', offmax
        CALL juDFT_error(TRIM(msg), &
          hint='sigma_z within one channel is that channel own overlap, so the block has '// &
               'to be the identity; a region left out of the integral moves it off', &
          calledby='melem_spin_channel_check')
      END IF
    END DO
  END SUBROUTINE melem_spin_channel_check

END MODULE m_melem_spin
