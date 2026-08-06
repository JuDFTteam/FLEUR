!--------------------------------------------------------------------------------
! Copyright (c) 2026 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------
!>  Provider of the spin (Pauli) operator matrix elements in the Bloch basis,
!>  on-site (b = 0, same k):
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
!>  spin up=1/dn=2) take their muffin-tin contribution from melem_spin_mt_block
!>  below: abc%cof contracted with radfun%integral, so u, udot and the local
!>  orbitals are all included, extended to cross spin. The local-to-global spin
!>  rotation is already applied when the abc coefficients are built, so no
!>  further rotation is applied here.
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
  PUBLIC :: melem_pauli_from_blocks, melem_spin_mt_block, melem_spin_sumrule
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

  !> Muffin-tin contribution to the four spin-block overlaps
  !>   o_ss'(m,n) += <phi^s_m | phi^s'_n>_MT ,   s,s' in {1,2}.
  !>
  !> Mirrors wannierlib_mmk0_sph (same abc%cof + radfun%integral contraction, so
  !> u / udot / local-orbital channels are all summed via the n_r index), crossing
  !> the two spinor components abc(:,1) and abc(:,2).
  !>
  !> NO explicit spin rotation here: the modern library calc_abc ALREADY applies
  !> the SU(2) local<->global rotation (ccchi = conjg(nococonv%umat)) and combines
  !> both spinor blocks, so abc%cof are the (local-frame) spin components. The old
  !> fleur-8.1 wann_mmk0_updown_sph used raw abcof + an explicit ccchi; that path is
  !> outdated. For a single global axis (Fe-FM_z, beta=0) local==global. A true
  !> noncollinear texture (Mn3Ir) would need a global-frame rotation downstream, but
  !> the spin sum rule (norm, |<sigma>|) is rotation-invariant, so this suffices here.
  SUBROUTINE melem_spin_mt_block(atoms, abc, radfun, o_uu, o_dd, o_ud, o_du)
    TYPE(t_atoms),  INTENT(IN) :: atoms
    TYPE(t_abc),    INTENT(IN) :: abc(:, :)        ! (2 spin, ntype) local-frame coeffs
    TYPE(t_radfun), INTENT(IN) :: radfun(:)        ! (ntype) : %integral(n_r,n_r2,l,1,1)
    COMPLEX, INTENT(INOUT) :: o_uu(:, :), o_dd(:, :), o_ud(:, :), o_du(:, :)   ! (nb,nb)

    INTEGER :: nb, i, j, ntyp, iat, l, ll1, mm, lm, n_r, n_r2
    COMPLEX :: loc(2, 2)

    ! The four spin blocks are indexed with the slots of their own components.
    INTEGER :: js1, js2
    !> The spin index comes first. An assumed-shape dummy accepts the transpose
    !> without complaint and reinterprets the two indices on the same memory, which
    !> is right only when there is a single atom type.
    IF (SIZE(abc, 1) /= 2 .OR. SIZE(abc, 2) /= atoms%ntype) &
      CALL judft_bug("melem_spin_mt_block: the coefficients must have shape (2,ntype)")

    js1 = radial_slot(radfun, 1); js2 = radial_slot(radfun, 2)
    nb = SIZE(o_uu, 1)
    DO j = 1, nb                       ! ket band
      DO i = 1, nb                     ! bra band
        loc = CMPLX(0.0, 0.0)
        DO ntyp = 1, atoms%ntype
          DO l = 0, atoms%lmax(ntyp)
            ll1 = l*(l + 1)
            DO mm = -l, l
              lm = ll1 + mm
              DO iat = 1, atoms%neq(ntyp)
                ! radial overlap carries the spin indices (jspins=2: up/down radials differ):
                !   loc(s1,s2) uses integral(n_r,n_r2,l,s1,s2)
                DO n_r = 1, abc(1, ntyp)%n_r(l)
                  DO n_r2 = 1, abc(1, ntyp)%n_r(l)
                    loc(1, 1) = loc(1, 1) + abc(1, ntyp)%cof(i,lm,n_r,iat)*CONJG(abc(1, ntyp)%cof(j,lm,n_r2,iat))*radfun(ntyp)%integral(n_r,n_r2,l,1,1)
                    loc(2, 2) = loc(2, 2) + abc(2, ntyp)%cof(i,lm,n_r,iat)*CONJG(abc(2, ntyp)%cof(j,lm,n_r2,iat))*radfun(ntyp)%integral(n_r,n_r2,l,js2,js2)
                    loc(1, 2) = loc(1, 2) + abc(1, ntyp)%cof(i,lm,n_r,iat)*CONJG(abc(2, ntyp)%cof(j,lm,n_r2,iat))*radfun(ntyp)%integral(n_r,n_r2,l,js1,js2)
                    loc(2, 1) = loc(2, 1) + abc(2, ntyp)%cof(i,lm,n_r,iat)*CONJG(abc(1, ntyp)%cof(j,lm,n_r2,iat))*radfun(ntyp)%integral(n_r,n_r2,l,js2,js1)
                  END DO
                END DO
              END DO
            END DO
          END DO
        END DO
        o_uu(i, j) = o_uu(i, j) + loc(1, 1)
        o_dd(i, j) = o_dd(i, j) + loc(2, 2)
        o_ud(i, j) = o_ud(i, j) + loc(1, 2)
        o_du(i, j) = o_du(i, j) + loc(2, 1)
      END DO
    END DO
  END SUBROUTINE melem_spin_mt_block

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
