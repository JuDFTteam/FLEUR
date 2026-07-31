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
!>  spin up=1/dn=2) are the sum of two pieces, each block-selected by spin:
!>    * interstitial : reuse wannierlib_mmnkb_int at b=0 (global frame, no rot.);
!>    * muffin-tin   : melem_spin_mt_block below, mirrored on the library
!>                     routine wannierlib_mmk0_sph (abc%cof + radfun%integral, so
!>                     u/udot/LO are all included), extended to cross spin. The
!>                     local<->global spin rotation is already applied by the
!>                     modern calc_abc, so no ccchi is re-applied here.
MODULE m_melem_spin
  USE m_juDFT
  USE m_constants, ONLY : ImagUnit, oUnit
  USE m_types_atoms
  USE m_types_abc
  USE m_types_radfun
  USE m_types_nococonv
  USE m_types_stars
  USE m_types_lapw
  USE m_types_mat
  USE m_wannierlib_mmkb_int, ONLY : wannierlib_mmnkb_int
  IMPLICIT NONE
  PRIVATE
  PUBLIC :: melem_spin_bloch, melem_pauli_from_blocks, melem_spin_mt_block, &
            melem_spin_sumrule, melem_spin_peratom
CONTAINS

  !> Bloch-basis spin matrices at one k, MT + interstitial. Returns s0(nb,nb,3)
  !> and prints the sum-rule check. abc(:,:) = (2 spin, ntype); zMat holds the
  !> full two-component spinor (spin-down block at row offset nv(1)+nlotot).
  SUBROUTINE melem_spin_bloch(atoms, abc, radfun, nococonv, stars, lapw, zMat, num_bands, ik, s0, l_check)
    TYPE(t_atoms),     INTENT(IN)  :: atoms
    TYPE(t_abc),       INTENT(IN)  :: abc(:, :)          ! (2, ntype)
    TYPE(t_radfun),    INTENT(IN)  :: radfun(:)          ! (ntype)
    TYPE(t_nococonv),  INTENT(IN)  :: nococonv           ! %alph(:), %beta(:) per type
    TYPE(t_stars),     INTENT(IN)  :: stars
    TYPE(t_lapw),      INTENT(IN)  :: lapw
    TYPE(t_mat),       INTENT(IN)  :: zMat               ! full spinor eigenvectors at k
    INTEGER,           INTENT(IN)  :: num_bands, ik
    COMPLEX,           INTENT(OUT) :: s0(:, :, :)        ! (num_bands, num_bands, 3)
    LOGICAL,           INTENT(IN)  :: l_check            ! print the spin sum-rule for this k

    COMPLEX, ALLOCATABLE :: oi_uu(:, :), oi_dd(:, :), oi_ud(:, :), oi_du(:, :)  ! interstitial blocks (global)
    COMPLEX, ALLOCATABLE :: c_uu(:, :), c_dd(:, :), c_ud(:, :), c_du(:, :)      ! total-trace copies for sum-rule
    COMPLEX, ALLOCATABLE :: spa(:, :, :, :)                                      ! per-atom MT spin (global)
    INTEGER :: nb, io_dn, gb(3), na

    nb    = num_bands
    io_dn = lapw%nv(1) + atoms%nlotot      ! row offset of the spin-down block in zMat
    gb    = 0                              ! b = 0 (same k, on-site)

    ALLOCATE(oi_uu(nb, nb), oi_dd(nb, nb), oi_ud(nb, nb), oi_du(nb, nb))
    oi_uu = CMPLX(0.0, 0.0); oi_dd = CMPLX(0.0, 0.0)
    oi_ud = CMPLX(0.0, 0.0); oi_du = CMPLX(0.0, 0.0)

    ! ---- interstitial part only (global frame, b=0), block-selected by spin offset ----
    !   o_ab: bra spin a -> ioff, ket spin b -> ioff_b ;  up offset 0, dn offset io_dn
    CALL wannierlib_mmnkb_int(stars, lapw, lapw, 1, 1, zMat, zMat, gb, oi_uu, ioff=0,     ioff_b=0)
    CALL wannierlib_mmnkb_int(stars, lapw, lapw, 1, 1, zMat, zMat, gb, oi_dd, ioff=io_dn, ioff_b=io_dn)
    CALL wannierlib_mmnkb_int(stars, lapw, lapw, 1, 1, zMat, zMat, gb, oi_ud, ioff=0,     ioff_b=io_dn)
    CALL wannierlib_mmnkb_int(stars, lapw, lapw, 1, 1, zMat, zMat, gb, oi_du, ioff=io_dn, ioff_b=0)

    ! ---- total spin (GLOBAL) = interstitial Pauli + sum_atom (MT per-atom, rotated to global) ----
    !   the MT part must be summed per-atom AFTER the local->global rotation, otherwise an
    !   AFM (beta_2=pi) would add two local-frame "+" moments and the total would not vanish.
    CALL melem_pauli_from_blocks(oi_uu, oi_dd, oi_ud, oi_du, s0)
    ALLOCATE(spa(nb, nb, 3, atoms%nat))
    CALL melem_spin_peratom(atoms, abc, radfun, nococonv, spa)
    DO na = 1, atoms%nat
      s0(:, :, :) = s0(:, :, :) + spa(:, :, :, na)
    END DO

    ! ---- sum rule: the norm o_uu+o_dd is a (frame-invariant) trace; rebuild the total
    !      diagonal from interstitial + MT-local just for the diagnostic print ----
    IF (l_check) THEN
      ALLOCATE(c_uu(nb, nb), c_dd(nb, nb), c_ud(nb, nb), c_du(nb, nb))
      c_uu = oi_uu; c_dd = oi_dd; c_ud = oi_ud; c_du = oi_du
      CALL melem_spin_mt_block(atoms, abc, radfun, c_uu, c_dd, c_ud, c_du)
      CALL melem_spin_sumrule(s0, c_uu, c_dd, ik, tol=1.0e-3)
      DEALLOCATE(c_uu, c_dd, c_ud, c_du)
    END IF

    DEALLOCATE(oi_uu, oi_dd, oi_ud, oi_du, spa)
  END SUBROUTINE melem_spin_bloch

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

    ! Radial spin slots. radfun%integral is allocated (.,.,.,jspins,jspins), so with a
    ! single radial set (jspins=1, e.g. l_soc=T/l_noco=F) all four spin blocks share
    ! slot 1; indexing 2 there ran past the array. Bound read from the array itself.
    INTEGER :: js1, js2
    js1 = 1; js2 = MERGE(1, 2, SIZE(radfun(1)%integral, 4) < 2)
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

  !> Per-atom (site-resolved) muffin-tin Pauli spin: spa(nb,nb,3,nat), rotated to
  !> the GLOBAL spin frame. Same MT contraction as melem_spin_mt_block but
  !> KEEPING the global atom index na instead of summing. The interstitial is NOT
  !> site-resolvable and is excluded here, so spa is the muffin-tin site spin (the
  !> physical local moment; in an AFM this is +M on one sublattice and -M on the
  !> other).
  !>
  !> The abc%cof are LOCAL-frame spinor coefficients (calc_abc combines the two
  !> spinor blocks in each atom's local quantization axis). For a collinear FM_z
  !> (beta=0) local==global, but for a noncollinear/AFM texture (beta_a /= 0) each
  !> atom's (sx,sy,sz) must be rotated local->global by R_z(alpha_a) R_y(beta_a)
  !> using the noco angles, or the AFM sublattice comes out with the wrong sign
  !> (both moments look "+"). We apply that rotation per atom here.
  SUBROUTINE melem_spin_peratom(atoms, abc, radfun, nococonv, spa)
    TYPE(t_atoms),    INTENT(IN)  :: atoms
    TYPE(t_abc),      INTENT(IN)  :: abc(:, :)        ! (2 spin, ntype)
    TYPE(t_radfun),   INTENT(IN)  :: radfun(:)
    TYPE(t_nococonv), INTENT(IN)  :: nococonv         ! %alph(:), %beta(:) per type
    COMPLEX,          INTENT(OUT) :: spa(:, :, :, :)  ! (nb,nb,3,nat): 1=sx 2=sy 3=sz per atom, GLOBAL frame

    INTEGER :: nb, i, j, ntyp, iat, na, l, ll1, mm, lm, n_r, n_r2
    COMPLEX :: loc(2, 2), cx, cy, cz
    REAL    :: ca, sa, cb, sb

    ! Radial spin slots. radfun%integral is allocated (.,.,.,jspins,jspins), so with a
    ! single radial set (jspins=1, e.g. l_soc=T/l_noco=F) all four spin blocks share
    ! slot 1; indexing 2 there ran past the array. Bound read from the array itself.
    INTEGER :: js1, js2
    js1 = 1; js2 = MERGE(1, 2, SIZE(radfun(1)%integral, 4) < 2)
    nb = SIZE(spa, 1)
    spa = CMPLX(0.0, 0.0)
    DO j = 1, nb
      DO i = 1, nb
        na = 0
        DO ntyp = 1, atoms%ntype
          DO iat = 1, atoms%neq(ntyp)
            na = na + 1
            loc = CMPLX(0.0, 0.0)
            ca = COS(nococonv%alph(ntyp)); sa = SIN(nococonv%alph(ntyp))
            cb = COS(nococonv%beta(ntyp)); sb = SIN(nococonv%beta(ntyp))
            DO l = 0, atoms%lmax(ntyp)
              ll1 = l*(l + 1)
              DO mm = -l, l
                lm = ll1 + mm
                DO n_r = 1, abc(1, ntyp)%n_r(l)
                  DO n_r2 = 1, abc(1, ntyp)%n_r(l)
                    loc(1,1) = loc(1,1) + abc(1, ntyp)%cof(i,lm,n_r,iat)*CONJG(abc(1, ntyp)%cof(j,lm,n_r2,iat))*radfun(ntyp)%integral(n_r,n_r2,l,1,1)
                    loc(2,2) = loc(2,2) + abc(2, ntyp)%cof(i,lm,n_r,iat)*CONJG(abc(2, ntyp)%cof(j,lm,n_r2,iat))*radfun(ntyp)%integral(n_r,n_r2,l,js2,js2)
                    loc(1,2) = loc(1,2) + abc(1, ntyp)%cof(i,lm,n_r,iat)*CONJG(abc(2, ntyp)%cof(j,lm,n_r2,iat))*radfun(ntyp)%integral(n_r,n_r2,l,js1,js2)
                    loc(2,1) = loc(2,1) + abc(2, ntyp)%cof(i,lm,n_r,iat)*CONJG(abc(1, ntyp)%cof(j,lm,n_r2,iat))*radfun(ntyp)%integral(n_r,n_r2,l,js2,js1)
                  END DO
                END DO
              END DO
            END DO
            ! local-frame Pauli components at this site ...
            cx = loc(1,2) + loc(2,1)                ! sigma_x (local)
            cy = -ImagUnit * (loc(1,2) - loc(2,1))  ! sigma_y (local)
            cz = loc(1,1) - loc(2,2)                ! sigma_z (local)
            ! ... rotated local->global by R_z(alpha) R_y(beta) (noco convention)
            spa(i,j,1,na) =  ca*cb*cx - sa*cy + ca*sb*cz
            spa(i,j,2,na) =  sa*cb*cx + ca*cy + sa*sb*cz
            spa(i,j,3,na) = -sb*cx           + cb*cz
          END DO
        END DO
      END DO
    END DO
  END SUBROUTINE melem_spin_peratom

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
