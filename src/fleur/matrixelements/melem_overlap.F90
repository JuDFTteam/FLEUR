!--------------------------------------------------------------------------------
! Copyright (c) 2026 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------
!>  Overlap of two sets of eigenvector coefficients over the interstitial at one
!>  k-point. The spin operators need it: a Pauli matrix is the identity in real
!>  space, so the interstitial part of <psi_m|sigma|psi_n> is the plain overlap of
!>  the spin components, and the four spin blocks are four such overlaps.
MODULE m_melem_overlap
   USE m_types_stars
   USE m_types_lapw
   USE m_types_mat
   USE m_types_atoms
   USE m_types_abc
   USE m_types_cell
   USE m_types_radfun
   USE m_melem_mmkb_int
   USE m_melem_mmkb_sph
   USE m_melem_mmkb_vac, ONLY: melem_mmkb_vac
   USE m_types_melem_vacabc, ONLY: t_melem_vacabc
   USE m_melem_ujugaunt, ONLY: melem_ujugaunt
   USE m_melem_check, ONLY: MELEM_CHECK_TOL
   USE m_constants, ONLY: oUnit
   USE m_judft
   IMPLICIT NONE
   PRIVATE
   PUBLIC :: melem_overlap_interstitial, melem_overlap_states
   PUBLIC :: melem_overlap_check_identity

CONTAINS

   !> ovl(m,n) = <psi^a_m | psi^b_n> over the interstitial, at one k for both sides.
   !>
   !> The two sides may come from one matrix or from two: ioff_a and ioff_b give the
   !> first row of each, which is how the components of a stacked spinor are picked
   !> apart, and lapw_a / lapw_b may differ, which is how the two channels of a
   !> collinear calculation are combined.
   !>
   !> The overlap of two plane waves over the interstitial is a Fourier component of
   !> the step function, Theta(G_n - G_m): ig maps that difference to its star and
   !> rgphs carries the phase.
   !>
   !> Real eigenvectors get their own branch. There the step-function factor is used
   !> without its phase and the intermediate is not conjugated, since both operations
   !> are no-ops on real data, and the result is real: it is returned in the real part.
   SUBROUTINE melem_overlap_interstitial(stars, lapw_a, lapw_b, zmat_a, zmat_b, &
                                         ioff_a, ioff_b, ovl)
      TYPE(t_stars), INTENT(IN)  :: stars
      TYPE(t_lapw),  INTENT(IN)  :: lapw_a, lapw_b
      TYPE(t_mat),   INTENT(IN)  :: zmat_a, zmat_b
      INTEGER,       INTENT(IN)  :: ioff_a, ioff_b
      COMPLEX,       INTENT(OUT) :: ovl(:,:)

      COMPLEX, ALLOCATABLE :: stepf(:,:), zstep(:,:)
      REAL,    ALLOCATABLE :: stepf_r(:,:), zstep_r(:,:), ovl_r(:,:)
      COMPLEX :: theta
      INTEGER :: nv_a, nv_b, nb, i, j, i1, i2, i3, in

      IF (zmat_a%l_real .NEQV. zmat_b%l_real) &
         CALL judft_bug("melem_overlap_interstitial: one side is real and the other complex")

      nv_a = lapw_a%nv(1)
      nv_b = lapw_b%nv(1)
      nb   = SIZE(ovl, 1)

      IF (zmat_a%l_real) THEN
         ALLOCATE(stepf_r(nv_b, nv_a), source=0.0)
         ALLOCATE(zstep_r(nv_a, nb), ovl_r(nb, nb))
      ELSE
         ALLOCATE(stepf(nv_b, nv_a), source=CMPLX(0.0, 0.0))
         ALLOCATE(zstep(nv_a, nb))
      END IF

      DO i = 1, nv_a
         DO j = 1, nv_b
            i1 = lapw_b%k1(j,1) - lapw_a%k1(i,1)
            i2 = lapw_b%k2(j,1) - lapw_a%k2(i,1)
            i3 = lapw_b%k3(j,1) - lapw_a%k3(i,1)
            in = stars%ig(i1, i2, i3)
            IF (in == 0) CYCLE
            theta = stars%rgphs(i1,i2,i3) * stars%ustep(in)
            IF (zmat_a%l_real) THEN
               stepf_r(j,i) = REAL(theta)
            ELSE
               stepf(j,i) = CONJG(theta)
            END IF
         END DO
      END DO

      IF (zmat_a%l_real) THEN
         CALL dgemm('T', 'N', nv_a, nb, nv_b, REAL(1.0), stepf_r, nv_b, &
                    zmat_b%data_r(1+ioff_b, 1), zmat_b%matsize1, REAL(0.0), zstep_r, nv_a)
         CALL dgemm('T', 'N', nb, nb, nv_a, REAL(1.0), &
                    zmat_a%data_r(1+ioff_a, 1), zmat_a%matsize1, &
                    zstep_r, nv_a, REAL(0.0), ovl_r, nb)
         ovl = CMPLX(ovl_r, 0.0)
      ELSE
         CALL zgemm('T', 'N', nv_a, nb, nv_b, CMPLX(1.0,0.0), stepf, nv_b, &
                    zmat_b%data_c(1+ioff_b, 1), zmat_b%matsize1, CMPLX(0.0,0.0), zstep, nv_a)
         zstep = CONJG(zstep)
         CALL zgemm('T', 'N', nb, nb, nv_a, CMPLX(1.0,0.0), &
                    zmat_a%data_c(1+ioff_a, 1), zmat_a%matsize1, &
                    zstep, nv_a, CMPLX(0.0,0.0), ovl, nb)
      END IF
   END SUBROUTINE melem_overlap_interstitial

   !> The overlap of two sets of states over the whole cell: muffin-tin part, interstitial
   !> part and, on a film, vacuum part, all into the same matrix.
   !>
   !> The vacuum expansions are optional and a bulk caller leaves them out. Handing them in
   !> unexpanded (nv2 == 0) is the same as leaving them out, so a caller that has them only
   !> sometimes needs one call site rather than two.
   !>
   !> The two sides are independent: each has its own k-point, its own basis and its own
   !> coefficients, and gb is the reciprocal lattice vector that brings the second one back
   !> to the first Brillouin zone. Nothing here requires them to be a k-point and one of its
   !> neighbours; a pair of neighbours of a third point works the same way, with gb the
   !> difference of the two vectors that folded them.
   !>
   !> Both halves accumulate, so ovl carries whatever the caller left in it: it is the
   !> caller that decides whether a slot starts at zero. And kdiff has to contain the
   !> difference the two sides span, or the muffin-tin half stops the run rather than
   !> guessing which table entry to use.
   SUBROUTINE melem_overlap_states(stars, atoms, lapw_a, lapw_b, zmat_a, zmat_b, &
                                   abc_a, abc_b, jspin_a, jspin_b, bkpt_a, bkpt_b, gb, &
                                   ujug, kdiff, nntot, ioff_a, ioff_b, ovl, vac_a, vac_b)
      TYPE(t_stars), INTENT(IN) :: stars
      TYPE(t_atoms), INTENT(IN) :: atoms
      TYPE(t_lapw),  INTENT(IN) :: lapw_a, lapw_b
      TYPE(t_mat),   INTENT(IN) :: zmat_a, zmat_b
      TYPE(t_abc),   INTENT(IN) :: abc_a(:), abc_b(:)
      INTEGER,       INTENT(IN) :: jspin_a, jspin_b
      REAL,          INTENT(IN) :: bkpt_a(3), bkpt_b(3)
      INTEGER,       INTENT(IN) :: gb(3)
      COMPLEX,       INTENT(IN) :: ujug(0:, 0:, :, :, :, :)
      REAL,          INTENT(IN) :: kdiff(:, :)
      INTEGER,       INTENT(IN) :: nntot, ioff_a, ioff_b
      COMPLEX,    INTENT(INOUT) :: ovl(:, :)
      TYPE(t_melem_vacabc), INTENT(IN), OPTIONAL :: vac_a, vac_b

      LOGICAL :: l_vac

      CALL melem_mmkb_int(stars, lapw_a, lapw_b, jspin_a, jspin_b, zmat_a, zmat_b, gb, ovl, &
                          ioff=ioff_a, ioff_b=ioff_b)
      CALL melem_mmkb_sph(atoms, abc_a, abc_b, bkpt_b, gb, bkpt_a, ujug, kdiff, nntot, ovl)

      !> Two statements because Fortran does not promise to stop at the first .AND.
      l_vac = .FALSE.
      IF (PRESENT(vac_a) .AND. PRESENT(vac_b)) l_vac = (vac_a%nv2 > 0 .AND. vac_b%nv2 > 0)
      IF (l_vac) CALL melem_mmkb_vac(vac_a, vac_b, gb, ovl)
   END SUBROUTINE melem_overlap_states

   !>  M_mn(k, k) has to be the identity. The eigenvectors are orthonormal over the whole
   !>  cell, so the regional halves of the overlap must add up to delta_mn -- and this is the
   !>  one invariant that exercises every region at once.
   !>
   !>  It is the only test that can tell a WRONG vacuum term from a missing one. On a film
   !>  the interpolated bands come out exact on the coarse mesh whatever the gauge, because
   !>  H_W is then a unitary rotation of the input spectrum, so they cannot see it. This can:
   !>  a region left out, or counted with the wrong measure, shows up here immediately.
   !>
   !>  It works in bulk too, with two halves instead of three, which is what makes it usable:
   !>  bulk is the control. If bulk is clean and a film is not, the vacuum term is at fault
   !>  rather than the check.
   !>
   !>  Warns rather than stops, like m_melem_check: a tolerance right for one basis is not
   !>  obviously right for the next.
   SUBROUTINE melem_overlap_check_identity(stars, atoms, cell, lapw, zmat, abc, radfun, &
                                           jspin_rad, ioff, nbnd, ik, vac, tol, l_ok)
      TYPE(t_stars), INTENT(IN) :: stars
      TYPE(t_atoms), INTENT(IN) :: atoms
      TYPE(t_cell),  INTENT(IN) :: cell
      TYPE(t_lapw),  INTENT(IN) :: lapw
      TYPE(t_mat),   INTENT(IN) :: zmat
      TYPE(t_abc),   INTENT(IN) :: abc(:)
      TYPE(t_radfun), INTENT(IN) :: radfun(:)
      INTEGER, INTENT(IN) :: jspin_rad, ioff, nbnd, ik
      TYPE(t_melem_vacabc), INTENT(IN), OPTIONAL :: vac
      REAL, INTENT(IN), OPTIONAL :: tol
      LOGICAL, INTENT(OUT), OPTIONAL :: l_ok

      REAL :: kdiff0(3, 1), t, dmax, omax, d
      COMPLEX, ALLOCATABLE :: ujug0(:, :, :, :, :, :), ovl(:, :)
      INTEGER :: i, j
      CHARACTER(LEN=12) :: reg

      t = MELEM_CHECK_TOL
      IF (PRESENT(tol)) t = tol
      reg = 'MT+int'
      !> Two statements: an unbuilt expansion may be passed in bulk, and Fortran does not
      !> promise to stop at the first .AND.
      IF (PRESENT(vac)) THEN
         IF (vac%nv2 > 0) reg = 'MT+int+vac'
      END IF

      !> b = 0 needs its own one-entry table: the muffin-tin half finds b in kdiff BY VALUE,
      !> and the neighbour table has no zero in it.
      kdiff0 = 0.0
      CALL melem_ujugaunt(atoms, cell, 1, kdiff0, radfun, radfun, jspin_rad, jspin_rad, &
                          .FALSE., 1, ujug0)

      ALLOCATE (ovl(nbnd, nbnd), source=CMPLX(0.0, 0.0))
      CALL melem_overlap_states(stars, atoms, lapw, lapw, zmat, zmat, abc, abc, &
                                jspin_rad, jspin_rad, lapw%bkpt, lapw%bkpt, [0, 0, 0], &
                                ujug0, kdiff0, 1, ioff, ioff, ovl, vac_a=vac, vac_b=vac)

      dmax = 0.0; omax = 0.0
      DO j = 1, nbnd
         DO i = 1, nbnd
            IF (i == j) THEN
               d = ABS(ovl(i, j) - CMPLX(1.0, 0.0)); dmax = MAX(dmax, d)
            ELSE
               d = ABS(ovl(i, j)); omax = MAX(omax, d)
            END IF
         END DO
      END DO

      IF (dmax > t .OR. omax > t) THEN
         WRITE (oUnit, '(a,i0,a,2(a,es12.4))') 'wannierlib overlap check [k=', ik, ', '// &
            TRIM(reg)//']: M(k,k) is not the identity', '  diagonal ', dmax, '  off-diagonal ', omax
      ELSE
         WRITE (oUnit, '(a,i0,a,2(a,es12.4))') 'wannierlib overlap check [k=', ik, ', '// &
            TRIM(reg)//']: M(k,k) = 1 ok', '  diagonal ', dmax, '  off-diagonal ', omax
      END IF
      IF (PRESENT(l_ok)) l_ok = (dmax <= t .AND. omax <= t)

      DEALLOCATE (ovl, ujug0)
   END SUBROUTINE melem_overlap_check_identity

END MODULE m_melem_overlap
