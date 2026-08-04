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
   USE m_judft
   IMPLICIT NONE
   PRIVATE
   PUBLIC :: melem_overlap_interstitial

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

END MODULE m_melem_overlap
