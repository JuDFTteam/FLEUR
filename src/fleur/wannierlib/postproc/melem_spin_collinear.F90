!--------------------------------------------------------------------------------
! Copyright (c) 2026 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------
!>  The combined 2N spin operator of a collinear jspins=2 calculation, in real space.
!>
!>  Its ingredient is the cross-spin overlap <up|dn> on the coarse mesh, which is a Bloch
!>  quantity and is produced with the other coarse matrices. What is left here is the part
!>  that cannot run until BOTH spin channels have been wannierised: rotating that overlap
!>  with both gauges, assembling the 2N Pauli matrices, and exporting them.
!>
!>  sigma_z is +/-1 on the diagonal by the orthonormality of each channel, so it is written
!>  down rather than computed, and the transverse components follow from the single rotated
!>  block.
MODULE m_melem_spin_collinear
   USE m_juDFT
   USE m_constants, ONLY: ImagUnit, oUnit
   USE m_types_cell
   USE m_types_kpts
   USE m_types_mpi
   USE m_melem_ft, ONLY: melem_ft_to_real_reduce
   IMPLICIT NONE
   PRIVATE

   PUBLIC :: melem_rspauli_collinear

CONTAINS

   SUBROUTINE melem_rspauli_collinear(num_wann, x0, v_ch, cell, kpts, distk, fmpi)
      INTEGER, INTENT(IN) :: num_wann
      !> The cross-spin overlap <up|dn> on this rank's k-slice, in the Bloch basis:
      !> (num_bands, num_bands, nk_loc), in ascending global-k order.
      COMPLEX, INTENT(IN) :: x0(:, :, :)
      !> The Wannier gauge V = u_opt.u_matrix of each spin channel over the whole mesh:
      !> (num_bands, num_wann, nkptf, channel). Both channels are needed at once, which is
      !> why this step waits until neither wannierization is still running.
      COMPLEX, INTENT(IN) :: v_ch(:, :, :, :)
      TYPE(t_cell), INTENT(IN) :: cell
      TYPE(t_kpts), INTENT(IN) :: kpts
      INTEGER, INTENT(IN) :: distk(:)
      TYPE(t_mpi), INTENT(IN) :: fmpi

      INTEGER :: nb, nw, n2, nkl, kl, gk, iu, irpt, i, j, kk, nrpts
      INTEGER, ALLOCATABLE :: gk_loc(:), irvec(:, :), ndegen(:)
      COMPLEX, ALLOCATABLE :: Xk(:, :), tmp(:, :)
      COMPLEX, ALLOCATABLE :: sig_loc(:, :, :, :), s1(:, :, :), sr(:, :, :, :)

      nb = SIZE(x0, 1); nw = num_wann; n2 = 2*nw

      !> The gauges index the same manifolds the overlap does, and a mismatch would be a
      !> silently wrong rotation rather than a failure.
      IF (SIZE(v_ch, 1) /= nb .OR. SIZE(v_ch, 2) /= nw .OR. SIZE(v_ch, 4) /= 2) &
         CALL juDFT_error("melem_rspauli_collinear: the gauges do not match the manifold", &
                          calledby="melem_rspauli_collinear")

      nkl = COUNT(distk == fmpi%irank); ALLOCATE (gk_loc(nkl)); j = 0
      DO i = 1, SIZE(distk)
         IF (distk(i) == fmpi%irank) THEN; j = j + 1; gk_loc(j) = i; END IF
      END DO
      IF (SIZE(x0, 3) < nkl) CALL juDFT_error( &
         "melem_rspauli_collinear: fewer overlap slices than k-points on this rank", &
         calledby="melem_rspauli_collinear")

      ALLOCATE (Xk(nw, nw), tmp(nb, nw))
      ALLOCATE (sig_loc(n2, n2, 3, MAX(1, nkl)), source=CMPLX(0.0, 0.0))

      DO kl = 1, nkl
         gk = gk_loc(kl)
         ! rotate to the WF gauge: X = V_up^dagger o_ud V_dn
         tmp = MATMUL(x0(:, :, kl), v_ch(:, :, gk, 2))
         Xk = MATMUL(CONJG(TRANSPOSE(v_ch(:, :, gk, 1))), tmp)
         ! assemble the 2N Pauli in the WF gauge (sigma_z block-diagonal +/- I by orthonormality)
         DO i = 1, nw
            sig_loc(i, i, 3, kl) = CMPLX(1.0, 0.0)
            sig_loc(nw + i, nw + i, 3, kl) = CMPLX(-1.0, 0.0)
         END DO
         sig_loc(1:nw, nw + 1:n2, 1, kl) = Xk
         sig_loc(nw + 1:n2, 1:nw, 1, kl) = CONJG(TRANSPOSE(Xk))
         sig_loc(1:nw, nw + 1:n2, 2, kl) = -ImagUnit*Xk
         sig_loc(nw + 1:n2, 1:nw, 2, kl) = ImagUnit*CONJG(TRANSPOSE(Xk))
      END DO
      DEALLOCATE (Xk, tmp)

      ! FT-reduce each of the 3 components (collective), rank 0 writes rspauli.1
      DO kk = 1, 3
         CALL melem_ft_to_real_reduce(cell, kpts, sig_loc(:, :, kk, :), gk_loc, fmpi%mpi_comm, s1, irvec, ndegen, nrpts)
         IF (kk == 1) ALLOCATE (sr(n2, n2, nrpts, 3))
         sr(:, :, :, kk) = s1; DEALLOCATE (s1)
      END DO
      IF (fmpi%irank == 0) THEN
         OPEN (newunit=iu, file='rspauli.1', status='replace')
         DO irpt = 1, nrpts
            DO j = 1, n2
               DO i = 1, n2
                  DO kk = 1, 3
                     WRITE (iu, '(i3,1x,i3,1x,i3,1x,i3,1x,i3,1x,i3,1x,f20.8,1x,f20.8)') &
                        irvec(1, irpt), irvec(2, irpt), irvec(3, irpt), i, j, kk, REAL(sr(i, j, irpt, kk)), AIMAG(sr(i, j, irpt, kk))
                  END DO
               END DO
            END DO
         END DO
         CLOSE (iu)
         WRITE (oUnit, '(a,i0,a)') 'wannierlib: wrote rspauli.1 (combined 2N collinear spin, ', nrpts, ' R-vectors, distributed FT)'
      END IF
      DEALLOCATE (sig_loc, gk_loc)
      IF (ALLOCATED(sr)) DEALLOCATE (sr)
      IF (ALLOCATED(irvec)) DEALLOCATE (irvec, ndegen)
   END SUBROUTINE melem_rspauli_collinear

END MODULE m_melem_spin_collinear
