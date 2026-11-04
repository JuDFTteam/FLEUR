!--------------------------------------------------------------------------------
! Copyright (c) 2026 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------

!> Hands the gauge Wannier90 returned to `results`, completed over the whole mesh.
!>
!> THE TWO FACTORS ARE KEPT APART, not multiplied, because they answer different questions:
!> u_opt fixes WHICH subspace was selected, u_mlwf fixes HOW it was rotated inside that
!> subspace. The total gauge is the product, in this order:
!>
!>       U(k) = u_opt(k) . u_mlwf(k)          (num_bands x num_wann)
!>
!> Multiplying them is one line for the reader; separating them again afterwards is not.
MODULE m_wannierlib_gauge_results
   USE m_juDFT
   USE m_constants, ONLY: oUnit
   USE m_types_kpts
   USE m_types_input
   USE m_types_mpi
#ifdef CPP_MPI
   USE mpi
#endif
   USE m_types_misc
   IMPLICIT NONE
   PRIVATE
   PUBLIC :: wannierlib_store_gauge
CONTAINS

   !> Copy the Wannier gauge into results, completed over the whole mesh.
   !>
   !> Nothing downstream of here is in the testset: the only consumer is
   !> interpolate_bandstructure, on the electron-phonon path in dfpt/superconductivity, and
   !> both of its call sites ask it not to write, so there is no file to compare. What the
   !> routine reports instead is the invariant that consumer depends on, so the mesh is
   !> checked where it is handed over rather than where it is used.
   !>
   !> The electron-phonon interpolation reads results%U_mat and results%U_dis; what run_w90
   !> returns is each rank's own k-points with the rest left at zero, which is what the
   !> operator pass expects, so the sum that completes the mesh is done on a copy and the
   !> arrays the caller holds are left exactly as they were.
   SUBROUTINE wannierlib_store_gauge(fmpi, distk, kpts, input, num_bands, num_wann, jspin, &
                                     u_matrix, u_opt, results)
      TYPE(t_mpi), INTENT(IN) :: fmpi
      INTEGER, INTENT(IN) :: distk(:)
      TYPE(t_kpts), INTENT(IN) :: kpts
      TYPE(t_input), INTENT(IN) :: input
      INTEGER, INTENT(IN) :: num_bands, num_wann, jspin
      COMPLEX, INTENT(IN) :: u_matrix(:, :, :)      ! (nw, nw, nk), this rank's k-points
      COMPLEX, INTENT(IN) :: u_opt(:, :, :)         ! (nb, nw, nk), this rank's k-points
      TYPE(t_results), INTENT(INOUT) :: results

      COMPLEX, ALLOCATABLE :: buf(:, :, :)
      COMPLEX, ALLOCATABLE :: vk(:, :), ovl(:, :)
      INTEGER :: ikpt, ierr, iw, nrow
      REAL :: dev

      IF (SIZE(u_matrix, 3) == kpts%nkptf) THEN
         buf = u_matrix
         DO ikpt = 1, kpts%nkptf
            IF (distk(ikpt) /= fmpi%irank) buf(:, :, ikpt) = CMPLX(0.0, 0.0)
         END DO
#ifdef CPP_MPI
         IF (fmpi%isize > 1) &
            CALL MPI_ALLREDUCE(MPI_IN_PLACE, buf, SIZE(buf), MPI_DOUBLE_COMPLEX, MPI_SUM, &
                               fmpi%mpi_comm, ierr)
#endif
         IF (.NOT. ALLOCATED(results%U_mat)) THEN
            ALLOCATE (results%U_mat(num_wann, num_wann, kpts%nkptf, input%jspins), &
                      source=CMPLX(0.0, 0.0))
         END IF
         results%U_mat(:, :, :, jspin) = buf
         DEALLOCATE (buf)
      END IF

      !> Only meaningful when the subspace was chosen, i.e. more bands than functions;
      !> with num_bands == num_wann u_opt is the projection itself and carries nothing.
      !> That is the condition the routine this replaces used, expressed as a test because
      !> here u_opt is always allocated.
      IF (num_bands > num_wann .AND. SIZE(u_opt, 3) == kpts%nkptf) THEN
         buf = u_opt
         DO ikpt = 1, kpts%nkptf
            IF (distk(ikpt) /= fmpi%irank) buf(:, :, ikpt) = CMPLX(0.0, 0.0)
         END DO
#ifdef CPP_MPI
         IF (fmpi%isize > 1) &
            CALL MPI_ALLREDUCE(MPI_IN_PLACE, buf, SIZE(buf), MPI_DOUBLE_COMPLEX, MPI_SUM, &
                               fmpi%mpi_comm, ierr)
#endif
         IF (.NOT. ALLOCATED(results%U_dis)) THEN
            ALLOCATE (results%U_dis(num_bands, num_wann, kpts%nkptf, input%jspins), &
                      source=CMPLX(0.0, 0.0))
         END IF
         results%U_dis(:, :, :, jspin) = buf
         DEALLOCATE (buf)
      END IF

      !> The one property a consumer can rely on, measured on the arrays it will read and
      !> not on the working copy above, which would only restate the sum just performed.
      !>
      !> V(k) = u_opt(k) . u_mlwf(k) has orthonormal columns, so V^dagger V is the num_wann
      !> identity at every k of the mesh. A k-point the reduction failed to collect is left
      !> at zero and gives the null matrix instead: a deviation of 1, ten orders of
      !> magnitude above the rounding a product of this size accumulates. That is the
      !> failure this guards against, because each rank is handed its own k-points and the
      !> rest as zeros, so a wrong distk or a missing rank is silent in every other way.
      IF (fmpi%irank == 0 .AND. ALLOCATED(results%U_mat)) THEN
         nrow = MERGE(num_bands, num_wann, num_bands > num_wann)
         ALLOCATE (vk(nrow, num_wann), ovl(num_wann, num_wann))
         dev = 0.0
         DO ikpt = 1, kpts%nkptf
            IF (num_bands > num_wann) THEN
               vk = MATMUL(results%U_dis(:, :, ikpt, jspin), results%U_mat(:, :, ikpt, jspin))
            ELSE
               vk = results%U_mat(:, :, ikpt, jspin)
            END IF
            ovl = MATMUL(TRANSPOSE(CONJG(vk)), vk)
            DO iw = 1, num_wann
               ovl(iw, iw) = ovl(iw, iw) - CMPLX(1.0, 0.0)
            END DO
            dev = MAX(dev, MAXVAL(ABS(ovl)))
         END DO
         DEALLOCATE (vk, ovl)
         WRITE (oUnit, '(a,i0,a,es12.5)') 'wannierlib gauge: spin ', jspin, &
            ', max |V^H V - 1| over the mesh = ', dev
      END IF
   END SUBROUTINE wannierlib_store_gauge
END MODULE m_wannierlib_gauge_results
