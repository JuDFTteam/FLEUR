!--------------------------------------------------------------------------------
! Copyright (c) 2026 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------
!>  The Wannier-gauge Hamiltonian on the coarse mesh, H_W(k).
!>
!>  Every interpolation driver needs it before it can interpolate anything -- the band
!>  driver to get the bands, the operator drivers to get the eigenvectors they project on,
!>  the velocity to differentiate. It is built here once, so a correction has a single
!>  place to reach.
!>
!>      eigval2(m,k) = sum_i eig_i(k) |u_opt(i,m,k)|^2        (the outer window, if any)
!>      H_W(i,j,k)   = sum_m eigval2(m,k) conjg(U(m,i,k)) U(m,j,k)
!>
!>  Without disentanglement u_opt is the identity and eigval2 is the input spectrum, so
!>  H_W is a unitary rotation of it and its eigenvalues ARE the ab-initio ones. That is
!>  what makes the interpolation exact on the mesh it was built from, which is what the
!>  WannFeBccInterp test asserts.
MODULE m_melem_hamk
   USE m_types_melem_manifold, ONLY: t_melem_manifold
   IMPLICIT NONE
   PRIVATE

   PUBLIC :: melem_build_hamk

CONTAINS

   SUBROUTINE melem_build_hamk(this, eig, u_matrix, u_opt, ham_k)
      TYPE(t_melem_manifold), INTENT(IN) :: this
      REAL,    INTENT(IN) :: eig(:, :)                       !< (num_bands, nk)
      COMPLEX, INTENT(IN) :: u_matrix(:, :, :)               !< (nw,nw,nk)  MLWF gauge
      COMPLEX, INTENT(IN) :: u_opt(:, :, :)                  !< (nb,nw,nk)  disentangled
      COMPLEX, ALLOCATABLE, INTENT(OUT) :: ham_k(:, :, :)    !< (nw,nw,nk)

      INTEGER :: nw, nb, nk, k, i, j, m, counter
      LOGICAL :: have_dis
      REAL, ALLOCATABLE :: eigval2(:, :), eigval_opt(:)

      nw = this%num_wann; nb = this%num_bands; nk = SIZE(u_matrix, 3)
      have_dis = (nb > nw)
      ALLOCATE (eigval2(nw, nk), source=0.0)
      IF (have_dis) THEN
         ALLOCATE (eigval_opt(nb))
         DO k = 1, nk
            counter = 0; eigval_opt = 0.0
            DO j = 1, nb
               IF (eig(j, k) >= this%dis_win_min .AND. eig(j, k) <= this%dis_win_max) THEN
                  counter = counter + 1; eigval_opt(counter) = eig(j, k)
               END IF
            END DO
            DO m = 1, nw
               DO i = 1, counter
                  eigval2(m, k) = eigval2(m, k) + eigval_opt(i) * ABS(u_opt(i, m, k))**2
               END DO
            END DO
         END DO
         DEALLOCATE (eigval_opt)
      ELSE
         eigval2(1:nw, :) = eig(1:nw, :)
      END IF
      ALLOCATE (ham_k(nw, nw, nk), source=CMPLX(0.0, 0.0))
      DO k = 1, nk
         DO j = 1, nw
            DO i = 1, nw
               DO m = 1, nw
                  ham_k(i, j, k) = ham_k(i, j, k) + eigval2(m, k) * CONJG(u_matrix(m, i, k)) * u_matrix(m, j, k)
               END DO
            END DO
         END DO
      END DO
      DEALLOCATE (eigval2)
   END SUBROUTINE melem_build_hamk

END MODULE m_melem_hamk
