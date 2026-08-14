!--------------------------------------------------------------------------------
! Copyright (c) 2026 Peter Gruenberg Institut, Forschungszentrum Juelich, Germany
! This file is part of FLEUR and available as free software under the conditions
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------

MODULE m_local_d_spin_projector_core
   USE m_intgr, ONLY: intgr0
   USE m_juDFT, ONLY: juDFT_error
   IMPLICIT NONE
   PRIVATE

   INTEGER, PARAMETER, PUBLIC :: local_d_l = 2
   INTEGER, PARAMETER, PUBLIC :: local_d_n_orbitals = 2*local_d_l + 1
   INTEGER, PARAMETER, PUBLIC :: local_d_n_spins = 2

   PUBLIC :: local_d_spin_contract_density
   PUBLIC :: local_d_spin_ordered_radial_metric

CONTAINS

   SUBROUTINE local_d_spin_ordered_radial_metric(radial_values, r0, dx, n_grid, metric)
      ! Construct the ordered radial metric used by the local d x spin density
      ! matrix. radial_values has dimensions
      !   (radial mesh, large/small component, radial order, local spin).
      ! Radial order 1 is u_l, order 2 is udot_l, and subsequent orders are
      ! local orbitals of the same l. The spin components are in the native
      ! local muffin-tin spin frame.
      !
      ! All ordered (i,j,s,s') combinations are integrated independently.
      ! In particular, no i <-> j symmetry is imposed inside a cross-spin
      ! block: S(i,j,s,s') is related to S(j,i,s',s), not generally to
      ! S(j,i,s,s').
      REAL, INTENT(IN)  :: radial_values(:, :, :, :)
      REAL, INTENT(IN)  :: r0, dx
      INTEGER, INTENT(IN) :: n_grid
      REAL, INTENT(OUT) :: metric(:, :, :, :)

      INTEGER :: i, j, ispin, jspin, n_radial
      REAL, ALLOCATABLE :: integrand(:)

      n_radial = SIZE(radial_values, 3)
      IF (n_grid < 2 .OR. n_grid > SIZE(radial_values, 1)) THEN
         CALL juDFT_error("Invalid radial mesh size for local d-spin projector", &
                          calledby="m_local_d_spin_projector_core")
      END IF
      IF (SIZE(radial_values, 2) /= 2) THEN
         CALL juDFT_error("Local d-spin projector requires large and small radial components", &
                          calledby="m_local_d_spin_projector_core")
      END IF
      IF (SIZE(radial_values, 4) /= local_d_n_spins) THEN
         CALL juDFT_error("Local d-spin projector requires two local spin components", &
                          calledby="m_local_d_spin_projector_core")
      END IF
      IF (n_radial < 1) THEN
         CALL juDFT_error("Local d-spin projector requires at least one radial function", &
                          calledby="m_local_d_spin_projector_core")
      END IF
      IF (SIZE(metric, 1) /= n_radial .OR. SIZE(metric, 2) /= n_radial .OR. &
          SIZE(metric, 3) /= local_d_n_spins .OR. SIZE(metric, 4) /= local_d_n_spins) THEN
         CALL juDFT_error("Ordered radial metric has inconsistent dimensions", &
                          calledby="m_local_d_spin_projector_core")
      END IF

      ALLOCATE(integrand(n_grid))
      DO ispin = 1, local_d_n_spins
         DO jspin = 1, local_d_n_spins
            DO i = 1, n_radial
               DO j = 1, n_radial
                  integrand = radial_values(1:n_grid, 1, i, ispin)*radial_values(1:n_grid, 1, j, jspin) &
                            + radial_values(1:n_grid, 2, i, ispin)*radial_values(1:n_grid, 2, j, jspin)
                  CALL intgr0(integrand, r0, dx, n_grid, metric(i, j, ispin, jspin))
               END DO
            END DO
         END DO
      END DO
   END SUBROUTINE local_d_spin_ordered_radial_metric

   SUBROUTINE local_d_spin_contract_density(coefficients, metric, rho_d)
      ! Contract one band's l=2 muffin-tin coefficients with the ordered
      ! radial metric. coefficients is indexed as (m, local spin, radial order)
      ! with m=-2,...,+2. rho_d is neither normalized nor symmetrized: its trace
      ! is the physical muffin-tin d weight of this band on the selected atom.
      COMPLEX, INTENT(IN)  :: coefficients(-local_d_l:, :, :)
      REAL,    INTENT(IN)  :: metric(:, :, :, :)
      COMPLEX, INTENT(OUT) :: rho_d(-local_d_l:local_d_l, local_d_n_spins, &
                                    -local_d_l:local_d_l, local_d_n_spins)

      INTEGER :: i, j, m, mp, ispin, jspin, n_radial

      n_radial = SIZE(coefficients, 3)
      IF (SIZE(coefficients, 1) /= local_d_n_orbitals .OR. &
          SIZE(coefficients, 2) /= local_d_n_spins) THEN
         CALL juDFT_error("Local d-spin coefficients must have dimensions (5,2,n_radial)", &
                          calledby="m_local_d_spin_projector_core")
      END IF
      IF (SIZE(metric, 1) /= n_radial .OR. SIZE(metric, 2) /= n_radial .OR. &
          SIZE(metric, 3) /= local_d_n_spins .OR. SIZE(metric, 4) /= local_d_n_spins) THEN
         CALL juDFT_error("Local d-spin metric does not match coefficient dimensions", &
                          calledby="m_local_d_spin_projector_core")
      END IF

      rho_d = CMPLX(0.0, 0.0)
      DO m = -local_d_l, local_d_l
         DO ispin = 1, local_d_n_spins
            DO mp = -local_d_l, local_d_l
               DO jspin = 1, local_d_n_spins
                  DO i = 1, n_radial
                     DO j = 1, n_radial
                        rho_d(m, ispin, mp, jspin) = rho_d(m, ispin, mp, jspin) &
                           + coefficients(m, ispin, i)*CONJG(coefficients(mp, jspin, j)) &
                           * metric(i, j, ispin, jspin)
                     END DO
                  END DO
               END DO
            END DO
         END DO
      END DO
   END SUBROUTINE local_d_spin_contract_density

END MODULE m_local_d_spin_projector_core
