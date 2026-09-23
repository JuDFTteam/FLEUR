!--------------------------------------------------------------------------------
! Copyright (c) 2026 Peter Gruenberg Institut, Forschungszentrum Juelich, Germany
! This file is part of FLEUR and available as free software under the conditions
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------

MODULE m_local_d_spin_projector
   USE m_juDFT, ONLY: juDFT_error
   USE m_local_d_spin_projector_core, ONLY: local_d_l, local_d_n_spins, &
                                            local_d_spin_contract_density, &
                                            local_d_spin_ordered_radial_metric
   USE m_types_abc, ONLY: t_abc
   USE m_types_atoms, ONLY: t_atoms
   USE m_types_radfun, ONLY: t_radfun
   IMPLICIT NONE
   PRIVATE

   PUBLIC :: local_d_spin_density_matrix

CONTAINS

   SUBROUTINE local_d_spin_density_matrix(abc_spin, radfun, atoms, itype, band, iatom_l, rho_d, ordered_metric)
      ! Return the unnormalized local d x spin reduced density matrix for one
      ! band and one atom. Orbital indices use FLEUR's complex Y_2m ordering
      ! m=-2,...,+2; spin indices 1,2 are the native local muffin-tin spin
      ! components produced by t_abc%calc_abc.
      !
      ! The ordered metric is rebuilt from radfun%r instead of using the cached
      ! cross-spin blocks of radfun%integral. The latter currently impose
      ! radial-index symmetry inside a fixed cross-spin block, which is not
      ! valid when the two spin components have different radial bases.
      TYPE(t_abc),    INTENT(IN) :: abc_spin(:)
      TYPE(t_radfun), INTENT(IN) :: radfun
      TYPE(t_atoms),  INTENT(IN) :: atoms
      INTEGER, INTENT(IN) :: itype, band, iatom_l
      COMPLEX, INTENT(OUT) :: rho_d(-local_d_l:local_d_l, local_d_n_spins, &
                                    -local_d_l:local_d_l, local_d_n_spins)
      REAL, ALLOCATABLE, OPTIONAL, INTENT(OUT) :: ordered_metric(:, :, :, :)

      COMPLEX, ALLOCATABLE :: coefficients(:, :, :)
      REAL, ALLOCATABLE :: metric(:, :, :, :)
      INTEGER :: i, lm, m, n_grid, n_radial, ispin

      IF (itype < 1 .OR. itype > atoms%ntype) THEN
         CALL juDFT_error("Invalid atom type in local d-spin projector", calledby="m_local_d_spin_projector")
      END IF
      IF (SIZE(abc_spin) /= local_d_n_spins) THEN
         CALL juDFT_error("Local d-spin projector requires two t_abc spin components", &
                          calledby="m_local_d_spin_projector")
      END IF
      IF (.NOT. ALLOCATED(radfun%n_r) .OR. .NOT. ALLOCATED(radfun%r)) THEN
         CALL juDFT_error("Radial functions are not initialized in local d-spin projector", &
                          calledby="m_local_d_spin_projector")
      END IF
      IF (LBOUND(radfun%n_r, 1) > local_d_l .OR. UBOUND(radfun%n_r, 1) < local_d_l) THEN
         CALL juDFT_error("Radial-function metadata does not contain l=2", calledby="m_local_d_spin_projector")
      END IF

      n_radial = radfun%n_r(local_d_l)
      n_grid = atoms%jri(itype)
      IF (n_radial < 1 .OR. n_radial > SIZE(radfun%r, 3)) THEN
         CALL juDFT_error("Invalid l=2 radial-function count in local d-spin projector", &
                          calledby="m_local_d_spin_projector")
      END IF
      IF (n_grid < 2 .OR. n_grid > SIZE(radfun%r, 1) .OR. SIZE(radfun%r, 2) /= 2 .OR. &
          UBOUND(radfun%r, 4) < local_d_l .OR. SIZE(radfun%r, 5) < local_d_n_spins) THEN
         CALL juDFT_error("Radial-function storage is inconsistent in local d-spin projector", &
                          calledby="m_local_d_spin_projector")
      END IF

      DO ispin = 1, local_d_n_spins
         IF (.NOT. ALLOCATED(abc_spin(ispin)%cof) .OR. .NOT. ALLOCATED(abc_spin(ispin)%n_r)) THEN
            CALL juDFT_error("abc coefficients are not initialized in local d-spin projector", &
                             calledby="m_local_d_spin_projector")
         END IF
         IF (band < 1 .OR. band > SIZE(abc_spin(ispin)%cof, 1)) THEN
            CALL juDFT_error("Invalid band in local d-spin projector", calledby="m_local_d_spin_projector")
         END IF
         IF (iatom_l < 1 .OR. iatom_l > SIZE(abc_spin(ispin)%cof, 4)) THEN
            CALL juDFT_error("Invalid equivalent-atom index in local d-spin projector", &
                             calledby="m_local_d_spin_projector")
         END IF
         IF (LBOUND(abc_spin(ispin)%n_r, 1) > local_d_l .OR. &
             UBOUND(abc_spin(ispin)%n_r, 1) < local_d_l .OR. &
             abc_spin(ispin)%n_r(local_d_l) /= n_radial .OR. &
             SIZE(abc_spin(ispin)%cof, 3) < n_radial) THEN
            CALL juDFT_error("abc and radial l=2 orders are inconsistent", calledby="m_local_d_spin_projector")
         END IF
         IF (LBOUND(abc_spin(ispin)%cof, 2) > local_d_l**2 .OR. &
             UBOUND(abc_spin(ispin)%cof, 2) < local_d_l*(local_d_l + 2)) THEN
            CALL juDFT_error("abc coefficients do not contain the complete l=2 block", &
                             calledby="m_local_d_spin_projector")
         END IF
      END DO

      ALLOCATE(coefficients(-local_d_l:local_d_l, local_d_n_spins, n_radial), SOURCE=CMPLX(0.0, 0.0))
      DO ispin = 1, local_d_n_spins
         DO m = -local_d_l, local_d_l
            lm = local_d_l*(local_d_l + 1) + m
            DO i = 1, n_radial
               coefficients(m, ispin, i) = abc_spin(ispin)%cof(band, lm, i, iatom_l)
            END DO
         END DO
      END DO

      ALLOCATE(metric(n_radial, n_radial, local_d_n_spins, local_d_n_spins))
      CALL local_d_spin_ordered_radial_metric( &
         radfun%r(1:n_grid, 1:2, 1:n_radial, local_d_l, 1:local_d_n_spins), &
         atoms%rmsh(1, itype), atoms%dx(itype), n_grid, metric)
      CALL local_d_spin_contract_density(coefficients, metric, rho_d)

      IF (PRESENT(ordered_metric)) THEN
         ALLOCATE(ordered_metric(n_radial, n_radial, local_d_n_spins, local_d_n_spins), SOURCE=metric)
      END IF
   END SUBROUTINE local_d_spin_density_matrix

END MODULE m_local_d_spin_projector
