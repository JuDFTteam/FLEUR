!--------------------------------------------------------------------------------
! Independent validation of the local d x spin density-matrix analysis layer.
! No jDOS transformation or production XAS/RIXS code is used here.
!--------------------------------------------------------------------------------

PROGRAM local_d_spin_analysis_test
   USE m_local_d_spin_analysis
   USE m_local_d_spin_projector_core, ONLY: local_d_l, local_d_n_orbitals, local_d_n_spins
   IMPLICIT NONE

   REAL, PARAMETER :: tolerance = 5.0e-12
   COMPLEX :: complex_to_real(-local_d_l:local_d_l, local_d_n_orbitals)
   COMPLEX :: effective_ml(local_d_n_t2g, local_d_n_t2g)
   COMPLEX :: jeff(local_d_n_jeff, local_d_n_jeff)
   COMPLEX :: rho_d(-local_d_l:local_d_l, local_d_n_spins, &
                    -local_d_l:local_d_l, local_d_n_spins)
   COMPLEX :: rho_real(local_d_n_orbitals, local_d_n_spins, &
                       local_d_n_orbitals, local_d_n_spins)
   COMPLEX :: rho_real_reference(local_d_n_orbitals, local_d_n_spins, &
                                 local_d_n_orbitals, local_d_n_spins)
   COMPLEX :: rho_t2g(local_d_n_t2g, local_d_n_spins, local_d_n_t2g, local_d_n_spins)
   COMPLEX :: rho_jeff(local_d_n_jeff, local_d_n_jeff)
   COMPLEX :: physical_l_complex(-local_d_l:local_d_l, -local_d_l:local_d_l, 3)
   COMPLEX :: physical_l_real(local_d_n_orbitals, local_d_n_orbitals, 3)
   COMPLEX :: expected_real_l(local_d_n_orbitals, local_d_n_orbitals, 3)
   COMPLEX :: physical_l_t2g(local_d_n_t2g, local_d_n_t2g, 3)
   COMPLEX :: effective_l_t2g(local_d_n_t2g, local_d_n_t2g, 3)
   COMPLEX :: spin_matrices(local_d_n_spins, local_d_n_spins, 3)
   COMPLEX :: expected_t2g_l(local_d_n_t2g, local_d_n_t2g, 3)
   COMPLEX :: state_t2g(local_d_n_jeff)
   COMPLEX :: j_matrices(local_d_n_jeff, local_d_n_jeff, 3)
   COMPLEX :: j_squared(local_d_n_jeff, local_d_n_jeff)
   COMPLEX :: transformed(local_d_n_jeff, local_d_n_jeff)
   REAL :: orbital_weights(local_d_n_orbitals)
   REAL :: mj_weights(local_d_n_jeff), j_weights(2)
   REAL :: normalized_mj(local_d_n_jeff), normalized_j(2)
   REAL :: spin_expectation(3), orbital_expectation(3)
   REAL :: total_weight, t2g_weight, eg_weight
   REAL :: max_error
   LOGICAL :: normalization_defined
   INTEGER :: failures, channel

   failures = 0
   max_error = 0.0

   CALL local_d_complex_to_real_unitary(complex_to_real)
   CALL check_bound("complex-to-real unitarity", &
                    unitary_error(complex_to_real), tolerance, failures, max_error)
   CALL check_complex("xy coefficient at m=-2", complex_to_real(-2, local_d_orb_xy), &
                      CMPLX(0.0, 1.0/SQRT(2.0)), tolerance, failures, max_error)
   CALL check_complex("xy coefficient at m=+2", complex_to_real(2, local_d_orb_xy), &
                      CMPLX(0.0, -1.0/SQRT(2.0)), tolerance, failures, max_error)
   CALL check_complex("yz coefficient at m=-1", complex_to_real(-1, local_d_orb_yz), &
                      CMPLX(0.0, 1.0/SQRT(2.0)), tolerance, failures, max_error)
   CALL check_complex("xz coefficient at m=+1", complex_to_real(1, local_d_orb_xz), &
                      CMPLX(-1.0/SQRT(2.0), 0.0), tolerance, failures, max_error)

   CALL initialize_real_density(rho_real_reference)
   CALL real_to_complex_density(rho_real_reference, complex_to_real, rho_d)
   CALL local_d_transform_to_real(rho_d, rho_real)
   CALL check_bound("full complex real-orbital density transformation", &
                    MAXVAL(ABS(rho_real - rho_real_reference)), tolerance, failures, max_error)
   CALL local_d_character_weights(rho_real, orbital_weights, t2g_weight, eg_weight, total_weight)
   CALL check_real("xy weight", orbital_weights(local_d_orb_xy), 0.3, tolerance, failures, max_error)
   CALL check_real("yz weight", orbital_weights(local_d_orb_yz), 0.4, tolerance, failures, max_error)
   CALL check_real("xz weight", orbital_weights(local_d_orb_xz), 0.5, tolerance, failures, max_error)
   CALL check_real("x2-y2 weight", orbital_weights(local_d_orb_x2y2), 0.6, tolerance, failures, max_error)
   CALL check_real("z2 weight", orbital_weights(local_d_orb_z2), 0.7, tolerance, failures, max_error)
   CALL check_real("t2g weight", t2g_weight, 1.2, tolerance, failures, max_error)
   CALL check_real("eg weight", eg_weight, 1.3, tolerance, failures, max_error)
   CALL check_real("t2g plus eg equals total d weight", t2g_weight + eg_weight, &
                   total_weight, tolerance, failures, max_error)

   CALL local_d_angular_momentum(physical_l_complex)
   CALL local_d_real_angular_momentum(physical_l_real, physical_l_t2g, effective_l_t2g)
   CALL initialize_expected_real_l(expected_real_l)
   CALL initialize_expected_t2g_l(expected_t2g_l)
   CALL check_bound("complex l=2 angular-momentum algebra", &
                    commutator_error(physical_l_complex(:, :, 1), physical_l_complex(:, :, 2), &
                                     CMPLX(0.0, 1.0)*physical_l_complex(:, :, 3)), &
                    tolerance, failures, max_error)
   CALL check_bound("physical-L action on all real cubic harmonics", &
                    MAXVAL(ABS(physical_l_real - expected_real_l)), &
                    tolerance, failures, max_error)
   CALL check_bound("projected physical-L matrices", MAXVAL(ABS(physical_l_t2g - expected_t2g_l)), &
                    tolerance, failures, max_error)
   CALL check_bound("projected physical-L reversed algebra", &
                    commutator_error(physical_l_t2g(:, :, 1), physical_l_t2g(:, :, 2), &
                                     CMPLX(0.0, -1.0)*physical_l_t2g(:, :, 3)), &
                    tolerance, failures, max_error)
   CALL check_bound("effective-l standard algebra", &
                    commutator_error(effective_l_t2g(:, :, 1), effective_l_t2g(:, :, 2), &
                                     CMPLX(0.0, 1.0)*effective_l_t2g(:, :, 3)), &
                    tolerance, failures, max_error)

   CALL local_d_effective_ml_unitary(effective_ml)
   CALL check_bound("effective-m_l unitarity", unitary_error(effective_ml), &
                    tolerance, failures, max_error)
   transformed = CMPLX(0.0, 0.0)
   transformed(1:3, 1:3) = MATMUL(CONJG(TRANSPOSE(effective_ml)), &
                                  MATMUL(effective_l_t2g(:, :, 3), effective_ml))
   CALL check_bound("effective-m_l l_eff,z eigenvalues", &
                    MAXVAL(ABS(transformed(1:3, 1:3) - diagonal_three([-1.0, 0.0, 1.0]))), &
                    tolerance, failures, max_error)
   transformed = CMPLX(0.0, 0.0)
   transformed(1:3, 1:3) = MATMUL(CONJG(TRANSPOSE(effective_ml)), &
      MATMUL(effective_l_t2g(:, :, 1) + CMPLX(0.0, 1.0)*effective_l_t2g(:, :, 2), effective_ml))
   CALL check_bound("effective-m_l ladder elements", &
                    MAXVAL(ABS(transformed(1:3, 1:3) - expected_ladder_plus())), &
                    tolerance, failures, max_error)
   transformed = CMPLX(0.0, 0.0)
   transformed(1:3, 1:3) = MATMUL(CONJG(TRANSPOSE(effective_ml)), &
                                  MATMUL(physical_l_t2g(:, :, 3), effective_ml))
   CALL check_bound("physical-L sign in effective-m_l basis", &
                    MAXVAL(ABS(transformed(1:3, 1:3) - diagonal_three([1.0, 0.0, -1.0]))), &
                    tolerance, failures, max_error)

   CALL local_d_spin_matrices(spin_matrices)
   CALL check_bound("spin commutator", &
                    commutator_error(spin_matrices(:, :, 1), spin_matrices(:, :, 2), &
                                     CMPLX(0.0, 1.0)*spin_matrices(:, :, 3)), &
                    tolerance, failures, max_error)
   CALL test_spin_expectations(failures, max_error)
   CALL test_orbital_expectations(effective_ml, physical_l_t2g, effective_l_t2g, failures, max_error)

   CALL local_d_jeff_unitary(jeff)
   CALL check_bound("j_eff unitarity", unitary_error(jeff), tolerance, failures, max_error)
   CALL build_total_j(effective_l_t2g, spin_matrices, j_matrices)
   j_squared = MATMUL(j_matrices(:, :, 1), j_matrices(:, :, 1)) &
             + MATMUL(j_matrices(:, :, 2), j_matrices(:, :, 2)) &
             + MATMUL(j_matrices(:, :, 3), j_matrices(:, :, 3))
   transformed = MATMUL(CONJG(TRANSPOSE(jeff)), MATMUL(j_matrices(:, :, 3), jeff))
   CALL check_bound("j_eff J_z eigenvalues", &
                    MAXVAL(ABS(transformed - diagonal_six([-0.5, 0.5, -1.5, -0.5, 0.5, 1.5]))), &
                    tolerance, failures, max_error)
   transformed = MATMUL(CONJG(TRANSPOSE(jeff)), MATMUL(j_squared, jeff))
   CALL check_bound("j_eff J_squared eigenvalues", &
                    MAXVAL(ABS(transformed - diagonal_six([0.75, 0.75, 3.75, 3.75, 3.75, 3.75]))), &
                    tolerance, failures, max_error)

   DO channel = 1, local_d_n_jeff
      state_t2g = jeff(:, channel)
      CALL flat_state_to_t2g_density(state_t2g, rho_t2g)
      CALL local_d_transform_to_jeff(rho_t2g, rho_jeff)
      CALL check_bound("pure ideal j_eff density", &
                       MAXVAL(ABS(rho_jeff - pure_channel_density(channel))), &
                       tolerance, failures, max_error)
      CALL local_d_jeff_weights(rho_t2g, mj_weights, j_weights, normalized_mj, normalized_j, &
                                normalization_defined)
      CALL check_bound("pure ideal j_eff m_j weights", &
                       MAXVAL(ABS(mj_weights - pure_channel_weights(channel))), &
                       tolerance, failures, max_error)
      IF (.NOT. normalization_defined) THEN
         WRITE (*, '(a,i0)') "FAIL: pure ideal j_eff normalization undefined for channel ", channel
         failures = failures + 1
      END IF
      CALL check_bound("pure ideal normalized j_eff weights", &
                       MAXVAL(ABS(normalized_mj - pure_channel_weights(channel))), &
                       tolerance, failures, max_error)
   END DO

   CALL initialize_mixed_t2g_eg_density(jeff(:, local_d_jeff_half_minus), rho_real)
   CALL local_d_character_weights(rho_real, orbital_weights, t2g_weight, eg_weight, total_weight)
   CALL local_d_extract_t2g(rho_real, rho_t2g)
   CALL local_d_jeff_weights(rho_t2g, mj_weights, j_weights, normalized_mj, normalized_j, &
                             normalization_defined)
   CALL check_real("mixed-state raw t2g weight", t2g_weight, 0.75, tolerance, failures, max_error)
   CALL check_real("mixed-state complementary eg weight", eg_weight, 0.25, tolerance, failures, max_error)
   CALL check_real("mixed-state raw j_eff=1/2 weight", j_weights(1), 0.75, tolerance, failures, max_error)
   CALL check_real("mixed-state raw j_eff=3/2 weight", j_weights(2), 0.0, tolerance, failures, max_error)
   CALL check_real("mixed-state t2g-normalized j_eff=1/2 weight", normalized_j(1), 1.0, &
                   tolerance, failures, max_error)
   CALL check_real("mixed-state full local d weight", total_weight, 1.0, tolerance, failures, max_error)

   rho_t2g = CMPLX(0.0, 0.0)
   CALL local_d_jeff_weights(rho_t2g, mj_weights, j_weights, normalized_mj, normalized_j, &
                             normalization_defined)
   IF (normalization_defined .OR. MAXVAL(ABS(normalized_mj)) > tolerance .OR. &
       MAXVAL(ABS(normalized_j)) > tolerance) THEN
      WRITE (*, '(a)') "FAIL: empty t2g normalization must be undefined and zero"
      failures = failures + 1
   ELSE
      WRITE (*, '(a)') "PASS: empty t2g normalization is explicitly undefined"
   END IF

   CALL test_manifold_additivity(jeff, failures, max_error)

   WRITE (*, '(a,es24.16)') "Maximum selected absolute error = ", max_error
   IF (failures /= 0) THEN
      WRITE (*, '(a,i0)') "LOCAL D-SPIN ANALYSIS TEST: FAIL; failures = ", failures
      ERROR STOP 1
   END IF
   WRITE (*, '(a)') "LOCAL D-SPIN ANALYSIS TEST: PASS"

CONTAINS

   SUBROUTINE initialize_real_density(density)
      COMPLEX, INTENT(OUT) :: density(local_d_n_orbitals, local_d_n_spins, &
                                      local_d_n_orbitals, local_d_n_spins)

      density = CMPLX(0.0, 0.0)
      density(local_d_orb_xy, 1, local_d_orb_xy, 1) = CMPLX(0.2, 0.0)
      density(local_d_orb_xy, 2, local_d_orb_xy, 2) = CMPLX(0.1, 0.0)
      density(local_d_orb_yz, 1, local_d_orb_yz, 1) = CMPLX(0.4, 0.0)
      density(local_d_orb_xz, 2, local_d_orb_xz, 2) = CMPLX(0.5, 0.0)
      density(local_d_orb_x2y2, 1, local_d_orb_x2y2, 1) = CMPLX(0.6, 0.0)
      density(local_d_orb_z2, 2, local_d_orb_z2, 2) = CMPLX(0.7, 0.0)
      density(local_d_orb_xy, 1, local_d_orb_yz, 2) = CMPLX(0.03, 0.04)
      density(local_d_orb_yz, 2, local_d_orb_xy, 1) = CMPLX(0.03, -0.04)
   END SUBROUTINE initialize_real_density

   SUBROUTINE real_to_complex_density(real_density, unitary, complex_density)
      COMPLEX, INTENT(IN) :: real_density(local_d_n_orbitals, local_d_n_spins, &
                                           local_d_n_orbitals, local_d_n_spins)
      COMPLEX, INTENT(IN) :: unitary(-local_d_l:local_d_l, local_d_n_orbitals)
      COMPLEX, INTENT(OUT) :: complex_density(-local_d_l:local_d_l, local_d_n_spins, &
                                              -local_d_l:local_d_l, local_d_n_spins)
      INTEGER :: m, mp, orbital, orbital_p, ispin, jspin

      complex_density = CMPLX(0.0, 0.0)
      DO m = -local_d_l, local_d_l
         DO ispin = 1, local_d_n_spins
            DO mp = -local_d_l, local_d_l
               DO jspin = 1, local_d_n_spins
                  DO orbital = 1, local_d_n_orbitals
                     DO orbital_p = 1, local_d_n_orbitals
                        complex_density(m, ispin, mp, jspin) = complex_density(m, ispin, mp, jspin) &
                           + unitary(m, orbital)*real_density(orbital, ispin, orbital_p, jspin) &
                           * CONJG(unitary(mp, orbital_p))
                     END DO
                  END DO
               END DO
            END DO
         END DO
      END DO
   END SUBROUTINE real_to_complex_density

   SUBROUTINE initialize_expected_t2g_l(expected)
      COMPLEX, INTENT(OUT) :: expected(local_d_n_t2g, local_d_n_t2g, 3)

      expected = CMPLX(0.0, 0.0)
      expected(1, 3, 1) = CMPLX(0.0, -1.0)
      expected(3, 1, 1) = CMPLX(0.0, 1.0)
      expected(1, 2, 2) = CMPLX(0.0, 1.0)
      expected(2, 1, 2) = CMPLX(0.0, -1.0)
      expected(2, 3, 3) = CMPLX(0.0, 1.0)
      expected(3, 2, 3) = CMPLX(0.0, -1.0)
   END SUBROUTINE initialize_expected_t2g_l

   SUBROUTINE initialize_expected_real_l(expected)
      ! Independent action of L on (xy,yz,xz,x2-y2,z2), derived by applying
      ! L_z|m>=m|m> and the l=2 ladder operators to the declared real states.
      COMPLEX, INTENT(OUT) :: expected(local_d_n_orbitals, local_d_n_orbitals, 3)
      REAL :: sqrt_three

      sqrt_three = SQRT(3.0)
      expected = CMPLX(0.0, 0.0)

      expected(local_d_orb_xy, local_d_orb_xz, 1) = CMPLX(0.0, -1.0)
      expected(local_d_orb_xz, local_d_orb_xy, 1) = CMPLX(0.0, 1.0)
      expected(local_d_orb_yz, local_d_orb_x2y2, 1) = CMPLX(0.0, -1.0)
      expected(local_d_orb_x2y2, local_d_orb_yz, 1) = CMPLX(0.0, 1.0)
      expected(local_d_orb_yz, local_d_orb_z2, 1) = CMPLX(0.0, -sqrt_three)
      expected(local_d_orb_z2, local_d_orb_yz, 1) = CMPLX(0.0, sqrt_three)

      expected(local_d_orb_xy, local_d_orb_yz, 2) = CMPLX(0.0, 1.0)
      expected(local_d_orb_yz, local_d_orb_xy, 2) = CMPLX(0.0, -1.0)
      expected(local_d_orb_xz, local_d_orb_x2y2, 2) = CMPLX(0.0, -1.0)
      expected(local_d_orb_x2y2, local_d_orb_xz, 2) = CMPLX(0.0, 1.0)
      expected(local_d_orb_xz, local_d_orb_z2, 2) = CMPLX(0.0, sqrt_three)
      expected(local_d_orb_z2, local_d_orb_xz, 2) = CMPLX(0.0, -sqrt_three)

      expected(local_d_orb_xy, local_d_orb_x2y2, 3) = CMPLX(0.0, 2.0)
      expected(local_d_orb_x2y2, local_d_orb_xy, 3) = CMPLX(0.0, -2.0)
      expected(local_d_orb_yz, local_d_orb_xz, 3) = CMPLX(0.0, 1.0)
      expected(local_d_orb_xz, local_d_orb_yz, 3) = CMPLX(0.0, -1.0)
   END SUBROUTINE initialize_expected_real_l

   SUBROUTINE test_spin_expectations(n_failures, largest_error)
      INTEGER, INTENT(INOUT) :: n_failures
      REAL, INTENT(INOUT) :: largest_error
      COMPLEX :: density(local_d_n_t2g, local_d_n_spins, local_d_n_t2g, local_d_n_spins)
      REAL :: expectation(3)

      density = CMPLX(0.0, 0.0)
      density(1, 1, 1, 1) = CMPLX(2.0, 0.0)
      CALL local_d_spin_expectation(density, expectation)
      CALL check_vector("native-frame pure spin up", expectation, [0.0, 0.0, 1.0], &
                        tolerance, n_failures, largest_error)

      density = CMPLX(0.0, 0.0)
      density(1, 2, 1, 2) = CMPLX(2.0, 0.0)
      CALL local_d_spin_expectation(density, expectation)
      CALL check_vector("native-frame pure spin down", expectation, [0.0, 0.0, -1.0], &
                        tolerance, n_failures, largest_error)

      density = CMPLX(0.0, 0.0)
      density(1, 1, 1, 1) = CMPLX(0.5, 0.0)
      density(1, 2, 1, 2) = CMPLX(0.5, 0.0)
      density(1, 1, 1, 2) = CMPLX(0.5, 0.0)
      density(1, 2, 1, 1) = CMPLX(0.5, 0.0)
      CALL local_d_spin_expectation(density, expectation)
      CALL check_vector("native-frame x spinor", expectation, [0.5, 0.0, 0.0], &
                        tolerance, n_failures, largest_error)

      density(1, 1, 1, 2) = CMPLX(0.0, -0.5)
      density(1, 2, 1, 1) = CMPLX(0.0, 0.5)
      CALL local_d_spin_expectation(density, expectation)
      CALL check_vector("native-frame y spinor", expectation, [0.0, 0.5, 0.0], &
                        tolerance, n_failures, largest_error)
   END SUBROUTINE test_spin_expectations

   SUBROUTINE test_orbital_expectations(ml_unitary, physical_l, effective_l, n_failures, largest_error)
      COMPLEX, INTENT(IN) :: ml_unitary(local_d_n_t2g, local_d_n_t2g)
      COMPLEX, INTENT(IN) :: physical_l(local_d_n_t2g, local_d_n_t2g, 3)
      COMPLEX, INTENT(IN) :: effective_l(local_d_n_t2g, local_d_n_t2g, 3)
      INTEGER, INTENT(INOUT) :: n_failures
      REAL, INTENT(INOUT) :: largest_error
      COMPLEX :: density(local_d_n_t2g, local_d_n_spins, local_d_n_t2g, local_d_n_spins)
      REAL :: expectation(3)
      INTEGER :: orbital, orbital_p

      density = CMPLX(0.0, 0.0)
      DO orbital = 1, local_d_n_t2g
         DO orbital_p = 1, local_d_n_t2g
            density(orbital, 1, orbital_p, 1) = ml_unitary(orbital, 3)*CONJG(ml_unitary(orbital_p, 3))
         END DO
      END DO
      CALL local_d_orbital_expectation(density, physical_l, expectation)
      CALL check_vector("m_eff=+1 physical-L expectation", expectation, [0.0, 0.0, -1.0], &
                        tolerance, n_failures, largest_error)
      CALL local_d_orbital_expectation(density, effective_l, expectation)
      CALL check_vector("m_eff=+1 effective-l expectation", expectation, [0.0, 0.0, 1.0], &
                        tolerance, n_failures, largest_error)
   END SUBROUTINE test_orbital_expectations

   SUBROUTINE build_total_j(effective_l, spin, total_j)
      COMPLEX, INTENT(IN) :: effective_l(local_d_n_t2g, local_d_n_t2g, 3)
      COMPLEX, INTENT(IN) :: spin(local_d_n_spins, local_d_n_spins, 3)
      COMPLEX, INTENT(OUT) :: total_j(local_d_n_jeff, local_d_n_jeff, 3)
      INTEGER :: direction, orbital, orbital_p, ispin, jspin, row, column

      total_j = CMPLX(0.0, 0.0)
      DO direction = 1, 3
         DO orbital = 1, local_d_n_t2g
            DO ispin = 1, local_d_n_spins
               row = (orbital - 1)*local_d_n_spins + ispin
               DO orbital_p = 1, local_d_n_t2g
                  DO jspin = 1, local_d_n_spins
                     column = (orbital_p - 1)*local_d_n_spins + jspin
                     IF (ispin == jspin) total_j(row, column, direction) = &
                        total_j(row, column, direction) + effective_l(orbital, orbital_p, direction)
                     IF (orbital == orbital_p) total_j(row, column, direction) = &
                        total_j(row, column, direction) + spin(ispin, jspin, direction)
                  END DO
               END DO
            END DO
         END DO
      END DO
   END SUBROUTINE build_total_j

   SUBROUTINE flat_state_to_t2g_density(state, density)
      COMPLEX, INTENT(IN) :: state(local_d_n_jeff)
      COMPLEX, INTENT(OUT) :: density(local_d_n_t2g, local_d_n_spins, &
                                      local_d_n_t2g, local_d_n_spins)
      INTEGER :: orbital, orbital_p, ispin, jspin, row, column

      DO orbital = 1, local_d_n_t2g
         DO ispin = 1, local_d_n_spins
            row = (orbital - 1)*local_d_n_spins + ispin
            DO orbital_p = 1, local_d_n_t2g
               DO jspin = 1, local_d_n_spins
                  column = (orbital_p - 1)*local_d_n_spins + jspin
                  density(orbital, ispin, orbital_p, jspin) = state(row)*CONJG(state(column))
               END DO
            END DO
         END DO
      END DO
   END SUBROUTINE flat_state_to_t2g_density

   SUBROUTINE initialize_mixed_t2g_eg_density(t2g_state, density)
      COMPLEX, INTENT(IN) :: t2g_state(local_d_n_jeff)
      COMPLEX, INTENT(OUT) :: density(local_d_n_orbitals, local_d_n_spins, &
                                      local_d_n_orbitals, local_d_n_spins)
      COMPLEX :: t2g_density(local_d_n_t2g, local_d_n_spins, local_d_n_t2g, local_d_n_spins)

      CALL flat_state_to_t2g_density(t2g_state, t2g_density)
      density = CMPLX(0.0, 0.0)
      density(1:local_d_n_t2g, :, 1:local_d_n_t2g, :) = 0.75*t2g_density
      density(local_d_orb_x2y2, 1, local_d_orb_x2y2, 1) = CMPLX(0.25, 0.0)
   END SUBROUTINE initialize_mixed_t2g_eg_density

   SUBROUTINE test_manifold_additivity(jeff_unitary, n_failures, largest_error)
      COMPLEX, INTENT(IN) :: jeff_unitary(local_d_n_jeff, local_d_n_jeff)
      INTEGER, INTENT(INOUT) :: n_failures
      REAL, INTENT(INOUT) :: largest_error
      COMPLEX :: density_a(local_d_n_t2g, local_d_n_spins, local_d_n_t2g, local_d_n_spins)
      COMPLEX :: density_b(local_d_n_t2g, local_d_n_spins, local_d_n_t2g, local_d_n_spins)
      REAL :: mj_a(local_d_n_jeff), mj_b(local_d_n_jeff), mj_sum(local_d_n_jeff)
      REAL :: j_a(2), j_b(2), j_sum(2)

      CALL flat_state_to_t2g_density(jeff_unitary(:, local_d_jeff_half_minus), density_a)
      CALL flat_state_to_t2g_density(jeff_unitary(:, local_d_jeff_three_plus_three), density_b)
      density_a = 0.4*density_a
      density_b = 0.6*density_b
      CALL local_d_jeff_weights(density_a, mj_a, j_a)
      CALL local_d_jeff_weights(density_b, mj_b, j_b)
      CALL local_d_jeff_weights(density_a + density_b, mj_sum, j_sum)
      CALL check_bound("manifold-summed m_j weight additivity", MAXVAL(ABS(mj_sum - mj_a - mj_b)), &
                       tolerance, n_failures, largest_error)
      CALL check_bound("manifold-summed j_eff weight additivity", MAXVAL(ABS(j_sum - j_a - j_b)), &
                       tolerance, n_failures, largest_error)
   END SUBROUTINE test_manifold_additivity

   FUNCTION pure_channel_density(channel) RESULT(density)
      INTEGER, INTENT(IN) :: channel
      COMPLEX :: density(local_d_n_jeff, local_d_n_jeff)

      density = CMPLX(0.0, 0.0)
      density(channel, channel) = CMPLX(1.0, 0.0)
   END FUNCTION pure_channel_density

   FUNCTION pure_channel_weights(channel) RESULT(weights)
      INTEGER, INTENT(IN) :: channel
      REAL :: weights(local_d_n_jeff)

      weights = 0.0
      weights(channel) = 1.0
   END FUNCTION pure_channel_weights

   FUNCTION diagonal_three(values) RESULT(matrix)
      REAL, INTENT(IN) :: values(3)
      COMPLEX :: matrix(3, 3)
      INTEGER :: i

      matrix = CMPLX(0.0, 0.0)
      DO i = 1, 3
         matrix(i, i) = CMPLX(values(i), 0.0)
      END DO
   END FUNCTION diagonal_three

   FUNCTION diagonal_six(values) RESULT(matrix)
      REAL, INTENT(IN) :: values(6)
      COMPLEX :: matrix(6, 6)
      INTEGER :: i

      matrix = CMPLX(0.0, 0.0)
      DO i = 1, 6
         matrix(i, i) = CMPLX(values(i), 0.0)
      END DO
   END FUNCTION diagonal_six

   FUNCTION expected_ladder_plus() RESULT(matrix)
      COMPLEX :: matrix(3, 3)

      matrix = CMPLX(0.0, 0.0)
      matrix(2, 1) = CMPLX(SQRT(2.0), 0.0)
      matrix(3, 2) = CMPLX(SQRT(2.0), 0.0)
   END FUNCTION expected_ladder_plus

   FUNCTION commutator_error(left, right, expected) RESULT(error)
      COMPLEX, INTENT(IN) :: left(:, :), right(:, :), expected(:, :)
      REAL :: error

      error = MAXVAL(ABS(MATMUL(left, right) - MATMUL(right, left) - expected))
   END FUNCTION commutator_error

   FUNCTION unitary_error(matrix) RESULT(error)
      COMPLEX, INTENT(IN) :: matrix(:, :)
      COMPLEX :: identity(SIZE(matrix, 2), SIZE(matrix, 2))
      REAL :: error
      INTEGER :: i

      identity = CMPLX(0.0, 0.0)
      DO i = 1, SIZE(identity, 1)
         identity(i, i) = CMPLX(1.0, 0.0)
      END DO
      error = MAXVAL(ABS(MATMUL(CONJG(TRANSPOSE(matrix)), matrix) - identity))
   END FUNCTION unitary_error

   SUBROUTINE check_bound(name, value, bound, n_failures, largest_error)
      CHARACTER(LEN=*), INTENT(IN) :: name
      REAL, INTENT(IN) :: value, bound
      INTEGER, INTENT(INOUT) :: n_failures
      REAL, INTENT(INOUT) :: largest_error

      largest_error = MAX(largest_error, value)
      IF (value > bound) THEN
         WRITE (*, '(a,a,a,es24.16,a,es24.16)') "FAIL: ", TRIM(name), " error ", value, " > ", bound
         n_failures = n_failures + 1
      ELSE
         WRITE (*, '(a,a,a,es24.16)') "PASS: ", TRIM(name), " error = ", value
      END IF
   END SUBROUTINE check_bound

   SUBROUTINE check_real(name, value, expected, bound, n_failures, largest_error)
      CHARACTER(LEN=*), INTENT(IN) :: name
      REAL, INTENT(IN) :: value, expected, bound
      INTEGER, INTENT(INOUT) :: n_failures
      REAL, INTENT(INOUT) :: largest_error

      CALL check_bound(name, ABS(value - expected), bound, n_failures, largest_error)
   END SUBROUTINE check_real

   SUBROUTINE check_complex(name, value, expected, bound, n_failures, largest_error)
      CHARACTER(LEN=*), INTENT(IN) :: name
      COMPLEX, INTENT(IN) :: value, expected
      REAL, INTENT(IN) :: bound
      INTEGER, INTENT(INOUT) :: n_failures
      REAL, INTENT(INOUT) :: largest_error

      CALL check_bound(name, ABS(value - expected), bound, n_failures, largest_error)
   END SUBROUTINE check_complex

   SUBROUTINE check_vector(name, value, expected, bound, n_failures, largest_error)
      CHARACTER(LEN=*), INTENT(IN) :: name
      REAL, INTENT(IN) :: value(3), expected(3), bound
      INTEGER, INTENT(INOUT) :: n_failures
      REAL, INTENT(INOUT) :: largest_error

      CALL check_bound(name, MAXVAL(ABS(value - expected)), bound, n_failures, largest_error)
   END SUBROUTINE check_vector

END PROGRAM local_d_spin_analysis_test
