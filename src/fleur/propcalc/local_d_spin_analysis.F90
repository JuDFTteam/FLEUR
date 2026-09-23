!--------------------------------------------------------------------------------
! Copyright (c) 2026 Peter Gruenberg Institut, Forschungszentrum Juelich, Germany
! This file is part of FLEUR and available as free software under the conditions
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------

MODULE m_local_d_spin_analysis
   USE m_juDFT, ONLY: juDFT_error
   USE m_local_d_spin_projector_core, ONLY: local_d_l, local_d_n_orbitals, local_d_n_spins
   IMPLICIT NONE
   PRIVATE

   INTEGER, PARAMETER, PUBLIC :: local_d_orb_xy = 1
   INTEGER, PARAMETER, PUBLIC :: local_d_orb_yz = 2
   INTEGER, PARAMETER, PUBLIC :: local_d_orb_xz = 3
   INTEGER, PARAMETER, PUBLIC :: local_d_orb_x2y2 = 4
   INTEGER, PARAMETER, PUBLIC :: local_d_orb_z2 = 5
   INTEGER, PARAMETER, PUBLIC :: local_d_n_t2g = 3
   INTEGER, PARAMETER, PUBLIC :: local_d_n_jeff = 6

   INTEGER, PARAMETER, PUBLIC :: local_d_jeff_half_minus = 1
   INTEGER, PARAMETER, PUBLIC :: local_d_jeff_half_plus = 2
   INTEGER, PARAMETER, PUBLIC :: local_d_jeff_three_minus_three = 3
   INTEGER, PARAMETER, PUBLIC :: local_d_jeff_three_minus_one = 4
   INTEGER, PARAMETER, PUBLIC :: local_d_jeff_three_plus_one = 5
   INTEGER, PARAMETER, PUBLIC :: local_d_jeff_three_plus_three = 6

   PUBLIC :: local_d_complex_to_real_unitary
   PUBLIC :: local_d_transform_to_real
   PUBLIC :: local_d_extract_t2g
   PUBLIC :: local_d_character_weights
   PUBLIC :: local_d_angular_momentum
   PUBLIC :: local_d_real_angular_momentum
   PUBLIC :: local_d_effective_ml_unitary
   PUBLIC :: local_d_spin_matrices
   PUBLIC :: local_d_spin_expectation
   PUBLIC :: local_d_orbital_expectation
   PUBLIC :: local_d_jeff_unitary
   PUBLIC :: local_d_transform_to_jeff
   PUBLIC :: local_d_jeff_weights

CONTAINS

   SUBROUTINE local_d_complex_to_real_unitary(unitary)
      ! Columns express the real cubic harmonics in FLEUR's complex Y_2m
      ! basis. Complex rows are ordered m=-2,-1,0,+1,+2; real columns are
      ! xy, yz, xz, x2-y2, z2.
      COMPLEX, INTENT(OUT) :: unitary(-local_d_l:local_d_l, local_d_n_orbitals)

      REAL :: inverse_sqrt_two

      inverse_sqrt_two = 1.0/SQRT(2.0)
      unitary = CMPLX(0.0, 0.0)
      unitary(-2, local_d_orb_xy) = CMPLX(0.0, inverse_sqrt_two)
      unitary( 2, local_d_orb_xy) = CMPLX(0.0, -inverse_sqrt_two)
      unitary(-1, local_d_orb_yz) = CMPLX(0.0, inverse_sqrt_two)
      unitary( 1, local_d_orb_yz) = CMPLX(0.0, inverse_sqrt_two)
      unitary(-1, local_d_orb_xz) = CMPLX(inverse_sqrt_two, 0.0)
      unitary( 1, local_d_orb_xz) = CMPLX(-inverse_sqrt_two, 0.0)
      unitary(-2, local_d_orb_x2y2) = CMPLX(inverse_sqrt_two, 0.0)
      unitary( 2, local_d_orb_x2y2) = CMPLX(inverse_sqrt_two, 0.0)
      unitary( 0, local_d_orb_z2) = CMPLX(1.0, 0.0)
   END SUBROUTINE local_d_complex_to_real_unitary

   SUBROUTINE local_d_transform_to_real(rho_d, rho_real)
      ! rho_real=(U^dagger x I_spin) rho_d (U x I_spin). Neither the full
      ! density matrix nor any subspace is normalized here.
      COMPLEX, INTENT(IN) :: rho_d(-local_d_l:, :, -local_d_l:, :)
      COMPLEX, INTENT(OUT) :: rho_real(local_d_n_orbitals, local_d_n_spins, &
                                       local_d_n_orbitals, local_d_n_spins)

      COMPLEX :: unitary(-local_d_l:local_d_l, local_d_n_orbitals)
      INTEGER :: m, mp, orbital, orbital_p, ispin, jspin

      CALL validate_d_spin_density(rho_d, local_d_n_orbitals, "complex d-spin density")
      CALL local_d_complex_to_real_unitary(unitary)
      rho_real = CMPLX(0.0, 0.0)
      DO orbital = 1, local_d_n_orbitals
         DO ispin = 1, local_d_n_spins
            DO orbital_p = 1, local_d_n_orbitals
               DO jspin = 1, local_d_n_spins
                  DO m = -local_d_l, local_d_l
                     DO mp = -local_d_l, local_d_l
                        rho_real(orbital, ispin, orbital_p, jspin) = &
                           rho_real(orbital, ispin, orbital_p, jspin) &
                           + CONJG(unitary(m, orbital))*rho_d(m, ispin, mp, jspin) &
                           * unitary(mp, orbital_p)
                     END DO
                  END DO
               END DO
            END DO
         END DO
      END DO
   END SUBROUTINE local_d_transform_to_real

   SUBROUTINE local_d_extract_t2g(rho_real, rho_t2g)
      ! t2g orbital order is xy, yz, xz. The complementary eg orbitals are
      ! x2-y2 and z2 and remain represented by the full real-basis density.
      COMPLEX, INTENT(IN) :: rho_real(:, :, :, :)
      COMPLEX, INTENT(OUT) :: rho_t2g(local_d_n_t2g, local_d_n_spins, &
                                      local_d_n_t2g, local_d_n_spins)

      CALL validate_d_spin_density(rho_real, local_d_n_orbitals, "real d-spin density")
      rho_t2g = rho_real(1:local_d_n_t2g, :, 1:local_d_n_t2g, :)
   END SUBROUTINE local_d_extract_t2g

   SUBROUTINE local_d_character_weights(rho_real, orbital_weights, t2g_weight, eg_weight, total_weight)
      COMPLEX, INTENT(IN) :: rho_real(:, :, :, :)
      REAL, INTENT(OUT) :: orbital_weights(local_d_n_orbitals)
      REAL, INTENT(OUT) :: t2g_weight, eg_weight, total_weight

      INTEGER :: orbital, ispin

      CALL validate_d_spin_density(rho_real, local_d_n_orbitals, "real d-spin density")
      orbital_weights = 0.0
      DO orbital = 1, local_d_n_orbitals
         DO ispin = 1, local_d_n_spins
            orbital_weights(orbital) = orbital_weights(orbital) &
                                     + REAL(rho_real(orbital, ispin, orbital, ispin))
         END DO
      END DO
      t2g_weight = SUM(orbital_weights(1:local_d_n_t2g))
      eg_weight = SUM(orbital_weights(local_d_n_t2g + 1:local_d_n_orbitals))
      total_weight = SUM(orbital_weights)
   END SUBROUTINE local_d_character_weights

   SUBROUTINE local_d_angular_momentum(angular_momentum)
      ! Construct L from the l=2 ladder algebra in the ordered complex basis.
      ! hbar=1 throughout this analysis module.
      COMPLEX, INTENT(OUT) :: angular_momentum(-local_d_l:local_d_l, &
                                               -local_d_l:local_d_l, 3)

      COMPLEX :: ladder_plus(-local_d_l:local_d_l, -local_d_l:local_d_l)
      COMPLEX :: ladder_minus(-local_d_l:local_d_l, -local_d_l:local_d_l)
      INTEGER :: m

      ladder_plus = CMPLX(0.0, 0.0)
      DO m = -local_d_l, local_d_l - 1
         ladder_plus(m + 1, m) = CMPLX(SQRT(REAL(local_d_l*(local_d_l + 1) - m*(m + 1))), 0.0)
      END DO
      ladder_minus = CONJG(TRANSPOSE(ladder_plus))
      angular_momentum(:, :, 1) = 0.5*(ladder_plus + ladder_minus)
      angular_momentum(:, :, 2) = CMPLX(0.0, -0.5)*(ladder_plus - ladder_minus)
      angular_momentum(:, :, 3) = CMPLX(0.0, 0.0)
      DO m = -local_d_l, local_d_l
         angular_momentum(m, m, 3) = CMPLX(REAL(m), 0.0)
      END DO
   END SUBROUTINE local_d_angular_momentum

   SUBROUTINE local_d_real_angular_momentum(physical_l_real, physical_l_t2g, effective_l_t2g)
      ! Projected physical L obeys [Lx,Ly]=-iLz in t2g. The effective
      ! l=1 operator is therefore defined as l_eff=-L_projected.
      COMPLEX, INTENT(OUT) :: physical_l_real(local_d_n_orbitals, local_d_n_orbitals, 3)
      COMPLEX, INTENT(OUT) :: physical_l_t2g(local_d_n_t2g, local_d_n_t2g, 3)
      COMPLEX, INTENT(OUT) :: effective_l_t2g(local_d_n_t2g, local_d_n_t2g, 3)

      COMPLEX :: angular_momentum(-local_d_l:local_d_l, -local_d_l:local_d_l, 3)
      COMPLEX :: unitary(-local_d_l:local_d_l, local_d_n_orbitals)
      INTEGER :: direction

      CALL local_d_angular_momentum(angular_momentum)
      CALL local_d_complex_to_real_unitary(unitary)
      DO direction = 1, 3
         physical_l_real(:, :, direction) = MATMUL(CONJG(TRANSPOSE(unitary)), &
                                                   MATMUL(angular_momentum(:, :, direction), unitary))
         physical_l_t2g(:, :, direction) = physical_l_real(1:local_d_n_t2g, 1:local_d_n_t2g, direction)
      END DO
      effective_l_t2g = -physical_l_t2g
   END SUBROUTINE local_d_real_angular_momentum

   SUBROUTINE local_d_effective_ml_unitary(unitary)
      ! Columns are m_eff=-1,0,+1 in the real t2g order xy,yz,xz.
      COMPLEX, INTENT(OUT) :: unitary(local_d_n_t2g, local_d_n_t2g)

      REAL :: inverse_sqrt_two

      inverse_sqrt_two = 1.0/SQRT(2.0)
      unitary = CMPLX(0.0, 0.0)
      unitary(local_d_orb_yz, 1) = CMPLX(inverse_sqrt_two, 0.0)
      unitary(local_d_orb_xz, 1) = CMPLX(0.0, -inverse_sqrt_two)
      unitary(local_d_orb_xy, 2) = CMPLX(1.0, 0.0)
      unitary(local_d_orb_yz, 3) = CMPLX(-inverse_sqrt_two, 0.0)
      unitary(local_d_orb_xz, 3) = CMPLX(0.0, -inverse_sqrt_two)
   END SUBROUTINE local_d_effective_ml_unitary

   SUBROUTINE local_d_spin_matrices(spin_matrices)
      ! Native local muffin-tin spin frame; no global-spin rotation is made.
      COMPLEX, INTENT(OUT) :: spin_matrices(local_d_n_spins, local_d_n_spins, 3)

      spin_matrices = CMPLX(0.0, 0.0)
      spin_matrices(1, 2, 1) = CMPLX(0.5, 0.0)
      spin_matrices(2, 1, 1) = CMPLX(0.5, 0.0)
      spin_matrices(1, 2, 2) = CMPLX(0.0, -0.5)
      spin_matrices(2, 1, 2) = CMPLX(0.0, 0.5)
      spin_matrices(1, 1, 3) = CMPLX(0.5, 0.0)
      spin_matrices(2, 2, 3) = CMPLX(-0.5, 0.0)
   END SUBROUTINE local_d_spin_matrices

   SUBROUTINE local_d_spin_expectation(rho, expectation)
      ! Works for either the full five-orbital density or an extracted orbital
      ! subspace. The result is raw (not divided by the subspace trace).
      COMPLEX, INTENT(IN) :: rho(:, :, :, :)
      REAL, INTENT(OUT) :: expectation(3)

      COMPLEX :: spin_matrices(local_d_n_spins, local_d_n_spins, 3)
      COMPLEX :: value
      INTEGER :: direction, orbital, ispin, jspin

      CALL validate_square_spin_density(rho, "spin expectation density")
      CALL local_d_spin_matrices(spin_matrices)
      expectation = 0.0
      DO direction = 1, 3
         value = CMPLX(0.0, 0.0)
         DO orbital = 1, SIZE(rho, 1)
            DO ispin = 1, local_d_n_spins
               DO jspin = 1, local_d_n_spins
                  value = value + rho(orbital, ispin, orbital, jspin)*spin_matrices(jspin, ispin, direction)
               END DO
            END DO
         END DO
         expectation(direction) = REAL(value)
      END DO
   END SUBROUTINE local_d_spin_expectation

   SUBROUTINE local_d_orbital_expectation(rho, orbital_matrices, expectation)
      ! Raw expectation of an orbital operator in the supplied orbital
      ! subspace. orbital_matrices can be physical L or effective l_eff.
      COMPLEX, INTENT(IN) :: rho(:, :, :, :)
      COMPLEX, INTENT(IN) :: orbital_matrices(:, :, :)
      REAL, INTENT(OUT) :: expectation(3)

      COMPLEX :: value
      INTEGER :: direction, orbital, orbital_p, ispin

      CALL validate_square_spin_density(rho, "orbital expectation density")
      IF (SIZE(orbital_matrices, 1) /= SIZE(rho, 1) .OR. &
          SIZE(orbital_matrices, 2) /= SIZE(rho, 1) .OR. SIZE(orbital_matrices, 3) /= 3) THEN
         CALL juDFT_error("Orbital operators do not match local density", calledby="m_local_d_spin_analysis")
      END IF
      expectation = 0.0
      DO direction = 1, 3
         value = CMPLX(0.0, 0.0)
         DO orbital = 1, SIZE(rho, 1)
            DO orbital_p = 1, SIZE(rho, 1)
               DO ispin = 1, local_d_n_spins
                  value = value + rho(orbital, ispin, orbital_p, ispin) &
                                * orbital_matrices(orbital_p, orbital, direction)
               END DO
            END DO
         END DO
         expectation(direction) = REAL(value)
      END DO
   END SUBROUTINE local_d_orbital_expectation

   SUBROUTINE local_d_jeff_unitary(unitary)
      ! Columns contain ideal |j_eff,m_j> states in the row basis
      ! xy-up,xy-down,yz-up,yz-down,xz-up,xz-down. Column order is
      ! 1/2:-1/2,+1/2, then 3/2:-3/2,-1/2,+1/2,+3/2.
      COMPLEX, INTENT(OUT) :: unitary(local_d_n_jeff, local_d_n_jeff)

      REAL :: inverse_sqrt_two, inverse_sqrt_three, inverse_sqrt_six, sqrt_two_thirds

      inverse_sqrt_two = 1.0/SQRT(2.0)
      inverse_sqrt_three = 1.0/SQRT(3.0)
      inverse_sqrt_six = 1.0/SQRT(6.0)
      sqrt_two_thirds = SQRT(2.0/3.0)
      unitary = CMPLX(0.0, 0.0)

      unitary(2, local_d_jeff_half_minus) = CMPLX(inverse_sqrt_three, 0.0)
      unitary(3, local_d_jeff_half_minus) = CMPLX(-inverse_sqrt_three, 0.0)
      unitary(5, local_d_jeff_half_minus) = CMPLX(0.0, inverse_sqrt_three)

      unitary(1, local_d_jeff_half_plus) = CMPLX(-inverse_sqrt_three, 0.0)
      unitary(4, local_d_jeff_half_plus) = CMPLX(-inverse_sqrt_three, 0.0)
      unitary(6, local_d_jeff_half_plus) = CMPLX(0.0, -inverse_sqrt_three)

      unitary(4, local_d_jeff_three_minus_three) = CMPLX(inverse_sqrt_two, 0.0)
      unitary(6, local_d_jeff_three_minus_three) = CMPLX(0.0, -inverse_sqrt_two)

      unitary(2, local_d_jeff_three_minus_one) = CMPLX(sqrt_two_thirds, 0.0)
      unitary(3, local_d_jeff_three_minus_one) = CMPLX(inverse_sqrt_six, 0.0)
      unitary(5, local_d_jeff_three_minus_one) = CMPLX(0.0, -inverse_sqrt_six)

      unitary(1, local_d_jeff_three_plus_one) = CMPLX(sqrt_two_thirds, 0.0)
      unitary(4, local_d_jeff_three_plus_one) = CMPLX(-inverse_sqrt_six, 0.0)
      unitary(6, local_d_jeff_three_plus_one) = CMPLX(0.0, -inverse_sqrt_six)

      unitary(3, local_d_jeff_three_plus_three) = CMPLX(-inverse_sqrt_two, 0.0)
      unitary(5, local_d_jeff_three_plus_three) = CMPLX(0.0, -inverse_sqrt_two)
   END SUBROUTINE local_d_jeff_unitary

   SUBROUTINE local_d_transform_to_jeff(rho_t2g, rho_jeff)
      COMPLEX, INTENT(IN) :: rho_t2g(:, :, :, :)
      COMPLEX, INTENT(OUT) :: rho_jeff(local_d_n_jeff, local_d_n_jeff)

      COMPLEX :: rho_flat(local_d_n_jeff, local_d_n_jeff)
      COMPLEX :: unitary(local_d_n_jeff, local_d_n_jeff)
      INTEGER :: orbital, orbital_p, ispin, jspin, row, column

      CALL validate_d_spin_density(rho_t2g, local_d_n_t2g, "t2g spin density")
      DO orbital = 1, local_d_n_t2g
         DO ispin = 1, local_d_n_spins
            row = (orbital - 1)*local_d_n_spins + ispin
            DO orbital_p = 1, local_d_n_t2g
               DO jspin = 1, local_d_n_spins
                  column = (orbital_p - 1)*local_d_n_spins + jspin
                  rho_flat(row, column) = rho_t2g(orbital, ispin, orbital_p, jspin)
               END DO
            END DO
         END DO
      END DO
      CALL local_d_jeff_unitary(unitary)
      rho_jeff = MATMUL(CONJG(TRANSPOSE(unitary)), MATMUL(rho_flat, unitary))
   END SUBROUTINE local_d_transform_to_jeff

   SUBROUTINE local_d_jeff_weights(rho_t2g, mj_weights, j_weights, normalized_mj_weights, &
                                   normalized_j_weights, normalization_defined)
      ! Raw weights retain the local t2g weight. Optional normalized weights
      ! divide by Tr(rho_t2g) only when that positive trace is numerically
      ! defined; otherwise they are returned as zero and the flag is false.
      COMPLEX, INTENT(IN) :: rho_t2g(:, :, :, :)
      REAL, INTENT(OUT) :: mj_weights(local_d_n_jeff)
      REAL, INTENT(OUT) :: j_weights(2)
      REAL, OPTIONAL, INTENT(OUT) :: normalized_mj_weights(local_d_n_jeff)
      REAL, OPTIONAL, INTENT(OUT) :: normalized_j_weights(2)
      LOGICAL, OPTIONAL, INTENT(OUT) :: normalization_defined

      COMPLEX :: rho_jeff(local_d_n_jeff, local_d_n_jeff)
      REAL :: threshold, t2g_weight
      LOGICAL :: is_defined
      INTEGER :: channel

      CALL local_d_transform_to_jeff(rho_t2g, rho_jeff)
      DO channel = 1, local_d_n_jeff
         mj_weights(channel) = REAL(rho_jeff(channel, channel))
      END DO
      j_weights(1) = SUM(mj_weights(1:2))
      j_weights(2) = SUM(mj_weights(3:6))
      t2g_weight = SUM(mj_weights)
      threshold = 100.0*EPSILON(t2g_weight)*MAX(1.0, ABS(t2g_weight))
      is_defined = t2g_weight > threshold
      IF (PRESENT(normalized_mj_weights)) THEN
         normalized_mj_weights = 0.0
         IF (is_defined) normalized_mj_weights = mj_weights/t2g_weight
      END IF
      IF (PRESENT(normalized_j_weights)) THEN
         normalized_j_weights = 0.0
         IF (is_defined) normalized_j_weights = j_weights/t2g_weight
      END IF
      IF (PRESENT(normalization_defined)) normalization_defined = is_defined
   END SUBROUTINE local_d_jeff_weights

   SUBROUTINE validate_d_spin_density(rho, n_orbitals, description)
      COMPLEX, INTENT(IN) :: rho(:, :, :, :)
      INTEGER, INTENT(IN) :: n_orbitals
      CHARACTER(LEN=*), INTENT(IN) :: description

      IF (SIZE(rho, 1) /= n_orbitals .OR. SIZE(rho, 2) /= local_d_n_spins .OR. &
          SIZE(rho, 3) /= n_orbitals .OR. SIZE(rho, 4) /= local_d_n_spins) THEN
         CALL juDFT_error(TRIM(description)//" has inconsistent dimensions", &
                          calledby="m_local_d_spin_analysis")
      END IF
   END SUBROUTINE validate_d_spin_density

   SUBROUTINE validate_square_spin_density(rho, description)
      COMPLEX, INTENT(IN) :: rho(:, :, :, :)
      CHARACTER(LEN=*), INTENT(IN) :: description

      IF (SIZE(rho, 1) < 1 .OR. SIZE(rho, 1) /= SIZE(rho, 3) .OR. &
          SIZE(rho, 2) /= local_d_n_spins .OR. SIZE(rho, 4) /= local_d_n_spins) THEN
         CALL juDFT_error(TRIM(description)//" has inconsistent dimensions", &
                          calledby="m_local_d_spin_analysis")
      END IF
   END SUBROUTINE validate_square_spin_density

END MODULE m_local_d_spin_analysis
