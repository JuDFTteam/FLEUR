!--------------------------------------------------------------------------------
! Copyright (c) 2026 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions 
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------
!--------------------------------------------------------------------------------

PROGRAM rixs_finiteq_test
   USE m_constants, ONLY: tpi_const
   USE m_rixs_momentum, ONLY: rixs_fold_rlu, rixs_phased_site_amplitude, rixs_site_phase
   USE m_rixs_spectrum, ONLY: rixs_add_finiteq_spinor_site_amplitudes, rixs_spinor_amplitude
   USE m_xas_amplitudes, ONLY: xas_emission_from_absorption
   IMPLICIT NONE

   CALL test_full_q_site_interference()
   CALL test_bloch_site_phase_invariances()
   CALL test_absorption_emission_conjugation()
   CALL test_one_site_qzero_kernel_limit()
   WRITE(*, '(a)') "finite-Q RIXS focused tests: PASS"

CONTAINS

   SUBROUTINE assert_close_complex(actual, expected, tolerance, message)
      COMPLEX, INTENT(IN) :: actual, expected
      REAL, INTENT(IN) :: tolerance
      CHARACTER(LEN=*), INTENT(IN) :: message
      IF (ABS(actual - expected) > tolerance) THEN
         WRITE(*, '(a,4es24.16)') TRIM(message)//": ", REAL(actual), AIMAG(actual), REAL(expected), AIMAG(expected)
         ERROR STOP 1
      END IF
   END SUBROUTINE assert_close_complex

   SUBROUTINE assert_close_real(actual, expected, tolerance, message)
      REAL, INTENT(IN) :: actual, expected, tolerance
      CHARACTER(LEN=*), INTENT(IN) :: message
      IF (ABS(actual - expected) > tolerance) THEN
         WRITE(*, '(a,2es24.16)') TRIM(message)//": ", actual, expected
         ERROR STOP 1
      END IF
   END SUBROUTINE assert_close_real

   SUBROUTINE test_full_q_site_interference()
      REAL :: q_full(3), q_reduced(3), tau(3, 4)
      COMPLEX :: sum_full, sum_reduced
      INTEGER :: i

      q_full = [1.0, 0.0, 1.0]
      q_reduced = rixs_fold_rlu(q_full)
      tau(:, 1) = [0.0, 0.0, 0.0]
      tau(:, 2) = [0.0, 0.0, 0.5]
      tau(:, 3) = [0.5, 0.5, 0.0]
      tau(:, 4) = [0.5, 0.5, 0.5]
      sum_full = CMPLX(0.0, 0.0)
      sum_reduced = CMPLX(0.0, 0.0)
      DO i = 1, 4
         CALL assert_close_real(ABS(rixs_site_phase(q_full, tau(:, i))), 1.0, 2.0e-6, "site phase magnitude")
         sum_full = sum_full + rixs_site_phase(q_full, tau(:, i))
         sum_reduced = sum_reduced + rixs_site_phase(q_reduced, tau(:, i))
      END DO
      CALL assert_close_complex(sum_full, CMPLX(0.0, 0.0), 5.0e-6, "full-Q four-Ir destructive sum")
      CALL assert_close_complex(sum_reduced, CMPLX(4.0, 0.0), 5.0e-6, "reduced-q four-Ir constructive sum")
   END SUBROUTINE test_full_q_site_interference

   SUBROUTINE test_bloch_site_phase_invariances()
      REAL :: q_full(3), k_v(3), k_n(3), reciprocal_shift(3), tau(3), lattice_translation(3), origin_shift(3)
      COMPLEX :: matrix_abs, matrix_emit, original, translated, shifted, expected_origin_phase
      COMPLEX :: missing_phase_translated, double_phase_original, double_phase_translated

      q_full = [0.5, -0.5, 1.0]
      k_v = [0.75, 0.25, 0.0]
      k_n = rixs_fold_rlu(k_v + q_full)
      reciprocal_shift = ANINT(k_v + q_full - k_n)
      tau = [0.13, 0.27, 0.41]
      lattice_translation = [1.0, 0.0, 0.0]
      origin_shift = [0.17, -0.11, 0.23]
      matrix_abs = CMPLX(0.7, -0.2)
      matrix_emit = CMPLX(-0.4, 0.9)

      original = rixs_phased_site_amplitude(q_full, tau, matrix_emit*matrix_abs)
      translated = rixs_phased_site_amplitude(q_full, tau + lattice_translation, &
         matrix_emit*EXP(CMPLX(0.0, tpi_const*DOT_PRODUCT(k_v, lattice_translation))) * &
         matrix_abs*EXP(CMPLX(0.0, -tpi_const*DOT_PRODUCT(k_n, lattice_translation))))
      CALL assert_close_complex(translated, original, 5.0e-6, "unit-cell representative invariance")

      shifted = rixs_phased_site_amplitude(q_full, tau + origin_shift, &
         matrix_emit*EXP(CMPLX(0.0, tpi_const*DOT_PRODUCT(k_v, origin_shift))) * &
         matrix_abs*EXP(CMPLX(0.0, -tpi_const*DOT_PRODUCT(k_n, origin_shift))))
      expected_origin_phase = EXP(CMPLX(0.0, tpi_const*DOT_PRODUCT(reciprocal_shift, origin_shift)))
      CALL assert_close_complex(shifted, expected_origin_phase*original, 5.0e-6, "overall-origin amplitude phase")
      CALL assert_close_real(ABS(shifted)**2, ABS(original)**2, 5.0e-6, "overall-origin intensity invariance")

      missing_phase_translated = matrix_emit*EXP(CMPLX(0.0, tpi_const*DOT_PRODUCT(k_v, lattice_translation))) * &
                                 matrix_abs*EXP(CMPLX(0.0, -tpi_const*DOT_PRODUCT(k_n, lattice_translation)))
      IF (ABS(missing_phase_translated - matrix_emit*matrix_abs) < 0.1) THEN
         ERROR STOP "representative test is insensitive to a missing photon phase"
      END IF
      double_phase_original = rixs_site_phase(q_full, tau)**2*matrix_emit*matrix_abs
      double_phase_translated = rixs_site_phase(q_full, tau + lattice_translation)**2 * &
         matrix_emit*EXP(CMPLX(0.0, tpi_const*DOT_PRODUCT(k_v, lattice_translation))) * &
         matrix_abs*EXP(CMPLX(0.0, -tpi_const*DOT_PRODUCT(k_n, lattice_translation)))
      IF (ABS(double_phase_translated - double_phase_original) < 0.1) THEN
         ERROR STOP "representative test is insensitive to a double-counted photon phase"
      END IF
   END SUBROUTINE test_bloch_site_phase_invariances

   SUBROUTINE test_absorption_emission_conjugation()
      COMPLEX :: absorption(2, 2), emission(2, 2), double_conjugated(2, 2)

      ! A genuinely complex synthetic result, as generated for a complex
      ! polarization.  Both production emission wrappers call this same helper.
      absorption = RESHAPE([CMPLX(1.0, 2.0), CMPLX(-0.5, 0.75), CMPLX(0.3, -1.1), CMPLX(-2.0, -0.4)], [2, 2])
      emission = xas_emission_from_absorption(absorption)
      IF (MAXVAL(ABS(emission - CONJG(absorption))) > 1.0e-7) THEN
         ERROR STOP "emission is not the Hermitian conjugate of absorption"
      END IF
      double_conjugated = xas_emission_from_absorption(CONJG(absorption))
      IF (MAXVAL(ABS(double_conjugated - emission)) < 0.1) THEN
         ERROR STOP "complex-polarization test cannot detect caller-side double conjugation"
      END IF
   END SUBROUTINE test_absorption_emission_conjugation

   SUBROUTINE test_one_site_qzero_kernel_limit()
      REAL :: eig_n(2)
      COMPLEX :: matrix_abs(2, 2), matrix_emit(2, 2), amplitude(2, 2), denominator, expected
      INTEGER :: band_v, band_n

      eig_n = [0.4, 0.7]
      matrix_abs = RESHAPE([CMPLX(0.2, 0.1), CMPLX(-0.3, 0.8), CMPLX(0.7, -0.4), CMPLX(0.1, 0.6)], [2, 2])
      matrix_emit = RESHAPE([CMPLX(-0.5, 0.2), CMPLX(0.9, 0.3), CMPLX(0.4, -0.7), CMPLX(-0.2, 0.1)], [2, 2])
      amplitude = CMPLX(0.0, 0.0)
      CALL rixs_add_finiteq_spinor_site_amplitudes(eig_n, -1.0, 1.8, 0.12, matrix_abs, matrix_emit, &
         rixs_site_phase([0.0, 0.0, 0.0], [0.37, 0.21, 0.44]), 1, 2, 1, 2, amplitude)
      DO band_v = 1, 2
         DO band_n = 1, 2
            denominator = CMPLX(1.8 - (eig_n(band_n) + 1.0), 0.12)
            expected = rixs_spinor_amplitude(matrix_emit(band_v, :), matrix_abs(band_n, :), denominator)
            CALL assert_close_complex(amplitude(band_v, band_n), expected, 2.0e-6, "one-site Q=0 kernel limit")
         END DO
      END DO
   END SUBROUTINE test_one_site_qzero_kernel_limit

END PROGRAM rixs_finiteq_test
