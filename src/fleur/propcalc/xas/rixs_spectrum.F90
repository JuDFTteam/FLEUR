!--------------------------------------------------------------------------------
! Copyright (c) 2026 Peter Gruenberg Institut, Forschungszentrum Juelich, Germany
! This file is part of FLEUR and available as free software under the conditions
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------

MODULE m_rixs_spectrum
   USE m_juDFT, ONLY: juDFT_error
   USE m_xas_spectrum, ONLY: xas_gaussian_broadening
   IMPLICIT NONE
   PRIVATE

   REAL, PARAMETER :: rixs_occ_tol = 1.0e-10

   PUBLIC :: rixs_accumulate_scalar_spin_trace_spectrum
   PUBLIC :: rixs_scalar_spin_trace_abs2

CONTAINS

   SUBROUTINE rixs_accumulate_scalar_spin_trace_spectrum(loss_grid, eig_band, occ_band, wk, core_energy, omega_in, gamma_core, &
                                                         eta_loss, matrix_abs_spin, matrix_emit_spin, intensity)
      ! In the guarded scalar jspins=1 prototype the bands are spin-degenerate.
      ! The final electron and valence-hole spin labels are orthogonal final
      ! states, so they are traced incoherently:
      !   sum_{sigma_v,sigma_n} |sum_mj M_emit(sigma_v) M_abs(sigma_n)/denom|^2.
      ! This replaces the old one-spin coherent contraction, which is not a
      ! valid scalar L2/L3 RIXS trace.
      REAL,    INTENT(IN)    :: loss_grid(:), eig_band(:), occ_band(:)
      REAL,    INTENT(IN)    :: wk, core_energy, omega_in, gamma_core, eta_loss
      COMPLEX, INTENT(IN)    :: matrix_abs_spin(:, :, :), matrix_emit_spin(:, :, :)
      REAL,    INTENT(INOUT) :: intensity(:)

      INTEGER :: i_grid, n_band, v_band
      REAL    :: f_v, vacancy_n, weight, loss_energy, gaussian, strength
      COMPLEX :: denominator

      CALL rixs_check_spin_trace_inputs(loss_grid, eig_band, occ_band, matrix_abs_spin, matrix_emit_spin, gamma_core, eta_loss, &
                                        intensity)

      DO n_band = 1, SIZE(eig_band)
         vacancy_n = 1.0 - occ_band(n_band)
         IF (vacancy_n <= rixs_occ_tol) CYCLE
         denominator = CMPLX(omega_in - (eig_band(n_band) - core_energy), gamma_core)
         DO v_band = 1, SIZE(eig_band)
            f_v = occ_band(v_band)
            IF (f_v <= rixs_occ_tol) CYCLE

            strength = rixs_scalar_spin_trace_abs2(matrix_abs_spin(n_band, :, :), matrix_emit_spin(v_band, :, :), denominator)
            IF (strength < TINY(strength)) CYCLE
            weight = wk*f_v*vacancy_n
            IF (weight <= 0.0) CYCLE
            loss_energy = eig_band(n_band) - eig_band(v_band)
            DO i_grid = 1, SIZE(loss_grid)
               gaussian = xas_gaussian_broadening(loss_grid(i_grid) - loss_energy, eta_loss)
               IF (gaussian == 0.0) CYCLE
               intensity(i_grid) = intensity(i_grid) + weight*strength*gaussian
            END DO
         END DO
      END DO
   END SUBROUTINE rixs_accumulate_scalar_spin_trace_spectrum

   REAL FUNCTION rixs_scalar_spin_trace_abs2(abs_by_mj_spin, emit_by_mj_spin, denominator) RESULT(strength)
      COMPLEX, INTENT(IN) :: abs_by_mj_spin(:, :), emit_by_mj_spin(:, :)
      COMPLEX, INTENT(IN) :: denominator

      INTEGER :: sigma_n, sigma_v, n_spin
      COMPLEX :: kh_amp

      IF (SIZE(abs_by_mj_spin, 1) /= SIZE(emit_by_mj_spin, 1)) THEN
         CALL juDFT_error("RIXS spin-trace mj dimensions differ", calledby="m_rixs_spectrum")
      END IF
      IF (SIZE(abs_by_mj_spin, 2) /= SIZE(emit_by_mj_spin, 2) .OR. SIZE(abs_by_mj_spin, 2) < 1) THEN
         CALL juDFT_error("RIXS spin-trace spin dimensions differ", calledby="m_rixs_spectrum")
      END IF

      strength = 0.0
      n_spin = SIZE(abs_by_mj_spin, 2)
      DO sigma_v = 1, n_spin
         DO sigma_n = 1, n_spin
            kh_amp = SUM(emit_by_mj_spin(:, sigma_v)*abs_by_mj_spin(:, sigma_n))/denominator
            strength = strength + ABS(kh_amp)**2
         END DO
      END DO
   END FUNCTION rixs_scalar_spin_trace_abs2

   SUBROUTINE rixs_check_spin_trace_inputs(loss_grid, eig_band, occ_band, matrix_abs_spin, matrix_emit_spin, gamma_core, &
                                           eta_loss, intensity)
      REAL,    INTENT(IN) :: loss_grid(:), eig_band(:), occ_band(:)
      COMPLEX, INTENT(IN) :: matrix_abs_spin(:, :, :), matrix_emit_spin(:, :, :)
      REAL,    INTENT(IN) :: gamma_core, eta_loss, intensity(:)

      IF (gamma_core <= 0.0) CALL juDFT_error("gammaCore must be positive in RIXS spin trace", calledby="m_rixs_spectrum")
      IF (eta_loss <= 0.0) CALL juDFT_error("etaLoss must be positive in RIXS spin trace", calledby="m_rixs_spectrum")
      IF (SIZE(loss_grid) /= SIZE(intensity)) THEN
         CALL juDFT_error("loss grid and intensity sizes differ in RIXS spin trace", calledby="m_rixs_spectrum")
      END IF
      IF (SIZE(eig_band) /= SIZE(occ_band)) THEN
         CALL juDFT_error("RIXS spin-trace eig_band and occ_band sizes differ", calledby="m_rixs_spectrum")
      END IF
      IF (SIZE(matrix_abs_spin, 1) < SIZE(eig_band) .OR. SIZE(matrix_emit_spin, 1) < SIZE(eig_band)) THEN
         CALL juDFT_error("RIXS spin-trace matrix band dimensions are too small", calledby="m_rixs_spectrum")
      END IF
      IF (SIZE(matrix_abs_spin, 2) /= SIZE(matrix_emit_spin, 2)) THEN
         CALL juDFT_error("RIXS spin-trace absorption/emission mj dimensions differ", calledby="m_rixs_spectrum")
      END IF
      IF (SIZE(matrix_abs_spin, 3) /= SIZE(matrix_emit_spin, 3) .OR. SIZE(matrix_abs_spin, 3) < 1) THEN
         CALL juDFT_error("RIXS spin-trace absorption/emission spin dimensions differ", calledby="m_rixs_spectrum")
      END IF
   END SUBROUTINE rixs_check_spin_trace_inputs

END MODULE m_rixs_spectrum
