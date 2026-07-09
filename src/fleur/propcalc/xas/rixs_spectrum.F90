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

   PUBLIC :: rixs_accumulate_direct_spectrum

CONTAINS

   SUBROUTINE rixs_accumulate_direct_spectrum(loss_grid, eig_band, occ_band, wk, core_energy, omega_in, gamma_core, &
                                              eta_loss, matrix_abs, matrix_emit, intensity)
      ! Direct same-k independent-particle RIXS prototype:
      ! I(Omega) += sum_vn wk*f_v*(1-f_n)
      !             * |sum_mj M_emit(v,mj) M_abs(n,mj)
      !                /(omega_in-(eig_n-core_energy)+i Gamma)|^2
      !             * g_eta_loss(Omega-(eig_n-eig_v)).
      REAL,    INTENT(IN)    :: loss_grid(:), eig_band(:), occ_band(:)
      REAL,    INTENT(IN)    :: wk, core_energy, omega_in, gamma_core, eta_loss
      COMPLEX, INTENT(IN)    :: matrix_abs(:, :), matrix_emit(:, :)
      REAL,    INTENT(INOUT) :: intensity(:)

      INTEGER :: i_grid, n_band, v_band
      REAL    :: f_v, vacancy_n, weight, loss_energy, gaussian, strength
      COMPLEX :: denominator, kh_amp

      CALL rixs_check_spectrum_inputs(loss_grid, eig_band, occ_band, matrix_abs, matrix_emit, gamma_core, eta_loss, intensity)

      DO n_band = 1, SIZE(eig_band)
         vacancy_n = 1.0 - occ_band(n_band)
         IF (vacancy_n <= rixs_occ_tol) CYCLE
         denominator = CMPLX(omega_in - (eig_band(n_band) - core_energy), gamma_core)
         DO v_band = 1, SIZE(eig_band)
            f_v = occ_band(v_band)
            IF (f_v <= rixs_occ_tol) CYCLE
            kh_amp = SUM(matrix_emit(v_band, :)*matrix_abs(n_band, :))/denominator
            strength = ABS(kh_amp)**2
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
   END SUBROUTINE rixs_accumulate_direct_spectrum

   SUBROUTINE rixs_check_spectrum_inputs(loss_grid, eig_band, occ_band, matrix_abs, matrix_emit, gamma_core, eta_loss, intensity)
      REAL,    INTENT(IN) :: loss_grid(:), eig_band(:), occ_band(:)
      COMPLEX, INTENT(IN) :: matrix_abs(:, :), matrix_emit(:, :)
      REAL,    INTENT(IN) :: gamma_core, eta_loss, intensity(:)

      IF (gamma_core <= 0.0) CALL juDFT_error("gammaCore must be positive in RIXS spectrum", calledby="m_rixs_spectrum")
      IF (eta_loss <= 0.0) CALL juDFT_error("etaLoss must be positive in RIXS spectrum", calledby="m_rixs_spectrum")
      IF (SIZE(loss_grid) /= SIZE(intensity)) THEN
         CALL juDFT_error("loss grid and intensity sizes differ in RIXS spectrum", calledby="m_rixs_spectrum")
      END IF
      IF (SIZE(eig_band) /= SIZE(occ_band)) THEN
         CALL juDFT_error("RIXS eig_band and occ_band sizes differ", calledby="m_rixs_spectrum")
      END IF
      IF (SIZE(matrix_abs, 1) < SIZE(eig_band) .OR. SIZE(matrix_emit, 1) < SIZE(eig_band)) THEN
         CALL juDFT_error("RIXS matrix band dimensions are too small", calledby="m_rixs_spectrum")
      END IF
      IF (SIZE(matrix_abs, 2) /= SIZE(matrix_emit, 2)) THEN
         CALL juDFT_error("RIXS absorption/emission mj dimensions differ", calledby="m_rixs_spectrum")
      END IF
   END SUBROUTINE rixs_check_spectrum_inputs

END MODULE m_rixs_spectrum
