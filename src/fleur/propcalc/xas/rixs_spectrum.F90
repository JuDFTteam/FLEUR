!--------------------------------------------------------------------------------
! Copyright (c) 2026 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------

MODULE m_rixs_spectrum
   USE m_juDFT, ONLY: juDFT_error
   USE m_xas_spectrum, ONLY: xas_gaussian_broadening
   IMPLICIT NONE
   PRIVATE

   REAL, PARAMETER :: rixs_occupation_tolerance = 1.0e-10

   PUBLIC :: rixs_occupation_tolerance
   PUBLIC :: rixs_accumulate_scalar_spin_trace_spectrum
   PUBLIC :: rixs_scalar_spin_trace_abs2
   PUBLIC :: rixs_accumulate_spinor_spectrum
   PUBLIC :: rixs_spinor_amplitude
   PUBLIC :: rixs_add_finiteq_spinor_site_amplitudes
   PUBLIC :: rixs_accumulate_finiteq_spinor_spectrum

CONTAINS

   SUBROUTINE rixs_accumulate_scalar_spin_trace_spectrum(loss_grid, eig_band, occ_band, wk, core_energy, omega_in, gamma_core, &
                                                         eta_loss, matrix_abs_spin, matrix_emit_spin, valence_band_min, &
                                                         valence_band_max, intermediate_band_min, intermediate_band_max, &
                                                         intensity)
      ! In the scalar jspins=1 branch the bands are spin-degenerate.
      ! The final electron and valence-hole spin labels are orthogonal final
      ! states, so they are traced incoherently:
      !   sum_{sigma_v,sigma_n} |sum_mj M_emit(sigma_v) M_abs(sigma_n)/denom|^2.
      ! This replaces the old one-spin coherent contraction, which is not a
      ! valid scalar L2/L3 RIXS trace.
      REAL,    INTENT(IN)    :: loss_grid(:), eig_band(:), occ_band(:)
      REAL,    INTENT(IN)    :: wk, core_energy, omega_in, gamma_core, eta_loss
      COMPLEX, INTENT(IN)    :: matrix_abs_spin(:, :, :), matrix_emit_spin(:, :, :)
      INTEGER, INTENT(IN)    :: valence_band_min, valence_band_max, intermediate_band_min, intermediate_band_max
      REAL,    INTENT(INOUT) :: intensity(:)

      INTEGER :: i_grid, n_band, v_band
      REAL    :: f_v, vacancy_n, weight, loss_energy, gaussian, strength
      COMPLEX :: denominator

      CALL rixs_check_spin_trace_inputs(loss_grid, eig_band, occ_band, matrix_abs_spin, matrix_emit_spin, gamma_core, eta_loss, &
                                        valence_band_min, valence_band_max, intermediate_band_min, intermediate_band_max, &
                                        intensity)

      DO n_band = intermediate_band_min, intermediate_band_max
         vacancy_n = 1.0 - occ_band(n_band)
         IF (vacancy_n <= rixs_occupation_tolerance) CYCLE
         denominator = CMPLX(omega_in - (eig_band(n_band) - core_energy), gamma_core)
         DO v_band = valence_band_min, valence_band_max
            f_v = occ_band(v_band)
            IF (f_v <= rixs_occupation_tolerance) CYCLE

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

   SUBROUTINE rixs_accumulate_spinor_spectrum(loss_grid, eig_band, occ_band, wk, core_energy, omega_in, gamma_core, &
                                              eta_loss, matrix_abs, matrix_emit, valence_band_min, valence_band_max, &
                                              intermediate_band_min, intermediate_band_max, intensity)
      ! First-variation noco bands already contain a coherent two-component
      ! spinor. The spin components are contracted inside the XAS matrix-element
      ! routines, so RIXS forms one coherent core-mj amplitude per (v,n,k,atom)
      ! and must not apply the scalar spin-degenerate S1 trace.
      REAL,    INTENT(IN)    :: loss_grid(:), eig_band(:), occ_band(:)
      REAL,    INTENT(IN)    :: wk, core_energy, omega_in, gamma_core, eta_loss
      COMPLEX, INTENT(IN)    :: matrix_abs(:, :), matrix_emit(:, :)
      INTEGER, INTENT(IN)    :: valence_band_min, valence_band_max, intermediate_band_min, intermediate_band_max
      REAL,    INTENT(INOUT) :: intensity(:)

      INTEGER :: i_grid, n_band, v_band
      REAL    :: f_v, vacancy_n, weight, loss_energy, gaussian, strength
      COMPLEX :: amplitude, denominator

      CALL rixs_check_spinor_inputs(loss_grid, eig_band, occ_band, matrix_abs, matrix_emit, gamma_core, eta_loss, &
                                    valence_band_min, valence_band_max, intermediate_band_min, intermediate_band_max, &
                                    intensity)

      DO v_band = valence_band_min, valence_band_max
         f_v = occ_band(v_band)
         IF (f_v <= rixs_occupation_tolerance) CYCLE
         DO n_band = intermediate_band_min, intermediate_band_max
            vacancy_n = 1.0 - occ_band(n_band)
            IF (vacancy_n <= rixs_occupation_tolerance) CYCLE

            denominator = CMPLX(omega_in - (eig_band(n_band) - core_energy), gamma_core)
            amplitude = rixs_spinor_amplitude(matrix_emit(v_band, :), matrix_abs(n_band, :), denominator)
            strength = ABS(amplitude)**2
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
   END SUBROUTINE rixs_accumulate_spinor_spectrum

   COMPLEX FUNCTION rixs_spinor_amplitude(emit_by_mj, abs_by_mj, denominator) RESULT(amplitude)
      COMPLEX, INTENT(IN) :: emit_by_mj(:), abs_by_mj(:), denominator

      amplitude = SUM(emit_by_mj*abs_by_mj)/denominator
   END FUNCTION rixs_spinor_amplitude

   SUBROUTINE rixs_add_finiteq_spinor_site_amplitudes(eig_intermediate, core_energy, omega_in, gamma_core, &
                                                       matrix_abs, matrix_emit, site_phase, valence_band_min, &
                                                       valence_band_max, intermediate_band_min, intermediate_band_max, &
                                                       amplitude_vn, site_partial_vn)
      ! Add one absorber site's complex contribution. site_phase must contain
      ! exp(+i 2*pi Q_full_rlu.tau_fractional), not a phase made from reduced q.
      ! No absolute square is taken here: all basis sites/types must be summed
      ! before the spectrum is accumulated.
      REAL,    INTENT(IN)    :: eig_intermediate(:), core_energy, omega_in, gamma_core
      COMPLEX, INTENT(IN)    :: matrix_abs(:, :), matrix_emit(:, :), site_phase
      INTEGER, INTENT(IN)    :: valence_band_min, valence_band_max, intermediate_band_min, intermediate_band_max
      COMPLEX, INTENT(INOUT) :: amplitude_vn(:, :)
      COMPLEX, OPTIONAL, INTENT(OUT) :: site_partial_vn(:, :)

      INTEGER :: band_v, band_n
      COMPLEX :: denominator, partial

      IF (gamma_core <= 0.0) CALL juDFT_error("gammaCore must be positive in finite-Q RIXS", calledby="m_rixs_spectrum")
      IF (SIZE(matrix_abs, 1) < SIZE(eig_intermediate) .OR. &
          SIZE(matrix_abs, 2) /= SIZE(matrix_emit, 2)) THEN
         CALL juDFT_error("Finite-Q RIXS site matrix dimensions are inconsistent", calledby="m_rixs_spectrum")
      END IF
      IF (SIZE(amplitude_vn, 1) < SIZE(matrix_emit, 1) .OR. &
          SIZE(amplitude_vn, 2) < SIZE(matrix_abs, 1)) THEN
         CALL juDFT_error("Finite-Q RIXS amplitude array is too small", calledby="m_rixs_spectrum")
      END IF
      IF (PRESENT(site_partial_vn)) THEN
         IF (ANY(SHAPE(site_partial_vn) /= SHAPE(amplitude_vn))) THEN
            CALL juDFT_error("Finite-Q RIXS site-partial array shape differs from amplitude array", &
                             calledby="m_rixs_spectrum")
         END IF
         site_partial_vn = CMPLX(0.0, 0.0)
      END IF

      DO band_n = intermediate_band_min, intermediate_band_max
         denominator = CMPLX(omega_in - (eig_intermediate(band_n) - core_energy), gamma_core)
         DO band_v = valence_band_min, valence_band_max
            partial = site_phase*rixs_spinor_amplitude(matrix_emit(band_v, :), matrix_abs(band_n, :), denominator)
            amplitude_vn(band_v, band_n) = amplitude_vn(band_v, band_n) + partial
            IF (PRESENT(site_partial_vn)) site_partial_vn(band_v, band_n) = partial
         END DO
      END DO
   END SUBROUTINE rixs_add_finiteq_spinor_site_amplitudes

   SUBROUTINE rixs_accumulate_finiteq_spinor_spectrum(loss_grid, eig_valence, occ_valence, eig_intermediate, &
                                                       occ_intermediate, wk, eta_loss, amplitude_vn, valence_band_min, &
                                                       valence_band_max, intermediate_band_min, intermediate_band_max, &
                                                       intensity)
      REAL,    INTENT(IN)    :: loss_grid(:), eig_valence(:), occ_valence(:), eig_intermediate(:), occ_intermediate(:)
      REAL,    INTENT(IN)    :: wk, eta_loss
      COMPLEX, INTENT(IN)    :: amplitude_vn(:, :)
      INTEGER, INTENT(IN)    :: valence_band_min, valence_band_max, intermediate_band_min, intermediate_band_max
      REAL,    INTENT(INOUT) :: intensity(:)

      INTEGER :: band_v, band_n, i_grid
      REAL :: f_v, vacancy_n, strength, weight, loss_energy, gaussian

      IF (eta_loss <= 0.0) CALL juDFT_error("etaLoss must be positive in finite-Q RIXS", calledby="m_rixs_spectrum")
      IF (SIZE(loss_grid) /= SIZE(intensity)) THEN
         CALL juDFT_error("Finite-Q RIXS loss grid and intensity sizes differ", calledby="m_rixs_spectrum")
      END IF
      IF (SIZE(eig_valence) /= SIZE(occ_valence) .OR. SIZE(eig_intermediate) /= SIZE(occ_intermediate)) THEN
         CALL juDFT_error("Finite-Q RIXS eigenvalue and occupation sizes differ", calledby="m_rixs_spectrum")
      END IF
      IF (SIZE(amplitude_vn, 1) < SIZE(eig_valence) .OR. SIZE(amplitude_vn, 2) < SIZE(eig_intermediate)) THEN
         CALL juDFT_error("Finite-Q RIXS amplitude dimensions are too small", calledby="m_rixs_spectrum")
      END IF

      DO band_v = valence_band_min, valence_band_max
         f_v = occ_valence(band_v)
         IF (f_v <= rixs_occupation_tolerance) CYCLE
         DO band_n = intermediate_band_min, intermediate_band_max
            vacancy_n = 1.0 - occ_intermediate(band_n)
            IF (vacancy_n <= rixs_occupation_tolerance) CYCLE
            strength = ABS(amplitude_vn(band_v, band_n))**2
            IF (strength < TINY(strength)) CYCLE
            weight = wk*f_v*vacancy_n
            IF (weight <= 0.0) CYCLE
            loss_energy = eig_intermediate(band_n) - eig_valence(band_v)
            DO i_grid = 1, SIZE(loss_grid)
               gaussian = xas_gaussian_broadening(loss_grid(i_grid) - loss_energy, eta_loss)
               IF (gaussian == 0.0) CYCLE
               intensity(i_grid) = intensity(i_grid) + weight*strength*gaussian
            END DO
         END DO
      END DO
   END SUBROUTINE rixs_accumulate_finiteq_spinor_spectrum

   SUBROUTINE rixs_check_spinor_inputs(loss_grid, eig_band, occ_band, matrix_abs, matrix_emit, gamma_core, eta_loss, &
                                       valence_band_min, valence_band_max, intermediate_band_min, intermediate_band_max, &
                                       intensity)
      REAL,    INTENT(IN) :: loss_grid(:), eig_band(:), occ_band(:)
      COMPLEX, INTENT(IN) :: matrix_abs(:, :), matrix_emit(:, :)
      REAL,    INTENT(IN) :: gamma_core, eta_loss, intensity(:)
      INTEGER, INTENT(IN) :: valence_band_min, valence_band_max, intermediate_band_min, intermediate_band_max

      IF (gamma_core <= 0.0) CALL juDFT_error("gammaCore must be positive in spinor RIXS", calledby="m_rixs_spectrum")
      IF (eta_loss <= 0.0) CALL juDFT_error("etaLoss must be positive in spinor RIXS", calledby="m_rixs_spectrum")
      IF (SIZE(loss_grid) /= SIZE(intensity)) THEN
         CALL juDFT_error("loss grid and intensity sizes differ in spinor RIXS", calledby="m_rixs_spectrum")
      END IF
      IF (SIZE(eig_band) /= SIZE(occ_band)) THEN
         CALL juDFT_error("spinor RIXS eig_band and occ_band sizes differ", calledby="m_rixs_spectrum")
      END IF
      IF (SIZE(matrix_abs, 1) < SIZE(eig_band) .OR. SIZE(matrix_emit, 1) < SIZE(eig_band)) THEN
         CALL juDFT_error("spinor RIXS matrix band dimensions are too small", calledby="m_rixs_spectrum")
      END IF
      IF (SIZE(matrix_abs, 2) /= SIZE(matrix_emit, 2) .OR. SIZE(matrix_abs, 2) < 1) THEN
         CALL juDFT_error("spinor RIXS absorption/emission mj dimensions differ", calledby="m_rixs_spectrum")
      END IF
      IF (valence_band_min < 1 .OR. valence_band_max > SIZE(eig_band) .OR. valence_band_min > valence_band_max) THEN
         CALL juDFT_error("spinor RIXS valence band loop bounds are invalid", calledby="m_rixs_spectrum")
      END IF
      IF (intermediate_band_min < 1 .OR. intermediate_band_max > SIZE(eig_band) .OR. &
          intermediate_band_min > intermediate_band_max) THEN
         CALL juDFT_error("spinor RIXS intermediate band loop bounds are invalid", calledby="m_rixs_spectrum")
      END IF
   END SUBROUTINE rixs_check_spinor_inputs

   SUBROUTINE rixs_check_spin_trace_inputs(loss_grid, eig_band, occ_band, matrix_abs_spin, matrix_emit_spin, gamma_core, &
                                           eta_loss, valence_band_min, valence_band_max, intermediate_band_min, &
                                           intermediate_band_max, intensity)
      REAL,    INTENT(IN) :: loss_grid(:), eig_band(:), occ_band(:)
      COMPLEX, INTENT(IN) :: matrix_abs_spin(:, :, :), matrix_emit_spin(:, :, :)
      REAL,    INTENT(IN) :: gamma_core, eta_loss, intensity(:)
      INTEGER, INTENT(IN) :: valence_band_min, valence_band_max, intermediate_band_min, intermediate_band_max

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
      IF (valence_band_min < 1 .OR. valence_band_max > SIZE(eig_band) .OR. valence_band_min > valence_band_max) THEN
         CALL juDFT_error("RIXS valence band loop bounds are invalid", calledby="m_rixs_spectrum")
      END IF
      IF (intermediate_band_min < 1 .OR. intermediate_band_max > SIZE(eig_band) .OR. &
          intermediate_band_min > intermediate_band_max) THEN
         CALL juDFT_error("RIXS intermediate band loop bounds are invalid", calledby="m_rixs_spectrum")
      END IF
   END SUBROUTINE rixs_check_spin_trace_inputs

END MODULE m_rixs_spectrum
