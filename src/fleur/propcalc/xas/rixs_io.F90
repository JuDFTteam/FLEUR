!--------------------------------------------------------------------------------
! Copyright (c) 2026 Peter Gruenberg Institut, Forschungszentrum Juelich, Germany
! This file is part of FLEUR and available as free software under the conditions
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------

MODULE m_rixs_io
   USE m_constants, ONLY: hartree_to_ev_const
   USE m_juDFT, ONLY: juDFT_error
   USE m_rixs_spectrum, ONLY: rixs_scalar_spin_trace_abs2
   USE m_types_xas, ONLY: t_xas
   USE m_xas_spectrum, ONLY: xas_gaussian_broadening
   IMPLICIT NONE
   PRIVATE

   PUBLIC :: rixs_write_spectrum_text
   PUBLIC :: rixs_open_contribution_table, rixs_write_contribution_rows, rixs_close_contribution_table
   PUBLIC :: rixs_print_setup_summary
   PUBLIC :: rixs_print_pair_summary
   PUBLIC :: rixs_print_contribution_check
   PUBLIC :: rixs_polarization_string
   PUBLIC :: rixs_energy_label

CONTAINS

   SUBROUTINE rixs_write_spectrum_text(filename, loss_grid, intensity)
      CHARACTER(LEN=*), INTENT(IN) :: filename
      REAL,             INTENT(IN) :: loss_grid(:), intensity(:)

      INTEGER :: io_unit, i_grid, io_status
      CHARACTER(LEN=256) :: io_message

      IF (SIZE(loss_grid) /= SIZE(intensity)) THEN
         CALL juDFT_error("loss_grid and intensity sizes differ in rixs_write_spectrum_text", calledby="m_rixs_io")
      END IF

      OPEN(NEWUNIT=io_unit, FILE=TRIM(filename), STATUS="replace", ACTION="write", FORM="formatted", &
           IOSTAT=io_status, IOMSG=io_message)
      IF (io_status /= 0) CALL juDFT_error("Cannot open RIXS spectrum text file: "//TRIM(io_message), calledby="m_rixs_io")

      WRITE(io_unit, '(a)', IOSTAT=io_status, IOMSG=io_message) "# FLEUR independent-particle RIXS spectrum"
      IF (io_status == 0) WRITE(io_unit, '(a)', IOSTAT=io_status, IOMSG=io_message) "# columns: loss_energy_Ha loss_energy_eV intensity"
      DO i_grid = 1, SIZE(loss_grid)
         IF (io_status /= 0) EXIT
         WRITE(io_unit, '(3es24.16e3)', IOSTAT=io_status, IOMSG=io_message) loss_grid(i_grid), &
            loss_grid(i_grid)*hartree_to_ev_const, intensity(i_grid)
      END DO
      CLOSE(io_unit)

      IF (io_status /= 0) CALL juDFT_error("Cannot write RIXS spectrum text file: "//TRIM(io_message), calledby="m_rixs_io")
   END SUBROUTINE rixs_write_spectrum_text

   SUBROUTINE rixs_open_contribution_table(filename, edge, absorber_z, pol_in, pol_out, omega_in, mpi_rank, io_unit)
      CHARACTER(LEN=*), INTENT(IN) :: filename, edge, pol_in, pol_out
      INTEGER,          INTENT(IN) :: absorber_z, mpi_rank
      REAL,             INTENT(IN) :: omega_in
      INTEGER,          INTENT(OUT) :: io_unit

      INTEGER :: io_status
      CHARACTER(LEN=256) :: io_message

      OPEN(NEWUNIT=io_unit, FILE=TRIM(filename), STATUS="replace", ACTION="write", FORM="formatted", &
           IOSTAT=io_status, IOMSG=io_message)
      IF (io_status /= 0) CALL juDFT_error("Cannot open RIXS contribution table: "//TRIM(io_message), calledby="m_rixs_io")

      WRITE(io_unit, '(a)') "# FLEUR independent-particle RIXS contribution table"
      WRITE(io_unit, '(a,a)') "# edge = ", TRIM(edge)
      WRITE(io_unit, '(a,i0)') "# absorberZ = ", absorber_z
      WRITE(io_unit, '(a,a)') "# incoming polarization = ", TRIM(pol_in)
      WRITE(io_unit, '(a,a)') "# outgoing polarization = ", TRIM(pol_out)
      WRITE(io_unit, '(a,es24.16e3,a,es24.16e3,a)') "# omegaIn = ", omega_in, " Ha (", &
                                                     omega_in*hartree_to_ev_const, " eV)"
      WRITE(io_unit, '(a,i0)') "# MPI rank = ", mpi_rank
      WRITE(io_unit, '(a)') "# This file contains only contributions evaluated by this MPI rank."
      WRITE(io_unit, '(a)') "# Scalar RIXS uses the spin-degenerate incoherent trace over final electron/hole spin labels."
      WRITE(io_unit, '(a)') "# amplitude_abs2 is sum_{sigma_v,sigma_n} |A_{sigma_v sigma_n}|^2."
      WRITE(io_unit, '(a)') "# amplitude_real and amplitude_imag are zero placeholders for scalar S1 production."
      WRITE(io_unit, '(a)') "# No single coherent complex amplitude represents the spin-traced contribution."
      WRITE(io_unit, '(a)') "# weighted_strength = k_weight * f_v * (1 - f_n) * amplitude_abs2."
      WRITE(io_unit, '(a)') "# columns:"
      WRITE(io_unit, '(a)') "# ikpt band_v band_n absorber_atom absorber_type "// &
                            "eps_v_Ha eps_n_Ha core_energy_Ha occupation_v occupation_n k_weight "// &
                            "loss_energy_Ha loss_energy_eV denominator_real denominator_imag denominator_abs2 "// &
                            "amplitude_real amplitude_imag amplitude_abs2 weighted_strength"
   END SUBROUTINE rixs_open_contribution_table

   SUBROUTINE rixs_write_contribution_rows(io_unit, ikpt, absorber_atom, absorber_type, eig_band, occupation, k_weight, &
                                           core_energy, omega_in, gamma_core, matrix_abs_spin, matrix_emit_spin, loss_grid, &
                                           eta_loss, valence_band_min, valence_band_max, intermediate_band_min, &
                                           intermediate_band_max, contribution_intensity)
      INTEGER, INTENT(IN) :: io_unit, ikpt, absorber_atom, absorber_type
      REAL,    INTENT(IN) :: eig_band(:), occupation(:), k_weight, core_energy, omega_in, gamma_core
      COMPLEX, INTENT(IN) :: matrix_abs_spin(:, :, :), matrix_emit_spin(:, :, :)
      REAL,    INTENT(IN), OPTIONAL :: loss_grid(:), eta_loss
      INTEGER, INTENT(IN) :: valence_band_min, valence_band_max, intermediate_band_min, intermediate_band_max
      REAL,    INTENT(INOUT), OPTIONAL :: contribution_intensity(:)

      INTEGER :: band_v, band_n, i_grid
      REAL :: occupation_v, occupation_n, vacancy_n, loss_energy, denominator_abs2, amplitude_abs2
      REAL :: weighted_strength, gaussian
      COMPLEX :: denominator
      LOGICAL :: l_accumulate_check

      IF (io_unit == -1) RETURN
      IF (SIZE(eig_band) /= SIZE(occupation) .OR. SIZE(matrix_abs_spin, 1) < SIZE(eig_band) .OR. &
          SIZE(matrix_emit_spin, 1) < SIZE(eig_band) .OR. SIZE(matrix_abs_spin, 2) /= SIZE(matrix_emit_spin, 2) .OR. &
          SIZE(matrix_abs_spin, 3) /= SIZE(matrix_emit_spin, 3)) THEN
         CALL juDFT_error("Inconsistent RIXS contribution-table dimensions", calledby="m_rixs_io")
      END IF
      IF (valence_band_min < 1 .OR. valence_band_max > SIZE(eig_band) .OR. valence_band_min > valence_band_max) THEN
         CALL juDFT_error("Invalid RIXS valence band bounds in contribution table", calledby="m_rixs_io")
      END IF
      IF (intermediate_band_min < 1 .OR. intermediate_band_max > SIZE(eig_band) .OR. &
          intermediate_band_min > intermediate_band_max) THEN
         CALL juDFT_error("Invalid RIXS intermediate band bounds in contribution table", calledby="m_rixs_io")
      END IF
      l_accumulate_check = PRESENT(loss_grid) .AND. PRESENT(eta_loss) .AND. PRESENT(contribution_intensity)
      IF (l_accumulate_check) THEN
         IF (SIZE(loss_grid) /= SIZE(contribution_intensity)) THEN
            CALL juDFT_error("RIXS contribution-check spectrum has inconsistent loss-grid size", calledby="m_rixs_io")
         END IF
         IF (eta_loss <= 0.0) CALL juDFT_error("etaLoss must be positive in RIXS contribution check", calledby="m_rixs_io")
      END IF

      DO band_n = intermediate_band_min, intermediate_band_max
         occupation_n = occupation(band_n)
         vacancy_n = 1.0 - occupation_n
         IF (vacancy_n <= 1.0e-10) CYCLE
         denominator = CMPLX(omega_in - (eig_band(band_n) - core_energy), gamma_core)
         denominator_abs2 = ABS(denominator)**2
         DO band_v = valence_band_min, valence_band_max
            occupation_v = occupation(band_v)
            IF (occupation_v <= 1.0e-10) CYCLE
            amplitude_abs2 = rixs_scalar_spin_trace_abs2(matrix_abs_spin(band_n, :, :), matrix_emit_spin(band_v, :, :), &
                                                         denominator)
            IF (amplitude_abs2 < TINY(amplitude_abs2)) CYCLE
            loss_energy = eig_band(band_n) - eig_band(band_v)
            weighted_strength = k_weight*occupation_v*vacancy_n*amplitude_abs2
            IF (l_accumulate_check) THEN
               DO i_grid = 1, SIZE(loss_grid)
                  gaussian = xas_gaussian_broadening(loss_grid(i_grid) - loss_energy, eta_loss)
                  IF (gaussian == 0.0) CYCLE
                  contribution_intensity(i_grid) = contribution_intensity(i_grid) + weighted_strength*gaussian
               END DO
            END IF
            WRITE(io_unit, '(5(i0,1x),15(es24.16e3,1x))') ikpt, band_v, band_n, absorber_atom, absorber_type, &
               eig_band(band_v), eig_band(band_n), core_energy, occupation_v, occupation_n, k_weight, &
               loss_energy, loss_energy*hartree_to_ev_const, REAL(denominator), AIMAG(denominator), denominator_abs2, &
               0.0, 0.0, amplitude_abs2, weighted_strength
         END DO
      END DO
   END SUBROUTINE rixs_write_contribution_rows

   SUBROUTINE rixs_close_contribution_table(io_unit)
      INTEGER, INTENT(INOUT) :: io_unit
      INTEGER :: io_status

      IF (io_unit == -1) RETURN
      CLOSE(io_unit, IOSTAT=io_status)
      io_unit = -1
      IF (io_status /= 0) CALL juDFT_error("Cannot close RIXS contribution table", calledby="m_rixs_io")
   END SUBROUTINE rixs_close_contribution_table

   SUBROUTINE rixs_print_setup_summary(rixs, l_noco, l_soc)
      TYPE(t_xas), INTENT(IN) :: rixs
      LOGICAL,     INTENT(IN) :: l_noco, l_soc

      WRITE(*, '(a)') " ---------- RIXS setup ------------------------------"
      WRITE(*, '(a,l1)') " RIXS enabled             : ", rixs%l_rixs
      WRITE(*, '(a,i0)') " Absorber Z               : ", rixs%rixs_absorber_z
      WRITE(*, '(a,a)') " Edge                     : ", TRIM(rixs%rixs_edge)
      WRITE(*, '(a,f12.6,a,f12.6,a)') " Incident energy omegaIn  : ", rixs%rixs_omega_in, " Ha (", &
                                      rixs%rixs_omega_in*hartree_to_ev_const, " eV)"
      WRITE(*, '(a,f12.6,a,f12.6,a)') " Core broadening Gamma    : ", rixs%rixs_gamma_core, " Ha (", &
                                      rixs%rixs_gamma_core*hartree_to_ev_const, " eV)"
      WRITE(*, '(a,f12.6,a,f12.6,a)') " Loss broadening etaLoss  : ", rixs%rixs_eta_loss, " Ha (", &
                                      rixs%rixs_eta_loss*hartree_to_ev_const, " eV)"
      WRITE(*, '(a,f12.6,a,f12.6,a,f12.6,a,f12.6,a)') " Loss window              : ", &
         rixs%rixs_loss_min, " ... ", rixs%rixs_loss_max, " Ha (", &
         rixs%rixs_loss_min*hartree_to_ev_const, " ... ", rixs%rixs_loss_max*hartree_to_ev_const, " eV)"
      WRITE(*, '(a,i0)') " Number of loss points    : ", rixs%rixs_n_loss
      WRITE(*, '(a,a)') " Incoming polarizations   : ", TRIM(rixs_polarization_string(rixs%rixs_in_polarizations))
      WRITE(*, '(a,a)') " Outgoing polarizations   : ", TRIM(rixs_polarization_string(rixs%rixs_out_polarizations))
      WRITE(*, '(a,a)') " Output prefix            : ", TRIM(rixs%rixs_output_prefix)
      WRITE(*, '(a)') " Approximation            : direct same-k independent-particle"
      IF (l_noco) THEN
         WRITE(*, '(a)') " Spinor treatment         : first-variation coherent core-mj amplitude"
      ELSE
         WRITE(*, '(a)') " Scalar spin treatment    : spin-degenerate S1 trace"
      END IF
      WRITE(*, '(a)') " Symmetry treatment       : full-k only, no star reconstruction"
      WRITE(*, '(a)') " Absorber-site sum        : incoherent local-site intensities"
      WRITE(*, '(a,a)') " Valence band window      : ", TRIM(rixs_band_window_string(rixs%l_rixs_valence_band_min, &
                                                                                         rixs%rixs_valence_band_min, &
                                                                                         rixs%l_rixs_valence_band_max, &
                                                                                         rixs%rixs_valence_band_max))
      WRITE(*, '(a,a)') " Intermediate band window : ", TRIM(rixs_band_window_string(rixs%l_rixs_intermediate_band_min, &
                                                                                         rixs%rixs_intermediate_band_min, &
                                                                                         rixs%l_rixs_intermediate_band_max, &
                                                                                         rixs%rixs_intermediate_band_max))
      IF (rixs%l_rixs_valence_band_max .OR. rixs%l_rixs_intermediate_band_max) THEN
         WRITE(*, '(a)') " NOTE: RIXS band-window upper bounds are clamped to the available bands at each k point."
      END IF
      WRITE(*, '(a,l1)') " Noco enabled             : ", l_noco
      WRITE(*, '(a,l1)') " SOC enabled              : ", l_soc
      WRITE(*, '(a,l1)') " Write contributions      : ", rixs%rixs_write_contributions
      IF (l_noco .AND. .NOT. rixs%rixs_write_contributions) THEN
         WRITE(*, '(a)') " Spinor contribution output: not requested"
      END IF
      IF (rixs%rixs_write_contributions) THEN
         WRITE(*, '(a)') " NOTE: RIXS contribution output may be very large."
         WRITE(*, '(a)') "       One rank-suffixed table is written per requested polarization pair and MPI rank."
      END IF
      WRITE(*, '(a)') " ----------------------------------------------------"
   END SUBROUTINE rixs_print_setup_summary

   SUBROUTINE rixs_print_pair_summary(pol_in, pol_out, loss_grid, intensity)
      CHARACTER(LEN=*), INTENT(IN) :: pol_in, pol_out
      REAL,             INTENT(IN) :: loss_grid(:), intensity(:)

      INTEGER :: i_max
      REAL :: d_omega, integrated

      IF (SIZE(loss_grid) /= SIZE(intensity)) THEN
         CALL juDFT_error("loss_grid and intensity sizes differ in RIXS summary", calledby="m_rixs_io")
      END IF
      d_omega = 0.0
      IF (SIZE(loss_grid) > 1) d_omega = (loss_grid(SIZE(loss_grid)) - loss_grid(1))/REAL(SIZE(loss_grid) - 1)
      integrated = SUM(intensity)*d_omega
      i_max = MAXLOC(intensity, DIM=1)

      WRITE(*, '(a,a,a,es18.10)') " RIXS ", TRIM(pol_in)//TRIM(pol_out), ": integrated intensity = ", integrated
      WRITE(*, '(a,a,a,es18.10)') " RIXS ", TRIM(pol_in)//TRIM(pol_out), ": max intensity        = ", intensity(i_max)
      WRITE(*, '(a,a,a,es18.10,a,es18.10,a)') " RIXS ", TRIM(pol_in)//TRIM(pol_out), ": max at loss          = ", &
         loss_grid(i_max), " Ha (", loss_grid(i_max)*hartree_to_ev_const, " eV)"
   END SUBROUTINE rixs_print_pair_summary

   SUBROUTINE rixs_print_contribution_check(pol_in, pol_out, spectrum, contribution_spectrum)
      CHARACTER(LEN=*), INTENT(IN) :: pol_in, pol_out
      REAL,             INTENT(IN) :: spectrum(:), contribution_spectrum(:)

      REAL, PARAMETER :: abs_tol = 1.0e-12, rel_tol = 1.0e-10
      REAL :: max_abs_diff, max_rel_diff, norm
      CHARACTER(LEN=4) :: status

      IF (SIZE(spectrum) /= SIZE(contribution_spectrum)) THEN
         CALL juDFT_error("RIXS contribution check spectra have different sizes", calledby="m_rixs_io")
      END IF

      max_abs_diff = MAXVAL(ABS(spectrum - contribution_spectrum))
      norm = MAX(MAXVAL(ABS(spectrum)), TINY(norm))
      max_rel_diff = max_abs_diff/norm
      IF (max_abs_diff < abs_tol .OR. max_rel_diff < rel_tol) THEN
         status = "PASS"
      ELSE
         status = "FAIL"
      END IF

      WRITE(*, '(a,a,a)') " RIXS ", TRIM(pol_in)//TRIM(pol_out), " contribution-spectrum check:"
      WRITE(*, '(a,es18.10)') "   max abs diff = ", max_abs_diff
      WRITE(*, '(a,es18.10)') "   max rel diff = ", max_rel_diff
      WRITE(*, '(a,a)') "   status       = ", TRIM(status)
   END SUBROUTINE rixs_print_contribution_check

   FUNCTION rixs_band_window_string(l_min, band_min, l_max, band_max) RESULT(label)
      LOGICAL, INTENT(IN) :: l_min, l_max
      INTEGER, INTENT(IN) :: band_min, band_max
      CHARACTER(LEN=32) :: label

      IF (.NOT. l_min .AND. .NOT. l_max) THEN
         label = "all"
      ELSE IF (l_min .AND. l_max) THEN
         WRITE(label, '(i0,a,i0)') band_min, " ... ", band_max
      ELSE IF (l_min) THEN
         WRITE(label, '(i0,a)') band_min, " ... all"
      ELSE
         WRITE(label, '(a,i0)') "1 ... ", band_max
      END IF
   END FUNCTION rixs_band_window_string

   FUNCTION rixs_polarization_string(polarizations) RESULT(pol_string)
      LOGICAL, INTENT(IN) :: polarizations(3)
      CHARACTER(LEN=16) :: pol_string

      pol_string = ""
      IF (polarizations(1)) pol_string = TRIM(pol_string)//" x"
      IF (polarizations(2)) pol_string = TRIM(pol_string)//" y"
      IF (polarizations(3)) pol_string = TRIM(pol_string)//" z"
      pol_string = ADJUSTL(pol_string)
   END FUNCTION rixs_polarization_string

   FUNCTION rixs_energy_label(energy) RESULT(label)
      REAL, INTENT(IN) :: energy
      CHARACTER(LEN=12) :: label
      INTEGER :: i_char

      WRITE(label, '(f12.6)') energy
      label = ADJUSTL(label)
      DO i_char = 1, LEN(label)
         IF (label(i_char:i_char) == ".") label(i_char:i_char) = "p"
         IF (label(i_char:i_char) == "-") label(i_char:i_char) = "m"
         IF (label(i_char:i_char) == " ") label(i_char:i_char) = ""
      END DO
   END FUNCTION rixs_energy_label

END MODULE m_rixs_io
