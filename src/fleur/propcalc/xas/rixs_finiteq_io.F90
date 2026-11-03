!--------------------------------------------------------------------------------
! Copyright (c) 2026 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------

MODULE m_rixs_finiteq_io
   USE m_constants, ONLY: hartree_to_ev_const
   USE m_juDFT, ONLY: juDFT_error
   USE m_rixs_spectrum, ONLY: rixs_occupation_tolerance
   USE m_xas_spectrum, ONLY: xas_gaussian_broadening
   IMPLICIT NONE
   PRIVATE

   PUBLIC :: rixs_open_finiteq_pair_table, rixs_open_finiteq_site_table
   PUBLIC :: rixs_write_finiteq_pair_rows, rixs_write_finiteq_site_rows
   PUBLIC :: rixs_finiteq_label

CONTAINS

   SUBROUTINE rixs_open_finiteq_pair_table(filename, edge, absorber_z, pol_in, pol_out, omega_in, q_full_rlu, q_reduced_rlu, &
                                           mpi_rank, io_unit)
      CHARACTER(LEN=*), INTENT(IN) :: filename, edge, pol_in, pol_out
      INTEGER, INTENT(IN) :: absorber_z, mpi_rank
      REAL, INTENT(IN) :: omega_in, q_full_rlu(3), q_reduced_rlu(3)
      INTEGER, INTENT(OUT) :: io_unit
      INTEGER :: io_status
      CHARACTER(LEN=256) :: io_message

      OPEN(NEWUNIT=io_unit, FILE=TRIM(filename), STATUS="replace", ACTION="write", IOSTAT=io_status, IOMSG=io_message)
      IF (io_status /= 0) CALL juDFT_error("Cannot open finite-Q RIXS pair table: "//TRIM(io_message), &
                                           calledby="m_rixs_finiteq_io")
      WRITE(io_unit, '(a)') "# FLEUR finite-Q first-variation spinor RIXS coherent pair table"
      WRITE(io_unit, '(a,a,1x,a,i0,1x,a,a,1x,a,a)') "# edge=", TRIM(edge), "absorberZ=", absorber_z, &
                                                    "incoming=", TRIM(pol_in), "outgoing=", TRIM(pol_out)
      WRITE(io_unit, '(a,es24.16e3,a,3es24.16e3,a,3es24.16e3)') "# omegaIn_Ha=", omega_in, &
         " Q_full_rlu=", q_full_rlu, " q_reduced_rlu=", q_reduced_rlu
      WRITE(io_unit, '(a,i0)') "# MPI rank = ", mpi_rank
      WRITE(io_unit, '(a)') "# Site amplitudes have already been summed coherently before amplitude_abs2 is formed."
      WRITE(io_unit, '(a)') "# reciprocal_shift is the integer vector with k_v + Q_full = k_n + reciprocal_shift."
      WRITE(io_unit, '(a)') "# columns: ikpt_v ikpt_n band_v band_n k_v(3) k_n(3) reciprocal_shift(3) "// &
         "eps_v_Ha eps_n_Ha occupation_v occupation_n k_weight loss_Ha loss_eV amplitude_real amplitude_imag "// &
         "amplitude_abs2 weighted_strength"
   END SUBROUTINE rixs_open_finiteq_pair_table

   SUBROUTINE rixs_open_finiteq_site_table(filename, edge, absorber_z, pol_in, pol_out, omega_in, q_full_rlu, q_reduced_rlu, &
                                           mpi_rank, io_unit)
      CHARACTER(LEN=*), INTENT(IN) :: filename, edge, pol_in, pol_out
      INTEGER, INTENT(IN) :: absorber_z, mpi_rank
      REAL, INTENT(IN) :: omega_in, q_full_rlu(3), q_reduced_rlu(3)
      INTEGER, INTENT(OUT) :: io_unit
      INTEGER :: io_status
      CHARACTER(LEN=256) :: io_message

      OPEN(NEWUNIT=io_unit, FILE=TRIM(filename), STATUS="replace", ACTION="write", IOSTAT=io_status, IOMSG=io_message)
      IF (io_status /= 0) CALL juDFT_error("Cannot open finite-Q RIXS site table: "//TRIM(io_message), &
                                           calledby="m_rixs_finiteq_io")
      WRITE(io_unit, '(a)') "# FLEUR finite-Q first-variation spinor RIXS site-partial amplitude table"
      WRITE(io_unit, '(a,a,1x,a,i0,1x,a,a,1x,a,a)') "# edge=", TRIM(edge), "absorberZ=", absorber_z, &
                                                    "incoming=", TRIM(pol_in), "outgoing=", TRIM(pol_out)
      WRITE(io_unit, '(a,es24.16e3,a,3es24.16e3,a,3es24.16e3)') "# omegaIn_Ha=", omega_in, &
         " Q_full_rlu=", q_full_rlu, " q_reduced_rlu=", q_reduced_rlu
      WRITE(io_unit, '(a,i0)') "# MPI rank = ", mpi_rank
      WRITE(io_unit, '(a)') "# Sum phased_partial as complex numbers over absorber sites/types to recover the pair amplitude."
      WRITE(io_unit, '(a)') "# Do not sum site abs2 values to form a spectrum."
      WRITE(io_unit, '(a)') "# columns: ikpt_v ikpt_n band_v band_n absorber_atom absorber_type tau_fractional(3) "// &
         "phase_real phase_imag denominator_real denominator_imag denominator_abs2 local_real local_imag "// &
         "phased_partial_real phased_partial_imag"
   END SUBROUTINE rixs_open_finiteq_site_table

   SUBROUTINE rixs_write_finiteq_site_rows(io_unit, ikpt_v, ikpt_n, absorber_atom, absorber_type, tau_fractional, &
                                           site_phase, eig_intermediate, occupation_valence, occupation_intermediate, &
                                           core_energy, omega_in, gamma_core, site_partial_vn, valence_band_min, &
                                           valence_band_max, intermediate_band_min, intermediate_band_max)
      INTEGER, INTENT(IN) :: io_unit, ikpt_v, ikpt_n, absorber_atom, absorber_type
      REAL, INTENT(IN) :: tau_fractional(3), eig_intermediate(:), occupation_valence(:), occupation_intermediate(:)
      REAL, INTENT(IN) :: core_energy, omega_in, gamma_core
      COMPLEX, INTENT(IN) :: site_phase, site_partial_vn(:, :)
      INTEGER, INTENT(IN) :: valence_band_min, valence_band_max, intermediate_band_min, intermediate_band_max
      INTEGER :: band_v, band_n, io_status
      COMPLEX :: denominator, local_amplitude
      CHARACTER(LEN=256) :: io_message

      IF (io_unit == -1) RETURN
      DO band_v = valence_band_min, valence_band_max
         IF (occupation_valence(band_v) <= rixs_occupation_tolerance) CYCLE
         DO band_n = intermediate_band_min, intermediate_band_max
            IF (1.0 - occupation_intermediate(band_n) <= rixs_occupation_tolerance) CYCLE
            denominator = CMPLX(omega_in - (eig_intermediate(band_n) - core_energy), gamma_core)
            local_amplitude = site_partial_vn(band_v, band_n)/site_phase
            WRITE(io_unit, '(6(i0,1x),12(es24.16e3,1x))', IOSTAT=io_status, IOMSG=io_message) &
               ikpt_v, ikpt_n, band_v, band_n, absorber_atom, absorber_type, tau_fractional, &
               REAL(site_phase), AIMAG(site_phase), REAL(denominator), AIMAG(denominator), ABS(denominator)**2, &
               REAL(local_amplitude), AIMAG(local_amplitude), REAL(site_partial_vn(band_v, band_n)), &
               AIMAG(site_partial_vn(band_v, band_n))
            IF (io_status /= 0) CALL juDFT_error("Cannot write finite-Q RIXS site row: "//TRIM(io_message), &
                                                 calledby="m_rixs_finiteq_io")
         END DO
      END DO
   END SUBROUTINE rixs_write_finiteq_site_rows

   SUBROUTINE rixs_write_finiteq_pair_rows(io_unit, ikpt_v, ikpt_n, k_v, k_n, reciprocal_shift, eig_valence, &
                                           occupation_valence, eig_intermediate, occupation_intermediate, k_weight, &
                                           amplitude_vn, loss_grid, eta_loss, valence_band_min, valence_band_max, &
                                           intermediate_band_min, intermediate_band_max, contribution_intensity)
      INTEGER, INTENT(IN) :: io_unit, ikpt_v, ikpt_n, reciprocal_shift(3)
      REAL, INTENT(IN) :: k_v(3), k_n(3), eig_valence(:), occupation_valence(:), eig_intermediate(:)
      REAL, INTENT(IN) :: occupation_intermediate(:), k_weight, loss_grid(:), eta_loss
      COMPLEX, INTENT(IN) :: amplitude_vn(:, :)
      INTEGER, INTENT(IN) :: valence_band_min, valence_band_max, intermediate_band_min, intermediate_band_max
      REAL, INTENT(INOUT) :: contribution_intensity(:)
      INTEGER :: band_v, band_n, i_grid, io_status
      REAL :: loss_energy, strength, weight, weighted_strength, gaussian
      CHARACTER(LEN=256) :: io_message

      IF (io_unit == -1) RETURN
      DO band_v = valence_band_min, valence_band_max
         IF (occupation_valence(band_v) <= rixs_occupation_tolerance) CYCLE
         DO band_n = intermediate_band_min, intermediate_band_max
            IF (1.0 - occupation_intermediate(band_n) <= rixs_occupation_tolerance) CYCLE
            strength = ABS(amplitude_vn(band_v, band_n))**2
            weight = k_weight*occupation_valence(band_v)*(1.0 - occupation_intermediate(band_n))
            loss_energy = eig_intermediate(band_n) - eig_valence(band_v)
            weighted_strength = weight*strength
            DO i_grid = 1, SIZE(loss_grid)
               gaussian = xas_gaussian_broadening(loss_grid(i_grid) - loss_energy, eta_loss)
               contribution_intensity(i_grid) = contribution_intensity(i_grid) + weighted_strength*gaussian
            END DO
            WRITE(io_unit, '(4(i0,1x),6(es24.16e3,1x),3(i0,1x),11(es24.16e3,1x))', &
                  IOSTAT=io_status, IOMSG=io_message) ikpt_v, ikpt_n, band_v, band_n, k_v, k_n, reciprocal_shift, &
               eig_valence(band_v), eig_intermediate(band_n), occupation_valence(band_v), &
               occupation_intermediate(band_n), k_weight, loss_energy, loss_energy*hartree_to_ev_const, &
               REAL(amplitude_vn(band_v, band_n)), AIMAG(amplitude_vn(band_v, band_n)), strength, weighted_strength
            IF (io_status /= 0) CALL juDFT_error("Cannot write finite-Q RIXS pair row: "//TRIM(io_message), &
                                                 calledby="m_rixs_finiteq_io")
         END DO
      END DO
   END SUBROUTINE rixs_write_finiteq_pair_rows

   FUNCTION rixs_finiteq_label(q_full_rlu) RESULT(label)
      REAL, INTENT(IN) :: q_full_rlu(3)
      CHARACTER(LEN=60) :: label
      CHARACTER(LEN=18) :: part
      INTEGER :: i, j

      label = "Q"
      DO i = 1, 3
         WRITE(part, '(f12.6)') q_full_rlu(i)
         part = ADJUSTL(part)
         DO j = 1, LEN_TRIM(part)
            IF (part(j:j) == ".") part(j:j) = "p"
            IF (part(j:j) == "-") part(j:j) = "m"
         END DO
         label = TRIM(label)//"_"//TRIM(part)
      END DO
   END FUNCTION rixs_finiteq_label

END MODULE m_rixs_finiteq_io
