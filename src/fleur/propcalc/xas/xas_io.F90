!--------------------------------------------------------------------------------
! Copyright (c) 2026 Peter Gruenberg Institut, Forschungszentrum Juelich, Germany
! This file is part of FLEUR and available as free software under the conditions
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------

MODULE m_xas_io
   USE m_juDFT, ONLY: juDFT_error
   USE m_xas_core, ONLY: t_xas_core_state
   USE m_xas_transitions, ONLY: t_xas_core_descriptor, t_xas_transition_record, &
                                xas_attach_l_channel_amplitudes, xas_attach_polarization_amplitude, &
                                xas_core_descriptor_from_state, xas_init_transition_record, xas_transition_absM2, &
                                xas_transition_l_reconstruction_error, xas_transition_total_weighted_strength
   IMPLICIT NONE
   PRIVATE

   PUBLIC :: xas_write_spectrum_text
   PUBLIC :: xas_open_transition_table, xas_write_transition_rows, xas_close_transition_table
   PUBLIC :: xas_l_channel_reconstruction_error

CONTAINS

   SUBROUTINE xas_write_spectrum_text(filename, energy_grid, intensity, energy_unit)
      CHARACTER(LEN=*), INTENT(IN) :: filename
      REAL,             INTENT(IN) :: energy_grid(:)
      REAL,             INTENT(IN) :: intensity(:)
      CHARACTER(LEN=*), INTENT(IN) :: energy_unit

      INTEGER :: io_unit, i_grid, io_status
      CHARACTER(LEN=256) :: io_message

      IF (SIZE(energy_grid) /= SIZE(intensity)) THEN
         CALL juDFT_error("energy_grid and intensity sizes differ in xas_write_spectrum_text", calledby="m_xas_io")
      END IF

      OPEN(NEWUNIT=io_unit, FILE=TRIM(filename), STATUS="replace", ACTION="write", FORM="formatted", &
           IOSTAT=io_status, IOMSG=io_message)
      IF (io_status /= 0) THEN
         CALL juDFT_error("Cannot open XAS spectrum text file: "//TRIM(io_message), calledby="m_xas_io")
      END IF

      WRITE(io_unit, '(a)', IOSTAT=io_status, IOMSG=io_message) "# XAS spectrum"
      IF (io_status == 0) WRITE(io_unit, '(a)', IOSTAT=io_status, IOMSG=io_message) "# energy_unit: "//TRIM(energy_unit)
      IF (io_status == 0) WRITE(io_unit, '(a,i0)', IOSTAT=io_status, IOMSG=io_message) "# n_grid: ", SIZE(energy_grid)
      IF (io_status == 0) WRITE(io_unit, '(a)', IOSTAT=io_status, IOMSG=io_message) "# intensity: arbitrary units"
      IF (io_status == 0) WRITE(io_unit, '(a)', IOSTAT=io_status, IOMSG=io_message) "# columns: energy intensity"

      DO i_grid = 1, SIZE(energy_grid)
         IF (io_status /= 0) EXIT
         WRITE(io_unit, '(2es24.16e3)', IOSTAT=io_status, IOMSG=io_message) energy_grid(i_grid), intensity(i_grid)
      END DO

      CLOSE(io_unit)

      IF (io_status /= 0) THEN
         CALL juDFT_error("Cannot write XAS spectrum text file: "//TRIM(io_message), calledby="m_xas_io")
      END IF
   END SUBROUTINE xas_write_spectrum_text

   SUBROUTINE xas_open_transition_table(filename, edge, absorber_z, polarization, mpi_rank, io_unit)
      CHARACTER(LEN=*), INTENT(IN) :: filename, edge, polarization
      INTEGER,          INTENT(IN) :: absorber_z, mpi_rank
      INTEGER,          INTENT(OUT) :: io_unit

      INTEGER :: io_status
      CHARACTER(LEN=256) :: io_message

      OPEN(NEWUNIT=io_unit, FILE=TRIM(filename), STATUS="replace", ACTION="write", FORM="formatted", &
           IOSTAT=io_status, IOMSG=io_message)
      IF (io_status /= 0) THEN
         CALL juDFT_error("Cannot open XAS transition table: "//TRIM(io_message), calledby="m_xas_io")
      END IF
      WRITE(io_unit, '(a)') "# FLEUR independent-particle XAS transition table"
      WRITE(io_unit, '(a,a)') "# edge = ", TRIM(edge)
      WRITE(io_unit, '(a,i0)') "# absorberZ = ", absorber_z
      WRITE(io_unit, '(a,a)') "# polarization = ", TRIM(polarization)
      WRITE(io_unit, '(a,i0)') "# MPI rank = ", mpi_rank
      WRITE(io_unit, '(a)') "# This file contains only transitions evaluated by this MPI rank."
      WRITE(io_unit, '(a)') "# columns:"
      WRITE(io_unit, '(a)') "# ikpt_full ikpt_parent star_index band absorber_atom absorber_type "// &
                            "transition_energy_Ha transition_energy_eV occupation k_weight one_minus_occ absM2 "// &
                            "total_weighted_strength"
   END SUBROUTINE xas_open_transition_table

   SUBROUTINE xas_write_transition_rows(io_unit, ikpt_full, ikpt_parent, star_index, bksym, absorber_atom, absorber_type, &
                                        spin_channel, core_state, pol_id, pol_label, eps_sph, eig_band, occupation, &
                                        k_weight, matrix, hartree_to_ev, final_l, matrix_lchan)
      INTEGER, INTENT(IN) :: io_unit, ikpt_full, ikpt_parent, star_index, bksym, absorber_atom, absorber_type, spin_channel
      TYPE(t_xas_core_state), INTENT(IN) :: core_state
      INTEGER, INTENT(IN) :: pol_id
      CHARACTER(LEN=*), INTENT(IN) :: pol_label
      COMPLEX, INTENT(IN) :: eps_sph(:)
      REAL,    INTENT(IN) :: eig_band(:), occupation(:), k_weight, hartree_to_ev
      COMPLEX, INTENT(IN) :: matrix(:, :)
      INTEGER, OPTIONAL, INTENT(IN) :: final_l(:)
      COMPLEX, OPTIONAL, INTENT(IN) :: matrix_lchan(:, :, :)

      TYPE(t_xas_core_descriptor) :: core
      TYPE(t_xas_transition_record) :: record
      INTEGER :: band
      REAL :: one_minus_occ, abs_m2, transition_energy

      IF (SIZE(eig_band) /= SIZE(occupation) .OR. SIZE(matrix, 1) < SIZE(eig_band)) THEN
         CALL juDFT_error("Inconsistent XAS transition-table dimensions", calledby="m_xas_io")
      END IF
      CALL xas_core_descriptor_from_state(core, core_state)
      DO band = 1, SIZE(eig_band)
         one_minus_occ = 1.0 - occupation(band)
         IF (one_minus_occ <= 1.0e-10) CYCLE
         CALL xas_init_transition_record(record, core, ikpt_parent, ikpt_full, star_index, bksym, band, &
                                         absorber_atom, absorber_type, spin_channel, eig_band(band), &
                                         occupation(band), k_weight)
         CALL xas_attach_polarization_amplitude(record, pol_id, pol_label, eps_sph, matrix, band)
         IF (PRESENT(final_l) .AND. PRESENT(matrix_lchan)) THEN
            CALL xas_attach_l_channel_amplitudes(record, final_l, matrix_lchan, band)
         END IF
         transition_energy = record%transition_energy
         abs_m2 = xas_transition_absM2(record)
         WRITE(io_unit, '(6(i0,1x),7(es24.16e3,1x))') ikpt_full, ikpt_parent, star_index, band, &
            absorber_atom, absorber_type, transition_energy, transition_energy*hartree_to_ev, record%occupation, &
            record%k_weight, record%one_minus_occ, abs_m2, xas_transition_total_weighted_strength(record)
      END DO
   END SUBROUTINE xas_write_transition_rows

   SUBROUTINE xas_l_channel_reconstruction_error(max_error, core_state, pol_id, pol_label, eps_sph, eig_band, occupation, &
                                                 k_weight, matrix, final_l, matrix_lchan, max_checks)
      REAL,                  INTENT(INOUT) :: max_error
      TYPE(t_xas_core_state), INTENT(IN)    :: core_state
      INTEGER,               INTENT(IN)    :: pol_id
      CHARACTER(LEN=*),      INTENT(IN)    :: pol_label
      COMPLEX,               INTENT(IN)    :: eps_sph(:)
      REAL,                  INTENT(IN)    :: eig_band(:), occupation(:), k_weight
      COMPLEX,               INTENT(IN)    :: matrix(:, :)
      INTEGER,               INTENT(IN)    :: final_l(:)
      COMPLEX,               INTENT(IN)    :: matrix_lchan(:, :, :)
      INTEGER, OPTIONAL,     INTENT(IN)    :: max_checks

      TYPE(t_xas_core_descriptor) :: core
      TYPE(t_xas_transition_record) :: record
      INTEGER :: band, n_checked, n_limit

      n_limit = HUGE(1)
      IF (PRESENT(max_checks)) n_limit = max_checks
      CALL xas_core_descriptor_from_state(core, core_state)
      n_checked = 0
      DO band = 1, SIZE(eig_band)
         IF (1.0 - occupation(band) <= 1.0e-10) CYCLE
         CALL xas_init_transition_record(record, core, 1, 1, 1, 1, band, 1, 1, 1, eig_band(band), &
                                         occupation(band), k_weight)
         CALL xas_attach_polarization_amplitude(record, pol_id, pol_label, eps_sph, matrix, band)
         CALL xas_attach_l_channel_amplitudes(record, final_l, matrix_lchan, band)
         max_error = MAX(max_error, xas_transition_l_reconstruction_error(record))
         n_checked = n_checked + 1
         IF (n_checked >= n_limit) EXIT
      END DO
   END SUBROUTINE xas_l_channel_reconstruction_error

   SUBROUTINE xas_close_transition_table(io_unit)
      INTEGER, INTENT(INOUT) :: io_unit
      INTEGER :: io_status

      IF (io_unit == -1) RETURN
      CLOSE(io_unit, IOSTAT=io_status)
      io_unit = -1
      IF (io_status /= 0) CALL juDFT_error("Cannot close XAS transition table", calledby="m_xas_io")
   END SUBROUTINE xas_close_transition_table

END MODULE m_xas_io
