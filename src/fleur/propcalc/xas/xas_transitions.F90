!--------------------------------------------------------------------------------
! Copyright (c) 2026 Peter Gruenberg Institut, Forschungszentrum Juelich, Germany
! This file is part of FLEUR and available as free software under the conditions
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------

MODULE m_xas_transitions
   USE m_juDFT, ONLY: juDFT_error
   USE m_xas_amplitudes, ONLY: t_xas_transition_amplitudes
   USE m_xas_core, ONLY: t_xas_core_state
   IMPLICIT NONE
   PRIVATE

   PUBLIC :: t_xas_core_descriptor
   PUBLIC :: t_xas_polarized_amplitudes
   PUBLIC :: t_xas_transition_record
   PUBLIC :: xas_core_descriptor_from_state
   PUBLIC :: xas_init_transition_record
   PUBLIC :: xas_attach_polarization_amplitude
   PUBLIC :: xas_attach_l_channel_amplitudes
   PUBLIC :: xas_allowed_final_l_channels
   PUBLIC :: xas_transition_absM2
   PUBLIC :: xas_transition_absM2_l
   PUBLIC :: xas_transition_l_reconstruction_error
   PUBLIC :: xas_transition_total_weighted_strength

   TYPE t_xas_core_descriptor
      CHARACTER(LEN=8) :: edge_name = ""
      INTEGER :: edge_id = 0
      INTEGER :: core_state_index = 1
      INTEGER :: n = 0
      INTEGER :: lc = -1
      INTEGER :: kappa = 0
      INTEGER :: twice_j = 0
      REAL    :: energy = 0.0
      INTEGER, ALLOCATABLE :: twice_mj(:)
   END TYPE t_xas_core_descriptor

   TYPE t_xas_polarized_amplitudes
      INTEGER :: pol_id = 0
      CHARACTER(LEN=1) :: pol_label = ""
      COMPLEX, ALLOCATABLE :: eps_sph(:)
      TYPE(t_xas_transition_amplitudes) :: amp
      ! Complex final-l amplitudes. Their incoherent strengths are useful
      ! diagnostics, but the full XAS observable uses the coherent total amp.
      INTEGER, ALLOCATABLE :: final_l(:)
      TYPE(t_xas_transition_amplitudes), ALLOCATABLE :: amp_l(:)
   END TYPE t_xas_polarized_amplitudes

   TYPE t_xas_transition_record
      ! Local/streaming transition metadata. This is intentionally not a
      ! global transition database: XAS writes or accumulates each record
      ! immediately, while later RIXS can reuse the complex amplitudes before
      ! any incoherent |M|^2 reduction.
      INTEGER :: ikpt_parent = 0
      INTEGER :: ikpt_full = 0
      INTEGER :: star_index = 1
      INTEGER :: bksym = 1
      INTEGER :: band = 0
      INTEGER :: absorber_atom = 0
      INTEGER :: absorber_type = 0
      INTEGER :: spin_channel = 0
      REAL    :: transition_energy = 0.0
      REAL    :: eig = 0.0
      REAL    :: core_energy = 0.0
      REAL    :: occupation = 0.0
      REAL    :: one_minus_occ = 0.0
      REAL    :: k_weight = 0.0
      TYPE(t_xas_core_descriptor) :: core
      TYPE(t_xas_polarized_amplitudes), ALLOCATABLE :: pol_amp(:)
   END TYPE t_xas_transition_record

CONTAINS

   SUBROUTINE xas_core_descriptor_from_state(core, core_state, core_state_index)
      TYPE(t_xas_core_descriptor), INTENT(OUT) :: core
      TYPE(t_xas_core_state),      INTENT(IN)  :: core_state
      INTEGER, OPTIONAL,           INTENT(IN)  :: core_state_index

      IF (.NOT. ALLOCATED(core_state%twice_mj)) THEN
         CALL juDFT_error("XAS core state has no twice_mj mapping", calledby="m_xas_transitions")
      END IF

      core%edge_name = core_state%edge_name
      core%edge_id = xas_edge_id(core_state%edge_name)
      IF (PRESENT(core_state_index)) core%core_state_index = core_state_index
      core%n = core_state%n
      core%lc = core_state%lc
      core%kappa = core_state%kappa
      core%twice_j = core_state%twice_j
      core%energy = core_state%energy
      ALLOCATE(core%twice_mj(SIZE(core_state%twice_mj)))
      core%twice_mj = core_state%twice_mj
   END SUBROUTINE xas_core_descriptor_from_state

   SUBROUTINE xas_init_transition_record(record, core, ikpt_parent, ikpt_full, star_index, bksym, band, &
                                         absorber_atom, absorber_type, spin_channel, eig, occupation, k_weight)
      TYPE(t_xas_transition_record), INTENT(OUT) :: record
      TYPE(t_xas_core_descriptor),   INTENT(IN)  :: core
      INTEGER,                       INTENT(IN)  :: ikpt_parent, ikpt_full, star_index, bksym, band
      INTEGER,                       INTENT(IN)  :: absorber_atom, absorber_type, spin_channel
      REAL,                          INTENT(IN)  :: eig, occupation, k_weight

      record%ikpt_parent = ikpt_parent
      record%ikpt_full = ikpt_full
      record%star_index = star_index
      record%bksym = bksym
      record%band = band
      record%absorber_atom = absorber_atom
      record%absorber_type = absorber_type
      record%spin_channel = spin_channel
      record%eig = eig
      record%core_energy = core%energy
      record%transition_energy = eig - core%energy
      record%occupation = occupation
      record%one_minus_occ = 1.0 - occupation
      record%k_weight = k_weight
      record%core = core
   END SUBROUTINE xas_init_transition_record

   SUBROUTINE xas_attach_polarization_amplitude(record, pol_id, pol_label, eps_sph, matrix, band)
      TYPE(t_xas_transition_record), INTENT(INOUT) :: record
      INTEGER,                       INTENT(IN)    :: pol_id, band
      CHARACTER(LEN=*),              INTENT(IN)    :: pol_label
      COMPLEX,                       INTENT(IN)    :: eps_sph(:)
      COMPLEX,                       INTENT(IN)    :: matrix(:, :)

      INTEGER :: n_pol

      IF (.NOT. ALLOCATED(record%core%twice_mj)) THEN
         CALL juDFT_error("XAS transition record has no core twice_mj mapping", calledby="m_xas_transitions")
      END IF
      IF (band /= record%band) THEN
         CALL juDFT_error("XAS transition record band mismatch", calledby="m_xas_transitions")
      END IF

      n_pol = 1
      IF (ALLOCATED(record%pol_amp)) DEALLOCATE(record%pol_amp)
      ALLOCATE(record%pol_amp(n_pol))
      record%pol_amp(n_pol)%pol_id = pol_id
      record%pol_amp(n_pol)%pol_label = ADJUSTL(pol_label)
      ALLOCATE(record%pol_amp(n_pol)%eps_sph(LBOUND(eps_sph, 1):UBOUND(eps_sph, 1)))
      record%pol_amp(n_pol)%eps_sph = eps_sph
      ! XAS uses sum_mj |M(m_j)|^2. RIXS needs the complex M(m_j)
      ! amplitudes and their physical twice_mj mapping before squaring.
      CALL record%pol_amp(n_pol)%amp%set_from_matrix_row(matrix, band, twice_mj=record%core%twice_mj)
   END SUBROUTINE xas_attach_polarization_amplitude

   SUBROUTINE xas_attach_l_channel_amplitudes(record, final_l, matrix_lchan, band, pol_index)
      TYPE(t_xas_transition_record), INTENT(INOUT) :: record
      INTEGER,                       INTENT(IN)    :: final_l(:), band
      COMPLEX,                       INTENT(IN)    :: matrix_lchan(:, :, :)
      INTEGER, OPTIONAL,             INTENT(IN)    :: pol_index

      INTEGER :: i_l, i_pol, n_channels

      IF (.NOT. ALLOCATED(record%pol_amp)) THEN
         CALL juDFT_error("XAS transition record has no polarization amplitudes for l-channel storage", &
                          calledby="m_xas_transitions")
      END IF
      IF (.NOT. ALLOCATED(record%core%twice_mj)) THEN
         CALL juDFT_error("XAS transition record has no core twice_mj mapping", calledby="m_xas_transitions")
      END IF
      IF (band /= record%band) THEN
         CALL juDFT_error("XAS transition l-channel band mismatch", calledby="m_xas_transitions")
      END IF
      IF (SIZE(matrix_lchan, 1) < band .OR. SIZE(matrix_lchan, 3) < SIZE(final_l)) THEN
         CALL juDFT_error("XAS l-channel matrix dimensions are too small", calledby="m_xas_transitions")
      END IF

      i_pol = 1
      IF (PRESENT(pol_index)) i_pol = pol_index
      IF (i_pol < 1 .OR. i_pol > SIZE(record%pol_amp)) THEN
         CALL juDFT_error("Invalid XAS transition polarization index for l channels", calledby="m_xas_transitions")
      END IF

      n_channels = SIZE(final_l)
      IF (ALLOCATED(record%pol_amp(i_pol)%final_l)) DEALLOCATE(record%pol_amp(i_pol)%final_l)
      IF (ALLOCATED(record%pol_amp(i_pol)%amp_l)) DEALLOCATE(record%pol_amp(i_pol)%amp_l)
      ALLOCATE(record%pol_amp(i_pol)%final_l(n_channels))
      ALLOCATE(record%pol_amp(i_pol)%amp_l(n_channels))
      record%pol_amp(i_pol)%final_l = final_l

      DO i_l = 1, n_channels
         CALL record%pol_amp(i_pol)%amp_l(i_l)%set_from_matrix_row(matrix_lchan(:, :, i_l), band, &
                                                                    twice_mj=record%core%twice_mj)
      END DO
   END SUBROUTINE xas_attach_l_channel_amplitudes

   SUBROUTINE xas_allowed_final_l_channels(core, lmax, n_channels, final_l)
      TYPE(t_xas_core_descriptor), INTENT(IN)  :: core
      INTEGER,                     INTENT(IN)  :: lmax
      INTEGER,                     INTENT(OUT) :: n_channels
      INTEGER,                     INTENT(OUT) :: final_l(:)

      n_channels = 0
      final_l = -1
      IF (core%lc - 1 >= 0 .AND. core%lc - 1 <= lmax) THEN
         n_channels = n_channels + 1
         IF (n_channels <= SIZE(final_l)) final_l(n_channels) = core%lc - 1
      END IF
      IF (core%lc + 1 <= lmax) THEN
         n_channels = n_channels + 1
         IF (n_channels <= SIZE(final_l)) final_l(n_channels) = core%lc + 1
      END IF
      IF (n_channels > SIZE(final_l)) THEN
         CALL juDFT_error("Too many XAS final-l channels", calledby="m_xas_transitions")
      END IF
   END SUBROUTINE xas_allowed_final_l_channels

   REAL FUNCTION xas_transition_absM2(record, pol_index) RESULT(abs_m2)
      TYPE(t_xas_transition_record), INTENT(IN) :: record
      INTEGER, OPTIONAL,             INTENT(IN) :: pol_index

      INTEGER :: i_pol

      IF (.NOT. ALLOCATED(record%pol_amp)) THEN
         CALL juDFT_error("XAS transition record has no polarization amplitudes", calledby="m_xas_transitions")
      END IF
      i_pol = 1
      IF (PRESENT(pol_index)) i_pol = pol_index
      IF (i_pol < 1 .OR. i_pol > SIZE(record%pol_amp)) THEN
         CALL juDFT_error("Invalid XAS transition polarization index", calledby="m_xas_transitions")
      END IF
      abs_m2 = record%pol_amp(i_pol)%amp%absM2
   END FUNCTION xas_transition_absM2

   REAL FUNCTION xas_transition_absM2_l(record, i_l, pol_index) RESULT(abs_m2)
      TYPE(t_xas_transition_record), INTENT(IN) :: record
      INTEGER,                       INTENT(IN) :: i_l
      INTEGER, OPTIONAL,             INTENT(IN) :: pol_index

      INTEGER :: i_pol

      i_pol = 1
      IF (PRESENT(pol_index)) i_pol = pol_index
      IF (.NOT. ALLOCATED(record%pol_amp)) THEN
         CALL juDFT_error("XAS transition record has no polarization amplitudes", calledby="m_xas_transitions")
      END IF
      IF (i_pol < 1 .OR. i_pol > SIZE(record%pol_amp)) THEN
         CALL juDFT_error("Invalid XAS transition polarization index", calledby="m_xas_transitions")
      END IF
      IF (.NOT. ALLOCATED(record%pol_amp(i_pol)%amp_l)) THEN
         CALL juDFT_error("XAS transition record has no l-channel amplitudes", calledby="m_xas_transitions")
      END IF
      IF (i_l < 1 .OR. i_l > SIZE(record%pol_amp(i_pol)%amp_l)) THEN
         CALL juDFT_error("Invalid XAS transition l-channel index", calledby="m_xas_transitions")
      END IF
      abs_m2 = record%pol_amp(i_pol)%amp_l(i_l)%absM2
   END FUNCTION xas_transition_absM2_l

   REAL FUNCTION xas_transition_l_reconstruction_error(record, pol_index) RESULT(max_error)
      TYPE(t_xas_transition_record), INTENT(IN) :: record
      INTEGER, OPTIONAL,             INTENT(IN) :: pol_index

      INTEGER :: i_pol, i_l, i_mj
      COMPLEX, ALLOCATABLE :: m_sum(:)

      i_pol = 1
      IF (PRESENT(pol_index)) i_pol = pol_index
      IF (.NOT. ALLOCATED(record%pol_amp)) THEN
         CALL juDFT_error("XAS transition record has no polarization amplitudes", calledby="m_xas_transitions")
      END IF
      IF (i_pol < 1 .OR. i_pol > SIZE(record%pol_amp)) THEN
         CALL juDFT_error("Invalid XAS transition polarization index", calledby="m_xas_transitions")
      END IF
      IF (.NOT. ALLOCATED(record%pol_amp(i_pol)%amp_l)) THEN
         CALL juDFT_error("XAS transition record has no l-channel amplitudes", calledby="m_xas_transitions")
      END IF

      ALLOCATE(m_sum(record%pol_amp(i_pol)%amp%n_mj), SOURCE=CMPLX(0.0, 0.0))
      DO i_l = 1, SIZE(record%pol_amp(i_pol)%amp_l)
         IF (record%pol_amp(i_pol)%amp_l(i_l)%n_mj /= SIZE(m_sum)) THEN
            CALL juDFT_error("XAS l-channel amplitude dimension mismatch", calledby="m_xas_transitions")
         END IF
         m_sum = m_sum + record%pol_amp(i_pol)%amp_l(i_l)%m
      END DO

      max_error = 0.0
      DO i_mj = 1, SIZE(m_sum)
         max_error = MAX(max_error, ABS(record%pol_amp(i_pol)%amp%m(i_mj) - m_sum(i_mj)))
      END DO
   END FUNCTION xas_transition_l_reconstruction_error

   REAL FUNCTION xas_transition_total_weighted_strength(record, pol_index) RESULT(strength)
      TYPE(t_xas_transition_record), INTENT(IN) :: record
      INTEGER, OPTIONAL,             INTENT(IN) :: pol_index

      strength = record%k_weight*record%one_minus_occ*xas_transition_absM2(record, pol_index)
   END FUNCTION xas_transition_total_weighted_strength

   INTEGER FUNCTION xas_edge_id(edge_name) RESULT(edge_id)
      CHARACTER(LEN=*), INTENT(IN) :: edge_name

      SELECT CASE (TRIM(ADJUSTL(edge_name)))
      CASE ("K", "k", "1s1/2", "1S1/2")
         edge_id = 1
      CASE ("L2", "l2", "2p1/2", "2P1/2")
         edge_id = 2
      CASE ("L3", "l3", "2p3/2", "2P3/2")
         edge_id = 3
      CASE DEFAULT
         edge_id = 0
      END SELECT
   END FUNCTION xas_edge_id

END MODULE m_xas_transitions
