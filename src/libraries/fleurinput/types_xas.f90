!--------------------------------------------------------------------------------
! Copyright (c) 2026 Peter Gruenberg Institut, Forschungszentrum Juelich, Germany
! This file is part of FLEUR and available as free software under the conditions
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------

MODULE m_types_xas
   USE m_juDFT
   USE m_types_fleurinput_base
   IMPLICIT NONE
   PRIVATE
   PUBLIC :: t_xas

   TYPE, EXTENDS(t_fleurinput_base) :: t_xas
      LOGICAL :: l_xas = .FALSE.
      ! Chemical selection for the absorber set. All atom types with this
      ! nuclear charge are summed; site-resolved selection is a later extension.
      INTEGER :: absorber_z = -1
      ! Supported production edges: K, L2, L3. Combined L23 is intentionally
      ! rejected so that L2 and L3 can be run and validated independently.
      CHARACTER(LEN=16) :: edge = "L3"
      REAL :: eta = 0.03
      ! If no explicit eMin/eMax is provided, the driver scans the available
      ! empty-state transitions and prints the resolved grid window at runtime.
      LOGICAL :: l_energy_window = .FALSE.
      REAL :: e_min = 0.0
      REAL :: e_max = 0.0
      INTEGER :: n_energy = 401
      LOGICAL :: polarizations(3) = [.TRUE., .TRUE., .TRUE.]
      CHARACTER(LEN=64) :: output_prefix = "xas"
      LOGICAL :: write_transitions = .FALSE.
      ! Independent-particle RIXS input shares this type with XAS because both
      ! paths use the same core-state and matrix-element postprocessing
      ! infrastructure.
      LOGICAL :: l_rixs = .FALSE.
      INTEGER :: rixs_absorber_z = -1
      CHARACTER(LEN=16) :: rixs_edge = "L3"
      REAL :: rixs_omega_in = 0.0
      REAL :: rixs_gamma_core = 0.0
      REAL :: rixs_loss_min = 0.0
      REAL :: rixs_loss_max = 0.0
      INTEGER :: rixs_n_loss = 0
      REAL :: rixs_eta_loss = 0.0
      LOGICAL :: rixs_in_polarizations(3) = [.FALSE., .FALSE., .FALSE.]
      LOGICAL :: rixs_out_polarizations(3) = [.FALSE., .FALSE., .FALSE.]
      CHARACTER(LEN=64) :: rixs_output_prefix = "rixs"
      LOGICAL :: rixs_write_contributions = .FALSE.
      LOGICAL :: rixs_write_state_character = .FALSE.
      INTEGER :: rixs_state_ligand_z = -1
      LOGICAL :: l_rixs_valence_band_min = .FALSE.
      LOGICAL :: l_rixs_valence_band_max = .FALSE.
      INTEGER :: rixs_valence_band_min = 1
      INTEGER :: rixs_valence_band_max = 1
      LOGICAL :: l_rixs_intermediate_band_min = .FALSE.
      LOGICAL :: l_rixs_intermediate_band_max = .FALSE.
      INTEGER :: rixs_intermediate_band_min = 1
      INTEGER :: rixs_intermediate_band_max = 1
   CONTAINS
      PROCEDURE :: read_xml => read_xml_xas
      PROCEDURE :: mpi_bc => mpi_bc_xas
   END TYPE t_xas

CONTAINS

   SUBROUTINE mpi_bc_xas(this, mpi_comm, irank)
      USE m_mpi_bc_tool
      CLASS(t_xas), INTENT(INOUT) :: this
      INTEGER,      INTENT(IN)    :: mpi_comm
      INTEGER,      INTENT(IN), OPTIONAL :: irank

      INTEGER :: rank

      IF (PRESENT(irank)) THEN
         rank = irank
      ELSE
         rank = 0
      END IF

      CALL mpi_bc(this%l_xas, rank, mpi_comm)
      CALL mpi_bc(this%absorber_z, rank, mpi_comm)
      CALL mpi_bc(rank, mpi_comm, this%edge)
      CALL mpi_bc(this%eta, rank, mpi_comm)
      CALL mpi_bc(this%l_energy_window, rank, mpi_comm)
      CALL mpi_bc(this%e_min, rank, mpi_comm)
      CALL mpi_bc(this%e_max, rank, mpi_comm)
      CALL mpi_bc(this%n_energy, rank, mpi_comm)
      CALL mpi_bc(this%polarizations(1), rank, mpi_comm)
      CALL mpi_bc(this%polarizations(2), rank, mpi_comm)
      CALL mpi_bc(this%polarizations(3), rank, mpi_comm)
      CALL mpi_bc(rank, mpi_comm, this%output_prefix)
      CALL mpi_bc(this%write_transitions, rank, mpi_comm)
      CALL mpi_bc(this%l_rixs, rank, mpi_comm)
      CALL mpi_bc(this%rixs_absorber_z, rank, mpi_comm)
      CALL mpi_bc(rank, mpi_comm, this%rixs_edge)
      CALL mpi_bc(this%rixs_omega_in, rank, mpi_comm)
      CALL mpi_bc(this%rixs_gamma_core, rank, mpi_comm)
      CALL mpi_bc(this%rixs_loss_min, rank, mpi_comm)
      CALL mpi_bc(this%rixs_loss_max, rank, mpi_comm)
      CALL mpi_bc(this%rixs_n_loss, rank, mpi_comm)
      CALL mpi_bc(this%rixs_eta_loss, rank, mpi_comm)
      CALL mpi_bc(this%rixs_in_polarizations(1), rank, mpi_comm)
      CALL mpi_bc(this%rixs_in_polarizations(2), rank, mpi_comm)
      CALL mpi_bc(this%rixs_in_polarizations(3), rank, mpi_comm)
      CALL mpi_bc(this%rixs_out_polarizations(1), rank, mpi_comm)
      CALL mpi_bc(this%rixs_out_polarizations(2), rank, mpi_comm)
      CALL mpi_bc(this%rixs_out_polarizations(3), rank, mpi_comm)
      CALL mpi_bc(rank, mpi_comm, this%rixs_output_prefix)
      CALL mpi_bc(this%rixs_write_contributions, rank, mpi_comm)
      CALL mpi_bc(this%rixs_write_state_character, rank, mpi_comm)
      CALL mpi_bc(this%rixs_state_ligand_z, rank, mpi_comm)
      CALL mpi_bc(this%l_rixs_valence_band_min, rank, mpi_comm)
      CALL mpi_bc(this%l_rixs_valence_band_max, rank, mpi_comm)
      CALL mpi_bc(this%rixs_valence_band_min, rank, mpi_comm)
      CALL mpi_bc(this%rixs_valence_band_max, rank, mpi_comm)
      CALL mpi_bc(this%l_rixs_intermediate_band_min, rank, mpi_comm)
      CALL mpi_bc(this%l_rixs_intermediate_band_max, rank, mpi_comm)
      CALL mpi_bc(this%rixs_intermediate_band_min, rank, mpi_comm)
      CALL mpi_bc(this%rixs_intermediate_band_max, rank, mpi_comm)
   END SUBROUTINE mpi_bc_xas

   SUBROUTINE read_xml_xas(this, xml)
      USE m_types_xml
      CLASS(t_xas), INTENT(INOUT) :: this
      TYPE(t_xml),  INTENT(INOUT) :: xml

      CHARACTER(LEN=255) :: value_string
      INTEGER :: number_nodes, num_tokens, i_token
      LOGICAL :: has_e_min, has_e_max

      this%l_xas = .FALSE.
      this%absorber_z = -1
      this%edge = "L3"
      this%eta = 0.03
      this%l_energy_window = .FALSE.
      this%e_min = 0.0
      this%e_max = 0.0
      this%n_energy = 401
      this%polarizations = [.TRUE., .TRUE., .TRUE.]
      this%output_prefix = "xas"
      this%write_transitions = .FALSE.
      CALL rixs_reset_defaults(this)

      number_nodes = xml%GetNumberOfNodes('/fleurInput/output/xas')
      IF (number_nodes > 1) CALL juDFT_error("Only one output/xas section is allowed.", calledby="m_types_xas")
      IF (number_nodes == 1) THEN

         this%l_xas = .TRUE.
         IF (xml%GetNumberOfNodes('/fleurInput/output/xas/@l_xas') == 1) THEN
            this%l_xas = evaluateFirstBoolOnly(xml%GetAttributeValue('/fleurInput/output/xas/@l_xas'))
         END IF
         IF (this%l_xas) THEN

            IF (xml%GetNumberOfNodes('/fleurInput/output/xas/@absorberZ') == 1) THEN
               this%absorber_z = evaluateFirstIntOnly(xml%GetAttributeValue('/fleurInput/output/xas/@absorberZ'))
            ELSE
               CALL juDFT_error("XAS is enabled but output/xas/@absorberZ is missing.", calledby="m_types_xas")
            END IF

            IF (xml%GetNumberOfNodes('/fleurInput/output/xas/@edge') == 1) THEN
               this%edge = TRIM(ADJUSTL(xml%GetAttributeValue('/fleurInput/output/xas/@edge')))
            END IF
            IF (xml%GetNumberOfNodes('/fleurInput/output/xas/@eta') == 1) THEN
               this%eta = evaluateFirstOnly(xml%GetAttributeValue('/fleurInput/output/xas/@eta'))
            END IF
            has_e_min = xml%GetNumberOfNodes('/fleurInput/output/xas/@eMin') == 1
            has_e_max = xml%GetNumberOfNodes('/fleurInput/output/xas/@eMax') == 1
            IF (has_e_min .NEQV. has_e_max) THEN
               CALL juDFT_error("XAS energy window requires both eMin and eMax or neither.", calledby="m_types_xas")
            END IF
            IF (has_e_min) THEN
               this%e_min = evaluateFirstOnly(xml%GetAttributeValue('/fleurInput/output/xas/@eMin'))
               this%e_max = evaluateFirstOnly(xml%GetAttributeValue('/fleurInput/output/xas/@eMax'))
               this%l_energy_window = .TRUE.
               IF (this%e_min >= this%e_max) CALL juDFT_error("XAS eMin must be smaller than eMax.", calledby="m_types_xas")
            END IF
            IF (xml%GetNumberOfNodes('/fleurInput/output/xas/@nEnergy') == 1) THEN
               this%n_energy = evaluateFirstIntOnly(xml%GetAttributeValue('/fleurInput/output/xas/@nEnergy'))
            END IF
            IF (this%n_energy < 2) CALL juDFT_error("XAS nEnergy must be at least 2.", calledby="m_types_xas")

            IF (xml%GetNumberOfNodes('/fleurInput/output/xas/@polarizations') == 1) THEN
               value_string = TRIM(ADJUSTL(xml%GetAttributeValue('/fleurInput/output/xas/@polarizations')))
               num_tokens = xml%countStringTokens(value_string)
               IF (num_tokens < 1) CALL juDFT_error("XAS polarizations list is empty.", calledby="m_types_xas")
               this%polarizations = .FALSE.
               DO i_token = 1, num_tokens
                  CALL xas_set_polarization(this%polarizations, xml%popFirstStringToken(value_string))
               END DO
            END IF
            IF (.NOT. ANY(this%polarizations)) CALL juDFT_error("No valid linear XAS polarizations selected.", calledby="m_types_xas")

            IF (xml%GetNumberOfNodes('/fleurInput/output/xas/@outputPrefix') == 1) THEN
               this%output_prefix = TRIM(ADJUSTL(xml%GetAttributeValue('/fleurInput/output/xas/@outputPrefix')))
            END IF
            IF (LEN_TRIM(this%output_prefix) == 0) CALL juDFT_error("XAS outputPrefix must not be empty.", calledby="m_types_xas")
            IF (xml%GetNumberOfNodes('/fleurInput/output/xas/@writeTransitions') == 1) THEN
               this%write_transitions = evaluateFirstBoolOnly(xml%GetAttributeValue('/fleurInput/output/xas/@writeTransitions'))
            END IF
            CALL xas_validate_edge_name(this%edge, "XAS")
         END IF
      END IF
      CALL read_xml_rixs(this, xml)
   END SUBROUTINE read_xml_xas

   SUBROUTINE rixs_reset_defaults(this)
      CLASS(t_xas), INTENT(INOUT) :: this

      this%l_rixs = .FALSE.
      this%rixs_absorber_z = -1
      this%rixs_edge = "L3"
      this%rixs_omega_in = 0.0
      this%rixs_gamma_core = 0.0
      this%rixs_loss_min = 0.0
      this%rixs_loss_max = 0.0
      this%rixs_n_loss = 0
      this%rixs_eta_loss = 0.0
      this%rixs_in_polarizations = [.FALSE., .FALSE., .FALSE.]
      this%rixs_out_polarizations = [.FALSE., .FALSE., .FALSE.]
      this%rixs_output_prefix = "rixs"
      this%rixs_write_contributions = .FALSE.
      this%rixs_write_state_character = .FALSE.
      this%rixs_state_ligand_z = -1
      this%l_rixs_valence_band_min = .FALSE.
      this%l_rixs_valence_band_max = .FALSE.
      this%rixs_valence_band_min = 1
      this%rixs_valence_band_max = 1
      this%l_rixs_intermediate_band_min = .FALSE.
      this%l_rixs_intermediate_band_max = .FALSE.
      this%rixs_intermediate_band_min = 1
      this%rixs_intermediate_band_max = 1
   END SUBROUTINE rixs_reset_defaults

   SUBROUTINE read_xml_rixs(this, xml)
      USE m_types_xml
      CLASS(t_xas), INTENT(INOUT) :: this
      TYPE(t_xml),  INTENT(INOUT) :: xml

      CHARACTER(LEN=255) :: value_string
      INTEGER :: number_nodes, num_tokens, i_token

      number_nodes = xml%GetNumberOfNodes('/fleurInput/output/rixs')
      IF (number_nodes == 0) RETURN
      IF (number_nodes > 1) CALL juDFT_error("Only one output/rixs section is allowed.", calledby="m_types_xas")

      this%l_rixs = .TRUE.
      IF (xml%GetNumberOfNodes('/fleurInput/output/rixs/@l_rixs') == 1) THEN
         this%l_rixs = evaluateFirstBoolOnly(xml%GetAttributeValue('/fleurInput/output/rixs/@l_rixs'))
      END IF
      IF (.NOT. this%l_rixs) RETURN

      IF (xml%GetNumberOfNodes('/fleurInput/output/rixs/@absorberZ') == 1) THEN
         this%rixs_absorber_z = evaluateFirstIntOnly(xml%GetAttributeValue('/fleurInput/output/rixs/@absorberZ'))
      ELSE
         CALL juDFT_error("RIXS is enabled but output/rixs/@absorberZ is missing.", calledby="m_types_xas")
      END IF
      IF (xml%GetNumberOfNodes('/fleurInput/output/rixs/@edge') == 1) THEN
         this%rixs_edge = TRIM(ADJUSTL(xml%GetAttributeValue('/fleurInput/output/rixs/@edge')))
      ELSE
         CALL juDFT_error("RIXS is enabled but output/rixs/@edge is missing.", calledby="m_types_xas")
      END IF
      CALL xas_validate_edge_name(this%rixs_edge, "RIXS")
      IF (xml%GetNumberOfNodes('/fleurInput/output/rixs/@omegaIn') == 1) THEN
         this%rixs_omega_in = evaluateFirstOnly(xml%GetAttributeValue('/fleurInput/output/rixs/@omegaIn'))
      ELSE
         CALL juDFT_error("RIXS is enabled but output/rixs/@omegaIn is missing.", calledby="m_types_xas")
      END IF
      IF (xml%GetNumberOfNodes('/fleurInput/output/rixs/@gammaCore') == 1) THEN
         this%rixs_gamma_core = evaluateFirstOnly(xml%GetAttributeValue('/fleurInput/output/rixs/@gammaCore'))
      ELSE
         CALL juDFT_error("RIXS is enabled but output/rixs/@gammaCore is missing.", calledby="m_types_xas")
      END IF
      IF (xml%GetNumberOfNodes('/fleurInput/output/rixs/@lossMin') == 1) THEN
         this%rixs_loss_min = evaluateFirstOnly(xml%GetAttributeValue('/fleurInput/output/rixs/@lossMin'))
      ELSE
         CALL juDFT_error("RIXS is enabled but output/rixs/@lossMin is missing.", calledby="m_types_xas")
      END IF
      IF (xml%GetNumberOfNodes('/fleurInput/output/rixs/@lossMax') == 1) THEN
         this%rixs_loss_max = evaluateFirstOnly(xml%GetAttributeValue('/fleurInput/output/rixs/@lossMax'))
      ELSE
         CALL juDFT_error("RIXS is enabled but output/rixs/@lossMax is missing.", calledby="m_types_xas")
      END IF
      IF (this%rixs_loss_min >= this%rixs_loss_max) CALL juDFT_error("RIXS lossMin must be smaller than lossMax.", calledby="m_types_xas")
      IF (xml%GetNumberOfNodes('/fleurInput/output/rixs/@nLoss') == 1) THEN
         this%rixs_n_loss = evaluateFirstIntOnly(xml%GetAttributeValue('/fleurInput/output/rixs/@nLoss'))
      ELSE
         CALL juDFT_error("RIXS is enabled but output/rixs/@nLoss is missing.", calledby="m_types_xas")
      END IF
      IF (this%rixs_n_loss < 2) CALL juDFT_error("RIXS nLoss must be at least 2.", calledby="m_types_xas")
      IF (xml%GetNumberOfNodes('/fleurInput/output/rixs/@etaLoss') == 1) THEN
         this%rixs_eta_loss = evaluateFirstOnly(xml%GetAttributeValue('/fleurInput/output/rixs/@etaLoss'))
      ELSE
         CALL juDFT_error("RIXS is enabled but output/rixs/@etaLoss is missing.", calledby="m_types_xas")
      END IF
      IF (this%rixs_gamma_core <= 0.0) CALL juDFT_error("RIXS gammaCore must be positive.", calledby="m_types_xas")
      IF (this%rixs_eta_loss <= 0.0) CALL juDFT_error("RIXS etaLoss must be positive.", calledby="m_types_xas")

      this%rixs_in_polarizations = .FALSE.
      IF (xml%GetNumberOfNodes('/fleurInput/output/rixs/@incomingPolarizations') == 1) THEN
         value_string = TRIM(ADJUSTL(xml%GetAttributeValue('/fleurInput/output/rixs/@incomingPolarizations')))
      ELSE
         value_string = "x"
      END IF
      num_tokens = xml%countStringTokens(value_string)
      IF (num_tokens < 1) CALL juDFT_error("RIXS incomingPolarizations list is empty.", calledby="m_types_xas")
      DO i_token = 1, num_tokens
         CALL xas_set_polarization(this%rixs_in_polarizations, xml%popFirstStringToken(value_string))
      END DO

      this%rixs_out_polarizations = .FALSE.
      IF (xml%GetNumberOfNodes('/fleurInput/output/rixs/@outgoingPolarizations') == 1) THEN
         value_string = TRIM(ADJUSTL(xml%GetAttributeValue('/fleurInput/output/rixs/@outgoingPolarizations')))
      ELSE
         value_string = "x"
      END IF
      num_tokens = xml%countStringTokens(value_string)
      IF (num_tokens < 1) CALL juDFT_error("RIXS outgoingPolarizations list is empty.", calledby="m_types_xas")
      DO i_token = 1, num_tokens
         CALL xas_set_polarization(this%rixs_out_polarizations, xml%popFirstStringToken(value_string))
      END DO
      IF (.NOT. ANY(this%rixs_in_polarizations)) CALL juDFT_error("No valid linear RIXS incoming polarizations selected.", calledby="m_types_xas")
      IF (.NOT. ANY(this%rixs_out_polarizations)) CALL juDFT_error("No valid linear RIXS outgoing polarizations selected.", calledby="m_types_xas")

      IF (xml%GetNumberOfNodes('/fleurInput/output/rixs/@outputPrefix') == 1) THEN
         this%rixs_output_prefix = TRIM(ADJUSTL(xml%GetAttributeValue('/fleurInput/output/rixs/@outputPrefix')))
      END IF
      IF (LEN_TRIM(this%rixs_output_prefix) == 0) CALL juDFT_error("RIXS outputPrefix must not be empty.", calledby="m_types_xas")
      IF (xml%GetNumberOfNodes('/fleurInput/output/rixs/@writeContributions') == 1) THEN
         this%rixs_write_contributions = evaluateFirstBoolOnly(xml%GetAttributeValue('/fleurInput/output/rixs/@writeContributions'))
      END IF
      IF (xml%GetNumberOfNodes('/fleurInput/output/rixs/@writeStateCharacter') == 1) THEN
         this%rixs_write_state_character = &
            evaluateFirstBoolOnly(xml%GetAttributeValue('/fleurInput/output/rixs/@writeStateCharacter'))
      END IF
      IF (xml%GetNumberOfNodes('/fleurInput/output/rixs/@stateLigandZ') == 1) THEN
         this%rixs_state_ligand_z = evaluateFirstIntOnly(xml%GetAttributeValue('/fleurInput/output/rixs/@stateLigandZ'))
      END IF
      IF (this%rixs_write_state_character .AND. this%rixs_state_ligand_z <= 0) THEN
         CALL juDFT_error("RIXS writeStateCharacter requires a positive stateLigandZ.",calledby="m_types_xas")
      END IF
      CALL rixs_read_band_windows(this, xml)
   END SUBROUTINE read_xml_rixs

   SUBROUTINE rixs_read_band_windows(this, xml)
      USE m_types_xml
      CLASS(t_xas), INTENT(INOUT) :: this
      TYPE(t_xml),  INTENT(INOUT) :: xml

      this%l_rixs_valence_band_min = xml%GetNumberOfNodes('/fleurInput/output/rixs/@valenceBandMin') == 1
      this%l_rixs_valence_band_max = xml%GetNumberOfNodes('/fleurInput/output/rixs/@valenceBandMax') == 1
      this%l_rixs_intermediate_band_min = xml%GetNumberOfNodes('/fleurInput/output/rixs/@intermediateBandMin') == 1
      this%l_rixs_intermediate_band_max = xml%GetNumberOfNodes('/fleurInput/output/rixs/@intermediateBandMax') == 1

      IF (this%l_rixs_valence_band_min) THEN
         this%rixs_valence_band_min = evaluateFirstIntOnly(xml%GetAttributeValue('/fleurInput/output/rixs/@valenceBandMin'))
      END IF
      IF (this%l_rixs_valence_band_max) THEN
         this%rixs_valence_band_max = evaluateFirstIntOnly(xml%GetAttributeValue('/fleurInput/output/rixs/@valenceBandMax'))
      END IF
      IF (this%l_rixs_intermediate_band_min) THEN
         this%rixs_intermediate_band_min = evaluateFirstIntOnly(xml%GetAttributeValue('/fleurInput/output/rixs/@intermediateBandMin'))
      END IF
      IF (this%l_rixs_intermediate_band_max) THEN
         this%rixs_intermediate_band_max = evaluateFirstIntOnly(xml%GetAttributeValue('/fleurInput/output/rixs/@intermediateBandMax'))
      END IF

      IF (this%l_rixs_valence_band_min) THEN
         IF (this%rixs_valence_band_min <= 0) CALL juDFT_error("RIXS valenceBandMin must be positive.", calledby="m_types_xas")
      END IF
      IF (this%l_rixs_valence_band_max) THEN
         IF (this%rixs_valence_band_max <= 0) CALL juDFT_error("RIXS valenceBandMax must be positive.", calledby="m_types_xas")
      END IF
      IF (this%l_rixs_intermediate_band_min) THEN
         IF (this%rixs_intermediate_band_min <= 0) CALL juDFT_error("RIXS intermediateBandMin must be positive.", calledby="m_types_xas")
      END IF
      IF (this%l_rixs_intermediate_band_max) THEN
         IF (this%rixs_intermediate_band_max <= 0) CALL juDFT_error("RIXS intermediateBandMax must be positive.", calledby="m_types_xas")
      END IF

      IF (this%l_rixs_valence_band_min .AND. this%l_rixs_valence_band_max) THEN
         IF (this%rixs_valence_band_min > this%rixs_valence_band_max) THEN
            CALL juDFT_error("RIXS valenceBandMin must not be larger than valenceBandMax.", calledby="m_types_xas")
         END IF
      END IF
      IF (this%l_rixs_intermediate_band_min .AND. this%l_rixs_intermediate_band_max) THEN
         IF (this%rixs_intermediate_band_min > this%rixs_intermediate_band_max) THEN
            CALL juDFT_error("RIXS intermediateBandMin must not be larger than intermediateBandMax.", calledby="m_types_xas")
         END IF
      END IF
   END SUBROUTINE rixs_read_band_windows

   SUBROUTINE xas_set_polarization(polarizations, token)
      LOGICAL,          INTENT(INOUT) :: polarizations(3)
      CHARACTER(LEN=*), INTENT(IN)    :: token

      SELECT CASE (TRIM(ADJUSTL(token)))
      CASE ("x", "X")
         polarizations(1) = .TRUE.
      CASE ("y", "Y")
         polarizations(2) = .TRUE.
      CASE ("z", "Z")
         polarizations(3) = .TRUE.
      CASE DEFAULT
         CALL juDFT_error("Unsupported XAS/RIXS polarization. Use linear x, y, z for now.", calledby="m_types_xas")
      END SELECT
   END SUBROUTINE xas_set_polarization

   SUBROUTINE xas_validate_edge_name(edge, context)
      CHARACTER(LEN=*), INTENT(IN) :: edge
      CHARACTER(LEN=*), OPTIONAL, INTENT(IN) :: context
      CHARACTER(LEN=16) :: label

      label = "XAS"
      IF (PRESENT(context)) label = TRIM(ADJUSTL(context))

      SELECT CASE (TRIM(ADJUSTL(edge)))
      CASE ("K", "k", "1s1/2", "1S1/2", "L2", "l2", "2p1/2", "2P1/2", "L3", "l3", "2p3/2", "2P3/2")
         RETURN
      CASE ("L23", "l23", "L2,3", "l2,3")
         CALL juDFT_error("Combined L23 edge is not implemented for "//TRIM(label)//"; run L2 and L3 separately.", calledby="m_types_xas")
      CASE DEFAULT
         CALL juDFT_error("Unsupported "//TRIM(label)//" edge. Supported edges are K, L2, and L3.", calledby="m_types_xas")
      END SELECT
   END SUBROUTINE xas_validate_edge_name

END MODULE m_types_xas
