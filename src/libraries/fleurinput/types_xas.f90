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

      number_nodes = xml%GetNumberOfNodes('/fleurInput/output/xas')
      IF (number_nodes == 0) RETURN
      IF (number_nodes > 1) CALL juDFT_error("Only one output/xas section is allowed.", calledby="m_types_xas")

      this%l_xas = .TRUE.
      IF (xml%GetNumberOfNodes('/fleurInput/output/xas/@l_xas') == 1) THEN
         this%l_xas = evaluateFirstBoolOnly(xml%GetAttributeValue('/fleurInput/output/xas/@l_xas'))
      END IF
      IF (.NOT. this%l_xas) RETURN

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
      CALL xas_validate_edge_name(this%edge)
   END SUBROUTINE read_xml_xas

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
         CALL juDFT_error("Unsupported XAS polarization. Use linear x, y, z for now.", calledby="m_types_xas")
      END SELECT
   END SUBROUTINE xas_set_polarization

   SUBROUTINE xas_validate_edge_name(edge)
      CHARACTER(LEN=*), INTENT(IN) :: edge

      SELECT CASE (TRIM(ADJUSTL(edge)))
      CASE ("K", "k", "1s1/2", "1S1/2", "L2", "l2", "2p1/2", "2P1/2", "L3", "l3", "2p3/2", "2P3/2")
         RETURN
      CASE ("L23", "l23", "L2,3", "l2,3")
         CALL juDFT_error("Combined L23 edge is not implemented yet; run L2 and L3 separately.", calledby="m_types_xas")
      CASE DEFAULT
         CALL juDFT_error("Unsupported XAS edge. Supported edges are K, L2, and L3.", calledby="m_types_xas")
      END SELECT
   END SUBROUTINE xas_validate_edge_name

END MODULE m_types_xas
