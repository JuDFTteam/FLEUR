!--------------------------------------------------------------------------------
! Copyright (c) 2020 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------

MODULE m_types_juPhon
   USE m_judft
   USE m_types_fleurinput_base
   IMPLICIT NONE
   PRIVATE

   TYPE, EXTENDS(t_fleurinput_base) :: t_juPhon
      LOGICAL :: l_jpCheck = .FALSE.
      LOGICAL :: l_jpTest  = .FALSE.
      LOGICAL :: l_potout  = .FALSE.
      LOGICAL :: l_eigout  = .FALSE.
      LOGICAL :: l_symTsh  = .FALSE.
      LOGICAL :: l_symTdm  = .FALSE.
      LOGICAL :: l_bfkq    = .FALSE.
      INTEGER :: jplmax    = 0
      REAL    :: kgqmax    = 0.0
      REAL    :: gqmax     = 0.0
      REAL    :: eps_pert  = 0.00001
      REAL    :: eDiffcut  = 1e-12
      REAL    :: qpt_ph(3)

   CONTAINS
      PROCEDURE :: read_xml => read_xml_juPhon
      PROCEDURE :: mpi_bc => mpi_bc_juPhon
   END TYPE t_juPhon

   PUBLIC t_juPhon

CONTAINS

   SUBROUTINE mpi_bc_juPhon(this, mpi_comm, irank)
      USE m_mpi_bc_tool

      CLASS(t_juPhon), INTENT(INOUT) ::this
      INTEGER, INTENT(IN)           :: mpi_comm
      INTEGER, INTENT(IN), OPTIONAL :: irank

      INTEGER :: rank

      IF (PRESENT(irank)) THEN
         rank = irank
      ELSE
         rank = 0
      END IF

      CALL mpi_bc(this%l_potout, rank, mpi_comm)
      CALL mpi_bc(this%l_eigout, rank, mpi_comm)

   END SUBROUTINE mpi_bc_juPhon

   SUBROUTINE read_xml_juPhon(this, xml)
      USE m_types_xml

      CLASS(t_juPhon), INTENT(INOUT) :: this
      TYPE(t_xml), INTENT(INOUT)     :: xml

      INTEGER::numberNodes
      CHARACTER(len=100) :: xPathA

      numberNodes = xml%GetNumberOfNodes('/fleurInput/output/juPhon')
      IF (numberNodes == 1) THEN
         this%l_jpCheck = evaluateFirstBoolOnly(xml%GetAttributeValue('/fleurInput/output/juPhon/@l_jpCheck'))
         this%l_jpTest  = evaluateFirstBoolOnly(xml%GetAttributeValue('/fleurInput/output/juPhon/@l_jpTest'))
         this%l_potout  = evaluateFirstBoolOnly(xml%GetAttributeValue('/fleurInput/output/juPhon/@l_potout'))
         this%l_eigout  = evaluateFirstBoolOnly(xml%GetAttributeValue('/fleurInput/output/juPhon/@l_eigout'))
         this%l_symTsh  = evaluateFirstBoolOnly(xml%GetAttributeValue('/fleurInput/output/juPhon/@l_symTsh'))
         this%l_symTdm  = evaluateFirstBoolOnly(xml%GetAttributeValue('/fleurInput/output/juPhon/@l_symTdm'))
         this%l_bfkq    = evaluateFirstBoolOnly(xml%GetAttributeValue('/fleurInput/output/juPhon/@l_bfkq'))
         this%jplmax    = evaluateFirstIntOnly(xml%GetAttributeValue('/fleurInput/output/juPhon/@l_jplmax'))
         this%kgqmax    = evaluateFirstOnly(xml%GetAttributeValue('/fleurInput/output/juPhon/@kgqmax'))
         this%gqmax     = evaluateFirstOnly(xml%GetAttributeValue('/fleurInput/output/juPhon/@gqmax'))
         this%eps_pert  = evaluateFirstOnly(xml%GetAttributeValue('/fleurInput/output/juPhon/@eps_pert'))
         this%eDiffcut  = evaluateFirstOnly(xml%GetAttributeValue('/fleurInput/output/juPhon/@eDiffcut'))
         this%qpt_ph(1) = 0.0
         this%qpt_ph(2) = 0.0
         this%qpt_ph(3) = 0.0    
      ENDIF

   END SUBROUTINE read_xml_juPhon

END MODULE m_types_juPhon
