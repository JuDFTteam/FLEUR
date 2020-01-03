!--------------------------------------------------------------------------------
! Copyright (c) 2016 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------

MODULE m_types_mpinp
   USE m_judft
   USE m_types_fleurinput_base
   IMPLICIT NONE
   PRIVATE

   TYPE, EXTENDS(t_fleurinput_base):: t_mpinp
      real                   :: g_cutoff = 0.0
      real                   :: linear_dep_tol = 0.0
   CONTAINS
      PROCEDURE :: read_xml => read_xml_mpinp
      PROCEDURE :: mpi_bc => mpi_bc_mpinp
   END TYPE t_mpinp
   PUBLIC t_mpinp

CONTAINS

   SUBROUTINE mpi_bc_mpinp(this, mpi_comm, irank)
      USE m_mpi_bc_tool
      CLASS(t_mpinp), INTENT(INOUT)::this
      INTEGER, INTENT(IN):: mpi_comm
      INTEGER, INTENT(IN), OPTIONAL::irank
      INTEGER ::rank
      IF (PRESENT(irank)) THEN
         rank = irank
      ELSE
         rank = 0
      END IF
      CALL mpi_bc(this%g_cutoff,rank,mpi_comm)
      CALL mpi_bc(this%linear_dep_tol,rank,mpi_comm)
   END SUBROUTINE mpi_bc_mpinp

   SUBROUTINE read_xml_mpinp(this, xml)
      USE m_types_xml
      CLASS(t_mpinp), INTENT(INout):: this
      TYPE(t_xml), INTENT(in)     :: xml

      INTEGER::numberNodes, ntype, itype
      CHARACTER(len=100)  :: xPathA
      CHARACTER(len=4), allocatable  :: xc_name

      ntype = xml%GetNumberOfNodes('/fleurInput/atomGroups/atomGroup')
      numberNodes = xml%GetNumberOfNodes('/fleurInput/calculationSetup/prodBasis')
      IF (numberNodes == 1) THEN
         this%g_cutoff=evaluateFirstOnly(xml%GetAttributeValue('/fleurInput/calculationSetup/prodBasis/@gcutm'))
         this%linear_dep_tol=evaluateFirstOnly(xml%GetAttributeValue('/fleurInput/calculationSetup/prodBasis/@tolerance'))
      ENDIF
   END SUBROUTINE read_xml_mpinp
END MODULE m_types_mpinp
