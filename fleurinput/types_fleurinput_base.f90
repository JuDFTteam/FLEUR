!--------------------------------------------------------------------------------
! Copyright (c) 2016 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------

MODULE m_types_fleurinput_base
  IMPLICIT NONE
  PRIVATE
  !This module defines an abstract datatype all fleurinput-datatypes should
  !implement

  TYPE,ABSTRACT:: t_fleurinput_base
   CONTAINS
     PROCEDURE :: read_xml
     PROCEDURE :: mpi_bc
  END TYPE t_fleurinput_base
  PUBLIC t_fleurinput_base

CONTAINS
  SUBROUTINE read_xml(this,xml)
    USE m_types_xml
    CLASS(t_fleurinput_base),INTENT(INOUT):: this
    TYPE(t_xml),INTENT(IN)              :: xml
  END SUBROUTINE read_xml
  SUBROUTINE mpi_bc(this,mpi_comm,irank)
    USE m_types_xml
    CLASS(t_fleurinput_base),INTENT(INOUT):: this
    INTEGER,INTENT(IN)                    :: mpi_comm
    INTEGER,INTENT(IN),OPTIONAL::irank
  END SUBROUTINE mpi_bc
  
END MODULE m_types_fleurinput_base

