!--------------------------------------------------------------------------------
! Copyright (c) 2016 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------

!> This module defines the abstract type each of the FLEUR setup datatypes should implement
MODULE m_types_fleur_setup
  IMPLICIT NONE
  TYPE,ABSTRACT:: t_fleursetup
     INTEGER:: version
   CONTAINS
     PROCEDURE(broadcast_type),DEFERRED,PASS(tt) :: broadcast
     PROCEDURE(write_type),DEFERRED,PASS(tt)     :: write
     PROCEDURE(read_type),DEFERRED,PASS(tt)      :: read
     PROCEDURE(read_xml),DEFERRED,PASS(tt)       :: read_xml
     GENERIC :: WRITE(formatted)=>write
     GENERIC :: READ(formatted)=>read
  END TYPE t_fleursetup

  INTERFACE
     SUBROUTINE broadcast_type(tt,mpi_comm,origin)
       IMPORT t_fleursetup
       CLASS(t_fleursetup),INTENT(INOUT):: tt
       INTEGER,INTENT(IN)               :: mpi_comm
       INTEGER,INTENT(IN),OPTIONAL      :: origin
     END SUBROUTINE broadcast_type
  END INTERFACE

  INTERFACE
     SUBROUTINE read_xml(tt)
       IMPORT t_fleursetup
       CLASS(t_fleursetup),INTENT(OUT):: tt
     END SUBROUTINE read_xml
  END INTERFACE
  
  INTERFACE
     SUBROUTINE  write_type(tt, unit, iotype, v_list, iostat, iomsg)
       IMPORT t_fleursetup
       CLASS(t_fleursetup),INTENT(in)  :: tt
       INTEGER, INTENT(IN)             :: unit
       CHARACTER(*), INTENT(IN)        :: iotype
       INTEGER, INTENT(IN)             :: v_list(:)
       INTEGER, INTENT(OUT)            :: iostat
       CHARACTER(*), INTENT(INOUT)     :: iomsg
     END SUBROUTINE write_type
  END INTERFACE

  INTERFACE
     SUBROUTINE  read_type(tt, unit, iotype, v_list, iostat, iomsg)
       IMPORT t_fleursetup
       CLASS(t_fleursetup),INTENT(INOUT) :: tt
       INTEGER, INTENT(IN)             :: unit
       CHARACTER(*), INTENT(IN)        :: iotype
       INTEGER, INTENT(IN)             :: v_list(:)
       INTEGER, INTENT(OUT)            :: iostat
       CHARACTER(*), INTENT(INOUT)     :: iomsg
     END SUBROUTINE read_type
  END INTERFACE
END MODULE m_types_fleur_setup
