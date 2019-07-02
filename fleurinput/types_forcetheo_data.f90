!--------------------------------------------------------------------------------
! Copyright (c) 2016 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------


!>This module defines the basic type for calculations of force-theorem type
!! Here only a dummy is defined that should be extended by a custom made data-type
!! The functionality is encoded into four functions/subroutines:
!! start: This routine is called in each SC-loop before the force theorem-loop
!! next_job: This function returns .true. if another job should be done, it also modifies its
!!           arguments appropriately to perform the calculation
!! eval: Here the calculation is done in this function, the results might be stored in
!!           MODULE variables, IF a .TRUE. is returned, the rest of the loop (charge generation)
!!           is skipped
!! postprocess: After the calculation here some IO etc. can be done
!!
!!
!! An example for a non-trivial force-theorem type extending this datatype can be found in
!! forcetheorem/mae.F90

MODULE m_types_forcetheo_data
  USE m_juDFT
  USE m_types_fleurinput_base
  IMPLICIT NONE
  PRIVATE
  PUBLIC:: t_forcetheo_data
  TYPE,EXTENDS(t_fleurinput_base) :: t_forcetheo_data
     INTEGER :: mode=0 !MAE=1,DMI=2,Jij=3,ssdisp=4
     REAL,ALLOCATABLE:: qvec(:,:) !DMI,Jij,ssdisp
     REAL,ALLOCATABLE:: theta(:)  !DMI,MAE,jij(1st only)
     REAL,ALLOCATABLE:: phi(:)    !MAE
    
   CONTAINS
     PROCEDURE :: read_xml=>read_xml_forcetheo_data
     PROCEDURE :: mpi_bc=>mpi_bc_forcetheo_data
  END TYPE t_forcetheo_data
CONTAINS
   SUBROUTINE mpi_bc_forcetheo_data(this,mpi_comm,irank)
    use m_mpi_bc_tool
    CLASS(t_forcetheo_data),INTENT(INOUT)::this
    integer,INTENT(IN):: mpi_comm
    INTEGER,INTENT(IN),OPTIONAL::irank
    INTEGER ::rank
    if (present(irank)) THEN
       rank=irank
    else
       rank=0
    end if
    CALL mpi_bc(this%mode,rank,mpi_comm)
    CALL mpi_bc(this%qvec,rank,mpi_comm)
    CALL mpi_bc(this%theta,rank,mpi_comm)
    CALL mpi_bc(this%phi ,rank,mpi_comm)

  END SUBROUTINE mpi_bc_forcetheo_data

  SUBROUTINE read_xml_forcetheo_data(this,xml)
    USE m_types_xml
    CLASS(t_forcetheo_data),INTENT(INOUT):: this
    TYPE(t_xml),INTENT(IN)             :: xml
    CHARACTER(len=200)::str
    
    IF (xml%GetNumberOfNodes('/fleurInput/forceTheorem/MAE')==1) THEN
       this%mode=1
       str=xml%GetAttributeValue('/fleurInput/forceTheorem/MAE/@theta')
       CALL evaluateList(this%theta,str)
       str=xml%GetAttributeValue('/fleurInput/forceTheorem/MAE/@phi')
       CALL evaluateList(this%phi,str)
    ENDIF
    IF (xml%GetNumberOfNodes('/fleurInput/forceTheorem/DMI')==1) THEN
       this%mode=2
       this%qvec=xml%read_q_list('/fleurInput/forceTheorem/DMI/qVectors')
       str=xml%GetAttributeValue('/fleurInput/forceTheorem/DMI/@theta')
       CALL evaluateList(this%theta,str)
       str=xml%GetAttributeValue('/fleurInput/forceTheorem/DMI/@phi')
       CALL evaluateList(this%phi,str)
    ENDIF
    IF (xml%GetNumberOfNodes('/fleurInput/forceTheorem/Jij')==1) THEN
       this%mode=3
       this%qvec=xml%read_q_list('/fleurInput/forceTheorem/Jij/qVectors')
       ALLOCATE(this%theta(1))
       this%theta(1)=evaluateFirstOnly(xml%GetAttributeValue('/fleurInput/forceTheorem/Jij/@thetaj'))      
    ENDIF
    
    IF (xml%GetNumberOfNodes('/fleurInput/forceTheorem/spinSpiralDispersion')==1) THEN
       this%mode=4
       this%qvec=xml%read_q_list('/fleurInput/forceTheorem/spinSpiralDispersion')
    ENDIF
  END SUBROUTINE read_xml_forcetheo_data

END MODULE m_types_forcetheo_data

