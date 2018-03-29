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

MODULE m_types_forcetheo
  TYPE :: t_forcetheo
     LOGICAL,PRIVATE :: firstloop
   CONTAINS
     PROCEDURE :: start   =>forcetheo_start
     PROCEDURE :: next_job=>forcetheo_next_job 
     PROCEDURE :: eval    =>forcetheo_eval
     PROCEDURE :: postprocess => forcetheo_postprocess
  END TYPE t_forcetheo

CONTAINS
  SUBROUTINE forcetheo_start(this)
    IMPLICIT NONE
    CLASS(t_forcetheo),INTENT(INOUT):: this
    this%firstloop=.TRUE.
  END SUBROUTINE forcetheo_start

  LOGICAL FUNCTION forcetheo_next_job(this,lastiter,noco)
    USE m_types_setup
    IMPLICIT NONE
    CLASS(t_forcetheo),INTENT(INOUT):: this
     LOGICAL,INTENT(IN)                  :: lastiter
    !Stuff that might be modified...
    TYPE(t_noco),INTENT(INOUT) :: noco
    forcetheo_next_job=this%firstloop
    this%firstloop=.FALSE.
  END FUNCTION forcetheo_next_job

  FUNCTION forcetheo_eval(this,results)RESULT(skip)
    USE m_types_misc
    IMPLICIT NONE
    CLASS(t_forcetheo),INTENT(INOUT):: this
    LOGICAL :: skip
    !Stuff that might be used...
    TYPE(t_results),INTENT(IN) :: results
    skip=.FALSE.
  END FUNCTION forcetheo_eval

  SUBROUTINE forcetheo_postprocess(this)
    IMPLICIT NONE
    CLASS(t_forcetheo),INTENT(INOUT):: this
  END SUBROUTINE forcetheo_postprocess

  
END MODULE m_types_forcetheo


!!$MODULE t_extends
!!$  USE m_types_forcetheo
!!$  TYPE,EXTENDS(t_forcetheo):: t_jij
!!$   CONTAINS
!!$     PROCEDURE :: next_job=>jij_next_job
!!$  END TYPE t_jij
!!$
!!$  TYPE,EXTENDS(t_forcetheo):: t_dmi
!!$   CONTAINS
!!$     PROCEDURE :: next_job=>dmi_next_job
!!$  END TYPE t_dmi
!!$CONTAINS
!!$  LOGICAL FUNCTION jij_next_job(this,mpi)
!!$    IMPLICIT NONE
!!$    CLASS(t_jij),INTENT(INOUT):: this
!!$    INTEGER :: mpi
!!$  END FUNCTION jij_next_job
!!$
!!$  LOGICAL FUNCTION dmi_next_job(this,mpi)
!!$    IMPLICIT NONE
!!$    CLASS(t_dmi),INTENT(INOUT):: this
!!$    INTEGER :: mpi
!!$  END FUNCTION dmi_next_job
!!$END MODULE t_extends
