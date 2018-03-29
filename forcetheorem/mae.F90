!--------------------------------------------------------------------------------
! Copyright (c) 2016 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------

MODULE m_types_mae
  USE m_judft
  USE m_types_forcetheo
  TYPE,EXTENDS(t_forcetheo) :: t_forcetheo_mae
     INTEGER :: directions_done
     REAL,ALLOCATABLE:: theta(:)
     REAL,ALLOCATABLE:: phi(:)
     REAL,ALLOCATABLE:: evsum(:)
   CONTAINS
     PROCEDURE :: start   =>mae_start
     PROCEDURE :: next_job=>mae_next_job 
     PROCEDURE :: eval    =>mae_eval
     PROCEDURE :: postprocess => mae_postprocess
     PROCEDURE :: init   => mae_init
  END TYPE t_forcetheo_mae

CONTAINS
  SUBROUTINE mae_init(this,theta_s,phi_s)
    USE m_calculator
    IMPLICIT NONE
    CLASS(t_forcetheo_mae),INTENT(INOUT):: this
    CHARACTER(len=*),INTENT(INOUT)      :: theta_s,phi_s

    CALL evaluateList(this%theta,theta_s)
    CALL evaluateList(this%phi,phi_s)

    IF (SIZE(this%phi).NE.SIZE(this%theta)) CALL &
         judft_error("Lists for theta/phi must have the same length in MAE force theorem calculations")
    ALLOCATE(this%evsum(SIZE(this%phi)))
    this%evsum=0
  END SUBROUTINE mae_init
    

  SUBROUTINE mae_start(this)
    IMPLICIT NONE
    CLASS(t_forcetheo_mae),INTENT(INOUT):: this
    this%directions_done=0
    CALL this%t_forcetheo%start() !call routine of basis type
  END SUBROUTINE  mae_start


  LOGICAL FUNCTION mae_next_job(this,lastiter,noco)
    USE m_types_setup
    USE m_xmlOutput
    IMPLICIT NONE
    CLASS(t_forcetheo_mae),INTENT(INOUT):: this
    LOGICAL,INTENT(IN)                  :: lastiter
    !Stuff that might be modified...
    TYPE(t_noco),INTENT(INOUT) :: noco
       IF (.NOT.lastiter) THEN
          mae_next_job=this%t_forcetheo%next_job(lastiter,noco)
          RETURN
       ENDIF
       !OK, now we start the MAE-loop
       this%directions_done=this%directions_done+1
       mae_next_job=(this%directions_done<=SIZE(this%phi)) !still angles to do

       noco%theta=this%theta(this%directions_done)
       noco%phi=this%phi(this%directions_done)
       IF (this%directions_done.NE.1) CALL closeXMLElement('Forcetheorem_Loop_MAE')
       CALL openXMLElementPoly('Forcetheorem_Loop_MAE',(/'No'/),(/this%directions_done/))
  END FUNCTION mae_next_job

  FUNCTION mae_eval(this,results)RESULT(skip)
    USE m_types_misc
    IMPLICIT NONE
    CLASS(t_forcetheo_mae),INTENT(INOUT):: this
    LOGICAL :: skip
    !Stuff that might be used...
    TYPE(t_results),INTENT(IN) :: results
       
    this%evsum(this%directions_done)=results%seigv
    skip=.TRUE.
  END FUNCTION  mae_eval

  SUBROUTINE mae_postprocess(this)
    USE m_xmlOutput
    IMPLICIT NONE
    CLASS(t_forcetheo_mae),INTENT(INOUT):: this

    !Locals
    INTEGER:: n
    CHARACTER(LEN=12):: attributes(3)
    !Now output the results
    call closeXMLElement('Forcetheorem_Loop_MAE')
    CALL openXMLElementPoly('Forcetheorem_MAE',(/'Angles'/),(/SIZE(this%evsum)/))
    DO n=1,SIZE(this%evsum)
       WRITE(attributes(1),'(f12.7)') this%theta(n)
       WRITE(attributes(2),'(f12.7)') this%phi(n)
       WRITE(attributes(3),'(f12.7)') this%evsum(n)     
       CALL writeXMLElementForm('Angle',(/'theta ','phi   ','ev-sum'/),attributes,&
                                   reshape((/5,3,6,12,12,12/),(/3,2/)))
    END DO
    CALL closeXMLElement('Forcetheorem_MAE')
  END SUBROUTINE mae_postprocess
END MODULE m_types_mae
