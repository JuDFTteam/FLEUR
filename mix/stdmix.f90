!--------------------------------------------------------------------------------
! Copyright (c) 2016 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------

MODULE m_stmix
  !
  !      straight mixing,     r.pentcheva, iff, 1996
  !
  !     sm   : input charge density of iteration m
  !     sm1  : input charge density of iteration m+1
  !     fsm  : output minus input charge densityof iteration m
  !
CONTAINS
  SUBROUTINE stmix(&
       &                 atoms,input,noco,&
       &                 fsm,fsm_mag,sm)
    USE m_types_mixvector
    USE m_types
    IMPLICIT NONE
    TYPE(t_input),INTENT(IN)        :: input
    TYPE(t_noco),INTENT(IN)         :: noco
    TYPE(t_atoms),INTENT(IN)        :: atoms
    TYPE(t_mixvector),INTENT(IN)    :: fsm,fsm_mag
    TYPE(t_mixvector),INTENT(INOUT) :: sm
    !     ..
    !     .. Local Scalars ..
    INTEGER imap
    REAL,PARAMETER:: tol_6=1.0e-6
    !     ..
    !
    sm = sm + input%alpha*fsm

    IF ( ABS(input%spinf-1.0e0).LE.tol_6 .OR. input%jspins.EQ.1 .or.input%imix.ne.0) THEN
       !  Done with     sm1 = sm + alpha * F(sm)
       !No spin
       RETURN
    ELSE
       sm = sm + input%alpha/2.0*(input%spinf-1.0)*fsm_mag
       !     -->perform simple mixing with the mixing parameters
       !        for charge and spin
       !
       !       sm1+/_ = (sm+/_) + alpha* F(sm)
       !                +/-0.5alpha(spinf-1)( F(sm+) + F(sm-) )
    END IF

  END SUBROUTINE stmix
END MODULE m_stmix
