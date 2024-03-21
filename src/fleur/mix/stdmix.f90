!--------------------------------------------------------------------------------
! Copyright (c) 2016 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------

MODULE m_stmix
  !! Simple or straight mixing
  
CONTAINS
  SUBROUTINE stmix(atoms,input,noco,fsm,fsm_mag,sm)
    !!Simple mixing
    USE m_types_mixvector
    USE m_types
    IMPLICIT NONE
    TYPE(t_input),INTENT(IN)        :: input  
    TYPE(t_noco),INTENT(IN)         :: noco
    TYPE(t_atoms),INTENT(IN)        :: atoms
    TYPE(t_mixvector),INTENT(IN)    :: fsm !!Difference of input and output density
    TYPE(t_mixvector),INTENT(IN)    :: fsm_mag !!Difference of input and output magnetisation
    TYPE(t_mixvector),INTENT(INOUT) :: sm          !!input: input-density, output: mixed density 
    !     ..
    !     .. Local Scalars ..
    REAL,PARAMETER:: tol_6=1.0e-6
    !     ..
    !
    sm = sm + input%alpha*fsm

    IF ( ABS(input%spinf-1.0).LE.tol_6 .OR. input%jspins.EQ.1 .or.input%imix.ne.0.or.noco%l_noco) RETURN
    !  Done with     sm1 = sm + alpha * F(sm)
    
    
    !Spin enhancement factor
    sm = sm + input%alpha/2.0*(input%spinf-1.0)*fsm_mag
       !     -->perform simple mixing with the mixing parameters
       !        for charge and spin
       !
       !       sm1+/_ = (sm+/_) + alpha* F(sm)
       !                +/-0.5alpha(spinf-1)( F(sm+) + F(sm-) )
       ! The F(sm+) and F(sm-) terms do not only include diagonal elements of the density matrices (as one could think) 
       ! but also contain off-diag. elements (jspins=3,4) of the density matrices in the fully noncolinear case. 
       ! Choosing a spinf>1 therefore might be helpful when it comes to converging noncolinear systems.
       ! DW: Actually for noco, this spinf>1 has the potential to break symmetry in the noco case, hence it is disabled. 
   
  END SUBROUTINE stmix
END MODULE m_stmix
