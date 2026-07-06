!--------------------------------------------------------------------------------
! Copyright (c) 2026 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions 
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------
      MODULE m_gf_PhaseMatrix 
      IMPLICIT NONE
      PRIVATE 
      COMPLEX,SAVE,ALLOCATABLE::P(:) 
      PUBLIC getLargePhaseMatrix,initPhaseMatrix,getPhaseMatrix 
      CONTAINS 
      !<--S: initPhaseMatrix(lapw,cell,gfinp,l_noco)                    
      SUBROUTINE initPhaseMatrix(jspin,lapw,cell,gfinp,l_noco) 
      ! INITIALIZES THE MODULE                                          
      ! The Phase factors are calculated and stored in P                
      USE m_gf_types 
      IMPLICIT NONE 
      INTEGER,INTENT(IN)     :: jspin 
      TYPE(t_lapw),INTENT(IN) :: lapw 
      TYPE(t_cell),INTENT(IN) :: cell 
      TYPE(t_gfinp),INTENT(IN) :: gfinp 
      LOGICAL, INTENT(IN)     :: l_noco 
      ! Locals                                                          
      INTEGER::n 
      REAL   ::dp(3),kg(3),dpr(3),kgr(3),s 
                                                                        
      IF (ALLOCATED(P)) THEN 
         IF (SIZE(P)/=lapw%nv2_tot) THEN 
            !The size of the phasematrix must be changed because        
            !2d-basis size is changed...                                
            DEALLOCATE(P) 
            ALLOCATE(P(lapw%nv2_tot)) 
         ENDIF 
      ELSE 
         ALLOCATE(P(lapw%nv2_tot)) 
      ENDIF 
      P(:) = CMPLX(0.0,0.0) 
      dp(3) = 0 
      kg(3) = 0 
      IF (l_noco) THEN 
         DO n = 1, Lapw%nv2_Tot/2 
            dp(1) = gfinp%dp1 
            dp(2) = gfinp%dp2 
            dpr=matmul(cell%amat,dp) 
!            CALL COTRA0(dp,dpr,Cell%amat)                              
                                           !k-dependence cancels!       
            kg(1) =lapw%kp%k1p(n,jspin) 
                                           !                            
            kg(2) =lapw%kp%k2p(n,jspin) 
            kgr=matmul(kg,cell%bmat) 
!           CALL COTRA3(kg,kgr,Cell%bmat)                               
            s = DOT_PRODUCT(kgr,dpr) 
            P(n) = EXP(CMPLX(0,s)) 
            P(n+lapw%nv2_tot/2) = EXP(CMPLX(0,s)) 
         ENDDO 
      ELSE 
         DO n = 1, Lapw%nv2_Tot 
            dp(1) = gfinp%dp1 
            dp(2) = gfinp%dp2 
            dpr=matmul(cell%amat,dp) 
            !CALL COTRA0(dp,dpr,Cell%amat)                              
                                           !k-dependence cancels!       
            kg(1) =lapw%kp%k1p(n,jspin) 
                                           !                            
            kg(2) =lapw%kp%k2p(n,jspin) 
            kgr=matmul(kg,cell%bmat) 
            !CALL COTRA3(kg,kgr,Cell%bmat)                              
            s = DOT_PRODUCT(kgr,dpr) 
            P(n) = EXP(CMPLX(0,s)) 
         ENDDO 
      ENDIF 
      END SUBROUTINE 
      !>                                                                
      FUNCTION getPhaseMatrix() RESULT(PM) 
      !returns a lapw%nv2_tot x lapw%nv2_tot Matrix of Phase factors    
      IMPLICIT NONE 
      COMPLEX::PM(SIZE(P),SIZE(P)) 
      INTEGER::i 
      PM=0.0 
      DO i=1,SIZE(P) 
         PM(i,i)=P(i) 
      ENDDO 
      END FUNCTION 
                                                                        
      FUNCTION getLargePhaseMatrix() RESULT(PM) 
      !returns a 2*lapw%nv2_tot x 2*lapw%nv2_tot Matrix of Phase factors
      IMPLICIT NONE 
      COMPLEX::PM(2*SIZE(P),2*SIZE(P)) 
      INTEGER::nv2,i 
      PM=0.0 
      nv2=SIZE(P) 
      DO i=1,nv2 
         PM(i,i)=P(i) 
         PM(i+nv2,i+nv2)=P(i) 
      ENDDO 
      END FUNCTION 
      END                                           
