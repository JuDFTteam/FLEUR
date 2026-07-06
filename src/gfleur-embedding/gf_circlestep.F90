!--------------------------------------------------------------------------------
! Copyright (c) 2026 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions 
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------
      MODULE m_gf_circlestep 
          IMPLICIT NONE
!****************************************************                   
!     Calculate the 2d-stepfunction that cuts out                       
!     circles in the embedding planes.                                  
!     Frank Freimuth, October 2007                                      
!****************************************************                   
      CONTAINS 
      !<--S: gf_circlestep                                              
                                                                        
      SUBROUTINE gf_circlestep(                                         &
     &                     bmat,                                        &
     &                     taual,rmt,amat33,c,                          &
     &                     mx1,mx2,                                     &
     &                     circlestep)                                  
                                                                        
      USE m_constants, ONLY: pimach 
      USE m_gf_bessel 
                                                                        
      IMPLICIT NONE 
                                                                        
      REAL, INTENT(IN)      :: bmat(3,3) 
      REAL, INTENT(IN)      :: taual(:,:) 
      REAL, INTENT(IN)      :: rmt(:) 
      REAL, INTENT(IN)      :: amat33 
      REAL, INTENT(IN)      :: c 
      INTEGER, INTENT(IN)   :: mx1,mx2 
      COMPLEX,INTENT(INOUT) :: circlestep(-mx1:mx1,-mx2:mx2,2) 
                                                                        
      INTEGER             :: n,k1,k2 
      REAL                :: zcoor,pi,gs 
      REAL                :: gvec(2),rg(2) 
      REAL                :: circleradius2 
      REAL                :: circleradius 
      REAL                :: kpara 
      COMPLEX             :: phase 
      REAL                :: argument 
                                                                        
      pi=pimach() 
      DO n = 1,SIZE(rmt) 
                                                      !cut right dividin
       IF(abs(c/2.0-taual(3,n)*amat33)<rmt(n))THEN 
         zcoor=c/2.0-taual(3,n)*amat33 
         circleradius2 = rmt(n)**2-zcoor**2 
         circleradius  = sqrt(circleradius2) 
         circlestep(0,0,2)=                                             &
     &      circlestep(0,0,2)+pi*circleradius2                          
         DO k2=-mx2,mx2 
            gvec(2)=k2 
            DO k1 = -mx1,mx1 
               gvec(1)=k1 
               IF (k1 == 0.AND.k2 == 0) CYCLE 
               rg  = MATMUL(gvec,bmat(1:2,1:2)) 
               argument=k1*taual(1,n)+k2*taual(2,n) 
               argument=2*pi*argument 
               phase=cmplx(cos(argument),sin(argument)) 
               kpara=sqrt( rg(1)*rg(1) + rg(2)*rg(2) ) 
               circlestep(k1,k2,2) = circlestep(k1,k2,2)+               &
     &              2.0*pi*circleradius/kpara*gf_bessel1(circleradius   &
     &              *kpara)*phase                                       
                  !k1                                                   
            ENDDO 
               !k2                                                      
         ENDDO 
              !right dividing plane is cut                              
       ENDIF 
                                                       !cut left dividin
       IF(abs(-c/2.0-taual(3,n)*amat33)<rmt(n))THEN 
         zcoor=-c/2.0-taual(3,n)*amat33 
         circleradius2 = rmt(n)**2-zcoor**2 
         circleradius  = sqrt(circleradius2) 
         circlestep(0,0,1)=                                             &
     &      circlestep(0,0,1)+pi*circleradius2                          
         DO k2=-mx2,mx2 
            gvec(2)=k2 
            DO k1 = -mx1,mx1 
               gvec(1)=k1 
               IF (k1 == 0.AND.k2 == 0) CYCLE 
               rg  = MATMUL(gvec,bmat(1:2,1:2)) 
               argument=k1*taual(1,n)+k2*taual(2,n) 
               argument=2*pi*argument 
               phase=cmplx(cos(argument),sin(argument)) 
               kpara=sqrt( rg(1)*rg(1) + rg(2)*rg(2) ) 
               circlestep(k1,k2,1) = circlestep(k1,k2,1)+               &
     &              2.0*pi*circleradius/kpara*gf_bessel1(circleradius   &
     &              *kpara)*phase                                       
                  !k1                                                   
            ENDDO 
               !k2                                                      
         ENDDO 
              !left dividing plane is cut                               
       ENDIF 
      ENDDO 
      END SUBROUTINE gf_circlestep 
                                                                        
      !>                                                                
                                                                        
      !<--S: gf_circlestep_embdesc                                      
                                                                        
      SUBROUTINE gf_circlestep_embdesc(                                 &
     &     bmat,                                                        &
     &     embdesc,                                                     &
     &     mx1,mx2,                                                     &
     &     circlestep)                                                  
                                                                        
      USE m_constants, ONLY: pimach 
      USE m_gf_embdesc 
      USE m_gf_bessel 
                                                                        
      IMPLICIT NONE 
                                                                        
      REAL, INTENT(IN)           :: bmat(3,3) 
      TYPE(t_embdesc),INTENT(IN) :: embdesc 
      INTEGER, INTENT(IN)        :: mx1,mx2 
      COMPLEX,INTENT(INOUT)      :: circlestep(-mx1:,-mx2:) 
                                                                        
      INTEGER             :: i,nn,n,k1,k2 
      REAL                :: pi,gs 
      REAL                :: gvec(2),rg(2) 
      REAL                :: circleradius2 
      REAL                :: circleradius 
      REAL                :: kpara 
      COMPLEX             :: phase 
      REAL                :: argument 
                                                                        
      pi=pimach() 
                  !in-out                                               
       DO i = 1,2 
         IF (i == 1) nn = embdesc%cut_atoms_in 
         IF (i == 2) nn = embdesc%cut_atoms_out 
         IF (nn == 0) CYCLE 
         DO n = 1,nn 
            circleradius2 = embdesc%atoms_rmt(n,i)**2                   &
     &           -embdesc%atoms_pos(3,n,i)**2                           
            circleradius = SQRT(circleradius2) 
            circlestep(0,0) =                                           &
     &           circlestep(0,0)+pi*circleradius2                       
            DO k2 =-mx2,mx2 
               DO k1 = -mx1,mx1 
                  IF (k1 == 0.AND.k2 == 0) CYCLE 
                  rg = MATMUL((/1.*k1,1.*k2/),bmat(1:2,1:2)) 
                  argument = rg(1)*embdesc%atoms_pos(1,n,i)+rg(2)       &
     &                 *embdesc%atoms_pos(2,n,i)                        
                  phase = CMPLX(COS(argument),SIN(argument)) 
                  kpara = SQRT( rg(1)*rg(1) + rg(2)*rg(2) ) 
                  circlestep(k1,k2) = circlestep(k1,k2)+                &
     &                 2.0*pi*circleradius/kpara*gf_bessel1(circleradius&
     &                 *kpara)*phase                                    
               ENDDO 
            ENDDO 
         ENDDO 
      ENDDO 
                                                                        
      END SUBROUTINE gf_circlestep_embdesc 
                                                                        
      !>                                                                
                                                                        
      END                                           
