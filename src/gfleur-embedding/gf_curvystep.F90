!--------------------------------------------------------------------------------
! Copyright (c) 2026 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions 
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------
      MODULE m_gf_curvystep 
          IMPLICIT NONE
!****************************************************                   
!     Calculate the 2d-stepfunction for the spherical                   
!     hats.                                                             
!     Frank Freimuth, October 2007                                      
!****************************************************                   
      CONTAINS 
      !<-- S: subroutine gf_curvystep                                   
      SUBROUTINE gf_curvystep(                                          &
     &                     bmat,                                        &
     &                     taual,rmt,amat33,c,                          &
     &                     mx1,mx2,mx3,                                 &
     &                     curvystep)                                   
                                                                        
      USE m_constants, ONLY: pimach 
      USE m_gf_curvypartfourier, ONLY: gf_curvypartfou 
                                                                        
      IMPLICIT NONE 
                                                                        
      REAL, INTENT(IN)      :: bmat(3,3) 
      REAL, INTENT(IN)      :: taual(:,:) 
      REAL, INTENT(IN)      :: rmt(:) 
      REAL, INTENT(IN)      :: amat33 
      REAL, INTENT(IN)      :: c 
      INTEGER, INTENT(IN)   :: mx1,mx2,mx3 
      COMPLEX,INTENT(INOUT) :: curvystep(-mx1:mx1,-mx2:mx2,-mx3:mx3,2) 
                                                                        
      INTEGER             :: n,k1,k2,k3 
      REAL                :: zcoor,pi,gs 
      REAL                :: gvec(3),rg(3) 
      REAL                :: circleradius2 
      REAL                :: circleradius 
      REAL                :: kpara,korthogo 
      COMPLEX             :: value 
      COMPLEX             :: phase 
      REAL                :: argument 
                                                                        
      pi=pimach() 
      DO n = 1,SIZE(rmt) 
                                                      !cut right dividin
       IF(abs(c/2.0-taual(3,n)*amat33)<rmt(n))THEN 
         zcoor=c/2.0-taual(3,n)*amat33 
         circleradius2 = rmt(n)**2-zcoor**2 
         circleradius  = sqrt(circleradius2) 
         DO k3=-mx3,mx3 
            gvec(3)=k3 
!            if(k3.eq.0)cycle                                           
            DO k2=-mx2,mx2 
              gvec(2)=k2 
              DO k1 = -mx1,mx1 
               gvec(1)=k1 
               rg  = MATMUL(gvec,bmat) 
               argument=k1*taual(1,n)+k2*taual(2,n) 
               argument=2*pi*argument 
               phase=cmplx(cos(argument),sin(argument)) 
               kpara=sqrt( rg(1)*rg(1) + rg(2)*rg(2) ) 
               korthogo=rg(3) 
               CALL gf_curvypartfou(                                    &
     &              rmt(n),circleradius,kpara,korthogo,                 &
     &              zcoor,                                              &
     &              value)                                              
               IF(zcoor<0.0) value=conjg(value) 
               curvystep(k1,k2,k3,2) = curvystep(k1,k2,k3,2)+           &
     &           2.0*pi*value*phase*exp(cmplx(0.0,-zcoor*korthogo))     
                    !k1                                                 
              ENDDO 
                  !k2                                                   
            ENDDO 
               !k3                                                      
         ENDDO 
              !right dividing plane is cut                              
       ENDIF 
                                                       !cut left dividin
       IF(abs(-c/2.0-taual(3,n)*amat33)<rmt(n))THEN 
         zcoor=-c/2.0-taual(3,n)*amat33 
         circleradius2 = rmt(n)**2-zcoor**2 
         circleradius  = sqrt(circleradius2) 
         DO k3=-mx3,mx3 
            gvec(3)=k3 
!            if(k3.eq.0)cycle                                           
            DO k2=-mx2,mx2 
              gvec(2)=k2 
              DO k1 = -mx1,mx1 
               gvec(1)=k1 
               rg  = MATMUL(gvec,bmat) 
               argument=k1*taual(1,n)+k2*taual(2,n) 
               argument=2*pi*argument 
               phase=cmplx(cos(argument),sin(argument)) 
               kpara=sqrt( rg(1)*rg(1) + rg(2)*rg(2) ) 
               korthogo=rg(3) 
               CALL gf_curvypartfou(                                    &
     &              rmt(n),circleradius,kpara,korthogo,                 &
     &              zcoor,                                              &
     &              value)                                              
               IF(zcoor<0.0) value=conjg(value) 
               curvystep(k1,k2,k3,1) = curvystep(k1,k2,k3,1)+           &
     &           2.0*pi*value*phase*exp(cmplx(0.0,-zcoor*korthogo))     
                    !k1                                                 
              ENDDO 
                  !k2                                                   
            ENDDO 
               !k3                                                      
         ENDDO 
              !left dividing plane is cut                               
       ENDIF 
      ENDDO 
      END SUBROUTINE gf_curvystep 
      !>                                                                
                                                                        
      !<-- S: subroutine gf_curvystep                                   
      SUBROUTINE gf_curvystep_embdesc(                                  &
     &                     bmat,embdesc,right,                          &
     &                     mx1,mx2,mx3,                                 &
     &                     curvystep)                                   
      USE m_constants 
      USE m_gf_embdesc 
      USE m_gf_curvypartfourier, ONLY: gf_curvypartfou 
                                                                        
      IMPLICIT NONE 
                                                                        
      REAL, INTENT(IN)           :: bmat(3,3) 
      TYPE(t_embdesc),INTENT(IN) :: embdesc 
      INTEGER, INTENT(IN)        :: mx1,mx2,mx3 
      LOGICAL,INTENT(IN)         :: right 
      COMPLEX,INTENT(INOUT)      :: curvystep(-mx1:,-mx2:,-mx3: ) 
                                                                        
      INTEGER             :: i,nn,n,k1,k2,k3 
      REAL                :: rg(3) 
      REAL                :: circleradius 
      REAL                :: kpara,korthogo 
      COMPLEX             :: value 
      COMPLEX             :: phase 
      REAL                :: argument,pi,sgn,sgnn 
                                                                        
      pi=pimach() 
      IF (right) THEN 
         sgn =-1.0 
      ELSE 
         sgn = 1.0 
      ENDIF 
                                                                        
                 !in-out                                                
      DO i = 1,2 
         IF (i == 1) THEN 
            sgnn = 1 
            nn = embdesc%cut_atoms_in 
         ELSEIF (i == 2) THEN 
            nn  = embdesc%cut_atoms_out 
            sgnn =-1 
         ENDIF 
         IF (nn == 0) CYCLE 
         DO n = 1,nn 
            circleradius = embdesc%atoms_rmt(n,i)**2-embdesc%atoms_pos(3&
     &           ,n,i)**2                                               
            circleradius = SQRT(circleradius) 
            DO k3 =-mx3,mx3 
               DO k2 =-mx2,mx2 
                  DO k1 = -mx1,mx1 
                     rg = MATMUL((/1.*k1,1.*k2,1.*k3/),bmat) 
                     argument = rg(1)*embdesc%atoms_pos(1,n,i)+rg(2     &
     &                    )*embdesc%atoms_pos(2,n,i)                    
                     phase = CMPLX(COS(argument),SIN(argument)) 
                     kpara = SQRT( rg(1)*rg(1) + rg(2)*rg(2) ) 
                     korthogo = sgnn*sgn*rg(3) 
                     CALL gf_curvypartfou(                              &
     &                    embdesc%atoms_rmt(n,i),circleradius,kpara     &
     &                    ,korthogo,embdesc%atoms_pos(3,n,i),value)     
!                     IF(i == 2) value = CONJG(value)                   
                     curvystep(k1,k2,k3) = curvystep(k1,k2,k3)+         &
     &                    2.0*pi*value*phase*EXP(CMPLX(0.0,             &
     &                    -embdesc%atoms_pos(3,n,i)*korthogo))          
                        !k1                                             
                  ENDDO 
                        !k2                                             
               ENDDO 
                        !k3                                             
            ENDDO 
               !n                                                       
         ENDDO 
               ! in-out                                                 
      ENDDO 
                                                                        
      END SUBROUTINE 
      !>                                                                
                                                                        
      END                                           
