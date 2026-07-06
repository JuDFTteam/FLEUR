!--------------------------------------------------------------------------------
! Copyright (c) 2026 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions 
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------
      MODULE m_gf_CBS_embtest 
          IMPLICIT NONE
!*****************************************************************      
! DESC:Try to improve embedding potential                               
!                          Daniel Wortmann, Wed Nov 19 04:40:09 2003    
!*****************************************************************      
      CONTAINS 
                                                                        
      !<-- S:gf_CBS_restricted(emb,ev,ew,curr)                          
      SUBROUTINE gf_CBS_restricted(emb,phi,dphi,ew,maxkappa) 
!-----------------------------------------------                        
!                                                                       
!           (last modified:09-06-26) D. Wortmann                        
!-----------------------------------------------                        
      USE m_gf_math 
      IMPLICIT NONE 
      !<-- Arguments                                                    
      COMPLEX,INTENT(OUT)    :: emb(:,:) 
      COMPLEX,INTENT(IN)     :: phi(:,:),dphi(:,:) 
      COMPLEX,INTENT(IN)     :: ew(:) 
      REAL   ,INTENT(IN)     :: maxkappa 
      !>                                                                
      !<-- Locals                                                       
      COMPLEX,ALLOCATABLE :: overlapp(:,:) 
      COMPLEX,ALLOCATABLE :: phi_m(:,:) 
      COMPLEX,ALLOCATABLE :: phi_p(:,:) 
      COMPLEX,ALLOCATABLE :: dphi_m(:,:) 
      INTEGER             :: n,nn,n2,m,mm 
      !>                                                                
      !<--select states for embedding potential                         
      n2 = SIZE(ew) 
      m  = COUNT(ABS(AIMAG(ew))<maxkappa) 
      WRITE(*,*) "Selected ",m," out of ",n2," for Sigma" 
      ALLOCATE(overlapp(m,m)) 
      ALLOCATE(phi_m(n2,m)) 
      ALLOCATE(phi_p(n2,m)) 
      ALLOCATE(dphi_m(n2,m)) 
      mm = 0 
      DO n = 1,SIZE(ew) 
         IF (ABS(AIMAG(ew(n)))<maxkappa) THEN 
            mm = mm+1 
            phi_m(:,mm) = phi(:,n) 
            dphi_m(:,mm) = dphi(:,n) 
         ENDIF 
      ENDDO 
      !>                                                                
      !<-- construct inverse overlapp matrix                            
      DO n = 1,m 
         DO nn = 1,m 
            overlapp(n,nn) = dot_PRODUCT(phi_m(:,n),phi_m(:,nn)) 
         ENDDO 
      ENDDO 
      overlapp = mat_inverse(overlapp) 
      !>                                                                
      !<-- construct projectors                                         
  !    DO n = 1,m                                                       
  !       phi_p(n,:) = priv_ortho(phi_m,n)                              
  !       PRINT*,dot_PRODUCT(phi_p(n,:),phi_m(:,m))                     
  !    ENDDO                                                            
      phi_p=MATMUL(phi_m,overlapp) 
      !>                                                                
      !<-- construct the embeding potential                             
      emb = 0 
      DO n = 1,m 
         emb = emb+direct_prod(dphi_m(:,n),phi_p(:,n)) 
      ENDDO 
      !>                                                                
                                                                        
      END SUBROUTINE 
      !>                                                                
                                                                        
      !<-- S:gf_CBS_embtest(emb,ev,ew,curr)                             
      SUBROUTINE gf_CBS_embtest(emb1,ev,ew,curr) 
!******************************************                             
!                                                                       
!                          D. Wortmann                                  
!******************************************                             
      IMPLICIT NONE 
      !<--Arguments                                                     
      COMPLEX,INTENT(INOUT) :: emb1(:,:) 
      COMPLEX,INTENT(IN)    ::ev(:,:),ew(:) 
      REAL,INTENT(IN)       ::curr(:) 
      !>                                                                
      !<--Locals                                                        
      INTEGER :: n1,n2,n 
      INTEGER :: nv2 
      LOGICAL :: ok 
      REAL    :: dcurr,dcurr2 
                                                !tmp storage of emb     
      COMPLEX :: emb(SIZE(emb1,1),SIZE(emb1,2)) 
                                                !projection vectors...  
      COMPLEX :: b(SIZE(emb1,1),SIZE(emb1,2)) 
      !>                                                                
      nv2 = size(curr) 
      emb=emb1 
      !<--Calculate Projectors for each eigenvector                     
      DO n=1,nv2 
         b(:,n)=priv_ortho(ev,n) 
      ENDDO 
      !>                                                                
                                                                        
      !<--Correct Diagonal elments!                                     
      ok=.FALSE. 
      DO WHILE(.NOT.OK) 
         ok=.TRUE. 
         DO n1=1,nv2 
            dcurr=2*AIMAG(dot_PRODUCT(ev(:,n1),MATMUL(emb,ev(:,n1)))) 
            IF (AIMAG(ew(n1))<1E-6) THEN 
               !Bloch state!                                            
               dcurr=dcurr-curr(n1) 
            ELSE 
               !Ev. state                                               
               dcurr=dcurr 
            ENDIF 
            IF (abs(dcurr)>1E-9) THEN 
               !Correct the embpot                                      
               ok=.FALSE. 
               emb=emb-CMPLX(0.0,dcurr/2.0)*direct_prod(b(:,n1),b(:,n1)) 
               dcurr2=2*AIMAG(dot_PRODUCT(ev(:,n1),MATMUL(emb,ev(:,n1)))&
     &              )                                                   
               IF(AIMAG(ew(n1))<1E-6) THEN 
                  !Bloch state!                                         
                  dcurr2=dcurr2-curr(n1) 
               ELSE 
                  !Ev. state                                            
                  dcurr2=dcurr2 
               ENDIF 
               IF (ABS(dcurr2)>ABS(dcurr)) STOP                         &
     &              'Error1 in gf_CBS_embtest'                          
            ENDIF 
         ENDDO 
      ENDDO 
      !>                                                                
      !<--Test off-diagonal elements                                    
      dcurr2=0.0 
      ok=.FALSE. 
      DO WHILE(.NOT.OK) 
         ok=.TRUE. 
         DO n1=1,nv2 
            DO n2=1,nv2 
               IF (n1==n2) CYCLE 
               dcurr=0.5*cmplx(0,1.)*(dot_PRODUCT(ev(:,n1),MATMUL(emb   &
     &              ,ev(:,n2)))-dot_PRODUCT(MATMUL(emb,ev(:,n1)),ev(:,2)&
     &              ))                                                  
               dcurr2=max(abs(dcurr),dcurr2) 
!               IF (abs(dcurr)>1E-9) THEN                               
!                  !Correct the embpot                                  
!                  ok=.FALSE.                                           
!                  emb=emb-CMPLX(0,dcurr)*direct_prod(b(:,n1),b(:,n2))  
!                  dcurr2=0.5*CMPLX(0,1.)*(dot_PRODUCT(ev(:,n1)         
!     +                 ,MATMUL(emb,ev(:,n2)))-dot_PRODUCT(MATMUL(emb   
!     +                 ,ev(:,n1)),ev(:,2)))                            
!                  write(*,*) 'a',dcurr,dcurr2                          
!                  IF (ABS(dcurr2)>ABS(dcurr)) STOP                     
!     +              'Error2 in gf_CBS_embtest'                         
!               ENDIF                                                   
            ENDDO 
         ENDDO 
      ENDDO 
      !>                                                                
      WRITE(*,*) 'Corrected Embedding potential:',dcurr2 
      WRITE(*,*) MAXVAL(ABS(emb-emb1)),SUM(ABS(emb-emb1))/SIZE(emb1) 
!      emb1=emb                                                         
      END SUBROUTINE 
      !>                                                                
                                                                        
      !<-- F:direct_prod(a,b)RESULT(c)                                  
      FUNCTION direct_prod(a,b)RESULT(c) 
!******************************************                             
!     Calculated direct product of two complex :: vectors               
!     D. Wortmann                                                       
!******************************************                             
      IMPLICIT NONE 
      !<--Arguments                                                     
      COMPLEX,INTENT(IN) :: a(:),b(:) 
      COMPLEX            :: c(SIZE(a),size(b)) 
      !>                                                                
                                                                        
      c=SPREAD(a,dim=2,ncopies=SIZE(b))*SPREAD(conjg(b),dim=1,ncopies   &
     &     =SIZE(a))                                                    
      END FUNCTION 
      !>                                                                
                                                                        
      !<-- F:priv_ortho(ev,n)Result(v)                                  
      FUNCTION priv_ortho(ev,n)RESULT(v) 
!******************************************                             
!     Return vector orthonormal to ev(:,n)                              
!     Uses QR-Factorization by lapack                                   
!                          D. Wortmann                                  
!******************************************                             
#include "cpp_double.h"
      IMPLICIT NONE 
      !<--Arguments                                                     
      COMPLEX,INTENT(IN) :: ev(:,:) 
      INTEGER,INTENT(IN) :: n 
      COMPLEX            :: v(SIZE(ev,1)) 
      !>                                                                
      !<--Locals                                                        
      COMPLEX :: mat(SIZE(ev,1),SIZE(ev,2)) 
      COMPLEX :: vec(Size(ev,1)),work(Size(ev,1)) 
      INTEGER :: i,ii,nv2,info 
      !>                                                                
      nv2=SIZE(ev,1) 
      !<--Generate Matrix such that last vector is vector n             
      ii=0 
      DO i=1,nv2 
         IF (i==n) CYCLE 
         ii=ii+1 
         mat(:,ii)=ev(:,i) 
      ENDDO 
      mat(:,nv2)=ev(:,n) 
      !>                                                                
      !<--Call LAPACK for QR-factorization                              
      CALL CPP_LAPACK_cgeqrf( nv2, nv2, mat, nv2, vec, WORK, nv2, INFO ) 
      CALL CPP_LAPACK_cungqr( nv2,nv2,nv2,mat,nv2,vec,WORK,nv2,INFO ) 
      !>                                                                
      !<--Normalize vector such that <O|ev>=1                           
      v=mat(:,nv2)/dot_product(mat(:,nv2),ev(:,n)) 
      !>                                                                
                                                                        
      END FUNCTION 
      !>                                                                
      END                                           
