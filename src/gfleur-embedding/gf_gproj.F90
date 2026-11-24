!--------------------------------------------------------------------------------
! Copyright (c) 2026 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions 
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------
      MODULE m_gf_proj 
      IMPLICIT NONE
!********************************************************************** 
!     This MODULE CONTAINS two subroutines to calculate the Projections 
!     of the Green's function onto the boundary surfaces                
!     gproj: Projection of the Green fct                                
!     dgproj: Projection of the dreivative of the Green fct             
!********************************************************************** 
      PRIVATE 
      PUBLIC::gf_gproj,gf_dgproj,gf_gprojnohelpregion 
      CONTAINS 
      !<--S:gf_GProj(r1,r2,Nv2,Nv,cell,k,l_noco,G,Gij)                  
      SUBROUTINE gf_GProj(                                              &
     &                 r1,r2,jspin,lapw,lapw_gf,cell,l_noco,l_sph,G,                  &
     &                 Gij)                                             
!********************************************************************** 
!     * This subroutine calculates projections of the Greenfunktion     
!     * gij = g(z1,z2)                                                  
!     * arguments:                                                      
!     * r1,r2  z-positions of planes r=0 for left side, r=1 for right si
!     *                                                                 
!     *                           Daniel Wortmann, Tokyo, 2001,         
!                                         cleaned       Juelich 2005    
!                                                                       
                                                                        
!********************************************************************** 
      USE m_gf_types 
                                                                        
      IMPLICIT NONE 
      !<--Arguments                                                     
      INTEGER, INTENT(IN)       :: r1, r2,jspin 
      TYPE(t_cell),INTENT(IN)   :: cell 
      TYPE(t_lapw),INTENT(IN)   :: lapw 
      TYPE(t_lapw_gf),INTENT(IN) :: lapw_gf
      LOGICAL, INTENT(IN)       :: l_noco,l_sph
      COMPLEX, INTENT(IN)       :: G(:,:) 
      COMPLEX, INTENT(OUT)      :: Gij(:,:) 
      !>                                                                
      !<--locals                                                        
      COMPLEX :: exp1, exp2, f1 
      INTEGER :: n1, n2 ,nv_tot
                                                                        
      !>                                                                
                                                                        
!!!   The surface projection of $G(\vec g,\vec g')$ given by:           
!!!   \[                                                                
!!!   G_{z,z'}(\vec g_parallel,\vec g'_parallel)=                       
!!!   /frac{1}{\tilde d}\sum_{g_z,g'_z)e^{i g_z z} e^{-i g'_z z'} G(\vec
!!!   g,\vec g')                                                        
!!!   \]                                                                
                                                                        
         Gij = CMPLX(0.0,0.0) 
                                                                        
         !factor is r1*2-1                                              
         nv_tot=lapw%nv_tot
         if (l_sph) nv_tot=lapw_gf%nv_tot_sphere
         exp1 = EXP(CMPLX(0.0,(r1*2-1)*cell%bmat(3,3)*cell%Z1)) 
         exp2 = EXP(CMPLX(0.0,-1.*(r2*2-1)*cell%bmat(3,3)*cell%Z1)) 
         IF (l_noco) THEN 
           DO n1 = 1, nv_tot/2
            f1 = 1/cell%amat(3,3)*exp1**lapw%k3(n1,jspin) 
            DO n2 = 1, nv_tot/2
               Gij(lapw%kp(n1,jspin),lapw%kp(n2,jspin))=            &
     &              Gij(lapw%kp(n1,jspin),lapw%kp(n2,jspin)) + f1   &
     &              *exp2**lapw%k3(n2,jspin)*G(n1,n2)                 
               Gij(lapw%kp(n1,jspin)+lapw_gf%nv2_tot/2,lapw%kp(n2,jspin&
     &              )) = Gij(lapw%kp(n1,jspin)+lapw_gf%nv2_tot/2         &
     &              ,lapw%kp(n2,jspin))+f1*exp2**lapw%k3(n2,jspin)  &
     &              *G(n1+nv_tot/2,n2)
               Gij(lapw%kp(n1,jspin),lapw%kp(n2,jspin)+lapw_gf%nv2_tot &
     &              /2) = Gij(lapw%kp(n1,jspin),lapw%kp(n2,jspin)   &
     &              +lapw_gf%nv2_tot/2)+f1*exp2**lapw%k3(n2,jspin)  *G(n1&
     &              ,n2+nv_tot/2)
               Gij(lapw%kp(n1,jspin)+lapw_gf%nv2_tot/2,lapw%kp(n2,jspin&
     &              )+lapw_gf%nv2_tot/2) = Gij(lapw%kp(n1,jspin)         &
     &              +lapw_gf%nv2_tot/2,lapw%kp(n2,jspin)+lapw_gf%nv2_tot/2) &
     &              +f1*exp2**lapw%k3(n2,jspin) *G(n1+nv_tot/2,n2&
     &              +nv_tot/2)
            ENDDO 
           ENDDO 
         ELSE 
           DO n1 = 1, nv_tot
            f1 = 1/cell%amat(3,3)*exp1**lapw%k3(n1,jspin) 
            DO n2 = 1, nv_tot
               Gij(lapw%kp(n1,jspin),lapw%kp(n2,jspin)) =           &
     &              Gij(lapw%kp(n1,jspin),lapw%kp(n2,jspin)) + f1   &
     &              *exp2**lapw%k3(n2,jspin)*G(n1,n2)                 
            ENDDO 
           ENDDO 
         ENDIF 
                                                                        
                                                                        
                                                                        
                                                                        
      END SUBROUTINE gf_gproj 
                                                                        
      SUBROUTINE gf_dgproj(G,Gij,Z1,Z2,Nmat,Nv2,Nv,Bz,Dt,diag,fact,     &
     &     Kp,K3)                                                       
!********************************************************************** 
!     * This subroutine calculates projections of the Derivative of the 
!     Greenfunktion  dgij=dg(z1,z2)                                     
!     * arguments:                                                      
!     * z1,z2  z-positions of planes                                    
!     * bz     reciprocal lattice vector in z-direction                 
!     * kp(nv) 2D-index of basis vector as generated by modified apws   
!     * k3(nv) z-index of basis vector                                  
!     *                                                                 
!     *                                                                 
!     *                           Daniel Wortmann, Tokyo, 2001          
!********************************************************************** 
      IMPLICIT NONE 
                                                                        
      INTEGER, INTENT(IN)  :: Nmat, Nv2, Nv 
      REAL, INTENT(IN)     :: Z1, Z2, Bz, Dt,diag,fact 
      INTEGER, INTENT(IN)  :: Kp(Nv), K3(Nv) 
      COMPLEX, INTENT(IN)  :: G(Nmat,Nmat) 
      COMPLEX, INTENT(OUT) :: Gij(Nv2,Nv2) 
                                                                        
!     locals                                                            
                                                                        
      COMPLEX :: exp1, exp2, f1 
      INTEGER :: n1, n2 
                                                                        
!     use f90-like initialisation of Gij                                
      Gij = CMPLX(0.0,0.0) 
!                                                                       
!     On the diagonal substract 2 because of delta-fcn                  
!                                                                       
      DO n1=1,nv2 
         Gij(n1,n1)=diag 
      ENDDO 
                                                                        
      exp1 = EXP(CMPLX(0.0,-Bz*Z1)) 
      exp2 = EXP(CMPLX(0.0,Bz*Z2)) 
      DO n1 = 1, Nv 
         f1 = 1/Dt*exp1**K3(n1) 
         DO n2 = 1, Nv 
            Gij(Kp(n1),Kp(n2)) = Gij(Kp(n1),Kp(n2))+fact*               &
     &         cmplx(0,-1.*BZ*k3(n2))*f1*exp2**K3(n2)*G(n1,n2)          
         ENDDO 
      ENDDO 
      END SUBROUTINE gf_dgproj 
                                                                        
!*******************************************************************    
!     l_nohelpregion=.true. => Project onto curvy boundary planes to    
!     circumvent the use of help regions.                               
!     Frank Freimuth, October 2007                                      
!*******************************************************************    
      SUBROUTINE gf_gprojnohelpregion(layer,cell,lapw,lapw_gf,l_noco,g,gij)
      USE m_gf_types 
      USE m_gf_curvy2dprojector 
      IMPLICIT NONE 
      INTEGER,INTENT(IN)        :: layer
      TYPE(t_cell),INTENT(IN)   :: cell 
      TYPE(t_lapw),INTENT(IN)   :: lapw 
      TYPE(t_lapw_gf),INTENT(IN) :: lapw_gf
      LOGICAL, INTENT(IN)       :: l_noco 
      COMPLEX, INTENT(IN)       :: G(:,:) 
      COMPLEX, INTENT(OUT)      :: Gij(:,:) 
                                                                        
      COMPLEX,ALLOCATABLE :: projtwodright(:,:) 
      COMPLEX,ALLOCATABLE :: gnoughtp(:,:) 
      INTEGER matsize 
                                                                        
         matsize=size(g,1) 
                                                                        
                                                                        
         CALL gf_curvy2dprojector(layer,cell,lapw,lapw_gf,.TRUE.) 
                                                                        
                                                                        
         ALLOCATE( projtwodright(2*lapw_gf%nv2_tot,matsize) ) 
         ALLOCATE(gnoughtp(matsize,2*lapw_gf%nv2_tot)) 
         CALL gf_get_curvy2dproject(lapw,lapw_gf,projtwodright) 
!         projtwodright(1:lapw_gf%nv2_tot,1:matsize)=                      
!     =        curvyproject(1:lapw_gf%nv2_tot,1:matsize,1)                 
!         projtwodright(lapw_gf%nv2_tot+1:2*lapw_gf%nv2_tot,1:matsize)=       
!     =        curvyproject(1:lapw_gf%nv2_tot,1:matsize,2)                 
         CALL zgemm('N','C',matsize,2*lapw_gf%nv2_tot,matsize,    &
     &                           cmplx(1.0,0.0),                        &
     &                           g,matsize,projtwodright,2*lapw_gf%nv2_tot,&
     &                           cmplx(0.0,0.0),gnoughtp,matsize)       
         CALL zgemm('N','N',2*lapw_gf%nv2_tot,2*lapw_gf%nv2_tot,     &
     &                   matsize,cmplx(1.0,0.0),projtwodright,          &
     &                   2*lapw_gf%nv2_tot,gnoughtp,matsize,cmplx(0.0,0.0),&
     &                   gij,2*lapw_gf%nv2_tot)                            
!         gij=matmul(projtwodright,gnoughtp)                            
         DEALLOCATE(projtwodright) 
         DEALLOCATE(gnoughtp) 
                                                                        
      END SUBROUTINE gf_gprojnohelpregion 
                                                                        
                                                                        
      END                                           
