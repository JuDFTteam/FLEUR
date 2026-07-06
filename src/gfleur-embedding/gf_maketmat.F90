!--------------------------------------------------------------------------------
! Copyright (c) 2026 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions 
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------
      MODULE m_gf_maketmat 
      use m_juDFT
          IMPLICIT NONE
      CONTAINS 
      SUBROUTINE gf_maketmat_neumann(                                   &
     &                                Nv2,Gij,                          &
     &                                tr)                               
!********************************************************************** 
!     This SUBROUTINE takes the projections of the Gmatrix as calculated
!     gf_gfcn and calculates the Transfer Matrix from that              
!                                                                       
!    WARNING, THIS SUBROUTINE USES $\Psi = (\psi,-\partial\psi)$        
!                                Daniel Wortmann, Tokyo, 2001           
!     (last modified: 04-09-07)                                         
!********************************************************************** 
!     For better performance: calls to blas library instead of matmul.  
!     Frank Freimuth, November 2007                                     
!********************************************************************** 
      USE m_gf_types,  ONLY:t_gfinp 
      USE m_gf_math,  ONLY:mat_inverse,zmat_product 
#ifdef CPP_CUOVLP                                                       
      USE m_gf_curvy2dprojector,ONLY:basisoverlaps 
#endif                                                                  
      IMPLICIT NONE 
      !<-- Arguments                                                    
      INTEGER,      INTENT(IN) :: nv2 
      COMPLEX,INTENT(IN)       :: Gij(2*nv2,2*nv2) 
      COMPLEX,INTENT(OUT)      :: TR(2*nv2,2*nv2) 
      !>                                                                
                                                                        
      !<--Locals                                                        
      COMPLEX,ALLOCATABLE :: G11(:,:),G12(:,:),G21(:,:),G22(:,:) 
      ALLOCATE(g11(nv2,nv2),g12(nv2,nv2),g21(nv2,nv2),g22(nv2,nv2)) 
      !>                                                                
                                                                        
      g11=gij(1:nv2,1:nv2) 
      g12=gij(1:nv2,nv2+1:2*nv2) 
      g21=gij(nv2+1:2*nv2,1:nv2) 
      g22=gij(nv2+1:2*nv2,nv2+1:2*nv2) 
                                                                        
#ifdef CPP_CUOVLP                                                       
         g11=matmul(g11,basisoverlaps(:,:)) 
         g12=matmul(g12,basisoverlaps(:,:)) 
         g21=matmul(g21,basisoverlaps(:,:)) 
         g22=matmul(g22,basisoverlaps(:,:)) 
#endif                                                                  
                                                                        
      g12 = mat_inverse(g12) 
                                                                        
      tr(1:nv2,1:nv2)=zmat_product(g22,g12) 
                                                                        
      TR(Nv2+1:2*Nv2,1:Nv2) =cmplx(- 2.0,0.0)*g12 
                                                                        
      TR(Nv2+1:2*Nv2,Nv2+1:2*Nv2)=zmat_product(g12,g11) 
                                                                        
      tr(1:nv2,nv2+1:2*nv2)=                                            &
     &          zmat_product(g22,tr(nv2+1:2*nv2,nv2+1:2*nv2))           
                                                                        
      TR(1:Nv2,Nv2+1:2*Nv2)=cmplx(0.5,0.0)*(g21-tr(1:nv2,nv2+1:2*nv2)) 
                                                                        
      DEALLOCATE(g11,g12,g21,g22) 
      END SUBROUTINE gf_maketmat_neumann 
                                                                        
      SUBROUTINE gf_maketmat2(layer,l,Nv2,rkp,TR,gfinp,en,pot_aux) 
!*****************************************************************      
! This subroutine calculates the inverse T-matrix for region II         
!*****************************************************************      
      USE m_gf_types,  ONLY:t_gfinp 
      USE m_gf_energies  ,ONLY:gf_z 
                                                                        
      IMPLICIT NONE 
      INTEGER,INTENT(IN)        :: layer 
      REAL,INTENT(IN)           :: l 
      INTEGER,      INTENT(IN)  :: nv2,en 
      TYPE(t_gfinp),INTENT(IN)  :: gfinp 
      COMPLEX,INTENT(INOUT)     :: TR(2*nv2,2*nv2) 
      REAL,          INTENT(IN) :: rkp(nv2) 
      COMPLEX,INTENT(IN)        :: pot_aux 
                                                                        
! locals                                                                
      INTEGER n 
      COMPLEX :: k 
                                                                        
      TR=0.0 
      IF (l>epsilon(0.0)) THEN 
                                                                        
         DO n=1,nv2 
            k = SQRT(2*(gf_z(en,layer)-pot_aux)-rkp(n)**2) 
            IF (ABS(k)<3*EPSILON(0.0)) THEN 
               k=CMPLX(0.0,3*EPSILON(0.0)) 
            ENDIF 
            TR(n,n) = cos(k*l) 
            TR(n,n+nv2) = -sin(k*l)/k 
            TR(n+nv2,n) = k*sin(k*l) 
            TR(n+nv2,n+nv2) = cos(k*l) 
         ENDDO 
      ELSE 
! zero sized region II, TR as identity matrix                           
         DO n=1,2*nv2 
            TR(n,n)=1.0 
         ENDDO 
      ENDIF 
                                                                        
      END SUBROUTINE gf_maketmat2 
                                                                        
      !<-- S:gf_maketmat(Nv2,Gij,DGij,TR,gfinp)                         
      SUBROUTINE gf_maketmat(Nv2,Gij,DGij,                              &
     &     TR,gfinp)                                                    
!********************************************************************** 
!     This SUBROUTINE takes the projections of the Gmatrix as calculated
!     gf_gfcn and calculates the Transfer Matrix from that              
!                                                                       
!    WARNING, THIS SUBROUTINE USES $\Psi = (\psi,-\partial\psi)$        
!                                Daniel Wortmann, Tokyo, 2001           
!     (last modified: 04-09-07)                                         
!********************************************************************** 
      USE m_gf_types,  ONLY:t_gfinp 
      USE m_gf_math,  ONLY:mat_inverse 
      IMPLICIT NONE 
      !<-- Arguments                                                    
      INTEGER,      INTENT(IN) :: nv2 
      TYPE(t_gfinp),INTENT(IN) :: gfinp 
      COMPLEX,INTENT(INOUT)    :: Gij(2*nv2,2*nv2,2),DGij(2*nv2,2*nv2,2) 
      COMPLEX,INTENT(INOUT)    :: TR(2*nv2,2*nv2,2) 
      !>                                                                
                                                                        
      !<--Locals                                                        
      COMPLEX,ALLOCATABLE :: G11(:,:),G12(:,:),G21(:,:),G22(:,:) 
      INTEGER             :: region 
                                                                        
      ALLOCATE(g11(nv2,nv2),g12(nv2,nv2),g21(nv2,nv2),g22(nv2,nv2)) 
      !>                                                                
                                                                        
                    !loop over both regions                             
      DO region=1,1 
                                                                        
!                                                                       
!     Either USE original matrices or swap 1<->2                        
!                                                                       
         IF (region==1) THEN 
            g11=gij(1:nv2,1:nv2,1) 
            g12=gij(1:nv2,nv2+1:2*nv2,1) 
            g21=gij(nv2+1:2*nv2,1:nv2,1) 
            g22=gij(nv2+1:2*nv2,nv2+1:2*nv2,1) 
         ELSE 
            gij(:,:,2)=-1.0*gij(:,:,2) 
            !if (gfinp%l_addemb) dgij(:,:,2)=-1.0*dgij(:,:,2)           
            g22=gij(1:nv2,1:nv2,2) 
            g21=gij(1:nv2,nv2+1:2*nv2,2) 
            g12=gij(nv2+1:2*nv2,1:nv2,2) 
            g11=gij(nv2+1:2*nv2,nv2+1:2*nv2,2) 
         ENDIF 
         !only the inverse of g12 is needed                             
         g12 = mat_inverse(g12) 
         !>                                                             
                                                                        
         !Now construct T!                                              
                                                                        
         IF (gfinp%l_addemb) THEN 
            !<-- Calculation of T with non-zero derivative of G         
                                                                        
             CALL juDFT_error("Calculating a T-matrix with addemb = T not supported",calledby="gf_maketmat.F90")
                                ! this is a bit more complicated!       
                                ! the normal derivative is not zero     
!            ALLOCATE ( A(Nv2,Nv2) )                                    
!            A=0.0                                                      
!            DO I=1,nv2                                                 
!               A(i,i)=1.0                                              
!            ENDDO                                                      
!            A=A+0.5*dg22-0.5*MATMUL(MATMUL(g22,g12),dg12)              
!            A=cmat_inverse(A)                                          
!            TR(1:Nv2,1:Nv2,region)=MATMUL(a,(MATMUL(g22,g12)-          
!     +           0.5*MATMUL(G22,MATMUL(g12,dg11))+0.5*dg21))           
!            TR(1:Nv2,Nv2+1:2*Nv2,region)=MATMUL(a,                     
!     +           (0.5*MATMUL(MATMUL(g22,g12),g11)-0.5*g21))            
!            TR(Nv2+1:2*Nv2,1:Nv2,region)=2.*g12-MATMUL(g12,dg11)+      
!     +           MATMUL(g12,MATMUL(dg12,TR(1:Nv2,1:Nv2,region)))       
!            TR(Nv2+1:2*Nv2,Nv2+1:2*Nv2,region)=                        
!     +           MATMUL(g12,MATMUL(dg12,TR(1:Nv2,Nv2+1:2*Nv2,region))) 
!     +           +MATMUL(g12,g11)                                      
!             CALL juDFT_error("CHECK sign+L<->R",calledby="gf_maketmat.F90")
         ELSE 
            TR(1:Nv2,1:Nv2,region) = MATMUL(g22,g12) 
            TR(Nv2+1:2*Nv2,1:Nv2,region) =- 2*g12 
            TR(1:Nv2,Nv2+1:2*Nv2,region)                                &
     &           =0.5*(g21-MATMUL(g22,MATMUL(g12,g11)))                 
            TR(Nv2+1:2*Nv2,Nv2+1:2*Nv2,region) =MATMUL(g12,g11) 
            !>                                                          
         ENDIF 
                                                                        
      ENDDO 
      !>                                                                
                                                                        
      DEALLOCATE(g11,g12,g21,g22) 
      END SUBROUTINE gf_maketmat 
                                                                        
      END                                           
