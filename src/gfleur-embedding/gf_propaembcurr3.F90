!--------------------------------------------------------------------------------
! Copyright (c) 2026 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions 
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------
      MODULE m_gf_propaembcurr3 
      IMPLICIT NONE
      PRIVATE 
      PUBLIC::gf_propaembcurr3 
      CONTAINS 
      SUBROUTINE gf_propaembcurr3(                                      &
     &           layers,nv2,en,nk,jspin,lapw,                           &
     &           bkpts,sym,cell,gfinp,mpi)
!************************************************                       
!     Calculate the current by propagating the                          
!     embedding potentials to the planes by which                       
!     the system is divided into subsystems.                            
!     Frank Freimuth, November 2007                                     
!************************************************                       
      USE m_gf_iotmat 
      USE m_gf_types 
      USE m_gf_embedding 
      USE m_gf_propaemb 
      USE m_gf_writetrans,ONLY:writetrans 
      USE m_gf_propagate_embpot 
                                                                        
      IMPLICIT NONE 
      TYPE(t_layers),INTENT(IN)::layers 
      type(t_mpi),intent(in)   ::mpi
      INTEGER,INTENT(IN)::nv2 
      INTEGER,INTENT(IN)::en 
      INTEGER,INTENT(IN)::nk 
      INTEGER,INTENT(IN)::jspin 
      TYPE(t_lapw),INTENT(IN)::lapw 
      REAL,INTENT(IN)::bkpts(:,:) 
      TYPE(t_sym),INTENT(IN)::sym 
      TYPE(t_cell),INTENT(IN)::cell 
      TYPE(t_gfinp),INTENT(IN)::gfinp 
                                                                        
      COMPLEX,ALLOCATABLE:: g1(:,:),g2(:,:) 
      COMPLEX,ALLOCATABLE:: embpot_right(:,:,:),embpot_left(:,:,:) 
      INTEGER            :: layer 
      REAL               :: j(layers%num_layers+1) 
                                                                        
      ALLOCATE( g1(lapw%nv2_tot,lapw%nv2_tot) ) 
      ALLOCATE( g2(lapw%nv2_tot,lapw%nv2_tot) ) 
      ALLOCATE( embpot_left (nv2,nv2,1:layers%num_layers+1) ) 
      ALLOCATE( embpot_right(nv2,nv2,0:layers%num_layers) ) 
                                                                        
      CALL gf_getemb2(g1,1,1,en,nk,jspin,lapw) 
      CALL gf_getemb2(g2,2,layers%num_layers,en,nk,jspin,lapw) 
                                                                        
      embpot_left(:,:,1)=cmplx(-2.0,0.0)*g1 
!propagate left embedding potential                                     
      DO layer=1,layers%num_layers 
         CALL gf_propagate_embpot_left(layer,en,nk,jspin,               &
     &                                 lapw,embpot_left(:,:,layer),     &
     &                                 gfinp%l_nohelpregion,            &
     &                                 embpot_left(:,:,layer+1))        
      ENDDO 
      embpot_left=cmplx(-1.0,0.0)*embpot_left 
                                                                        
      embpot_right(:,:,layers%num_layers)=cmplx(2.0,0.0)*g2 
!propagate right embedding potential                                    
      DO layer=layers%num_layers,1,-1 
         CALL gf_propagate_embpot_right(layer,en,nk,jspin,              &
     &                                 lapw,embpot_right(:,:,layer),    &
     &                                 gfinp%l_nohelpregion,            &
     &                                 embpot_right(:,:,layer-1))       
      ENDDO 
                                                                        
      DO layer=1,layers%num_layers+1 
       CALL gf_landauer1plane(                                          &
     &                     embpot_left(:,:,layer),                      &
     &                     embpot_right(:,:,layer-1),                   &
     &                     j(layer)  )                                  
                                                                        
      ENDDO 
                                                                        
      CALL writetrans(en,nk,jspin,bkpts,sym,cell,3,j,mpi)
                                                                        
      END SUBROUTINE gf_propaembcurr3 
                                                                        
      SUBROUTINE gf_Landauer1Plane(                                     &
     &                             g1,g2,                               &
     &                             j)                                   
!c********************************************************************* 
!c     subroutine to calculate the current from two embedding           
!c     potentials on the same plane                                     
!c                                                                      
!c                                      Daniel Wortmann                 
!c********************************************************************* 
      USE m_gf_math 
      IMPLICIT NONE 
      COMPLEX,INTENT(IN)  :: g1(:,:),g2(:,:) 
      REAL,INTENT(OUT)    :: j 
                                                                        
      COMPLEX             :: G(size(g1,1),size(g1,1)) 
      COMPLEX :: A(SIZE(g1,1),SIZE(g1,1)),B(SIZE(g1,1),SIZE(g1,1)) 
                                                                        
      G=mat_inverse(G1+G2) 
      IF(.FALSE.)THEN 
        j = 2.0*real(trace(matmul(matmul(imag2d(G1),G),                 &
     &     matmul(imag2d(g2)                                            &
     &     ,TRANSPOSE(CONJG(G))))))                                     
                                                                        
      ELSE 
        ! longer version                                                
        ! Re(G1 G12 G2 G21^* - G1 G12 G2^* G21^*)                       
        A=matmul(matmul(g1,g),g2) 
        A=matmul(A,transpose(conjg(g))) 
        B=matmul(matmul(g1,g),transpose(conjg(g2))) 
        B=matmul(B,transpose(conjg(g))) 
                                                                        
        j=-2.*REAL(trace(A-B)) 
      ENDIF 
                                                                        
      END SUBROUTINE gf_landauer1plane 
      END                                           
