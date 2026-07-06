!--------------------------------------------------------------------------------
! Copyright (c) 2026 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions 
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------
      MODULE m_gf_propaembcurr 
          IMPLICIT NONE
      CONTAINS 
      SUBROUTINE gf_propaembcurr(                                       &
     &           layers,nv2,en,nk,jspin,lapw,                           &
     &           bkpts,sym,cell,mpi)
!************************************************                       
!     Calculate the current by propagating the                          
!     embedding potentials to the planes by which                       
!     the system is divided into subsystems.                            
!     Frank Freimuth, November 2007                                     
!************************************************                       
      USE m_gf_read_tmat 
      USE m_gf_types 
      USE m_gf_embedding 
      USE m_gf_propaemb 
                                                                        
      IMPLICIT NONE 
      TYPE(t_layers),INTENT(IN)::layers 
      type(t_mpi),intent(in)::mpi

      INTEGER,INTENT(IN)::nv2 
      INTEGER,INTENT(IN)::en 
      INTEGER,INTENT(IN)::nk 
      INTEGER,INTENT(IN)::jspin 
      TYPE(t_lapw),INTENT(IN)::lapw 
      REAL,INTENT(IN)::bkpts(:,:) 
      TYPE(t_sym),INTENT(IN)::sym 
      TYPE(t_cell),INTENT(IN)::cell 
                                                                        
      COMPLEX,ALLOCATABLE:: tmat(:,:) 
      COMPLEX,ALLOCATABLE:: g1(:,:),g2(:,:) 
      COMPLEX,ALLOCATABLE:: embpot_in(:,:),embpot_out(:,:) 
                                                                        
      ALLOCATE( tmat(2*nv2,2*nv2) ) 
      ALLOCATE( g1(lapw%nv2_tot,lapw%nv2_tot) ) 
      ALLOCATE( g2(lapw%nv2_tot,lapw%nv2_tot) ) 
      ALLOCATE( embpot_in(nv2,nv2) ) 
      ALLOCATE( embpot_out(nv2,nv2) ) 
                                                                        
      CALL gf_read_tmat(                                                &
     &                          layers,nv2,en,nk,jspin,                 &
     &                          lapw,                                   &
     &                          tmat)                                   
                                                                        
      CALL gf_getemb(g1,g2,1,en,nk,jspin,lapw) 
      embpot_in=-2*g1 
      g2=2*g2 
      CALL gf_propaemb(.FALSE.,nv2,embpot_in,tmat,embpot_out) 
      embpot_out=-embpot_out 
      CALL gf_landauer1plane(en,nk,jspin,bkpts,sym,cell,                &
     &                        nv2,embpot_out,g2,mpi)
                                                                        
      DEALLOCATE( tmat ) 
      END SUBROUTINE gf_propaembcurr 
                                                                        
      SUBROUTINE gf_Landauer1Plane(en,nk,jspin,bkpts,sym,cell,nv2,g1,g2,mpi)
!c********************************************************************* 
!c     subroutine to calculate the current from two embedding           
!c     potentials on the same plane                                     
!c                                                                      
!c                                      Daniel Wortmann                 
!c********************************************************************* 
      USE m_gf_writetrans,ONLY:writetrans 
      USE m_gf_math 
      USE m_gf_types 
      IMPLICIT NONE
      type(t_mpi),intent(in)::mpi

      INTEGER,INTENT(IN)  :: en 
      INTEGER,INTENT(IN)  :: nk 
      INTEGER,INTENT(IN)  :: jspin 
      REAL,INTENT(IN)     :: bkpts(:,:) 
      TYPE(t_sym),INTENT(IN)            :: sym 
      TYPE(t_cell),INTENT(IN)           :: cell 
      INTEGER,INTENT(IN)  :: nv2 
      COMPLEX,INTENT(IN)  :: g1(:,:),g2(:,:) 
                                                                        
      REAL                :: j1(3+nv2) 
      COMPLEX             :: G(size(g1,1),size(g1,1)) 
      COMPLEX :: A(SIZE(g1,1),SIZE(g1,1)),B(SIZE(g1,1),SIZE(g1,1)) 
                                                                        
      G=mat_inverse(G1+G2) 
      j1(1) = 2.0*real(trace(matmul(matmul(imag2d(G1),G),               &
     &     matmul(imag2d(g2)                                            &
     &     ,TRANSPOSE(CONJG(G))))))                                     
                                                                        
      ! longer version                                                  
      ! Re(G1 G12 G2 G21^* - G1 G12 G2^* G21^*)                         
      A=matmul(matmul(g1,g),g2) 
      A=matmul(A,transpose(conjg(g))) 
      B=matmul(matmul(g1,g),transpose(conjg(g2))) 
      B=matmul(B,transpose(conjg(g))) 
                                                                        
      j1(1)=-2.*REAL(trace(A-B)) 
                                                                        
      CALL writetrans(en,nk,jspin,bkpts,sym,cell,3,j1,mpi)
                                                                        
                                                                        
      END SUBROUTINE gf_landauer1plane 
      END                                           
