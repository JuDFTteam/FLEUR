!--------------------------------------------------------------------------------
! Copyright (c) 2026 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions 
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------
      MODULE m_gf_read_tmat 
          IMPLICIT NONE
      CONTAINS 
      SUBROUTINE gf_read_tmat(                                          &
     &                       layers,nv2,en,nk,jspin,lapw,               &
     &                       tmat)                                      
!******************************************                             
!    Read the T-matrices of the subsystems                              
!    and calculate the T-matrix of the                                  
!    composed system.                                                   
!    Frank Freimuth, November 2007                                      
!******************************************                             
      USE m_gf_types 
      USE m_gf_iotmat 
      IMPLICIT NONE 
      TYPE(t_layers),INTENT(IN)::layers 
      TYPE(t_lapw_gf),INTENT(IN)::lapw 
      INTEGER,INTENT(IN):: nv2 
      INTEGER,INTENT(IN):: en 
      INTEGER,INTENT(IN):: nk 
      INTEGER,INTENT(IN):: jspin 
      COMPLEX,INTENT(OUT):: tmat(:,:) 
                                                                        
      COMPLEX, ALLOCATABLE :: tr(:,:,:),t_temp(:,:) 
      INTEGER layer_ind 
                                                                        
      ALLOCATE( tr(2*Nv2,2*Nv2,2) ) 
      ALLOCATE( t_temp(2*nv2,2*nv2) ) 
      CALL gf_READ_tmat2(1,en,nk,jspin,lapw,tr(:,:,1)) 
                                                                        
      DO layer_ind=2,layers%num_layers 
         CALL gf_READ_tmat2(layer_ind,en,nk,jspin,lapw,tr(:,:,2)) 
        CALL zgemm('N','N',2*nv2,2*nv2,                        &
     &      2*nv2,cmplx(1.0,0.0),tr(:,:,2),2*nv2,                       &
     &      tr(:,:,1),2*nv2,cmplx(0.0,0.0),t_temp,                      &
     &             2*nv2)                                               
           tr(:,:,1)=t_temp 
      ENDDO 
      DEALLOCATE( t_temp ) 
      tmat(:,:)=tr(:,:,1) 
                                                                        
      END SUBROUTINE gf_read_tmat 
      END                                           
