!--------------------------------------------------------------------------------
! Copyright (c) 2026 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions 
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------
      MODULE m_gf_iotmat 
      IMPLICIT NONE
!-----------------------------------------------                        
! DESC:Handle IO of t-matrix, option to keep it memory for all energies 
!                 Daniel Wortmann, (07-12-14)                           
!-----------------------------------------------                        
      PRIVATE 
      COMPLEX,ALLOCATABLE,SAVE :: tmat_storage(:,:,:,:) 
      INTEGER,SAVE             :: nk_storage,jspin_storage 
      PUBLIC gf_READ_tmat2,gf_WRITE_tmat,init_tmat_storage 
                                                                        
      CONTAINS 
      !<--S: gf_READ_tmat2(layer,en,nk,jspin,lapw_gf,tmat)                 
                                                                        
      SUBROUTINE gf_READ_tmat2(                                         &
     &     layer,en,nk,jspin,lapw_gf,                                      &
     &     tmat)                                                        
!******************************************                             
!    Read the T-matrices of the subsystems                              
!    and calculate the T-matrix of the                                  
!    composed system.                                                   
!    Frank Freimuth, November 2007                                      
!******************************************                             
      USE m_gf_types 
      USE m_gf_io2dmat 
      USE m_gf_tmatregularize 
      use m_juDFT 
      IMPLICIT NONE 
      INTEGER,INTENT(IN)      :: layer 
      TYPE(t_lapw_gf),INTENT(IN) :: lapw_gf 
      INTEGER,INTENT(IN)      :: en 
      INTEGER,INTENT(IN)      :: nk 
      INTEGER,INTENT(IN)      :: jspin 
      COMPLEX,INTENT(OUT)     :: tmat(:,:) 
      CALL timestart("IO: reading tmat") 
      IF (ALLOCATED(tmat_storage)) THEN 
         IF (nk /= nk_storage) CALL                                     &
     &        juDFT_error("io_tmat:Wrong nk in storage")                  
         IF (jspin /= jspin_storage) CALL                               &
     &        juDFT_error("io_tmat:Wrong nk in storage")                  
         tmat = tmat_storage(:SIZE(tmat,1),:SIZE(tmat,2),en,layer) 
      ELSE 
         IF(.NOT.(gf_read2dmat(                                         &
     &        IO2D_TMAT,layer,0,en,nk,jspin,lapw_gf,tmat )))               &
     &        CALL juDFT_error("tmat_layer: Data not found in file")      
      ENDIF 
      CALL gf_tmatregularize( tmat ) 
      CALL timestop("IO: reading tmat") 
      END SUBROUTINE gf_READ_tmat2 
                                                                        
      !>                                                                
      !<-- S: gf_write_tmat(layer,en,nk,jspin,lapw_gf,tmat)                
                                                                        
      SUBROUTINE gf_write_tmat(layer,en,nk,jspin,lapw_gf,tmat) 
!-----------------------------------------------                        
!  write t-matrix to storage                                            
!           (last modified: 2004-00-00) D. Wortmann                     
!-----------------------------------------------                        
      USE m_gf_types 
      USE m_gf_io2dmat 
      IMPLICIT NONE 
      !<--Arguments                                                     
      INTEGER,INTENT(IN)      :: layer 
      TYPE(t_lapw_gf),INTENT(IN) :: lapw_gf 
      INTEGER,INTENT(IN)      :: en 
      INTEGER,INTENT(IN)      :: nk 
      INTEGER,INTENT(IN)      :: jspin 
      COMPLEX,INTENT(IN)      :: tmat(:,:) 
                                                                        
      !>                                                                
      CALL timestart("IO: writing tmat") 
      IF (ALLOCATED(tmat_storage)) THEN 
         nk_storage = nk 
         jspin_storage=jspin 
         tmat_storage(:SIZE(tmat,1),:SIZE(tmat,2),en,layer) = tmat 
      ELSE 
      CALL gf_write2dmat(IO2D_TMAT,                                     &
     &     layer,0,en,nk,jspin,lapw_gf,Tmat(:,:))                          
      ENDIF 
      CALL timestop("IO: writing tmat") 
      END SUBROUTINE 
                                                                        
      !>                                                                
                                                                        
      !<-- S: init_tmat_storage(en,layer,nv2)                           
      SUBROUTINE init_tmat_storage(l_noco,energies,layers,nv2) 
!-----------------------------------------------                        
!allocate storage for t-matrix                                          
!           (last modified: 2004-00-00) D. Wortmann                     
!-----------------------------------------------                        
      IMPLICIT NONE 
      !<--Arguments                                                     
      LOGICAL,INTENT(IN)     :: l_noco 
      INTEGER,INTENT(IN)     :: energies,layers,nv2 
      !>                                                                
      IF (l_noco) THEN 
         ALLOCATE(tmat_storage(4*nv2,4*nv2,energies,layers)) 
                                                                        
      ELSE 
         ALLOCATE(tmat_storage(2*nv2,2*nv2,energies,layers)) 
      ENDIF 
      WRITE(*,*) "Keep T-matrix in memory" 
      END SUBROUTINE 
      !>                                                                
                                                                        
      END                                           
