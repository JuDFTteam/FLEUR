!--------------------------------------------------------------------------------
! Copyright (c) 2026 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions 
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------
      MODULE m_gf_propagate_embpot 
      use m_juDFT
          IMPLICIT NONE
      CONTAINS 
      SUBROUTINE gf_propagate_embpot_right(                             &
     &                                    layer,en,nk,jspin,            &
     &                                    lapw,embpot_in,l_nohelpregion,&
     &                                    embpot_out)                   
!*********************************************                          
!     Propagate right embedding potential.                              
!     Frank Freimuth, November 2007                                     
!*********************************************                          
      USE m_gf_iotmat 
      USE m_gf_types 
      USE m_gf_propaemb 
      IMPLICIT NONE 
      INTEGER, INTENT(IN)     :: layer 
      INTEGER, INTENT(IN)     :: en,nk,jspin 
      TYPE(t_lapw),INTENT(IN) :: lapw 
      COMPLEX, INTENT(IN)     :: embpot_in (:,:) 
      LOGICAL, INTENT(IN)     :: l_nohelpregion 
      COMPLEX, INTENT(OUT)    :: embpot_out(:,:) 
                                                                        
      COMPLEX,ALLOCATABLE:: tmat(:,:) 
                                                                        
      ALLOCATE( tmat(2*lapw%nv2_tot,2*lapw%nv2_tot) ) 
                                                                        
      CALL gf_read_tmat2(                                               &
     &                   layer,en,nk,jspin,                             &
     &                   lapw,                                          &
     &                   tmat)                                          
      CALL gf_propaemb(.TRUE.,lapw%nv2_tot,embpot_in,tmat,embpot_out) 
                                                                        
      IF(.NOT.l_nohelpregion)THEN 
!          CALL juDFT_error("not yet implemented",calledby="gf_propagate_embpot.F90")
      ENDIF 
      DEALLOCATE(tmat) 
      END SUBROUTINE gf_propagate_embpot_right 
                                                                        
                                                                        
      SUBROUTINE gf_propagate_embpot_left(                              &
     &                                    layer,en,nk,jspin,            &
     &                                    lapw,embpot_in,l_nohelpregion,&
     &                                    embpot_out)                   
!*********************************************                          
!     Propagate left embedding potential.                               
!     Frank Freimuth, November 2007                                     
!*********************************************                          
      USE m_gf_iotmat 
      USE m_gf_types 
      USE m_gf_propaemb 
      IMPLICIT NONE 
      INTEGER, INTENT(IN)     :: layer 
      INTEGER, INTENT(IN)     :: en,nk,jspin 
      TYPE(t_lapw),INTENT(IN) :: lapw 
      COMPLEX, INTENT(IN)     :: embpot_in (:,:) 
      LOGICAL, INTENT(IN)     :: l_nohelpregion 
      COMPLEX, INTENT(OUT)    :: embpot_out(:,:) 
                                                                        
      COMPLEX,ALLOCATABLE:: tmat(:,:) 
                                                                        
      ALLOCATE( tmat(2*lapw%nv2_tot,2*lapw%nv2_tot) ) 
                                                                        
      CALL gf_read_tmat2(                                               &
     &                   layer,en,nk,jspin,                             &
     &                   lapw,                                          &
     &                   tmat)                                          
      CALL gf_propaemb(.FALSE.,lapw%nv2_tot,embpot_in,tmat,embpot_out) 
                                                                        
      IF(.NOT.l_nohelpregion)THEN 
!          CALL juDFT_error("not yet implemented",calledby="gf_propagate_embpot.F90")
      ENDIF 
      DEALLOCATE(tmat) 
      END SUBROUTINE gf_propagate_embpot_left 
      END                                           
