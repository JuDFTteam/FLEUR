!--------------------------------------------------------------------------------
! Copyright (c) 2026 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions 
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------
      MODULE m_gf_gfleur_propagate 
      use m_juDFT
      IMPLICIT NONE
!*************************************************                      
!     Propagate the embedding potentials of the                         
!     composed system to the subsystems.                                
!     Might require a lot of disk space!                                
!     Frank Freimuth, November 2007                                     
!*************************************************                      
      CONTAINS 
      SUBROUTINE gf_gfleur_propagate(layers,fmpi,lapw,lapw_gf,gfinp,nk,jspin,ld &
     &     ,bk)
      USE m_gf_types
      USE m_gf_iotmat
      USE m_gf_propagate_embpot
      USE m_gf_embedding
      USE m_gf_io2dmat
      USE m_gf_energies
      USE m_gf_current_single
      USE m_gf_tmat
      USE m_gf_math
      IMPLICIT NONE
      TYPE(t_layers),INTENT(IN) :: layers
      TYPE(t_gfmpi),INTENT(IN)  :: fmpi
      TYPE(t_lapw),INTENT(IN)   :: lapw
      TYPE(t_lapw_gf),INTENT(IN) :: lapw_gf
      TYPE(t_embinp),INTENT(IN)  :: gfinp
      INTEGER,INTENT(IN)        :: nk
      INTEGER,INTENT(IN)        :: jspin
      TYPE(t_gflayer),INTENT(IN) :: ld(:)
      REAL   ,INTENT(IN)        :: bk(:,:)
                                                                        
      COMPLEX,ALLOCATABLE :: g1(:,:),g2(:,:) 
      COMPLEX,ALLOCATABLE :: embpot_right(:,:,:),embpot_left(:,:,:) 
      INTEGER             :: layer 
      LOGICAL             :: l_notused 
      INTEGER             :: key,en,nv2,en_loop 
      print *, "Propagate with rank:",fmpi%fmpi%irank
      key = IO2D_EMB 
      ALLOCATE( g1(lapw_gf%nv2_tot,lapw_gf%nv2_tot) ) 
      ALLOCATE( g2(lapw_gf%nv2_tot,lapw_gf%nv2_tot) ) 
      nv2=lapw_gf%nv2_tot 
      ALLOCATE( embpot_left (nv2,nv2,1:layers%num_layers) ) 
      ALLOCATE( embpot_right(nv2,nv2,1:layers%num_layers) ) 
      CALL gf_io2dmat_memorysync()
#ifdef CPP_TMAT_DIRECT                                                  
      DO layer = 1,layers%num_layers 
         REWIND(700+layer) 
      ENDDO 
#endif                                                                  
      IF (gf_bias_layer()<layers%num_layers+1) THEN 
         WRITE(*,*) "Generating complex BS" 
         DO en_loop = 1,fmpi%ke_ENperPE 
            en = fmpi%ke_energies(en_loop) 
            CALL GF_TMAT(                                               &
     &           gf_bias_layer(),.FALSE.,layers,en,nk,jspin,            &
     &           ld(1)%fi%sym,ld(gf_bias_layer())%fi%cell,fmpi,lapw,     &
     &           lapw_gf,bk,gfinp,CMPLX(0.0,0.0),on_the_fly = .TRUE.)
         ENDDO
      ENDIF 
                                                                        
      CALL timestart("gf_propagate_left")
      DO en_loop = 1,fmpi%ke_ENperPE
         en = fmpi%ke_energies(en_loop) 
                                                                        
         CALL gf_GETEMB2(G1,1,1,en,nk,jspin,lapw,lapw_gf) 

         embpot_left(:,:,1) = CMPLX(-2.0,0.0)*g1
         !propagate left embedding potential

        DO layer=1,layers%num_layers-1 
           IF (layer == gf_bias_layer()-1) THEN 
              CALL gf_GETEMB2(embpot_left(:,:,layer+1),gf_bias_layer(),1&
     &             ,en,nk,jspin,lapw,lapw_gf)
               embpot_left(:,:,layer+1)=CMPLX(-2.0,0.0)*embpot_left(:,:,layer+1)
!              g1 = embpot_left(:,:,layer)-imag2d(embpot_left(:,:,layer)
!              WRITE(*,*) layer,en                                      
!              CALL gf_propagate_embpot_left(layer,en,nk,jspin,         
!     >             lapw,g1,                                            
!     >             gfinp%l_nohelpregion,                               
!     >             embpot_left(:,:,layer+1))                           
           ELSE 
           !!$omp parallel sections default(shared)
           !!$omp section
              CALL gf_propagate_embpot_left(layer,en,nk,jspin,          &
     &             lapw,lapw_gf,embpot_left(:,:,layer),                         &
     &             gfinp%l_nohelpregion,                                &
     &             embpot_left(:,:,layer+1))
           !!$omp section
!              CALL gf_write2dmat(IO2D_EMB,layer,1,en,nk,jspin,lapw_gf,        &
!     &                         cmplx(-0.5,0.0)*embpot_left(:,:,layer))
           !!$omp end parallel sections
           ENDIF 
        ENDDO 
!        CALL gf_write2dmat(IO2D_EMB,layers%num_layers,1,en,nk,jspin,lapw_gf,        &
!     &                    cmplx(-0.5,0.0)*embpot_left(:,:,layers%num_layers))

        embpot_left=cmplx(-0.5,0.0)*embpot_left
        IF (gfinp%curr<16) THEN
          DO layer = 1,layers%num_layers
             CALL gf_write2dmat(IO2D_EMB,layer,1,en,nk,jspin,lapw_gf,  &
     &                         embpot_left(:,:,layer))
          ENDDO
        ENDIF

        if (gfinp%l_surface.or.gfinp%curr>15) then
          CALL gf_propagate_embpot_left(layers%num_layers,en,nk,jspin,          &
     &             lapw,lapw_gf,-2.*embpot_left(:,:,layers%num_layers),                         &
     &             gfinp%l_nohelpregion,                                &
     &              embpot_left(:,:,1))
     		  IF (gfinp%curr<16) THEN
     		        CALL gf_write2dmat(IO2D_EMB,layers%num_layers+1,1,en,nk,jspin,lapw_gf,        &
     &                         -0.5*embpot_left(:,:,1))
              ELSE
                 CALL gf_current_single(layers,lapw,lapw_gf,ld(1)%fi%cell,ld(1)%fi%sym,fmpi,bk,nk,en,jspin,        &
     &                         -0.5*embpot_left(:,:,1))
              ENDIF
        endif
      ENDDO 
      CALL timestop("gf_propagate_left")
      IF (gfinp%curr>15) return
#ifdef CPP_TMAT_DIRECT                                                    
      DO layer = 1,layers%num_layers 
         REWIND(700+layer) 
      ENDDO 
#endif
      CALL timestart("gf_propagate_right")
      DO en_loop = 1,fmpi%ke_ENperPE 
         en = fmpi%ke_energies(en_loop) 
         CALL gf_GETEMB2(G2,2,layers%num_layers,en,nk,jspin,lapw,lapw_gf)


         embpot_right(:,:,layers%num_layers) = CMPLX(2.0,0.0)*g2 
!propagate right embedding potential                                    

        DO layer=layers%num_layers,2,-1 
           IF (layer == gf_bias_layer()+1) THEN 
              CALL gf_GETEMB2(embpot_right(:,:,layer-1),gf_bias_layer() &
     &             ,2,en,nk,jspin,lapw,lapw_gf)
               embpot_right(:,:,layer-1)=CMPLX(2.0,0.0)*embpot_right(:,:,layer-1)
!              g1 = embpot_right(:,:,layer)-imag2d(embpot_right(:,:,laye
!     +             ))                                                  
!              CALL gf_propagate_embpot_right(layer,en,nk,jspin,        
!     >             lapw,g1,                                            
!     >             gfinp%l_nohelpregion,                               
!     >             embpot_right(:,:,layer-1))                          
           ELSE 
           !!$omp parallel sections default(shared)
           !!$omp section
              CALL gf_propagate_embpot_right(layer,en,nk,jspin,         &
     &             lapw,lapw_gf,embpot_right(:,:,layer),                        &
     &             gfinp%l_nohelpregion,                                &
     &             embpot_right(:,:,layer-1))
            !!$omp section
     !         CALL gf_write2dmat(IO2D_EMB,layer,2,en,nk,jspin,lapw_gf,        &
     !&          cmplx(0.5,0.0)*embpot_right(:,:,layer))
            !!$omp end parallel sections
           ENDIF 
        ENDDO 
        !CALL gf_write2dmat(IO2D_EMB,1,2,en,nk,jspin,lapw_gf,        &
     !&          cmplx(0.5,0.0)*embpot_right(:,:,1))
        embpot_right=cmplx(0.5,0.0)*embpot_right
        DO layer = 1,layers%num_layers
           CALL gf_write2dmat(IO2D_EMB,layer,2,en,nk,jspin,lapw_gf,   &
     &          embpot_right(:,:,layer))
        enddo
      ENDDO 
      CALL timestop("gf_propagate_right")
      CALL gf_io2dmat_memorysync()
      END SUBROUTINE gf_gfleur_propagate 
      END                                           
