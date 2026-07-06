!--------------------------------------------------------------------------------
! Copyright (c) 2026 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions 
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------
      MODULE m_gf_tmat 
          IMPLICIT NONE
      CONTAINS 
      SUBROUTINE gf_tmat_layers(                                        &
     &                  layer,layers,en,nk,jspin,                       &
     &                  mpi,lapw,lapw_gf,gfinp,pot_aux,gij,                     &
     &                  tmat)                                           
!********************************************************************** 
!     Support for dividing the system into subsystems ("layers"):       
!     Calculate the T-matrix of subsystem with index                    
!     "layer". T-matrix may be written to file.                         
!     Frank Freimuth, November 2007                                     
!********************************************************************** 
      USE m_gf_types
      USE m_gf_maketmat,  ONLY: gf_maketmat_neumann,gf_maketmat2 
      USE m_gf_iotmat 
      USE m_gf_math 
      !>                                                                
      IMPLICIT NONE 
      !<-- Arguments                                                    
                                                                        
      INTEGER,      INTENT(IN)             :: layer 
      TYPE(t_layers),INTENT(IN)            :: layers 
      INTEGER,      INTENT(IN)             :: nk,en,jspin 
      TYPE(t_gfmpi),  INTENT(IN)             :: mpi 
      TYPE(t_embinp),INTENT(IN)             :: gfinp 
      TYPE(t_lapw), INTENT(IN)             :: lapw 
      TYPE(t_lapw_gf),INTENT(IN) :: lapw_gf
      COMPLEX,      INTENT(IN)             :: pot_aux 
      COMPLEX,      INTENT(IN)             :: gij(:,:) 
      COMPLEX,      INTENT(INOUT),OPTIONAL :: tmat(:,:) 
                                                                        
      !>                                                                
      !<-- Locals                                                       
      COMPLEX, ALLOCATABLE :: t1(:,:),T2(:,:) 
      REAL                 :: helpthick 

                                                                        
      ALLOCATE( t1(2*lapw_gf%nv2_tot,2*lapw_gf%nv2_tot) ) 
      CALL gf_maketmat_neumann(                                         &
     &                    lapw_gf%nv2_tot,gij,                             &
     &                    t1)                                           
      IF( (.NOT.gfinp%l_nohelpregion) .AND. (layer/=1) )THEN 
            helpthick=0.5*(layers%d(layer-1)-layers%c(layer-1)+         &
     &                      layers%d(layer)  -layers%c(layer)   )       
            ALLOCATE( t2(2*lapw_gf%nv2_tot,2*lapw_gf%nv2_tot) ) 
            CALL gf_maketmat2(layer,helpthick,lapw_gf%nv2_tot,lapw_gf%rkp, &
     &                         t2(:,:),gfinp,en,pot_aux)                
            t1=zmat_product(t1,t2) 
            DEALLOCATE( t2 ) 
      ENDIF 
      CALL gf_write_tmat(layer,en,nk,jspin,lapw_gf,T1(:,:)) 
                                                                        
      IF(present(tmat))THEN 
         tmat=t1 
      ENDIF 
      DEALLOCATE( t1 ) 
      END SUBROUTINE gf_tmat_layers 
                                                                        
      !<-- S:GF_TMAT(en,nk,jspin,sym,cell,mpi,kp,bk,gfinp,gij,dgij)     
      SUBROUTINE GF_TMAT(                                               &
     &     layer,l_layers,layers,en,nk,jspin,sym,cell,                  &
     &     mpi,lapw,lapw_gf,bk,gfinp,pot_aux,                                   &
     &     gij,dgij,on_the_fly)                                         
!********************************************************************** 
!     This subroutine is called for all t-matrix related calculations   
!     - gf_maketmat is called to actually construct the t-matrix        
!     - gf_CBS is called for CBS calculations incl. the construction of 
!        embedding potentials                                           
!     *                                                                 
!     *                           Daniel Wortmann, Tokyo, 2001          
!     (last modified: 04-09-07)                                         
!********************************************************************** 
!     Support for dividing the system into subsystems ("layers"):       
!     l_layers=.true. => read in the T-matrices of the subsystems and   
!     multiply them in order to calculate CBS and/or embedding potential
!     from their product corresponding to the T-matrix of the composed  
!     system.                                                           
!     Frank Freimuth, November 2007                                     
!********************************************************************** 
                                                                        
      !<-- Use                                                          
      use m_juDFT 
      USE m_gf_types
      USE m_gf_CBS,       ONLY: gf_CBS 
      USE m_gf_maketmat,  ONLY: gf_maketmat_neumann,gf_maketmat2 
      USE m_gf_io2dmat 
      USE m_gf_PhaseMatrix,ONLY:getLargePhaseMatrix 
      USE m_gf_rt_coefs, ONLY:gf_wave_match 
      USE m_gf_curvy2dprojector 
      USE m_gf_tmatregularize 
      USE m_gf_math,ONLY:zmat_product 
      USE m_gf_READ_tmat 
      USE m_gf_iotmat 
      USE m_gf_energies
      !>                                                                
      IMPLICIT NONE 
      !<-- Arguments                                                    
      INTEGER,INTENT(IN)       :: layer 
      LOGICAL, INTENT(IN)      :: l_layers 
      TYPE(t_layers),INTENT(IN):: layers 
      INTEGER,      INTENT(IN) :: nk,en,jspin 
      REAL,         INTENT(IN) :: bk(:,:) 
      TYPE(t_gfmpi),  INTENT(IN) :: mpi 
      TYPE(t_embinp),INTENT(IN) :: gfinp 
      TYPE(t_lapw),INTENT(IN)  :: lapw 
      TYPE(t_lapw_gf),INTENT(IN) :: lapw_gf
      TYPE(t_cell),INTENT(IN)  :: cell 
      TYPE(t_sym),INTENT(IN)   :: sym 
      COMPLEX,INTENT(IN)       :: pot_aux 
      COMPLEX,    INTENT(INOUT),OPTIONAL :: gij(:,:,:) 
      COMPLEX,INTENT(INOUT)    ,OPTIONAL :: dgij(:,:,:) 
      LOGICAL,INTENT(IN),OPTIONAL        :: on_the_fly 
                                                                        
      !>                                                                
      !<-- Locals                                                       
      INTEGER                  :: layer_ind 
      COMPLEX, ALLOCATABLE :: tr(:,:,:),T(:,:) 
      INTEGER i,j 
      REAL time1,time2 
      REAL helpthick 
                                                                        
                                                                        
      ALLOCATE( tr(2*lapw_gf%nv2_tot,2*lapw_gf%nv2_tot,2) ) 
      ALLOCATE(  T(2*lapw_gf%nv2_tot,2*lapw_gf%nv2_tot) ) 
                       !read in T-matrices of subsystems                
      IF(l_layers)THEN 
         CALL gf_READ_tmat(layers,lapw_gf%nv2_tot,en,nk,jspin,lapw_gf,tr(:,:,1&
     &        ))                                                        
         IF(.NOT.gfinp%l_nohelpregion)THEN 
             helpthick=0.5*( layers%d(1)-layers%c(1)+                   &
     &         layers%d(layers%num_layers)-layers%c(layers%num_layers)) 
             CALL gf_maketmat2(layer,helpthick,lapw_gf%nv2_tot,lapw_gf%rkp,&
     &         TR(:,:,2),gfinp,en,pot_aux)                              
         ENDIF 
      ELSE 
        ! T-matrix might be read-in from a file                         
         IF (.NOT.PRESENT(gij)) THEN 
            CALL gf_READ_tmat2(                                         &
     &           layer,en,nk,jspin,lapw_gf,                             &
     &           TR(:,:,1))                                             
         ELSE 
            IF(.NOT.(gf_read2dmat(IO2D_TMAT,layer,0,en,nk,jspin,lapw_gf &
     &           ,TR(:,:,1))).AND. gfinp%l_nohelpregion) THEN           
               !<-- T-Matrix must be calculated
               CALL timestart("gf_maketmat")
               CALL gf_maketmat_neumann(                                &
     &              lapw_gf%nv2_tot,Gij,                                   &
     &              TR(:,:,1))                                          
               IF(.NOT.gfinp%l_nohelpregion)THEN 
                  helpthick = layers%d(1)-layers%c(1) 
                  CALL gf_maketmat2(layer,helpthick,lapw_gf%nv2_tot        &
     &                 ,lapw_gf%rkp,TR(:,:,2),gfinp,en,pot_aux)         
               ENDIF
               CALL timestop("gf_maketmat")
               !<-- T-Matrices can be written to a file                 
               CALL gf_write2dmat(IO2D_TMAT,1,0,en,nk,jspin,lapw_gf,TR(:,: &
     &              ,1))                                                
               IF(.NOT.gfinp%l_nohelpregion)THEN 
                  CALL gf_write2dmat(IO2D_TMAT,2,0,en,nk,jspin,         &
     &                 lapw_gf,TR(:,:,2))                                  
               ENDIF 
               !>                                                       
            ENDIF 
            !>                                                          
         ENDIF 
      ENDIF 

      t(:,:)=zmat_product(getLargePhaseMatrix(),tr(:,:,1)) 
      IF(gfinp%l_nohelpregion)THEN 
         tr(:,:,2)=cmplx(0.0,0.0) 
         DO i=1,2*lapw_gf%nv2_tot 
            tr(i,i,2)=cmplx(1.0,0.0) 
         ENDDO 
      ELSE 
         t(:,:)=zmat_product(t(:,:),tr(:,:,2)) 
      ENDIF 
                                                                        
      !<-- Calculate CBS
      CALL timestart("gf_CBS")
      IF ( gfinp%l_cbs) THEN 
         CALL gf_CBS(layers,0,nk,bk(:2,nk),en,jspin,gfinp,lapw,lapw_gf, &
     &        cell, t(:,:),sym%mrot,mpi,TR(:,:,1),TR(:,:,2))            
      ELSE 
         IF (PRESENT(on_the_fly)) THEN 
            IF (on_the_fly)                                             &
     &           CALL gf_CBS(layers,layer,nk,bk(:2,nk),en,jspin,gfinp   &
     &           ,lapw,lapw_gf,cell, t(:,:),sym%mrot,mpi,TR(:,:,1),TR(:,:,2))   
         ENDIF 
      ENDIF 
      CALL timestop("gf_CBS")
      !>                                                                
                                                                        
      !<-- Do a wavefunction matching using the T-matrix                
      IF (gfinp%curr ==-2) THEN 
         CALL gf_wave_match(lapw_gf%nv2_tot,en,nk,jspin,sym,cell,lapw,   &
     &        lapw_gf,Tr(:,:,1),gfinp,bk,3,mpi)
      ENDIF 
      !>                                                                
                                                                        
      DEALLOCATE(tr,T) 
                                                                        
      END SUBROUTINE gf_tmat 
      !>                                                                
      END                                           
