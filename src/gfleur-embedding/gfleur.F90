!--------------------------------------------------------------------------------
! Copyright (c) 2026 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions 
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------
       PROGRAM gfleur 
!********************************************************************** 
!     This PROGRAM implements the embedding Green FUNCTION formalism    
!     on the basis of the FLEUR code                                    
!                                                                       
!     Daniel Wortmann 2000-2004                                         
!********************************************************************** 
!     Modifications in the call tree that allow the decomposition of    
!     the problem into smaller problems related to subsystems ("layers")
!     Frank Freimuth, November 2007                                     
!     CALL tree:                                                        
!     gfleur                                                            
!     |                                                                 
!     |-gf_setup                                                        
!     |                                                                 
!     (Self-Consistency loop it)                                        
!     it  |                                                             
!     it  |-gf_vgen                                                     
!     it  |                                                             
!     it  (Spin loop jspin)                                             
!     it  jspin  |                                                      
!     it  jspin  (kpoint-loop nk)                                       
!     it  jspin  nk   |                                                 
!     it  jspin  nk   (layer-loop layer)                                
!     it  jspin  nk    layer                                            
!     it  jspin  nk    layer                                            
!     it  jspin  nk    layer                                            
!*********************************************************************  
                                                                        
                                                                        
      !<-- Uses

      USE m_gf_vacuum_pot
      USE m_gf_vacuum_hs
      use m_gf_vacuum_charge
      USE m_gf_boundaryembpot 
      USE m_gf_charge 
      USE m_gf_optional 
      use m_juDFT
      USE m_fleur_interface,ONLY:fleur_calcenpara
      USE m_gf_energies,ONLY:gf_noen,gf_weightz,direction 
      USE m_gf_types 
      USE m_gf_PhaseMatrix,ONLY:initPhaseMatrix 
      USE m_gf_cdnskval 
      USE m_gf_gencdn 
      USE m_gf_vgen 
      USE m_gf_hs1lapw
      USE m_gf_dos
      USE m_gf_ab_coef 
      USE m_gf_io2dmat,ONLY:io2dmat_finalize => finalize 
      USE m_gf_setup 
      USE m_gf_mix 
      !USE m_gf_vacuum
      USE m_gf_apws 
      USE m_gf_version 
      USE m_gf_nognofinalize 
      USE m_gf_gmatfromtmat 
      USE m_gf_curvy2dprojector,ONLY: gf_curvy2dealloc 
      USE m_gf_gfleur_basic 
      USE m_gf_gfleur_propagate 
      USE m_gf_gfleur_compose 
      USE m_gf_tlmplminit 
      USE m_fleur_hsetup 
      USE m_gf_makepseudocharge 
      USE m_gf_iodop 
      USE m_gf_intcoul
      USE m_gf_totalmix
      USE m_gf_writetrans
      use m_gf_bandbending
      use m_gf_embpot_postprocess
      use m_gf_CBS,ONLY:outCBS_openfile,outCBS_closefile
      IMPLICIT NONE 
#ifdef CPP_MPI                                                          
      INCLUDE 'mpif.h' 
#endif                                                                  
                                                                        
      !>                                                                
      !<-- Vars for setup                                               
                                                                        
      INTEGER :: Jspins 
      TYPE(t_atoms),ALLOCATABLE    :: atoms(:) 
      TYPE(t_sphhar),ALLOCATABLE   :: sphhar(:) 
      TYPE(t_sym)                  :: sym 
      TYPE(t_cell),ALLOCATABLE     :: cell(:) 
      TYPE(t_xcpot),ALLOCATABLE    :: xcpot(:) 
      TYPE(t_mix)                  :: mix 
      TYPE(t_kpts)                 :: kpts 
      TYPE(t_enpara),ALLOCATABLE   :: enpara(:) 
      TYPE(t_soc)                  :: soc 
      TYPE(t_noco),ALLOCATABLE     :: noco(:) 
      TYPE(t_fermi)                :: fermi 
      TYPE(t_stars),ALLOCATABLE    :: stars(:) 
      TYPE(t_gfinp)                :: gfinp 
      TYPE(t_mpi)                  :: mpi 
      TYPE(t_lapw),ALLOCATABLE     :: lapw(:) 
      TYPE(t_layers)               :: layers 
      TYPE(t_potential),ALLOCATABLE:: potential(:) 
      TYPE(t_charge),ALLOCATABLE   :: charge(:) 
      TYPE(t_tlmplm),ALLOCATABLE   :: tlmplm_DATA(:) 
      TYPE(t_raddata),ALLOCATABLE  :: raddata(:) 
                                                                        
      !>                                                                
      !<-- integers for loops                                           
                                                                        
      INTEGER  :: it 
      INTEGER  :: layer_loop,nk_loop 
                                   !k-point loop                        
      INTEGER  :: nk 
                                     !spin                              
      INTEGER  :: jspin,nspins 
!     keeps track of DATA                                               
      LOGICAL  :: newspin 
      INTEGER  :: chargelayers 
                                                                        
      !>                                                                
      !<-- Mpi                                                          
      INTEGER               :: ierr 
      !>                                                                
      !<-- For timekeeping                                              
                                                                        
      CHARACTER(LEN = 80)      :: time_string 
                                                                        
      !>                                                                
      !<-- lapw a+b coef for construction of charge density             
                                                                        
      COMPLEX,ALLOCATABLE :: lapw_ab(:,:,:,:) 
                                                                        
      !>                                                                
      !<-- charge density                                               
                                                                        
                                      !max no of eigenstates per k-point
      INTEGER                :: neigd 
      LOGICAL                :: first_state,first_g 
      REAL   ,ALLOCATABLE    :: qtot_nuc(:),qtot_el(:) 
                                                                        
      !>                                                                
      !<-- File testing                                                 
                                                                        
      LOGICAL  :: l_exist_cdn
                                                                        
      !>                                                                
      !<-- potential                                                    
                                                                        
                                                !potential in aux-volume
      COMPLEX                   :: pot_aux(2,3) 
                                                !vacuum potential in    
      REAL                      :: vac_pot 
                                               !surface calculation     
                                                                        
      !>                                                                
      !<-- the matrices!!                                               
                                                                        
      COMPLEX, ALLOCATABLE :: g(:,:) 
      COMPLEX,ALLOCATABLE  :: gij(:,:,:),dgij(:,:,:) 
      COMPLEX,ALLOCATABLE  :: G_sum(:,:) 
      INTEGER              :: num_layers 
      INTEGER              :: layer 
                                                 !LDA+U                 
      COMPLEX,ALLOCATABLE  :: vs_mmp(:,:,:,:,:) 
      INTEGER, PARAMETER :: LMAXB = 3 
      !<-- init gf-part of FLEUR                                        
                                                                        
#ifdef CPP_MPI                                                          
      CALL MPI_INIT(ierr) 
#endif                                                                  
                      !Say hello                                        
      CALL gf_hello() 
      CALL gf_setup(layers,jspins,atoms,mpi,gfinp,cell,                 &
     &     sym,xcpot,mix,lapw,                                          &
     &     noco ,soc,fermi,stars,sphhar,enpara,kpts)                    
      CALL gf_OPTIONAL(jspins,cell(1),gfinp,lapw(1),noco(1),kpts,sym) 
                                                                        
      !<-- Allocate potential                                           
      ALLOCATE(potential(layers%num_layers)) 
      IF(gfinp%l_gmat)THEN 
         DO layer = 1,layers%num_layers 
            ALLOCATE( potential(layer)%vr(MAXVAL(Atoms(layer)%jri),     &
     &           0:MAXVAL(Sphhar(layer)%nlh),                           &
     &           Atoms(layer)%ntype,Jspins))                            
            IF (noco(1)%l_noco) THEN 
               ALLOCATE( potential(layer)%vpw(stars(layer)%nq3,3) ) 
            ELSE 
               ALLOCATE( potential(layer)%vpw(stars(layer)%nq3,jspins) ) 
            ENDIF 
            potential(layer)%vpw = 0.0 
         ENDDO 
      ENDIF 
      !>                                                                
      !<-- Allocate charge                                              
      ALLOCATE( charge(layers%num_layers) ) 
      IF(gfinp%l_charge)THEN 
         ALLOCATE( qtot_el(0:layers%num_layers)                         &
     &        ,qtot_nuc(0:layers%num_layers))                           
         DO layer_loop = 1,mpi%kl_layerPerPE 
            layer = mpi%kl_layers(layer_loop) 
            ALLOCATE( charge(layer)%rho_NEW(                            &
     &           MAXVAL(atoms(layer)%jri),                              &
     &           0:MAXVAL(sphhar(layer)%nlh),                           &
     &           atoms(layer)%ntype,jspins) )                           
            ALLOCATE( charge(layer)%pwd_NEW(                            &
     &           stars(layer)%nq3,jspins) )                             
            ALLOCATE( charge(layer)%qmtl_new(                           &
     &           0:MAXVAL(atoms(layer)%lmax),                           &
     &           atoms(layer)%ntype))                                   
            charge(layer)%qmtl_new = 0.0 
            charge(layer)%pwd_new  = 0.0 
            charge(layer)%rho_new  = 0.0 
         ENDDO 
      ENDIF 
      !>                                                                
      !<--SC-loop                                                       
                                                                        
                                      !self consistency loop
      CPP_juDFT_timestart("iteration loop")
      DO it = 1,mix%iter 
         WRITE(*,*) "Starting iteration: ",it," of ",MAX(1,mix%iter) 
         !<-- Generate total Hartree potential                          
         INQUIRE(FILE ="gf_cdn.hdf",EXIST= l_exist_cdn)
         if (.not.l_exist_cdn.and.gfinp%l_charge)then
              write(*,*) "WARNING, no gf_cdn but l_charge=T"
              write(6,*)  "WARNING, no gf_cdn but l_charge=T"
              if (mix%iter>1) call juDFT_error("No self-consistency withount gf_cdn.hdf")
         endif
         IF ((gfinp%l_charge.AND.l_exist_cdn).OR.(gfinp%l_CBS.AND.l_exist_cdn)) THEN
            CPP_juDFT_timestart("gf_makepseudocharge") 
            DO layer_loop = 1,mpi%kl_LayerperPE 
               layer  = mpi%kl_layers(layer_loop) 
                                         !only call makepseudocharge wit
               IF (mpi%k_kpts(1)==1) THEN 
                  CALL gf_makepseudocharge(layer,jspins,atoms(layer)    &
     &              ,stars(layer),sphhar(layer),cell(layer),sym         &
     &                 ,layers,noco(1),mpi)
               ENDIF 
            ENDDO 
            CPP_juDFT_timestop("gf_makepseudocharge")
         ENDIF 
         IF (gfinp%l_charge.AND.l_exist_cdn) THEN

            CPP_juDFT_timestart("gf_makeintcoulpot") 
            DO layer = 1,layers%num_layers 
               potential(layer)%vpw = 0.0 
               potential(layer)%vr= 0.0 
            ENDDO 
#ifdef CPP_MPI                                                          
               CALL MPI_BARRIER(MPI_COMM_WORLD,ierr) 
#endif                                                                  
               IF (mpi%pe0) THEN
                    CALL gf_makeintcoulpot(jspins,layers,stars,mpi,    &
     &           gfinp,potential,vac_pot,atoms,cell,sym)
               endif
            CPP_juDFT_timestop("gf_makeintcoulpot") 
#ifdef CPP_MPI                                                          
               CALL MPI_BARRIER(MPI_COMM_WORLD,ierr) 
#endif                                                                  
                                                                        
         ENDIF 
         !>                                                             
         !<-- Generate Potential                                        
         IF (gfinp%l_gmat) THEN 
            CPP_juDFT_timestart("gf_vgen") 
            WRITE(6,*) "Potential Generation" 
            CALL gf_vgen(layers,gfinp,                                  &
     &           .TRUE.,                                                &
     &           stars,sphhar,sym,cell,                                 &
     &           atoms,mpi,xcpot,mix,                                   &
     &           jspins,noco,enpara,                                           &
     &           pot_aux,potential)                                     
                         ! this is the vchk-switch
            if(gfinp%l_surface.and.l_exist_cdn)call gf_vacuum_totalpot(vac_pot,stars(1),xcpot(1),cell(1),&
                   jspins,noco(1)%l_noco,mpi%self_subcom)
            CPP_juDFT_timestop("gf_vgen") 
            WRITE(6,*) "Potential Generation done" 
            !if (gfinp%l_surface) call gf_vacuum_check_vacpot( &
            !potential(layers%num_layers)%vpw,stars(layers%num_layers),cell(layers%num_layers)%amat(3,3),vac_pot)


            IF (gfinp%l_cbs.AND.l_exist_cdn.and.mpi%isize==1) CALL gf_boundaryembpot(layers  &
     &           ,stars,cell,atoms,gfinp)



         !>                                                             

         !adjust energy parameters
         DO layer=1,layers%num_layers
            CALL fleur_calcenpara(jspins,mpi,enpara(layer)        &
     &                 ,atoms(layer),potential(layer)%vr)
         ENDDO
         ENDIF
         !<-- Spin loop                                                 
         nspins = jspins 
                                     !only one spin for IEC             
         IF (gfinp%l_IEC) nspins = 1 
         IF ( noco(1)%l_noco ) nspins = 1 
                               !spin loop                               
         DO jspin  =  1,nspins 
            !<--Setup for spinloop                                      
            IF (gfinp%l_dos)                                            &
     &           CALL gf_dos_init(layers,gfinp,atoms,cell(1),sym,lapw(1)&
     &           ,jspins,kpts%nkpts,potential,enpara,noco(1)%l_noco)
                                                                        
            !>                                                          
            IF (gfinp%l_CBS) THEN
                if (noco(1)%l_noco) THEN
                     CALL outCBS_openfile(jspin,jspins,kpts,sum(layers%c),2*lapw(1)%nv2d,gfinp%l_hdfio)
                else
                     CALL outCBS_openfile(jspin,jspins,kpts,sum(layers%c),2*lapw(1)%nv2d,gfinp%l_hdfio)
                endif
            endif

            IF ( ANY(atoms(:)%n_u>0) ) THEN 
               CPP_juDFT_timestart("fleur_usetup") 
               IF (.NOT.ALLOCATED(vs_mmp)) ALLOCATE(vs_mmp(-LMAXB:LMAXB,&
     &              -LMAXB:LMAXB,MAXVAL(atoms(:)%n_u),Jspins,           &
     &              layers%num_layers))                                 
               DO layer = 1,layers%num_layers 
                  CALL fleur_usetup(jspins,atoms(layer),mpi,            &
     &                 sphhar(layer),enpara(layer),                     &
     &                 potential(layer)%vr,vs_mmp(:,:,:,:,layer),layer)
               ENDDO 
               CPP_juDFT_timestop("fleur_usetup") 
            ELSE 
               IF (.NOT.ALLOCATED(vs_mmp)) ALLOCATE(vs_mmp(-LMAXB:LMAXB,&
     &              -LMAXB:LMAXB,1,Jspins,layers%num_layers))           
            ENDIF 
            !>                                                          
            !<-- k-loop                                                 
!            WRITE(*,*) "Parallel(kperPE)", mpi%irank, mpi%k_kperPE     
!            WRITE(*,*) "Parallel(kl_layerpe)", mpi%irank,              
!     $           mpi%kl_layerperPE                                     
                                        !k-loop                         
            DO nk_loop = 1,mpi%k_kperPE 
               nk=mpi%k_kpts(nk_loop) 
               DO layer_loop=1,mpi%kl_layerperPE 
                  layer = mpi%kl_layers(layer_loop)
	              IF (mpi%pe0) call priv_layer_progress('G0 for layer:',layer,layer==1,layer_loop==mpi%kl_layerperPE)

                  !<-- Setup of Projectors,APWs& Phasematrix            
!                   WRITE(*,"(i0,a,i0,a,i0,a,i0,a,i0,a,i0)") mpi%irank
!     $                 ,"k-loop:",nk_loop,"->",nk," Layerloop:"
!     $                 ,layer_loop,"->",layer                          
                  IF(gfinp%l_nohelpregion)THEN 
                     CALL gf_curvy2dealloc() 
                  ENDIF 
                  CPP_juDFT_timestart("gf_apws") 
                  CALL gf_APWS(Jspins,jspin,kpts%bk(:,nk),mpi%pe0,      &
     &                 noco(layer),gfinp,Cell(layer),lapw(layer),layer) 
                                                                        
                  !IF (gfinp%l_surface.AND.layer == layers%num_layers)   &
                  !     THEN
                  !   CALL gf_vacuum_init(jspin,cell(layer)              &
                  !        ,kpts,lapw(layer),gfinp,vac_pot,stars(layer),noco(layers%num_layers)%l_noco)
                  !   CALL gf_vacuum_Poti_init()
                  !ENDIF
                  CALL initPhaseMatrix(jspin,lapw(layer),cell(layer),   &
     &                           gfinp,noco(layer)%l_noco)              
                  CPP_juDFT_timestop("gf_apws") 
                  !>                                                    
                  !<-- Calculation of G_0                               
                  IF(gfinp%l_gmat)THEN 
                     CPP_juDFT_timestart("gf_gfleur_basic")
                     CALL gf_gfleur_basic(jspin,jspins,layer,layers,nk  &
     &                    ,gfinp,atoms(layer),sphhar(layer),stars(layer)&
     &                    ,sym,cell(layer),mpi,soc,noco(layer),kpts%bk(:&
     &                    ,:),lapw(layer),enpara(layer),potential(layer &
     &                    )%vpw,potential(layer)%vr,vs_mmp(:,:,:,:,layer&
     &                    ),tlmplm_DATA(layer),raddata(layer),pot_aux)
                     CPP_juDFT_timestop("gf_gfleur_basic")
                  ENDIF 
                  !>                                                    
                     !layer                                             
               ENDDO 
               IF (gfinp%l_surface) THEN
                 CPP_juDFT_timestart("gf_vacuum_diagonalize")
                 call gf_vacuum_diagonalize(mpi,lapw(mpi%kl_layers(1)),  &
                          stars(mpi%kl_layers(1)),cell(mpi%kl_layers(1)),&
                          sym,kpts,enpara(mpi%kl_layers(1)),jspins,nk,noco(mpi%kl_layers(1))%l_noco)
                 CPP_juDFT_timestop("gf_vacuum_diagonalize")
               endif

               !<-- Propagation of SIGMA                                
#ifdef CPP_MPI                                                          
               CALL MPI_BARRIER(MPI_COMM_WORLD,ierr) 
#endif                                                                  
               IF(gfinp%l_gmat.and.(((gfinp%l_dos.or.gfinp%l_charge).AND.layers%num_layers>1).OR.          &
     &              (gfinp%curr >3)))THEN
                  CPP_juDFT_timestart("gf_gfleur_propagate") 
                  CALL gf_gfleur_propagate(layers,mpi                   &
     &                 ,lapw(mpi%kl_layers(1)),gfinp,nk,jspin,sym,cell  &
     &                 ,kpts%bk)                                        
                  CPP_juDFT_timestop("gf_gfleur_propagate") 
               ENDIF 
               !>                                                       
#ifdef CPP_MPI                                                          
               CALL MPI_BARRIER(MPI_COMM_WORLD,ierr) 
#endif
               if (gfinp%curr>15) cycle !nk loop
               !<--Second layerloop for G_0->G                          
               chargelayers = 1;layer = 1 
               IF(gfinp%l_charge.or.gfinp%l_dos)chargelayers = layers%num_layers
               DO layer_loop = 1,MIN(chargelayers,mpi%kl_layerperPE) 
                  IF (gfinp%l_gmat.and.((gfinp%l_dos.or.gfinp%l_charge).OR.mpi%kl_layers(layer_loop) == 1)) &
     &                 THEN
                     layer = mpi%kl_layers(layer_loop) 
                     IF (.not.(gfinp%l_charge.or.gfinp%l_doslayer(layer))) cycle
                     CPP_juDFT_timestart("gf_gfleur_compose")
                     IF (mpi%pe0) call priv_layer_progress('G for layer:',layer,layer==1,layer_loop==mpi%kl_layerperPE)

                     CALL gf_gfleur_compose(layer,noco(layer),gfinp     &
     &                    ,layers,nk,jspin,sym,cell(layer),mpi          &
     &                    ,lapw(layer),kpts,pot_aux,charge(layer)       &
     &                    ,atoms(layer),stars(layer),sphhar(layer)      &
     &                    ,potential(layer)%vr,enpara(layer))
                     CPP_juDFT_timestop("gf_gfleur_compose")
                  ENDIF 
                     !layer                                             
               ENDDO 
               if (gfinp%l_surface.and.mpi%pe0.and.gfinp%l_charge.and.gfinp%l_gmat) then
               		call gf_vacuum_makecharge(stars(1),lapw(1),sym,cell(1),kpts,enpara(1),layers,noco(1)%l_noco,jspin,jspins,nk)
               endif
               !>                                                       
               WRITE(*,"(a,i5,a,i5,a,i1,a,i1,a,i2,a,i2,a)") "Finished: k:(" ,nk,"/"&
     &              ,kpts%Nkpts,")  Spin: (",jspin,"/",jspins,") Iter: (",it,"/",mix%iter,")"
               if (gfinp%l_embspectral) CALL gf_embpot_spectral(jspin,nk,gfinp,layers,lapw(1),cell(1),kpts)      !postprocessing in m_gf_embpot_postproc


            ENDDO !nk
            !>                                                          
            !<-- Construct new charge                                   
            neigd = 0 
            IF (gfinp%l_charge.and.gfinp%l_gmat) THEN
               qtot_el = 0.0 
               qtot_nuc = 0.0 
               DO layer_loop = 1,mpi%kl_layerperPE 
                  layer = mpi%kl_layers(layer_loop) 
                  CPP_juDFT_timestart("gf_gencdn") 
                                                                        
                  CALL gf_gencdn(layer,jspin,.FALSE.,jspins,            &
     &                 gfinp,atoms(layer),cell(layer),sym,kpts          &
     &                 ,stars(layer),sphhar(layer),mpi,enpara(layer)    &
     &                 ,potential(layer)%vr,charge(layer)%pwd_NEW       &
     &                 ,charge(layer)%rho_NEW,charge(layer)%qmtl_NEW    &
     &                 ,neigd,soc,noco(1)%l_noco,qtot_el(1:),qtot_nuc(1:))
                                                                        
                  CPP_juDFT_timestop("gf_gencdn") 
                  charge(layer)%qmtl_NEW = 0.0 
                     !layers                                            
               ENDDO 
            ENDIF 
            !>                                                          
            IF (gfinp%l_CBS) CALL outCBS_closefile(noco(1)%l_noco.or.(jspin==jspins))
               !spin loop                                               
         ENDDO 
         !>                                                             
         !<-- mix charge                                                
         IF (gfinp%l_charge.and.l_exist_cdn) THEN
            CPP_juDFT_timestart("gf_mix") 
            CALL gf_charge_sum(gfinp%l_surface,qtot_nuc,qtot_el)
            if (gfinp%l_surface) then
                      IF (mpi%pe0) WRITE(6,"('Vacuum charge:',f9.4)")  gf_vacuum_totalcharge(jspins,cell(1),stars(1))
                      qtot_el(0)=qtot_el(0)+gf_vacuum_totalcharge(jspins,cell(1),stars(1))
            endif
            IF (mpi%pe0) WRITE(6,"('Total charge neutrality:',f9.4,'+',f9.4,'=',f9.4)")          &
     &       -1*qtot_el(0),qtot_nuc(0), -qtot_el(0)+qtot_nuc(0)
            DO layer_loop = 1,mpi%kl_LayerperPE 
               layer  = mpi%kl_layers(layer_loop) 
               IF (mpi%k_kpts(1) == 1) THEN 
                  WHERE( ABS(qtot_el)<1E-10) 
                     qtot_el = 10E-10 
                  END WHERE 
                  CALL gf_charge(jspins,mpi,stars(layer),atoms(layer)   &
     &                 ,cell(layer),mix,gfinp,noco(layer),sphhar(layer) &
     &                 ,sym,charge(layer)%pwd_NEW(:,:),charge(layer     &
     &                 )%rho_NEW(:,:,:,:),layer,qtot_nuc/qtot_el)       
               ENDIF 
               charge(layer)%pwd_NEW   = 0.0 
               charge(layer)%rho_NEW   = 0.0 
               charge(layer)%qmtl_new  = 0.0 
            ENDDO 
            IF (.NOT.gfinp%l_potmix.AND.gfinp%l_charge.and.l_exist_cdn) THEN
               IF (gfinp%l_totalmix) THEN 
                  CALL gf_totalmix(jspins,layers,atoms,sphhar,stars,mix &
     &                 ,gfinp,mpi)
               ELSE 
                  CALL gf_mix(atoms,cell,sphhar                         &
     &                 ,stars,gfinp,noco,mpi,sym,enpara,mix,jspins,layers)
               ENDIF 
               !if (gfinp%l_surface.and.mpi%pe0) call gf_vacuum_copycharge(mpi%self_subcom,mix%alpha)
               CPP_juDFT_timestop("gf_mix") 
            ENDIF 
         ENDIF 
         !>                                                             
#ifdef CPP_MPI                                                          
         CALL MPI_barrier(MPI_COMM_WORLD,ierr) 
         IF (mpi%pe0) WRITE(*,*) "Finished IT:",it 
#endif                                                                  
            !self consistency loop                                      
      ENDDO
      CPP_juDFT_timestop("iteration loop")

      if (layers%num_layers==1.and.gfinp%l_CBS)   &
        call gf_bandbending(lapw(1),jspins,kpts%nkpts)

      !>
      CPP_juDFT_timestart("gf_writedos")
      IF (gfinp%l_dos) THEN
        IF (noco(1)%l_noco) THEN
         CALL gf_writedos(mpi,layers,atoms,cell,gfinp,3,kpts%nkpts)
        ELSE
         CALL gf_writedos(mpi,layers,atoms,cell,gfinp,jspins,kpts%nkpts)
        ENDIF
       ENDIF
      CPP_juDFT_timestop("gf_writedos")
      CPP_juDFT_timestart("io2dmat_finalize")
      CALL io2dmat_finalize() 
      call gf_writetrans_hdfclose()
      CPP_juDFT_timestop("io2dmat_finalize")
      if (mpi%pe0) CALL writetimes()
      CALL juDFT_end("GFLEUR done")
      contains
      subroutine priv_layer_progress(text,layer,first,last)
      implicit none
      !write out some information on the progress using nonadvancing IO
      character(len=*),intent(in)::text
      integer,intent(in)         ::layer
      logical,intent(in)         ::first,last

      if (first) write(*,"(A)",advance='no') text
      if (last) then
           write(*,"(i3)") layer
      else
           write(*,"(i3)",advance='no') layer
      endif
      end subroutine
      END                                           
