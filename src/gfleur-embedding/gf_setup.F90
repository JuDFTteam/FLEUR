!--------------------------------------------------------------------------------
! Copyright (c) 2026 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions 
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------
      MODULE m_gf_setup 
      use m_juDFT
      IMPLICIT NONE
      PRIVATE 
      PUBLIC :: gf_setup 
      CONTAINS 
      !<-- S: gf_setup                                                  
                                                                        
      SUBROUTINE gf_setup(layers,jspins,atoms,mpi,gfinp,cell,           &
     &      sym,xcpot,mix,lapw,                                         &
     &      noco,soc,fermi,stars,sphhar,enpara,kpts)                    
      use m_gf_stepsanaly, only: gf_calcstepsanaly, gf_initstepsanaly
      USE m_gf_types 
      USE m_gf_energies, gf_energies_init => init 
      USE m_gf_iosetup 
#include "juDFT_env.h" 
      USE m_gf_inp2plot 
      USE m_gf_iogfinp,ONLY:gf_rgfinp 
      USE m_gf_kpts,ONLY:gf_loadkpts 
      USE m_gf_stars,ONLY:gf_stargen 
      USE m_gf_io2dmat 
      USE m_hdf_tools,ONLY:hdf_init 
      USE m_gf_apws,ONLY:gf_apws_dim 
      USE m_gf_inp_describe 
      USE m_gf_out 
      USE m_gf_embedding 
      USE m_gf_enpara 
      USE m_gf_writetrans
!      use m_gf_curvy2dprojector,only:gf_curvy2dinit                    
      USE m_gf_flippot 
      USE m_gf_kptsp1x1 
      USE m_fleur_INTERFACE,ONLY:fleur_renpara,fleur_mksphhar           &
     &     ,fleur_prp_qfft,fleur_prp_xcfft
      USE m_gf_iotmat,ONLY:init_tmat_storage 
      USE  m_gf_mpi_groups 
                                                                        
      !>                                                                
      IMPLICIT NONE 
      !<--Arguments                                                     
                                                                        
      INTEGER,INTENT(OUT)        ::jspins 
      TYPE(t_atoms),ALLOCATABLE,INTENT(OUT)  ::atoms(:) 
      TYPE(t_mpi),INTENT(OUT)    ::mpi 
      TYPE(t_gfinp),INTENT(OUT)  ::gfinp 
      TYPE(t_cell),ALLOCATABLE,INTENT(OUT)   ::cell(:) 
      TYPE(t_sym),INTENT(OUT)    ::sym 
      TYPE(t_xcpot),INTENT(OUT),ALLOCATABLE  ::xcpot(:) 
      TYPE(t_mix),INTENT(OUT)    ::mix 
      TYPE(t_lapw),ALLOCATABLE,INTENT(OUT)   ::lapw(:) 
      TYPE(t_soc),INTENT(OUT)    ::soc 
      TYPE(t_noco),ALLOCATABLE,INTENT(OUT)   ::noco(:) 
      TYPE(t_fermi),INTENT(OUT)  ::fermi 
      TYPE(t_stars),ALLOCATABLE,INTENT(OUT)  ::stars(:) 
      TYPE(t_sphhar),ALLOCATABLE,INTENT(OUT) ::sphhar(:) 
      TYPE(t_kpts),INTENT(OUT)   ::kpts 
      TYPE(t_enpara),ALLOCATABLE,INTENT(OUT) ::enpara(:) 
      TYPE(t_layers),INTENT(OUT) ::layers 
      CHARACTER(LEN=10) filename 
      INTEGER layerind 
                                                                        
      !>                                                                
      !<--Locals                                                        
                                                                        
      LOGICAL:: lexist,l_tmatmemory 
#ifdef CPP_MPI                                                          
      INCLUDE 'mpif.h' 
      INTEGER ierr,irank 
      CALL MPI_COMM_RANK (MPI_COMM_WORLD,irank,ierr) 
      mpi%pe0=irank==0 
#else                                                                   
      mpi%pe0=.TRUE. 
#endif                                                                  
      CALL hdf_init() 
#ifdef CPP_MPI                                                          
      CALL MPI_BARRIER(mpi_comm_world,ierr) 
#endif                                                                  
                                                                        
      !>                                                                
                                                                        
      !Start the timing                                                 
      CPP_juDFT_timestart("Setup") 
      !Read most of atoms,cell and sym from potential file              
      CALL gf_out_newheader('SETUP of the Calculation') 
      CPP_juDFT_timestart("reading setup")
      CALL gf_read_layerinfo(layers) 
      ALLOCATE( atoms(layers%num_layers)) 
      ALLOCATE( cell(layers%num_layers)) 
      ALLOCATE( lapw(layers%num_layers)) 
      ALLOCATE( stars(layers%num_layers)) 
      ALLOCATE( sphhar(layers%num_layers)) 
      ALLOCATE( enpara(layers%num_layers)) 
      ALLOCATE( noco(layers%num_layers)) 
      ALLOCATE( xcpot(layers%num_layers)) 
      ALLOCATE( gfinp%napw(layers%num_layers)) 
      DO layerind=1,layers%num_layers 
         IF(layerind>1)THEN 
           DEALLOCATE(sym%mrot) 
           DEALLOCATE(sym%invtab) 
           DEALLOCATE(sym%tau) 
          ENDIF 
         CALL gf_readatt(layerind,jspins,                               &
     &     atoms(layerind),cell(layerind),sym)
          cell(layerind)%c=layers%c(layerind)
      ENDDO 
      !Read the gf_inp-file                                             
      CALL gf_rgfinp(layers,gfinp,atoms,xcpot,mix,lapw,soc,fermi,       &
     &                stars,noco)
      if (jspins==1.and.noco(1)%l_noco) then
           write(*,*) "WARNING, in a noco-calculation jspin is set to 2"
           jspins=2
      endif
      gfinp%l_nogno = (.NOT.gfinp%l_fullgreen .AND.(gfinp%l_dos.or.gfinp%l_charge))
                                                                        
      IF (mpi%pe0) CALL gf_inp_describe(                                &
     &                  layers,jspins,atoms,cell,gfinp)                 
!      IF (gfinp%napw == 0 ) THEN   ! region 1 analytically             
!         gfinp%c1 = gfinp%d/2.0                                        
!         gfinp%c2 = gfinp%d/2.0                                        
!      ENDIF                                                            
      ! read the non-collinear input                                    
!      CALL gf_rnocoinp(atoms,noco)                                     
      !Read the kpts                                                    
      CALL gf_out_newheader('INPUT of kpts&enpara') 
      CALL gf_loadkpts(cell(1),sym,mpi%pe0,kpts) 
      !Generate a new kpts-set for p(1x1) unit cell if requested        
!      IF (gfinp%kpts /="none") THEN                                    
!         CALL gf_kptsp1x1(cell,sym,kpts,gfinp)                         
!      ENDIF                                                            
      !Load enparas                                                     
      CALL gf_loadenpara(jspins,layers,atoms,enpara)

! artificial setup                                                      
#ifdef CPP_DEBUG_WITH_STEP                                              
!      INQUIRE(FILE='step_para',EXIST=lexist)                           
!      IF (lexist) THEN                                                 
!      OPEN(123,FILE='step_para',STATUS='old')                          
!      READ(123,*) cell%z1                                              
!      READ(123,*) gfinp%d                                              
!      CLOSE(123)                                                       
!      ENDIF                                                            
#endif                                                                  
                                                                        
                                                                        
                                                                        
      !read energies                                                    
      CALL gf_out_newheader('Energy Grid for GFleur') 
      CALL gf_energies_init(mpi%pe0) 
      !Setup MPI from gf_mpi_groups                                     
      CALL gf_setup_mpi_groups(mpi,layers,kpts,gf_noen()) 
      CPP_juDFT_timestop("reading setup")
      !Generate the stars                                               
      CALL gf_out_newheader('INFORMATION on stars') 
      IF (mpi%pe0) WRITE(6,'("2-D stars",/)') 
      stars(:)%gmax=stars(:)%gmax_inp 
      CPP_juDFT_timestart("Generation of stars") 
      DO layerind=1,layers%num_layers 
         CALL gf_stargen(sym,cell(layerind)%bmat,                       &
     &               stars(layerind)%gmax,stars(layerind))              
      ENDDO 
      CPP_juDFT_timestop("Generation of stars") 
                                                                        
      !Generate the step functions for potential                        
      CPP_juDFT_timestart("gf_step_pot") 
      DO layerind=1,layers%num_layers 
!        IF (mpi%pe0) THEN                                              
!         CALL gf_step_pot(layerind,stars(layerind))                    
!        ENDIF                                                          
      ENDDO 
      CPP_juDFT_timestop("gf_step_pot") 
#ifdef CPP_MPI                                                          
      CALL MPI_BARRIER(mpi_comm_world,ierr) 
#endif                                                                  
      CPP_juDFT_timestart("fft-box")
      !Generate the fft-box for the GGA                                 
      DO layerind=1,layers%num_layers 
         CALL fleur_prp_xcfft(stars(layerind),lapw(layerind),           &
     &                       cell(layerind),xcpot(layerind))            
      ENDDO 
      !Generate the spherical Harmonics                                 
      CPP_juDFT_timestop("fft-box")
      CALL gf_out_newheader('INFORMATION on sphercial harmonics') 
      IF(gfinp%l_gmat)THEN 
        CPP_juDFT_timestart("generate sphhar")

       DO layerind=1,layers%num_layers 
         CALL fleur_mksphhar(                                           &
     &         atoms(layerind),cell(layerind),sym,sphhar(layerind))     
       ENDDO
       CPP_juDFT_timestop("generate sphhar")
      ENDIF 
                                                                        
      !Generate lapw-dimensions                                         
      lapw(1:layers%num_layers)%rkmax=                                  &
     &        lapw(1:layers%num_layers)%rkmax_inp                       
      DO layerind=1,layers%num_layers 
         CALL gf_apws_dim(jspins,                                       &
     &     noco(layerind)%qss,                                          &
     &     kpts,gfinp,cell(layerind),                                   &
     &     (/MINVAL(atoms(layerind)%rmt*atoms(layerind)%lmax)           &
     &     ,MAXVAL(atoms(layerind)%rmt*atoms(layerind)%lmax)/)          &
     &     ,lapw(layerind),layerind)                                    
         lapw(layerind)%nbasfcn = lapw(layerind)%nvd +                  &
     &                     atoms(layerind)%nlotot                       
      IF (noco(layerind)%l_noco)  lapw(layerind)%nbasfcn =              &
     &               2*lapw(layerind)%nbasfcn                           
      ENDDO 
                                                                        
      DO layerind=1,layers%num_layers 
        !Generate FFT-box for charge density                            
        CALL fleur_prp_qfft(stars(layerind),noco(layerind),             &
     &                   lapw(layerind),cell(layerind))                 
      ENDDO 
                                                                        
      !<--Setup of stepfunction                                         
                                                                        
      CPP_juDFT_timestart("Setup of step-functions") 
      DO layerind=1,layers%num_layers 
         CALL gf_initstepsanaly(stars(layerind),gfinp%napw(layerind)) 
         IF (mpi%pe0) THEN 
            CALL gf_calcstepsanaly(layerind,            &
     &           cell(layerind),                                        &
     &           stars(layerind),gfinp%napw(layerind))                  
         ENDIF 
      ENDDO 
      CPP_juDFT_timestop("Setup of step-functions") 
                                                                        
!      CPP_juDFT_timestart("gf_curvy2dinit")                              
!      if(gfinp%l_nohelpregion)then                                     
!       do layerind=1,layers%num_layers                                 
!          IF (mpi%pe0) CALL gf_curvy2dinit(layerind                    
!     $         ,gfinp%napw(layerind),stars(layerind),cell(layerind),   
!     $         .TRUE.)                                                 
!       enddo                                                           
!      endif                                                            
!      CPP_juDFT_timestop("gf_curvy2dinit")                               
                                                                        
      !>                                                                
      !<-- Flip the spins if "spinflip" file exists                     
                                                                        
      INQUIRE(FILE ='gf_spinflip',EXIST = lexist) 
      IF (lexist) THEN 
          CALL juDFT_error("think about it",calledby="gf_setup.F90")
!         CALL gf_flipPot(jspins,atoms,stars,sphhar,enpara,gfinp        
!     +        ,noco%l_noco)                                            
      ENDIF 
                                                                        
      !>                                                                
      !OPEN 2D-file
      CPP_juDFT_timestart("init io2dmat")

      CALL init(kpts,gfinp,lapw(1),sym,layers,jspins,                   &
     &           noco(1)%l_noco                                         &
     &     ,mpi,l_tmatmemory)
      CPP_juDFT_timestop("init io2dmat")

                                !this is from gf_io2dmat                

      IF (l_tmatmemory) THEN 
         CALL init_tmat_storage(noco(1)%l_noco,gf_noen()                &
     &        ,layers%num_layers,SIZE(lapw(1)%global2dList,2))          
      ENDIF 
! check the embedding potential                                         
      IF ( gfinp%l_addemb ) THEN 
          CALL gf_checkemb(layers,cell(1)%amat)
      ENDIF 
      !write an xyz-file for plotting                                   
!      CALL gf_inp2plot(atoms,cell,gfinp)                               
                                                                        
!if scattering states for any energy are calculated open eig-file       
!      IF (gfinp%l_charge.and.any(direction.ne.0)) THEN                 
!         CALL openeig(lapw%nbasfcn,lapw%nv2d*gf_noEn(),kpts%nkpts,jspin
!     +        ,MAXVAL(atoms%lmax),atoms%nlod,atoms%ntype,create=.TRUE. 
!     +        ,readonly=.FALSE.)                                       
!      ENDIF                                                            
                                                                        
      !init debugging-module                                            
#ifdef CPP_DEBUG                                                        
      !CALL gf_debug_init(sphhar,atoms,cell,stars,sym,mix)
#endif                                                                  
                                                                        
      !do all kinds of sanity checks on input                           
      CALL priv_setup_test(gfinp,atoms(1),mpi,mix,noco(1)) 

      IF (gfinp%l_hdfio) THEN
         !open file for transmission data if needed
         if (gfinp%curr.ne.0) call gf_writetrans_hdfopen(mpi,kpts,sym,layers%num_layers,jspins)
      ENDIF
                                                                        
      CALL gf_out_newheader('GF-SETUP done') 
      CPP_juDFT_timestop("Setup") 
                                                                        
      END SUBROUTINE 
                                                                        
      !>                                                                
      ! all further subroutines are private                             
                                                                        
      !<-- S: priv_setup_test(gfinp,mix,lapw,kpts,enpara)               
      SUBROUTINE priv_setup_test(gfinp,atoms,mpi,mix,noco) 
!-----------------------------------------------                        
!   Test the input for correct combination of switches                  
!   More tests should be added here!!!                                  
!           (last modified: 2004-00-00) D. Wortmann                     
!-----------------------------------------------                        
      USE m_gf_types 
      use m_juDFT 
      IMPLICIT NONE 
      !<--Arguments                                                     
                                                                        
      TYPE(t_gfinp),INTENT(IN)  :: gfinp 
      TYPE(t_atoms),INTENT(IN)  :: atoms 
      TYPE(t_mpi),INTENT(IN)    :: mpi 
      TYPE(t_mix),INTENT(INOUT) :: mix 
      TYPE(t_noco),INTENT(IN)   :: noco 
                                                                        
      !>                                                                
      !<-- Locals                                                       
      LOGICAL             :: test 
      LOGICAL             :: ok 
      ok = .TRUE. 
      !>                                                                
      !<-- Things that only should be warned about                      
                                                                        
      IF (.NOT.gfinp%l_charge.AND.mix%iter>0) THEN 
         WRITE(6,*)                                                     &
     &      "WARNING, using iter>0 in non-selfconsistent calculation"   
      ENDIF 
                                                                        
      !>                                                                
                                                                        
      !<-- Tests that should disappear somtime :-)                      
                                                                        
      IF (gfinp%npw /= 0)                                               &
     &      CALL juDFT_error("Region II can only be treated analytically",calledby="gf_setup.F90")
                                                                        
      !>                                                                
      !<-- Test for CBS-mode                                            
                                                                        
      IF (gfinp%l_CBS) THEN 
         mix%iter=1 
         IF (.NOT.gfinp%l_gmat.AND..NOT.gfinp%l_tmat) ok = Error(       &
     &        "Check the tmat,gmat-switches for CBS-mode")              
         IF (gfinp%l_charge) ok = Error                                 &
     &        ("charge-generation cannot be used in CBS-MODE")          
         IF (gfinp%l_addemb) ok = Error                                 &
     &        ("l_addemb cannot be used in CBS-mode")                   
!         INQUIRE(file ="gf_cdn.hdf",exist = test)                      
!         IF (test) ok = Error("Found gf_cdn.hdf-file in CBS-mode")     
!         IF (MINVAL(atoms%lmax)<11) ok = Error("Too small lmax for CBS"
      ENDIF 
                                                                        
      !>                                                                
      !<-- Test for sc                                                  
      IF (gfinp%l_charge) THEN 
         IF (.NOT.gfinp%l_addemb) ok       = Error(                     &
     &        "Addemb must be specified for generation of charge")      
         !IF (MAXVAL(atoms%lmax)>10) ok = Error(                         &
     !&        "Too large lmax for charge generation")
      ENDIF 
      !>
      !Test for current
      IF (gfinp%curr .ne. 0) THEN
          IF (.NOT.gfinp%l_addemb) ok       = Error(                     &
     &        "Addemb must be specified for calculation of conductances")
      ENDIF
      !<-- unsorted misc tests                                          
                                                                        
      IF (gfinp%napw(1) == 0.AND.noco%l_noco) ok = error(               &
     &     "No noco for a pure vacuum calculation")                     
      IF (gfinp%l_IEC.AND.gfinp%l_addemb) ok = error(                   &
     &     "in IEC mode you must set addemb = F!")
                                                                        
      !>                                                                
                                                                        
      IF (.NOT.ok) CALL juDFT_error("GF-setup:setup_test") 
                                                                        
                                                                        
      CONTAINS 
      !Error Message output!                                            
      FUNCTION ERROR(string) 
      CHARACTER(LEN=*),INTENT(IN) :: String 
      LOGICAL              :: ERROR 
      WRITE(6,*) "Setup Error in gf_setup detected" 
      WRITE(6,*) String 
      WRITE(*,*) "Setup Error in gf_setup detected" 
      WRITE(*,*) String 
                      !indicates presence of Error :-), can changed to  
      ERROR = .FALSE. 
                      !ignore all tests                                 
      END FUNCTION 
      END SUBROUTINE 
      !>                                                                
      END                                           
