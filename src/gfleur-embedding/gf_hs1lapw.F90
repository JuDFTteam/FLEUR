!--------------------------------------------------------------------------------
! Copyright (c) 2026 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions 
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------
      MODULE m_gf_hs1lapw 
      use m_juDFT 
      IMPLICIT NONE
      PRIVATE 

      PUBLIC gf_hs1lapw 
      CONTAINS 
      !<-- S: gf_hs1lapw()                                              
                                                                        
      SUBROUTINE gf_hs1lapw(                                            &
     &     jspin,jspins,layer,layers,nk,                                &
     &     gfinp,atoms,sphhar,stars,sym,cell,mpi,soc,noco,              &
     &     bk,lapw,enpara,vpw,vr,vs_mmp,tlmplm_DATA,raddata)            
!*****************************************************************      
!     DESC:This subroutine generates the FLAPW-Hamiltonian and the      
!     Overlapp                                                          
!     Matrix by calling the appropriate FLEUR-subroutines               
!     Daniel Wortmann, Sat Mar  2 13:46:30 2002                         
!*****************************************************************      

      USE m_gf_hsdata 
      use m_gf_types
      USE gf_hsint,ONLY    : hsintstep
#include "juDFT_env.h" 
                            !interface to fleur                         
      USE m_fleur_hsetup 
      USE m_gf_spectrum,ONLY: gf_spectrum_clean 
      use m_gf_stepsanaly, only: gf_initstepsanaly, gf_stepf_nohelpregion
!      use m_gf_overlap                                                 
      IMPLICIT NONE 
      !<--Arguments                                                     
                                                                        
      INTEGER,INTENT(IN)          :: jspin,jspins 
      INTEGER                     :: layer 
      TYPE(t_layers),INTENT(IN)   :: layers 
      INTEGER,INTENT(IN)          :: nk 
      !     Type                                                        
      TYPE(t_gfinp), INTENT(IN)   :: gfinp 
      TYPE(t_atoms),INTENT(IN)    :: atoms 
      TYPE(t_sphhar),INTENT(IN)   :: sphhar 
      TYPE(t_stars),INTENT(IN)    :: stars 
      TYPE(t_sym),INTENT(IN)      :: sym 
      TYPE(t_cell),INTENT(IN)     :: cell 
      TYPE(t_mpi),INTENT(IN)      :: mpi 
      TYPE(t_soc),INTENT(IN)      :: soc 
      TYPE(t_noco),INTENT(IN)     :: noco 
                                                    ! Kpt               
      REAL, INTENT(IN)            :: Bk(3) 
      TYPE(t_lapw),INTENT(INOUT)  :: lapw 
      TYPE(t_enpara),INTENT(INOUT) :: enpara 
      !     .. Array Arguments ..                                       
                                                                        
      REAL, INTENT(INOUT)         :: vr(:,0:,:,:) 
      COMPLEX,INTENT(INOUT)       :: vpw(:,:) 
      COMPLEX,INTENT(IN)          :: vs_mmp(-3:,-3:,:,:) 
      TYPE(t_tlmplm),INTENT(INOUT)   :: tlmplm_DATA 
      TYPE(t_raddata),INTENT(INOUT)  :: raddata 
                                                                        
      !>                                                                
      !<-- Locals                                                       
                                                                        
      ! for LDA+u                                                       
      INTEGER, PARAMETER :: LMAXB = 3 
      ! dimensions                                                      
      INTEGER            ::matsize,err
                                                                        
      !Large data fields                                                
      COMPLEX,ALLOCATABLE    :: ustep(:,:,:) 
      COMPLEX, POINTER :: h(:), s(:) 
      INTEGER                  :: i 
      write(6,*) "Old nbasfcn:",lapw%nbasfcn
      lapw%nbasfcn=lapw%nv_tot+atoms%nlotot
      if (noco%l_noco) lapw%nbasfcn=lapw%nbasfcn+atoms%nlotot
      write(6,*) "New nbasfcn:",lapw%nbasfcn
!#ifdef CPP_MPI                                                         
!      Matsize = Lapw%nbasfcn*(Lapw%nbasfcn+1+mpi%n_size)/2/mpi%n_size  
!#else                                                                  
      Matsize = Lapw%nbasfcn*(Lapw%nbasfcn+1)/2 
!#endif                                                                 
                                                                        
      !>                                                                
                                                                        
                                                                        
      ALLOCATE(ustep(-stars%mx1:stars%mx1,-stars%mx2:stars%mx2,         &
     &        -2*gfinp%napw(layer):2*gfinp%napw(layer)),stat=err)
      if(err/=0) then
        write(6,*) "Step function:",layer
        write(6,*) stars%mx1,stars%mx2,gfinp%napw(layer)
        call juDFT_error("Could not allocate ustep",calledby="gf_hs1lapw")
      endif
      IF(.NOT.gfinp%l_nohelpregion)THEN 
         CALL juDFT_error("l_nohelpregion required",calledby              &
     &        ="gf_hs1lapw.F90")
      ELSE 
         CALL gf_initstepsanaly(stars,gfinp%napw(layer)) 
         CALL gf_stepf_nohelpregion(layer,stars%mx1       &
     &        ,stars%mx2,gfinp%napw(layer),ustep)                       
      ENDIF 
                                                                        
      !>                                                                
                                                                        
      !>                                                                
      !<-- This has to be done for every kpt and every spin             
      ! Now the core part of the routine follows                        
      ! It sets up the H and S matrix and is allways executed           

      CALL  gf_spectrum_clean() 
      ! Assign pointers to storage                                      
      CALL gf_assignHS(matsize,jspin,nk,layer,h,s) 
      !<-- interstitial contribution to H+S                             
      IF ( gfinp%l_solwil ) THEN 
         !<--soler wiliams                                              
                                                                  !!")  
         CALL juDFT_error("Soler Wiliams Basis not supported anymore !")

         CPP_juDFT_timestart("gf_hsint") 
!         CALL hsintsolwil(Lapw,matsize,Atoms,                          
!     >        jspin,                                                   
!     +        CELL,                                                    
!     &        bk,                                                      
!     +        stars%mx1,Stars%mx2,                                     
!     &        vpw,                                                     
!     &        gfinp,h,s)                                               
         lapw%nmat = lapw%nv(jspin) 
         IF (noco%L_noco) lapw%nmat = lapw%nv(1)+lapw%nv(2) 
         CPP_juDFT_timestop("gf_hsint") 
         !>                                                             
      ELSE 
         !<--standard setup                                             
         CPP_juDFT_timestart("gf_hsint") 
         CALL hsintstep(                                                &
     &        matsize,lapw,stars,                                       &
     &        gfinp%napw(layer),ustep,jspin,mpi,bk,cell%bbmat,          &
     &        atoms%nlotot                                              &
     &        ,vpw,noco%l_noco,h,s)                                     
         CPP_juDFT_timestop("gf_hsint") 
                                                                        
         !>                                                             
      ENDIF 
      !>                                                                
      !<-- MT-contribution to H and S                                   
      i=(lapw%nv(1)*(lapw%nv(1)+1))/2+1
      CPP_juDFT_timestart("hssphn") 
                                                                        
      CALL fleur_hssphn(jspins,jspin,atoms,sphhar,sym,cell              &
     &     ,mpi,soc,noco,lapw,enpara,bk,vr,vs_mmp,                      &
     &     tlmplm_data,raddata,                                         &
     &     h,s)                                                         
      CPP_juDFT_timestop("hssphn") 
                                                                        
                                                                        
      !>                                                                
                                                                        
      !<-- calculate restriced overlapp matrix                          
!      IF (gfinp%l_overlap) CALL gf_init_overlap(jspin,jspins           
!     +     ,gfinp%napw(layer),bk,ustep,lapw,stars,cell,sym,mpi         
!     +     ,raddata,atoms)                                             
                                                                        
      !<-- write all data to files                                      
                                                                        
                                                                        
#ifdef CPP_MPI                                                          
!      Matsize = Lapw%nbasfcn*(Lapw%nbasfcn+1)/2                        
!      CALL gf_collectHS(mpi,lapw%nbasfcn,matsize)                      
#endif                                                                  
      IF ((gfinp%l_charge.or.gfinp%l_dos.OR.gfinp%l_writehs.OR..NOT.gfinp%l_spectral).AND. &
     &     mpi%n_rank == 0) CALL gf_writeHS(gfinp%l_savemem)
                                                                        
      !>                                                                
                                                                        
      !>                                                                
      END SUBROUTINE 
                                                                        
      !>                                                                
                                                                        
      END                                           
