!--------------------------------------------------------------------------------
! Copyright (c) 2026 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions 
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------
      MODULE m_gf_gfleur_basic 
      use m_juDFT
      IMPLICIT NONE
      CONTAINS 
      SUBROUTINE gf_gfleur_basic(                                       &
     &     jspin,jspins,layer,layers,nk,                                &
     &     gfinp,atoms,sphhar,stars,sym,cell,mpi,soc,noco,              &
     &     bk,lapw,enpara,vpw,vr,vs_mmp,tlmplm_DATA,raddata,            &
     &     pot_aux)                                                     
!*************************************************                      
!     Set up of the Hamiltonian (always).                               
!     In the spectral mode: Diagonalization of the                      
!     Hamiltonian.                                                      
!     If the number of layers is larger than one                        
!     the transfer-matrix of the layer is computed                      
!     in addition.                                                      
!     Frank Freimuth, November 2007                                     
!*************************************************                      
      USE m_gf_types 
      USE m_gf_hs1lapw 
      USE m_gf_spectrum 
      USE m_gf_get_spectrum 
#include "juDFT_env.h" 
      USE m_gf_gfcn 
      USE m_gf_tmat, ONLY  : gf_tmat_layers 
      use m_juDFT 
      USE m_gf_energies,ONLY:gf_noen 
      USE m_gf_io2dmat
      USE m_gf_get_spectrum 
                                                                        
      IMPLICIT NONE 
      INTEGER,INTENT(IN)          :: jspin,jspins 
      INTEGER,INTENT(IN)          :: layer 
      TYPE(t_layers),INTENT(IN)   :: layers 
      INTEGER,INTENT(IN)          :: nk 
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
      REAL, INTENT(IN)            :: Bk(:,:) 
      TYPE(t_lapw),INTENT(INOUT)  :: lapw 
      TYPE(t_enpara),INTENT(INOUT) :: enpara 
      COMPLEX,INTENT(INOUT)       :: vpw(:,:) 
      REAL, INTENT(INOUT)         :: vr(:,0:,:,:) 
      COMPLEX,INTENT(IN)          :: vs_mmp(-3:,-3:,:,:) 
      TYPE(t_tlmplm),INTENT(INOUT)   :: tlmplm_DATA 
      TYPE(t_raddata),INTENT(INOUT)  :: raddata 
      COMPLEX,INTENT(INOUT)          :: pot_aux(2,2) 
                                                                        
                                      !not used in this subroutine      
      COMPLEX, ALLOCATABLE :: g(:,:) 
      COMPLEX, ALLOCATABLE :: gij(:,:,:),gij2(:,:,:),unity(:,:)
      INTEGER en 
      LOGICAL l_writespectrum 
      CPP_juDFT_timestart("gf_hs1lapw")
      CALL gf_hs1lapw(                                                  &
     &     jspin,jspins,layer,layers,nk,                                &
     &     gfinp,atoms,sphhar,stars,sym,cell,mpi,soc,noco,              &
     &     bk(1:3,nk),lapw,enpara,vpw,vr,vs_mmp,tlmplm_DATA,raddata)
      CPP_juDFT_timestop("gf_hs1lapw")
      IF(gfinp%l_spectral)THEN 
         l_writespectrum=.FALSE. 
         !IF(gfinp%l_charge)l_writespectrum=.TRUE.
         CPP_juDFT_timestart("gf_get_spectrum") 
         CALL gf_spectrum_clean() 
         CALL gf_get_spectrum(layer,jspin,gfinp,cell,lapw,.FALSE.,      &
     &        gfinp%l_fullgreen,gfinp%l_nogno,gfinp%l_nohelpregion,     &
     &        l_writespectrum,nk,.false.)
         CPP_juDFT_timestop("gf_get_spectrum") 
      ENDIF 

      !calculate the transfer matrix
      IF(layers%num_layers>1)THEN 
         ALLOCATE( gij(2*lapw%nv2_tot,2*lapw%nv2_tot,1) ) 
         !ALLOCATE( gij2(2*lapw%nv2_tot,2*lapw%nv2_tot,1) )
         !ALLOCATE( unity(lapw%nv2_tot,lapw%nv2_tot))
         !UNITY=0.0
         !DO en=1,lapw%nv2_tot
         !  unity(en,en)=1.0
         !ENDDO
         !energy loop
         DO en = 1,gf_noen() 
         IF (gfinp%napw(layer)>-1) THEN 
         !<-- Standard Calculation of retarded GF                       
                                                                        
            CPP_juDFT_timestart("gf_gfcn") 
            CALL gf_GFCN(.TRUE.,layer,en,nk,jspin,cell,lapw,gfinp       &
     &           ,noco%l_noco,sym%invs,.false.,g,gij,real_energy=gfinp%l_CBS.or.gfinp%curr.ne.0)

            !CALL gf_GFCN(.TRUE.,layer,en,nk,jspin,cell,lapw,gfinp       &
            !     ,noco%l_noco,sym%invs,.false.,g,gij2,real_energy=.true.)

            !gij(1:lapw%nv2_tot,1:lapw%nv2_tot,1)=gij(1:lapw%nv2_tot,1:lapw%nv2_tot,1)
            !gij(1:lapw%nv2_tot,lapw%nv2_tot+1:2*lapw%nv2_tot,1)=gij2(1:lapw%nv2_tot,lapw%nv2_tot+1:2*lapw%nv2_tot,1)
            !gij(lapw%nv2_tot+1:2*lapw%nv2_tot,lapw%nv2_tot+1:2*lapw%nv2_tot,1)=gij(lapw%nv2_tot+1:2*lapw%nv2_tot,lapw%nv2_tot+1:2*lapw%nv2_tot,1)
            !gij(lapw%nv2_tot+1:2*lapw%nv2_tot,1:lapw%nv2_tot,1)=gij2(lapw%nv2_tot+1:2*lapw%nv2_tot,1:lapw%nv2_tot,1)

            CPP_juDFT_timestop("gf_gfcn") 
         ELSE 
         !<-- Vacuum calculation
            CALL juDFT_error("BUG:vacuum in gfleur_basic")
            !CALL gf_vacuum(gfinp,lapw                                   &
     !&                       ,en,nk,jspin)
                        !>                                              
         ENDIF 
         CPP_juDFT_timestart("gf_tmat") 
         IF (pot_aux(1,jspin) /= pot_aux(2,jspin)) CALL                 &
     &      juDFT_error                                                   &
     &      ("Cannot treat different aux_potentials")                   
                                                                        
         CALL gf_TMAT_layers(                                           &
     &           layer,layers,en,nk,jspin,                              &
     &           mpi,lapw,gfinp,pot_aux(1,jspin),gij(:,:,1))            
         CPP_juDFT_timestop("gf_tmat") 
             !energy loop                                               
       ENDDO 
       DEALLOCATE( gij)
       !DEALLOCATE(gij2,unity)
            !calculate transfer matrix                                  
      ENDIF 
      END SUBROUTINE gf_gfleur_basic 
      END                                           
