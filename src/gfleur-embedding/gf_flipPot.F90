!--------------------------------------------------------------------------------
! Copyright (c) 2026 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions 
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------
      MODULE m_gf_flipPot 
      use m_juDFT
      use m_juDFT 
      IMPLICIT NONE
!-----------------------------------------------                        
! DESC: Small subroutine that read the potential/charge                 
! for embedding calculations and flips spin for                         
!    z>0                                                                
!                 Daniel Wortmann, (04-08-17)                           
!-----------------------------------------------                        
      CONTAINS 
      !<-- S: gf_flipPot(jspins,atoms,stars,sphhar,enpara,l_noco)       
      SUBROUTINE gf_flipPot(jspins,atoms,stars,sphhar,enpara,gfinp      &
     &     ,noco,mpi)
!-----------------------------------------------                        
!           (last modified: 04-08-17) D. Wortmann                       
!-----------------------------------------------                        
      USE m_gf_types 
      USE m_gf_iodop 
!      USE m_fleur_enpara                                               
      IMPLICIT NONE 
      !<--Arguments                                                     
                                                                        
      INTEGER                          layer 
      INTEGER,INTENT(IN)            :: jspins 
      TYPE(t_atoms),INTENT(IN)      :: atoms 
      TYPE(t_stars),INTENT(IN)      :: stars 
      TYPE(t_sphhar),INTENT(IN)     :: sphhar 
      TYPE(t_enpara),INTENT(INOUT)  :: enpara 
      TYPE(t_embinp),INTENT(IN)      :: gfinp 
      TYPE(t_gfmpi),INTENT(IN)        :: mpi 
      TYPE(t_noco),INTENT(IN)       :: noco
                                                                        
      !>                                                                
      !<-- Locals                                                       
      INTEGER           :: n,na 
      LOGICAL           :: l_exist 
      COMPLEX           :: vpw(stars%ng3,jspins) 
      REAL              :: vr(MAXVAL(Atoms%jri),0:MAXVAL(Sphhar%nlh),   &
     &     Atoms%ntype,Jspins)                                          
      REAL,ALLOCATABLE  :: el_tmp(:),ello_tmp(:) 
      !>                                                                
                                                                        
      IF(jspins/= 2) CALL juDFT_error("No Spinflip for jspins = 1")
      IF(noco%l_noco) CALL juDFT_error("No Spinflip for noco")
                                                                        
                                                                        
      !<-- the potential                                                
                                                                        
      CALL gf_loddop(gf_potfile,layer,jspins,atoms,stars,sphhar,        &
     &     vr,vpw,noco,.FALSE.)
      CALL local_flip(atoms,stars,vpw,vr) 
      CALL gf_wrtdop(GF_POTFILE,layer,jspins,                           &
     &     gfinp,atoms,stars,sphhar,                                    &
     &     vr,vpw,noco%l_noco,mpi%self_subcom)
                                                                        
      !>                                                                
      !<-- The charge                                                   
                                                                        
      INQUIRE(FILE ="gf_cdn.hdf",EXIST = l_exist) 
      IF (l_exist) THEN 
         CALL gf_loddop(gf_cdnFILE,jspins,layer,atoms,stars,sphhar,     &
     &        vr,vpw,noco,.FALSE.)
         CALL local_flip(atoms,stars,vpw,vr) 
         CALL gf_wrtdop(GF_CDNFILE,jspins,layer,                        &
     &        gfinp,atoms,stars,sphhar,                                 &
     &        vr,vpw,noco%l_noco,mpi%self_subcom)
      ENDIF 
                                                                        
      !>                                                                
      !<-- The energy -parameters....                                   

      ALLOCATE(el_tmp(SIZE(enpara%el0,1)),ello_tmp(SIZE(enpara%ello0,1))) 
      na = 1
      DO n = 1,atoms%ntype 
         !flip up/down spin if positive position                        
         IF (atoms%pos(3,na)>0) THEN 
            el_tmp = enpara%el0(:,n,1) 
            enpara%el0(:,n,1) = enpara%el0(:,n,2) 
            enpara%el0(:,n,2) = el_tmp 
            ello_tmp = enpara%ello0(:,n,1) 
            enpara%ello0(:,n,1) = enpara%ello0(:,n,2) 
            enpara%ello0(:,n,2) = ello_tmp 
         ENDIF 
         !if at zero, average spin!                                     
         IF (atoms%pos(3,na)>0) THEN 
            enpara%el0(:,n,1) = (enpara%el0(:,n,2)+enpara%el0(:,n,1))/2 
            enpara%el0(:,n,2) = enpara%el0(:,n,1) 
            enpara%ello0(:,n,1) = (enpara%ello0(:,n,2)+enpara%ello0(:,n,1))&
     &           /2                                                     
            enpara%ello0(:,n,2) = enpara%ello0(:,n,1) 
         ENDIF 
         na = na+atoms%neq(n) 
      ENDDO 
      DEALLOCATE(el_tmp,ello_tmp) 
!      CALL fleur_writeenpara(atoms,jspins,enpara)                      
       CALL juDFT_error("writing of enpara no longer implemented",calledby="gf_flipPot.F90")

      !>                                                                
      CALL juDFT_end("Potential was spin-flipped for z>0") 
                                                                        
      CONTAINS 
      !<-- S: local_flip(atoms,stars,vpw,vr)                            
      SUBROUTINE local_flip(atoms,stars,vpw,vr) 
!-----------------------------------------------                        
!  performs the actual flipping.....                                    
!           (last modified: 2004-00-00) D. Wortmann                     
!-----------------------------------------------                        
      USE m_gf_types 
      USE m_gf_rsmesh 
                                                                        
      IMPLICIT NONE 
      !<--Arguments                                                     
      TYPE(t_atoms),INTENT(IN) :: atoms 
      TYPE(t_stars),INTENT(IN) :: stars 
      COMPLEX,INTENT(INOUT)      :: vpw(:,:) 
      REAL   ,INTENT(INOUT)      :: vr(:,:,:,:) 
      !>                                                                
      !<-- Locals                                                       
      INTEGER             :: na,n 
      TYPE(t_rsdata)      :: up,down 
      REAL,ALLOCATABLE    ::vr_tmp(:,:),tmpxy(:,:) 
                                                                        
      !>                                                                
                                                                        
      !<-- convert MT-part                                              
      ALLOCATE(vr_tmp(SIZE(vr,1),SIZE(vr,2))) 
      na = 1 
      DO n = 1,atoms%ntype 
         !flip up/down spin if positive position                        
         IF (atoms%pos(3,na)>0) THEN 
            vr_tmp = vr(:,:,n,1) 
            vr(:,:,n,1) = vr(:,:,n,2) 
            vr(:,:,n,2) = vr_tmp 
         ENDIF 
         !if at zero, average spin!                                     
         IF (atoms%pos(3,na)>0) THEN 
            vr(:,:,n,1) = (vr(:,:,n,2)+vr(:,:,n,1))/2 
            vr(:,:,n,2) = vr(:,:,n,1) 
         ENDIF 
         na = na+atoms%neq(n) 
      ENDDO 
      !>                                                                
      !<-- convert INT-part                                             
                                                                        
      CALL gf_rsforstars(stars,up) 
      CALL gf_rsforstars(stars,down) 
      CALL gf_pwtors(vpw(:,1),stars,up) 
      CALL gf_pwtors(vpw(:,2),stars,down) 
      !now flip data                                                    
      ALLOCATE(tmpxy(up%grid(1),up%grid(2))) 
      DO n = 1,up%grid(3)/2 
         tmpxy = up%data(:,:,n) 
         up%data(:,:,n) = down%data(:,:,n) 
         down%data(:,:,n) = tmpxy 
      ENDDO 
      CALL gf_rstopw(up,stars,vpw(:,1)) 
      CALL gf_rstopw(down,stars,vpw(:,2)) 
                                                                        
      !>                                                                
      DEALLOCATE(up%data,down%data,tmpxy,vr_tmp) 
                                                                        
                                                                        
      END SUBROUTINE 
      !>                                                                
                                                                        
      END SUBROUTINE 
      !>                                                                
      END                                           
