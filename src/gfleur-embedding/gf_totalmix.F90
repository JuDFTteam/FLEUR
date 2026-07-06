!--------------------------------------------------------------------------------
! Copyright (c) 2026 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions 
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------
      MODULE m_gf_totalmix 
      IMPLICIT NONE
!-------------------------------------------------------------          
!Module provides an interface to the mixing part of FLEUR               
!as used in the GF-code                                                 
!------------------------------------------------------------           
      PRIVATE 
      PUBLIC gf_totalmix 
      CONTAINS 
                                                                        
      !<-- S:gf_totalmix(qtot_nuc,qtot_el)                              
      SUBROUTINE gf_totalmix(jspins,layers,atoms,sphhar,stars,mix,gfinp,mpi)
!-----------------------------------------------                        
!                                                                       
!           (last modified:09-11-30) D. Wortmann                        
!-----------------------------------------------                        
      USE m_gf_iodop 
      USE m_gf_types 
      USE m_gf_math 
      use m_juDFT 
      IMPLICIT NONE 
      !<-- Arguments                                                    
      INTEGER,INTENT(IN)        :: jspins 
      TYPE(t_layers),INTENT(IN) :: layers 
      TYPE(t_atoms),INTENT(IN)  :: atoms(:) 
      TYPE(t_sphhar),INTENT(IN) :: sphhar(:) 
      TYPE(t_stars),INTENT(IN)  :: stars(:) 
      TYPE(t_gfmix),INTENT(IN)    :: mix 
      TYPE(t_embinp),INTENT(IN)  :: gfinp 
      TYPE(t_gfmpi),intent(in)    :: mpi
      !>                                                                
      !<-- Locals                                                       
      INTEGER,PARAMETER   :: maxiter = 100 
      LOGICAL             :: input 
      INTEGER             :: layer,iter,l 
      REAL                :: charge 
      REAL                :: in(layers%num_layers,maxiter) 
      REAL                :: out(layers%num_layers,maxiter) 
      REAL                :: new(layers%num_layers) 
      REAL   ,ALLOCATABLE :: vr(:,:,:,:), vr_diff(:,:,:,:) 
      COMPLEX,ALLOCATABLE :: vpw(:,:),vpw_diff(:,:) 
      !>                                                                
                                                                        
      !<-- read layerresolved charges                                   
      WRITE(*,*) "In gf_totalmix" 
      OPEN(99, FILE = "gf_totalmix") 
      layer = 0;input = .TRUE. 
      iter = 1 
      DO 
         READ(99,"(i5,f30.20)",END = 100) l,charge 
         IF (l<layer) THEN 
            input=.NOT.input 
            IF (input) THEN 
               iter = iter+1 
               IF (iter>maxiter) CALL juDFT_error                         &
     &              ("gf_totalmix more than maxiter iterations")        
            ENDIF 
         ENDIF 
         layer = l 
         IF (input) THEN 
            in(layer,iter) = charge 
         ELSE 
            out(layer,iter) = charge 
         ENDIF 
      ENDDO 
  100 CLOSE(99) 
      !>                                                                
      !<-- Do a mixing                                                  
      IF (iter>1) THEN 
!         new = simple_broyden(in(:,:iter),in(:,:iter)-out(:,:iter))    
         new  = (1-mix%alpha)*in(:,iter)+mix%alpha*out(:,iter) 
      ELSE 
         new = (1-mix%alpha)*in(:,1)+mix%alpha*out(:,1) 
      ENDIF 
      !>                                                                
      !<-- construct new charges                                        
      DO layer = 1,layers%num_layers 
         ALLOCATE( vr(MAXVAL(Atoms(layer)%jri),0:MAXVAL(Sphhar(layer    &
     &        )%nlh),Atoms(layer)%ntype,Jspins))                        
         ALLOCATE( vpw(stars(layer)%ng3,jspins) ) 
         ALLOCATE( vr_diff(MAXVAL(Atoms(layer)%jri)                     &
     &        ,0:MAXVAL(Sphhar(layer)%nlh),Atoms(layer)%ntype,Jspins))  
         ALLOCATE( vpw_diff(stars(layer)%ng3,jspins) ) 
                                                                        
         CALL gf_loddop(GF_CDNSTARTFILE,layer,jspins,                   &
     &        atoms(layer),stars(layer),sphhar(layer),                  &
     &        vr,vpw)
         CALL gf_loddop(GF_CDNDIFFFILE,layer,jspins,                    &
     &        atoms(layer),stars(layer),sphhar(layer),                  &
     &        vr_diff,vpw_diff)
         vr = vr+(new(layer)-in(layer,1))*vr_diff 
         vpw = vpw+(new(layer)-in(layer,1))*vpw_diff 
         CALL gf_renamepot(GF_CDNFILE,mpi%iodop_subcom,layer)
         CALL gf_wrtdop(GF_CDNFILE,layer,jspins,                        &
     &        gfinp,atoms(layer),stars(layer),sphhar(layer),            &
     &        vr,vpw,.FALSE.,1)                                         
         DEALLOCATE(vpw,vr,vr_diff,vpw_diff) 
      ENDDO 
      !>                                                                
      END SUBROUTINE 
      !>                                                                
                                                                        
                                                                        
                                                                        
                                                                        
                                                                        
      END                                           
