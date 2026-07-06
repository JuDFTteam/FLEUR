!--------------------------------------------------------------------------------
! Copyright (c) 2026 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions 
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------
      MODULE m_gf_makeCoulBoundary 
          IMPLICIT NONE
      CONTAINS 
                                                                        
      !<-- S: gf_makeCoulBoundary(layer,layers,cell,stars)              
      SUBROUTINE gf_makeCoulBoundary(layer,layers,cell,stars) 
!-----------------------------------------------                        
!   Use the Projectors to construct the boundary values                 
!   of the Coulomb potential and save this to embpot file               
!           (last modified: 07-12-12) D. Wortmann                       
!-----------------------------------------------                        
      USE m_gf_types 
      USE  m_gf_potproject,ONLY: gf_makePotProjector 
      USE m_gf_iodop 
      USE hdf5 
      USE m_gf_io2dmat 
      USE m_hdf_tools 
      USE m_gf_ioCoulBoundary 
      IMPLICIT NONE 
      !<--Arguments                                                     
      INTEGER,INTENT(IN)        :: layer 
      TYPE(t_layers),INTENT(IN) :: layers 
      TYPE(t_cell),INTENT(IN)   :: cell 
      TYPE(t_stars),INTENT(IN)  :: stars 
      !>                                                                
      !<-- Locals                                                       
      INTEGER(hid_t)      :: fid 
      COMPLEX,ALLOCATABLE :: curvyproj(:,:,:) 
      COMPLEX,ALLOCATABLE :: vpw(:),qpw(:) 
      COMPLEX,ALLOCATABLE :: v_bound(:) 
      INTEGER             :: status 
      !>                                                                
                                                                        
      ALLOCATE(curvyproj(stars%nq2,stars%nq3,2),STAT = status) 
      IF (status /= 0) THEN 
         WRITE(*,*) "WARNING, could not generate boundary potential" 
         RETURN 
      ENDIF 
      ALLOCATE(vpw(stars%nq3),qpw(stars%nq3)) 
      ALLOCATE(v_bound(stars%nq2)) 
      CALL gf_lodcoul(GF_POTFILE,layer,vpw = vpw(:)) 
      CALL gf_makePotProjector(layer,layers,cell,stars,.FALSE.          &
     &     ,curvyproj)                                                  
                                                                        
      IF (layer == 1) THEN 
         WRITE(6,*) "Save Coulomb potential for left boundary" 
         !left side                                                     
         v_bound = MATMUL(curvyproj(:,:,1),vpw(:)) 
         CALL gf_SAVEcoulboundary(1,1,v_bound) 
      ENDIF 
      IF (layer == layers%num_layers) THEN 
         !right side                                                    
         WRITE(6,*) "Save Coulomb potential for right boundary" 
         v_bound = MATMUL(curvyproj(:,:,2),vpw(:)) 
         CALL gf_SAVEcoulboundary(layer,2,v_bound) 
      ENDIF 
      DEALLOCATE(curvyproj,vpw, v_bound) 
                                                                        
      END SUBROUTINE 
      !>                                                                
                                                                        
      END                                           
