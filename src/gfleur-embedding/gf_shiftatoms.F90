!--------------------------------------------------------------------------------
! Copyright (c) 2026 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions 
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------
      MODULE m_gf_shiftatoms 
          IMPLICIT NONE
      CONTAINS 
      !<-- S:                                                           
      SUBROUTINE gf_shiftatoms(atpos,ucell) 
!-----------------------------------------------                        
!   shift all the atomic positions in 2D-Unit cell such that they are in
!   positive (x,y) range                                                
!           (last modified: 2004-00-00) D. Wortmann                     
!-----------------------------------------------                        
      IMPLICIT NONE 
      !<--Arguments                                                     
      REAL,INTENT(INOUT)    :: atpos(:,:) 
      REAL   ,INTENT(IN)    :: ucell(2,2) 
      !>                                                                
      !<-- Locals                                                       
      INTEGER             :: n 
      REAL                :: uinv(2,2) 
      REAL                :: pos(2) 
      !>                                                                
                                                                        
      !<-- inverse of ucell                                             
      uinv=ucell 
      uinv(1,1)=-1.*ucell(2,2) 
      uinv(2,2)=-1.*ucell(1,1) 
      uinv=uinv/(ucell(1,2)*ucell(2,1)-ucell(1,1)*ucell(2,2)) 
      !>                                                                
                                                                        
      DO n = 1,SIZE(atpos,2) 
         !to internal coords                                            
         pos = MATMUL(uinv,atpos(:,n)) 
         DO WHILE (ANY(pos(:)<0)) 
            pos(:) = pos(:)+(/1.0,1.0/) 
         ENDDO 
         !subtract vectors if possible                                  
         DO WHILE (pos(1)>1.-epsilon(0.)) 
            pos(1) = pos(1)-1.0 
         ENDDO 
         DO WHILE (pos(2)>1.-epsilon(0.)) 
            pos(2) = pos(2)-1.0 
         ENDDO 
         !Go back to real space                                         
         atpos(:,n)=matmul(ucell,pos) 
      ENDDO 
      END SUBROUTINE 
      !>                                                                
      END                                           
