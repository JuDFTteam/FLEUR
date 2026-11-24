!--------------------------------------------------------------------------------
! Copyright (c) 2026 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions 
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------
      MODULE m_gf_povrayColors 
          IMPLICIT NONE
      CONTAINS 
      !<-- F: getColors                                                 
      FUNCTION getColors() RESULT(Col) 
!-----------------------------------------------                        
! Colors for atoms                                                      
!             (last modified: 2004-00-00) D. Wortmann                   
!-----------------------------------------------                        
      IMPLICIT NONE 
      !<--Arguments                                                     
      REAL                :: Col(3,103) 
      !>                                                                
      !<-- Locals                                                       
      INTEGER             :: rgb(3,103) 
      !>                                                                
                                                                        
      !Set all colors to default                                        
      rgb(1,:)=190 
      rgb(2,:)=190 
      rgb(3,:)=190 
                                                                        
      !Set some colors !Please add more                                 
      rgb(:,1)=(/255,20,147/) 
      rgb(:,2)=(/250,235,215/) 
      rgb(:,3) = (/255,192,203/) 
      rgb(:,4) = (/178,34,34/) 
      rgb(:,5) = (/34,139,34/) 
      rgb(:,6) = (/0,255,0/) 
      rgb(:,7) = (/112,128,144/) 
      rgb(:,8) = (/0,191,255/) 
      rgb(:,9) = (/255,0,0/) 
      rgb(:,10) = (/218,165,32/) 
      rgb(:,11) = (/255,105,180/) 
      rgb(:,12) = (/0,0,255/) 
      rgb(:,13) = (/34,139,34/) 
      rgb(:,14) = (/190,190,190/) 
      rgb(:,15) = (/218,165,32/) 
      rgb(:,16) = (/255,165,0/) 
      rgb(:,17) = (/255,255,0/) 
      rgb(:,18) = (/0,255,0/) 
      rgb(:,19) = (/255,192,203/) 
      rgb(:,20) = (/255,20,14/) 
      rgb(:,21) = (/128,128,128/) 
                                                                        
      !Calculate rgb values                                             
                                                                        
      col=real(rgb)/255.0 
                                                                        
      END FUNCTION 
      !>                                                                
                                                                        
      END                                           
