!--------------------------------------------------------------------------------
! Copyright (c) 2026 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions 
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------
      MODULE m_gf_OPTIONAL 
      use m_juDFT
      IMPLICIT NONE
!-----------------------------------------------                        
! DESC:Collection of tools and service routines                         
!                 Daniel Wortmann, (05-08-17)                           
!-----------------------------------------------                        
      PRIVATE 
      PUBLIC gf_OPTIONAL 
      CONTAINS 
      !<-- S: gf_optional(jspins,cell,gfinp,lapw,noco,kpts,sym)         
      SUBROUTINE gf_optional(jspins,atoms,input,nococonv,cell,gfinp,lapw,lapw_gf,noco,kpts,sym) 
!-----------------------------------------------                        
!    Checks for various files and calls corresponding routines          
!           (last modified: 05-08-17) D. Wortmann                       
!-----------------------------------------------                        
      USE m_gf_kptsp1x1 
      USE m_gf_types 
      USE m_gf_embdos 
      IMPLICIT NONE 
      !<--Arguments                                                     
      INTEGER,INTENT(IN)         :: jspins 
      TYPE(t_lapw),INTENT(INOUT) :: lapw 
      TYPE(t_embinp),INTENT(IN)   :: gfinp 
      TYPE(t_cell),INTENT(IN)    :: cell 
      TYPE(t_atoms),INTENT(IN)   :: atoms
      TYPE(t_input),INTENT(IN)   :: input
      TYPE(t_nococonv),INTENT(IN):: nococonv
      TYPE(t_lapw_gf),INTENT(INOUT) :: lapw_gf
      TYPE(t_noco),INTENT(IN)    :: noco 
      TYPE(t_kpts),INTENT(IN)    :: kpts 
      TYPE(t_sym),INTENT(IN)     :: sym 
      !>                                                                
      !<-- Locals                                                       
      LOGICAL             :: l_exist 
      !>                                                                
                                                                        
      !A p1x1 kpts-file should be generated                             
      INQUIRE(FILE='gf_pot.kpts.hdf',EXIST=l_exist) 
      IF (l_exist) THEN 
         CALL gf_kptsp1x1(cell,sym,kpts,gfinp) 
          CALL juDFT_end("kpts-generated")
      ENDIF 
      INQUIRE(FILE ='embdos_inp',EXIST = l_exist) 
      IF (l_exist) THEN 
         CALL gf_embdos(jspins,atoms,input,sym,nococonv,lapw,lapw_gf,kpts,cell,noco,gfinp)
          CALL juDFT_end("embdos")
      ENDIF 
      END SUBROUTINE 
      !>                                                                
                                                                        
      END                                           
