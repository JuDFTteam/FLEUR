!--------------------------------------------------------------------------------
! Copyright (c) 2026 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions 
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------
      MODULE m_gf_tlmplminit 
      use m_juDFT
          IMPLICIT NONE
      CONTAINS 
      SUBROUTINE gf_tlmplminit(jspins,atoms,noco,tlmplm_data,raddata) 
!-----------------------------------------------                        
!   allocate large arrays for radial data and for the                   
!   tlmplm stuff for later use                                          
!           (last modified: 05-04-08) D. Wortmann                       
!-----------------------------------------------                        
      USE m_gf_types 
      USE m_fleur_hsetup 
      IMPLICIT NONE 
      !<--Arguments                                                     
      INTEGER,INTENT(IN)           :: jspins 
      TYPE(t_atoms),INTENT(IN)     :: atoms 
      TYPE(t_noco),INTENT(IN)      :: noco 
      TYPE(t_tlmplm),INTENT(INOUT) :: tlmplm_DATA 
      TYPE(t_raddata),INTENT(INOUT) :: raddata 
      !>                                                                
      !<-- Locals                                                       
      INTEGER             :: j,err,lmd,lmplmd 
      lmd = MAXVAL(Atoms%lmax0)*(MAXVAL(Atoms%lmax0)+2) 
      lmplmd = (lmd*(lmd+3))/2 
                                                                        
      !>                                                                
      !<-- Allocate raddata                                             
                                                                        
      IF (.NOT.ALLOCATED(raddata%duds)) THEN 
         ALLOCATE(raddata%ulos(atoms%nlod,atoms%ntype,jspins),          &
     &        raddata%dulos(atoms%nlod,atoms%ntype,jspins),             &
     &        raddata%uulon(atoms%nlod,atoms%ntype,jspins),             &
     &        raddata%dulon(atoms%nlod,atoms%ntype,jspins),             &
     &        raddata%uloulopn(atoms%nlod,atoms%nlod,atoms%ntype,jspins)&
     &        ,raddata%us(0:MAXVAL(atoms%lmax0),atoms%ntype,jspins)     &
     &        ,raddata%uds(0:MAXVAL(atoms%lmax0),atoms%ntype,jspins)    &
     &        ,raddata%dus(0:MAXVAL(atoms%lmax0),atoms%ntype,jspins)    &
     &        ,raddata%duds(0:MAXVAL(atoms%lmax0),atoms%ntype,jspins)   &
     &        ,raddata%ddn(0:MAXVAL(atoms%lmax0),atoms%ntype,jspins),   &
     &        STAT = err)                                               
         IF (err /= 0) THEN 
            WRITE (6,*)                                                 &
     &              'gf_hs1lapw: an error occured during allocation of' 
            WRITE (6,*) 'the raddata-arrays: ',err 
            CALL juDFT_error                                               &
     &           ("Error occured during allocation of the raddata-arrays"&
     &           ,calledby ="gf_tlmplminit.F90")  
         ENDIF 
      ENDIF 
                                                                        
      !>                                                                
      !<-- Allocate storage for the  t(l'm',lm) matrices                
                                                                        
      IF (.NOT.ALLOCATED(tlmplm_data%tuu)) THEN 
         j = 1 ; IF (noco%l_noco) j = 2 
         ALLOCATE(tlmplm_data%tuu(0:lmplmd,atoms%ntype,j),STAT = err) 
         ALLOCATE(tlmplm_data%tud(0:lmplmd,atoms%ntype,j),STAT = err) 
         ALLOCATE(tlmplm_data%tdd(0:lmplmd,atoms%ntype,j),STAT = err) 
         ALLOCATE(tlmplm_data%tdu(0:lmplmd,atoms%ntype,j),STAT = err) 
         ALLOCATE(tlmplm_data%tdulo(0:lmd,-atoms%llod:atoms%llod        &
     &        ,SUM(atoms%nlo),j),STAT = err)                            
         ALLOCATE(tlmplm_data%tuulo(0:lmd,-atoms%llod:atoms%llod        &
     &        ,SUM(atoms%nlo),j),STAT = err)                            
         ALLOCATE(tlmplm_data%tuloulo(-atoms%llod:atoms%llod,           &
     &        -atoms%llod:atoms%llod,SUM(atoms%nlo*(atoms%nlo+1)/2),j   &
     &        ),STAT = err)                                             
         ALLOCATE(tlmplm_data%ind(0:lmd,0:lmd,atoms%ntype,j),STAT =     &
     &        err )                                                     
         IF (err /= 0) THEN 
            WRITE (6,*)                                                 &
     &           'gf_hs1lapw: an error occured during allocation of'    
            WRITE (6,*) 'the tuu, tdd etc.s: ',err 
            CALL juDFT_error                                              &
     &           ("Error occured during allocation of the tuu, tdd etc" &
     &           ,calledby ="gf_tlmplminit.F90")
         ENDIF 
      ENDIF 
      tlmplm_data%tdulo(:,:,:,:) = CMPLX(0.0,0.0) 
      tlmplm_data%tuulo(:,:,:,:) = CMPLX(0.0,0.0) 
      tlmplm_data%tuloulo(:,:,:,:) = CMPLX(0.0,0.0) 
                                                                        
                                                                        
      !>                                                                
      END SUBROUTINE gf_tlmplminit 
      !>                                                                
      END                                           
