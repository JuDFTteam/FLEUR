!--------------------------------------------------------------------------------
! Copyright (c) 2026 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions 
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------
      MODULE m_gf_landauer 
          IMPLICIT NONE
      CONTAINS 
      SUBROUTINE gf_landauer(nv2,en,nk,jspin,cell,sym,lapw,bkpts,gij,mpi)
!*********************************************************************  
!     Subroutine to calculate current within the energy loop            
!     calls subroutines from gf_current                                 
!                                           Daniel Wortmann             
!*********************************************************************  
      USE m_gf_embedding,ONLY:gf_getemb 
      USE m_gf_types 
      USE m_gf_writetrans,ONLY:writetrans 
      USE m_gf_current 
#ifdef CPP_CUOVLP                                                       
      USE m_gf_curvy2dprojector,ONLY:basisoverlaps 
#endif                                                                  
      IMPLICIT NONE 
      INTEGER ,INTENT(IN)               :: nv2,en,nk,jspin 
      type(t_mpi),intent(in)::mpi

      TYPE(t_sym),INTENT(IN)            :: sym 
      TYPE(t_cell),INTENT(IN)           :: cell 
      TYPE(t_lapw),INTENT(IN)           :: lapw 
      REAL,INTENT(IN)                   :: bkpts(:,:) 
      COMPLEX,DIMENSION(:,:),INTENT(IN) :: gij 
                                                                        
      COMPLEX,ALLOCATABLE :: G1(:,:),G2(:,:) 
      REAL                :: j(3+nv2) 
      LOGICAL l_nohelpregion 
                                                                        
      ALLOCATE(G1(nv2,nv2),g2(nv2,nv2)) 
                                                                        
      !Load embedding potentials                                        
                                                                        
      CALL gf_getemb(g1,g2,1,en,nk,jspin,lapw) 
      !calculate with landauer                                          
#ifdef CPP_CUOVLP                                                       
      l_nohelpregion=.TRUE. 
      IF(l_nohelpregion)THEN 
         g1=matmul(basisoverlaps(:,:,1),g1) 
         g2=matmul(basisoverlaps(:,:,2),g2) 
      ENDIF 
#endif                                                                  
      CALL gf_landauer2plane(gij(1:nv2,nv2+1:2*nv2),g1,g2,j(1),j(2)) 
                                                                        
      !channel decomposition                                            
!      CALL gf_channels(g1,g2,gij(1:nv2,nv2+1:2*nv2),j(4:),j(3),lapw)   
                                                                        
      CALL writetrans(en,nk,jspin,bkpts,sym,cell,3,j(1:2),mpi)
                                                                        
      DEALLOCATE(g1,g2) 
                                                                        
      END SUBROUTINE 
                                                                        
                                                                        
      END                                           
