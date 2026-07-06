!--------------------------------------------------------------------------------
! Copyright (c) 2026 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions 
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------
!******************************************************************     
      MODULE m_gf_addembnhr 
          IMPLICIT NONE
      CONTAINS 
      SUBROUTINE gf_addembnhr(                                          &
     &     cell,en,nk,jspin,layer,lapw,lapw_gf,                                 &
     &     g)                                                           
      USE m_gf_embedding 
      USE m_gf_types 
      USE m_gf_curvy2dprojector 
                                                                        
      IMPLICIT NONE 
      INTEGER,INTENT(IN) :: en,nk,jspin,layer 
      TYPE(t_lapw),INTENT(IN)::lapw 
      TYPE(t_lapw_gf),INTENT(IN) :: lapw_gf
      TYPE(t_cell),INTENT(IN)     :: cell 
      COMPLEX,INTENT(INOUT)::g(:,:) 
                                                                        
      COMPLEX,ALLOCATABLE::projtwodright(:,:) 
      COMPLEX,ALLOCATABLE::g1(:,:),g2(:,:) 
      COMPLEX,ALLOCATABLE::projtwodleft(:,:) 
      COMPLEX,ALLOCATABLE::bigg(:,:) 
      COMPLEX,ALLOCATABLE::gnoughtp(:,:) 
      INTEGER matsize 
                                                                        
      matsize=lapw_gf%nmat_sphere
      ALLOCATE( g1(lapw_gf%nv2_tot,lapw_gf%nv2_tot) ) 
      ALLOCATE( g2(lapw_gf%nv2_tot,lapw_gf%nv2_tot) ) 
      ALLOCATE( bigg(2*lapw_gf%nv2_tot,2*lapw_gf%nv2_tot) ) 
      ALLOCATE( gnoughtp(matsize,2*lapw_gf%nv2_tot)) 
      CALL gf_getemb(g1,g2,layer,en,nk,jspin,lapw,lapw_gf) 
      bigg(:,:)=cmplx(0.0,0.0) 
      bigg(1:lapw_gf%nv2_tot,1:lapw_gf%nv2_tot)=                              &
     &        g1(1:lapw_gf%nv2_tot,1:lapw_gf%nv2_tot)                         
      bigg(lapw_gf%nv2_tot+1:2*lapw_gf%nv2_tot,                               &
     &     lapw_gf%nv2_tot+1:2*lapw_gf%nv2_tot)=                              &
     &        g2(1:lapw_gf%nv2_tot,1:lapw_gf%nv2_tot)                         
      DEALLOCATE(g1,g2) 
      CALL gf_curvy2dprojector(layer,cell,lapw,lapw_gf,.TRUE.) 
                                                                        
                                                                        
      ALLOCATE( projtwodright(matsize,2*lapw_gf%nv2_tot) ) 
      CALL gf_get_curvy2dproject(lapw,lapw_gf,projtwodright,.TRUE.) 
!      projtwodright(1:matsize,1:lapw_gf%nv2_tot)=                         
!     =    transpose(conjg(curvyproject(1:lapw_gf%nv2_tot,1:matsize,1)))   
!      projtwodright(1:matsize,lapw_gf%nv2_tot+1:2*lapw_gf%nv2_tot)=          
!     =    transpose(conjg(curvyproject(1:lapw_gf%nv2_tot,1:matsize,2)))   
      CALL zgemm('N','N',matsize,2*lapw_gf%nv2_tot,2*lapw_gf%nv2_tot,&
     &                    cmplx(1.0,0.0),projtwodright,                 &
     &                    matsize,bigg,2*lapw_gf%nv2_tot,                  &
     &                    cmplx(0.0,0.0),gnoughtp,matsize)              
      DEALLOCATE(bigg) 
                                                                        
      CALL zgemm('N','C',matsize,matsize,2*lapw_gf%nv2_tot,       &
     &                   cmplx(1.0,0.0),gnoughtp,                       &
     &                   matsize,projtwodright,matsize,cmplx(1.0,0.0),  &
     &                   g,matsize)                                     
                                                                        
      DEALLOCATE(projtwodright) 
      DEALLOCATE(gnoughtp) 
      END SUBROUTINE gf_addembnhr 
      END                                           
