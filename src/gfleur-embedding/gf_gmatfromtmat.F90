!--------------------------------------------------------------------------------
! Copyright (c) 2026 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions 
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------
      MODULE m_gf_gmatfromtmat 
      use m_juDFT
          IMPLICIT NONE
!***************************************************                    
!     Obtain the surface projected Green's function                     
!     from the transfer matrix. The Green's function                    
!     is the one with van Neumann boundary conditions.                  
!     Frank Freimuth, September 2007                                    
!***************************************************                    
      CONTAINS 
      SUBROUTINE gf_gmatfromtmat(                                       &
     &                           nv2,en,nk,jspin,lapw,                  &
     &                           l_embed,cell,                          &
     &                           gij,                                   &
     &                           tmat)                                  
      USE m_gf_math, ONLY: mat_inverse 
      USE m_gf_io2dmat 
      USE m_gf_types, ONLY: t_lapw, t_cell 
      USE m_gf_projembed 
                                                                        
      IMPLICIT NONE 
      INTEGER, INTENT(IN)           :: nv2,en,nk,jspin 
      TYPE(t_lapw),INTENT(IN)       :: lapw 
      LOGICAL, INTENT(IN)           :: l_embed 
      TYPE(t_cell),INTENT(IN)       :: cell 
      COMPLEX, INTENT(OUT)          :: gij(2*nv2,2*nv2) 
      COMPLEX, OPTIONAL, INTENT(IN) :: tmat(2*nv2,2*nv2) 
                                                                        
      COMPLEX, ALLOCATABLE :: tmat11(:,:) 
      COMPLEX, ALLOCATABLE :: tmat22(:,:) 
      COMPLEX, ALLOCATABLE :: tmat12(:,:) 
      COMPLEX, ALLOCATABLE :: tmat21(:,:) 
      COMPLEX, ALLOCATABLE :: tmat_tmp(:,:) 
                                                                        
      ALLOCATE( tmat11(nv2,nv2) ) 
      ALLOCATE( tmat22(nv2,nv2) ) 
      ALLOCATE( tmat12(nv2,nv2) ) 
      ALLOCATE( tmat21(nv2,nv2) ) 
                                                                        
      IF(present(tmat))THEN 
         tmat11(1:nv2,1:nv2)=tmat(1:nv2,1:nv2) 
         tmat12(1:nv2,1:nv2)=tmat(1:nv2,nv2+1:2*nv2) 
         tmat21(1:nv2,1:nv2)=tmat(nv2+1:2*nv2,1:nv2) 
         tmat22(1:nv2,1:nv2)=tmat(nv2+1:2*nv2,nv2+1:2*nv2) 
      ELSE 
         ALLOCATE( tmat_tmp(2*nv2,2*nv2) ) 
! read transfer-matrix from file                                        
         IF(.NOT.(gf_read2dmat(                                         &
     &         IO2D_TMAT,1,0,en,nk,jspin,lapw,tmat_tmp)))               &
     &          CALL juDFT_error("tmat1",calledby="gf_gmatfromtmat.F90")
         tmat11(1:nv2,1:nv2)=tmat_tmp(1:nv2,1:nv2) 
         tmat12(1:nv2,1:nv2)=tmat_tmp(1:nv2,nv2+1:2*nv2) 
         tmat21(1:nv2,1:nv2)=tmat_tmp(nv2+1:2*nv2,1:nv2) 
         tmat22(1:nv2,1:nv2)=tmat_tmp(nv2+1:2*nv2,nv2+1:2*nv2) 
         DEALLOCATE(tmat_tmp) 
      ENDIF 
      gij(1:nv2,nv2+1:2*nv2)=-2.0*mat_inverse(tmat21) 
      gij(1:nv2,1:nv2)=matmul(gij(1:nv2,nv2+1:2*nv2),tmat22) 
      gij(nv2+1:2*nv2,nv2+1:2*nv2)=matmul(tmat11,gij(1:nv2,nv2+1:2*nv2)) 
      gij(nv2+1:2*nv2,1:nv2)=2.0*tmat12+matmul(tmat11,gij(1:nv2,1:nv2)) 
      DEALLOCATE(tmat11,tmat12,tmat21,tmat22) 
                                                                        
      IF(l_embed)THEN 
           CALL gf_projembed(                                           &
     &           lapw,en,nk,jspin,cell,.FALSE.,                         &
     &           .FALSE.,.FALSE.,1,                                     &
     &           gij(:,:))                                              
      ENDIF 
                                                                        
      END SUBROUTINE 
      END                                           
