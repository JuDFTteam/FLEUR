!--------------------------------------------------------------------------------
! Copyright (c) 2026 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions 
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------
      MODULE m_gf_propaemb 
      USE m_constants, ONLY: oUnit
      USE m_gf_math
          IMPLICIT NONE
      PRIVATE
      PUBLIC :: gf_propaemb
      !ill-conditioned propagations reported so far
      INTEGER, SAVE :: nwarn = 0
      CONTAINS 
      SUBROUTINE gf_propaemb(                                           &
     &                     l_invs,nv2,embpot_in,tmat,                   &
     &                     embpot_out)                                  
!****************************************************                   
!     Propagate the embedding potential using the                       
!     transfer matrix "tmat".                                           
!     l_invs=.true. => use the inverse of "tmat".                       
!     Frank Freimuth, November 2007                                     
!****************************************************                   
      IMPLICIT NONE 
      LOGICAL,INTENT(IN)::l_invs 
      INTEGER,INTENT(IN)::nv2 
      COMPLEX,INTENT(IN)::embpot_in(:,:) 
      COMPLEX,INTENT(IN)::tmat(:,:) 
      COMPLEX,INTENT(OUT)::embpot_out(:,:) 
                                                                        
      COMPLEX,ALLOCATABLE:: t11(:,:),t22(:,:),t12(:,:),t21(:,:) 
      COMPLEX,ALLOCATABLE:: tmat_tmp(:,:),ew(:)
      REAL :: anorm,cond
                                                                        
      ALLOCATE( t11(nv2,nv2) ) 
      ALLOCATE( t12(nv2,nv2) ) 
      ALLOCATE( t21(nv2,nv2) ) 
      ALLOCATE( t22(nv2,nv2) ) 
      ALLOCATE( ew(nv2))
      IF(l_invs)THEN 
         ALLOCATE( tmat_tmp(2*nv2,2*nv2) ) 
         tmat_tmp=mat_inverse(tmat) 
         t11(1:nv2,1:nv2)=tmat_tmp(1:nv2,1:nv2) 
         t12(1:nv2,1:nv2)=tmat_tmp(1:nv2,nv2+1:2*nv2) 
         t21(1:nv2,1:nv2)=tmat_tmp(nv2+1:2*nv2,1:nv2) 
         t22(1:nv2,1:nv2)=tmat_tmp(nv2+1:2*nv2,nv2+1:2*nv2) 
         DEALLOCATE( tmat_tmp ) 
      ELSE 
         t11(1:nv2,1:nv2)=tmat(1:nv2,1:nv2) 
         t12(1:nv2,1:nv2)=tmat(1:nv2,nv2+1:2*nv2) 
         t21(1:nv2,1:nv2)=tmat(nv2+1:2*nv2,1:nv2) 
         t22(1:nv2,1:nv2)=tmat(nv2+1:2*nv2,nv2+1:2*nv2) 
      ENDIF 
!      t12=matmul(t12,embpot_in)                                        
      CALL cblas_matmul_sq(t12,embpot_in) 
!      t22=matmul(t22,embpot_in)                                        
      CALL cblas_matmul_sq(t22,embpot_in) 
      t11=t11+t12 
      anorm=MAXVAL(SUM(ABS(t11),DIM=1))
      t11=mat_inverse(t11)
      !1-norm condition of t11+t12*Sigma: thick layers lose the evanescent states
      cond=anorm*MAXVAL(SUM(ABS(t11),DIM=1))
      IF (cond>1.0E12.AND.nwarn<10) THEN
         nwarn=nwarn+1
         WRITE(oUnit,"(a,es10.2,a)") "WARNING: propagation of the embedding potential ill-conditioned, cond=", &
     &        cond," (thinner layers or a smaller |k+G| cutoff help)"
      ENDIF

      t21=t21+t22
!      embpot_out=matmul(t11,t21)                                       
      CALL cblas_matmul_sq2(t21,t11,embpot_out) 
      END SUBROUTINE gf_propaemb 
      END                                           
