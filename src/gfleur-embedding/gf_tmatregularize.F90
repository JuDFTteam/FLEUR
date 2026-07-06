!--------------------------------------------------------------------------------
! Copyright (c) 2026 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions 
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------
      MODULE m_gf_tmatregularize 
          IMPLICIT NONE
#include "cpp_double.h"                                                 
      CONTAINS 
      SUBROUTINE gf_tmatregularize(tmat) 
      USE m_gf_math 
      IMPLICIT NONE 
      COMPLEX, INTENT(INOUT):: tmat(:,:) 
                                                                        
      COMPLEX,ALLOCATABLE :: ev(:,:),ev_inv(:,:) 
      COMPLEX,ALLOCATABLE :: ew(:) 
      COMPLEX,ALLOCATABLE :: t_temp(:,:) 
      INTEGER i,nv2 
      COMPLEX lambda 
      REAL cutoff 
      LOGICAL l_exist 
                                                                        
                                                                        
      cutoff=0.4E-1 
                                                                        
      INQUIRE(FILE='regularize',EXIST=l_exist) 
      IF(.NOT.l_exist)THEN 
         RETURN 
      ELSE 
         OPEN(888,FILE='regularize') 
         READ(888,*)cutoff 
         PRINT*,"cutoff=",cutoff 
         CLOSE(888) 
      ENDIF 
                                                                        
                                                                        
      nv2=size(tmat,1) 
      ALLOCATE( ev_inv(nv2,nv2) ) 
      ALLOCATE( ev(nv2,nv2) ) 
      ALLOCATE( ew(nv2) ) 
      ALLOCATE( t_temp(nv2,nv2) ) 
      t_temp = Tmat 
      CALL eigenvalues(T_temp,ew,ev) 
      ev_inv=mat_inverse(ev) 
      DO i=1,nv2 
         lambda=ew(i) 
         IF(abs(lambda)<cutoff)THEN 
            lambda=cmplx(cutoff,0.0) 
         ELSEIF(abs(lambda)>(1.0/cutoff))THEN 
            lambda=cmplx(1.0/cutoff,0.0) 
         ENDIF 
         ev(:,i)=lambda*ev(:,i) 
      ENDDO 
!      tmat=matmul(ev,ev_inv)                                           
      CALL CPP_BLAS_cgemm('N','N',nv2,nv2,                              &
     &      nv2,cmplx(1.0,0.0),ev,nv2,                                  &
     &    ev_inv,nv2,cmplx(0.0,0.0),tmat,nv2)                           
      DEALLOCATE(ev,ew,t_temp) 
      END SUBROUTINE 
      END                                           
