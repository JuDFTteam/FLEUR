!--------------------------------------------------------------------------------
! Copyright (c) 2026 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions 
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------
      MODULE m_gf_mkposdef 
      use m_juDFT
          IMPLICIT NONE
#include "realcomplex.h"                                                
#include "cpp_double.h"
      CONTAINS 
      !<-- S:priv_mkposdef(matsize,l_inversion,over,over_real,uhumatrix_
      SUBROUTINE gf_mkposdef(matsize,over)
!-----------------------------------------------                        
!in some cases the overlap matrix is not positive definite              
!the following piece of code helps here                                 
!Calculate the eigenvalues e(j) below toler and the corresponding       
!eigenvectors v(j) of the overlap matrix.                               
!Shift the eigenvalues by adding (toler-e(j))*v(j)*v^H(j) to            
!the overlap matrix.                                                    
!     Frank Freimuth, January-February 2007                             
!           (last modified: 07-02-22) D. Wortmann                       
!-----------------------------------------------                        
      USE m_gf_spectrum 
      IMPLICIT NONE 
      !<--Arguments                                                     
      INTEGER,INTENT(IN)     :: matsize 
      COMPLEX,INTENT(INOUT)  :: over(:,:) 
      !>                                                                
      !<-- Locals                                                       
      COMPLEX,ALLOCATABLE :: work(:) 
      REAL,ALLOCATABLE    :: rwork(:) 
      INTEGER,ALLOCATABLE :: iwork(:) 
      INTEGER,ALLOCATABLE :: ifail(:) 
      COMPLEX,ALLOCATABLE :: overeigvec1(:,:) 
      REAL                :: toler,abstol,vl,vu 
      INTEGER             :: lwork,liwork,lrwork,info,ne,ind1 
                                                                        
      !>                                                                
!diagonalize overlap matrix                                             
      PRINT*,"make metric positive definite" 
      lwork = 2*matsize
      ALLOCATE(work(lwork))
      liwork = 5*matsize
      ALLOCATE(iwork(liwork))
      lrwork = 7*matsize
      ALLOCATE(rwork(lrwork))
      ALLOCATE(ifail(matsize))
                         !eigenvalues below toler will be shifted       
      toler = 1.0E-12
      vl =-1.0E12
      vu = toler
      abstol = 2*TINY(abstol)
      ALLOCATE(overeigvec1(matsize,matsize))
      uhumatrix(:,:) = over(:,:)

      CALL CPP_LAPACK_cheevx('V','V','U',matsize,uhumatrix,       &
     &        matsize,                                                  &
     &        vl,vu,1,matsize,abstol,ne,uhueigval,overeigvec1,matsize,  &
     &        work,lwork,rwork,iwork,ifail,info)                        
      IF(info /= 0) CALL juDFT_error("overdiag: cheevx",calledby="gf_mkposdef.F90")
      DEALLOCATE(iwork)
      DEALLOCATE(work)
      DEALLOCATE(rwork)
      DEALLOCATE(ifail)
!remove problematic eigenvalues                                         
      PRINT*,"eigenvalues that will be removed:"
      DO ind1 = 1,ne
          PRINT*,uhueigval(ind1)
          uhueigval(ind1) = toler-uhueigval(ind1)
      ENDDO
      DO ind1 = 1,ne
          CALL CPP_BLAS_cher('U',matsize,uhueigval(ind1),        &
     &           overeigvec1(1,ind1),1,over,matsize)                    
      ENDDO
      DEALLOCATE(overeigvec1)
      END SUBROUTINE gf_mkposdef 
      !>                                                                
      END                                           
