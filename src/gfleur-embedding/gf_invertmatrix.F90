!--------------------------------------------------------------------------------
! Copyright (c) 2026 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions 
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------
      MODULE m_gf_invertmatrix 
!----------------------------------------------------------------       
!  Module contains subroutines to invert large matrices                 
!                                                                       
!   gf_invertmatrix : Public frontend for gf_fleur                      
!   gf_inv_cmat     : inverts a large general complex matrix            
!   gf_inv_rmat     : inverts a large general real matrix               
!   gf_inv_hmat     : inverts a large hermitian matrix                  
!   gf_inv_Smat     : inverts a large symmetric matrix                  
!                                                                       
!---------------------------------------------------------------        
      use m_juDFT 
      IMPLICIT NONE
      CONTAINS 
                                                                        
      !<-- S: gf_invertmatrix(A,l_inversion,l_addemb,l_realenergy)      
                                                                        
      SUBROUTINE gf_invertmatrix(A,l_inversion,l_addemb,l_realenergy) 
!-----------------------------------------------                        
!       chooses correct inversion routine for GFLEUR                    
!           (last modified: 2004-00-00) D. Wortmann                     
!-----------------------------------------------                        
      IMPLICIT NONE 
      !<--Arguments                                                     
      COMPLEX,INTENT(INOUT) :: A(:,:) 
      LOGICAL,INTENT(IN)    :: l_inversion,l_addemb,l_realenergy 
      !>                                                                
      !<-- locals                                                       
      REAL, ALLOCATABLE :: B(:,:) 
      REAL ,PARAMETER   :: eps = 1.E-9 
      LOGICAL           :: isreal,issym 
      INTEGER           :: i 
      !>                                                                
                                                                        
      !<-- tests                                                        
      isreal  = .FALSE.;issym = .FALSE. 
      IF (l_realenergy) THEN 
         !if no real energy we have a general complex matrix!           
         !Otherwise the matrix should be hermitian(symmetric)           
         !However, this symmetry might be broken by the embedding       
         !potential and therefore we test it for column 3 of the matrix 
         IF (l_inversion.AND..NOT.l_addemb) THEN 
                                       !real and symmetric!             
            isreal=.TRUE.;issym=.TRUE. 
         ELSEIF(l_inversion) THEN 
            !We should have a real and symmetric matrix, but better test
            IF (ALL(ABS(AIMAG(A(3,:)))<eps)) isreal=.TRUE. 
            IF (ALL(AIMAG(A(3,:)) == AIMAG(A(:,3)))) issym = .TRUE. 
         ELSE 
            !no inversion symmetry -> no real matrix                    
            !but should be still hermitian                              
!            IF (ALL(abs(CONJG(A(3,:))-A(:,3)).lt.1.E-12)) then         
!               issym=.TRUE.                                            
!            else                                                       
!               do i=1,size(A,1)                                        
!                  if(conjg(a(3,i)).ne.a(i,3))then                      
!                     print*,"i=",i                                     
!                     print*,"a(i,3)=",a(i,3)                           
!                     print*,"conjg(a(3,i))=",conjg(a(3,i))             
!                  endif                                                
!               enddo                                                   
!            endif                                                      
         ENDIF 
      ENDIF 
      !>                                                                
                                                                        
      IF (isreal) THEN 
         !<-- we have a real-matrix!                                    
         ALLOCATE(B(SIZE(A,1),SIZE(A,2))) 
         B = REAL(A) 
         IF (issym) THEN 
                                !real-symmtric                          
            CALL gf_inv_smat(B) 
         ELSE 
                                !general real                           
            CALL gf_inv_rmat(B) 
         ENDIF 
         A = CMPLX(B,0.0) 
         DEALLOCATE(B) 
         !>                                                             
      ELSE 
         !<-- complex-matrix                                            
         IF (issym) THEN 
            CALL gf_inv_hmat(A) 
         ELSE 
            CALL gf_inv_cmat(A) 
         ENDIF 
         !>                                                             
      ENDIF 
                                                                        
      END SUBROUTINE 
                                                                        
      !>                                                                
                                                                        
      !<--S: gf_inv_cmat(A)                                             
      SUBROUTINE gf_inv_cmat(A) 
!***********************************************************************
!     Subroutine encapsules the LAPACK calls needed to invert a large ma
!                                                                       
!     Daniel Wortmann, Juelich 2001                                     
!***********************************************************************
      IMPLICIT NONE 
#include "cpp_double.h"                                                 
      COMPLEX,INTENT(INOUT)::A(:,:) 
!     LOCALS                                                            
      INTEGER             :: nmat 
      INTEGER             :: lwork,info,ilaenv 
      INTEGER             :: ipiv(size(A,1)) 
      COMPLEX,ALLOCATABLE :: work(:) 
                                                                        
      EXTERNAL ilaenv 
                                                                        
      nmat = SIZE(A,1) 
      !<-- Before starting determine best size of lapack work-array!    
                                                                        
#ifndef CPP_APC                                                         
      !on PC no ilaenv-call is available...                             
      lwork=ilaenv(1,'cgetri',' ',Nmat,Nmat,-1,-1) 
      lwork=Nmat*lwork 
      IF (lwork<0)                                                   &
     &     CALL juDFT_error                                               &
     &     ('gf_gfcn: failed to determine LAPACK work-size')            
#else                                                                   
      lwork=Nmat 
#endif                                                                  
      ALLOCATE(work(lwork),STAT=info) 
      IF (info/=0) THEN 
         WRITE(6,*)                                                     &
     & 'WARNING, failed to ALLOCATE LAPACK workspace in gf_invertmatrix'
         WRITE(*,*)                                                     &
     &        'WARNING, failed to ALLOCATE LAPACK workspace'            
         WRITE(6,*)                                                     &
     &        'Using default value'                                     
         lwork=Nmat 
         ALLOCATE(work(lwork)) 
      ENDIF 
                                                                        
      !>                                                                
                                                                        
                                                                        
      !<-- Now inverte the Matrix!                                      
                                                                        
      WRITE(*,*) 'Non-Hermitian complex matrix inversion' 
      CALL CPP_LAPACK_cgetrf(Nmat,Nmat,A,Nmat,ipiv,info) 
      IF ( info /= 0 ) CALL juDFT_error('cgetrf (inversion of H-zS)') 
      CALL CPP_LAPACK_cgetri(Nmat,A,Nmat,ipiv,work,lwork,info) 
      IF ( info /= 0 ) CALL juDFT_error('cgetri (inversion of H-zS)') 
                                                                        
      !>                                                                
      DEALLOCATE(work) 
                                                                        
      END SUBROUTINE 
      !>                                                                
      !<--S: gf_inv_hmat(A)                                             
      SUBROUTINE gf_inv_hmat(A) 
!***********************************************************************
!     Subroutine encapsulates the LAPACK calls needed to invert a large 
!                                                                       
!     Daniel Wortmann, Juelich 2001                                     
!***********************************************************************
      IMPLICIT NONE 
#include "cpp_double.h"                                                 
      COMPLEX,INTENT(INOUT)::A(:,:) 
!     LOCALS                                                            
      INTEGER             :: nmat,i,j 
      INTEGER             :: lwork,info,ilaenv 
      INTEGER             :: ipiv(size(A,1)) 
      COMPLEX,ALLOCATABLE :: work(:) 
                                                                        
      EXTERNAL ilaenv 
                                                                        
      nmat = SIZE(A,1) 
      !<-- Before starting determine best size of lapack work-array!    
#ifndef CPP_APC                                                         
      !on PC no ilaenv-call is available...                             
      lwork=ilaenv(1,'cgetri',' ',Nmat,Nmat,-1,-1) 
      lwork=Nmat*lwork 
      IF (lwork<0)                                                   &
     &     CALL juDFT_error                                               &
     &     ('gf_gfcn: failed to determine LAPACK work-size')            
#else                                                                   
      lwork=Nmat 
#endif                                                                  
      ALLOCATE(work(lwork),STAT=info) 
      IF (info/=0) THEN 
         WRITE(6,*)                                                     &
     & 'WARNING, failed to ALLOCATE LAPACK workspace in gf_invertmatrix'
         WRITE(*,*)                                                     &
     &        'WARNING, failed to ALLOCATE LAPACK workspace'            
         WRITE(6,*)                                                     &
     &        'Using default value'                                     
         lwork=Nmat 
         ALLOCATE(work(lwork)) 
      ENDIF 
      !>                                                                

                                                                        
      !<-- Now invert the Matrix!                                       
      WRITE(*,*) 'Hermitian matrix inversion' 
      CALL CPP_LAPACK_chetrf('L',Nmat,A,Nmat,ipiv,work,lwork,info) 
      IF ( info /= 0 ) CALL juDFT_error('chetrf (inversion of H-zS)') 
      CALL CPP_LAPACK_chetri('L',Nmat,A,Nmat,ipiv,work(1:Nmat),info) 
      IF ( info /= 0 ) CALL juDFT_error('chetri (inversion of H-zS)') 
!     chetri returns only the lower diagonal part, so we set the upper  
!     diag...                                                           
      DO i = 1,nmat 
         DO j = i+1,nmat 
            A(i,j) = CONJG(A(j,i)) 
         ENDDO 
      ENDDO 
      !>                                                                
      DEALLOCATE(work) 
                                                                        
      END SUBROUTINE 
      !>                                                                
      !<--S: gf_inv_rmat(A)                                             
      SUBROUTINE gf_inv_rmat(A) 
!***********************************************************************
!     Subroutine encapsules the LAPACK calls needed to invert a large ma
!                                                                       
!     Daniel Wortmann, Juelich 2001                                     
!***********************************************************************
      IMPLICIT NONE 
#include "cpp_double.h"                                                 
!     Arguments                                                         
      REAL,INTENT(INOUT)::A(:,:) 
!     LOCALS                                                            
      INTEGER      :: nmat 
      INTEGER      :: lwork,info,ilaenv 
      INTEGER      :: ipiv(SIZE(A,1)) 
      REAL,ALLOCATABLE::work(:) 
                                                                        
      EXTERNAL ilaenv 
      nmat = SIZE(A,1) 
      !<-- Before starting determine best size of lapack work-array!    
#ifndef CPP_APC                                                         
      !on PC no ilaenv-call is available...                             
      lwork=ilaenv(1,'ssytrf',' ',Nmat,Nmat,-1,-1) 
      lwork=Nmat*lwork 
      IF (lwork<0)                                                   &
     &     CALL juDFT_error                                               &
     &     ('gf_gfcn: failed to determine LAPACK work-size')            
#else                                                                   
      lwork=Nmat 
#endif                                                                  
      ALLOCATE(work(lwork),STAT=info) 
      IF (info/=0) THEN 
         WRITE(6,*)                                                     &
     & 'WARNING, failed to ALLOCATE LAPACK workspace in gf_invertmatrix'
         WRITE(*,*)                                                     &
     &        'WARNING, failed to ALLOCATE LAPACK workspace'            
         WRITE(6,*)                                                     &
     &        'Using default value'                                     
         lwork=Nmat 
         ALLOCATE(work(lwork)) 
      ENDIF 
      !>                                                                
                                                                        
      !<-- Now invert the Matrix!                                       
                                                                        
      WRITE(*,*) 'Non-Hermitian real matrix inversion' 
      CALL CPP_LAPACK_sgetrf(Nmat,Nmat,A,Nmat,ipiv,info) 
      IF ( info /= 0 ) CALL juDFT_error('cgetrf (inversion of H-zS)') 
      CALL CPP_LAPACK_sgetri(Nmat,A,Nmat,ipiv,work,lwork,info) 
      IF ( info /= 0 ) CALL juDFT_error('cgetri (inversion of H-zS)') 
                                                                        
      !>                                                                
                                                                        
      DEALLOCATE(work) 
                                                                        
      END SUBROUTINE 
      !>                                                                
      !<--S: gf_inv_smat(A)                                             
      SUBROUTINE gf_inv_smat(A) 
!***********************************************************************
!     Subroutine encapsules the LAPACK calls needed to invert a large ma
!                                                                       
!     Daniel Wortmann, Juelich 2001                                     
!***********************************************************************
      IMPLICIT NONE 
#include "cpp_double.h"                                                 
!     Arguments                                                         
      REAL,INTENT(INOUT)::A(:,:) 
!     LOCALS                                                            
      INTEGER      :: nmat,i,j 
      INTEGER      :: lwork,info,ilaenv 
      INTEGER      :: ipiv(SIZE(A,1)) 
      REAL,ALLOCATABLE::work(:) 
                                                                        
      EXTERNAL ilaenv 
      nmat = SIZE(A,1) 
      !<-- Before starting determine best size of lapack work-array!    
#ifndef CPP_APC                                                         
      !on PC no ilaenv-call is available...                             
      lwork=ilaenv(1,'ssytrf',' ',Nmat,Nmat,-1,-1) 
      lwork=Nmat*lwork 
      IF (lwork<0)                                                   &
     &     CALL juDFT_error                                               &
     &     ('gf_gfcn: failed to determine LAPACK work-size')            
#else                                                                   
      lwork=Nmat 
#endif                                                                  
      ALLOCATE(work(lwork),STAT=info) 
      IF (info/=0) THEN 
         WRITE(6,*)                                                     &
     & 'WARNING, failed to ALLOCATE LAPACK workspace in gf_invertmatrix'
         WRITE(*,*)                                                     &
     &        'WARNING, failed to ALLOCATE LAPACK workspace'            
         WRITE(6,*)                                                     &
     &        'Using default value'                                     
         lwork=Nmat 
         ALLOCATE(work(lwork)) 
      ENDIF 
      !>                                                                
                                                                        
      !<-- Now inverte the Matrix!                                      
      WRITE(*,*) "Real-symmetric matrix inversion" 
      CALL CPP_LAPACK_ssytrf('L',Nmat,A,Nmat,ipiv,work,lwork,info) 
      IF ( info /= 0 ) CALL juDFT_error('ssytrf (inversion of H-zS)') 
      CALL CPP_LAPACK_ssytri('L',Nmat,A,Nmat,ipiv,work(1:Nmat),info) 
      IF ( info/=0 ) CALL juDFT_error('ssytri (inversion of H-zS)') 
      DO i = 1,nmat 
         DO j = i+1,nmat 
            A(i,j) = A(j,i) 
         ENDDO 
      ENDDO 
      !>                                                                
      DEALLOCATE(work) 
                                                                        
      END SUBROUTINE 
      !>                                                                
                                                                        
                                                                        
      END                                           
