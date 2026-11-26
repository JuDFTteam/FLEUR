!--------------------------------------------------------------------------------
! Copyright (c) 2026 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions 
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------
      MODULE m_gf_calcspectrum2 
          IMPLICIT NONE
      CONTAINS 
      SUBROUTINE gf_calcspectrum2(l_inversion,layer,lapw                &
     &     ,jspin)                                                      
!-----------------------------------------------                        
! solve generalized eigenvalue problem for spectral representation      
! the hamiltonian&overlapp are taken from gf_hsdata                     
! The vector of eigenvalues: uhueigval.                                 
! The matrix of eigenvectors: uhumatrix.                                
! The left surface projected uhumatrix: uhuprojone.                     
! The right surface projected uhumatrix: uhuprojtwo.                    
! Frank Freimuth, January-February 2007                                 
!           (last modified: 07-02-22) D. Wortmann                       
!-----------------------------------------------                        
      USE m_gf_hsdata 
      USE m_gf_types 
      use m_juDFT
#ifdef CPP_WANNIER                                                      
      USE m_gf_selfenergy 
#endif                                                                  
      USE m_gf_spectrum 
      USE m_gf_mkposdef 
      IMPLICIT NONE 
      !<-- Arguments                                                    
                                                                        
      INTEGER,INTENT(IN)        :: jspin,layer 
      TYPE(t_lapw),INTENT(IN)   :: lapw 
      LOGICAL,INTENT(IN)        :: l_inversion 
                                                                        
      !>                                                                
      !<-- Locals                                                       
      INTEGER                  :: matsize 
      LOGICAL                  :: l_firstrun,l_posdef 
      REAL,ALLOCATABLE         :: over_REAL(:,:) 
      REAL,ALLOCATABLE         :: uhumatrix_REAL(:,:) 
      COMPLEX,ALLOCATABLE      :: over(:) 
      COMPLEX,ALLOCATABLE      :: ham(:) 
      INTEGER                  :: in,i,j,n1,n2,info 
      COMPLEX,POINTER :: hp(:,:),sp(:,:) 
      COMPLEX,ALLOCATABLE      :: work(:) 
      REAL,ALLOCATABLE         :: rwork(:) 
      INTEGER,ALLOCATABLE      :: iwork(:),ifail(:) 
      INTEGER                  :: lwork,lrwork,liwork 
      COMPLEX                  :: exp1,f1 
      LOGICAL                  :: l_ldauwan 
      REAL                     :: lb,ub 
      INTEGER                  :: iu,ne,ind1 
      REAL                     :: toler,abstol 
      !>                                                                
                                                                        
      INQUIRE(FILE='ldauwan',EXIST=l_ldauwan) 
      matsize = lapw%nmat 
      IF(ALLOCATED(uhuprojtwo)) CALL juDFT_error("uhu, no uhu") 
      IF(ALLOCATED(uhueigval))  CALL juDFT_error("uhueigval, no uhu") 
      IF(ALLOCATED(uhumatrix))  CALL juDFT_error("uhumatrix, no uhu") 
      CALL gf_getHSaddr(layer,hp,sp) 
!***************************************************************        
!        Solve the generalized eigenvalue problem.                      
!***************************************************************        
!if posdef if present, the eigenvalues of the overlap-matrix            
!are calculated and those which are too small are made bigger           
      abstol=2.0*TINY(abstol) 
      INQUIRE(FILE ='posdef',EXIST = l_posdef) 
      ALLOCATE(over((matsize*(matsize+1))/2)) 
      CALL priv_pack_upper(sp,matsize,over) 
      info=1 
      IF (l_posdef) GOTO 200
      CALL timestart("Diagonalization")
      CALL zpptrf('U',matsize,over,info) 
      IF(info/=0)CALL priv_pack_upper(sp,matsize,over) 
      CALL timestop("Diagonalization")
  200 CONTINUE 
      IF(info/=0)THEN 
         CALL timestart("Make Posdef")
         WRITE(*,*) "make metrix positive definite" 
         ALLOCATE(uhumatrix(matsize,matsize)) 
         ALLOCATE(uhueigval(matsize)) 
         ALLOCATE ( work(2*matsize),iwork(5*matsize),ifail(matsize) ) 
         ALLOCATE( rwork(8*matsize)) 
                         !eigenvalues below toler will be shifted       
         toler = 1.0E-12 
         lb =-1.0E12 
         ub = toler 
         CALL zhpevx('V','V','U',matsize,over,lb,ub,1,iu,    &
     &             abstol,                                              &
     &             ne,uhueigval,uhumatrix,matsize,                      &
     &             work,rwork,iwork,ifail,info)                         
         IF(info/=0) CALL juDFT_error('chpevx-posdef') 
         DEALLOCATE(rwork,work,iwork,ifail) 
         CALL priv_pack_upper(sp,matsize,over) 
  !remove problematic eigenvalues                                       
         DO ind1 = 1,ne 
            uhueigval(ind1) = toler-uhueigval(ind1) 
         ENDDO 
         DO ind1 = 1,ne
            CALL zhpr('U',matsize,uhueigval(ind1),        &
     &           uhumatrix(1,ind1),1,over)                              
         ENDDO 
         CALL zpptrf('U',matsize,over,info) 
         IF(info/=0)CALL juDFT_error('cpptrf-second time')
         CALL timestop("Make Posdef")
      ENDIF 
      CALL timestart("Diagonalization")
      !DEALLOCATE(sp)
      ALLOCATE(ham((matsize*(matsize+1))/2)) 
      CALL priv_pack_upper(hp,matsize,ham) 
      !DEALLOCATE(hp)
      CALL zhpgst(1,'U',matsize,ham,over,info) 
      IF(info/=0)CALL juDFT_error('chpgst') 
                                                                        
      IF(.NOT.allocated(uhumatrix))ALLOCATE(uhumatrix(matsize,matsize)) 
      IF(.NOT.allocated(uhueigval))ALLOCATE(uhueigval(matsize)) 
      ALLOCATE ( work(2*matsize),iwork(5*matsize),ifail(matsize) ) 
      ALLOCATE( rwork(8*matsize)) 
      CALL zhpevx('V','A','U',matsize,ham,lb,ub,1,iu,abstol, &
     &             ne,uhueigval,uhumatrix,matsize,                      &
     &             work,rwork,iwork,ifail,info)                         
      IF(info/=0) CALL juDFT_error('chpevx') 
      DEALLOCATE(rwork,work,iwork,ifail) 
      DEALLOCATE(ham) 
                                                                        
      CALL ztptrs('U','N','N',matsize,matsize,               &
     &                   over,uhumatrix,matsize,info)                   
      IF(info/=0) CALL juDFT_error('ctptrs') 
      DEALLOCATE(over) 
      CALL timestop("Diagonalization")
      END SUBROUTINE 
      SUBROUTINE priv_pack_upper(full,n,packed)
      !pack the upper triangle of the FULL Hermitian matrix (standard
      !convention, see gf_hsdata) into LAPACK 'U' packed storage
      IMPLICIT NONE
      COMPLEX,INTENT(IN)  :: full(:,:)
      INTEGER,INTENT(IN)  :: n
      COMPLEX,INTENT(OUT) :: packed(:)
      INTEGER :: i,j,in
      in=0
      DO j=1,n
         DO i=1,j
            in=in+1
            packed(in)=full(i,j)
         ENDDO
      ENDDO
      END SUBROUTINE priv_pack_upper

      END                                           
