!------------------------------------------------                       
! Copyright (c) 2026 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions 
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------
!     wannprojone,wannprojtwo -> projection onto Wannier functions      
!     uhueigval -> vector of eigenvalues                                
!     Frank Freimuth, November 2007                                     
!------------------------------------------------                       
      MODULE m_gf_spectrum 
      use m_juDFT
      IMPLICIT NONE
      REAL,   ALLOCATABLE,SAVE :: uhueigval(:) 
      COMPLEX,ALLOCATABLE,SAVE :: uhuprojone(:,:) 
      COMPLEX,ALLOCATABLE,SAVE :: uhuprojtwo(:,:) 
      COMPLEX,ALLOCATABLE,SAVE :: uhumatrix(:,:) 
      COMPLEX,ALLOCATABLE,SAVE :: gright(:,:) 
      COMPLEX,ALLOCATABLE,SAVE :: gleft(:,:) 
#ifdef CPP_WANNIER                                                      
      COMPLEX,ALLOCATABLE,SAVE :: wgw(:,:) 
      COMPLEX,ALLOCATABLE,SAVE :: wannprojone(:,:) 
      COMPLEX,ALLOCATABLE,SAVE :: wannprojtwo(:,:) 
#endif                                                                  
      CONTAINS 
      !<-- S: gf_spectrum_clean()                                       
      SUBROUTINE  gf_spectrum_clean() 
!-----------------------------------------------                        
! deallocate module data                                                
!           (last modified: 2004-00-00) D. Wortmann                     
!-----------------------------------------------                        
      IMPLICIT NONE 
      IF (ALLOCATED(uhueigval))  DEALLOCATE(uhueigval) 
      IF (ALLOCATED(uhuprojone)) DEALLOCATE(uhuprojone) 
      IF (ALLOCATED(uhuprojtwo)) DEALLOCATE(uhuprojtwo) 
      IF (ALLOCATED(uhumatrix))  DEALLOCATE(uhumatrix) 
#ifdef CPP_WANNIER                                                      
      IF (ALLOCATED(wannprojone))DEALLOCATE(wannprojone) 
      IF (ALLOCATED(wannprojtwo))DEALLOCATE(wannprojtwo) 
#endif                                                                  
      END SUBROUTINE 
      !>                                                                
                                                                        
      SUBROUTINE gf_uhuproj(layer,jspin,lapw,lapw_gf,cell,l_nohelpregion,l_sph)
!****************************************************************       
!        Construct the surface projections of eigenvector-matrices.     
!        The surface projections are called uhuprojone and              
!        uhuprojtwo.                                                    
!        Frank Freimuth, November 2007                                  
!****************************************************************       
      USE m_gf_types 
      USE m_gf_curvy2dprojector 
      IMPLICIT NONE 
      INTEGER, INTENT(IN)       :: layer 
      INTEGER, INTENT(IN)       :: jspin 
      TYPE(t_lapw),INTENT(IN)   :: lapw
      TYPE(t_lapw_gf),INTENT(IN):: lapw_gf
      TYPE(t_cell),INTENT(IN)   :: cell 
      LOGICAL, INTENT(IN)       :: l_nohelpregion,l_sph
                                                                        
      COMPLEX                   :: exp1,f1 
      INTEGER matsize 
      INTEGER n1,n2 
      COMPLEX,ALLOCATABLE       :: uhuprojone_tmp(:,:) 
                                                                        
      if (l_sph) then
        matsize=lapw_gf%nmat_sphere
      else
        matsize=lapw%nmat
      endif
                                                                        
!      open(777,file='uhumatrix')                                       
!      write(777,*)uhumatrix(:,:)                                       
!      close(777)                                                       
                                                                        
      ALLOCATE( uhuprojone(2*lapw_gf%nv2_tot,matsize)) 
      IF(.NOT.l_nohelpregion)THEN 
       uhuprojone=cmplx(0.0,0.0) 
       exp1=exp(cmplx(0.0,-1.0*cell%bmat(3,3)*cell%z1)) 
       DO n2=1,matsize 
         DO n1=1,lapw%nv_tot 
          f1=1/sqrt(cell%amat(3,3))*exp1**lapw%k3(n1,jspin) 
          uhuprojone(lapw%kp(n1,jspin),n2)=                           &
     &     uhuprojone(lapw%kp(n1,jspin),n2)+f1*uhumatrix(n1,n2)       
               !n1                                                      
         ENDDO 
             !n2                                                        
       ENDDO 
       exp1=exp(cmplx(0.0,1.0*cell%bmat(3,3)*cell%z1)) 
       DO n2=1,matsize 
        DO n1=1,lapw%nv_tot 
         f1=1/sqrt(cell%amat(3,3))*exp1**lapw%k3(n1,jspin) 
         uhuprojone(lapw_gf%nv2_tot+lapw%kp(n1,jspin),n2)=               &
     &       uhuprojone(lapw_gf%nv2_tot+lapw%kp(n1,jspin),n2)+           &
     &       f1*uhumatrix(n1,n2)                                        
              !n1                                                       
        ENDDO 
             !n2                                                        
       ENDDO 
      ELSE 
                                                                        
         CALL gf_curvy2dprojector(layer,cell,lapw,lapw_gf,.TRUE.) 
                                                                        
!         CALL gf_curvy2dproject(lapw,uhumatrix,uhuprojone)             
        ALLOCATE( uhuprojone_tmp(2*lapw_gf%nv2_tot,matsize)) 
        CALL gf_get_curvy2dproject(lapw,lapw_gf,uhuprojone_tmp) 
!        uhuprojone_tmp(1:lapw_gf%nv2_tot,1:matsize) =                     
!     =              curvyproject(1:lapw_gf%nv2_tot,1:matsize,1)           
!        uhuprojone_tmp(lapw_gf%nv2_tot+1:2*lapw_gf%nv2_tot,1:matsize) =      
!     =              curvyproject(1:lapw_gf%nv2_tot,1:matsize,2)           
!        uhuprojone = matmul(uhuprojone_tmp,uhumatrix)                  
        CALL zgemm('N','N',2*lapw_gf%nv2_tot,matsize,matsize,              &
     &             cmplx(1.0,0.0),uhuprojone_tmp,2*lapw_gf%nv2_tot,        &
     &             uhumatrix,matsize,                                   &
     &             cmplx(0.0,0.0),uhuprojone,2*lapw_gf%nv2_tot)            
        DEALLOCATE( uhuprojone_tmp) 
      ENDIF 
      ALLOCATE(uhuprojtwo(matsize,2*lapw_gf%nv2_tot)) 
      uhuprojtwo=transpose(conjg(uhuprojone)) 
      END SUBROUTINE gf_uhuproj 
#ifdef CPP_WANNIER                                                      
      SUBROUTINE gf_wannprojspec() 
!************************************************************           
!     Calculate the projections of eigenvector matrix onto              
!     the Wannier function subspace.                                    
!     Frank Freimuth, November 2007                                     
!************************************************************           
      USE m_gf_selfenergy,ONLY:projwan 
      IMPLICIT NONE 
      INTEGER totwann 
      INTEGER matsize 
                                                                        
      IF(.NOT.allocated(projwan)) CALL juDFT_error("projwan",calledby="gf_spectrum.F90")
      matsize=size(projwan,1) 
      totwann=size(projwan,2) 
      ALLOCATE( wannprojone(totwann,matsize)) 
      CALL zgemm('C','N',totwann,matsize,matsize,                       &
     &             cmplx(1.0,0.0),projwan,matsize,                      &
     &             uhumatrix,matsize,                                   &
     &             cmplx(0.0,0.0),wannprojone,totwann)                  
      ALLOCATE( wannprojtwo(matsize,totwann)) 
      wannprojtwo=transpose(conjg(wannprojone)) 
      END SUBROUTINE gf_wannprojspec 
#endif                                                                  
      SUBROUTINE gf_read_spectrum(layer,nk)
      IMPLICIT NONE
!************************************************************           
!     Read uhumatrix, uhuprojone and uhuprojtwo from files.             
!     Frank Freimuth, November 2007                                     
!************************************************************           
      INTEGER, INTENT(IN) :: layer,nk 
                                                                        
      INTEGER nmat,nv2 
      CHARACTER(LEN=14)    :: uhuname 
                                                                        
      CALL gf_spectrum_clean() 
      WRITE(uhuname,777)layer,nk 
  777 FORMAT('uhu_',i3.3,"_",i6.6) 
      OPEN(222,FILE=uhuname,FORM='unformatted') 
      READ(222)nmat,nv2 
      ALLOCATE( uhuprojone(2*nv2,nmat) ) 
      ALLOCATE( uhuprojtwo(nmat,2*nv2) ) 
      ALLOCATE( uhumatrix(nmat,nmat) ) 
      ALLOCATE( uhueigval(nmat) ) 
      READ(222) uhuprojone 
      READ(222) uhuprojtwo 
      READ(222) uhueigval 
      READ(222) uhumatrix 
      CLOSE(222) 
      END SUBROUTINE gf_read_spectrum
                                                                        
      END                                           
