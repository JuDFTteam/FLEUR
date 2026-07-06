!--------------------------------------------------------------------------------
! Copyright (c) 2026 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions 
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------
      MODULE m_gf_get_spectrum 
          IMPLICIT NONE
      CONTAINS 
      SUBROUTINE gf_get_spectrum(                                       &
     &     layer,jspin,gfinp,cell,lapw,                                 &
     &     l_inversion,l_fullgreen,l_nogno,l_nohelpregion,              &
     &     l_writespectrum,nk,l_sph)
!***************************************************
!     Diagonalize the Hamiltonian and calculate                         
!     projections of the eigenvectors on various                        
!     subspaces.                                                        
!     Frank Freimuth, November 2007                                     
!***************************************************                    
      USE m_gf_types 
      USE m_gf_energies,ONLY:gf_z 
      use m_juDFT 
      USE m_gf_spectrum 
      USE m_gf_calcspectrum 
      USE m_gf_calcspectrum2 
!      use m_gf_overlap                                                 
      IMPLICIT NONE 
      INTEGER,INTENT(IN)        :: layer 
      INTEGER,INTENT(IN)        :: jspin 
      TYPE(t_gfinp), INTENT(IN) :: gfinp 
      TYPE(t_cell),INTENT(IN)   :: cell 
      TYPE(t_lapw),INTENT(IN)   :: lapw 
      LOGICAL,INTENT(IN)        :: l_inversion,l_fullgreen,l_nogno 
      LOGICAL,INTENT(IN)        :: l_nohelpregion 
      LOGICAL,INTENT(IN)        :: l_writespectrum,l_sph
      INTEGER,INTENT(IN)        :: nk 
                                                                        
      COMPLEX                   :: cmplxone,z,enedeno 
      COMPLEX,ALLOCATABLE       :: g_temp(:,:),g_temp2(:,:) 
      INTEGER                   :: matsize,j,i,info 
      LOGICAL                   :: l_ldauwan 
      INTEGER                   :: totwann 
      LOGICAL                   :: l_spectrum2 
      CHARACTER(LEN=14)          :: uhuname 
                                                                        
      INTEGER             :: n 
      REAL   ,ALLOCATABLE :: r_n(:),f_n(:) 
                                                                        
      REAL time1,time2 
      matsize = lapw%nmat 
      cmplxone = CMPLX(1.0,0.0) 
      INQUIRE(FILE='spectrum2',EXIST=l_spectrum2) 
                                                                        
      IF(.NOT.ALLOCATED(uhuprojone))THEN 
         IF(l_spectrum2)THEN 
            !CALL gf_calcspectrum2(                                      &
     !&           l_inversion,layer,lapw,                                &
     !&                          jspin)
            CALL gf_calcspectrum_simple(layer,lapw,jspin,l_sph)
         ELSE 
            CALL gf_calcspectrum(layer,lapw,jspin,l_sph)
         ENDIF 
         CALL gf_uhuproj(layer,jspin,lapw,cell,l_nohelpregion,l_sph)
#ifdef CPP_WANNIER                                                      
         IF(gfinp%l_addselfen)THEN 
            CALL gf_wannprojspec() 
         ENDIF 
#endif                                                                  
         IF(l_writespectrum)THEN
            CALL juDFT_error("BUG l_writespectrum not implemted",calledby="gf_get_spectrum")
            WRITE(uhuname,777)layer,nk 
  777       FORMAT('uhu_',i3.3,'_',i6.6) 
            OPEN(222,FILE=uhuname,FORM='unformatted') 
            WRITE(222)lapw%nmat,lapw%nv2_tot 
            WRITE(222)uhuprojone 
            WRITE(222)uhuprojtwo 
            WRITE(222)uhueigval 
            WRITE(222)uhumatrix 
            CLOSE(222) 
         ENDIF 
         !<-- norm of states                                            
!         IF (gfinp%l_overlap) THEN                                     
!            ALLOCATE(f_n(SIZE(uhumatrix,1)))                           
!            ALLOCATE(f_r(SIZE(uhumatrix,1)))                           
!            f_n = gf_overlap_vecnorm(uhumatrix,.TRUE.)                 
!            r_n = gf_overlap_vecnorm(uhumatrix,.FALSE.)                
!             write(6,*) "Eigenvalues of H"                             
!             WRITE(6,*)                                                
!     +            "n    ew   full_norm  restricted_norm"               
!            DO n = 1,SIZE(f_n)                                         
!               WRITE(6,"(i4,1x,5(f15.9,1x))") n,uhueigval(n),f_n(n)    
!     +        ,r_n(n)                                                  
!            ENDDO                                                      
!            DEALLOCATE(f_n,r_n)                                        
!         ENDIF                                                         
         IF(.NOT.l_fullgreen.AND..NOT.l_nogno.AND.                      &
     &             .NOT.gfinp%l_addselfen) DEALLOCATE(uhumatrix)        
            !construct uhumatrix,uhueigval,uhuprojone,uhuprojtwo        
      ENDIF 
                                                                        
      END SUBROUTINE gf_get_spectrum 
      END                                           
