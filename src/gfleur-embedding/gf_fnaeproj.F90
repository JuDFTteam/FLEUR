!--------------------------------------------------------------------------------
! Copyright (c) 2026 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions 
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------
      MODULE m_gf_fnaeproj 
      IMPLICIT NONE
!------------------------------------------------------------           
!     Calculate Surface Projection of Green's function from the         
!     spectral representation of Hamiltonian.                           
!     Optionally, the 3D Green's Function is calculated additionally.   
!     This is controlled by the logical l_fullgreen.                    
!     In the nogno-mode, the surface projections are                    
!     applied only to one side of the Green's function at one           
!     time.                                                             
!     The Green's function calculated here is in any case the           
!     one without embedding potential. An additional embedding          
!     potential may be added later on by the subroutine                 
!     gf_projembed.                                                     
!     Frank Freimuth, January - February 2007                           
!------------------------------------------------------------           
#include "realcomplex.h"                                                
#include "cpp_double.h"
      PRIVATE 
      PUBLIC gf_fnaeproj 
                                                                        
      CONTAINS 
      SUBROUTINE gf_fnaeproj(layer,                                     &
     &     jspin,gfinp,cell,lapw,en,l_sph,                                    &
     &     l_inversion,l_fullgreen,l_nogno,l_nohelpregion,              &
     &     gij,g,real_energy)
!*****************************************************                  
!     Calculate the "free" Green's function, i.e. the                   
!     Green's function without embedding potential.                     
!     gij: both indices from 2D basis                                   
!     g  : 3D Green's function                                          
!     gleft,gright: one index from 2D basis, the other                  
!     from 3D basis                                                     
!     Frank Freimuth, January - February 2007                           
!*****************************************************                  
      USE m_gf_types 
      USE m_gf_energies,ONLY:gf_z 
      use m_juDFT 
      USE m_gf_spectrum 
#include "juDFT_env.h"
      IMPLICIT NONE 
      INTEGER,INTENT(IN)        :: jspin,layer 
      TYPE(t_gfinp), INTENT(IN) :: gfinp 
      TYPE(t_cell),INTENT(IN)   :: cell 
      TYPE(t_lapw),INTENT(IN)   :: lapw 
      INTEGER,INTENT(IN)        :: en 
      LOGICAL,INTENT(IN)        :: l_inversion,l_fullgreen,l_nogno,l_sph
      LOGICAL,INTENT(IN)        :: l_nohelpregion 
      COMPLEX,INTENT(OUT)       :: gij(2*lapw%nv2_tot,2*lapw%nv2_tot) 
      COMPLEX,INTENT(OUT)       :: g(:,:) 
      LOGICAL,INTENT(IN),OPTIONAL:: real_energy
                                                                        
      COMPLEX                   :: cmplxone,z,enedeno 
      COMPLEX,ALLOCATABLE       :: g_temp(:,:),g_temp2(:,:) 
      INTEGER                   :: matsize,matsize1,j,i,info
      LOGICAL                   :: l_ldauwan
      INTEGER                   :: totwann 
      LOGICAL                   :: l_spectrum2 
      REAL                      :: smooth
                                                                        

      matsize = lapw%nmat 
      if (l_sph) matsize=lapw%nmat_sphere

      cmplxone = CMPLX(1.0,0.0) 
                                                                        
!*******************************************************************    
!     Calculate g using the uhuproj-matrices prepared above.            
!     g is the surface projected Green function without embedding.      
!*******************************************************************    
                                                                        
      IF(.NOT.( ALLOCATED(uhuprojone) .AND. ALLOCATED(uhuprojtwo) .AND. &
     &     ALLOCATED(uhueigval) )) CALL juDFT_error("fnae: no uhueigval") 
                                                                        
      z = gf_z(en,layer)
      smooth=0.0
      if (present(real_energy)) THEN
         if (real_energy) THEN
                  smooth=imag(z)
                  z=real(z)
         ENDIF
      ENDIF
                                                                        
      ALLOCATE(g_temp(2*lapw%nv2_tot,matsize),STAT=info) 
      IF(info/=0) CALL juDFT_error("g_temp1") 
      DO j = 1,matsize
         IF (smooth.ne.0.0) THEN
              enedeno = cmplxone/(z-CMPLX(uhueigval(j),smooth))
              enedeno=real(enedeno)
         ELSE
              enedeno = cmplxone/(z-CMPLX(uhueigval(j),0.0))
         ENDIF
         DO i = 1,2*lapw%nv2_tot 
            g_temp(i,j) = uhuprojone(i,j)*enedeno 
               !i                                                       
         ENDDO 
               !j                                                       
      ENDDO 
                                                                        
                                                                        
      CALL CPP_BLAS_cgemm('N','N',                                 &
     &      2*lapw%nv2_tot,2*lapw%nv2_tot,matsize                       &
     &     ,CMPLX(1.0,0.0),g_temp,2*lapw%nv2_tot,uhuprojtwo,matsize     &
     &     ,CMPLX(0.0,0.0),gij,2*lapw%nv2_tot)                          
                                                                        
      DEALLOCATE(g_temp) 
!********************************************************************** 
!     Calculate wgw using the wannproj-matrices prepared above.         
!     wgw is the Wannier-projected non-interacting Green function       
!     without embedding.                                                
!********************************************************************** 
#ifdef CPP_WANNIER                                                      
      IF(gfinp%l_addselfen)THEN 
       IF(.NOT.( ALLOCATED(wannprojone) .AND. ALLOCATED(wannprojtwo)    &
     &    )) CALL juDFT_error("fnae: no wannproj")                        
                                                                        
       z = gf_z(en,layer) 
       totwann=size(wannprojone,1) 
       ALLOCATE(g_temp(totwann,matsize),STAT=info) 
       IF(info/=0) CALL juDFT_error("g_temp1") 
       ALLOCATE(wgw(totwann,totwann),STAT=info) 
       IF(info/=0) CALL juDFT_error("g_temp1") 
       DO j = 1,matsize 
         IF (abs(z-uhueigval(j))<3*abs(smooth)) THEN
              enedeno = cmplxone/(z-CMPLX(uhueigval(j),smooth))
         ELSE
              enedeno = cmplxone/(z-CMPLX(uhueigval(j),0.0))
         ENDIF
         DO i = 1,totwann 
            g_temp(i,j) = wannprojone(i,j)*enedeno 
               !i                                                       
         ENDDO 
                !j                                                      
       ENDDO 
                                                                        
       CALL CPP_BLAS_cgemm('N','N',totwann,totwann,matsize,        &
     &     CMPLX(1.0,0.0),g_temp,totwann,wannprojtwo,matsize,           &
     &     CMPLX(0.0,0.0),wgw,totwann)                                  
                                                                        
       DEALLOCATE(g_temp) 
            !l_addselfen                                                
      ENDIF 
#endif                                                                  
!***********************************************************************
!     Optionally calculate also the 3D Green function without embedding.
!     This is done for l_fullgreen = .true.                             
!***********************************************************************
      IF(l_fullgreen)THEN
         !now we only calculate G up to smaller cutoff nv_sphere
         matsize1=lapw%nmat_sphere

         CPP_juDFT_timestart("gf_fnae_full")
         IF(.NOT.ALLOCATED(uhumatrix)) CALL juDFT_error                   &
     &        ("fullgreen, no uhumat")                                  
         IF(SIZE(g,1)/=matsize1)CALL juDFT_error("size(g,1)")
         IF(SIZE(g,2) /= matsize1)CALL juDFT_error("size(g,2)")
         ALLOCATE(g_temp(matsize1,matsize1),STAT=info)
         ALLOCATE(g_temp2(matsize1,matsize1),STAT=info)
         CPP_juDFT_timestart("g_temp")
         !Copy out part of the uhumatrix

         g_temp2(:lapw%nv_tot_sphere,:)=uhumatrix(:lapw%nv_tot_sphere,:matsize1)   !lapw part of vectors
         g_temp2(lapw%nv_tot_sphere+1:,:)=uhumatrix(lapw%nv_tot_sphere+1:,:matsize1)  !lo-part of vectors

         IF(info/=0) CALL juDFT_error("g_temp2")

         DO j = 1,matsize1
            IF (abs(z-uhueigval(j))<3*abs(smooth)) THEN
              enedeno = cmplxone/(z-CMPLX(uhueigval(j),smooth))
            ELSE
              enedeno = cmplxone/(z-CMPLX(uhueigval(j),0.0))
            ENDIF
            DO i = 1,matsize1
               g_temp(i,j) = g_temp2(i,j)*enedeno
                  !i                                                    
            ENDDO 
                  !j                                                    
         ENDDO 
         CPP_juDFT_timestop("g_temp")
         CPP_juDFT_timestart("dgemm")
         CALL CPP_BLAS_cgemm('N','C',matsize1,matsize1,matsize1,      &
     &        CMPLX(1.0,0.0),g_temp,matsize1,g_temp2,                  &
     &        matsize1,CMPLX(0.0,0.0),g,matsize1)
         CPP_juDFT_timestop("dgemm")
         DEALLOCATE(g_temp,g_temp2)
         CPP_juDFT_timestop("gf_fnae_full")
             !l_fullgreen                                               
      ENDIF 
                                                                        
!***********************************************************************
!     Optionally calculate the (reduced) Green's functions with one 3D-i
!     and one 2D-index.                                                 
!     gright: first index 3D and second index 2D                        
!     gleft : first index 2D and second index 3D                        
!***********************************************************************
      IF(l_nogno.OR.gfinp%l_addselfen)THEN
         CPP_juDFT_timestart("gf_fnae_red")
         g(:,:)=cmplx(0.0,0.0) 
         IF(allocated(gright))DEALLOCATE(gright)
         IF(allocated(gleft))DEALLOCATE(gleft)
         ALLOCATE(gright(matsize,2*lapw%nv2_tot))
         ALLOCATE(gleft(2*lapw%nv2_tot,matsize))
                                                                        
         ALLOCATE(g_temp(2*lapw%nv2_tot,matsize)) 
         ALLOCATE(g_temp2(2*lapw%nv2_tot,matsize)) 
                                                                        
         g_temp2=conjg(uhuprojone) 
         DO j=1,matsize 
           IF (abs(z-uhueigval(j))<3*abs(smooth)) THEN
              enedeno = cmplxone/(z-CMPLX(uhueigval(j),0.0))
           ELSE
              enedeno = cmplxone/(z-CMPLX(uhueigval(j),0.0))
           ENDIF

          DO i=1,2*lapw%nv2_tot 
           g_temp(i,j)=g_temp2(i,j)*enedeno 
                !i                                                      
          ENDDO 
          DO i=1,2*lapw%nv2_tot 
           gleft(i,j)=uhuprojone(i,j)*enedeno 
                !i                                                      
          ENDDO 
               !j                                                       
         ENDDO 
                                                                        
         gright=transpose(g_temp) 
         DEALLOCATE(g_temp2) 
         DEALLOCATE(g_temp) 
         CPP_juDFT_timestop("gf_fnae_red")
      ENDIF 
                                                                        
      END SUBROUTINE gf_fnaeproj 
      END                                           
