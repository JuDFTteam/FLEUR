!--------------------------------------------------------------------------------
! Copyright (c) 2026 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions 
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------
      MODULE m_gf_stepsanaly 
      use m_juDFT
      IMPLICIT NONE
                                                                        
      PRIVATE 
      INTEGER,SAVE   :: s_DIM(3),s_mx3,s_mx2,s_mx1 
      REAL,PARAMETER :: cutoff = 2.1 
      PUBLIC         :: gf_calcstepsanaly,gf_gspaceconvolve 
      PUBLIC         :: gf_stepf_nohelpregion,gf_initstepsanaly 
      PUBLIC         :: gf_steps_pot_convolve 
      CONTAINS 
                                                                        
      SUBROUTINE gf_stepf_nohelpregion(layer,                  &
     &     mx1,mx2,napw,                                                &
     &     ustep)                                                       
!******************************************************                 
!     Generate the stepfunction for the kinetic energy.                 
!     Needed for the version without help region.                       
!     Frank Freimuth, October 2007                                      
!******************************************************                 
      IMPLICIT NONE 
      INTEGER,INTENT(IN)  :: layer
      INTEGER,INTENT(IN)  :: mx1,mx2,napw 
      COMPLEX,INTENT(OUT) :: ustep(-mx1:,-mx2:,-2*napw:) 
                                                                        
      COMPLEX,ALLOCATABLE :: stepkomplett(:,:,:) 
      INTEGER             :: k1,k2,k3 
                                                                        
!Obtain the stepfunctions                                               
      ALLOCATE(stepkomplett(-s_mx1:s_mx1,-s_mx2:s_mx2,-s_mx3:s_mx3)) 
      CALL priv_loadstep(layer,stepkomplett = stepkomplett)
                                                                        
      DO k3 =-2*napw,2*napw 
         DO k2 =-mx2,mx2 
            DO k1 =-mx1,mx1 
               ustep(k1,k2,k3) = stepkomplett(k1,k2,k3) 
            ENDDO 
         ENDDO 
      ENDDO 
      DEALLOCATE(stepkomplett) 
                                                                        
      END SUBROUTINE gf_stepf_nohelpregion 
                                                                        
      SUBROUTINE gf_initstepsanaly(stars,napw) 
!***************************************************************        
!     Calculate the correct dimensioning for the step functions.        
!     Frank Freimuth, November 2007                                     
!***************************************************************        
      USE m_gf_types 
      use m_juDFT 
      IMPLICIT NONE 
      !<-- Arguments                                                    
      TYPE(t_stars),INTENT(IN) :: stars 
      INTEGER,INTENT(IN)       :: napw 
      INTEGER                  :: mx3 
                                                                        
      !<-- Check if 3rd dimension is sufficient for zylinder basis      
      IF (2*napw>cutoff*stars%mx3) THEN 
         WRITE(*,*) "Generating Step-functions with larger z-grid" 
         WRITE(*,*) INT(cutoff*stars%mx3)," -> ",2*napw 
         mx3 = FLOOR(2*napw/cutoff)+2 
      ELSE 
         mx3 = stars%mx3 
      ENDIF 
      !>                                                                
                                                                        
!Dimensions of step-function arrays.                                    
      s_mx1 = INT(cutoff*stars%mx1) 
      s_mx2 = INT(cutoff*stars%mx2) 
      s_mx3 = INT(cutoff*mx3) 
      s_DIM(1) = 2*s_mx1+1 
      s_DIM(2) = 2*s_mx2+1 
      s_DIM(3) = 2*s_mx3+1 
                                                                        
                                                                        
      END SUBROUTINE gf_initstepsanaly 
                                                                        
      SUBROUTINE gf_calcstepsanaly(layer,cell,stars,napw)
!********************************************************************   
!     Generate the stepfunctions needed in GFLEUR                       
!     (last modified: 06-10-09) D. Wortmann                             
!********************************************************************   
!     Modifications to perform the convolutions analytically            
!     without need of FFT and to                                        
!     calculate the stepfunctions analytically.                         
!     Frank Freimuth, April 2007                                        
!********************************************************************   
      USE hdf5 
      USE m_hdf_tools 
      USE m_gf_types 
      use m_juDFT 
      IMPLICIT NONE 
      !<--Arguments                                                     
      INTEGER,INTENT(IN)       :: layer
      TYPE(t_stars),INTENT(IN) :: stars 
      TYPE(t_cell),INTENT(IN)  :: cell 
      INTEGER,INTENT(IN)       :: napw 
      !>                                                                
      !<-- Locals                                                       
                                                                        
      INTEGER             :: ncutat,newat,nr,nl,n 
      REAL, ALLOCATABLE   :: cutpos(:,:),newpos(:,:),newrmt(:) 
      REAL, ALLOCATABLE   :: cutrmt(:),posL(:,:),posR(:,:),rmtL(:) 
      REAL, ALLOCATABLE   :: rmtR(:) 
      REAL                :: c,time1,time2,time3 
      LOGICAL             :: success 
      COMPLEX,ALLOCATABLE :: step_tmp(:,:,:) 
                                                                        
!Step function for potential (and kinetic energy in "nohelpregion"-mode)
      COMPLEX,ALLOCATABLE :: stepkomplett(:,:,:) 
!Step function for physical volume                                      
      COMPLEX,ALLOCATABLE :: stepphys(:,:,:) 
!Step functions for help regions (to apply potential shifts)            
                                                  !right                
      COMPLEX,ALLOCATABLE :: stephelpreg_r(:,:,:) 
                                                  !left                 
      COMPLEX,ALLOCATABLE :: stephelpreg_l(:,:,:) 
                                                                        
      INTEGER(HID_T)      :: fid,gid,varid,ffid 
      INTEGER             :: hdferr,mx3 
      !>                                                                
                                                                        
                                                                        
      !>                                                                
                                                                        
      !<-- Check if step functions are present                          
      CALL priv_loadstep(layer,success = success)
      IF (success) RETURN 
                                                                        
      WRITE(*,*) "Generating Step-function:",layer 
      !>                                                                
      !<-- read all information needed                                  
      CALL io_hdfopen('gf_setup.hdf',H5F_ACC_RDONLY_F,ffid,hdferr) 
      CALL io_gopen(ffid,io_layername(layer),fid,hdferr) 
      CALL io_gopen (fid,'stepinfo',gid,hdferr) 
      IF (hdferr <0 )  CALL juDFT_error                                   &
     &     ('Information about step functions not found')               
                                                                        
                                                                        
      CALL io_READ_att(gid, "c",c) 
      CALL io_READ_att(gid,"ncutat",ncutat) 
      CALL io_READ_att(gid,"nnewat",newat) 
      ALLOCATE(newpos(3,newat),cutpos(3,ncutat)) 
      ALLOCATE(newrmt(newat),cutrmt(ncutat)) 
      CALL io_READ_att(gid,"newpos",newpos) 
      CALL io_READ_att(gid,"newrmt",newrmt) 
      IF (ncutat>0) THEN 
         CALL io_READ_att(gid,"cutpos",cutpos) 
         CALL io_READ_att(gid,"cutrmt",cutrmt) 
      ENDIF 
      CALL io_gclose (gid, hdferr) 
      CALL io_gclose (fid, hdferr) 
      CALL io_hdfclose(ffid,hdferr) 
      !>                                                                
                                                                        
      !<-- generate the step-functions                                  
      ALLOCATE(stepkomplett(-s_mx1:s_mx1,-s_mx2:s_mx2,-s_mx3:s_mx3)) 
      ALLOCATE(stepphys(-s_mx1:s_mx1,-s_mx2:s_mx2,-s_mx3:s_mx3)) 
                                                                        
      stepkomplett=cmplx(0.0,0.0) 
      stepkomplett(0,0,0)=cmplx(1.0,0.0) 
      ALLOCATE(stephelpreg_r(-s_mx1:s_mx1,-s_mx2:s_mx2,-s_mx3:s_mx3)) 
      ALLOCATE(stephelpreg_l(-s_mx1:s_mx1,-s_mx2:s_mx2,-s_mx3:s_mx3)) 
      ALLOCATE(step_tmp(-s_mx1:s_mx1,-s_mx2:s_mx2,-s_mx3:s_mx3)) 
                                                                        
                                                                        
!*************************************************                      
! Step functions for interstitial regions.                              
!*************************************************                      
      CALL cpu_time(time1) 
      CALL priv_generate_accastep(cell%bmat(3,3),                       &
     &  cell%amat(3,3),c,step_tmp)                                      
      stepkomplett=stepkomplett-step_tmp 
      CALL priv_generate_cdstepright(cell%bmat(3,3),                    &
     &  cell%amat(3,3),c,cell%z1,step_tmp)                              
      stephelpreg_r=step_tmp 
      stephelpreg_l=conjg(step_tmp) 
      CALL cpu_time(time2) 
      time2=time2-time1 
                                                                        
!****************************************************************       
! Some of the MT-spheres in [-c/2,c/2] might cut the acca volume.       
! The common part of the volumes would be double counted then.          
! Calculate the "double counting" here to subtract it again later.      
!****************************************************************       
      CALL cpu_time(time1) 
      CALL priv_generate_ballpartstep_right(cell%omtil,cell%bmat,c,     &
     &                   cell%amat(3,3),newpos,newrmt,step_tmp)         
      stepkomplett=stepkomplett+step_tmp 
      stephelpreg_r=stephelpreg_r-step_tmp 
      CALL priv_generate_ballpartstep_left(cell%omtil,cell%bmat,c,      &
     &                   cell%amat(3,3),newpos,newrmt,step_tmp)         
      stepkomplett=stepkomplett+step_tmp 
      stephelpreg_l=stephelpreg_l-step_tmp 
      CALL cpu_time(time2) 
      time2=time2-time1 
                                                                        
!**************************************************************         
! Step functions due to atoms outside the embedding region.             
!**************************************************************         
      !<-- Sort cutatoms according to the sides                         
      IF (ncutat>0) THEN 
         nr = COUNT(cutpos(3,:)>0) 
         nl = COUNT(cutpos(3,:)<0) 
         ALLOCATE(posL(3,nL),rmtL(nl),posR(3,nr),rmtR(nR)) 
         DO n = 1,ncutat 
            IF (cutpos(3,n)<0) THEN 
               posL(:,nl) = cutpos(:,n) 
               rmtL(nl) = cutrmt(n) 
               nl = nl-1 
            ELSEIF (cutpos(3,n)>0) THEN 
               posR(:,nR) = cutpos(:,n) 
               rmtR(nR) = cutrmt(n) 
               nR = nR-1 
            ELSE 
               CALL juDFT_error("Atom cut at z = 0 from outside?") 
            ENDIF 
         ENDDO 
         IF (ABS(nl)+ABS(nr) /= 0) CALL juDFT_error(                      &
     &        "Misscount ingf_stepfunction")                            
         CALL cpu_time(time1) 
         CALL priv_generate_ballpartstep_right(cell%omtil,cell%bmat,c,  &
     &                  cell%amat(3,3),posR,rmtR,step_tmp)              
         stepkomplett=stepkomplett-step_tmp 
         stephelpreg_r=stephelpreg_r+step_tmp 
         CALL priv_generate_ballpartstep_left(cell%omtil,cell%bmat,c,   &
     &                  cell%amat(3,3),posL,rmtL,step_tmp)              
         stepkomplett=stepkomplett-step_tmp 
         stephelpreg_l=stephelpreg_l+step_tmp 
         CALL cpu_time(time2) 
         time2=time2-time1 
                                                                        
      ENDIF 
      !>                                                                
                                                                        
!***************************************************************        
! Step function for MT-spheres which are in [-c/2,c/2].                 
! Full spheres are considered here, irrespective of whether             
! they cut the acca volume.                                             
!***************************************************************        
      CALL cpu_time(time1) 
      CALL priv_generate_MTstep(cell%omtil,cell%bmat,newpos,            &
     &                                        newrmt,step_tmp)          
      CALL cpu_time(time2) 
      time2=time2-time1 
      stepphys = stepkomplett 
      stepkomplett=stepkomplett-step_tmp 
                                                                        
                                                                        
      CALL priv_savesteps(layer,stepkomplett=stepkomplett,              &
     &     stephelpreg_r = stephelpreg_r,                               &
     &     stephelpreg_l = stephelpreg_l,stepphys=stepphys)             
      DEALLOCATE(stepkomplett,stephelpreg_r,stephelpreg_l,stepphys) 
      DEALLOCATE(step_tmp) 
                                                                        
                                                                        
                                                                        
      IF (ncutat>0) DEALLOCATE(posR,rmtR,rmtL,posL) 
      DEALLOCATE(newpos,newrmt,cutpos,cutrmt) 
      END SUBROUTINE 
                                                                        
      !<--S: gf_gspaceconvolve(layer,stars,pot_aux,vpw,vpw_w)
      SUBROUTINE gf_gspaceconvolve(layer,stars,pot_aux,vpw     &
     &     ,vpw_w)                                                      
!********************************************************************   
!     Calculate the convolution in g-space without making use of fft.   
!     Frank Freimuth, April 2007                                        
!********************************************************************   
      USE m_gf_types 
      USE hdf5 
      USE m_hdf_tools 
      USE m_gf_fft_singleton 
      IMPLICIT NONE 
      !<--Arguments                                                     
      INTEGER,INTENT(IN)       :: layer
      TYPE(t_stars),INTENT(IN) :: stars 
      REAL   ,INTENT(IN)       :: pot_aux 
                                            !Planewave potential        
      COMPLEX,INTENT(IN)       :: vpw(:) 
                                            !Warped planewave potential 
      COMPLEX,INTENT(OUT)      :: vpw_w(:) 
      !>                                                                
      !<-- Locals                                                       
      COMPLEX,ALLOCATABLE :: v(:,:,:) 
                                                                        
      COMPLEX,ALLOCATABLE :: stephelpreg_r(:,:,:) 
      COMPLEX,ALLOCATABLE :: stepkomplett(:,:,:) 
      COMPLEX,ALLOCATABLE :: fftstepkomplett(:,:,:) 
      COMPLEX,ALLOCATABLE :: fftstephelp(:,:,:) 
      COMPLEX,ALLOCATABLE :: kstep(:,:,:) 
      COMPLEX,ALLOCATABLE :: vpw_w_temp(:,:,:) 
      INTEGER             :: i,k1,k2,k3,in,g1,g2,g3,h1,h2,h3 
      COMPLEX  potential 
      REAL time1,time2 
                                                                        
      !>                                                                
                                                                        
!Obtain the stepfunctions                                               
      ALLOCATE(stephelpreg_r(-s_mx1:s_mx1,-s_mx2:s_mx2,-s_mx3:s_mx3)) 
      ALLOCATE(stepkomplett(-s_mx1:s_mx1,-s_mx2:s_mx2,-s_mx3:s_mx3)) 
      CALL priv_loadstep(layer,stepkomplett = stepkomplett,    &
     &     stephelpreg_r = stephelpreg_r)                               
                                                                        
                                                                        
                                                                        
      CALL cpu_time(time1) 
                                                                        
!************************************************                       
! straigth forward convolution: takes much longer                       
!************************************************                       
      IF(.FALSE.)THEN 
                                                                        
!Expand the stars into planewaves.                                      
      ALLOCATE(V(-stars%mx1:stars%mx1,-stars%mx2:stars%mx2,             &
     &              -stars%mx3:stars%mx3))                              
      v = 0 
      DO k1 = -stars%mx1,stars%mx1 
       DO k2 = -stars%mx2,stars%mx2 
        DO k3 = -stars%mx3,stars%mx3 
         in = stars%ig(k1,k2,k3) 
         IF (in == 0) CYCLE 
         V(k1,k2,k3) = vpw(in) 
        ENDDO 
       ENDDO 
      ENDDO 
                                                                        
                                                                        
!Perform the convolution                                                
      ALLOCATE(vpw_w_temp(-stars%mx1:stars%mx1,-stars%mx2:stars%mx2,    &
     &              -stars%mx3:stars%mx3))                              
      vpw_w_temp=cmplx(0.0,0.0) 
                                                                        
                                                                        
      DO k3=-stars%mx3,stars%mx3 
       DO g3=-stars%mx3,stars%mx3 
        h3=g3-k3 
        IF(ABS(h3)>s_mx3)CYCLE 
        DO k2=-stars%mx2,stars%mx2 
         DO g2=-stars%mx2,stars%mx2 
          h2=g2-k2 
          IF(abs(h2)>s_mx2)CYCLE 
          DO k1=-stars%mx1,stars%mx1 
           DO g1=-stars%mx1,stars%mx1 
            h1=g1-k1 
            IF(abs(h1)>s_mx1)CYCLE 
!            if(stars%ig(k1,k2,k3)==0)cycle                             
!            if(stars%ig(g1,g2,g3)==0)cycle                             
            vpw_w_temp(g1,g2,g3)=vpw_w_temp(g1,g2,g3)                   &
     &             +v(k1,k2,k3)*stepkomplett(h1,h2,h3)                  
           ENDDO 
          ENDDO 
         ENDDO 
        ENDDO 
       ENDDO 
      ENDDO 
                                                                        
!Transform the planewaves back to stars                                 
      vpw_w=cmplx(0.0,0.0) 
      DO k1 = -stars%mx1,stars%mx1 
         DO k2 = -stars%mx2,stars%mx2 
            DO k3 = -stars%mx3,stars%mx3 
               in = stars%ig(k1,k2,k3) 
               IF (in == 0) CYCLE 
               vpw_w(in) = Vpw_w_temp(k1,k2,k3) 
            ENDDO 
         ENDDO 
      ENDDO 
                                                                        
                                                                        
!*******************************************                            
! Use fft: faster                                                       
!*******************************************                            
           !use fft                                                     
      ELSE 
                                                                        
                                                                        
      ALLOCATE(V(1:s_dim(1),1:s_dim(2),                                 &
     &                      1:s_dim(3)))                                
      v = 0 
      DO k3=1,s_dim(3) 
       DO k2=1,s_dim(2) 
        DO k1=1,s_dim(1) 
              g1=k1-1 
              g2=k2-1 
              g3=k3-1 
              IF(g1>(s_mx1+0.01)) g1=g1-s_dim(1) 
              IF(g2>(s_mx2+0.01)) g2=g2-s_dim(2) 
              IF(g3>(s_mx3+0.01)) g3=g3-s_dim(3) 
              IF(abs(g1)>stars%mx1)CYCLE 
              IF(abs(g2)>stars%mx2)CYCLE 
              IF(abs(g3)>stars%mx3)CYCLE 
         in = stars%ig(g1,g2,g3) 
         IF (in == 0) CYCLE 
         V(k1,k2,k3) = vpw(in) 
        ENDDO 
       ENDDO 
      ENDDO 
                                                                        
                                                                        
         ALLOCATE(fftstepkomplett(1:s_dim(1),1:s_dim(2),                &
     &                                          1:s_dim(3)))            
         fftstepkomplett=cmplx(0.0,0.0) 
         DO k3=1,s_dim(3) 
          DO k2=1,s_dim(2) 
           DO k1=1,s_dim(1) 
              g1=k1-1 
              g2=k2-1 
              g3=k3-1 
              IF(g1>(s_mx1+0.01)) g1=g1-s_dim(1) 
              IF(g2>(s_mx2+0.01)) g2=g2-s_dim(2) 
              IF(g3>(s_mx3+0.01)) g3=g3-s_dim(3) 
              fftstepkomplett(k1,k2,k3)=stepkomplett(g1,g2,g3) 
           ENDDO 
          ENDDO 
         ENDDO 
                                                                        
                                                                        
      v=fft(v) 
      fftstepkomplett=fft(fftstepkomplett) 
                                                                        
      v=v*fftstepkomplett 
      V = fft(V,inv=.TRUE.)/size(V) 
                                                                        
                                     !potential shift of right side     
      IF(abs(pot_aux)>1.E-10)THEN 
         ALLOCATE(fftstephelp(1:s_dim(1),1:s_dim(2),                    &
     &                                          1:s_dim(3)))            
         fftstephelp=cmplx(0.0,0.0) 
         DO k3=1,s_dim(3) 
          DO k2=1,s_dim(2) 
           DO k1=1,s_dim(1) 
              g1=k1-1 
              g2=k2-1 
              g3=k3-1 
              IF(g1>(s_mx1+0.01)) g1=g1-s_dim(1) 
              IF(g2>(s_mx2+0.01)) g2=g2-s_dim(2) 
              IF(g3>(s_mx3+0.01)) g3=g3-s_dim(3) 
              fftstephelp(k1,k2,k3)=stephelpreg_r(g1,g2,g3) 
           ENDDO 
          ENDDO 
         ENDDO 
         fftstephelp=fftstephelp*cmplx(pot_aux,0.0) 
         v=v+fftstephelp 
      ENDIF 
                                                                        
                                                                        
      vpw_w=cmplx(0.0,0.0) 
         DO k3=1,s_dim(3) 
          DO k2=1,s_dim(2) 
           DO k1=1,s_dim(1) 
              g1=k1-1 
              g2=k2-1 
              g3=k3-1 
              IF(g1>(s_mx1+0.01)) g1=g1-s_dim(1) 
              IF(g2>(s_mx2+0.01)) g2=g2-s_dim(2) 
              IF(g3>(s_mx3+0.01)) g3=g3-s_dim(3) 
              IF(abs(g1)>stars%mx1)CYCLE 
              IF(abs(g2)>stars%mx2)CYCLE 
              IF(abs(g3)>stars%mx3)CYCLE 
              in = stars%ig(g1,g2,g3) 
              IF (in == 0) CYCLE 
!              vpw_w(in)=v(k1,k2,k3)                                    
              vpw_w(in)=v(k1,k2,k3)+vpw_w(in) 
           ENDDO 
          ENDDO 
         ENDDO 
         DO in=1,stars%nq3 
            vpw_w(in)=vpw_w(in)/stars%nstr(in) 
         ENDDO 
      ENDIF 
                                                                        
                                                                        
                                                                        
                                                                        
         CALL cpu_time(time2) 
         time2=time2-time1 
                                                                        

      DEALLOCATE(v,stepkomplett) 
                                                                        
                                                                        
                                                                        
                                                                        
                                                                        
      END SUBROUTINE gf_gspaceconvolve 
      !>                                                                
      !<-- S: gf_steps_pot_convolve(layer,stars,pot_aux,vpw,vpw
      SUBROUTINE gf_steps_pot_convolve(layer,stars,vpw,vpw_w)
!********************************************************************   
!     Calculate the convolution in g-space without making use of fft.   
!     Frank Freimuth, April 2007                                        
!********************************************************************   
      USE m_gf_types 
      USE m_gf_fft_singleton 
      IMPLICIT NONE 
      !<--Arguments                                                     
      INTEGER,INTENT(IN)       :: layer 

      TYPE(t_stars),INTENT(IN) :: stars 
                                            !Planewave potential/charge 
      COMPLEX,INTENT(IN)       :: vpw(:) 
                                             !Warped planewave/charge po
      COMPLEX,INTENT(OUT)      :: vpw_w(:,:) 
      !>                                                                
      !<-- Locals                                                       
      COMPLEX,ALLOCATABLE :: v(:,:,:) 
      COMPLEX,ALLOCATABLE :: stepphys(:,:,:) 
      COMPLEX,ALLOCATABLE :: fftstep(:,:,:) 
      INTEGER             :: i,k1,k2,k3,in,g1,g2,g3,h1,h2,h3,g3i
                                                                        
      !>                                                                
                                                                        
!Obtain the stepfunctions                                               
      ALLOCATE(stepphys(-s_mx1:s_mx1,-s_mx2:s_mx2,-s_mx3:s_mx3)) 
      CALL priv_loadstep(layer,stepphys = stepphys)
                                                                        
      ALLOCATE(V(1:s_DIM(1),1:s_DIM(2),                                 &
     &                      1:s_dim(3)))                                
      v = 0 
      !<-- put charge/potential on regular grid                         
      DO k3 = 1,s_DIM(3) 
         DO k2 = 1,s_DIM(2) 
            DO k1 = 1,s_DIM(1) 
               g1 = k1-1 
               g2 = k2-1 
               g3 = k3-1 
               IF(g1>(s_mx1+0.01)) g1 = g1-s_DIM(1) 
               IF(g2>(s_mx2+0.01)) g2 = g2-s_DIM(2) 
               IF(g3>(s_mx3+0.01)) g3 = g3-s_DIM(3) 
               IF(ABS(g1)>stars%mx1)CYCLE 
               IF(ABS(g2)>stars%mx2)CYCLE 
               IF(ABS(g3)>stars%mx3)CYCLE 
               in = stars%ig(g1,g2,g3) 
               IF (in   == 0) CYCLE 
               V(k1,k2,k3)    = vpw(in) 
            ENDDO 
         ENDDO 
      ENDDO 
      !>                                                                
      !<-- Put step-function on grid                                    
      ALLOCATE(fftstep(1:s_DIM(1),1:s_DIM(2),1:s_DIM(3))) 
      fftstep = CMPLX(0.0,0.0) 
      DO k3 = 1,s_DIM(3) 
         DO k2 = 1,s_DIM(2) 
            DO k1 = 1,s_DIM(1) 
               g1 = k1-1 
               g2 = k2-1 
               g3 = k3-1 
               IF(g1>(s_mx1+0.01)) g1 = g1-s_DIM(1) 
               IF(g2>(s_mx2+0.01)) g2 = g2-s_DIM(2) 
               IF(g3>(s_mx3+0.01)) g3 = g3-s_DIM(3) 
               print *, "Smooth step"
               IF(ABS(g1)>stars%mx1)CYCLE
               IF(ABS(g2)>stars%mx2)CYCLE
               IF(ABS(g3)>stars%mx3)CYCLE
               in = stars%ig(g1,g2,g3)
               IF (in   == 0) CYCLE
               fftstep(k1,k2,k3) = stepphys(g1,g2,g3) 
            ENDDO 
         ENDDO 
      ENDDO 
      !>                                                                
      v = fft(v) 
      fftstep = fft(fftstep) 
                                                                        
      v = v*fftstep 
      DEALLOCATE(fftstep) 
      V = fft(V,inv=.TRUE.)/size(V) 
                                                                        
      !<-- Get charge/potential from grid                               
      vpw_w = CMPLX(0.0,0.0) 
      DO k3 = 1,s_DIM(3) 
         IF (k3-1 <= stars%mx3) THEN
            g3 = k3-1
            g3i= k3
         ELSE 
            g3 = k3-s_DIM(3)-1
            g3i = 2*stars%mx3+g3+2
         ENDIF
         IF (ABS(g3)>stars%mx3) CYCLE

         DO k2 = 1,s_DIM(2) 
            g2 = k2-1 
            IF(g2>(s_mx2+0.01)) g2 = g2-s_DIM(2) 
            IF(ABS(g2)>stars%mx2)CYCLE 
            DO k1 = 1,s_DIM(1) 
               g1 = k1-1 
               IF(g1>(s_mx1+0.01)) g1 = g1-s_DIM(1) 
               IF(ABS(g1)>stars%mx1)CYCLE 
               in = stars%ig(g1,g2,g3)
               IF (in<1) CYCLE
               in = stars%ig(g1,g2,0) 
               vpw_w(g3i,stars%ig2(in)) = v(k1,k2,k3)+vpw_w(g3i           &
     &              ,stars%ig2(in))                                     
            ENDDO 
         ENDDO 
      ENDDO 
      DO in = 1,stars%nq2 
         vpw_w(:,in) = vpw_w(:,in)/stars%nstr2(in) 
      ENDDO 
      !>                                                                
      DEALLOCATE(v,stepphys) 
                                                                        
                                                                        
                                                                        
                                                                        
                                                                        
      END SUBROUTINE 
      !>                                                                
      !<-- S: priv_loadstep(layer,stepkomplett,stephelpreg_r,st
                                                                        
      SUBROUTINE priv_loadstep(layer,stepkomplett,stephelpreg_r&
     &     ,stephelpreg_l,stepphys,stepMTOUT,success)                   
!-----------------------------------------------                        
!   Read stepfunctions from file                                        
!   -only those arguments present will be read                          
!   -success can be used to test if stepfunctions                       
!       with correct dims are present                                   
!        D. Wortmann                                                    
!                                                                       
!   Modification for the g-space-stepfunctions.                         
!   Frank Freimuth, April 2007                                          
!-----------------------------------------------                        
      USE hdf5 
      USE m_hdf_tools 
      use m_juDFT 
      IMPLICIT NONE 
                                                                        
      !<-- Arguments                                                    
      INTEGER,INTENT(IN) :: layer 

      COMPLEX,OPTIONAL,INTENT(OUT) :: stepkomplett                      &
     &                                (-s_mx1:,-s_mx2:,-s_mx3:)         
      COMPLEX,OPTIONAL,INTENT(OUT) :: stephelpreg_r                     &
     &                                (-s_mx1:,-s_mx2:,-s_mx3:)         
      COMPLEX,OPTIONAL,INTENT(OUT) :: stephelpreg_l                     &
     &                                (-s_mx1:,-s_mx2:,-s_mx3:)         
      COMPLEX,OPTIONAL,INTENT(OUT) :: stepphys                          &
     &                                (-s_mx1:,-s_mx2:,-s_mx3:)         
      COMPLEX,OPTIONAL,INTENT(OUT) :: stepMTOUT(-s_mx1:,-s_mx2:,-s_mx3:) 
      LOGICAL,OPTIONAL,INTENT(OUT) :: success 
      !>                                                                
                                                                        
      !<-- Locals                                                       
      INTEGER(HID_T)        :: fid,gid,varid,ffid,access_prp 
      INTEGER               :: hdferr,hdfdim(3),d(3)
      REAL,ALLOCATABLE:: tmp_r(:,:,:) 
      REAL,ALLOCATABLE:: tmp_i(:,:,:) 
      !>                                                                

      access_prp=H5P_DEFAULT_f 

      !<-- Check if file gf_setup.hdf contains group stepsanaly         
                                                                        
      CALL io_hdfopen ('gf_setup.hdf',H5F_ACC_RDONLY_F, ffid, hdferr     &
     &     ,access_prp)                                                 
      CALL io_gopen(ffid,io_layername(layer),fid,hdferr) 
      IF(.NOT.io_groupexists(fid,'stepsanaly')) THEN 
         IF (.NOT.PRESENT(success)) CALL juDFT_error                      &
     &        ("No Stepfunction found1")                                
         success = .FALSE. 
         CALL io_gclose(fid,hdferr) 
         CALL io_hdfclose (ffid, hdferr) 
         RETURN 
      ENDIF 
      !>                                                                
!Group stepsanaly is there. Let's see if its dimensions are correct.    
      CALL io_gopen (fid,'stepsanaly',gid,hdferr) 
      CALL io_READ_att(gid, "dimensions" ,hdfdim) 
      IF (  ANY(hdfdim < s_dim)  ) THEN 
         IF (.NOT.PRESENT(success)) THEN 
            WRITE(*,*) "Layer:",layer 
            WRITE(*,*) "hdfdim:",hdfdim 
            WRITE(*,*) "s_dim:",s_dim 
                                                                        
            CALL juDFT_error                                              &
     &        ("Not enough data in Stepfunction found")                 
         ENDIF 
         success = .FALSE. 
         CALL io_gclose (gid, hdferr) 
         CALL io_gclose (fid, hdferr) 
         CALL io_hdfclose (ffid, hdferr) 
         RETURN 
      ENDIF 
                                                                        
!At this point it can be assumed that the correct stepfunctions are in t
      IF (PRESENT(success)) success = .TRUE. 
      d=(hdfdim-1)/2
      ALLOCATE(tmp_r(-d(1):d(1),-d(2):d(2),-d(3):d(3)))
      ALLOCATE(tmp_i(-d(1):d(1),-d(2):d(2),-d(3):d(3)))
                                                                        
      !<--read the stepfunctions needed                                 
      IF (PRESENT(stepkomplett)) THEN 
         CALL io_dopen(gid, 'stepkomplett', varid, hdferr) 
         CALL io_READ(varid, (/1,1,1,1/),                               &
     &         (/hdfdim(1),hdfdim(2),hdfdim(3),1/),tmp_r)
         CALL io_READ(varid, (/1,1,1,2/),                               &
     &         (/hdfdim(1),hdfdim(2),hdfdim(3),1/),tmp_i)
         stepkomplett = cmplx(tmp_r(-s_mx1:s_mx1,-s_mx2:s_mx2,-s_mx3:s_mx3),tmp_i(-s_mx1:s_mx1,-s_mx2:s_mx2,-s_mx3:s_mx3))
         CALL io_dclose(varid, hdferr) 
      ENDIF 
      IF (PRESENT(stephelpreg_r)) THEN 
         CALL io_dopen(gid, 'stephelpreg_r', varid, hdferr) 
         CALL io_READ(varid, (/1,1,1,1/),                               &
     &         (/hdfdim(1),hdfdim(2),hdfdim(3),1/),tmp_r)
         CALL io_READ(varid, (/1,1,1,2/),                               &
     &         (/hdfdim(1),hdfdim(2),hdfdim(3),1/),tmp_i)
         stephelpreg_r = cmplx(tmp_r(-s_mx1:s_mx1,-s_mx2:s_mx2,-s_mx3:s_mx3),tmp_i(-s_mx1:s_mx1,-s_mx2:s_mx2,-s_mx3:s_mx3))
         CALL io_dclose(varid, hdferr) 
      ENDIF 
      IF (PRESENT(stephelpreg_l)) THEN 
         CALL io_dopen(gid, 'stephelpreg_l', varid, hdferr) 
                  CALL io_READ(varid, (/1,1,1,1/),                               &
     &         (/hdfdim(1),hdfdim(2),hdfdim(3),1/),tmp_r)
         CALL io_READ(varid, (/1,1,1,2/),                               &
     &         (/hdfdim(1),hdfdim(2),hdfdim(3),1/),tmp_i)
         stephelpreg_l = cmplx(tmp_r(-s_mx1:s_mx1,-s_mx2:s_mx2,-s_mx3:s_mx3),tmp_i(-s_mx1:s_mx1,-s_mx2:s_mx2,-s_mx3:s_mx3))
         CALL io_dclose(varid, hdferr) 
      ENDIF 
      IF (PRESENT(stepphys)) THEN 
         CALL io_dopen(gid, 'stepphys', varid, hdferr)
                  CALL io_READ(varid, (/1,1,1,1/),                               &
     &         (/hdfdim(1),hdfdim(2),hdfdim(3),1/),tmp_r)
         CALL io_READ(varid, (/1,1,1,2/),                               &
     &         (/hdfdim(1),hdfdim(2),hdfdim(3),1/),tmp_i)
         stepphys = cmplx(tmp_r(-s_mx1:s_mx1,-s_mx2:s_mx2,-s_mx3:s_mx3),tmp_i(-s_mx1:s_mx1,-s_mx2:s_mx2,-s_mx3:s_mx3))
         CALL io_dclose(varid, hdferr) 
      ENDIF 
      IF (PRESENT(stepMTOUT)) THEN 
         CALL io_dopen(gid, 'stepMTOUT', varid, hdferr) 
         CALL io_READ(varid, (/1,1,1,1/),                               &
     &         (/hdfdim(1),hdfdim(2),hdfdim(3),1/),tmp_r)
         CALL io_READ(varid, (/1,1,1,2/),                               &
     &         (/hdfdim(1),hdfdim(2),hdfdim(3),1/),tmp_i)
         stepmtout = cmplx(tmp_r(-s_mx1:s_mx1,-s_mx2:s_mx2,-s_mx3:s_mx3),tmp_i(-s_mx1:s_mx1,-s_mx2:s_mx2,-s_mx3:s_mx3))
         CALL io_dclose(varid, hdferr) 
      ENDIF 
      !>                                                                
                                                                        
      CALL io_gclose (gid, hdferr) 
      CALL io_gclose (fid, hdferr) 
      CALL io_hdfclose (ffid,hdferr) 
      DEALLOCATE(tmp_r,tmp_i) 
      END SUBROUTINE priv_loadstep 
                                                                        
      !>                                                                
                                                                        
                                                                        
                                                                        
      !<-- S: priv_savesteps(layer,stepkomplett,stephelpreg_r,stephelpre
      SUBROUTINE priv_savesteps(layer,stepkomplett,stephelpreg_r,       &
     &     stephelpreg_l,stepphys,stepMTOUT)                            
!-----------------------------------------------                        
!           (last modified: 2004-00-00) D. Wortmann                     
!                                                                       
!   Modification for the g-space-stepfunctions.                         
!   Frank Freimuth, April 2007                                          
!-----------------------------------------------                        
      USE hdf5 
      USE m_hdf_tools 
      use m_juDFT 
      IMPLICIT NONE 
      !<-- Arguments                                                    
      INTEGER,INTENT(IN)::layer 
      COMPLEX,OPTIONAL,INTENT(IN) :: stepkomplett                       &
     &                               (-s_mx1:,-s_mx2:,-s_mx3:)          
      COMPLEX,OPTIONAL,INTENT(IN) :: stephelpreg_r                      &
     &                               (-s_mx1:,-s_mx2:,-s_mx3:)          
      COMPLEX,OPTIONAL,INTENT(IN) :: stephelpreg_l                      &
     &                               (-s_mx1:,-s_mx2:,-s_mx3:)          
      COMPLEX,OPTIONAL,INTENT(IN) :: stepphys                           &
     &                               (-s_mx1:,-s_mx2:,-s_mx3:)          
      COMPLEX,OPTIONAL,INTENT(IN) :: stepMTOUT(-s_mx1:,-s_mx2:,-s_mx3:) 
                                                                        
      !>                                                                
      !<-- Locals                                                       
      INTEGER(HID_T) :: fid,gid,varid,ffid 
      INTEGER        :: hdferr 
                                                                        
      !>                                                                
      CALL io_hdfopen ('gf_setup.hdf',H5F_ACC_RDWR_F,ffid, hdferr) 
      CALL io_gopen(ffid,io_layername(layer),fid,hdferr) 
      IF(io_groupexists(fid,'stepsanaly')) THEN 
         CALL io_gdelete(fid,'stepsanaly',hdferr) 
      ENDIF 
      CALL io_gcreate(fid,'stepsanaly',gid,hdferr) 
      CALL io_WRITE_att(gid, "dimensions" ,(/s_dim(1),s_dim(2),         &
     &                                          s_dim(3)/))             
                                                                        
                                                                        
      IF (PRESENT(stepkomplett)) THEN 
       CALL io_createvar(gid,"stepkomplett",H5T_NATIVE_DOUBLE,          &
     &    (/s_dim(1),s_dim(2),s_dim(3),2/),varid)                       
       CALL io_WRITE(varid,(/1,1,1,1/),(/s_dim(1),s_dim(2),s_dim(3),1/),&
     &                        real(stepkomplett(-s_mx1:s_mx1,           &
     &            -s_mx2:s_mx2,-s_mx3:s_mx3))                       )   
       CALL io_WRITE(varid,(/1,1,1,2/),(/s_dim(1),s_dim(2),s_dim(3),1/),&
     &                        imag(stepkomplett(-s_mx1:s_mx1,           &
     &            -s_mx2:s_mx2,-s_mx3:s_mx3))                       )   
       CALL io_dclose(varid,hdferr) 
      ENDIF 
                                                                        
      IF (PRESENT(stephelpreg_r)) THEN 
       CALL io_createvar(gid,"stephelpreg_r",H5T_NATIVE_DOUBLE,         &
     &    (/s_dim(1),s_dim(2),s_dim(3),2/),varid)                       
       CALL io_WRITE(varid,(/1,1,1,1/),(/s_dim(1),s_dim(2),s_dim(3),1/),&
     &                        REAL(stephelpreg_r(-s_mx1:s_mx1,          &
     &            -s_mx2:s_mx2,-s_mx3:s_mx3))                       )   
       CALL io_WRITE(varid,(/1,1,1,2/),(/s_dim(1),s_dim(2),s_dim(3),1/),&
     &                        imag(stephelpreg_r(-s_mx1:s_mx1,          &
     &            -s_mx2:s_mx2,-s_mx3:s_mx3))                       )   
       CALL io_dclose(varid,hdferr) 
      ENDIF 
                                                                        
      IF (PRESENT(stephelpreg_l)) THEN 
       CALL io_createvar(gid,"stephelpreg_l",H5T_NATIVE_DOUBLE,         &
     &    (/s_dim(1),s_dim(2),s_dim(3),2/),varid)                       
       CALL io_WRITE(varid,(/1,1,1,1/),(/s_dim(1),s_dim(2),s_dim(3),1/),&
     &                        REAL(stephelpreg_l))                      
       CALL io_WRITE(varid,(/1,1,1,2/),(/s_dim(1),s_dim(2),s_dim(3),1/),&
     &                        imag(stephelpreg_l))                      
       CALL io_dclose(varid,hdferr) 
      ENDIF 
                                                                        
      IF (PRESENT(stepMTOUT)) THEN 
       CALL io_createvar(gid,"stepMTOUT",H5T_NATIVE_DOUBLE,             &
     &    (/s_dim(1),s_dim(2),s_dim(3),2/),varid)                       
       CALL io_WRITE(varid,(/1,1,1,1/),(/s_dim(1),s_dim(2),s_dim(3),1/),&
     &                        REAL(stepMTOUT))                          
       CALL io_WRITE(varid,(/1,1,1,2/),(/s_dim(1),s_dim(2),s_dim(3),1/),&
     &                        imag(stepMTOUT))                          
       CALL io_dclose(varid,hdferr) 
      ENDIF 
                                                                        
      IF (PRESENT(stepphys)) THEN 
       CALL io_createvar(gid,"stepphys",H5T_NATIVE_DOUBLE,              &
     &    (/s_dim(1),s_dim(2),s_dim(3),2/),varid)                       
       CALL io_WRITE(varid,(/1,1,1,1/),(/s_dim(1),s_dim(2),s_dim(3),1/),&
     &                        REAL(stepphys))                           
       CALL io_WRITE(varid,(/1,1,1,2/),(/s_dim(1),s_dim(2),s_dim(3),1/),&
     &                        imag(stepphys))                           
       CALL io_dclose(varid,hdferr) 
      ENDIF 
                                                                        
                                                                        
      !close group&file                                                 
      CALL io_gclose (gid, hdferr) 
      CALL io_gclose (fid, hdferr) 
      CALL io_hdfclose(ffid,hdferr) 
      END SUBROUTINE priv_savesteps 
      !>                                                                
                                                                        
                                                                        
      SUBROUTINE priv_generate_MTstep(omtil,bmat,taual,rmt,             &
     &            stepmtfui)                                            
!-----------------------------------------------                        
!  No equivalent atoms used here, size(taual,2) and size(rmt) must be   
!  the same!                                                            
!           (last modified: 06-10-09) D. Wortmann                       
!------------------------------------------------                       
!  Modified version: Calculate the Fourier transform                    
!  of the stepfunction which equals one in the MT-spheres               
!  and zero outside of the MTs.                                         
!  Frank Freimuth, April 2007                                           
!-----------------------------------------------                        
      USE m_constants,ONLY:PIMACH 
      IMPLICIT NONE 
      !<-- Arguments                                                    
      REAL ,INTENT(IN)       :: omtil,bmat(3,3) 
      REAL ,INTENT(IN)       :: taual(:,:) 
      REAL ,INTENT(IN)       :: rmt(:) 
                                                                        
      COMPLEX,INTENT(OUT)    :: stepmtfui(-s_mx1:,-s_mx2:,-s_mx3:) 
      !>                                                                
      !<-- Locals                                                       
      INTEGER             :: n,k1,k2,k3 
      REAL                :: c,rg(3),gs,th,gvec(3),tpi 
      !>                                                                
                                                                        
                                                                        
                                                                        
      tpi=2.*pimach() 
                                                                        
      stepmtfui = cmplx(0.0,0.0) 
                                                                        
      DO  n = 1,SIZE(rmt) 
         stepmtfui(0,0,0) = stepmtfui(0,0,0)+                           &
     &               CMPLX(2.0*tpi/3.*rmt(n)**3./omtil,0.0)             
      ENDDO 
                                                                        
      !  Loop over all atoms                                            
      ! No equivalent atoms used here!                                  
      DO n = 1,SIZE(rmt) 
         c = 3.*2.0*tpi/3.*rmt(n)**3/omtil 
         DO k1 = -s_mx1,s_mx1 
            gvec(1)=k1 
            DO k2 = -s_mx2,s_mx2 
               gvec(2)=k2 
               DO k3 = -s_mx3,s_mx3 
                  IF (k1 == 0.AND.k2 == 0.AND.k3 == 0) CYCLE 
                  gvec(3)=k3 
                  rg  = MATMUL(gvec,bmat) 
                  gs = SQRT(dot_PRODUCT(rg,rg))*rmt(n) 
                  th =-tpi*dot_PRODUCT(gvec,taual(:,n)) 
                  stepmtfui(k1,k2,k3) =                                 &
     &                 stepmtfui(k1,k2,k3)+(c*(SIN(gs)/gs-COS(gs))/     &
     &                 (gs*gs))*CMPLX(COS(th),SIN(th))                  
               ENDDO 
            ENDDO 
         ENDDO 
      ENDDO 
                                                                        
      END SUBROUTINE priv_generate_mtstep 
                                                                        
                                                                        
      SUBROUTINE priv_generate_accastep(bmat33,amat33,c,stepacca) 
!----------------------------------------------------------------       
!  Calculate the Fourier transform of the stepfunction which is         
!  zero in the volume [-c/2,c/2] and unity everywhere else.             
!             Frank Freimuth, April 2007                                
!----------------------------------------------------------------       
      IMPLICIT NONE 
      !<-- Arguments                                                    
      REAL   ,INTENT(IN)     :: amat33,bmat33 
      REAL ,INTENT(IN)       :: c 
                                                                        
      COMPLEX,INTENT(OUT)    :: stepacca(-s_mx1:,-s_mx2:,-s_mx3:) 
                                                                        
      !>                                                                
      !<--Locals                                                        
      INTEGER             :: n 
      REAL                :: gvec,z1,z2 
      !>                                                                
                                                                        
      z1=c/2.0 
      z2=amat33-z1 
      stepacca=cmplx(0.0,0.0) 
      stepacca(0,0,0) = CMPLX((amat33-c) /amat33,0.0) 
      !g_parallel == 0                                                  
!      DO n = -s_mx3,s_mx3                                              
!         if(n==0) cycle                                                
!         gvec = bmat33*n                                               
!         stepacca(0,0,n) = CMPLX(0.0,1/amat33/gvec)* (                 
!     $      EXP(CMPLX(0.,-gvec*z2))-                                   
!     $      EXP(CMPLX(0.,-gvec*z1))                    )               
!      ENDDO                                                            
                                                                        
      !g_parallel == 0                                                  
      DO n = -s_mx3,s_mx3 
         IF(n==0) CYCLE 
         gvec = bmat33*n 
         stepacca(0,0,n)=cmplx(-2.0/amat33/gvec*sin(gvec*z1),0.0) 
      ENDDO 
                                                                        
                                                                        
      END SUBROUTINE priv_generate_accastep 
                                                                        
                                                                        
      SUBROUTINE priv_generate_cdstepright(bmat33,amat33,c,d,           &
     &                                     cdstepright)                 
!----------------------------------------------------------------       
!  Calculate the Fourier transform of the stepfunction which is         
!  one in the volume [c/2,d/2] and zero everywhere else.                
!             Frank Freimuth, April 2007                                
!----------------------------------------------------------------       
      IMPLICIT NONE 
      !<-- Arguments                                                    
      REAL   ,INTENT(IN)     :: amat33,bmat33 
      REAL ,INTENT(IN)       :: c 
      REAL, INTENT(IN)       :: d 
                                                                        
      COMPLEX,INTENT(OUT)    :: cdstepright(-s_mx1:,-s_mx2:,-s_mx3:) 
                                                                        
      !>                                                                
      !<--Locals                                                        
      INTEGER             :: n 
      REAL                :: gvec,z1,z2 
      !>                                                                
                                                                        
      z1=c/2.0 
      z2=d/2.0 
      cdstepright=cmplx(0.0,0.0) 
      cdstepright(0,0,0) = CMPLX((d-c)/amat33/2.0,0.0) 
      !g_parallel == 0                                                  
      DO n = -s_mx3,s_mx3 
         IF(n==0) CYCLE 
         gvec = bmat33*n 
         cdstepright(0,0,n) = CMPLX(0.0,1/amat33/gvec)* (               &
     &      EXP(CMPLX(0.,-gvec*z2))-                                    &
     &      EXP(CMPLX(0.,-gvec*z1))                    )                
      ENDDO 
                                                                        
      END SUBROUTINE priv_generate_cdstepright 
                                                                        
                                                                        
      SUBROUTINE priv_generate_ballpartstep(omtil,bmat,c,amat33,taual,  &
     &            rmt,step)                                             
!-----------------------------------------------                        
!  Calculate the Fourier transform                                      
!  of the stepfunction which equals one in a part of                    
!  the MT-spheres and which is zero otherwise.                          
!  Frank Freimuth, April 2007                                           
!-----------------------------------------------                        
      USE m_constants,ONLY:PIMACH 
      USE m_gf_ballpartfourier,ONLY:ballpartfou
      IMPLICIT NONE 
      !<-- Arguments                                                    
      REAL ,INTENT(IN)     :: omtil,bmat(3,3) 
      REAL, INTENT(IN)     :: c,amat33 
      REAL ,INTENT(IN)     :: taual(:,:) 
      REAL ,INTENT(IN)     :: rmt(:) 
                                                                        
      COMPLEX,INTENT(OUT)  :: step(-s_mx1:,-s_mx2:,-s_mx3:) 
      !>                                                                
      !<-- Locals                                                       
      INTEGER             :: n,k1,k2,k3 
      REAL                :: rg(3),gs,th,gvec(3),tpi,pi 
      REAL zcoor,kpara,korthogo 
      COMPLEX value 
      !>                                                                
                                                                        
                                                                        
      pi=pimach() 
      tpi=2.*pimach() 
                                                                        
      step = cmplx(0.0,0.0) 
                                                                        
      DO n = 1,SIZE(rmt) 
                                                      !cut right dividin
       IF(abs(c/2.0-taual(3,n)*amat33)<rmt(n))THEN 
         zcoor=c/2.0-taual(3,n)*amat33 
         IF(zcoor>=0)THEN 
            step(0,0,0)=step(0,0,0)+pi*(  2.0/3.0*(rmt(n))**3-          &
     &               (rmt(n))**2*zcoor  +   zcoor**3/3.0     )/omtil    
         ELSE 
            step(0,0,0)=step(0,0,0)+pi*(  2.0/3.0*(rmt(n))**3+          &
     &               (rmt(n))**2*zcoor  -   zcoor**3/3.0     )/omtil    
         ENDIF 
         DO k1 = -s_mx1,s_mx1 
            gvec(1)=k1 
            DO k2 = -s_mx2,s_mx2 
               gvec(2)=k2 
               DO k3 = -s_mx3,s_mx3 
                  IF (k1 == 0.AND.k2 == 0.AND.k3 == 0) CYCLE 
                  gvec(3)=k3 
                  rg  = MATMUL(gvec,bmat) 
                  gs = SQRT(dot_PRODUCT(rg,rg))*rmt(n) 
                  th =-tpi*dot_PRODUCT(gvec,taual(:,n)) 
                  kpara=rg(3) 
                  korthogo=sqrt( rg(1)*rg(1) + rg(2)*rg(2) ) 
                  IF(zcoor>=0)THEN 
                    CALL ballpartfou(rmt(n),zcoor,kpara,korthogo,value) 
                    step(k1,k2,k3) = step(k1,k2,k3)+                    &
     &                       conjg(value)*cmplx(cos(th),sin(th))/omtil  
                  ELSE 
                    CALL ballpartfou(rmt(n),-zcoor,kpara,korthogo,value) 
                    step(k1,k2,k3) = step(k1,k2,k3)+                    &
     &                             value *cmplx(cos(th),sin(th))/omtil  
                  ENDIF 
               ENDDO 
            ENDDO 
         ENDDO 
              !right dividing plane is cut                              
       ENDIF 
                                                       !cut left dividin
       IF(abs(-c/2.0-taual(3,n)*amat33)<rmt(n))THEN 
         zcoor=-c/2.0-taual(3,n)*amat33 
         IF(zcoor>=0)THEN 
            step(0,0,0)=step(0,0,0)+pi*(  2.0/3.0*(rmt(n))**3-          &
     &               (rmt(n))**2*zcoor  +   zcoor**3/3.0     )/omtil    
         ELSE 
            step(0,0,0)=step(0,0,0)+pi*(  2.0/3.0*(rmt(n))**3+          &
     &               (rmt(n))**2*zcoor  -   zcoor**3/3.0     )/omtil    
         ENDIF 
         DO k1 = -s_mx1,s_mx1 
            gvec(1)=k1 
            DO k2 = -s_mx2,s_mx2 
               gvec(2)=k2 
               DO k3 = -s_mx3,s_mx3 
                  IF (k1 == 0.AND.k2 == 0.AND.k3 == 0) CYCLE 
                  gvec(3)=k3 
                  rg  = MATMUL(gvec,bmat) 
                  gs = SQRT(dot_PRODUCT(rg,rg))*rmt(n) 
                  th =-tpi*dot_PRODUCT(gvec,taual(:,n)) 
                  kpara=rg(3) 
                  korthogo=sqrt( rg(1)*rg(1) + rg(2)*rg(2) ) 
                  IF(zcoor>=0)THEN 
                    CALL ballpartfou(rmt(n),zcoor,kpara,korthogo,value) 
                    step(k1,k2,k3) = step(k1,k2,k3)+                    &
     &                       conjg(value)*cmplx(cos(th),sin(th))/omtil  
                  ELSE 
                    CALL ballpartfou(rmt(n),-zcoor,kpara,korthogo,value) 
                    step(k1,k2,k3) = step(k1,k2,k3)+                    &
     &                             value *cmplx(cos(th),sin(th))/omtil  
                  ENDIF 
               ENDDO 
            ENDDO 
         ENDDO 
              !left dividing plane is cut                               
       ENDIF 
            !loop over atoms                                            
      ENDDO 
                                                                        
      END SUBROUTINE priv_generate_ballpartstep 
                                                                        
                                                                        
      SUBROUTINE priv_generate_ballpartstep_right(                      &
     &            omtil,bmat,c,amat33,taual,                            &
     &            rmt,step)                                             
!-----------------------------------------------                        
!  Calculate the Fourier transform                                      
!  of the stepfunction which equals one in a part of                    
!  the MT-spheres and which is zero otherwise.                          
!  Frank Freimuth, April 2007                                           
!-----------------------------------------------                        
      USE m_constants,ONLY:PIMACH 
      USE m_gf_ballpartfourier,ONLY:ballpartfou
      IMPLICIT NONE 
      !<-- Arguments                                                    
      REAL ,INTENT(IN)     :: omtil,bmat(3,3) 
      REAL, INTENT(IN)     :: c,amat33 
      REAL ,INTENT(IN)     :: taual(:,:) 
      REAL ,INTENT(IN)     :: rmt(:) 
                                                                        
      COMPLEX,INTENT(OUT)  :: step(-s_mx1:,-s_mx2:,-s_mx3:) 
      !>                                                                
      !<-- Locals                                                       
      INTEGER             :: n,k1,k2,k3 
      REAL                :: rg(3),gs,th,gvec(3),tpi,pi 
      REAL zcoor,kpara,korthogo 
      COMPLEX value 
      !>                                                                
                                                                        
                                                                        
      pi=pimach() 
      tpi=2.*pimach() 
                                                                        
      step = cmplx(0.0,0.0) 
                                                                        
      DO n = 1,SIZE(rmt) 
                                                      !cut right dividin
       IF(abs(c/2.0-taual(3,n)*amat33)<rmt(n))THEN 
         zcoor=c/2.0-taual(3,n)*amat33 
         IF(zcoor>=0)THEN 
            step(0,0,0)=step(0,0,0)+pi*(  2.0/3.0*(rmt(n))**3-          &
     &               (rmt(n))**2*zcoor  +   zcoor**3/3.0     )/omtil    
         ELSE 
            step(0,0,0)=step(0,0,0)+pi*(  2.0/3.0*(rmt(n))**3+          &
     &               (rmt(n))**2*zcoor  -   zcoor**3/3.0     )/omtil    
         ENDIF 
         DO k1 = -s_mx1,s_mx1 
            gvec(1)=k1 
            DO k2 = -s_mx2,s_mx2 
               gvec(2)=k2 
               DO k3 = -s_mx3,s_mx3 
                  IF (k1 == 0.AND.k2 == 0.AND.k3 == 0) CYCLE 
                  gvec(3)=k3 
                  rg  = MATMUL(gvec,bmat) 
                  gs = SQRT(dot_PRODUCT(rg,rg))*rmt(n) 
                  th =-tpi*dot_PRODUCT(gvec,taual(:,n)) 
                  kpara=rg(3) 
                  korthogo=sqrt( rg(1)*rg(1) + rg(2)*rg(2) ) 
                  IF(zcoor>=0)THEN 
                    CALL ballpartfou(rmt(n),zcoor,kpara,korthogo,value) 
                    step(k1,k2,k3) = step(k1,k2,k3)+                    &
     &                       conjg(value)*cmplx(cos(th),sin(th))/omtil  
                  ELSE 
                    CALL ballpartfou(rmt(n),-zcoor,kpara,korthogo,value) 
                    step(k1,k2,k3) = step(k1,k2,k3)+                    &
     &                             value *cmplx(cos(th),sin(th))/omtil  
                  ENDIF 
               ENDDO 
            ENDDO 
         ENDDO 
              !right dividing plane is cut                              
       ENDIF 
            !loop over atoms                                            
      ENDDO 
                                                                        
      END SUBROUTINE priv_generate_ballpartstep_right 
                                                                        
      SUBROUTINE priv_generate_ballpartstep_left(                       &
     &            omtil,bmat,c,amat33,taual,                            &
     &            rmt,step)                                             
!-----------------------------------------------                        
!  Calculate the Fourier transform                                      
!  of the stepfunction which equals one in a part of                    
!  the MT-spheres and which is zero otherwise.                          
!  Frank Freimuth, April 2007                                           
!-----------------------------------------------                        
      USE m_constants,ONLY:PIMACH 
      USE m_gf_ballpartfourier,ONLY:ballpartfou
      IMPLICIT NONE 
      !<-- Arguments                                                    
      REAL ,INTENT(IN)     :: omtil,bmat(3,3) 
      REAL, INTENT(IN)     :: c,amat33 
      REAL ,INTENT(IN)     :: taual(:,:) 
      REAL ,INTENT(IN)     :: rmt(:) 
                                                                        
      COMPLEX,INTENT(OUT)  :: step(-s_mx1:,-s_mx2:,-s_mx3:) 
      !>                                                                
      !<-- Locals                                                       
      INTEGER             :: n,k1,k2,k3 
      REAL                :: rg(3),gs,th,gvec(3),tpi,pi 
      REAL zcoor,kpara,korthogo 
      COMPLEX value 
      !>                                                                
                                                                        
                                                                        
      pi=pimach() 
      tpi=2.*pimach() 
                                                                        
      step = cmplx(0.0,0.0) 
                                                                        
      DO n = 1,SIZE(rmt) 
                                                       !cut left dividin
       IF(abs(-c/2.0-taual(3,n)*amat33)<rmt(n))THEN 
         zcoor=-c/2.0-taual(3,n)*amat33 
         IF(zcoor>=0)THEN 
            step(0,0,0)=step(0,0,0)+pi*(  2.0/3.0*(rmt(n))**3-          &
     &               (rmt(n))**2*zcoor  +   zcoor**3/3.0     )/omtil    
         ELSE 
            step(0,0,0)=step(0,0,0)+pi*(  2.0/3.0*(rmt(n))**3+          &
     &               (rmt(n))**2*zcoor  -   zcoor**3/3.0     )/omtil    
         ENDIF 
         DO k1 = -s_mx1,s_mx1 
            gvec(1)=k1 
            DO k2 = -s_mx2,s_mx2 
               gvec(2)=k2 
               DO k3 = -s_mx3,s_mx3 
                  IF (k1 == 0.AND.k2 == 0.AND.k3 == 0) CYCLE 
                  gvec(3)=k3 
                  rg  = MATMUL(gvec,bmat) 
                  gs = SQRT(dot_PRODUCT(rg,rg))*rmt(n) 
                  th =-tpi*dot_PRODUCT(gvec,taual(:,n)) 
                  kpara=rg(3) 
                  korthogo=sqrt( rg(1)*rg(1) + rg(2)*rg(2) ) 
                  IF(zcoor>=0)THEN 
                    CALL ballpartfou(rmt(n),zcoor,kpara,korthogo,value) 
                    step(k1,k2,k3) = step(k1,k2,k3)+                    &
     &                       conjg(value)*cmplx(cos(th),sin(th))/omtil  
                  ELSE 
                    CALL ballpartfou(rmt(n),-zcoor,kpara,korthogo,value) 
                    step(k1,k2,k3) = step(k1,k2,k3)+                    &
     &                             value *cmplx(cos(th),sin(th))/omtil  
                  ENDIF 
               ENDDO 
            ENDDO 
         ENDDO 
              !left dividing plane is cut                               
       ENDIF 
            !loop over atoms                                            
      ENDDO 
                                                                        
      END SUBROUTINE priv_generate_ballpartstep_left 
                                                                        
      END                                           
