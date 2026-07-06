!--------------------------------------------------------------------------------
! Copyright (c) 2026 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions 
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------
      MODULE m_gf_rsmesh 
      use m_juDFT
!*****************************************************************      
! DESC:This module contains subroutines and datatypes to deal with      
! quantities on a real space mesh.                                      
! Public:                                                               
!      TYPE rsdata:contains the data on RS-grid                         
!      SUBROUTINE: gf_rstors interpolates data onto new RS-grid         
!      SUBROUTINE: gf_rsForStars generates a RS-datatype correctly      
!                  dimensionized for the stars                          
!      SUBROUTINE: gf_rstopw transforms data from rs to pw              
!      SUBROUTINE: gf_pwtors transforms data from pw to rs              
!      SUBROUTINE: gf_rswrite writes rs-data to hdf5-file               
!      SUBROUTINE: gf_rsread reads rs-data from hdf5-file               
!                          Daniel Wortmann, Fri Oct 18 10:53:49 2002    
!   Last Modified: $Author: wortmann $ $Date: 2008/04/14 15:44:13 $     
!*****************************************************************      
      USE m_gf_fft,ONLY:gf_fft3d,GF_FFT_TO_G_SPACE,GF_FFT_TO_REAL_SPACE 
      use m_juDFT 
      IMPLICIT NONE
      PRIVATE 
      TYPE t_rsdata 
        REAL         :: box(3) 
                                !Origin in Z-direction                  
        REAL         :: z_zero 
        INTEGER      :: grid(3) 
        REAL,POINTER :: DATA(:,:,:) 
      END TYPE 
                                                                        
      INTERFACE gf_rsread 
         MODULE PROCEDURE gf_rsread_filename,gf_rsread_fid 
      END INTERFACE 
                                                                        
      INTERFACE gf_rswrite 
         MODULE PROCEDURE gf_rswrite_filename,gf_rswrite_fid 
      END INTERFACE 
                                                                        
                                                                        
      PUBLIC::t_rsdata,gf_rstors,gf_rsForStars,gf_pwtors,gf_rstopw      &
     &     ,gf_rswrite,gf_rsread,gf_rssmoothz,gf_rscutbox,gf_rsputbox   
                                                                        
      CONTAINS 
      !<--S: gf_rstors(old,new,trans)                                   
      SUBROUTINE gf_rstors(old,new,trans) 
!*****************************************************************      
! DESC: Old rsdata is transformed to new rsdata. The vector trans       
! specifies the translation of the origin                               
!                                                                       
!                        Daniel Wortmann, Fri Oct 18 10:53:49 2002      
!*****************************************************************      
      IMPLICIT NONE 
!     Arguments                                                         
      TYPE(t_rsdata),INTENT(IN) ::old 
      TYPE(t_rsdata),INTENT(INOUT)::new 
      REAL,INTENT(IN)           ::trans(3) 
                                                                        
      INTEGER::xi,yi,zi 
      REAL   ::xr,yr,zr 
                                                                        
      DO xi=1,new%grid(1) 
         !Get real position                                             
         IF (xi<=new%grid(1)/2) THEN 
            xr=(xi-1.)/new%grid(1)*new%box(1)+trans(1) 
         ELSE 
            xr=(xi-new%grid(1)-1.)/new%grid(1)*new%box(1)+trans(1) 
         ENDIF 
         !transform to old coord. system                                
         DO WHILE (xr<-old%box(1)/2); xr=xr+old%box(1);ENDDO 
                                                            !move to -ol
         DO WHILE (xr>old%box(1)/2) ;xr=xr-old%box(1);ENDDO 
                                      !Transfrom to index               
         xr=xr/old%box(1)*old%grid(1) 
                                        !swap negative coord            
         IF (xr<0.0)  xr=xr+old%grid(1) 
         xr=xr+1 
                                                      !should not happen
         IF (floor(xr)>old%grid(1)) xr=xr-old%grid(1) 
         DO yi=1,new%grid(2) 
            !Get real position                                          
            IF (yi<=new%grid(2)/2) THEN 
               yr=(yi-1.)/new%grid(2)*new%box(2)+trans(2) 
            ELSE 
               yr=(yi-new%grid(2)-1.)/new%grid(2)*new%box(2)+trans(2) 
            ENDIF 
            !transform to old coord. system                             
            DO WHILE (yr<-old%box(2)/2); yr=yr+old%box(2);ENDDO 
                                                               !move to 
            DO WHILE (yr>old%box(2)/2) ;yr=yr-old%box(2);ENDDO 
                                                            !Transfrom t
            yr=yr/old%box(2)*old%grid(2) 
                                             !swap negative coord       
            IF (yr<0.0)  yr=yr+old%grid(2) 
            yr=yr+1 
            IF (floor(yr)>old%grid(2)) yr=yr-old%grid(2) 
            DO zi=1,new%grid(3) 
               !Get real position                                       
               IF (zi<=new%grid(3)/2) THEN 
                  zr=(zi-1.)/new%grid(3)*new%box(3)+trans(3) 
               ELSE 
                  zr=(zi-new%grid(3)-1.)/new%grid(3)*new%box(3)         &
     &                 +trans(3)                                        
               ENDIF 
               !transform to old coord. system                          
               DO WHILE (zr<-old%box(3)/2); zr=zr+old%box(3);ENDDO 
                                                                  !move 
               DO WHILE (zr>old%box(3)/2) ;zr=zr-old%box(3);ENDDO 
                                                               !Transfro
               zr=zr/old%box(3)*old%grid(3) 
                                                           !swap negativ
               IF (zr<0.0)  zr=zr+old%grid(3) 
               zr=zr+1 
               IF (floor(zr)>old%grid(3)) zr=zr-old%grid(3) 
               !write(*,'(3i5,3f12.5)')xi,yi,zi,xr,yr,zr                
               new%data(xi,yi,zi)=interpolate(old%data,xr,yr,zr) 
            ENDDO 
         ENDDO 
      ENDDO 
      END SUBROUTINE 
      !>                                                                
      !<--S:  gf_rsForStars(stars,rs)                                   

      SUBROUTINE gf_rsForStars(stars,rs) 
!*****************************************************************      
! DESC: Sets Box and grid of rs to value corresponding to stars         
!*****************************************************************      
      USE m_constants,ONLY:PIMACH 
      USE m_gf_types
      IMPLICIT NONE 
!     Arguments                                                         
      TYPE(t_stars),INTENT(IN)   ::stars 
      TYPE(t_rsdata),INTENT(OUT) ::rs 
                                                                        
      !calculate size of box from stars                                 
      rs%box(1)=stars%sk3(stars%ig(1,0,0)) 
      rs%box(2)=stars%sk3(stars%ig(0,1,0)) 
      rs%box(3)=stars%sk3(stars%ig(0,0,1)) 
      rs%box=2*PIMACH()/rs%box 
                                                                        
      !number of gid-points for FFT!                                    
                                                                        
      rs%grid(1)=stars%mx1*3 
      rs%grid(2)=stars%mx2*3 
      rs%grid(3)=stars%mx3*3 
                                                                        
      rs%z_zero=0 
      !allocate storage                                                 
                                                                        
      !if (ALLOCATED(rs%data)) DEALLOCATE( rs%data)                     
      ALLOCATE(rs%data(rs%grid(1),rs%grid(2),rs%grid(3))) 
      END SUBROUTINE 

      !>                                                                
      !<--S: gf_rstopw(rs,stars,pw)                                     

      SUBROUTINE gf_rstopw(rs,stars,pw) 
!*****************************************************************      
! DESC: Transforms data to rez. Space                                   
!*****************************************************************      
      USE m_gf_types
      IMPLICIT NONE 
!     Arguments                                                         
      TYPE(t_rsdata),INTENT(IN) ::rs 
      TYPE(t_stars),INTENT(IN)  ::stars 
      COMPLEX,INTENT(OUT)       ::pw(:) 
                                                                        
      REAL,ALLOCATABLE::afft(:) 
      !check dim                                                        
      IF (size(pw)/=stars%nq3.OR.                                     &
     &     stars%mx1*3/=rs%grid(1).OR.stars%mx2*3/=rs%grid(2)         &
     &     .OR.stars%mx3*3/=rs%grid(3))                               &
     &      CALL juDFT_error("Invalid Dimensions in gf_rstopw",calledby="gf_rsmesh.F90")
                                                                        
      !Allocate DATASPACE                                               
      ALLOCATE(afft(rs%grid(1)*rs%grid(2)*rs%grid(3))) 
      afft=reshape(rs%data,(/rs%grid(1)*rs%grid(2)*rs%grid(3)/)) 
                                                                        
      CALL gf_fft3d(afft,pw,stars,GF_FFT_TO_G_SPACE) 
                                                                        
      DEALLOCATE(afft) 
      END SUBROUTINE 

      !>                                                                
      !<--S: gf_pwtors(pw,stars,rs)                                     

      SUBROUTINE gf_pwtors(pw,stars,rs) 
!*****************************************************************      
! DESC: Transforms data to rez. Space                                   
!*****************************************************************      
      USE m_gf_types
      IMPLICIT NONE 
!     Arguments                                                         
      TYPE(t_rsdata),INTENT(INOUT) ::rs 
      TYPE(t_stars),INTENT(IN)   ::stars 
                                          !is INTENT(IN) only           
      COMPLEX,INTENT(INOUT)      :: pw(:) 
                                                                        
      REAL,ALLOCATABLE::afft(:) 
      !check dim                                                        
      IF (size(pw)/=stars%nq3.OR.                                     &
     &     stars%mx1*3/=rs%grid(1).OR.stars%mx2*3/=rs%grid(2)       &
     &     .OR.stars%mx3*3/=rs%grid(3))                               &
     &      CALL juDFT_error("Invalid Dimensions in gf_pwtors",calledby="gf_rsmesh.F90")
                                                                        
      !Allocate DATASPACE                                               
      ALLOCATE(afft(rs%grid(1)*rs%grid(2)*rs%grid(3))) 
      CALL gf_fft3d(afft,pw,stars,GF_FFT_TO_REAL_SPACE) 
      rs%data=reshape(afft,(/rs%grid(1),rs%grid(2),rs%grid(3)/)) 
      DEALLOCATE(afft) 
      END SUBROUTINE 

      !>                                                                
      !<--S: gf_rswrite_filename(rs,filename,dataname)                  
      SUBROUTINE gf_rswrite_filename(rs,filename,dataname) 
!*****************************************************************      
! DESC: Write the rs-data to the hdf5-file specified by filename        
!       The dataset and its dims will be named by dataname              
!                                                                       
!                        Daniel Wortmann, Fri Oct 18 10:53:49 2002      
!*****************************************************************      
      USE hdf5 
      USE m_hdf_tools 
      IMPLICIT NONE 
!     Arguments                                                         
      TYPE(t_rsdata),INTENT(IN) ::rs 
      CHARACTER*(*),INTENT(IN)  ::filename,dataname 
                                                                        
      INTEGER(HID_T)::fid 
      INTEGER       ::hdferr 
      LOGICAL:: f_exist 
                                                                        
      INQUIRE(FILE=filename//'.hdf',EXIST=f_exist) 
      IF (f_exist) THEN 
         CALL io_hdfopen (filename//'.hdf', H5F_ACC_RDWR_F, fid,hdferr) 
      ELSE 
         CALL h5fcreate_f (filename//'.hdf', H5F_ACC_TRUNC_F,fid,hdferr) 
      ENDIF 
      CALL gf_rswrite_fid(rs,fid,dataname) 
      CALL io_hdfclose(fid,hdferr) 
      END SUBROUTINE 
      !>                                                                
      !<--S: gf_rswrite_fid(rs,fid,dataname)                            
                                                                        
      SUBROUTINE gf_rswrite_fid(rs,fid,dataname) 
!*****************************************************************      
! DESC: Write the rs-data to the hdf5-file specified by filename        
!       The dataset and its dims will be named by dataname              
!                                                                       
!                        Daniel Wortmann, Fri Oct 18 10:53:49 2002      
!*****************************************************************      
      USE hdf5 
      USE m_hdf_tools 
      IMPLICIT NONE 
!     Arguments                                                         
      TYPE(t_rsdata),INTENT(IN) ::rs 
      INTEGER(HID_T),INTENT(IN) ::fid 
      CHARACTER*(*),INTENT(IN)  ::dataname 
                                                                        
      INTEGER(HID_T)::varid 
      INTEGER       ::hdferr 
                                                                        
      IF (io_dataexists(fid,dataname)) THEN 
         CALL io_dopen(fid,dataname,varid,hdferr) 
      ELSE 
         ! Write the dims                                               
         CALL io_createvar(fid,dataname,H5T_NATIVE_DOUBLE,ABS(rs%grid)  &
     &        ,varid)                                                   
         CALL io_WRITE_att(varid,'GRID',rs%grid) 
         CALL io_WRITE_att(varid,'BOX',rs%box) 
         CALL io_WRITE_Att(varid,'Z_ZERO',rs%z_zero) 
      ENDIF 
      !write data                                                       
      CALL io_write(varid,(/1,1,1/),abs(rs%grid),rs%data) 
      !close                                                            
      CALL io_dclose(varid,hdferr) 
      END SUBROUTINE 
                                                                        
      !>                                                                
      !<--S: gf_rssmoothz(rs,zstart,zstop)                              
      SUBROUTINE gf_rssmoothz(rs,zstart,zstop) 
!*****************************************************************      
! DESC: Smoothes the z-dependend data between zstart and zstop by using
!       polynomial 4-th order instead of the real data                  
!                                                                       
!                        Daniel Wortmann                                
!*****************************************************************      
      IMPLICIT NONE 
      INTEGER,INTENT(IN)             ::zstart,zstop
      TYPE(t_rsdata),INTENT(INOUT)   ::rs 
                                                                        
      INTEGER::x,y,zi 
      REAL   ::a1,a2,a3,a4,z 
      REAL   ::z1,z2,z3,z4 
                                                                        
                                                 !No smoothing because r
      IF (zstart<2.OR.zstop>rs%grid(3)-2) RETURN 
                                                 !to small!!            
      z1=zstart-1 
      z2=zstart 
      z3=zstop 
      z4=zstop+1 
      DO x=1,rs%grid(1) 
         DO y=1,rs%grid(2) 
            a1=rs%data(x,y,zstart-1) 
            a2=rs%data(x,y,zstart) 
            a3=rs%data(x,y,zstop) 
            a4=rs%data(x,y,zstop+1) 
            DO zi=zstart,zstop 
               z=1.0*zi 
               rs%data(x,y,zi)=                                         &
     &              (z-z2)*(z-z3)*(z-z4)/(z1-z2)/(z1-z3)/(z1-z4)*a1     &
     &              +(z-z1)*(z-z3)*(z-z4)/(z2-z1)/(z2-z3)/(z2-z4)*a2    &
     &              +(z-z1)*(z-z2)*(z-z4)/(z3-z1)/(z3-z2)/(z3-z4)*a3    &
     &              +(z-z1)*(z-z2)*(z-z3)/(z4-z1)/(z4-z2)/(z4-z3)*a4    
            ENDDO 
         ENDDO 
      ENDDO 
                                                                        
      END SUBROUTINE 
      !>                                                                
      !<--S: gf_rsread_filename(rs,filename,dataname)                   
      SUBROUTINE gf_rsread_filename(rs,filename,dataname) 
!*****************************************************************      
! DESC: Reads the rs-data from the hdf5-file specified by filename      
!       The dataset and its dims will be named by dataname              
!       THE RS-DATA SHOULD NOT BE ALLOCATED!                            
!                                                                       
!                        Daniel Wortmann, Fri Oct 18 10:53:49 2002      
!*****************************************************************      
      USE hdf5 
      USE m_hdf_tools 
      IMPLICIT NONE 
!     Arguments                                                         
      TYPE(t_rsdata),INTENT(OUT) ::rs 
      CHARACTER*(*),INTENT(IN)  ::filename,dataname 
                                                                        
      INTEGER(HID_T)::fid 
      INTEGER       ::hdferr 
      CALL io_hdfopen (filename//'.hdf', H5F_ACC_RDONLY_F,fid,hdferr) 
      CALL gf_rsread_fid(rs,fid,dataname) 
      CALL io_hdfclose(fid,hdferr) 
      END SUBROUTINE 
      !>                                                                
      !<--S: gf_rsread_fid(rs,filename,dataname)                        
      SUBROUTINE gf_rsread_fid(rs,fid,dataname) 
!*****************************************************************      
! DESC: Reads the rs-data from the hdf5-file specified by filename      
!       The dataset and its dims will be named by dataname              
!       THE RS-DATA SHOULD NOT BE ALLOCATED!                            
!                                                                       
!                        Daniel Wortmann, Fri Oct 18 10:53:49 2002      
!*****************************************************************      
      USE hdf5 
      USE m_hdf_tools 
      IMPLICIT NONE 
!     Arguments                                                         
      TYPE(t_rsdata),INTENT(OUT) :: rs 
      INTEGER(HID_T),INTENT(IN)  :: fid 
      CHARACTER*(*),INTENT(IN)   :: dataname 
                                                                        
      INTEGER(HID_T)::varid 
      INTEGER       ::hdferr 
      ! Get the dims                                                    
      CALL io_dopen(fid, dataname, varid, hdferr) 
      CALL io_read_att(varid,'GRID',rs%grid) 
      !Get the Box                                                      
      CALL io_read_att(varid,'BOX',rs%box) 
      CALL io_read_att(varid,'Z_ZERO',rs%z_zero) 
      !Read the Variable                                                
      ALLOCATE(rs%data(rs%grid(1),rs%grid(2),ABS(rs%grid(3)))) 
      CALL io_read(varid,(/1,1,1/),abs(rs%grid),rs%data) 
      !done                                                             
      CALL io_dclose(varid,hdferr) 
      END SUBROUTINE 
      !>                                                                
      !<--F: priv_interpolate(data,x,y,z) result(v)                     
      FUNCTION interpolate(data,x,y,z) RESULT(v) 
!*****************************************************************      
! DESC: private function to interpolate data                            
!*****************************************************************      
      IMPLICIT NONE 
!     Arguments                                                         
      REAL,INTENT(IN)::data(:,:,:),x,y,z 
      REAL::v 
                                                                        
                                                                        
      REAL::t,u,v1,v2 
      INTEGER::x1,x2,y1,y2,z1,z2 
      ! First Do two dim interpolations in the xy plane, than           
      ! 1D interpolation in z-direction                                 
                                                                        
      t=(x-floor(x)) 
      u=(y-floor(y)) 
                                                                        
      x1=floor(x) 
      IF (x1<size(data,1)) THEN 
         x2=x1+1 
      ELSE 
         x2=1 
      ENDIF 
      y1=floor(y) 
      IF (y1<size(data,2)) THEN 
         y2=y1+1 
      ELSE 
         y2=1 
      ENDIF 
      z1=floor(z) 
      IF (z1<size(data,3)) THEN 
         z2=z1+1 
      ELSE 
         z2=1 
      ENDIF 
                                                                        
      v1=(1-t)*(1-u)*data(x1,y1,z1)+t*(1-u)*data(x2,y1,z1)+t*u*data(x2  &
     &     ,y2,z1)*(1-t)*u*data(x1,y2,z1)                               
      v2=(1-t)*(1-u)*data(x1,y1,z2)+t*(1-u)*data(x2,y1,z2)+t*u*data(x2  &
     &     ,y2,z2)*(1-t)*u*data(x1,y2,z2)                               
                                                                        
      t=z-floor(z) 
                                                                        
      v=t*v2+(1-t)*v1 
                                                                        
                                                                        
                                                                        
      END FUNCTION 
      !>                                                                
      !<-- S: gf_rscutbox(rs_in,zmin,zmax,rs_out)                       
      SUBROUTINE gf_rscutbox(rs_in,zmin,zmax,rs_out) 
!-----------------------------------------------                        
!    Constructs a subbox containing the data between zmin and zmax      
!    the origin of the box will be at zmin, e.g. the z_zero will        
!    be set to the distance of the first point to zmin                  
!           (last modified: 05-05-11) D. Wortmann                       
!-----------------------------------------------                        
      IMPLICIT NONE 
      !<--Arguments                                                     
      TYPE(t_rsdata),INTENT(IN)  :: rs_in 
      REAL   ,INTENT(IN)         :: zmin,zmax 
      TYPE(t_rsdata),INTENT(OUT) :: rs_out 
      !>                                                                
      !<-- Locals                                                       
      INTEGER             :: i1,i2 
                                                                        
      !>                                                                
      rs_out%grid(1:2)=rs_in%grid(1:2) 
      rs_out%box(1:2)=rs_out%box(1:2) 
                                                                        
      IF (zmax<zmin) CALL juDFT_error("gf_rscutbox: zmax<zmin") 
      IF (zmin<0) THEN 
         IF (zmax>0) CALL juDFT_error("gf_rscutbox cannot cut this") 
         i1 = FLOOR((rs_in%box(3)+zmin)/rs_in%box(3)*rs_in%grid(3)) 
         i2 = FLOOR((rs_in%box(3)+zmax)/rs_in%box(3)*rs_in%grid(3))+1 
         rs_out%z_zero = i1*rs_in%box(3)/rs_in%grid(3)-(zmin+rs_in%box(3&
     &        ))                                                        
      ELSE 
         i1 = FLOOR(zmin/rs_in%box(3)*rs_in%grid(3)) 
         i2 = FLOOR(zmax/rs_in%box(3)*rs_in%grid(3))+1 
         rs_out%z_zero = i1*rs_in%box(3)/rs_in%grid(3)-zmin 
      ENDIF 
      rs_out%grid(3)=i2-i1+1 
      rs_out%box(3)=rs_out%grid(3)*rs_in%box(3)/rs_in%grid(3) 
                                                                        
      ALLOCATE(rs_out%data(rs_out%grid(1),rs_out%grid(2),rs_out%grid(3))&
     &     )                                                            
                                                                        
      rs_out%data(:,:,:) = rs_in%data(:,:,i1:i2) 
                                                                        
                                                                        
      END SUBROUTINE 
      !>                                                                
      !<-- S: gf_rsputbox(rs_b,zmin,rs)                                 
      SUBROUTINE gf_rsputbox(rs_b,zmin,rs) 
!-----------------------------------------------                        
!     puts the data of rs_box into rs                                   
!           (last modified: 2004-00-00) D. Wortmann                     
!-----------------------------------------------                        
      IMPLICIT NONE 
      !<--Arguments                                                     
      TYPE(t_rsdata),INTENT(IN)    :: rs_b 
      REAL                         :: zmin 
      TYPE(t_rsdata),INTENT(INOUT) :: rs 
      !>                                                                
      !<-- Locals                                                       
      INTEGER             :: i1,i2,nz,n 
      REAL                :: z,dist,z0,zmesh(rs_b%grid(3)) 
                                                                        
      !>                                                                
      IF (zmin<0) THEN 
         i1 = FLOOR((rs%box(3)+zmin)/rs%box(3)*rs%grid(3)) 
         i2 = FLOOR((rs%box(3)+MIN(0.0,zmin+rs_b%box(3)))/rs%box(3)     &
     &        *rs%grid(3))+1                                            
         i2 = MAX(i2,rs%grid(3)) 
         z0 = rs%box(3)+zmin 
      ELSE 
         i1 = FLOOR(zmin/rs%box(3)*rs%grid(3)) 
         i2 = FLOOR((zmin+rs_b%box(3))/rs%box(3)*rs%grid(3))+1 
         i2 = MIN(i2,int(rs%grid(3)/2)) 
         z0 = zmin 
      ENDIF 
      DO n = 1,rs_b%grid(3) 
         zmesh(n) = rs_b%z_zero+(n-1)*rs_b%box(3)/rs_b%grid(3) 
      ENDDO 
      DO n = i1,i2 
         z = n*rs%box(3)/rs%grid(3)-z0 
         DO nz = 1,rs_b%grid(3) 
            IF (zmesh(nz)>=z) EXIT 
         ENDDO 
         IF (nz == 1) THEN 
            rs%data(:,:,n)=rs_b%data(:,:,1) 
         ELSE IF (nz>rs_b%grid(3)) THEN 
            rs%data(:,:,n)=rs_b%data(:,:,rs_b%grid(3)) 
         ELSE 
            dist = (z-zmesh(nz-1))/(zmesh(nz)-zmesh(nz-1)) 
            rs%data(:,:,n) = rs_b%data(:,:,nz-1)+dist*(rs_b%data(:,:,nz)&
     &           - rs_b%data(:,:,nz-1))                                 
         ENDIF 
      ENDDO 
                                                                        
      END SUBROUTINE 
      !>                                                                
      !<-- S: gf_rsmergers(rs1,rs2,step,rs3)                            
      SUBROUTINE  gf_rsmergers(rs1,rs2,step,rs3) 
!-----------------------------------------------                        
!                                                                       
!           (last modified: 2004-00-00) D. Wortmann                     
!-----------------------------------------------                        
      IMPLICIT NONE 
      !<--Arguments                                                     
      TYPE(t_rsdata),INTENT(IN)  :: rs1,rs2 
      LOGICAL,INTENT(IN)         :: step(:) 
      TYPE(t_rsdata),INTENT(OUT) :: rs3 
      !>                                                                
      IF (ANY(rs1%grid-rs2%grid/=0).OR.ANY(ABS(rs1%box-rs2%box)>1E-6))  &
     &     CALL juDFT_error("Different rs-data in gf_rsmergers")          
      IF (SIZE(rs1%data) /= SIZE(step))                                 &
     &     CALL juDFT_error("Wrong dimension of step in gf_rsmergers")
      rs3%grid = rs1%grid 
      rs3%box  = rs3%box 
      rs3%z_zero = rs1%z_zero 
      ALLOCATE(rs3%data(rs3%grid(1),rs3%grid(2),rs3%grid(3))) 
      rs3%data = MERGE(rs1%data,rs2%data,RESHAPE(step,rs3%grid)) 
                                                                        
      END SUBROUTINE 
      !>                                                                
      END                                           
