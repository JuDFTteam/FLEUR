!--------------------------------------------------------------------------------
! Copyright (c) 2016 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------

      MODULE m_hdf_tools6 
!-----------------------------------------------                        
!     major rewrite of hdf_tools                                        
!     this module contains subroutines that create a                    
!     variable for an array and fill it with data at once               
!-----------------------------------------------                        
      USE hdf5 
      PRIVATE 
      INTERFACE io_WRITE_var 
      MODULE PROCEDURE io_WRITE_var_REAL1,io_WRITE_var_REAL2            &
     &     ,io_WRITE_var_REAL3,io_WRITE_var_REAL4,io_WRITE_var_integer1 &
     &     ,io_WRITE_var_integer2,io_WRITE_var_integer3                 &
     &     ,io_WRITE_var_integer4                                       
      END INTERFACE 
      INTERFACE io_READ_var 
      MODULE PROCEDURE io_READ_var_REAL1,io_READ_var_REAL2              &
     &     ,io_READ_var_REAL3,io_READ_var_REAL4,io_READ_var_integer1    &
     &     ,io_READ_var_integer2,io_READ_var_integer3                   &
     &     ,io_READ_var_integer4                                        
      END INTERFACE 
      PUBLIC :: io_WRITE_var,io_read_var 
      CONTAINS 
                                                                        
!-----------------------------------------------                        
! real-routines for writing                                             
!           (last modified: 07-07-31) D. Wortmann                       
!-----------------------------------------------                        
      !<-- S:io_WRITE_var_REAL1                                         
                                                                        
      SUBROUTINE io_WRITE_var_REAL1(gid,name,var,transprop) 
                                                                        
      USE m_hdf_tools4,ONLY:hdf_err,rkind,io_createvar 
      USE m_hdf_tools2,ONLY:io_WRITE 
      USE m_hdf_tools3,ONLY:io_dataexists 
                                                                        
      IMPLICIT NONE 
      !<-- Arguments                                                    
      INTEGER(HID_T),INTENT(IN) :: gid 
      CHARACTER*(*)              :: name 
      REAL(rkind),INTENT(IN)     :: var(:) 
      INTEGER(HID_T),INTENT(IN),OPTIONAL :: transprop 
      !>                                                                
      INTEGER(HID_T) :: varid 
      INTEGER        :: hdferr 
                                                                        
      IF (io_dataexists(gid,name)) call hdf_err                         &
     &     ("Variable could not be created:"//name)                     
      CALL io_createvar(gid,name, H5T_NATIVE_DOUBLE,SHAPE(var),varid) 
      IF (PRESENT(transprop)) THEN 
         CALL io_WRITE(varid,(/1/),SHAPE(var),var,transprop) 
      ELSE 
         CALL io_WRITE(varid,(/1/),SHAPE(var),var) 
      ENDIF 
                                                                        
      CALL h5dclose_f(varid,hdferr) 
                                                                        
      END SUBROUTINE 
                                                                        
      !>                                                                
      !<-- S:io_WRITE_var_REAL2                                         
                                                                        
      SUBROUTINE io_WRITE_var_REAL2(gid,name,var,transprop) 
                                                                        
      USE m_hdf_tools4,ONLY:hdf_err,rkind,io_createvar 
      USE m_hdf_tools2,ONLY:io_WRITE 
      USE m_hdf_tools3,ONLY:io_dataexists 
                                                                        
      IMPLICIT NONE 
      !<-- Arguments                                                    
      INTEGER(HID_T),INTENT(IN) :: gid 
      CHARACTER*(*)              :: name 
      REAL(rkind),INTENT(IN)     :: var(:,:) 
      INTEGER(HID_T),INTENT(IN),OPTIONAL :: transprop 
      !>                                                                
      INTEGER(HID_T) :: varid 
      INTEGER        :: hdferr 
                                                                        
      IF (io_dataexists(gid,name)) call hdf_err                         &
     &     ("Variable could not be created:"//name)                     
      CALL io_createvar(gid,name, H5T_NATIVE_DOUBLE,SHAPE(var),varid) 
      IF (PRESENT(transprop)) THEN 
         CALL io_WRITE(varid,(/1,1/),SHAPE(var),var,transprop) 
      ELSE 
         CALL io_WRITE(varid,(/1,1/),SHAPE(var),var) 
      ENDIF 
      CALL h5dclose_f(varid,hdferr) 
                                                                        
      END SUBROUTINE 
                                                                        
      !>                                                                
      !<-- S:io_WRITE_var_REAL3                                         
                                                                        
      SUBROUTINE io_WRITE_var_REAL3(gid,name,var,transprop) 
                                                                        
      use m_hdf_tools4,ONLY:hdf_err,rkind,io_createvar 
      use m_hdf_tools2,ONLY:io_WRITE 
      use m_hdf_tools3,ONLY:io_dataexists 
                                                                        
      IMPLICIT NONE 
      !<-- Arguments                                                    
      INTEGER(HID_T),INTENT(IN) :: gid 
      CHARACTER*(*)              :: name 
      REAL(rkind),INTENT(IN)     :: var(:,:,:) 
      INTEGER(HID_T),INTENT(IN),OPTIONAL :: transprop 
      !>                                                                
      INTEGER(HID_T) :: varid 
      INTEGER        :: hdferr 
                                                                        
      IF (io_dataexists(gid,name)) call hdf_err                         &
     &     ("Variable could not be created:"//name)                     
      CALL io_createvar(gid,name, H5T_NATIVE_DOUBLE,SHAPE(var),varid) 
                                                                        
      IF (PRESENT(transprop)) THEN 
         CALL io_WRITE(varid,(/1,1,1/),SHAPE(var),var,transprop) 
      ELSE 
         CALL io_WRITE(varid,(/1,1,1/),SHAPE(var),var) 
      ENDIF 
      CALL h5dclose_f(varid,hdferr) 
                                                                        
      END SUBROUTINE 
                                                                        
      !>                                                                
      !<-- S:io_WRITE_var_REAL4                                         
                                                                        
      SUBROUTINE io_WRITE_var_REAL4(gid,name,var,transprop) 
                                                                        
      use m_hdf_tools4,ONLY:hdf_err,rkind,io_createvar 
      use m_hdf_tools2,ONLY:io_WRITE 
      use m_hdf_tools3,ONLY:io_dataexists 
                                                                        
      IMPLICIT NONE 
      !<-- Arguments                                                    
      INTEGER(HID_T),INTENT(IN) :: gid 
      CHARACTER*(*)              :: name 
      REAL(rkind),INTENT(IN)     :: var(:,:,:,:) 
      INTEGER(HID_T),INTENT(IN),OPTIONAL :: transprop 
      !>                                                                
      INTEGER(HID_T) :: varid 
      INTEGER        :: hdferr 
                                                                        
      IF (io_dataexists(gid,name)) call hdf_err                         &
     &     ("Variable could not be created:"//name)                     
      CALL io_createvar(gid,name, H5T_NATIVE_DOUBLE,SHAPE(var),varid) 
      IF (PRESENT(transprop)) THEN 
         CALL io_WRITE(varid,(/1,1,1,1/),SHAPE(var),var,transprop) 
      ELSE 
         CALL io_WRITE(varid,(/1,1,1,1/),SHAPE(var),var) 
      ENDIF 
      CALL h5dclose_f(varid,hdferr) 
                                                                        
      END SUBROUTINE 
                                                                        
      !>                                                                
                                                                        
!-----------------------------------------------                        
! integer-routines for writing                                          
!           (last modified: 07-07-31) D. Wortmann                       
!-----------------------------------------------                        
      !<-- S:io_WRITE_var_INTEGER1                                      
                                                                        
      SUBROUTINE io_WRITE_var_integer1(gid,name,var,transprop) 
                                                                        
      use m_hdf_tools4,ONLY:hdf_err,rkind,io_createvar 
      use m_hdf_tools2,ONLY:io_WRITE 
      use m_hdf_tools3,ONLY:io_dataexists 
                                                                        
      IMPLICIT NONE 
      !<-- Arguments                                                    
      INTEGER(HID_T),INTENT(IN)  :: gid 
      CHARACTER*(*)              :: name 
      INTEGER,INTENT(IN)     :: var(:) 
      INTEGER(HID_T),INTENT(IN),OPTIONAL::transprop 
                                                                        
      !>                                                                
      INTEGER(HID_T) :: varid 
      INTEGER        :: hdferr 
                                                                        
      IF (io_dataexists(gid,name)) call hdf_err                         &
     &     ("Variable could not be created:"//name)                     
      CALL io_createvar(gid,name, H5T_NATIVE_INTEGER,SHAPE(var),varid) 
      IF (PRESENT(transprop)) THEN 
         CALL io_WRITE(varid,(/1/),SHAPE(var),var,transprop) 
      ELSE 
         CALL io_WRITE(varid,(/1/),SHAPE(var),var) 
      ENDIF 
      CALL h5dclose_f(varid,hdferr) 
                                                                        
      END SUBROUTINE 
                                                                        
      !>                                                                
      !<-- S:io_WRITE_var_INTEGER2                                      
                                                                        
      SUBROUTINE io_WRITE_var_INTEGER2(gid,name,var,transprop) 
                                                                        
      use m_hdf_tools4,ONLY:hdf_err,rkind,io_createvar 
      use m_hdf_tools2,ONLY:io_WRITE 
      use m_hdf_tools3,ONLY:io_dataexists 
                                                                        
      IMPLICIT NONE 
      !<-- Arguments                                                    
      INTEGER(HID_T),INTENT(IN)  :: gid 
      CHARACTER*(*)              :: name 
      INTEGER,INTENT(IN)     :: var(:,:) 
      INTEGER(HID_T),INTENT(IN),OPTIONAL::transprop 
      !>                                                                
      INTEGER(HID_T) :: varid 
      INTEGER        :: hdferr 
                                                                        
      IF (io_dataexists(gid,name)) call hdf_err                         &
     &     ("Variable could not be created:"//name)                     
      CALL io_createvar(gid,name, H5T_NATIVE_INTEGER,SHAPE(var),varid) 
      IF (PRESENT(transprop)) THEN 
         CALL io_WRITE(varid,(/1,1/),SHAPE(var),var,transprop) 
      ELSE 
         CALL io_WRITE(varid,(/1,1/),SHAPE(var),var) 
      ENDIF 
      CALL h5dclose_f(varid,hdferr) 
                                                                        
      END SUBROUTINE 
                                                                        
      !>                                                                
      !<-- S:io_WRITE_var_INTEGER3                                      
                                                                        
      SUBROUTINE io_WRITE_var_INTEGER3(gid,name,var,transprop) 
                                                                        
      use m_hdf_tools4,ONLY:hdf_err,rkind,io_createvar 
      use m_hdf_tools2,ONLY:io_WRITE 
      use m_hdf_tools3,ONLY:io_dataexists 
                                                                        
      IMPLICIT NONE 
      !<-- Arguments                                                    
      INTEGER(HID_T),INTENT(IN) :: gid 
      CHARACTER*(*)              :: name 
      INTEGER,INTENT(IN)     :: var(:,:,:) 
      INTEGER(HID_T),INTENT(IN),OPTIONAL :: transprop 
      !>                                                                
      INTEGER(HID_T) :: varid 
      INTEGER        :: hdferr 
                                                                        
      IF (io_dataexists(gid,name)) call hdf_err                         &
     &     ("Variable could not be created:"//name)                     
      CALL io_createvar(gid,name, H5T_NATIVE_INTEGER,SHAPE(var),varid) 
      IF (PRESENT(transprop)) THEN 
         CALL io_WRITE(varid,(/1,1,1/),SHAPE(var),var,transprop) 
      ELSE 
         CALL io_WRITE(varid,(/1,1,1/),SHAPE(var),var) 
      ENDIF 
      CALL h5dclose_f(varid,hdferr) 
                                                                        
      END SUBROUTINE 
                                                                        
      !>                                                                
      !<-- S:io_WRITE_var_INTEGER4                                      
                                                                        
      SUBROUTINE io_WRITE_var_INTEGER4(gid,name,var,transprop) 
                                                                        
      use m_hdf_tools4,ONLY:hdf_err,rkind,io_createvar 
      use m_hdf_tools2,ONLY:io_WRITE 
      use m_hdf_tools3,ONLY:io_dataexists 
                                                                        
      IMPLICIT NONE 
      !<-- Arguments                                                    
      INTEGER(HID_T),INTENT(IN) :: gid 
      CHARACTER*(*)              :: name 
      INTEGER,INTENT(IN)     :: var(:,:,:,:) 
      INTEGER(HID_T),INTENT(IN),OPTIONAL :: transprop 
      !>                                                                
      INTEGER(HID_T) :: varid 
      INTEGER        :: hdferr 
                                                                        
      IF (io_dataexists(gid,name)) call hdf_err                         &
     &     ("Variable could not be created:"//name)                     
      CALL io_createvar(gid,name, H5T_NATIVE_INTEGER,SHAPE(var),varid) 
      IF (PRESENT(transprop)) THEN 
         CALL io_WRITE(varid,(/1,1,1,1/),SHAPE(var),var,transprop) 
      ELSE 
         CALL io_WRITE(varid,(/1,1,1,1/),SHAPE(var),var) 
      ENDIF 
      CALL h5dclose_f(varid,hdferr) 
                                                                        
      END SUBROUTINE 
                                                                        
      !>                                                                
                                                                        
!-----------------------------------------------                        
!  real :: routines for reading                                         
!           (last modified: 07-07-31) D. Wortmann                       
!-----------------------------------------------                        
                                                                        
      !<-- S:io_read_var_real1                                          
                                                                        
      SUBROUTINE io_READ_var_real1(gid,name,var,transprop) 
      use m_hdf_tools4,ONLY:hdf_err,rkind 
      use m_hdf_tools2,ONLY:io_read 
      use m_hdf_tools3,ONLY:io_dataexists 
      IMPLICIT NONE 
      !<--Arguments                                                     
      INTEGER(HID_T),INTENT(IN)  :: gid 
      CHARACTER*(*)              :: name 
      REAL(rkind),INTENT(OUT)    :: var(:) 
      INTEGER(HID_T),INTENT(IN),OPTIONAL::transprop 
      !>                                                                
      !<-- Locals                                                       
      INTEGER(HID_T) :: varid,fspace 
      INTEGER        :: hdferr 
      INTEGER(HSIZE_T):: dims(7) 
      INTEGER(HSIZE_T):: maxdims(7) 
      !>                                                                
      IF (.NOT.io_dataexists(gid,name)) CALL hdf_err                    &
     &     ("Variable could not be read:"//name)                        
                                                                        
      CALL h5dopen_f(gid, name, varid, hdferr) 
      CALL h5dget_space_f(varid,fspace,hdferr) 
      CALL H5Sget_simple_extent_dims_f(fspace,dims,maxdims,hdferr) 
      IF (ANY(dims(:1)>SHAPE(var))) CALL hdf_err                        &
     &     ("Not enough data in:"//name)                                
      CALL h5sclose_f(fspace,hdferr) 
      IF (PRESENT(transprop)) THEN 
         CALL io_read(varid,(/1/),SHAPE(var),var,transprop) 
      ELSE 
         CALL io_read(varid,(/1/),SHAPE(var),var) 
      ENDIF 
      CALL h5dclose_f(varid,hdferr) 
                                                                        
      END SUBROUTINE 
                                                                        
      !>                                                                
      !<-- S:io_read_var_real2                                          
                                                                        
      SUBROUTINE io_READ_var_real2(gid,name,var,transprop) 
      use m_hdf_tools4,ONLY:hdf_err,rkind 
      use m_hdf_tools2,ONLY:io_read 
      use m_hdf_tools3,ONLY:io_dataexists 
      IMPLICIT NONE 
      !<--Arguments                                                     
      INTEGER(HID_T),INTENT(IN)  :: gid 
      CHARACTER*(*)              :: name 
      REAL(rkind),INTENT(OUT)    :: var(:,:) 
      INTEGER(HID_T),INTENT(IN),OPTIONAL::transprop 
      !>                                                                
      !<-- Locals                                                       
      INTEGER(HID_T) :: varid,fspace 
      INTEGER        :: hdferr 
      INTEGER(HSIZE_T):: dims(7) 
      INTEGER(HSIZE_T):: maxdims(7) 
      !>                                                                
      IF (.NOT.io_dataexists(gid,name)) CALL hdf_err                    &
     &     ("Variable could not be read:"//name)                        
                                                                        
      CALL h5dopen_f(gid, name, varid, hdferr) 
      CALL h5dget_space_f(varid,fspace,hdferr) 
      CALL H5Sget_simple_extent_dims_f(fspace,dims,maxdims,hdferr) 
      IF (ANY(dims(:2)>SHAPE(var))) CALL hdf_err                        &
     &     ("Not enough data in:"//name)                                
      CALL h5sclose_f(fspace,hdferr) 
      IF (PRESENT(transprop)) THEN 
         CALL io_read(varid,(/1,1/),SHAPE(var),var,transprop) 
      ELSE 
         CALL io_read(varid,(/1,1/),SHAPE(var),var) 
      ENDIF 
                                                                        
      CALL h5dclose_f(varid,hdferr) 
                                                                        
      END SUBROUTINE 
                                                                        
      !>                                                                
      !<-- S:io_read_var_real3                                          
                                                                        
      SUBROUTINE io_READ_var_real3(gid,name,var,transprop) 
      use m_hdf_tools4,ONLY:hdf_err,rkind 
      use m_hdf_tools2,ONLY:io_read 
      use m_hdf_tools3,ONLY:io_dataexists 
      IMPLICIT NONE 
      !<--Arguments                                                     
      INTEGER(HID_T),INTENT(IN)  :: gid 
      CHARACTER*(*)              :: name 
      REAL(rkind),INTENT(OUT)    :: var(:,:,:) 
      INTEGER(HID_T),INTENT(IN),OPTIONAL::transprop 
      !>                                                                
      !<-- Locals                                                       
      INTEGER(HID_T) :: varid,fspace 
      INTEGER        :: hdferr 
      INTEGER(HSIZE_T):: dims(7) 
      INTEGER(HSIZE_T):: maxdims(7) 
      !>                                                                
      IF (.NOT.io_dataexists(gid,name)) CALL hdf_err                    &
     &     ("Variable could not be read:"//name)                        
                                                                        
      CALL h5dopen_f(gid, name, varid, hdferr) 
      CALL h5dget_space_f(varid,fspace,hdferr) 
      CALL H5Sget_simple_extent_dims_f(fspace,dims,maxdims,hdferr) 
      IF (ANY(dims(:3)>SHAPE(var))) CALL hdf_err                        &
     &     ("Not enough data in:"//name)                                
      CALL h5sclose_f(fspace,hdferr) 
      IF (PRESENT(transprop)) THEN 
         CALL io_read(varid,(/1,1,1/),SHAPE(var),var,transprop) 
      ELSE 
         CALL io_read(varid,(/1,1,1/),SHAPE(var),var) 
      ENDIF 
                                                                        
                                                                        
      CALL h5dclose_f(varid,hdferr) 
                                                                        
      END SUBROUTINE 
                                                                        
      !>                                                                
      !<-- S:io_read_var_real4                                          
                                                                        
      SUBROUTINE io_READ_var_real4(gid,name,var,transprop) 
      use m_hdf_tools4,ONLY:hdf_err,rkind 
      use m_hdf_tools2,ONLY:io_read 
      use m_hdf_tools3,ONLY:io_dataexists 
      IMPLICIT NONE 
      !<--Arguments                                                     
      INTEGER(HID_T),INTENT(IN)  :: gid 
      CHARACTER*(*)              :: name 
      REAL(rkind),INTENT(OUT)    :: var(:,:,:,:) 
      INTEGER(HID_T),INTENT(IN),OPTIONAL::transprop 
      !>                                                                
      !<-- Locals                                                       
      INTEGER(HID_T) :: varid,fspace 
      INTEGER        :: hdferr 
      INTEGER(HSIZE_T):: dims(7) 
      INTEGER(HSIZE_T):: maxdims(7) 
      !>                                                                
      IF (.NOT.io_dataexists(gid,name)) CALL hdf_err                    &
     &     ("Variable could not be read:"//name)                        
                                                                        
      CALL h5dopen_f(gid, name, varid, hdferr) 
      CALL h5dget_space_f(varid,fspace,hdferr) 
      CALL H5Sget_simple_extent_dims_f(fspace,dims,maxdims,hdferr) 
      IF (ANY(dims(:4)>SHAPE(var))) CALL hdf_err                        &
     &     ("Not enough data in:"//name)                                
      CALL h5sclose_f(fspace,hdferr) 
            IF (PRESENT(transprop)) THEN 
         CALL io_read(varid,(/1,1,1,1/),SHAPE(var),var,transprop) 
      ELSE 
         CALL io_read(varid,(/1,1,1,1/),SHAPE(var),var) 
      ENDIF 
                                                                        
                                                                        
      CALL h5dclose_f(varid,hdferr) 
                                                                        
      END SUBROUTINE 
                                                                        
      !>                                                                
                                                                        
!-----------------------------------------------                        
!  integer :: routines for reading                                      
!           (last modified: 07-07-31) D. Wortmann                       
!-----------------------------------------------                        
                                                                        
      !<-- S:io_read_var_integer1                                       
                                                                        
      SUBROUTINE io_READ_var_integer1(gid,name,var,transprop) 
      use m_hdf_tools4,ONLY:hdf_err,rkind 
      use m_hdf_tools2,ONLY:io_read 
      use m_hdf_tools3,ONLY:io_dataexists 
      IMPLICIT NONE 
      !<--Arguments                                                     
      INTEGER(HID_T),INTENT(IN)  :: gid 
      CHARACTER*(*)              :: name 
      integer,INTENT(OUT)    :: var(:) 
      INTEGER(HID_T),INTENT(IN),OPTIONAL::transprop 
      !>                                                                
      !<-- Locals                                                       
      INTEGER(HID_T) :: varid,fspace 
      INTEGER        :: hdferr 
      INTEGER(HSIZE_T):: dims(7) 
      INTEGER(HSIZE_T):: maxdims(7) 
      !>                                                                
      IF (.NOT.io_dataexists(gid,name)) CALL hdf_err                    &
     &     ("Variable could not be read:"//name)                        
                                                                        
      CALL h5dopen_f(gid, name, varid, hdferr) 
      CALL h5dget_space_f(varid,fspace,hdferr) 
      CALL H5Sget_simple_extent_dims_f(fspace,dims,maxdims,hdferr) 
      IF (ANY(dims(:1)>SHAPE(var))) CALL hdf_err                        &
     &     ("Not enough data in:"//name)                                
      CALL h5sclose_f(fspace,hdferr) 
      IF (PRESENT(transprop)) THEN 
         CALL io_read(varid,(/1/),SHAPE(var),var,transprop) 
      ELSE 
         CALL io_read(varid,(/1/),SHAPE(var),var) 
      ENDIF 
                                                                        
                                                                        
      CALL h5dclose_f(varid,hdferr) 
                                                                        
      END SUBROUTINE 
                                                                        
      !>                                                                
      !<-- S:io_read_var_integer2                                       
                                                                        
      SUBROUTINE io_READ_var_integer2(gid,name,var,transprop) 
      use m_hdf_tools4,ONLY:hdf_err,rkind 
      use m_hdf_tools2,ONLY:io_read 
      use m_hdf_tools3,ONLY:io_dataexists 
      IMPLICIT NONE 
      !<--Arguments                                                     
      INTEGER(HID_T),INTENT(IN)  :: gid 
      CHARACTER*(*)              :: name 
      integer,INTENT(OUT)    :: var(:,:) 
      INTEGER(HID_T),INTENT(IN),OPTIONAL::transprop 
      !>                                                                
      !<-- Locals                                                       
      INTEGER(HID_T) :: varid,fspace 
      INTEGER        :: hdferr 
      INTEGER(HSIZE_T):: dims(7) 
      INTEGER(HSIZE_T):: maxdims(7) 
      !>                                                                
      IF (.NOT.io_dataexists(gid,name)) CALL hdf_err                    &
     &     ("Variable could not be read:"//name)                        
                                                                        
      CALL h5dopen_f(gid, name, varid, hdferr) 
      CALL h5dget_space_f(varid,fspace,hdferr) 
      CALL H5Sget_simple_extent_dims_f(fspace,dims,maxdims,hdferr) 
      IF (ANY(dims(:2)>SHAPE(var))) CALL hdf_err                        &
     &     ("Not enough data in:"//name)                                
      CALL h5sclose_f(fspace,hdferr) 
      IF (PRESENT(transprop)) THEN 
         CALL io_read(varid,(/1,1/),SHAPE(var),var,transprop) 
      ELSE 
         CALL io_read(varid,(/1,1/),SHAPE(var),var) 
      ENDIF 
                                                                        
                                                                        
      CALL h5dclose_f(varid,hdferr) 
                                                                        
      END SUBROUTINE 
                                                                        
      !>                                                                
      !<-- S:io_read_var_integer3                                       
                                                                        
      SUBROUTINE io_READ_var_integer3(gid,name,var,transprop) 
      use m_hdf_tools4,ONLY:hdf_err,rkind 
      use m_hdf_tools2,ONLY:io_read 
      use m_hdf_tools3,ONLY:io_dataexists 
      IMPLICIT NONE 
      !<--Arguments                                                     
      INTEGER(HID_T),INTENT(IN)  :: gid 
      CHARACTER*(*)              :: name 
      integer,INTENT(OUT)    :: var(:,:,:) 
      INTEGER(HID_T),INTENT(IN),OPTIONAL::transprop 
      !>                                                                
      !<-- Locals                                                       
      INTEGER(HID_T) :: varid,fspace 
      INTEGER        :: hdferr 
      INTEGER(HSIZE_T):: dims(7) 
      INTEGER(HSIZE_T):: maxdims(7) 
      !>                                                                
      IF (.NOT.io_dataexists(gid,name)) CALL hdf_err                    &
     &     ("Variable could not be read:"//name)                        
                                                                        
      CALL h5dopen_f(gid, name, varid, hdferr) 
      CALL h5dget_space_f(varid,fspace,hdferr) 
      CALL H5Sget_simple_extent_dims_f(fspace,dims,maxdims,hdferr) 
      IF (ANY(dims(:3)>SHAPE(var))) CALL hdf_err                        &
     &     ("Not enough data in:"//name)                                
      CALL h5sclose_f(fspace,hdferr) 
      IF (PRESENT(transprop)) THEN 
         CALL io_read(varid,(/1,1,1/),SHAPE(var),var,transprop) 
      ELSE 
         CALL io_read(varid,(/1,1,1/),SHAPE(var),var) 
      ENDIF 
                                                                        
                                                                        
      CALL h5dclose_f(varid,hdferr) 
                                                                        
      END SUBROUTINE 
                                                                        
      !>                                                                
      !<-- S:io_read_var_integer4                                       
      SUBROUTINE io_READ_var_integer4(gid,name,var,transprop) 
      use m_hdf_tools4,ONLY:hdf_err,rkind 
      use m_hdf_tools2,ONLY:io_read 
      use m_hdf_tools3,ONLY:io_dataexists 
      IMPLICIT NONE 
      !<--Arguments                                                     
      INTEGER(HID_T),INTENT(IN)  :: gid 
      CHARACTER*(*)              :: name 
      integer,INTENT(OUT)    :: var(:,:,:,:) 
      INTEGER(HID_T),INTENT(IN),OPTIONAL::transprop 
      !>                                                                
      !<-- Locals                                                       
      INTEGER(HID_T) :: varid,fspace 
      INTEGER        :: hdferr 
      INTEGER(HSIZE_T):: dims(7) 
      INTEGER(HSIZE_T):: maxdims(7) 
      !>                                                                
      IF (.NOT.io_dataexists(gid,name)) CALL hdf_err                    &
     &     ("Variable could not be read:"//name)                        
                                                                        
      CALL h5dopen_f(gid, name, varid, hdferr) 
      CALL h5dget_space_f(varid,fspace,hdferr) 
      CALL H5Sget_simple_extent_dims_f(fspace,dims,maxdims,hdferr) 
      IF (ANY(dims(:4)>SHAPE(var))) CALL hdf_err                        &
     &     ("Not enough data in:"//name)                                
      CALL h5sclose_f(fspace,hdferr) 
      IF (PRESENT(transprop)) THEN 
         CALL io_read(varid,(/1,1,1,1/),SHAPE(var),var,transprop) 
      ELSE 
         CALL io_read(varid,(/1,1,1,1/),SHAPE(var),var) 
      ENDIF 
                                                                        
      CALL h5dclose_f(varid,hdferr) 
                                                                        
      END SUBROUTINE 
      !>                                                                
                                                                        
      END                                           
