      MODULE m_hdf_tools1 
!-----------------------------------------------                        
!    major rewrite of hdf_tools!                                        
!    The part contains the io_[read/write]_att part                     
!    public is only the generic interface                               
!                                                                       
!-----------------------------------------------                        
      USE m_hdf_tools4 
      !PRIVATE                                                          
      !<-- definitions of interfaces                                    
      INTERFACE io_write_att 
      MODULE PROCEDURE io_write_attreal0,io_write_attreal1              &
     &     ,io_write_attreal2,io_write_attreal3                         
      MODULE PROCEDURE io_write_attint0,io_write_attint1                &
     &     ,io_write_attint2,io_write_attint3                           
      MODULE PROCEDURE io_write_attlog0 
      MODULE PROCEDURE io_write_attchar0 
      END INTERFACE 
      INTERFACE io_read_att 
      MODULE PROCEDURE io_read_attreal0,io_read_attreal1                &
     &     ,io_read_attreal2,io_read_attreal3                           
      MODULE PROCEDURE io_read_attint0,io_read_attint1                  &
     &     ,io_read_attint2,io_read_attint3                             
      MODULE PROCEDURE io_read_attlog0 
      MODULE PROCEDURE io_read_attchar0 
      END INTERFACE 
      !>                                                                
      public:: io_read_att,io_write_att 
      CONTAINS 
      !<-- implementations of the io_write/read_att subroutines         
                                                                        
!*****************************************************************      
!                                                                       
!     The following subroutines READ or WRITE a CHARACTER*-attributes   
!                                                                       
!                                                                       
!*****************************************************************      
      SUBROUTINE io_write_attchar0(did,name,DATA) 
!*****************************************************************      
      USE hdf5 
                                                                        
      IMPLICIT NONE 
      INTEGER(HID_T),INTENT(IN)  ::did 
      CHARACTER,INTENT(IN)       ::name*(*) 
      CHARACTER,INTENT(IN)         ::DATA*(*) 
      !locals                                                           
      INTEGER(HSIZE_T)::dims(7) 
      INTEGER(HID_t)::atid,sid 
      INTEGER       ::hdferr 
      dims=(/LEN(DATA),0,0,0,0,0,0/) 
      CALL h5screate_simple_f(1,dims,sid,hdferr) 
      CALL h5acreate_f(did, name,H5T_NATIVE_CHARACTER,sid,atid,hdferr) 
      CALL h5awrite_f(atid,H5T_NATIVE_CHARACTER,DATA,dims, hdferr) 
      CALL h5aclose_f(atid,hdferr) 
      CALL h5sclose_f(sid,hdferr) 
      CALL io_check('io_write_attchar0:'//name//'-'//DATA,hdferr) 
      END SUBROUTINE 
!*****************************************************************      
      SUBROUTINE io_read_attchar0(did,name,DATA) 
!*****************************************************************      
      USE hdf5 
                                                                        
      IMPLICIT NONE 
      INTEGER(HID_T),INTENT(IN)  ::did 
      CHARACTER,INTENT(IN)       ::name*(*) 
      CHARACTER,INTENT(OUT)      ::DATA*(*) 
                                                                        
      !locals                                                           
      INTEGER(HSIZE_T)::dims(7) 
      INTEGER(HID_t)::atid 
      INTEGER       ::hdferr 
      dims=(/LEN(DATA),0,0,0,0,0,0/) 
      CALL h5aopen_name_f(did, name, atid, hdferr) 
      CALL h5aread_f(atid,H5T_NATIVE_CHARACTER,DATA,dims, hdferr) 
      CALL h5aclose_f(atid,hdferr) 
      CALL io_check('io_write_attchar0:'//name//'-'//DATA,hdferr) 
      END SUBROUTINE 
                                                                        
                                                                        
!*****************************************************************      
!                                                                       
!     The following subroutines READ or WRITE LOGICAL-attributes        
!                                                                       
!                                                                       
!*****************************************************************      
      SUBROUTINE io_write_attlog0(did,name,DATA) 
!*****************************************************************      
      USE hdf5 
                                                                        
      IMPLICIT NONE 
      INTEGER(HID_T),INTENT(IN)  ::did 
      CHARACTER,INTENT(IN)       ::name*(*) 
      LOGICAL,INTENT(IN)         ::DATA 
      !locals                                                           
      INTEGER(HSIZE_T)::dims(7) 
      INTEGER(HID_t)::atid,sid 
      INTEGER       ::hdferr,dum 
      dims=(/1,0,0,0,0,0,0/) 
      IF (DATA) THEN 
         dum=1 
      ELSE 
         dum=0 
      ENDIF 
      CALL h5screate_simple_f(1,dims,sid,hdferr) 
      CALL h5acreate_f(did, name,H5T_NATIVE_INTEGER,sid,atid,hdferr) 
      CALL h5awrite_f(atid,H5T_NATIVE_INTEGER,dum,dims, hdferr) 
      CALL h5aclose_f(atid,hdferr) 
      CALL h5sclose_f(sid,hdferr) 
      CALL io_check('io_write_attlog0'//name,hdferr) 
      END SUBROUTINE 
!*****************************************************************      
      SUBROUTINE io_read_attlog0(did,name,DATA) 
!*****************************************************************      
      USE hdf5 
                                                                        
      IMPLICIT NONE 
      INTEGER(HID_T),INTENT(IN)  ::did 
      CHARACTER,INTENT(IN)       ::name*(*) 
      LOGICAL,INTENT(OUT)        ::DATA 
                                                                        
      !locals                                                           
      INTEGER(HSIZE_T)::dims(7) 
      INTEGER(HID_t)::atid 
      INTEGER       ::hdferr,dum 
      dims=(/1,0,0,0,0,0,0/) 
      CALL h5aopen_name_f(did, name, atid, hdferr) 
      CALL h5aread_f(atid,H5T_NATIVE_INTEGER,dum,dims, hdferr) 
      DATA=(dum==1) 
      CALL h5aclose_f(atid,hdferr) 
      CALL io_check('io_read_attlog0'//name,hdferr) 
      END SUBROUTINE 
                                                                        
                                                                        
                                                                        
!*****************************************************************      
!                                                                       
!     The following subroutines READ or WRITE REAL-attributes           
!                                                                       
!                                                                       
!*****************************************************************      
      SUBROUTINE io_write_attreal0(did,name,DATA) 
!*****************************************************************      
      USE hdf5 
                                                                        
      IMPLICIT NONE 
      INTEGER(HID_T),INTENT(IN)  ::did 
      CHARACTER,INTENT(IN)       ::name*(*) 
      REAL(rkind),INTENT(IN)    ::DATA 
                                                                        
      CALL io_write_attreal1(did,name,(/DATA/)) 
      END SUBROUTINE 
!*****************************************************************      
      SUBROUTINE io_write_attreal1(did,name,DATA) 
!*****************************************************************      
      USE hdf5 
                                                                        
      IMPLICIT NONE 
      INTEGER(HID_T),INTENT(IN)  ::did 
      CHARACTER,INTENT(IN)       ::name*(*) 
      REAL(rkind),INTENT(IN)    ::DATA(:) 
      !locals                                                           
      INTEGER(HSIZE_T)::dims(7) 
      INTEGER(HID_t)::atid,sid 
      INTEGER       ::hdferr 
      dims=(/SIZE(DATA),0,0,0,0,0,0/) 
      CALL h5screate_simple_f(1,dims,sid,hdferr) 
      CALL h5acreate_f(did, name,H5T_NATIVE_DOUBLE,sid,atid,hdferr) 
      CALL h5awrite_f(atid,H5T_NATIVE_DOUBLE,DATA,dims, hdferr) 
      CALL h5aclose_f(atid,hdferr) 
      CALL h5sclose_f(sid,hdferr) 
      CALL io_check('io_write_attreal1'//name,hdferr) 
      END SUBROUTINE 
!*****************************************************************      
      SUBROUTINE io_write_attreal2(did,name,DATA) 
!*****************************************************************      
      USE hdf5 
                                                                        
      IMPLICIT NONE 
      INTEGER(HID_T),INTENT(IN)  ::did 
      CHARACTER,INTENT(IN)       ::name*(*) 
      REAL(rkind),INTENT(IN)    ::DATA(:,:) 
      !locals                                                           
      INTEGER(HSIZE_T)::dims(7) 
      INTEGER(HID_t)::atid,sid 
      INTEGER       ::hdferr 
      dims=(/SIZE(DATA,1),SIZE(DATA,2),0,0,0,0,0/) 
      CALL h5screate_simple_f(2,dims,sid,hdferr) 
      CALL h5acreate_f(did, name,H5T_NATIVE_DOUBLE,sid,atid,hdferr) 
      CALL h5awrite_f(atid,H5T_NATIVE_DOUBLE,DATA,dims, hdferr) 
      CALL h5aclose_f(atid,hdferr) 
      CALL h5sclose_f(sid,hdferr) 
      CALL io_check('io_write_attreal2'//name,hdferr) 
      END SUBROUTINE 
!*****************************************************************      
      SUBROUTINE io_write_attreal3(did,name,DATA) 
!*****************************************************************      
      USE hdf5 
                                                                        
      IMPLICIT NONE 
      INTEGER(HID_T),INTENT(IN)  ::did 
      CHARACTER,INTENT(IN)       ::name*(*) 
      REAL(rkind),INTENT(IN)    ::DATA(:,:,:) 
      !locals                                                           
      INTEGER(HSIZE_T)::dims(7) 
      INTEGER(HID_t)::atid,sid 
      INTEGER       ::hdferr 
      dims=(/SIZE(DATA,1),SIZE(DATA,2),SIZE(DATA,3),0,0,0,0/) 
      CALL h5screate_simple_f(3,dims,sid,hdferr) 
      CALL h5acreate_f(did, name,H5T_NATIVE_DOUBLE,sid,atid,hdferr) 
      CALL h5awrite_f(atid,H5T_NATIVE_DOUBLE,DATA,dims, hdferr) 
      CALL h5aclose_f(atid,hdferr) 
      CALL h5sclose_f(sid,hdferr) 
      CALL io_check('io_write_attreal3'//name,hdferr) 
      END SUBROUTINE 
!*****************************************************************      
      SUBROUTINE io_read_attreal0(did,name,DATA) 
!*****************************************************************      
      USE hdf5 
                                                                        
      IMPLICIT NONE 
      INTEGER(HID_T),INTENT(IN)  ::did 
      CHARACTER,INTENT(IN)       ::name*(*) 
      REAL(rkind),INTENT(OUT)    ::DATA 
      !locals                                                           
      INTEGER(HSIZE_T)::dims(7) 
      INTEGER(HID_t)::atid 
      INTEGER       ::hdferr 
      REAL          ::buf(1) 
      dims=(/1,0,0,0,0,0,0/) 
      CALL h5aopen_name_f(did, name, atid, hdferr) 
      CALL io_check('io_read_attreal0-1: '//name,hdferr) 
      CALL h5aread_f(atid,H5T_NATIVE_DOUBLE,buf,dims, hdferr) 
      CALL io_check('io_read_attreal0-2: '//name,hdferr) 
      DATA=buf(1) 
      CALL h5aclose_f(atid,hdferr) 
      CALL io_check('io_read_attreal0-3: '//name,hdferr) 
      END SUBROUTINE 
!*****************************************************************      
      SUBROUTINE io_read_attreal1(did,name,DATA) 
!*****************************************************************      
      USE hdf5 
                                                                        
      IMPLICIT NONE 
      INTEGER(HID_T),INTENT(IN)  ::did 
      CHARACTER,INTENT(IN)       ::name*(*) 
      REAL(rkind),INTENT(OUT)    ::DATA(:) 
                                                                        
      !locals                                                           
      INTEGER(HSIZE_T)::dims(7) 
      INTEGER(HID_t)::atid 
      INTEGER       ::hdferr 
      dims=(/SIZE(DATA,1),0,0,0,0,0,0/) 
      CALL h5aopen_name_f(did, name, atid, hdferr) 
      CALL io_check('io_read_attreal1-1: '//name,hdferr) 
      CALL h5aread_f(atid,H5T_NATIVE_DOUBLE,DATA,dims, hdferr) 
      CALL io_check('io_read_attreal1-2: '//name,hdferr) 
      CALL h5aclose_f(atid,hdferr) 
      CALL io_check('io_read_attreal1-3: '//name,hdferr) 
      END SUBROUTINE 
!*****************************************************************      
      SUBROUTINE io_read_attreal2(did,name,DATA) 
!*****************************************************************      
      USE hdf5 
                                                                        
      IMPLICIT NONE 
      INTEGER(HID_T),INTENT(IN)  ::did 
      CHARACTER,INTENT(IN)       ::name*(*) 
      REAL(rkind),INTENT(OUT)    ::DATA(:,:) 
                                                                        
      !locals                                                           
      INTEGER(HSIZE_T)::dims(7) 
      INTEGER(HID_t)::atid 
      INTEGER       ::hdferr 
                                                                        
      dims=(/SIZE(DATA,1),SIZE(DATA,2),0,0,0,0,0/) 
      CALL h5aopen_name_f(did, name, atid, hdferr) 
      CALL io_check('io_real_attreal2-1: '//name,hdferr) 
      CALL h5aread_f(atid,H5T_NATIVE_DOUBLE,DATA,dims,hdferr) 
      CALL io_check('io_real_attreal2-2: '//name,hdferr) 
      CALL h5aclose_f(atid,hdferr) 
      CALL io_check('io_real_attreal2-3: '//name,hdferr) 
      END SUBROUTINE 
!*****************************************************************      
      SUBROUTINE io_read_attreal3(did,name,DATA) 
!*****************************************************************      
      USE hdf5 
                                                                        
      IMPLICIT NONE 
      INTEGER(HID_T),INTENT(IN)  ::did 
      CHARACTER,INTENT(IN)       ::name*(*) 
      REAL(rkind),INTENT(OUT)    ::DATA(:,:,:) 
                                                                        
      !locals                                                           
      INTEGER(HSIZE_T)::dims(7) 
      INTEGER(HID_t)::atid 
      INTEGER       ::hdferr 
      dims=(/SIZE(DATA,1),SIZE(DATA,2),SIZE(DATA,3),0,0,0,0/) 
      CALL h5aopen_name_f(did, name, atid, hdferr) 
      CALL io_check('io_read_attreal3-1: '//name,hdferr) 
      CALL h5aread_f(atid,H5T_NATIVE_DOUBLE,DATA,dims, hdferr) 
      CALL io_check('io_read_attreal3-2: '//name,hdferr) 
      CALL h5aclose_f(atid,hdferr) 
      CALL io_check('io_read_attreal3-3: '//name,hdferr) 
      END SUBROUTINE 
!*****************************************************************      
!                                                                       
!     The following subroutines READ or WRITE INTEGER-attributes        
!                                                                       
!                                                                       
!*****************************************************************      
      SUBROUTINE io_write_attint0(did,name,DATA) 
!*****************************************************************      
      USE hdf5 
                                                                        
      IMPLICIT NONE 
      INTEGER(HID_T),INTENT(IN)  ::did 
      CHARACTER,INTENT(IN)       ::name*(*) 
      INTEGER,INTENT(IN)    ::DATA 
      !locals                                                           
      INTEGER(HSIZE_T)::dims(7) 
      INTEGER(HID_t)::atid,sid 
      INTEGER       ::hdferr 
      dims=(/1,0,0,0,0,0,0/) 
      CALL h5screate_simple_f(1,dims,sid,hdferr) 
      CALL h5acreate_f(did, name,H5T_NATIVE_INTEGER,sid,atid,hdferr) 
      CALL h5awrite_f(atid,H5T_NATIVE_INTEGER,DATA,dims, hdferr) 
      CALL h5aclose_f(atid,hdferr) 
      CALL h5sclose_f(sid,hdferr) 
      CALL io_check('io_write_attint0'//name,hdferr) 
      END SUBROUTINE 
!*****************************************************************      
      SUBROUTINE io_write_attint1(did,name,DATA) 
!*****************************************************************      
      USE hdf5 
                                                                        
      IMPLICIT NONE 
      INTEGER(HID_T),INTENT(IN)  ::did 
      CHARACTER,INTENT(IN)       ::name*(*) 
      INTEGER,INTENT(IN)    ::DATA(:) 
      !locals                                                           
      INTEGER(HSIZE_T)::dims(7) 
      INTEGER(HID_t)::atid,sid 
      INTEGER       ::hdferr 
      dims=(/SIZE(DATA),0,0,0,0,0,0/) 
      CALL h5screate_simple_f(1,dims,sid,hdferr) 
      CALL h5acreate_f(did, name,H5T_NATIVE_INTEGER,sid,atid,hdferr) 
      CALL h5awrite_f(atid,H5T_NATIVE_INTEGER,DATA,dims, hdferr) 
      CALL h5aclose_f(atid,hdferr) 
      CALL h5sclose_f(sid,hdferr) 
      CALL io_check('io_write_attint1'//name,hdferr) 
      END SUBROUTINE 
!*****************************************************************      
      SUBROUTINE io_write_attint2(did,name,DATA) 
!*****************************************************************      
      USE hdf5 
                                                                        
      IMPLICIT NONE 
      INTEGER(HID_T),INTENT(IN)  ::did 
      CHARACTER,INTENT(IN)       ::name*(*) 
      INTEGER,INTENT(IN)    ::DATA(:,:) 
      !locals                                                           
      INTEGER(HSIZE_T)::dims(7) 
      INTEGER(HID_t)::atid,sid 
      INTEGER       ::hdferr 
      dims=(/SIZE(DATA,1),SIZE(DATA,2),0,0,0,0,0/) 
      CALL h5screate_simple_f(2,dims,sid,hdferr) 
      CALL h5acreate_f(did, name,H5T_NATIVE_INTEGER,sid,atid,hdferr) 
      CALL h5awrite_f(atid,H5T_NATIVE_INTEGER,DATA,dims, hdferr) 
      CALL h5aclose_f(atid,hdferr) 
      CALL h5sclose_f(sid,hdferr) 
      CALL io_check('io_write_attint2'//name,hdferr) 
      END SUBROUTINE 
!*****************************************************************      
      SUBROUTINE io_write_attint3(did,name,DATA) 
!*****************************************************************      
      USE hdf5 
                                                                        
      IMPLICIT NONE 
      INTEGER(HID_T),INTENT(IN)  ::did 
      CHARACTER,INTENT(IN)       ::name*(*) 
      INTEGER,INTENT(IN)    ::DATA(:,:,:) 
      !locals                                                           
      INTEGER(HSIZE_T)::dims(7) 
      INTEGER(HID_t)::atid,sid 
      INTEGER       ::hdferr 
      dims=(/SIZE(DATA,1),SIZE(DATA,2),SIZE(DATA,3),0,0,0,0/) 
      CALL h5screate_simple_f(3,dims,sid,hdferr) 
      CALL h5acreate_f(did, name,H5T_NATIVE_INTEGER,sid,atid,hdferr) 
      CALL h5awrite_f(atid,H5T_NATIVE_INTEGER,DATA,dims, hdferr) 
      CALL h5aclose_f(atid,hdferr) 
      CALL h5sclose_f(sid,hdferr) 
      CALL io_check('io_write_attint3'//name,hdferr) 
      END SUBROUTINE 
!*****************************************************************      
      SUBROUTINE io_read_attint0(did,name,DATA) 
!*****************************************************************      
      USE hdf5 
                                                                        
      IMPLICIT NONE 
      INTEGER(HID_T),INTENT(IN)  ::did 
      CHARACTER,INTENT(IN)       ::name*(*) 
      INTEGER,INTENT(OUT)    ::DATA 
      !locals                                                           
      INTEGER(HSIZE_T)::dims(7) 
      INTEGER(HID_t)::atid 
      INTEGER       ::hdferr 
      dims=(/1,0,0,0,0,0,0/) 
      CALL h5aopen_name_f(did, name, atid, hdferr) 
      CALL io_check('io_read_attint0-1: '//name,hdferr) 
      CALL h5aread_f(atid,H5T_NATIVE_INTEGER,DATA,dims, hdferr) 
      CALL io_check('io_read_attint0-2: '//name,hdferr) 
      CALL h5aclose_f(atid,hdferr) 
      CALL io_check('io_read_attint0-3: '//name,hdferr) 
      END SUBROUTINE 
!*****************************************************************      
      SUBROUTINE io_read_attint1(did,name,DATA) 
!*****************************************************************      
      USE hdf5 
                                                                        
      IMPLICIT NONE 
      INTEGER(HID_T),INTENT(IN)  ::did 
      CHARACTER,INTENT(IN)       ::name*(*) 
      INTEGER,INTENT(OUT)    ::DATA(:) 
                                                                        
      !locals                                                           
      INTEGER(HSIZE_T)::dims(7) 
      INTEGER(HID_t)::atid 
      INTEGER       ::hdferr 
      dims=(/SIZE(DATA,1),0,0,0,0,0,0/) 
      CALL h5aopen_name_f(did, name, atid, hdferr) 
      CALL io_check('io_read_attint1-1: '//name,hdferr) 
      CALL h5aread_f(atid,H5T_NATIVE_INTEGER,DATA,dims, hdferr) 
      CALL io_check('io_read_attint1-2: '//name,hdferr) 
      CALL h5aclose_f(atid,hdferr) 
      CALL io_check('io_read_attint1-3: '//name,hdferr) 
      END SUBROUTINE 
!*****************************************************************      
      SUBROUTINE io_read_attint2(did,name,DATA) 
!*****************************************************************      
      USE hdf5 
                                                                        
      IMPLICIT NONE 
      INTEGER(HID_T),INTENT(IN)  ::did 
      CHARACTER,INTENT(IN)       ::name*(*) 
      INTEGER,INTENT(OUT)    ::DATA(:,:) 
                                                                        
      !locals                                                           
      INTEGER(HSIZE_T)::dims(7) 
      INTEGER(HID_t)::atid 
      INTEGER       ::hdferr 
                                                                        
      dims=(/SIZE(DATA,1),SIZE(DATA,2),0,0,0,0,0/) 
      CALL h5aopen_name_f(did, name, atid, hdferr) 
      CALL io_check('io_read_attint2-1: '//name,hdferr) 
      CALL h5aread_f(atid,H5T_NATIVE_INTEGER,DATA,dims,hdferr) 
      CALL io_check('io_read_attint2-2: '//name,hdferr) 
      CALL h5aclose_f(atid,hdferr) 
      CALL io_check('io_read_attint2-3: '//name,hdferr) 
      END SUBROUTINE 
!*****************************************************************      
      SUBROUTINE io_read_attint3(did,name,DATA) 
!*****************************************************************      
      USE hdf5 
                                                                        
      IMPLICIT NONE 
      INTEGER(HID_T),INTENT(IN)  ::did 
      CHARACTER,INTENT(IN)       ::name*(*) 
      INTEGER,INTENT(OUT) ::DATA(:,:,:) 
                                                                        
      !locals                                                           
      INTEGER(HSIZE_T)::dims(7) 
      INTEGER(HID_t)::atid 
      INTEGER       ::hdferr 
      dims=(/SIZE(DATA,1),SIZE(DATA,2),SIZE(DATA,3),0,0,0,0/) 
      CALL h5aopen_name_f(did, name, atid, hdferr) 
      CALL io_check('io_read_attint3-1: '//name,hdferr) 
      CALL h5aread_f(atid,H5T_NATIVE_INTEGER,DATA,dims, hdferr) 
      CALL io_check('io_read_attint3-2: '//name,hdferr) 
      CALL h5aclose_f(atid,hdferr) 
      CALL io_check('io_read_attint3-3: '//name,hdferr) 
      END SUBROUTINE 
                                                                        
      !>                                                                
      END                                           
