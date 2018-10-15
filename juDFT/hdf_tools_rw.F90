!--------------------------------------------------------------------------------
! Copyright (c) 2016 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------

      MODULE m_hdf_tools2 
#include "juDFT_env.h"
!-----------------------------------------------                        
!     major rewrite of hdf_tools                                        
!     this module contains only the                                     
!     io_[read/write](id,start,ncount,data,transprop) to access data
!                                                                       
!     only the generic interface is public                              
!-----------------------------------------------                        
      USE m_hdf_tools4 
      !PRIVATE                                                          
      !<--definitions of interfaces                                     
                                                                        
      INTERFACE io_read 
      MODULE PROCEDURE io_read_real0,io_read_real1,io_read_real2        &
     &     ,io_read_real3,io_read_real4,io_read_real5,io_read_real6     
      MODULE PROCEDURE io_read_integer0,io_read_integer1                &
     &     ,io_read_integer2,io_read_integer3,io_read_integer4          &
     &     ,io_read_integer5,io_read_integer6                           
      MODULE PROCEDURE io_read_complex0,io_read_complex1                &
     &     ,io_read_complex2,io_read_complex3,io_read_complex4          &
     &     ,io_read_complex5                                            
      END INTERFACE 
      INTERFACE io_write 
      MODULE PROCEDURE io_write_real0,io_write_real1,io_write_real2     &
     &     ,io_write_real3,io_write_real4,io_write_real5,io_write_real6 
      MODULE PROCEDURE io_write_integer0,io_write_integer1              &
     &     ,io_write_integer2,io_write_integer3,io_write_integer4       &
     &     ,io_write_integer5,io_write_integer6                         
      MODULE PROCEDURE io_write_complex0,io_write_complex1              &
     &     ,io_write_complex2,io_write_complex3,io_write_complex4       &
     &     ,io_write_complex5                                           
      END INTERFACE 
                                                                        
      !>                                                                
      PUBLIC:: io_read,io_write 
      CONTAINS 
      !<-- implementations of the io_write/read subroutines             
                                                                        
!*****************************************************************      
!                                                                       
!     The following subroutines READ or WRITE int values                
!     from hdf-file                                                     
!                                                                       
!*****************************************************************      
      SUBROUTINE io_read_real0(did,start,ncount,DATA,transprop)
!*****************************************************************      
      USE hdf5 
                                                                        
      IMPLICIT NONE 
      INTEGER(HID_T),INTENT(IN)  ::did 
      INTEGER,INTENT(IN)         ::start(:),                            &
     &     ncount(:)
      REAL(rkind),INTENT(OUT)    ::DATA 
      INTEGER(HID_T),INTENT(IN),OPTIONAL::transprop 
      !locals                                                           
      INTEGER(HSIZE_t)::dims(7),foffset(SIZE(start)),fncount(SIZE(ncount))
      INTEGER(HID_t)::trans,fspace,memspace 
      INTEGER       ::hdferr,s,n 
      IF (.NOT.PRESENT(transprop)) THEN 
         trans=gettransprop() 
      ELSE 
         trans=transprop 
      ENDIF 
                       !                                                
      foffset=start-1 
                       !                                                
      fncount=ncount
                               !                                        
      dims=(/1,0,0,0,0,0,0/) 
      !check if size of ncount is ok!
      s=1 
      trans=H5p_DEFAULT_F 
      DO n=1,SIZE(ncount)
         IF (ncount(n)>0) s=s*ncount(n)
      ENDDO 
      !DO IO                                                            
      CALL h5dget_space_f(did,fspace,hdferr) 
      CALL h5sselect_hyperslab_f(fspace,H5S_SELECT_SET_F,foffset,fncount &
     &     ,hdferr)                                                     
      CALL h5screate_simple_f(1,dims,memspace,hdferr) 
      CALL h5dread_f(did,H5T_NATIVE_DOUBLE,DATA,dims,hdferr,memspace    &
     &     ,fspace,trans)                                               

      CALL io_check("io_read_real0 !",hdferr)
      CALL h5sclose_f(memspace,hdferr) 
      CALL h5sclose_f(fspace,hdferr) 
      CALL cleartransprop(trans) 

      CALL io_check("io_read_real0 !",hdferr)
      END SUBROUTINE 
                                                                        
                                                                        
      SUBROUTINE io_read_real1(did,start,ncount,DATA,transprop)
!*****************************************************************      
      USE hdf5 
                                                                        
      IMPLICIT NONE 
      INTEGER(HID_T),INTENT(IN)  ::did 
      INTEGER,INTENT(IN)         ::start(:),                            &
     &     ncount(:)
      REAL(rkind),INTENT(OUT)    ::DATA(:) 
      INTEGER(HID_T),INTENT(IN),OPTIONAL::transprop 
      !locals                                                           
      INTEGER(HSIZE_t)::dims(7),foffset(SIZE(start)),fncount(SIZE(ncount))
      INTEGER(HID_t)::trans,fspace,memspace 
      INTEGER       ::hdferr,s,n 
      IF (.NOT.PRESENT(transprop)) THEN 
         trans=gettransprop() 
      ELSE 
         trans=transprop 
      ENDIF 
                       !                                                
      foffset=start-1 
                       !                                                
      fncount=ncount
                                          !                             
      dims=(/SIZE(DATA,1),0,0,0,0,0,0/) 
      !check if size of ncount is ok!
      s=1 
      DO n=1,SIZE(ncount)
         IF (ncount(n)>0) s=s*ncount(n)
      ENDDO 
      IF (s/=SIZE(DATA)) CALL hdf_err('mismatch of sizes') 
      !DO IO                                                            
      CALL h5dget_space_f(did,fspace,hdferr) 
      CALL h5sselect_hyperslab_f(fspace,H5S_SELECT_SET_F,foffset,fncount &
     &     ,hdferr)                                                     
      CALL h5screate_simple_f(1,dims,memspace,hdferr) 
      CALL h5dread_f(did,H5T_NATIVE_DOUBLE,DATA,dims,hdferr,memspace    &
     &     ,fspace,trans)                                               

      CALL io_check("io_read_real1 !",hdferr)
      CALL h5sclose_f(memspace,hdferr) 
      CALL h5sclose_f(fspace,hdferr) 
      CALL cleartransprop(trans) 

      CALL io_check("io_read_real1 !",hdferr)
                                                                        
      END SUBROUTINE 
                                                                        
      SUBROUTINE io_read_real2(did,start,ncount,DATA,transprop)
!*****************************************************************      
      USE hdf5 
      IMPLICIT NONE 
                                                                        
      INTEGER(HID_T), INTENT(IN)  :: did 
      INTEGER,        INTENT(IN)  :: start(:),ncount(:)
      REAL,           INTENT(OUT) :: DATA(:,:) 
      INTEGER(HID_T), INTENT(IN), OPTIONAL :: transprop 
! locals                                                                
      INTEGER(HSIZE_t) :: dims(7),foffset(SIZE(start)) 
      INTEGER(HSIZE_t) :: fncount(SIZE(ncount))
      INTEGER(HID_t)   :: trans,fspace,memspace 
      INTEGER          :: hdferr,s,n 
                                                                        
       IF (.NOT.PRESENT(transprop)) THEN 
         trans=gettransprop() 
      ELSE 
         trans=transprop 
      ENDIF 
                                                                        
                       !                                                
      foffset=start-1 
                       !                                                
      fncount=ncount
                                                     !                  
      dims=(/SIZE(DATA,1),SIZE(DATA,2),0,0,0,0,0/) 
      !check if size of ncount is ok!
                               !read nothing                            
      IF (ANY(ncount<1)) RETURN
      s=1 
      DO n=1,SIZE(ncount)
         IF (ncount(n)>0) s=s*ncount(n)
      ENDDO 
      IF (s.NE.SIZE(DATA)) CALL hdf_err('mismatch of sizes') 
      !DO IO                                                            
      CALL h5dget_space_f(did,fspace,hdferr) 
      CALL h5sselect_hyperslab_f(fspace,H5S_SELECT_SET_F,foffset,fncount &
     &     ,hdferr)                                                     
      CALL h5screate_simple_f(2,dims,memspace,hdferr) 
      CALL h5dread_f(did,H5T_NATIVE_DOUBLE,DATA,dims,hdferr,memspace    &
     &     ,fspace,trans)                                               

      CALL io_check("io_read_real2 !",hdferr)
      CALL h5sclose_f(memspace,hdferr) 
      CALL h5sclose_f(fspace,hdferr) 
      CALL cleartransprop(trans) 

      CALL io_check("io_read_real2 !",hdferr)
                                                                        
      END SUBROUTINE 
                                                                        
      SUBROUTINE io_read_real3(did,start,ncount,DATA,transprop)
!*****************************************************************      
      USE hdf5 
      IMPLICIT NONE 
      INTEGER(HID_T),INTENT(IN)  ::did 
      INTEGER,INTENT(IN)         ::start(:),                            &
     &     ncount(:)
      REAL(rkind),INTENT(OUT)    ::DATA(:,:,:) 
      INTEGER(HID_T),INTENT(IN),OPTIONAL::transprop 
      !locals                                                           
      INTEGER(HSIZE_t)::dims(7),foffset(SIZE(start)),fncount(SIZE(ncount))
      INTEGER(HID_t)::trans,fspace,memspace 
      INTEGER       ::hdferr,s,n 
                                                                        
      IF (.NOT.PRESENT(transprop)) THEN 
         trans=gettransprop() 
      ELSE 
         trans=transprop 
      ENDIF 
                       !                                                
      foffset=start-1 
                       !                                                
      fncount=ncount
                                                                !       
      dims=(/SIZE(DATA,1),SIZE(DATA,2),SIZE(DATA,3),0,0,0,0/) 
      !check if size of ncount is ok!
      s=1 
                               !read nothing                            
      IF (ANY(ncount<1)) RETURN
      DO n=1,SIZE(ncount)
         IF (ncount(n)>0) s=s*ncount(n)
      ENDDO 
      IF (s.NE.SIZE(DATA)) CALL hdf_err("mismatch of sizes") 
      !DO IO                                                            
      CALL h5dget_space_f(did,fspace,hdferr) 
      CALL h5sselect_hyperslab_f(fspace,H5S_SELECT_SET_F,foffset,fncount &
     &     ,hdferr)                                                     
      CALL h5screate_simple_f(3,dims,memspace,hdferr) 
      CALL h5dread_f(did,H5T_NATIVE_DOUBLE,DATA,dims,hdferr,memspace    &
     &     ,fspace,trans)                                               

      CALL io_check("io_read_real3 !",hdferr)
      CALL h5sclose_f(memspace,hdferr) 
      CALL h5sclose_f(fspace,hdferr) 
      CALL cleartransprop(trans) 

      CALL io_check("io_read_real3 !",hdferr)
                                                                        
      END SUBROUTINE 
                                                                        
      SUBROUTINE io_read_real4(did,start,ncount,DATA,transprop)
!*****************************************************************      
      USE hdf5 
      IMPLICIT NONE 
                                                                        
      INTEGER(HID_T),INTENT(IN)  :: did 
      INTEGER,       INTENT(IN)  :: start(:),ncount(:)
      REAL,          INTENT(OUT) :: DATA(:,:,:,:) 
      INTEGER(HID_T),INTENT(IN),OPTIONAL :: transprop 
! locals                                                                
      INTEGER(HSIZE_t)::dims(7),foffset(SIZE(start)),fncount(SIZE(ncount))
      INTEGER(HID_t)::trans,fspace,memspace 
      INTEGER       ::hdferr,s,n 
                                                                        
      IF (.NOT.PRESENT(transprop)) THEN 
         trans=gettransprop() 
      ELSE 
         trans=transprop 
      ENDIF 
                       !                                                
      foffset=start-1 
                       !                                                
      fncount=ncount
      dims=(/SIZE(DATA,1),SIZE(DATA,2),SIZE(DATA,3),SIZE(DATA,4),0,0,0/) 
      !check if size of ncount is ok!
                               !read nothing                            
      IF (ANY(ncount<1)) RETURN
      s=1 
      DO n=1,SIZE(ncount)
         IF (ncount(n)>0) s=s*ncount(n)
      ENDDO 
      IF (s.NE.SIZE(DATA)) CALL hdf_err("mismatch of sizes") 
      !DO IO                                                            
      CALL h5dget_space_f(did,fspace,hdferr) 
      CALL h5sselect_hyperslab_f(fspace,H5S_SELECT_SET_F,foffset,fncount &
     &     ,hdferr)                                                     
      CALL h5screate_simple_f(4,dims,memspace,hdferr) 
      CALL h5dread_f(did,H5T_NATIVE_DOUBLE,DATA,dims,hdferr,memspace    &
     &     ,fspace,trans)                                               

      CALL io_check("io_read_real4 !",hdferr)
      CALL h5sclose_f(memspace,hdferr) 
      CALL h5sclose_f(fspace,hdferr) 
      CALL cleartransprop(trans) 

      CALL io_check("io_read_real4 !",hdferr)
      END SUBROUTINE 
                                                                        
      SUBROUTINE io_read_real5(did,start,ncount,DATA,transprop)
!*****************************************************************      
      USE hdf5 
                                                                        
      IMPLICIT NONE 
      INTEGER(HID_T),INTENT(IN)  ::did 
      INTEGER,INTENT(IN)         ::start(:),                            &
     &     ncount(:)
      REAL(rkind),INTENT(OUT)    ::DATA(:,:,:,:,:) 
      INTEGER(HID_T),INTENT(IN),OPTIONAL::transprop 
      !locals                                                           
      INTEGER(HSIZE_t)::dims(7),foffset(SIZE(start)),fncount(SIZE(ncount))
      INTEGER(HID_t)::trans,fspace,memspace 
      INTEGER       ::hdferr,s,n 
       IF (.NOT.PRESENT(transprop)) THEN 
         trans=gettransprop() 
      ELSE 
         trans=transprop 
      ENDIF 
                       !                                                
      foffset=start-1 
                       !                                                
      fncount=ncount
                                                                 !      
      dims=(/SIZE(DATA,1),SIZE(DATA,2),SIZE(DATA,3),SIZE(DATA,4)        &
     &     ,SIZE(DATA,5),0,0/)                                          
      !check if size of ncount is ok!
                               !read nothing                            
      IF (ANY(ncount<1)) RETURN
      s=1 
      DO n=1,SIZE(ncount)
         IF (ncount(n)>0) s=s*ncount(n)
      ENDDO 
      IF (s.NE.SIZE(DATA)) CALL hdf_err("mismatch of sizes") 
      !DO IO                                                            
      CALL h5dget_space_f(did,fspace,hdferr) 
      CALL h5sselect_hyperslab_f(fspace,H5S_SELECT_SET_F,foffset,fncount &
     &     ,hdferr)                                                     
      CALL h5screate_simple_f(5,dims,memspace,hdferr) 
      CALL h5dread_f(did,H5T_NATIVE_DOUBLE,DATA,dims,hdferr,memspace    &
     &     ,fspace,trans)                                               

      CALL io_check("io_read_real5 !",hdferr)
      CALL h5sclose_f(memspace,hdferr) 
      CALL h5sclose_f(fspace,hdferr) 
      CALL cleartransprop(trans) 

      CALL io_check("io_read_real5 !",hdferr)
                                                                        
      END SUBROUTINE 
                                                                        
      SUBROUTINE io_read_real6(did,start,ncount,DATA,transprop)
!*****************************************************************      
      USE hdf5 
                                                                        
      IMPLICIT NONE 
      INTEGER(HID_T),INTENT(IN)  ::did 
      INTEGER,INTENT(IN)         ::start(:),                            &
     &     ncount(:)
      REAL(rkind),INTENT(OUT)    ::DATA(:,:,:,:,:,:) 
      INTEGER(HID_T),INTENT(IN),OPTIONAL::transprop 
      !locals                                                           
      INTEGER(HSIZE_t)::dims(7),foffset(SIZE(start)),fncount(SIZE(ncount))
      INTEGER(HID_t)::trans,fspace,memspace 
      INTEGER       ::hdferr,s,n 
      IF (.NOT.PRESENT(transprop)) THEN 
         trans=gettransprop() 
      ELSE 
         trans=transprop 
      ENDIF 
                       !                                                
      foffset=start-1 
                       !                                                
      fncount=ncount
                                                                    !   
      dims(:)=(/SIZE(DATA,1),SIZE(DATA,2),SIZE(DATA,3),SIZE(DATA,4)     &
     &     ,SIZE(DATA,5),SIZE(DATA,6),0/)                               
      !check if size of ncount is ok!
                               !read nothing                            
      IF (ANY(ncount<1)) RETURN
      s=1 
      DO n=1,SIZE(ncount)
         IF (ncount(n)>0) s=s*ncount(n)
      ENDDO 
      IF (s.NE.SIZE(DATA)) CALL hdf_err("mismatch of sizes") 
      !DO IO                                                            
      CALL h5dget_space_f(did,fspace,hdferr) 
      CALL h5sselect_hyperslab_f(fspace,H5S_SELECT_SET_F,foffset,fncount &
     &     ,hdferr)                                                     
      CALL h5screate_simple_f(6,dims,memspace,hdferr) 
      CALL h5dread_f(did,H5T_NATIVE_DOUBLE,DATA,dims,hdferr,memspace    &
     &     ,fspace,trans)                                               

      CALL io_check("io_read_real6 !",hdferr)
      CALL h5sclose_f(memspace,hdferr) 
      CALL h5sclose_f(fspace,hdferr) 
      CALL cleartransprop(trans) 

      CALL io_check("io_read_real6 !",hdferr)
                                                                        
      END SUBROUTINE 
                                                                        
!---------------------------------------------------------------------- 
                                                                        
      SUBROUTINE io_write_real0(did,start,ncount,DATA,transprop)
!*****************************************************************      
      USE hdf5 
                                                                        
      IMPLICIT NONE 
                                                                        
!arguments                                                              
      INTEGER(HID_T), INTENT(IN)  :: did 
      INTEGER,        INTENT(IN)  :: start(:),ncount(:)
      REAL(rkind),    INTENT(IN)  :: DATA 
      INTEGER(HID_T), INTENT(IN), OPTIONAL :: transprop 
!locals                                                                 
      INTEGER(HSIZE_t) :: dims(7) 
      INTEGER(HSIZE_t) :: foffset(SIZE(start)),fncount(SIZE(ncount))
      INTEGER(HID_t)   :: trans,fspace,memspace 
      INTEGER          :: hdferr,s,n 
                                                                        
                                ! write a single real                   
      dims = (/1,0,0,0,0,0,0/) 
                                                                        
      IF (.NOT.PRESENT(transprop)) THEN 
         trans = gettransprop() 
      ELSE 
         trans = transprop 
      ENDIF 
      foffset = start-1 
      fncount  = ncount
                               ! write nothing                          
      IF (ANY(ncount<1)) RETURN
                               ! check if size of ncount is ok
      s=1 
      DO n = 1, SIZE(ncount)
         IF (ncount(n)>0) s=s*ncount(n)
      ENDDO 
      IF (s.NE.1) CALL hdf_err("mismatch of sizes") 
!do I/O                                                                 
      CALL h5dget_space_f(                                              &
     &                    did,fspace,hdferr)                            
      CALL h5sselect_hyperslab_f(                                       &
     &                           fspace,H5S_SELECT_SET_F,               &
     &                           foffset,fncount,hdferr)
      CALL h5screate_simple_f(                                          &
     &                        1,dims,memspace,hdferr)                   
      CALL h5dwrite_f(                                                  &
     &                did,H5T_NATIVE_DOUBLE,DATA,dims,hdferr,           &
     &                memspace,fspace,trans)                            
      CALL h5sclose_f(                                                  &
     &                memspace,hdferr)                                  
      CALL h5sclose_f(                                                  &
     &                fspace,hdferr)                                    
      CALL cleartransprop(trans) 

      CALL io_check("io_write_real0 !",hdferr)
      END SUBROUTINE io_write_real0 
                                                                        
!---------------------------------------------------------------------- 
                                                                        
      SUBROUTINE io_write_real1(did,start,ncount,DATA,transprop)
!*****************************************************************      
      USE hdf5 
                                                                        
      IMPLICIT NONE 
                                                                        
!arguments                                                              
      INTEGER(HID_T), INTENT(IN)  :: did 
      INTEGER,        INTENT(IN)  :: start(:),ncount(:)
      REAL(rkind),    INTENT(IN)  :: DATA(:) 
      INTEGER(HID_T), INTENT(IN), OPTIONAL :: transprop 
!locals                                                                 
      INTEGER(HSIZE_t) :: dims(7) 
      INTEGER(HSIZE_t) :: foffset(SIZE(start)),fncount(SIZE(ncount))
      INTEGER(HID_t)   :: trans,fspace,memspace 
      INTEGER          :: hdferr,s,n 
                                                                        
                                            ! write 1-dim array         
      dims = (/SIZE(DATA,1),0,0,0,0,0,0/) 
                                                                        
      IF (.NOT.PRESENT(transprop)) THEN 
         trans = gettransprop() 
      ELSE 
         trans = transprop 
      ENDIF 
      foffset = start-1 
      fncount  = ncount
                               ! write nothing                          
      IF (ANY(ncount<1)) RETURN
                               ! check if size of ncount is ok
      s=1 
      DO n = 1, SIZE(ncount)
         IF (ncount(n)>0) s=s*ncount(n)
      ENDDO 
      IF (s.NE.SIZE(DATA)) CALL hdf_err("Mismatch of sizes")
                                                                        
!do I/O                                                                 
      CALL h5dget_space_f(                                              &
     &                    did,fspace,hdferr)                            
      CALL h5sselect_hyperslab_f(                                       &
     &                           fspace,H5S_SELECT_SET_F,               &
     &                           foffset,fncount,hdferr)
      CALL h5screate_simple_f(                                          &
     &                        1,dims,memspace,hdferr)                   
      CALL h5dwrite_f(                                                  &
     &                did,H5T_NATIVE_DOUBLE,DATA,dims,hdferr,           &
     &                memspace,fspace,trans)                            
      CALL h5sclose_f(                                                  &
     &                memspace,hdferr)                                  
      CALL h5sclose_f(                                                  &
     &                fspace,hdferr)                                    
      CALL cleartransprop(trans) 

      CALL io_check("io_write_real1 !",hdferr)
                                                                        
      END SUBROUTINE io_write_real1 
                                                                        
!---------------------------------------------------------------------- 
                                                                        
      SUBROUTINE io_write_real2(did,start,ncount,DATA,transprop)
!*****************************************************************      
      USE hdf5 
                                                                        
      IMPLICIT NONE 
      INTEGER(HID_T),INTENT(IN)  ::did 
      INTEGER,INTENT(IN)         ::start(:),                            &
     &     ncount(:)
      REAL(rkind),INTENT(IN)    ::DATA(:,:) 
      INTEGER(HID_T),INTENT(IN),OPTIONAL::transprop 
      !locals                                                           
      INTEGER(HSIZE_t)::dims(7),foffset(SIZE(start)),fncount(SIZE(ncount))
      INTEGER(HID_t)::trans,fspace,memspace 
      INTEGER       ::hdferr,s,n 
      IF (.NOT.PRESENT(transprop)) THEN 
         trans=gettransprop() 
      ELSE 
         trans=transprop 
      ENDIF 
                       !                                                
      foffset=start-1 
                       !                                                
      fncount=ncount
                                                     !                  
      dims=(/SIZE(DATA,1),SIZE(DATA,2),0,0,0,0,0/) 
      !check if size of ncount is ok!
                               !read nothing                            
      IF (ANY(ncount<1)) RETURN
      s=1 
      DO n=1,SIZE(ncount)
         IF (ncount(n)>0) s=s*ncount(n)
      ENDDO 
      IF (s.NE.SIZE(DATA)) CALL hdf_err("mismatch of sizes") 
      !DO IO                                                            
      CALL h5dget_space_f(did,fspace,hdferr) 
      CALL h5sselect_hyperslab_f(fspace,H5S_SELECT_SET_F,foffset,fncount &
     &     ,hdferr)                                                     
      CALL h5screate_simple_f(2,dims,memspace,hdferr) 
      CALL h5dwrite_f(did,H5T_NATIVE_DOUBLE,DATA,dims,hdferr,memspace   &
     &     ,fspace,trans)                                               
      CALL h5sclose_f(memspace,hdferr) 
      CALL h5sclose_f(fspace,hdferr) 
      CALL cleartransprop(trans) 

      CALL io_check("io_write_real2 !",hdferr)
                                                                        
      END SUBROUTINE 
                                                                        
      SUBROUTINE io_write_real3(did,start,ncount,DATA,transprop)
!*****************************************************************      
      USE hdf5 
      IMPLICIT NONE 
                                                                        
      INTEGER(HID_T),INTENT(IN) :: did 
      INTEGER,       INTENT(IN) :: start(:),ncount(:)
      REAL,          INTENT(IN) :: DATA(:,:,:) 
      INTEGER(HID_T),INTENT(IN),OPTIONAL :: transprop 
!  locals                                                               
      INTEGER(HSIZE_t)::dims(7),foffset(SIZE(start)),fncount(SIZE(ncount))
      INTEGER(HID_t)::trans,fspace,memspace 
      INTEGER       ::hdferr,s,n 
      IF (.NOT.PRESENT(transprop)) THEN 
         trans=gettransprop() 
      ELSE 
         trans=transprop 
      ENDIF 
                       !                                                
      foffset=start-1 
                       !                                                
      fncount=ncount
                                                                !       
      dims=(/SIZE(DATA,1),SIZE(DATA,2),SIZE(DATA,3),0,0,0,0/) 
      !check if size of ncount is ok!
                               !read nothing                            
      IF (ANY(ncount<1)) RETURN
      s=1 
      DO n=1,SIZE(ncount)
         IF (ncount(n)>0) s=s*ncount(n)
      ENDDO 
      IF (s.NE.SIZE(DATA)) CALL hdf_err("mismatch of sizes") 
      !DO IO                                                            
      CALL h5dget_space_f(did,fspace,hdferr) 
      CALL h5sselect_hyperslab_f(fspace,H5S_SELECT_SET_F,foffset,fncount &
     &     ,hdferr)                                                     
      CALL h5screate_simple_f(3,dims,memspace,hdferr) 
      CALL h5dwrite_f(did,H5T_NATIVE_DOUBLE,DATA,dims,hdferr,memspace   &
     &     ,fspace,trans)                                               
      CALL h5sclose_f(memspace,hdferr) 
      CALL h5sclose_f(fspace,hdferr) 
      CALL cleartransprop(trans) 

      CALL io_check("io_write_real3 !",hdferr)
      END SUBROUTINE 
                                                                        
      SUBROUTINE io_write_real4(did,start,ncount,DATA,transprop)
!*****************************************************************      
      USE hdf5 
                                                                        
      IMPLICIT NONE 
      INTEGER(HID_T),INTENT(IN)  ::did 
      INTEGER,INTENT(IN)         ::start(:),                            &
     &     ncount(:)
      REAL(rkind),INTENT(IN)    ::DATA(:,:,:,:) 
      INTEGER(HID_T),INTENT(IN),OPTIONAL::transprop 
      !locals                                                           
      INTEGER(HSIZE_t)::dims(7),foffset(SIZE(start)),fncount(SIZE(ncount))
      INTEGER(HID_t)::trans,fspace,memspace 
      INTEGER       ::hdferr,s,n 
      IF (.NOT.PRESENT(transprop)) THEN 
         trans=gettransprop() 
      ELSE 
         trans=transprop 
      ENDIF 
                       !                                                
      foffset=start-1 
                       !                                                
      fncount=ncount
      dims=(/SIZE(DATA,1),SIZE(DATA,2),SIZE(DATA,3),SIZE(DATA,4),0,0,0/) 
      !check if size of ncount is ok!
                               !read nothing                            
      IF (ANY(ncount<1)) RETURN
      s=1 
      DO n=1,SIZE(ncount)
         IF (ncount(n)>0) s=s*ncount(n)
      ENDDO
      IF (s.NE.SIZE(DATA)) CALL hdf_err("mismatch of sizes") 
      !DO IO                                                            
      CALL h5dget_space_f(did,fspace,hdferr) 
      CALL h5sselect_hyperslab_f(fspace,H5S_SELECT_SET_F,foffset,fncount &
     &     ,hdferr)                                                     
      CALL h5screate_simple_f(4,dims,memspace,hdferr) 
      CALL h5dwrite_f(did,H5T_NATIVE_DOUBLE,DATA,dims,hdferr,memspace   &
     &     ,fspace,trans)                                               
      CALL h5sclose_f(memspace,hdferr) 
      CALL h5sclose_f(fspace,hdferr) 
      CALL cleartransprop(trans) 

      CALL io_check("io_write_real4 !",hdferr)
                                                                        
      END SUBROUTINE 
                                                                        
      SUBROUTINE io_write_real5(did,start,ncount,DATA,transprop)
!*****************************************************************      
      USE hdf5 
                                                                        
      IMPLICIT NONE 
      INTEGER(HID_T),INTENT(IN)  ::did 
      INTEGER,INTENT(IN)         ::start(:),                            &
     &     ncount(:)
      REAL(rkind),INTENT(IN)    ::DATA(:,:,:,:,:) 
      INTEGER(HID_T),INTENT(IN),OPTIONAL::transprop 
      !locals                                                           
      INTEGER(HSIZE_t)::dims(7),foffset(SIZE(start)),fncount(SIZE(ncount))
      INTEGER(HID_t)::trans,fspace,memspace 
      INTEGER       ::hdferr,s,n 
      IF (.NOT.PRESENT(transprop)) THEN 
         trans=gettransprop() 
      ELSE 
         trans=transprop 
      ENDIF 
                       !                                                
      foffset=start-1 
                       !                                                
      fncount=ncount
                                                                 !      
      dims=(/SIZE(DATA,1),SIZE(DATA,2),SIZE(DATA,3),SIZE(DATA,4)        &
     &     ,SIZE(DATA,5),0,0/)                                          
      !check if size of ncount is ok!
                               !read nothing                            
      IF (ANY(ncount<1)) RETURN
      s=1 
      DO n=1,SIZE(ncount)
         IF (ncount(n)>0) s=s*ncount(n)
      ENDDO 
      IF (s.NE.SIZE(DATA)) CALL hdf_err("mismatch of sizes") 
      !DO IO                                                            
      CALL h5dget_space_f(did,fspace,hdferr) 
      CALL h5sselect_hyperslab_f(fspace,H5S_SELECT_SET_F,foffset,fncount &
     &     ,hdferr)                                                     
      CALL h5screate_simple_f(5,dims,memspace,hdferr) 
      CALL h5dwrite_f(did,H5T_NATIVE_DOUBLE,DATA,dims,hdferr,memspace   &
     &     ,fspace,trans)                                               
      CALL h5sclose_f(memspace,hdferr) 
      CALL h5sclose_f(fspace,hdferr) 
      CALL cleartransprop(trans) 

      CALL io_check("io_write_real5 !",hdferr)
      END SUBROUTINE 
                                                                        
      SUBROUTINE io_write_real6(did,start,ncount,DATA,transprop)
!*****************************************************************      
      USE hdf5 
                                                                        
      IMPLICIT NONE 
      INTEGER(HID_T),INTENT(IN)  ::did 
      INTEGER,INTENT(IN)         ::start(:),                            &
     &     ncount(:)
      REAL(rkind),INTENT(IN)    ::DATA(:,:,:,:,:,:) 
      INTEGER(HID_T),INTENT(IN),OPTIONAL::transprop 
      !locals                                                           
      INTEGER(HSIZE_t)::dims(7),foffset(SIZE(start)),fncount(SIZE(ncount))
      INTEGER(HID_t)::trans,fspace,memspace 
      INTEGER       ::hdferr,s,n 
      IF (.NOT.PRESENT(transprop)) THEN 
         trans=gettransprop() 
      ELSE 
         trans=transprop 
      ENDIF 
                       !                                                
      foffset=start-1 
                       !                                                
      fncount=ncount
                                                                    !   
      dims(:)=(/SIZE(DATA,1),SIZE(DATA,2),SIZE(DATA,3),SIZE(DATA,4)     &
     &     ,SIZE(DATA,5),SIZE(DATA,6),0/)                               
      !check if size of ncount is ok!
                               !read nothing                            
      IF (ANY(ncount<1)) RETURN
      s=1 
      DO n=1,SIZE(ncount)
         IF (ncount(n)>0) s=s*ncount(n)
      ENDDO 
      IF (s.NE.SIZE(DATA)) CALL hdf_err("mismatch of sizes") 
      !DO IO                                                            
      CALL h5dget_space_f(did,fspace,hdferr) 
      CALL h5sselect_hyperslab_f(fspace,H5S_SELECT_SET_F,foffset,fncount &
     &     ,hdferr)                                                     
      CALL h5screate_simple_f(6,dims,memspace,hdferr) 
      CALL h5dwrite_f(did,H5T_NATIVE_DOUBLE,DATA,dims,hdferr,memspace   &
     &     ,fspace,trans)                                               
      CALL h5sclose_f(memspace,hdferr) 
      CALL h5sclose_f(fspace,hdferr) 
      CALL cleartransprop(trans) 

      CALL io_check("io_write_real6 !",hdferr)
                                                                        
      END SUBROUTINE 
                                                                        
!*****************************************************************      
!                                                                       
!     The following subroutines READ or WRITE INTEGER values            
!     from hdf-file                                                     
!                                                                       
!*****************************************************************      
      SUBROUTINE io_read_integer0(did,start,ncount,DATA,transprop)
!*****************************************************************      
      USE hdf5 
                                                                        
      IMPLICIT NONE 
      INTEGER(HID_T),INTENT(IN)  ::did 
      INTEGER,INTENT(IN)         ::start(:),                            &
     &     ncount(:)
      INTEGER,INTENT(OUT)    ::DATA 
                                                                        
      INTEGER(HID_T),INTENT(IN),OPTIONAL::transprop 
      !locals                                                           
      INTEGER(HSIZE_t)::dims(7),foffset(SIZE(start)),fncount(SIZE(ncount))
      INTEGER(HID_t)::trans,fspace,memspace 
      INTEGER       ::hdferr,s,n 
      IF (.NOT.PRESENT(transprop)) THEN 
         trans=gettransprop() 
      ELSE 
         trans=transprop 
      ENDIF 
                       !                                                
      foffset=start-1 
                       !                                                
      fncount=ncount
                               !                                        
      dims=(/1,0,0,0,0,0,0/) 
      !check if size of ncount is ok!
                               !read nothing                            
      IF (ANY(ncount<1)) RETURN
      s=1 
      DO n=1,SIZE(ncount)
         IF (ncount(n)>0) s=s*ncount(n)
      ENDDO 
      !DO IO                                                            
      CALL h5dget_space_f(did,fspace,hdferr) 
      CALL h5sselect_hyperslab_f(fspace,H5S_SELECT_SET_F,foffset,fncount &
     &     ,hdferr)                                                     
      CALL h5screate_simple_f(1,dims,memspace,hdferr) 
      CALL h5dread_f(did,H5T_NATIVE_INTEGER,DATA,dims,hdferr,memspace   &
     &     ,fspace,trans)                                               

      CALL io_check("io_read_integer0 !",hdferr)
      CALL h5sclose_f(memspace,hdferr) 
      CALL h5sclose_f(fspace,hdferr) 
      CALL cleartransprop(trans) 

      CALL io_check("io_read_integer0 !",hdferr)
      END SUBROUTINE 
                                                                        
      SUBROUTINE io_read_integer1(did,start,ncount,DATA,transprop)
!*****************************************************************      
      USE hdf5 
                                                                        
      IMPLICIT NONE 
      INTEGER(HID_T),INTENT(IN)  ::did 
      INTEGER,INTENT(IN)         ::start(:),                            &
     &     ncount(:)
      INTEGER,INTENT(OUT)    ::DATA(:) 
      INTEGER(HID_T),INTENT(IN),OPTIONAL::transprop 
      !locals                                                           
      INTEGER(HSIZE_t)::dims(7),foffset(SIZE(start)),fncount(SIZE(ncount))
      INTEGER(HID_t)::trans,fspace,memspace 
      INTEGER       ::hdferr,s,n 
      IF (.NOT.PRESENT(transprop)) THEN 
         trans=gettransprop() 
      ELSE 
         trans=transprop 
      ENDIF 
                       !                                                
      foffset=start-1 
                       !                                                
      fncount=ncount
                                          !                             
      dims=(/SIZE(DATA,1),0,0,0,0,0,0/) 
      !check if size of ncount is ok!
                               !read nothing                            
      IF (ANY(ncount<1)) RETURN
      s=1 
      DO n=1,SIZE(ncount)
         IF (ncount(n)>0) s=s*ncount(n)
      ENDDO 
      IF (s.NE.SIZE(DATA)) CALL hdf_err("mismatch of sizes") 
      !DO IO                                                            
      CALL h5dget_space_f(did,fspace,hdferr) 
      CALL h5sselect_hyperslab_f(fspace,H5S_SELECT_SET_F,foffset,fncount &
     &     ,hdferr)                                                     
      CALL h5screate_simple_f(1,dims,memspace,hdferr) 
      CALL h5dread_f(did,H5T_NATIVE_INTEGER,DATA,dims,hdferr,memspace   &
     &     ,fspace,trans)                                               

      CALL io_check("io_read_integer1 !",hdferr)
      CALL h5sclose_f(memspace,hdferr) 
      CALL h5sclose_f(fspace,hdferr) 
      CALL cleartransprop(trans) 

      CALL io_check("io_read_integer1 !",hdferr)
      END SUBROUTINE 
                                                                        
      SUBROUTINE io_read_integer2(did,start,ncount,DATA,transprop)
!*****************************************************************      
      USE hdf5 
                                                                        
      IMPLICIT NONE 
      INTEGER(HID_T),INTENT(IN)  ::did 
      INTEGER,INTENT(IN)         ::start(:),                            &
     &     ncount(:)
      INTEGER,INTENT(OUT)    ::DATA(:,:) 
      INTEGER(HID_T),INTENT(IN),OPTIONAL::transprop 
      !locals                                                           
      INTEGER(HSIZE_t)::dims(7),foffset(SIZE(start)),fncount(SIZE(ncount))
      INTEGER(HID_t)::trans,fspace,memspace 
      INTEGER       ::hdferr,s,n 
      IF (.NOT.PRESENT(transprop)) THEN 
         trans=gettransprop() 
      ELSE 
         trans=transprop 
      ENDIF 
                       !                                                
      foffset=start-1 
                       !                                                
      fncount=ncount
                                                     !                  
      dims=(/SIZE(DATA,1),SIZE(DATA,2),0,0,0,0,0/) 
      !check if size of ncount is ok!
                               !read nothing                            
      IF (ANY(ncount<1)) RETURN
      s=1 
      DO n=1,SIZE(ncount)
         IF (ncount(n)>0) s=s*ncount(n)
      ENDDO 
      IF (s.NE.SIZE(DATA)) CALL hdf_err("mismatch of sizes") 
      !DO IO                                                            
      CALL h5dget_space_f(did,fspace,hdferr) 
      CALL h5sselect_hyperslab_f(fspace,H5S_SELECT_SET_F,foffset,fncount &
     &     ,hdferr)                                                     
      CALL h5screate_simple_f(2,dims,memspace,hdferr) 
      CALL h5dread_f(did,H5T_NATIVE_INTEGER,DATA,dims,hdferr,memspace   &
     &     ,fspace,trans)                                               

      CALL io_check("io_read_integer2 !",hdferr)
      CALL h5sclose_f(memspace,hdferr) 
      CALL h5sclose_f(fspace,hdferr) 
      CALL cleartransprop(trans) 

      CALL io_check("io_read_integer2 !",hdferr)
      END SUBROUTINE 
                                                                        
      SUBROUTINE io_read_integer3(did,start,ncount,DATA,transprop)
!*****************************************************************      
      USE hdf5 
                                                                        
      IMPLICIT NONE 
      INTEGER(HID_T),INTENT(IN)  ::did 
      INTEGER,INTENT(IN)         ::start(:),                            &
     &     ncount(:)
      INTEGER,INTENT(OUT)    ::DATA(:,:,:) 
      INTEGER(HID_T),INTENT(IN),OPTIONAL::transprop 
      !locals                                                           
      INTEGER(HSIZE_t)::dims(7),foffset(SIZE(start)),fncount(SIZE(ncount))
      INTEGER(HID_t)::trans,fspace,memspace 
      INTEGER       ::hdferr,s,n 
      IF (.NOT.PRESENT(transprop)) THEN 
         trans=gettransprop() 
      ELSE 
         trans=transprop 
      ENDIF 
                       !                                                
      foffset=start-1 
                       !                                                
      fncount=ncount
                                                                !       
      dims=(/SIZE(DATA,1),SIZE(DATA,2),SIZE(DATA,3),0,0,0,0/) 
      !check if size of ncount is ok!
                               !read nothing                            
      IF (ANY(ncount<1)) RETURN
      s=1 
      DO n=1,SIZE(ncount)
         IF (ncount(n)>0) s=s*ncount(n)
      ENDDO 
      IF (s.NE.SIZE(DATA)) CALL hdf_err("mismatch of sizes") 
      !DO IO                                                            
      CALL h5dget_space_f(did,fspace,hdferr) 
      CALL h5sselect_hyperslab_f(fspace,H5S_SELECT_SET_F,foffset,fncount &
     &     ,hdferr)                                                     
      CALL h5screate_simple_f(3,dims,memspace,hdferr) 
      CALL h5dread_f(did,H5T_NATIVE_INTEGER,DATA,dims,hdferr,memspace   &
     &     ,fspace,trans)                                               

      CALL io_check("io_read_integer3 !",hdferr)
      CALL h5sclose_f(memspace,hdferr) 
      CALL h5sclose_f(fspace,hdferr) 
      CALL cleartransprop(trans) 

      CALL io_check("io_read_integer3 !",hdferr)
      END SUBROUTINE 
                                                                        
      SUBROUTINE io_read_integer4(did,start,ncount,DATA,transprop)
!*****************************************************************      
      USE hdf5 
                                                                        
      IMPLICIT NONE 
      INTEGER(HID_T),INTENT(IN)  ::did 
      INTEGER,INTENT(IN)         ::start(:),                            &
     &     ncount(:)
      INTEGER,INTENT(OUT)    ::DATA(:,:,:,:) 
      INTEGER(HID_T),INTENT(IN),OPTIONAL::transprop 
      !locals                                                           
      INTEGER(HSIZE_t)::dims(7),foffset(SIZE(start)),fncount(SIZE(ncount))
      INTEGER(HID_t)::trans,fspace,memspace 
      INTEGER       ::hdferr,s,n 
      IF (.NOT.PRESENT(transprop)) THEN 
         trans=gettransprop() 
      ELSE 
         trans=transprop 
      ENDIF 
                       !                                                
      foffset=start-1 
                       !                                                
      fncount=ncount
      dims=(/SIZE(DATA,1),SIZE(DATA,2),SIZE(DATA,3),SIZE(DATA,4),0,0,0/) 
      !check if size of ncount is ok!
                               !read nothing                            
      IF (ANY(ncount<1)) RETURN
      s=1 
      DO n=1,SIZE(ncount)
         IF (ncount(n)>0) s=s*ncount(n)
      ENDDO 
      IF (s.NE.SIZE(DATA)) CALL hdf_err("mismatch of sizes") 
      !DO IO                                                            
      CALL h5dget_space_f(did,fspace,hdferr) 
      CALL h5sselect_hyperslab_f(fspace,H5S_SELECT_SET_F,foffset,fncount &
     &     ,hdferr)                                                     
      CALL h5screate_simple_f(4,dims,memspace,hdferr) 
      CALL h5dread_f(did,H5T_NATIVE_INTEGER,DATA,dims,hdferr,memspace   &
     &     ,fspace,trans)                                               

      CALL io_check("io_read_integer4 !",hdferr)
      CALL h5sclose_f(memspace,hdferr) 
      CALL h5sclose_f(fspace,hdferr) 
      CALL cleartransprop(trans) 

      CALL io_check("io_read_integer4  !",hdferr)
                                                                        
      END SUBROUTINE 
                                                                        
      SUBROUTINE io_read_integer5(did,start,ncount,DATA,transprop)
!*****************************************************************      
      USE hdf5 
                                                                        
      IMPLICIT NONE 
      INTEGER(HID_T),INTENT(IN)  ::did 
      INTEGER,INTENT(IN)         ::start(:),                            &
     &     ncount(:)
      INTEGER,INTENT(OUT)    ::DATA(:,:,:,:,:) 
      INTEGER(HID_T),INTENT(IN),OPTIONAL::transprop 
      !locals                                                           
      INTEGER(HSIZE_t)::dims(7),foffset(SIZE(start)),fncount(SIZE(ncount))
      INTEGER(HID_t)::trans,fspace,memspace 
      INTEGER       ::hdferr,s,n 
      IF (.NOT.PRESENT(transprop)) THEN 
         trans=gettransprop() 
      ELSE 
         trans=transprop 
      ENDIF 
                       !                                                
      foffset=start-1 
                       !                                                
      fncount=ncount
                                                                 !      
      dims=(/SIZE(DATA,1),SIZE(DATA,2),SIZE(DATA,3),SIZE(DATA,4)        &
     &     ,SIZE(DATA,5),0,0/)                                          
      !check if size of ncount is ok!
                               !read nothing                            
      IF (ANY(ncount<1)) RETURN
      s=1 
      DO n=1,SIZE(ncount)
         IF (ncount(n)>0) s=s*ncount(n)
      ENDDO 
      IF (s.NE.SIZE(DATA)) CALL hdf_err("mismatch of sizes") 
      !DO IO                                                            
      CALL h5dget_space_f(did,fspace,hdferr) 
      CALL h5sselect_hyperslab_f(fspace,H5S_SELECT_SET_F,foffset,fncount &
     &     ,hdferr)                                                     
      CALL h5screate_simple_f(5,dims,memspace,hdferr) 
      CALL h5dread_f(did,H5T_NATIVE_INTEGER,DATA,dims,hdferr,memspace   &
     &     ,fspace,trans)                                               

      CALL io_check("io_read_integer5  !",hdferr)
      CALL h5sclose_f(memspace,hdferr) 
      CALL h5sclose_f(fspace,hdferr) 
      CALL cleartransprop(trans) 

      CALL io_check("io_read_integer5 !",hdferr)
      END SUBROUTINE 
                                                                        
      SUBROUTINE io_read_integer6(did,start,ncount,DATA,transprop)
!*****************************************************************      
      USE hdf5 
                                                                        
      IMPLICIT NONE 
      INTEGER(HID_T),INTENT(IN)  ::did 
      INTEGER,INTENT(IN)         ::start(:),                            &
     &     ncount(:)
      INTEGER,INTENT(OUT)    ::DATA(:,:,:,:,:,:) 
      INTEGER(HID_T),INTENT(IN),OPTIONAL::transprop 
      !locals                                                           
      INTEGER(HSIZE_t)::dims(7),foffset(SIZE(start)),fncount(SIZE(ncount))
      INTEGER(HID_t)::trans,fspace,memspace 
      INTEGER       ::hdferr,s,n 
      IF (.NOT.PRESENT(transprop)) THEN 
         trans=gettransprop() 
      ELSE 
         trans=transprop 
      ENDIF 
                       !                                                
      foffset=start-1 
                       !                                                
      fncount=ncount
                                                                    !   
      dims(:)=(/SIZE(DATA,1),SIZE(DATA,2),SIZE(DATA,3),SIZE(DATA,4)     &
     &     ,SIZE(DATA,5),SIZE(DATA,6),0/)                               
      !check if size of ncount is ok!
                               !read nothing                            
      IF (ANY(ncount<1)) RETURN
      s=1 
      DO n=1,SIZE(ncount)
         IF (ncount(n)>0) s=s*ncount(n)
      ENDDO 
      IF (s.NE.SIZE(DATA)) CALL hdf_err("mismatch of sizes") 
      !DO IO                                                            
      CALL h5dget_space_f(did,fspace,hdferr) 
      CALL h5sselect_hyperslab_f(fspace,H5S_SELECT_SET_F,foffset,fncount &
     &     ,hdferr)                                                     
      CALL h5screate_simple_f(6,dims,memspace,hdferr) 
      CALL h5dread_f(did,H5T_NATIVE_INTEGER,DATA,dims,hdferr,memspace   &
     &     ,fspace,trans)                                               

      CALL io_check("io_read_integer6 !",hdferr)
      CALL h5sclose_f(memspace,hdferr) 
      CALL h5sclose_f(fspace,hdferr) 
      CALL cleartransprop(trans) 

      CALL io_check("io_read_integer6 !",hdferr)
      END SUBROUTINE 
!---------------------------------------------------------------------- 
      SUBROUTINE io_write_integer0(did,start,ncount,DATA,transprop)
!*****************************************************************      
      USE hdf5 
                                                                        
      IMPLICIT NONE 
                                                                        
! arguments                                                             
                                                       ! setid          
      INTEGER(HID_T),INTENT(IN)  :: did 
      INTEGER,       INTENT(IN)  :: start(:),ncount(:)
      INTEGER,       INTENT(IN)  :: DATA 
      INTEGER(HID_T),INTENT(IN),OPTIONAL:: transprop 
                                                                        
! locals                                                                
      INTEGER(HSIZE_T) :: dims(7) 
      INTEGER(HSIZE_T) :: foffset(SIZE(start)),fncount(SIZE(ncount))
      INTEGER(HID_T)   :: trans,fspace,memspace 
      INTEGER          :: hdferr,s,n 
                                                                        
                               ! write a single integer                 
      dims = (/1,0,0,0,0,0,0/) 
                                                                        
      IF (.NOT.PRESENT(transprop)) THEN 
                                         ! set default                  
         trans = gettransprop() 
                                         ! transfer-properties          
      ELSE 
         trans = transprop 
      ENDIF 
      foffset = start-1 
      fncount  = ncount
                               ! write nothing                          
      IF (ANY(ncount<1)) RETURN
                               ! check if size of ncount is ok
      s=1 
      DO n = 1, SIZE(ncount)
         IF (ncount(n)>0) s=s*ncount(n)
      ENDDO 
      IF (s.NE.1) CALL hdf_err("mismatch of sizes") 
                                                                        
!do I/O                                                                 
      CALL h5dget_space_f(                                              &
     &                    did,fspace,hdferr)                            
      CALL h5sselect_hyperslab_f(                                       &
     &                           fspace,H5S_SELECT_SET_F,               &
     &                           foffset,fncount,hdferr)
      CALL h5screate_simple_f(                                          &
     &                        1,dims,memspace,hdferr)                   
      CALL h5dwrite_f(                                                  &
     &                did,H5T_NATIVE_INTEGER,DATA,dims,hdferr,          &
     &                memspace,fspace,trans)                            
      CALL h5sclose_f(                                                  &
     &                memspace,hdferr)                                  
      CALL h5sclose_f(                                                  &
     &                fspace,hdferr)                                    
      CALL cleartransprop(trans) 
      CALL io_check("io_write_integer0:",hdferr) 
                                                                        
      END SUBROUTINE  io_write_integer0 
                                                                        
!---------------------------------------------------------------------- 
      SUBROUTINE io_write_integer1(did,start,ncount,DATA,transprop)
!*****************************************************************      
      USE hdf5 
                                                                        
      IMPLICIT NONE 
                                                                        
!arguments                                                              
                                                       ! setid          
      INTEGER(HID_T),INTENT(IN)  :: did 
      INTEGER,       INTENT(IN)  :: start(:),ncount(:)
                                                       ! 1-dim          
      INTEGER,       INTENT(IN)  :: DATA(:) 
      INTEGER(HID_T),INTENT(IN),OPTIONAL:: transprop 
                                                                        
! locals                                                                
      INTEGER(HSIZE_T) :: dims(7) 
      INTEGER(HSIZE_T) :: foffset(SIZE(start)),fncount(SIZE(ncount))
      INTEGER(HID_T)   :: trans,fspace,memspace 
      INTEGER          :: hdferr,s,n 
                                                                        
                                            ! write 1-dim array         
      dims = (/SIZE(DATA,1),0,0,0,0,0,0/) 
                                                                        
      IF (.NOT.PRESENT(transprop)) THEN 
         trans = gettransprop() 
      ELSE 
         trans = transprop 
      ENDIF 
      foffset = start-1 
      fncount  = ncount
                               ! write nothing                          
      IF (ANY(ncount<1)) RETURN
                               ! check if size of ncount is ok
      s=1 
      DO n = 1, SIZE(ncount)
         IF (ncount(n)>0) s=s*ncount(n)
      ENDDO 
      IF (s.NE.SIZE(DATA)) CALL hdf_err("mismatch of sizes") 
                                                                        
!do I/O                                                                 
      CALL h5dget_space_f(                                              &
     &                    did,fspace,hdferr)                            
      CALL h5sselect_hyperslab_f(                                       &
     &                           fspace,H5S_SELECT_SET_F,               &
     &                           foffset,fncount,hdferr)
      CALL h5screate_simple_f(                                          &
     &                        1,dims,memspace,hdferr)                   
      CALL h5dwrite_f(                                                  &
     &                did,H5T_NATIVE_INTEGER,DATA,dims,hdferr,          &
     &                memspace,fspace,trans)                            
      CALL h5sclose_f(                                                  &
     &                memspace,hdferr)                                  
      CALL h5sclose_f(                                                  &
     &                fspace,hdferr)                                    
      CALL cleartransprop(trans) 
      CALL io_check("io_write_integer1:",hdferr) 
                                                                        
      END SUBROUTINE  io_write_integer1 
                                                                        
!---------------------------------------------------------------------- 
                                                                        
      SUBROUTINE io_write_integer2(did,start,ncount,DATA,transprop)
!*****************************************************************      
      USE hdf5 
                                                                        
      IMPLICIT NONE 
      INTEGER(HID_T),INTENT(IN)  ::did 
      INTEGER,INTENT(IN)         ::start(:),                            &
     &     ncount(:)
      INTEGER,INTENT(IN)    ::DATA(:,:) 
      INTEGER(HID_T),INTENT(IN),OPTIONAL::transprop 
      !locals                                                           
      INTEGER(HSIZE_t)::dims(7),foffset(SIZE(start)),fncount(SIZE(ncount))
      INTEGER(HID_t)::trans,fspace,memspace 
      INTEGER       ::hdferr,s,n 
      IF (.NOT.PRESENT(transprop)) THEN 
         trans=gettransprop() 
      ELSE 
         trans=transprop 
      ENDIF 
                       !                                                
      foffset=start-1 
                       !                                                
      fncount=ncount
                                                     !                  
      dims=(/SIZE(DATA,1),SIZE(DATA,2),0,0,0,0,0/) 
      !check if size of ncount is ok!
                               !read nothing                            
      IF (ANY(ncount<1)) RETURN
      s=1 
      DO n=1,SIZE(ncount)
         IF (ncount(n)>0) s=s*ncount(n)
      ENDDO 
      IF (s.NE.SIZE(DATA)) CALL hdf_err("mismatch of sizes") 
      !DO IO                                                            
      CALL h5dget_space_f(did,fspace,hdferr) 
      CALL h5sselect_hyperslab_f(fspace,H5S_SELECT_SET_F,foffset,fncount &
     &     ,hdferr)                                                     
      CALL h5screate_simple_f(2,dims,memspace,hdferr) 
      CALL h5dwrite_f(did,H5T_NATIVE_INTEGER,DATA,dims,hdferr,memspace  &
     &     ,fspace,trans)                                               
      CALL h5sclose_f(memspace,hdferr) 
      CALL h5sclose_f(fspace,hdferr) 
      CALL cleartransprop(trans) 
      CALL io_check("io_write_integer2:",hdferr) 
                                                                        
      END SUBROUTINE 
                                                                        
      SUBROUTINE io_write_integer3(did,start,ncount,DATA,transprop)
!*****************************************************************      
      USE hdf5 
                                                                        
      IMPLICIT NONE 
      INTEGER(HID_T),INTENT(IN)  ::did 
      INTEGER,INTENT(IN)         ::start(:),                            &
     &     ncount(:)
      INTEGER,INTENT(IN)    ::DATA(:,:,:) 
      INTEGER(HID_T),INTENT(IN),OPTIONAL::transprop 
      !locals                                                           
      INTEGER(HSIZE_t)::dims(7),foffset(SIZE(start)),fncount(SIZE(ncount))
      INTEGER(HID_t)::trans,fspace,memspace 
      INTEGER       ::hdferr,s,n 
      IF (.NOT.PRESENT(transprop)) THEN 
         trans=gettransprop() 
      ELSE 
         trans=transprop 
      ENDIF 
                       !                                                
      foffset=start-1 
                       !                                                
      fncount=ncount
                                                                !       
      dims=(/SIZE(DATA,1),SIZE(DATA,2),SIZE(DATA,3),0,0,0,0/) 
      !check if size of ncount is ok!
                               !read nothing                            
      IF (ANY(ncount<1)) RETURN
      s=1 
      DO n=1,SIZE(ncount)
         IF (ncount(n)>0) s=s*ncount(n)
      ENDDO 
      IF (s.NE.SIZE(DATA)) CALL hdf_err("mismatch of sizes") 
      !DO IO                                                            
      CALL h5dget_space_f(did,fspace,hdferr) 
      CALL h5sselect_hyperslab_f(fspace,H5S_SELECT_SET_F,foffset,fncount &
     &     ,hdferr)                                                     
      CALL h5screate_simple_f(3,dims,memspace,hdferr) 
      CALL h5dwrite_f(did,H5T_NATIVE_INTEGER,DATA,dims,hdferr,memspace  &
     &     ,fspace,trans)                                               
      CALL h5sclose_f(memspace,hdferr) 
      CALL h5sclose_f(fspace,hdferr) 
      CALL cleartransprop(trans) 
      CALL io_check("io_write_integer3:",hdferr) 
      END SUBROUTINE 
                                                                        
      SUBROUTINE io_write_integer4(did,start,ncount,DATA,transprop)
!*****************************************************************      
      USE hdf5 
                                                                        
      IMPLICIT NONE 
      INTEGER(HID_T),INTENT(IN)  ::did 
      INTEGER,INTENT(IN)         ::start(:),                            &
     &     ncount(:)
      INTEGER,INTENT(IN)    ::DATA(:,:,:,:) 
      INTEGER(HID_T),INTENT(IN),OPTIONAL::transprop 
      !locals                                                           
      INTEGER(HSIZE_t)::dims(7),foffset(SIZE(start)),fncount(SIZE(ncount))
      INTEGER(HID_t)::trans,fspace,memspace 
      INTEGER       ::hdferr,s,n 
      IF (.NOT.PRESENT(transprop)) THEN 
         trans=gettransprop() 
      ELSE 
         trans=transprop 
      ENDIF 
                       !                                                
      foffset=start-1 
                       !                                                
      fncount=ncount
      dims=(/SIZE(DATA,1),SIZE(DATA,2),SIZE(DATA,3),SIZE(DATA,4),0,0,0/) 
      !check if size of ncount is ok!
                               !read nothing                            
      IF (ANY(ncount<1)) RETURN
      s=1 
      DO n=1,SIZE(ncount)
         IF (ncount(n)>0) s=s*ncount(n)
      ENDDO 
      IF (s.NE.SIZE(DATA)) CALL hdf_err("mismatch of sizes") 
      !DO IO                                                            
      CALL h5dget_space_f(did,fspace,hdferr) 
      CALL h5sselect_hyperslab_f(fspace,H5S_SELECT_SET_F,foffset,fncount &
     &     ,hdferr)                                                     
      CALL h5screate_simple_f(4,dims,memspace,hdferr) 
      CALL h5dwrite_f(did,H5T_NATIVE_INTEGER,DATA,dims,hdferr,memspace  &
     &     ,fspace,trans)                                               
      CALL h5sclose_f(memspace,hdferr) 
      CALL h5sclose_f(fspace,hdferr) 
      CALL cleartransprop(trans) 
      CALL io_check("io_write_integer4:",hdferr) 
      END SUBROUTINE 
                                                                        
      SUBROUTINE io_write_integer5(did,start,ncount,DATA,transprop)
!*****************************************************************      
      USE hdf5 
                                                                        
      IMPLICIT NONE 
      INTEGER(HID_T),INTENT(IN)  ::did 
      INTEGER,INTENT(IN)         ::start(:),                            &
     &     ncount(:)
      INTEGER,INTENT(IN)    ::DATA(:,:,:,:,:) 
      INTEGER(HID_T),INTENT(IN),OPTIONAL::transprop 
      !locals                                                           
      INTEGER(HSIZE_t)::dims(7),foffset(SIZE(start)),fncount(SIZE(ncount))
      INTEGER(HID_t)::trans,fspace,memspace 
      INTEGER       ::hdferr,s,n 
      IF (.NOT.PRESENT(transprop)) THEN 
         trans=gettransprop() 
      ELSE 
         trans=transprop 
      ENDIF 
                       !                                                
      foffset=start-1 
                       !                                                
      fncount=ncount
                                                                 !      
      dims=(/SIZE(DATA,1),SIZE(DATA,2),SIZE(DATA,3),SIZE(DATA,4)        &
     &     ,SIZE(DATA,5),0,0/)                                          
      !check if size of ncount is ok!
                               !read nothing                            
      IF (ANY(ncount<1)) RETURN
      s=1 
      DO n=1,SIZE(ncount)
         IF (ncount(n)>0) s=s*ncount(n)
      ENDDO 
      IF (s.NE.SIZE(DATA)) CALL hdf_err("mismatch of sizes") 
      !DO IO                                                            
      CALL h5dget_space_f(did,fspace,hdferr) 
      CALL h5sselect_hyperslab_f(fspace,H5S_SELECT_SET_F,foffset,fncount &
     &     ,hdferr)                                                     
      CALL h5screate_simple_f(5,dims,memspace,hdferr) 
      CALL h5dwrite_f(did,H5T_NATIVE_INTEGER,DATA,dims,hdferr,memspace  &
     &     ,fspace,trans)                                               
      CALL h5sclose_f(memspace,hdferr) 
      CALL h5sclose_f(fspace,hdferr) 
      CALL cleartransprop(trans) 
      CALL io_check("io_write_integer5:",hdferr) 
      END SUBROUTINE 
                                                                        
      SUBROUTINE io_write_integer6(did,start,ncount,DATA,transprop)
!*****************************************************************      
      USE hdf5 
                                                                        
      IMPLICIT NONE 
      INTEGER(HID_T),INTENT(IN)  ::did 
      INTEGER,INTENT(IN)         ::start(:),                            &
     &     ncount(:)
      INTEGER,INTENT(IN)    ::DATA(:,:,:,:,:,:) 
      INTEGER(HID_T),INTENT(IN),OPTIONAL::transprop 
      !locals                                                           
      INTEGER(HSIZE_t)::dims(7),foffset(SIZE(start)),fncount(SIZE(ncount))
      INTEGER(HID_t)::trans,fspace,memspace 
      INTEGER       ::hdferr,s,n 
      IF (.NOT.PRESENT(transprop)) THEN 
         trans=gettransprop() 
      ELSE 
         trans=transprop 
      ENDIF 
                       !                                                
      foffset=start-1 
                       !                                                
      fncount=ncount
                                                                    !   
      dims(:)=(/SIZE(DATA,1),SIZE(DATA,2),SIZE(DATA,3),SIZE(DATA,4)     &
     &     ,SIZE(DATA,5),SIZE(DATA,6),0/)                               
      !check if size of ncount is ok!
                               !read nothing                            
      IF (ANY(ncount<1)) RETURN
      s=1 
      DO n=1,SIZE(ncount)
         IF (ncount(n)>0) s=s*ncount(n)
      ENDDO 
      IF (s.NE.SIZE(DATA)) CALL hdf_err("mismatch of sizes") 
      !DO IO                                                            
      CALL h5dget_space_f(did,fspace,hdferr) 
      CALL h5sselect_hyperslab_f(fspace,H5S_SELECT_SET_F,foffset,fncount &
     &     ,hdferr)                                                     
      CALL h5screate_simple_f(6,dims,memspace,hdferr) 
      CALL h5dwrite_f(did,H5T_NATIVE_INTEGER,DATA,dims,hdferr,memspace  &
     &     ,fspace,trans)                                               
      CALL h5sclose_f(memspace,hdferr) 
      CALL h5sclose_f(fspace,hdferr) 
      CALL cleartransprop(trans) 
      CALL io_check("io_write_integer6:",hdferr) 
      END SUBROUTINE 
!*****************************************************************      
!                                                                       
!     The following subroutines READ or WRITE COMPLEX values            
!     from hdf-file                                                     
!                                                                       
!*****************************************************************      
                                                                        
      SUBROUTINE io_read_complex0(did,start,ncount,DATA,trans)
!*****************************************************************      
      USE hdf5 
                                                                        
      IMPLICIT NONE 
      INTEGER(HID_T),INTENT(IN)  ::did 
      INTEGER,INTENT(IN)         ::start(:),                            &
     &     ncount(:)
      COMPLEX(rkind),INTENT(OUT)    ::DATA 
      INTEGER(HID_T),INTENT(IN),OPTIONAL::trans 
      !locals                                                           
      INTEGER::foffset(SIZE(start)),fncount(SIZE(start))
      REAL          ::a(2)
                     !                                                  
      foffset=start
      fncount=start
      !DO 2 calls to read real values                                   
      WHERE(start<0) 
         foffset=1
         fncount=2
      ENDWHERE 
      CALL io_read_real1(did,foffset,fncount,a,trans)
      DATA=CMPLX(a(1),a(2))
      END SUBROUTINE 
      SUBROUTINE io_read_complex1(did,start,ncount,DATA,trans)
!*****************************************************************      
      USE hdf5 
                                                                        
      IMPLICIT NONE 
      INTEGER(HID_T),INTENT(IN)  ::did 
      INTEGER,INTENT(IN)         ::start(:),                            &
     &     ncount(:)
      COMPLEX(rkind),INTENT(OUT)    ::DATA(:) 
      INTEGER(HID_T),INTENT(IN),OPTIONAL::trans 
      !locals                                                           
      INTEGER::foffset(SIZE(start)),fncount(SIZE(start))
      REAL,ALLOCATABLE::a(:,:)

      foffset=start
      fncount=ncount
      WHERE(start<0) 
         foffset=1
         fncount=2
      ENDWHERE 

      if (start(1)<0) then
          ALLOCATE(A(2,SIZE(DATA,1)))
          CALL io_read_real2(did,foffset,fncount,a,trans)
          DATA=CMPLX(a(1,:),a(2,:))
      else
          ALLOCATE(A(SIZE(DATA,1),2))
          CALL io_read_real2(did,foffset,fncount,a,trans)
          DATA=CMPLX(a(:,1),a(:,2))
      endif
      DEALLOCATE(a)
      END SUBROUTINE 
!---------------------------------------------------------------------- 
                                                                        
      SUBROUTINE io_read_complex2(did,start,ncount,DATA,trans)
!*****************************************************************      
      USE hdf5 
      IMPLICIT NONE 
! arguments                                                             
      INTEGER(HID_T),INTENT(IN)  :: did 
      INTEGER,       INTENT(IN)  :: start(:),ncount(:)
      COMPLEX,       INTENT(OUT) :: DATA(:,:) 
      INTEGER(HID_T),INTENT(IN),OPTIONAL :: trans 
! locals                                                                
      INTEGER           :: foffset(SIZE(start)),fncount(size(start))
      REAL, ALLOCATABLE :: a(:,:,:)
                                                                        
                     !                                                  
      foffset=start 
      fncount=ncount
      !DO 2 calls to read real values                                   
      WHERE(start<0) 
         foffset=1 
         fncount=2
      ENDWHERE
      if (count(start<0).ne.1) CPP_error("Wrong no of negatives")
      if (start(1)<0) then
          ALLOCATE(A(2,SIZE(DATA,1),size(data,2)))
          CALL io_read_real3(did,foffset,fncount,a,trans)
          DATA=CMPLX(a(1,:,:),a(2,:,:))
      else if (start(2)<0) then
          ALLOCATE(A(SIZE(DATA,1),2,size(data,2)))
          CALL io_read_real3(did,foffset,fncount,a,trans)
          DATA=CMPLX(a(:,1,:),a(:,2,:))
      else
          ALLOCATE(A(SIZE(DATA,1),size(data,2),2))
          CALL io_read_real3(did,foffset,fncount,a,trans)
          DATA=CMPLX(a(:,:,1),a(:,:,2))
      endif
      DEALLOCATE(a)
      END SUBROUTINE 
      SUBROUTINE io_read_complex3(did,start,ncount,DATA,trans)
!*****************************************************************      
      USE hdf5 
                                                                        
      IMPLICIT NONE 
      INTEGER(HID_T),INTENT(IN)  ::did 
      INTEGER,INTENT(IN)         ::start(:),                            &
     &     ncount(:)
      COMPLEX(rkind),INTENT(OUT)    ::DATA(:,:,:) 
      INTEGER(HID_T),INTENT(IN),OPTIONAL::trans 
      !locals                                                           
      INTEGER::foffset(SIZE(start)),fncount(size(start))
      REAL,ALLOCATABLE::a(:,:,:,:)

                     !                                                  
      foffset=start 
      fncount=ncount
      !DO 2 calls to read real values                                   
      WHERE(start<0) 
         foffset=1 
         fncount=2
      ENDWHERE 
      if (count(start<0).ne.1) CPP_error("Wrong no of negatives")
      if (start(1)<0) then
          ALLOCATE(A(2,SIZE(DATA,1),size(data,2),size(data,3)))
          CALL io_read_real4(did,foffset,fncount,a,trans)
          DATA=CMPLX(a(1,:,:,:),a(2,:,:,:))
      else if (start(4)<0) then
          ALLOCATE(A(SIZE(DATA,1),size(data,2),size(data,3),2))
          CALL io_read_real4(did,foffset,fncount,a,trans)
          DATA=CMPLX(a(:,:,:,1),a(:,:,:,2))
      else
          CPP_error("Wrong position of negative")
      endif
      DEALLOCATE(a)
      END SUBROUTINE 
      SUBROUTINE io_read_complex4(did,start,ncount,DATA,trans)
!*****************************************************************      
      USE hdf5 
                                                                        
      IMPLICIT NONE 
      INTEGER(HID_T),INTENT(IN)  ::did 
      INTEGER,INTENT(IN)         ::start(:),                            &
     &     ncount(:)
      COMPLEX(rkind),INTENT(OUT)    ::DATA(:,:,:,:) 
      INTEGER(HID_T),INTENT(IN),OPTIONAL::trans 
      !locals                                                           
      INTEGER::foffset(SIZE(start)),fncount(SIZE(start))
      REAL,ALLOCATABLE::a(:,:,:,:,:)

                     !                                                  
      foffset=start 
      fncount=ncount
      !DO 2 calls to read real values                                   
      WHERE(start<0) 
         foffset=1
         fncount=2
      ENDWHERE 
      if (count(start<0)/=1) CPP_error("Wrong no of negatives")
      if (start(1)<0) then
          ALLOCATE(A(2,SIZE(DATA,1),size(data,2),size(data,3),size(data,4)))
          CALL io_read_real5(did,foffset,fncount,a,trans)
          DATA=CMPLX(a(1,:,:,:,:),a(2,:,:,:,:))
      else if (start(5)<0) then
          ALLOCATE(A(SIZE(DATA,1),size(data,2),size(data,3),size(data,4),2))
          CALL io_read_real5(did,foffset,fncount,a,trans)
          DATA=CMPLX(a(:,:,:,:,1),a(:,:,:,:,2))
      else
          CPP_error("Wrong position of negative")
      endif
      DEALLOCATE(a)
      END SUBROUTINE 
      SUBROUTINE io_read_complex5(did,start,ncount,DATA,trans)
!*****************************************************************      
      USE hdf5 
                                                                        
      IMPLICIT NONE 
      INTEGER(HID_T),INTENT(IN)  ::did 
      INTEGER,INTENT(IN)         ::start(:),                            &
     &     ncount(:)
      COMPLEX(rkind),INTENT(OUT)    ::DATA(:,:,:,:,:) 
      INTEGER(HID_T),INTENT(IN),OPTIONAL::trans 
      !locals
      INTEGER::foffset(SIZE(start)),fncount(SIZE(start))
      REAL,ALLOCATABLE::a(:,:,:,:,:,:)

                     !                                                  
      foffset=start 
      fncount=ncount
      !DO 2 calls to read real values                                   
      WHERE(start<0) 
         foffset=1
         fncount=2
      ENDWHERE 
      if (count(start<0)/=1) CPP_error("Wrong no of negatives")
      if (start(1)<0) then
          ALLOCATE(A(2,SIZE(DATA,1),size(data,2),size(data,3),size(data,4),size(data,5)))
          CALL io_read_real6(did,foffset,fncount,a,trans)
          DATA=CMPLX(a(1,:,:,:,:,:),a(2,:,:,:,:,:))
      else if (start(6)<0) then
          ALLOCATE(A(SIZE(DATA,1),size(data,2),size(data,3),size(data,4),size(data,5),2))
          CALL io_read_real6(did,foffset,fncount,a,trans)
          DATA=CMPLX(a(:,:,:,:,:,1),a(:,:,:,:,:,2))
      else
          CPP_error("Wrong position of negative")
      endif
      DEALLOCATE(a)
      END SUBROUTINE 
      SUBROUTINE io_write_complex0(did,start,ncount,DATA,trans)
!*****************************************************************      
      USE hdf5 
                                                                        
      IMPLICIT NONE 
      INTEGER(HID_T),INTENT(IN)  ::did 
      INTEGER,INTENT(IN)         ::start(:),                            &
     &     ncount(:)
      COMPLEX(rkind),INTENT(IN)    ::DATA 
      INTEGER(HID_T),INTENT(IN),OPTIONAL::trans 
      !locals
      real   ::a(2)
      INTEGER::foffset(SIZE(start)),fncount(SIZE(start))
      fncount=ncount
      foffset=start 

      !DO 2 calls to read real values                                   
      WHERE(start<0) 
         foffset=1
         fncount=2
      ENDWHERE 
      a(1)=real(data)
      a(2)=aimag(data)
      CALL io_write_real1(did,foffset,fncount,a,trans)

      END SUBROUTINE 
      SUBROUTINE io_write_complex1(did,start,ncount,DATA,trans)
!*****************************************************************      
      USE hdf5 
                                                                        
      IMPLICIT NONE 
      INTEGER(HID_T),INTENT(IN)  ::did 
      INTEGER,INTENT(IN)         ::start(:),                            &
     &     ncount(:)
      COMPLEX(rkind),INTENT(IN)    ::DATA(:) 
      INTEGER(HID_T),INTENT(IN),OPTIONAL::trans 
      !locals                                                           
      INTEGER::foffset(SIZE(start)),fncount(SIZE(start))
      real::a(2,size(data))
      fncount=ncount                                !
      foffset=start 
      a(1,:)=real(data)
      a(2,:)=aimag(data)
      !DO 2 calls to read real values                                   
      WHERE(start<0) 
         foffset=1
         fncount=2
      ENDWHERE 
      CALL io_write_real2(did,foffset,fncount,a,trans)
      END SUBROUTINE 
      SUBROUTINE io_write_complex2(did,start,ncount,DATA,trans)
!*****************************************************************      
      USE hdf5 
                                                                        
      IMPLICIT NONE 
      INTEGER(HID_T),INTENT(IN)  ::did 
      INTEGER,INTENT(IN)         ::start(:),                            &
     &     ncount(:)
      COMPLEX(rkind),INTENT(IN)    ::DATA(:,:) 
      INTEGER(HID_T),INTENT(IN),OPTIONAL::trans 
      !locals                                                           
      INTEGER::foffset(SIZE(start)),fncount(SIZE(start))
      real::a(2,size(data,1),size(data,2))
      fncount=ncount
                     !                                                  
      foffset=start 
      !DO 2 calls to read real values                                   
      WHERE(start<0) 
         foffset=1 
         fncount=2
      ENDWHERE 
      a(1,:,:)=real(data)
      a(2,:,:)=aimag(data)
      CALL io_write_real3(did,foffset,fncount,a,trans)
      END SUBROUTINE 

      SUBROUTINE io_write_complex3(did,start,ncount,DATA,trans)
!*****************************************************************      
      USE hdf5 
                                                                        
      IMPLICIT NONE 
      INTEGER(HID_T),INTENT(IN)  ::did 
      INTEGER,INTENT(IN)         ::start(:),                            &
     &     ncount(:)
      COMPLEX(rkind),INTENT(IN)    ::DATA(:,:,:) 
      INTEGER(HID_T),INTENT(IN),OPTIONAL::trans 
      !locals                                                           
      INTEGER::foffset(SIZE(start)),fncount(SIZE(start))
      real::a(2,size(data,1),size(data,2),size(data,3))

      a(1,:,:,:)=real(data)
      a(2,:,:,:)=aimag(data)

      fncount=ncount
                     !                                                  
      foffset=start 
      !DO 2 calls to read real values                                   
      WHERE(start<0) 
         foffset=1
         fncount=2
      ENDWHERE 
      CALL io_write_real4(did,foffset,fncount,a,trans)
      END SUBROUTINE 
      SUBROUTINE io_write_complex4(did,start,ncount,DATA,trans)
!*****************************************************************      
      USE hdf5 
                                                                        
      IMPLICIT NONE 
      INTEGER(HID_T),INTENT(IN)  ::did 
      INTEGER,INTENT(IN)         ::start(:),                            &
     &     ncount(:)
      COMPLEX(rkind),INTENT(IN)    ::DATA(:,:,:,:) 
      INTEGER(HID_T),INTENT(IN),OPTIONAL::trans 
      !locals                                                           
      INTEGER::foffset(SIZE(start)),fncount(SIZE(start))
      real::a(2,size(data,1),size(data,2),size(data,3),size(data,4))

      a(1,:,:,:,:)=real(data)
      a(2,:,:,:,:)=aimag(data)

      fncount=ncount
                     !                                                  
      foffset=start 
      !DO 2 calls to read real values                                   
      WHERE(start<0) 
         foffset=1
         fncount=2
      ENDWHERE 
      CALL io_write_real5(did,foffset,fncount,a,trans)
      END SUBROUTINE 
      SUBROUTINE io_write_complex5(did,start,ncount,DATA,trans)
!*****************************************************************      
      USE hdf5 
                                                                        
      IMPLICIT NONE 
      INTEGER(HID_T),INTENT(IN)  ::did 
      INTEGER,INTENT(IN)         ::start(:),                            &
     &     ncount(:)
      COMPLEX(rkind),INTENT(IN)    ::DATA(:,:,:,:,:) 
      INTEGER(HID_T),INTENT(IN),OPTIONAL::trans 
      !locals                                                           
      INTEGER::foffset(SIZE(start)),fncount(SIZE(start))
      real::a(2,size(data,1),size(data,2),size(data,3),size(data,4),size(data,5))

      a(1,:,:,:,:,:)=real(data)
      a(2,:,:,:,:,:)=aimag(data)

      fncount=ncount
                     !                                                  
      foffset=start 
      !DO 2 calls to read real values                                   
      WHERE(start<0) 
         foffset=1
         fncount=2
      ENDWHERE 
      CALL io_write_real6(did,foffset,fncount,a,trans)
      END SUBROUTINE 
                                                                        
      !>                                                                
      END                                           
