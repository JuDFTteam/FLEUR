!--------------------------------------------------------------------------------
! Copyright (c) 2016 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------

      MODULE m_hdf_tools5 
!-----------------------------------------------                        
!     major rewrite of hdf_tools                                        
!     this module contains only the                                     
!     IO-with stride                                                    
!-----------------------------------------------                        
      USE m_hdf_tools4 
      PRIVATE 
      PUBLIC io_write_real1s,io_write_real2s,io_write_real3s 
      CONTAINS 
      !<--subroutines for IO with stride                                
      SUBROUTINE io_write_real3s(did,start,count,data,stride,transprop) 
!*****************************************************************      
      USE hdf5 
      IMPLICIT NONE 
                                                                        
      INTEGER(HID_T),INTENT(IN) :: did 
      INTEGER,       INTENT(IN) :: start(:),count(:),stride(:) 
      REAL,          INTENT(IN) :: data(:,:,:) 
      INTEGER(HID_T),INTENT(IN),OPTIONAL :: transprop 
!  locals                                                               
      INTEGER(HSIZE_t) :: dims(7),foffset(size(start)) 
      INTEGER(HSIZE_t) :: fcount(size(count)) 
      INTEGER(HSIZE_t) :: fstride(SIZE(stride)) 
      INTEGER(HID_t)::trans,fspace,memspace 
      INTEGER       ::hdferr,s,n 
                                                                        
      IF (.not.present(transprop)) THEN 
         trans=gettransprop() 
      ELSE 
         trans=transprop 
      ENDIF 
                       !                                                
      foffset=start-1 
      fstride = stride 
                       !                                                
      fcount=count 
                                                                !       
      dims=(/size(data,1),size(data,2),size(data,3),0,0,0,0/) 
! check if size of count is ok!                                         
                               !read nothing                            
      if (any(count<1)) return 
      s=1 
      DO n=1,size(count) 
         IF (count(n)>0) s=s*count(n) 
      ENDDO 
      IF (s.ne.size(data)) CALL hdf_err("Missmatch of sizes") 
! DO IO                                                                 
      CALL h5dget_space_f(                                              &
     &                    did,fspace,hdferr)                            
      CALL h5sselect_hyperslab_f(                                       &
     &                           fspace,H5S_SELECT_SET_F,               &
     &                           foffset,fcount,                        &
     &                           hdferr,fstride)                        
      CALL h5screate_simple_f(                                          &
     &                        3,dims,                                   &
     &                        memspace,hdferr)                          
      CALL h5dwrite_f(                                                  &
     &                did,H5T_NATIVE_DOUBLE,                            &
     &                data,dims,hdferr,                                 &
     &                memspace,fspace,trans)                            
      CALL h5sclose_f(                                                  &
     &                memspace,hdferr)                                  
      CALL h5sclose_f(                                                  &
     &                fspace,hdferr)                                    
      CALL cleartransprop(trans) 
                                                                        
      END SUBROUTINE io_write_real3s 
                                                                        
      SUBROUTINE io_write_real1s(                                       &
     &                           did,start,count,DATA,stride,transprop) 
!*****************************************************************      
      USE hdf5 
      IMPLICIT NONE 
                                                                        
!arguments                                                              
      INTEGER(HID_T), INTENT(IN)  :: did 
      INTEGER,        INTENT(IN)  :: start(:),COUNT(:),stride(:) 
      REAL(rkind),    INTENT(IN)  :: DATA(:) 
      INTEGER(HID_T), INTENT(IN), OPTIONAL :: transprop 
!locals                                                                 
      INTEGER(HSIZE_t) :: dims(7) 
      INTEGER(HSIZE_t) :: foffset(SIZE(start)),fcount(SIZE(count)) 
      INTEGER(HSIZE_t) :: fstride(SIZE(stride)) 
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
      fcount  = count 
      fstride = stride 
                               ! write nothing                          
      IF (ANY(count<1)) RETURN 
                               ! check if size of count is ok           
      s=1 
      DO n = 1, SIZE(count) 
         IF (COUNT(n)>0) s=s*COUNT(n) 
      ENDDO 
      IF (s.NE.SIZE(DATA)) CALL hdf_err("Missmatch of sizes") 
                                                                        
!do I/O                                                                 
      CALL h5dget_space_f(                                              &
     &                    did,                                          &
     &                    fspace,hdferr)                                
                                           ! dataset_id                 
                                           ! dataspace_id & error       
      CALL h5sselect_hyperslab_f(                                       &
     &                           fspace,H5S_SELECT_SET_F,               &
     &                           foffset,fcount,                        &
     &                           hdferr,fstride)                        
                                                           ! dataspace_i
                                                           ! starting po
                                                           ! error (out)
      CALL h5screate_simple_f(                                          &
     &                        1,dims,                                   &
     &                        memspace,hdferr)                          
                                                 ! rank & dimensions of 
                                                 ! memoryspace identifie
      CALL h5dwrite_f(                                                  &
     &                did,H5T_NATIVE_DOUBLE,                            &
     &                DATA,dims,hdferr,                                 &
     &                memspace,fspace,trans)                            
                                                 ! dataset_id, datatype_
                                                 ! data & dimensions, er
                                                 ! memoryspace_id, file-
                                                 ! Transfer property lis
      CALL h5sclose_f(                                                  &
     &                memspace,hdferr)                                  
      CALL h5sclose_f(                                                  &
     &                fspace,hdferr)                                    
      CALL cleartransprop(trans) 
                                                                        
      END SUBROUTINE io_write_real1s 
                                                                        
!---------------------------------------------------------------------- 
                                                                        
      SUBROUTINE io_write_real2s(                                       &
     &                           did,start,count,DATA,stride,transprop) 
!*****************************************************************      
      USE hdf5 
      IMPLICIT NONE 
                                                                        
!arguments                                                              
      INTEGER(HID_T), INTENT(IN)  :: did 
      INTEGER,        INTENT(IN)  :: start(:),COUNT(:),stride(:) 
      REAL(rkind),    INTENT(IN)  :: DATA(:,:) 
      INTEGER(HID_T), INTENT(IN), OPTIONAL :: transprop 
!locals                                                                 
      INTEGER(HSIZE_t) :: dims(7) 
      INTEGER(HSIZE_t) :: foffset(SIZE(start)),fcount(SIZE(count)) 
      INTEGER(HSIZE_t) :: fstride(SIZE(stride)) 
      INTEGER(HID_t)   :: trans,fspace,memspace 
      INTEGER          :: hdferr,s,n 
                                                                        
                                                       ! write 2-dim arr
      dims = (/SIZE(DATA,1),SIZE(DATA,2),0,0,0,0,0/) 
                                                                        
      IF (.NOT.PRESENT(transprop)) THEN 
         trans = gettransprop() 
      ELSE 
         trans = transprop 
      ENDIF 
      foffset = start-1 
      fcount  = count 
      fstride = stride 
                               ! write nothing                          
      IF (ANY(count<1)) RETURN 
                               ! check if size of count is ok           
      s=1 
      DO n = 1, SIZE(count) 
         IF (COUNT(n)>0) s=s*COUNT(n) 
      ENDDO 
      IF (s.NE.SIZE(DATA)) CALL hdf_err("Missmatch of sizes") 
                                                                        
!     write(*,*) 'foffset',foffset                                      
!     write(*,*) 'fcount', fcount                                       
!     write(*,*) 'fstride',fstride                                      
!do I/O                                                                 
      CALL h5dget_space_f(                                              &
     &                    did,                                          &
     &                    fspace,hdferr)                                
                                           ! dataset_id                 
                                           ! dataspace_id & error       
      CALL h5sselect_hyperslab_f(                                       &
     &                           fspace,H5S_SELECT_SET_F,               &
     &                           foffset,fcount,                        &
     &                           hdferr,fstride)                        
                                                           ! dataspace_i
                                                           ! starting po
                                                           ! error (out)
      CALL h5screate_simple_f(                                          &
     &                        2,dims,                                   &
     &                        memspace,hdferr)                          
                                                 ! rank & dimensions of 
                                                 ! memoryspace identifie
      CALL h5dwrite_f(                                                  &
     &                did,H5T_NATIVE_DOUBLE,                            &
     &                DATA,dims,hdferr,                                 &
     &                memspace,fspace,trans)                            
                                                 ! dataset_id, datatype_
                                                 ! data & dimensions, er
                                                 ! memoryspace_id, file-
                                                 ! Transfer property lis
      CALL h5sclose_f(                                                  &
     &                memspace,hdferr)                                  
      CALL h5sclose_f(                                                  &
     &                fspace,hdferr)                                    
      CALL cleartransprop(trans) 
                                                                        
      END SUBROUTINE io_write_real2s 
!---------------------------------------------------------------------- 
      !>                                                                
      END                                           
