!--------------------------------------------------------------------------------
! Copyright (c) 2016 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------

      MODULE m_hdf_tools3 
!-----------------------------------------------                        
!     major rewrite of hdf_tools                                        
!     this module contains only the                                     
!     subroutines                                                       
!                                                                       
!     io_[data/att/group]exists(id,name) to check for objects           
!                                                                       
!-----------------------------------------------                        
      CONTAINS 
      !<-- functions to test for objects                                
      !<-- F: io_dataexists(gid,name)RESULT(exists)                     
                                                                        
      FUNCTION io_dataexists(gid,name)RESULT(exist) 
!******************************************                             
!     checks if dataset called 'name' exists at position gid            
!                          D. Wortmann                                  
!******************************************                             
      USE hdf5 
      IMPLICIT NONE 
      !<--Arguments                                                     
      INTEGER(HID_T),INTENT(IN)  ::gid 
      CHARACTER,INTENT(IN)       ::name*(*) 
      LOGICAL                    ::exist 
      !>                                                                
      !<--Locals                                                        
      INTEGER(HID_t)::testid 
      INTEGER       ::hdferr 
      !>
      testid = 0
      hdferr = 0
      CALL h5eclear_f(hdferr) 
                                    !No automatic error checking!       
      CALL h5eset_auto_f(0, hdferr) 
      CALL h5dopen_f(gid,name,testid,hdferr) 
      exist=(hdferr.EQ.0) 
      CALL h5dclose_f(testid,hdferr) 
      CALL h5eclear_f(hdferr) 
                                    !Resume automatic error checking!   
      CALL h5eset_auto_f(1, hdferr) 
      END FUNCTION 
                                                                        
      !>                                                                
                                                                        
      !<-- F: io_attexists(gid,name)RESULT(exists)                      
      FUNCTION io_attexists(gid,name)RESULT(exist) 
!******************************************                             
!     checks if attribute called 'name' exists at position gid          
!                          D. Wortmann                                  
!******************************************                             
      USE hdf5 
      IMPLICIT NONE 
      !<--Arguments                                                     
      INTEGER(HID_T),INTENT(IN)  ::gid 
      CHARACTER,INTENT(IN)       ::name*(*) 
      LOGICAL                    ::exist 
      !>                                                                
      !<--Locals                                                        
      INTEGER(HID_t)::testid 
      INTEGER       ::hdferr 
      !>
      testid = 0
      hdferr = 0
      CALL h5eclear_f(hdferr) 
                                    !No automatic error checking!       
      CALL h5eset_auto_f(0, hdferr) 
      CALL h5aopen_name_f(gid,name,testid,hdferr) 
      exist=(hdferr.EQ.0) 
      CALL h5aclose_f(testid,hdferr) 
      CALL h5eclear_f(hdferr) 
                                    !Resume automatic error checking!   
      CALL h5eset_auto_f(1, hdferr) 
      END FUNCTION 
      !>                                                                
                                                                        
      !<-- F: io_groupexists(gid,name)RESULT(exists)                    
      FUNCTION io_groupexists(gid,name)RESULT(exist) 
!******************************************                             
!     checks if group called 'name' exists at position gid              
!                          D. Wortmann                                  
!******************************************                             
      USE hdf5 
      IMPLICIT NONE 
      !<--Arguments                                                     
      INTEGER(HID_T),INTENT(IN)  ::gid 
      CHARACTER,INTENT(IN)       ::name*(*) 
      LOGICAL                    ::exist 
      !>                                                                
      !<--Locals                                                        
      INTEGER(HID_t)::testid 
      INTEGER       ::hdferr 
      !>
      testid = 0
      hdferr = 0
      CALL h5eclear_f(hdferr) 
                                    !No automatic error checking!       
      CALL h5eset_auto_f(0, hdferr) 
      CALL h5gopen_f(gid,name,testid,hdferr) 
      exist=(hdferr.EQ.0) 
      CALL h5gclose_f(testid,hdferr) 
      CALL h5eclear_f(hdferr) 
                                    !Resume automatic error checking!   
      CALL h5eset_auto_f(1, hdferr) 
      END FUNCTION 
      !>                                                                
      !>                                                                
      END                                           
