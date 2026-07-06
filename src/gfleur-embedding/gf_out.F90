!--------------------------------------------------------------------------------
! Copyright (c) 2026 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions 
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------
      MODULE m_gf_out 
      IMPLICIT NONE
!*****************************************************************      
! DESC:Module provides subroutines to handle the Output to unit 6 and 16
!                          Daniel Wortmann, Tue May  6 17:33:58 2003    
!*****************************************************************      
      PRIVATE 
                           !opens the out&inf file                      
      PUBLIC gf_out_create 
                           !closes both files again                     
      PUBLIC gf_out_close 
                               !writes a new section header             
      PUBLIC gf_out_newheader 
      CONTAINS 
                                                                        
      SUBROUTINE gf_out_create() 
      IMPLICIT NONE 
      !Open real files if not on Pe0                                    
      IF (isPe0()) THEN 
         OPEN(16,FILE="inf",FORM="FORMATTED",ACTION="WRITE",STATUS="REPLACE")
#ifndef  __TOS_BGQ__
      !Do not open out-file on BlueGene
              OPEN(6,FILE="out",FORM="FORMATTED",ACTION="WRITE",STATUS="REPLACE")
#endif
                                                                        
      ELSE 
      !Open Scratch files!                                              
         CLOSE(6,ERR=111) 
         CLOSE(16,ERR=111) 
  111    OPEN(16,FILE="inf.parallel",FORM="FORMATTED")
         OPEN(6,FILE="out.parallel",FORM="FORMATTED")
      ENDIF 
                                                                        
      END SUBROUTINE 
                                                                        
      SUBROUTINE gf_out_close() 
      IMPLICIT NONE 
      CLOSE(6) 
      CLOSE(16) 
      END SUBROUTINE 
                                                                        
                                                                        
      SUBROUTINE gf_out_newheader(name) 
      IMPLICIT NONE 
      CHARACTER*(*),INTENT(IN)::name 
      IF (.NOT.isPe0()) RETURN 
      WRITE(6,999) name 
      WRITE(16,999) name 
  999 FORMAT (/,/,'<***',/,10x,a,/,'***>') 
      END SUBROUTINE 
                                                                        
      FUNCTION isPe0() 
!*****************************************************************      
!  private function to check for pe=0                                   
!*****************************************************************      
      IMPLICIT NONE 
      LOGICAL::isPe0 
#ifdef CPP_MPI                                                          
      INCLUDE 'mpif.h' 
      INTEGER ierr,irank 
      CALL MPI_COMM_RANK (MPI_COMM_WORLD,irank,ierr) 
      isPe0=irank==0 
#else                                                                   
      isPe0=.TRUE. 
#endif                                                                  
      END FUNCTION 
      END                                           
