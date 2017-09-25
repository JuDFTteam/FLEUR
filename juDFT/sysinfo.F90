!--------------------------------------------------------------------------------
! Copyright (c) 2016 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------


!To generate Memory usage info on LINUX we need to know the pagesize
!a default of 4096 bytes is assumed
#ifndef CPP_PAGESIZE
#define CPP_PAGESIZE 4096
#endif
 MODULE m_judft_sysinfo
      IMPLICIT NONE
    CONTAINS
  


    !Dump status file into out file (should be done only at end of run)
    SUBROUTINE print_memory_info()
      IMPLICIT NONE

      CHARACTER(LEN=1024):: line
      INTEGER            :: err
      OPEN(99,FILE="/proc/self/status",ERR=999)
      DO
         READ(99,"(a)",ERR=999,END=999) line
         WRITE(6,*) trim(line)
      ENDDO
 999  CLOSE(99,IOSTAT=err)
    END SUBROUTINE print_memory_info

    !Read mstat to find out current memory usage
    FUNCTION memory_usage_string()
      IMPLICIT NONE
      CHARACTER(len=10):: memory_usage_string
      INTEGER:: fid=543
      INTEGER*8:: idum
      LOGICAL:: firstcall=.TRUE.
      LOGICAL:: available=.FALSE.
      
      IF (firstcall) THEN
         firstcall=.FALSE.
         OPEN(fid,FILE="/proc/self/statm",status='old',action='read',iostat=idum)
         available=(idum==0)
      ENDIF
      IF (available) THEN
         REWIND(fid)
         READ(fid,*) idum
         WRITE(memory_usage_string,"(f8.3,a)") (CPP_PAGESIZE/(1024.*1024.*1024.))*idum,"GB"
      ELSE
         memory_usage_string=""
      ENDIF
    END FUNCTION memory_usage_string
    
   
      SUBROUTINE checkstack()
        CHARACTER(LEN=10):: l1,l2,l3,l4
        INTEGER          :: err
        LOGICAL          :: unlimited
        unlimited=.TRUE.  !set to true by default. 
        !If /proc/self/limits does not exist
        !or parsing fails no warning is issued
        OPEN(99,FILE="/proc/self/limits",ERR=999)
        DO
           READ(99,*,ERR=999,END=999) l1,l2,l3,l4
           IF (ALL((/INDEX(l1,"Max"),INDEX(l2,"stack"),INDEX(l3,"size")/)==1)) THEN
              unlimited=INDEX(l4,"unlim")==1
           ENDIF
        ENDDO
999     CLOSE(99,IOSTAT=err)
        IF (.NOT.unlimited) THEN
           WRITE(*,*) "*********** WARNING! ************"
           WRITE(*,*) "Your stacksize seems to be small"
           WRITE(*,*) "FLEUR might crash without further"
       WRITE(*,*) "notice. Try 'ulimit -s unlimited'"
       WRITE(*,*) "*********** WARNING! ************"
    ENDIF
  END SUBROUTINE checkstack
END MODULE m_judft_sysinfo

!program test
!  use m_judft_sysinfo
!  call checkstack()

!end program test
