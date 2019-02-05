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
  SUBROUTINE print_memory_info(io,maxmem)
    IMPLICIT NONE
    INTEGER,INTENT(in)          :: io
    LOGICAL,INTENT(IN),OPTIONAL :: maxmem
    CHARACTER(LEN=1024):: line
    INTEGER            :: err,irank
#ifdef CPP_MPI
    INCLUDE 'mpif.h'
    CALL MPI_COMM_RANK(MPI_COMM_WORLD,irank,err)
#else
    irank=0
#endif
    WRITE(io,"(a,i0,a,a)") "Rank:",irank," used ",TRIM(memory_usage_string(maxmem))
  END SUBROUTINE print_memory_info
  
  !Read mstat to find out current memory usage
  FUNCTION memory_usage_string(maxmem)
    IMPLICIT NONE
    LOGICAL,INTENT(IN),OPTIONAL :: maxmem
    CHARACTER(len=100):: memory_usage_string
    INTEGER:: fid=543
    INTEGER*8:: idum
    LOGICAL:: firstcall=.TRUE.
    LOGICAL:: available=.FALSE.
    CHARACTER(len=40)::line
    
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

    IF (PRESENT(maxmem)) THEN
       IF (maxmem) THEN
          !Try to find maximal memory usage
          OPEN(544,FILE="/proc/self/status",status='old',action='read',iostat=idum)
          IF (idum==0) THEN
             DO
                READ(544,'(a)',iostat=idum) line
                IF (idum.NE.0) EXIT
                IF (INDEX(line,"VmPeak:")>0) THEN
                   memory_usage_string=TRIM(memory_usage_string)//"/"//TRIM(ADJUSTL(line(8:)))
                   CLOSE(544)
                   RETURN
                ENDIF
             ENDDO
             CLOSE(544)
          END IF
       END IF
    END IF
  END FUNCTION memory_usage_string
    
   
  SUBROUTINE checkstack()
    CHARACTER(LEN=10):: l1,l2,l3,l4
    INTEGER          :: err
    LOGICAL          :: unlimited
#ifdef CPP_MPI
    include 'mpif.h'
    INTEGER:: ierr,irank
#endif    
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
999 CLOSE(99,IOSTAT=err)
#ifdef CPP_MPI
    CALL MPI_COMM_RANK(MPI_COMM_WORLD,irank,ierr)
    IF (irank.NE.0) THEN
       IF (.NOT.unlimited) WRITE(*,*)"Warning, stacksize limited at PE:",irank
       RETURN
    ENDIF
#endif    
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
