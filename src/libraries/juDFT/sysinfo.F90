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
  
integer function num_openmp()
!$  use omp_lib
num_openmp=0
!$ num_openmp=omp_get_max_threads()
end function  

integer function num_gpu()
#ifdef _OPENACC
   use openacc
#endif
   num_gpu=0
#ifdef _OPENACC
   num_gpu=acc_get_num_devices(acc_device_nvidia)
#endif
end function

function uname()
   CHARACTER(:), allocatable::uname
   integer::iostat
   character(len=200)::tmp
   CALL EXECUTE_COMMAND_LINE("uname -a >uname.log")
   OPEN(1234,FILE="uname.log",status='old',action='read',iostat=iostat)
   if (iostat==0) THEN 
      read(1234,'(a)') tmp
      close(1234,status='delete')
      uname=trim(tmp)
   else 
      uname="unkown"
   endif 
end function   
  
  !Dump status file into out file (should be done only at end of run)
  SUBROUTINE print_memory_info(io,maxmem)
#ifdef CPP_MPI
    USE mpi
#endif
    IMPLICIT NONE
    INTEGER,INTENT(in)          :: io
    LOGICAL,INTENT(IN),OPTIONAL :: maxmem
    CHARACTER(LEN=1024):: line
    INTEGER            :: err,irank
    LOGICAL :: l_mpi=.FALSE.
#ifdef CPP_MPI
    CALL mpi_initialized(l_mpi,err)
    IF (l_mpi) THEN
       CALL MPI_COMM_RANK(MPI_COMM_WORLD,irank,err)
    ELSE
       irank=0
    ENDIF
#else
    irank=0
#endif
    if(irank == 0) WRITE(io,"(a,i0,a,a)") "Rank:",irank," used ",TRIM(memory_usage_string(maxmem))
  END SUBROUTINE print_memory_info
  
  !Read mstat to find out current memory usage
  FUNCTION memory_usage_string(maxmem)
   use iso_c_binding
    IMPLICIT NONE
    LOGICAL,INTENT(IN),OPTIONAL :: maxmem
    CHARACTER(len=100):: memory_usage_string
    INTEGER:: fid=543
    INTEGER*8:: idum,rss
    LOGICAL:: firstcall=.TRUE.
    LOGICAL:: available=.FALSE.
    CHARACTER(len=40)::line
   
#if CPP_GPU_CUDA
   interface
   function gpu_mem_usage() bind(c)
      use, intrinsic :: iso_c_binding
      real(kind=c_double) :: gpu_mem_usage
   end function gpu_mem_usage
   end interface    
#endif

    IF (firstcall) THEN
       firstcall=.FALSE.
       OPEN(fid,FILE="/proc/self/statm",status='old',action='read',iostat=idum)
       available=(idum==0)
    ENDIF
    IF (available) THEN
       REWIND(fid)
       READ(fid,*) idum,rss
       WRITE(memory_usage_string,"(2f8.3,a)") (CPP_PAGESIZE/(1024.*1024.*1024.))*idum,(CPP_PAGESIZE/(1024.*1024.*1024.))*rss," GB"
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

#ifdef CPP_GPU_CUDA
   write(line,"(f10.3)") gpu_mem_usage()
   memory_usage_string=TRIM(memory_usage_string)//" GPU: "//trim(line)//"GB"
#endif      
  END FUNCTION memory_usage_string
    
   
  SUBROUTINE checkstack()
#ifdef CPP_MPI
    USE mpi
#endif
    CHARACTER(LEN=10):: l1,l2,l3,l4
    INTEGER          :: err
    LOGICAL          :: unlimited,l_mpi
#ifdef CPP_MPI
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
    CALL mpi_initialized(l_mpi,ierr)
    irank=0
    if (l_mpi) CALL MPI_COMM_RANK(MPI_COMM_WORLD,irank,ierr)
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
