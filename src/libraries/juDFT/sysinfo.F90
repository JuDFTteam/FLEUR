!--------------------------------------------------------------------------------
! Copyright (c) 2026 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------


MODULE m_judft_sysinfo
  IMPLICIT NONE


CONTAINS

subroutine fp_error_check(onoff)
   LOGICAL,INTENT(IN) :: onoff
   Integer::errorStatus
#if defined(CPP_DEBUG)
      interface
         function startFPErrorDetection() bind(C, name="startFPErrorDetection")
            use iso_c_binding
            INTEGER(c_int) startFPErrorDetection
         end function startFPErrorDetection
      end interface

      interface
         function stopFPErrorDetection() bind(C, name="stopFPErrorDetection")
            use iso_c_binding
            INTEGER(c_int) stopFPErrorDetection
         end function stopFPErrorDetection
      end interface

   if (onoff) then
      errorStatus=startFPErrorDetection()
   else
      errorStatus=stopFPErrorDetection()
   endif
#endif
end subroutine fp_error_check

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
  
  !Print memory info to unit io. With maxmem=.TRUE. additionally report the
  !peak resident set size (high-water mark) of this process.
  !If l_par=.TRUE., the peak is in addition aggregated over all MPI ranks (max
  !and sum), so that the worst-case rank (drives OOM) and the total job
  !footprint are both visible.
  !Be carefull not to set l_par=.TRUE. if you are calling only from one MPI
  !rank: the reduction would then never be matched and would block forever.
  !That is why it defaults to .FALSE. - a caller that forgets it loses a line
  !of output, it does not hang.
  SUBROUTINE print_memory_info(io,maxmem,l_par)
#ifdef CPP_MPI
    USE mpi
#endif
    USE iso_c_binding
    IMPLICIT NONE
    INTEGER,INTENT(in)          :: io
    LOGICAL,INTENT(IN),OPTIONAL :: maxmem
    LOGICAL,INTENT(IN),OPTIONAL :: l_par
    INTEGER            :: err,irank,isize
    LOGICAL            :: l_mpi,l_par_local
    REAL(c_double)     :: loc_peak,max_peak,sum_peak
    REAL(c_double),PARAMETER :: gb=1024._c_double**3
    CHARACTER(len=40)  :: l1,l2

    INTERFACE
       FUNCTION peak_rss_bytes() bind(c)
          USE, INTRINSIC :: iso_c_binding
          REAL(kind=c_double) :: peak_rss_bytes
       END FUNCTION peak_rss_bytes
    END INTERFACE

    irank=0
    isize=1
    l_mpi=.FALSE.
    l_par_local=.FALSE.
    IF (PRESENT(l_par)) l_par_local=l_par
#ifdef CPP_MPI
    CALL mpi_initialized(l_mpi,err)
    IF (l_mpi) THEN
       CALL MPI_COMM_RANK(MPI_COMM_WORLD,irank,err)
       CALL MPI_COMM_SIZE(MPI_COMM_WORLD,isize,err)
    ENDIF
#endif

    loc_peak=peak_rss_bytes()
    max_peak=loc_peak
    sum_peak=loc_peak
#ifdef CPP_MPI
    IF (l_mpi.AND.l_par_local) THEN
       CALL MPI_REDUCE(loc_peak,max_peak,1,MPI_DOUBLE_PRECISION,MPI_MAX,0,MPI_COMM_WORLD,err)
       CALL MPI_REDUCE(loc_peak,sum_peak,1,MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_WORLD,err)
    ENDIF
#endif

    !With l_par all PEs are here and rank 0 reports for the whole job. Without
    !it every PE that gets here reports its own numbers: on the error path the
    !reporting PE is the first one with an error message, not necessarily rank 0.
    IF (irank==0.OR..NOT.l_par_local) THEN
       WRITE(io,"(a,i0,2a)") "Memory (rank ",irank,"): ",TRIM(memory_usage_string(maxmem))
       IF (PRESENT(maxmem)) THEN
          IF (maxmem.AND.l_par_local.AND.max_peak>0._c_double) THEN
             WRITE(l1,"(f12.3)") max_peak/gb
             WRITE(l2,"(f12.3)") sum_peak/gb
             WRITE(io,"(a,i0,a,a,a,a,a)") "Peak resident set over ",isize, &
                  " rank(s): max ",TRIM(ADJUSTL(l1))," GB, total ",TRIM(ADJUSTL(l2))," GB"
          END IF
       END IF
    END IF
  END SUBROUTINE print_memory_info
  
  !Return a human-readable string with the current resident set size of this
  !process and, with maxmem=.TRUE., the peak resident set size (high-water mark).
  !Both come from the portable C helpers in mem_usage.c (getrusage / proc / mach),
  !so a real value is produced on Linux, macOS and BSD rather than only on Linux.
  FUNCTION memory_usage_string(maxmem)
    USE iso_c_binding
    IMPLICIT NONE
    LOGICAL,INTENT(IN),OPTIONAL :: maxmem
    CHARACTER(len=100):: memory_usage_string
    REAL(c_double)    :: cur,pk
    REAL(c_double),PARAMETER :: gb=1024._c_double**3
    CHARACTER(len=40) :: line

    INTERFACE
       FUNCTION peak_rss_bytes() bind(c)
          USE, INTRINSIC :: iso_c_binding
          REAL(kind=c_double) :: peak_rss_bytes
       END FUNCTION peak_rss_bytes
       FUNCTION current_rss_bytes() bind(c)
          USE, INTRINSIC :: iso_c_binding
          REAL(kind=c_double) :: current_rss_bytes
       END FUNCTION current_rss_bytes
    END INTERFACE
#if CPP_GPU_CUDA
   interface
   function gpu_mem_usage() bind(c)
      use, intrinsic :: iso_c_binding
      real(kind=c_double) :: gpu_mem_usage
   end function gpu_mem_usage
   end interface
#endif

    cur=current_rss_bytes()
    pk =peak_rss_bytes()

    memory_usage_string=""
    IF (cur>0._c_double) THEN
       WRITE(line,"(f9.3,a)") cur/gb," GB"
       memory_usage_string=TRIM(ADJUSTL(line))
    ELSE IF (pk>0._c_double) THEN
       !current RSS not available on this platform - report peak as best estimate
       WRITE(line,"(f9.3,a)") pk/gb," GB"
       memory_usage_string=TRIM(ADJUSTL(line))
    ENDIF

    IF (PRESENT(maxmem)) THEN
       IF (maxmem.AND.cur>0._c_double.AND.pk>0._c_double) THEN
          WRITE(line,"(f9.3,a)") pk/gb," GB"
          memory_usage_string=TRIM(memory_usage_string)//" (peak "//TRIM(ADJUSTL(line))//")"
       ENDIF
    ENDIF

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
