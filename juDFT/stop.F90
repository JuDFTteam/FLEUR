!--------------------------------------------------------------------------------
! Copyright (c) 2016 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------

MODULE m_juDFT_stop
  !-----------------------------------------------
  !    module to terminate Calculation, should be used instead
  !    of a simple STOP
  !
  !    error(message,calledby,hint,no,warning)
  !         message  : message string
  !         calledby : subroutine in which error occurs(optional)
  !         hint     : string with more information (optional)
  !         no       : error number (optional)
  !         warning  : logical indicating a warning message (optional)
  !
  !    warn(message,calledby,hint,no)
  !         shortcut for calling error with warning=.true.
  !
  !    juDFT_end(message)
  !         call this to terminate without error
  !
  !   IF the file "JUDFT_WARN_ONLY" is not present, warnings will lead to errors.
  !
  !   If the file "JUDFT_TRACE" is present, a stacktrace will be generated
  !   on some compilers
  !
  !
  !                 Daniel Wortmann (2010)
  !-----------------------------------------------
  USE m_judft_time
  USE m_judft_sysinfo
  USE m_judft_args
  IMPLICIT NONE
  PRIVATE
  PUBLIC juDFT_error,juDFT_warn,juDFT_end,judft_file_readable
CONTAINS

  SUBROUTINE judfT_file_readable(filename,warning)
    IMPLICIT NONE
    CHARACTER(len=*),INTENT(IN):: filename
    LOGICAL,INTENT(IN),OPTIONAL:: warning
    LOGICAL  :: l_exist

    INQUIRE(file=filename,exist=l_exist)
    IF (.not.l_exist) CALL judft_error("File not readable:"//filename,hint="You tried to read a file that is not present",warning=warning)
  END SUBROUTINE judfT_file_readable

  SUBROUTINE juDFT_error(message,calledby,hint,no,warning,file,line)
    USE m_judft_usage
    USE m_xmloutput
    IMPLICIT NONE
    CHARACTER*(*),INTENT(IN)          :: message
    CHARACTER*(*),OPTIONAL,INTENT(IN) :: calledby,hint
    INTEGER,OPTIONAL,INTENT(IN)       :: no
    LOGICAL,OPTIONAL,INTENT(IN)       :: warning
    CHARACTER*(*),OPTIONAL,INTENT(IN) :: file
    INTEGER,INTENT(IN),OPTIONAL       :: line

    LOGICAL                       :: callstop,warn,first_pe
    INTEGER                       :: isize,irank,e,i
    CHARACTER(len=100),ALLOCATABLE::message_list(:)
#ifdef CPP_MPI
    include 'mpif.h'
    LOGICAL :: first_parallel
    CALL MPI_COMM_RANK(MPI_COMM_WORLD,irank,e)
    CALL MPI_COMM_SIZE(MPI_COMM_WORLD,isize,e)
#else
    first_pe=.TRUE.
    isize=1
    irank=0
#endif
    warn = .FALSE.
    IF (PRESENT(warning)) warn = warning
    IF (warn) THEN
       !check if we stop nevertheless
       IF (judft_was_argument("-warn_only")) THEN
          callstop=.false.
       ELSE
          INQUIRE(FILE ="JUDFT_WARN_ONLY",EXIST= callstop)
          callstop  = .NOT.callstop
       ENDIF
    ELSE
       callstop = .TRUE.
    ENDIF
    
#ifdef CPP_MPI
    CALL collect_messages(message,message_list,first_pe)
#endif
    
    IF (first_pe) THEN
       IF (.NOT.warn) THEN
          WRITE(0,'(a)') "**************juDFT-Error*****************"
       ELSE
          WRITE(0,'(a)') "************juDFT-Warning*****************"
       ENDIF
       WRITE(0,"(3a)") "Error message:",message
       IF (PRESENT(calledby)) THEN
          WRITE(0,"(3a)") "Error occurred in subroutine:",calledby
       ENDIF
       IF (PRESENT(hint)) THEN
          WRITE(0,"(3a)") "Hint:",hint
       ENDIF
       IF (PRESENT(no)) THEN
          WRITE(0,"(1a,i0)") "Error number:",no
       ENDIF
       IF (PRESENT(file)) THEN
          IF (PRESENT(line)) THEN
             WRITE(0,"(3a,i0)") "Source:",file,":",line
          ELSE
             WRITE(0,"(3a)") "Source:",file
          ENDIF
       ENDIF
#ifdef CPP_MPI
       WRITE(0,'(a,i0,a,i0)') "Error from PE:",irank,"/",isize
       first_parallel=.TRUE.
       DO i=0,isize-1
          IF (i==irank) CYCLE
          IF (LEN_TRIM(message_list(i))>1)THEN
             IF (first_parallel) THEN
                WRITE(0,'(2a)') "Other PEs with error messages:"
                first_parallel=.FALSE.
             END IF
             WRITE(0,'(a,i4,2a)') "  ",i,"-",message_list(i)
          END IF
       END DO
#endif
       WRITE(0,'(2a)') "*****************************************"
       
       IF (.NOT.warn) CALL juDFT_time_lastlocation()
       IF (callstop.and.warn) WRITE(0,'(a)')"Warnings not ignored. Touch 'JUDFT_WARN_ONLY' to make the warning nonfatal"
       IF (callstop) THEN
          CALL writetimes()
          CALL print_memory_info(0,.TRUE.)
          IF (irank==0) THEN
             !Error on PE0 write info to out and out.xml
             WRITE(6,*) "***************ERROR***************"
             WRITE(6,*) message
             WRITE(6,*) "***************ERROR***************"
             CALL writeXMLElement('ERROR',(/"Message"/),(/message/))
             !try closing the out file
             CLOSE(6,iostat=e)
             !Try closing the xml-out
             CALL endXMLOutput()
          ENDIF
       END IF
    ELSE
       CALL priv_wait(2.0)
    ENDIF

    IF (callstop) THEN
       CALL add_usage_data("Error",message)
       !$OMP MASTER
       CALL send_usage_data()
       !$OMP END MASTER
       CALL juDFT_STOP()
    ENDIF
  END SUBROUTINE juDFT_error

  SUBROUTINE juDFT_warn(message,calledby,hint,no,file,line)
    IMPLICIT NONE
    CHARACTER*(*),INTENT(IN)          :: message
    CHARACTER*(*),OPTIONAL,INTENT(IN) :: calledby,hint
    INTEGER,OPTIONAL,INTENT(IN)       :: no
    CHARACTER*(*),OPTIONAL,INTENT(IN) :: file
    INTEGER,INTENT(IN),OPTIONAL       :: line

    CALL juDFT_error(message,calledby,hint,no,warning = .TRUE.,file=file,line=line)

  END SUBROUTINE juDFT_warn

  SUBROUTINE juDFT_END(message, irank, l_endXML)
    ! If irank is present every mpi process has to call this routine.
    ! Otherwise only a single mpi process is allowed to call the routine.
    USE m_xmlOutput
    USE m_judft_usage
    IMPLICIT NONE
#ifdef CPP_MPI
    INCLUDE 'mpif.h'
    INTEGER :: ierr
#endif
    CHARACTER*(*), INTENT(IN)      :: message
    INTEGER, OPTIONAL, INTENT(IN)  :: irank
    LOGICAL, OPTIONAL, INTENT(IN)  :: l_endXML

    LOGICAL l_endXML_local

    l_endXML_local = .TRUE.
    IF(PRESENT(l_endXML)) THEN
       l_endXML_local = l_endXML
    END IF

    IF(l_endXML_local) THEN
       IF(PRESENT(irank)) THEN
          IF (irank.EQ.0) CALL endXMLOutput()
       ELSE
          ! It is assumed that this is the only mpi process calling this routine.
          CALL endXMLOutput()
       END IF
    END IF
    IF (TRIM(message)=="") STOP ! simple stop if no end message is given
    WRITE(0,*)
    WRITE(0,*) "*****************************************"
    WRITE(0,*) "Run finished successfully"
    WRITE(0,*) "Stop message:"
    WRITE(0,*) "  ",message
    WRITE(0,*) "*****************************************"
    CALL writetimes()
    CALL print_memory_info(0,.true.)
    CALL send_usage_data()
#ifdef CPP_MPI
    IF(PRESENT(irank)) THEN
       write (*,*) "Going into post send barrier"
       CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)
       CALL MPI_FINALIZE(ierr)
    ELSE
       CALL juDFT_STOP(0)
    END IF
#endif
    STOP 'OK'
  END SUBROUTINE juDFT_END

  !this is a private subroutine that stops the calculations
  !different compilers might have to be added here
  SUBROUTINE juDFT_stop(errorCode)
#ifdef __INTEL_COMPILER
    USE ifcore
#endif
#ifdef CPP_MPI
    INCLUDE 'mpif.h'
#endif
    INTEGER, OPTIONAL, INTENT(IN)  :: errorCode
    INTEGER :: error
    LOGICAL :: calltrace
    LOGICAL, allocatable :: a(:)
#ifdef CPP_MPI
    INTEGER :: ierr
#endif
    error = 0

    IF(PRESENT(errorCode)) THEN
       error = errorCode
    END IF
    !Now try to generate a stack-trace if requested
    INQUIRE(FILE="JUDFT_TRACE",EXIST=calltrace)
    IF (judft_was_argument("-trace")) calltrace=.TRUE.
    IF (error.EQ.1) calltrace = .TRUE.
    IF (calltrace) THEN
#ifdef __INTEL_COMPILER
       CALL tracebackqq(USER_EXIT_CODE=-1) !return after traceback
#elif (defined(CPP_AIX)&&!defined(__PGI))
       CALL xl__trbk()
#endif
       ! cause an error, so that the compiler generates a stacktrace
       DEALLOCATE(a)
    ENDIF

#if defined(CPP_MPI)
    IF(error.EQ.0) THEN
       WRITE(0,*) ""
       WRITE(0,*) "Terminating all MPI processes."
       WRITE(0,*) "Note: This is a normal procedure."
       WRITE(0,*) "      Error messages in the following lines can be ignored."
       WRITE(0,*) ""
    END IF
    CALL MPI_ABORT(MPI_COMM_WORLD,error,ierr)
#endif
    STOP 'juDFT-STOPPED'
  END SUBROUTINE juDFT_stop


#ifdef CPP_MPI
  SUBROUTINE collect_messages(mymessage,message_list,first_pe)
    !This routine collects all error messages from all PE into an array
    !As not all PE might call this routine we use non-blocking communication
    !The integer first_pe is true if this pe is the one with the lowest rank among
    !all having error messages to report
    IMPLICIT NONE
    CHARACTER(len=*),INTENT(IN)                              :: mymessage
    CHARACTER(len=100),ASYNCHRONOUS,ALLOCATABLE,INTENT(OUT)  :: message_list(:)
    LOGICAL,INTENT(OUT)                                      :: first_pe
    INCLUDE 'mpif.h'
    INTEGER:: irank,isize,ierr,i
    LOGICAL:: completed
    INTEGER,ALLOCATABLE::ihandle(:)
    CHARACTER(len=100):: message
    REAL :: t1,t2
    
    message=mymessage
    
    CALL MPI_COMM_RANK(MPI_COMM_WORLD,irank,ierr)
    CALL MPI_COMM_SIZE(MPI_COMM_WORLD,isize,ierr)
    ALLOCATE(message_list(0:isize-1))
    ALLOCATE(ihandle(0:isize-1))
    message_list=""
    !Announce that I have a message to all PE
    DO i=0,isize-1
       CALL MPI_isend(message,100,MPI_CHARACTER,i,999,MPI_COMM_WORLD,ihandle(0),ierr)
    ENDDO
    !Collect all message
    DO i=0,isize-1
       CALL MPI_irecv(message_list(i),100,MPI_CHARACTER,i,999,MPI_COMM_WORLD,ihandle(i),ierr)
    ENDDO
    !Wait for 2 seconds
    CALL priv_wait(2.0)
    !Check if any PE with a lower rank also reports an error
    first_pe=.TRUE.
    DO i=0,isize-1
       CALL MPI_TEST(ihandle(i),completed,MPI_STATUS_IGNORE,ierr)
       IF (i<irank) first_pe=first_pe.AND..NOT.completed
    ENDDO
  END SUBROUTINE collect_messages
#endif

  SUBROUTINE priv_wait(sec)
    !Simple routine to wait for sec-seconds
    IMPLICIT NONE
    REAL::sec,t1,t2
    CALL cpu_TIME(t1)
    t2=t1
    DO WHILE(t2-t1<sec)
       CALL cpu_TIME(t2)
    ENDDO
  END SUBROUTINE priv_wait
END MODULE m_juDFT_stop


