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
  use m_juDFT_logging
  use m_juDFT_string,only:int2str
#ifdef CPP_MPI
  USE mpi
#endif
  IMPLICIT NONE
  PRIVATE
  CHARACTER(len=5),PARAMETER:: name="FLEUR"
  PUBLIC juDFT_error,juDFT_warn,juDFT_end,judft_file_readable, juDFT_BUG
CONTAINS

  SUBROUTINE judfT_file_readable(filename,warning)
    IMPLICIT NONE
    CHARACTER(len=*),INTENT(IN):: filename
    LOGICAL,INTENT(IN),OPTIONAL:: warning
    LOGICAL  :: l_exist

    INQUIRE(file=filename,exist=l_exist)
    IF (.not.l_exist) CALL judft_error("File not readable:"//filename,hint="You tried to read a file that is not present",warning=warning)
  END SUBROUTINE judfT_file_readable

  SUBROUTINE juDFT_BUG(message,calledby,hint,no,file,line)
   IMPLICIT NONE
   CHARACTER*(*),INTENT(IN)          :: message
   CHARACTER*(*),OPTIONAL,INTENT(IN) :: calledby,hint
   INTEGER,OPTIONAL,INTENT(IN)       :: no
   CHARACTER*(*),OPTIONAL,INTENT(IN) :: file
   INTEGER,INTENT(IN),OPTIONAL       :: line

   CALL juDFT_error(message,calledby,hint,no,bug = .TRUE.,file=file,line=line)
  END SUBROUTINE 

  SUBROUTINE juDFT_error(message,calledby,hint,no,warning,bug,file,line)

    USE iso_fortran_env ! for "output_unit"
    USE m_juDFT_internalParams
    USE m_judft_usage
    use m_juDFT_string
    USE m_judft_xmloutput
    IMPLICIT NONE
    CHARACTER*(*),INTENT(IN)          :: message
    CHARACTER*(*),OPTIONAL,INTENT(IN) :: calledby,hint
    INTEGER,OPTIONAL,INTENT(IN)       :: no
    LOGICAL,OPTIONAL,INTENT(IN)       :: warning,bug
    CHARACTER*(*),OPTIONAL,INTENT(IN) :: file
    INTEGER,INTENT(IN),OPTIONAL       :: line

    LOGICAL                       :: callstop,warn,first_pe
    LOGICAL                       :: l_mpi=.FALSE.
    INTEGER                       :: isize,irank,e,i
    CHARACTER(len=100),ALLOCATABLE::message_list(:)
    

   !For logging
    integer:: log_level=logmode_error
    type(t_log_message) :: log
#ifdef CPP_MPI
    LOGICAL :: first_parallel
    CALL MPI_INITIALIZED(l_mpi,e)
    IF (l_mpi) THEN
       CALL MPI_COMM_RANK(MPI_COMM_WORLD,irank,e)
       CALL MPI_COMM_SIZE(MPI_COMM_WORLD,isize,e)
    ELSE
       first_pe=.TRUE.
       isize=1
       irank=0
    ENDIF
#else
    first_pe=.TRUE.
    isize=1
    irank=0
#endif
    warn = .FALSE.
    IF (PRESENT(warning).and..not.present(bug)) warn = warning
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
    if (l_mpi) CALL collect_messages(message,message_list,first_pe)
#endif

    IF (first_pe) THEN
       IF (.NOT.warn) THEN
          IF (present(bug)) THEN
            WRITE(*,'(a)') "**************"//name//"-BUG*****************"
            log_level=logmode_bug
         ELSE  
            WRITE(*,'(a)') "**************"//name//"-Error*****************"
            log_level=logmode_error
          ENDIF   
       ELSE
         log_level=logmode_warning
         WRITE(*,'(a)') "************"//name//"-Warning*****************"
       ENDIF
       WRITE(*,"(3a)") "Error message: ",message
       call log%add("message",message)
       IF (PRESENT(calledby)) THEN
          WRITE(*,"(3a)") "Error occurred in subroutine: ",calledby
          call log%add("subroutine",calledby)
       ENDIF
       IF (PRESENT(hint)) THEN
          WRITE(*,"(3a)") "Hint: ",hint
          call log%add("hint",hint)
       ENDIF
       IF (PRESENT(no)) THEN
          WRITE(*,"(1a,i0)") "Error number: ",no
          call log%add("No",int2str(no))
       ENDIF
       IF (PRESENT(file)) THEN
          IF (PRESENT(line)) THEN
             WRITE(*,"(3a,i0)") "Source: ",file,":",line
             call log%add("Source",file//":"//int2str(line))
          ELSE
             WRITE(*,"(3a)") "Source: ",file
             call log%add("Source",file)
          ENDIF
       ENDIF
       IF (PRESENT(bug)) THEN
         write(0,*) "This is considered a BUG in "//name//". Please report it."
       ENDIF  
#ifdef CPP_MPI
       IF (l_mpi) THEN
          WRITE(*,'(a,i0,a,i0)') "Error from PE:",irank,"/",isize
          first_parallel=.TRUE.
          DO i=0,isize-1
             IF (i==irank) CYCLE
             IF (LEN_TRIM(message_list(i))>1)THEN
                IF (first_parallel) THEN
                   WRITE(*,'(2a)') "Other PEs with error messages:"
                   first_parallel=.FALSE.
                END IF
                WRITE(*,'(a,i4,2a)') "  ",i,"-",message_list(i)
             END IF
          END DO
       END IF
#endif
       WRITE(*,'(2a)') "*****************************************"

       IF (.NOT.warn) CALL juDFT_time_lastlocation(log)
       IF (callstop.and.warn) WRITE(*,'(a)')"Warnings not ignored. To make the warning nonfatal create a file 'JUDFT_WARN_ONLY' in the working directory or start FLEUR with the -warn_only command line option."
       IF (callstop) THEN
          CALL writetimes()
          CALL print_memory_info(output_unit,.TRUE.)
          IF (irank==0) THEN
             !Error on PE0 write info to out and out.xml
             WRITE(juDFT_outUnit,*) "***************ERROR***************"
             WRITE(juDFT_outUnit,*) message
             WRITE(juDFT_outUnit,*) "***************ERROR***************"
             !try closing the out file
             CLOSE(juDFT_outUnit,iostat=e)
             !Try closing the xml-out
             CALL endXMLOutput(errmsg=message)
          ENDIF
       END IF
    ELSE
       CALL priv_wait(2.0)
    ENDIF

    call log%report(log_level)

    IF (callstop) THEN
       CALL add_usage_data("Error",replace_text(message, new_line('A'), " "))
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
    USE iso_fortran_env ! for "output_unit"
    USE m_judft_xmlOutput
    USE m_judft_usage
    IMPLICIT NONE
    CHARACTER*(*), INTENT(IN)      :: message
    INTEGER, OPTIONAL, INTENT(IN)  :: irank
    LOGICAL, OPTIONAL, INTENT(IN)  :: l_endXML

    LOGICAL l_endXML_local, is_root
    LOGICAL :: l_mpi=.false.

    type(t_log_message)::log
#ifdef CPP_MPI
    INTEGER :: ierr
    CALL MPI_INITIALIZED(l_mpi,ierr)
#endif
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

    if(present(irank)) then
       is_root = (irank == 0)
    else
       is_root = .True.
    endif

    IF(is_root) THEN
       WRITE(*,*)
       WRITE(*,*) "********************************************************************"
       WRITE(*,*) "Run finished successfully"
       WRITE(*,*) "Stop message:"
       WRITE(*,*) "  ",message
       WRITE(*,*) "********************************************************************"
       WRITE(*,*) "If you publish work with contributions from FLEUR calculations,"
       WRITE(*,*) "please cite:"
       WRITE(*,*) ""
       WRITE(*,*) "  - The FLEUR project: https://www.flapw.de"
       WRITE(*,*) ""
       WRITE(*,*) "  - D. Wortmann et al., FLEUR, Zenodo, DOI: 10.5281/zenodo.7576163"
       WRITE(*,*) ""
       WRITE(*,*) "Please also consult on the website"
       WRITE(*,*) "  User Guide -> Reference -> References"
       WRITE(*,*) "for more information on relevant papers and example Bibtex entries."
       WRITE(*,*) "********************************************************************"
       FLUSH(output_unit)
      
    ENDIF

    !logging
    call log%add("Success",message)
    call log%report(logmode_status)

    CALL writetimes()
    CALL print_memory_info(output_unit,.true.)
    CALL send_usage_data()
#ifdef CPP_MPI
    IF (l_mpi) THEN
       IF(PRESENT(irank)) THEN
          CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)
          CALL MPI_ERRHANDLER_SET(MPI_COMM_WORLD,MPI_ERRORS_RETURN,ierr)
          CALL MPI_FINALIZE(ierr)
       ELSE
          CALL juDFT_STOP(0)
       END IF
    ENDIF
#endif
    if(is_root) then
       STOP
    else
       STOP
    endif
  END SUBROUTINE juDFT_END

  !this is a private subroutine that stops the calculations
  !different compilers might have to be added here
  SUBROUTINE juDFT_stop(errorCode)
#ifdef __INTEL_COMPILER
    USE ifcore
#endif
    INTEGER, OPTIONAL, INTENT(IN)  :: errorCode
    INTEGER :: error
    LOGICAL :: calltrace
    LOGICAL,ALLOCATABLE::a(:)
    LOGICAL :: l_mpi=.FALSE.

    !logging
    type(t_log_message):: log
#ifdef CPP_MPI
    INTEGER :: ierr
    CALL mpi_initialized(l_mpi,ierr)
#endif
    error = 1
    IF(PRESENT(errorCode)) THEN
       error = errorCode
    END IF

    !finalize logging
    call log%add("ExitCode",int2str(error))
    call log%report(logmode_status)
    
    call log_stop()
    
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
    ENDIF

#if defined(CPP_MPI)
    IF (l_mpi) THEN
       IF(error.EQ.0) THEN
          WRITE(*,*) ""
          WRITE(*,*) "Terminating all MPI processes."
          WRITE(*,*) "Note: This is a normal procedure."
          WRITE(*,*) "      Error messages in the following lines can be ignored."
          WRITE(*,*) ""
       END IF
       CALL MPI_ERRHANDLER_SET(MPI_COMM_WORLD,MPI_ERRORS_RETURN,ierr)
       CALL MPI_ABORT(MPI_COMM_WORLD,error,ierr)
    ENDIF
#endif
    IF (error.EQ.0) THEN
       STOP
    END IF
    STOP 1
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
