!--------------------------------------------------------------------------------
! Copyright (c) 2026 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
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
  !Length of the error messages exchanged between the PEs
  INTEGER,PARAMETER         :: MESSAGE_LENGTH=100
  !Time (in seconds) spent waiting for the error messages of the other PEs.
  !An error that is actually reported terminates the run, so it is worth
  !waiting for the messages of the other PEs. An ignored warning is not: it
  !does not terminate anything and may be issued many times in one run, so
  !the wait would simply be added to the runtime over and over again.
  REAL,PARAMETER            :: WAIT_ERROR=2.0
  REAL,PARAMETER            :: WAIT_WARNING=0.2
#ifdef CPP_MPI
  !One-sided (RMA) window used by collect_messages to gather the error messages
  !of all PEs. Since juDFT_error might be called by an arbitrary subset of the
  !PEs only, no matching communication call on the other PEs can be assumed.
  !Therefore the window is created once, collectively, in
  !juDFT_init_errormessages and afterwards only used with passive target
  !synchronization.
  INTEGER,SAVE                                         :: errmsg_win
  LOGICAL,SAVE                                         :: l_errmsg_win=.FALSE.
  CHARACTER(len=MESSAGE_LENGTH),ALLOCATABLE,SAVE,TARGET:: errmsg_buffer(:)
#endif
  PUBLIC juDFT_error,juDFT_warn,juDFT_end,judft_file_readable, juDFT_BUG
  PUBLIC juDFT_init_errormessages,juDFT_free_errormessages
CONTAINS

  SUBROUTINE juDFT_init_errormessages()
    !Create the RMA window used by collect_messages to gather the error
    !messages of all PEs.
    !This routine is collective and has to be called by all PEs of
    !MPI_COMM_WORLD. It is called from juDFT_init directly after MPI_INIT.
    IMPLICIT NONE
#ifdef CPP_MPI
    INTEGER                        :: isize,ierr
    LOGICAL                        :: l_mpi
    INTEGER(KIND=MPI_ADDRESS_KIND) :: winsize

    IF (l_errmsg_win) RETURN !already initialized
    CALL MPI_INITIALIZED(l_mpi,ierr)
    IF (.NOT.l_mpi) RETURN

    CALL MPI_COMM_SIZE(MPI_COMM_WORLD,isize,ierr)
    ALLOCATE(errmsg_buffer(0:isize-1))
    errmsg_buffer=""
    winsize=INT(isize,MPI_ADDRESS_KIND)*INT(MESSAGE_LENGTH,MPI_ADDRESS_KIND)
    CALL MPI_WIN_CREATE(errmsg_buffer,winsize,1,MPI_INFO_NULL,MPI_COMM_WORLD,errmsg_win,ierr)
    IF (ierr.NE.MPI_SUCCESS) THEN
       !No window available. This is not fatal, collect_messages will then
       !simply not report the messages of the other PEs.
       DEALLOCATE(errmsg_buffer)
       RETURN
    ENDIF
    !Never abort while reporting an error, failing RMA calls are simply skipped
    CALL MPI_WIN_SET_ERRHANDLER(errmsg_win,MPI_ERRORS_RETURN,ierr)
    l_errmsg_win=.TRUE.
#endif
  END SUBROUTINE juDFT_init_errormessages

  SUBROUTINE juDFT_free_errormessages()
    !Free the RMA window again. This is collective as well and has to be
    !called by all PEs before MPI_FINALIZE.
    IMPLICIT NONE
#ifdef CPP_MPI
    INTEGER :: ierr

    IF (.NOT.l_errmsg_win) RETURN
    l_errmsg_win=.FALSE.
    CALL MPI_WIN_FREE(errmsg_win,ierr)
    DEALLOCATE(errmsg_buffer)
#endif
  END SUBROUTINE juDFT_free_errormessages

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
    CHARACTER(len=MESSAGE_LENGTH),ALLOCATABLE::message_list(:)
    

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
    if (l_mpi) CALL collect_messages(message,message_list,first_pe,callstop)
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
       !Give the reporting PE time to write its message first
       CALL priv_wait(MERGE(WAIT_ERROR,WAIT_WARNING,callstop))
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
    call log_stop()
    CALL writetimes()
    CALL print_memory_info(output_unit,.true.)
    CALL send_usage_data()
#ifdef CPP_MPI
    IF (l_mpi) THEN
       IF(PRESENT(irank)) THEN
          CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)
          !All PEs are here, hence the collective free of the RMA window is safe
          CALL juDFT_free_errormessages()
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
  SUBROUTINE collect_messages(mymessage,message_list,first_pe,l_callstop)
    !This routine collects the error messages of all PEs into an array.
    !As not all PEs might call this routine, neither collective communication
    !nor matching point-to-point calls can be used. Instead every PE writes its
    !message with one-sided communication (RMA) into the window of all PEs.
    !In contrast to the non-blocking send/recv used before, all communication is
    !guaranteed to be complete when MPI_WIN_UNLOCK returns. Hence no MPI
    !operation can access the local buffers after this routine has returned
    !(which corrupted the heap whenever juDFT_error did return, e.g. for
    !warnings with JUDFT_WARN_ONLY).
    !first_pe is true if this PE is the one with the lowest rank among all PEs
    !having an error message to report. l_callstop indicates that the caller
    !terminates the calculation; if it does not, the message is removed from the
    !windows again so that it is not reported by a later call, and the wait for
    !the other PEs is shortened.
    IMPLICIT NONE
    CHARACTER(len=*),INTENT(IN)                           :: mymessage
    CHARACTER(len=MESSAGE_LENGTH),ALLOCATABLE,INTENT(OUT) :: message_list(:)
    LOGICAL,INTENT(OUT)                                   :: first_pe
    LOGICAL,INTENT(IN)                                    :: l_callstop
    INTEGER                        :: irank,isize,ierr,i
    INTEGER(KIND=MPI_ADDRESS_KIND) :: disp
    LOGICAL                        :: l_flag
    REAL                           :: t1,t2,waittime
    CHARACTER(len=MESSAGE_LENGTH)  :: message

    CALL MPI_COMM_RANK(MPI_COMM_WORLD,irank,ierr)
    CALL MPI_COMM_SIZE(MPI_COMM_WORLD,isize,ierr)
    ALLOCATE(message_list(0:isize-1))
    message_list=""
    first_pe=.TRUE.
    !Without a window we can only report our own message
    IF (.NOT.l_errmsg_win) RETURN

    message=mymessage
    disp=INT(irank,MPI_ADDRESS_KIND)*INT(MESSAGE_LENGTH,MPI_ADDRESS_KIND)
    !Announce my message by writing it into my slot in the window of all PEs.
    !Only a single lock is held at a time, so this can not deadlock even if
    !several PEs report an error simultaneously.
    DO i=0,isize-1
       CALL MPI_WIN_LOCK(MPI_LOCK_EXCLUSIVE,i,0,errmsg_win,ierr)
       IF (ierr.NE.MPI_SUCCESS) CYCLE
       CALL MPI_PUT(message,MESSAGE_LENGTH,MPI_CHARACTER,i,disp,MESSAGE_LENGTH,&
                    MPI_CHARACTER,errmsg_win,ierr)
       !After the unlock the put is complete, i.e. 'message' can be reused
       CALL MPI_WIN_UNLOCK(i,errmsg_win,ierr)
    ENDDO

    !Wait to give the other PEs the chance to report as well. Only an error
    !that is really reported is worth the full wait, see the comment at
    !WAIT_ERROR/WAIT_WARNING above.
    !MPI_IPROBE is called in between to ensure progress of the one-sided
    !communication also with MPI libraries without asynchronous progress.
    waittime=MERGE(WAIT_ERROR,WAIT_WARNING,l_callstop)
    CALL cpu_TIME(t1)
    t2=t1
    DO WHILE(t2-t1<waittime)
       CALL MPI_IPROBE(MPI_ANY_SOURCE,MPI_ANY_TAG,MPI_COMM_WORLD,l_flag,MPI_STATUS_IGNORE,ierr)
       CALL cpu_TIME(t2)
    ENDDO

    !Now look at the messages collected in my own window
    CALL MPI_WIN_LOCK(MPI_LOCK_EXCLUSIVE,irank,0,errmsg_win,ierr)
    IF (ierr==MPI_SUCCESS) THEN
       message_list=errmsg_buffer
       !Consume the messages so that a later call does not report them again
       errmsg_buffer=""
       CALL MPI_WIN_UNLOCK(irank,errmsg_win,ierr)
    ENDIF

    !Check if any PE with a lower rank also reports an error
    DO i=0,irank-1
       IF (LEN_TRIM(message_list(i))>0) first_pe=.FALSE.
    ENDDO

    !If we do not stop here (a warning in warn-only mode) the message has to be
    !removed from the windows again. Otherwise it would show up as a stale
    !message of "another PE" in a later error report.
    IF (.NOT.l_callstop) THEN
       message=""
       DO i=0,isize-1
          CALL MPI_WIN_LOCK(MPI_LOCK_EXCLUSIVE,i,0,errmsg_win,ierr)
          IF (ierr.NE.MPI_SUCCESS) CYCLE
          CALL MPI_PUT(message,MESSAGE_LENGTH,MPI_CHARACTER,i,disp,MESSAGE_LENGTH,&
                       MPI_CHARACTER,errmsg_win,ierr)
          CALL MPI_WIN_UNLOCK(i,errmsg_win,ierr)
       ENDDO
    ENDIF
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
