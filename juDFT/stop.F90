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
      IMPLICIT NONE
      PRIVATE
      PUBLIC juDFT_error,juDFT_warn,juDFT_end
      CONTAINS


      SUBROUTINE juDFT_error(message,calledby,hint,no,warning,file,line)
      IMPLICIT NONE
      CHARACTER*(*),INTENT(IN)          :: message
      CHARACTER*(*),OPTIONAL,INTENT(IN) :: calledby,hint
      INTEGER,OPTIONAL,INTENT(IN)       :: no
      LOGICAL,OPTIONAL,INTENT(IN)       :: warning
      CHARACTER*(*),OPTIONAL,INTENT(IN) :: file
      INTEGER,INTENT(IN),OPTIONAL       :: line

      LOGICAL :: callstop,warn
      CHARACTER(LEN=4)::PE
      !store all output in variable for single call to write in MPI case
      CHARACTER(len=300)::text(10)
      INTEGER           ::linenr,n
#ifdef CPP_MPI
      include 'mpif.h'
      INTEGER :: irank,e
      CALL MPI_COMM_RANK(MPI_COMM_WORLD,irank,e)
      WRITE(PE,"(i4)") irank
#else
      PE="****"
#endif
      warn = .FALSE.
      IF (PRESENT(warning)) warn = warning
      IF (warn) THEN
         !check if we stop nevertheless
         INQUIRE(FILE ="JUDFT_WARN_ONLY",EXIST= callstop)
         callstop  = .NOT.callstop
      ELSE
         callstop = .TRUE.
      ENDIF

      IF (.NOT.warn) THEN
         WRITE(text(1),*) PE,"**************juDFT-Error*****************"
      ELSE
         WRITE(text(1),*) PE,"************juDFT-Warning*****************"
      ENDIF
      WRITE(text(2),"(3a)") PE,"Error message:",message
      linenr=3
      IF (PRESENT(calledby)) THEN
         WRITE(text(3),"(3a)") PE,"Error occurred in subroutine:",calledby
         linenr=4
      ENDIF
      IF (PRESENT(hint)) THEN
         WRITE(text(linenr),"(3a)") PE,"Hint:",hint
         linenr=linenr+1
      ENDIF
      IF (PRESENT(no)) THEN
         WRITE(text(linenr),"(2a,i0)") PE,"Error number:",no
         linenr=linenr+1
      ENDIF
      IF (present(file)) THEN
          if (present(line)) THEN
                write(text(linenr),"(4a,i0)") PE,"Source:",file,":",line
          ELSE
                write(text(linenr),"(3a)") PE,"Source:",file
          ENDIF
          linenr=linenr+1
      ENDIF
      WRITE(text(linenr),*) PE,"*****************************************"

      if (.not.warn) CALL juDFT_time_lastlocation(PE)

      IF (callstop) THEN
         IF (warn) THEN
               linenr=linenr+1
               WRITE(text(linenr),'(a)')"Warnings not ignored. Touch 'JUDFT_WARN_ONLY' to make the warning nonfatal"
         ENDIF
         write(0,"(10(a,/))") (trim(text(n)),n=1,linenr)
         CALL juDFT_STOP()
      ENDIF
      write(0,"(10(a,/))") (trim(text(n)),n=1,linenr)
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

      SUBROUTINE juDFT_END(message, irank)
      USE m_xmlOutput
      IMPLICIT NONE
#ifdef CPP_MPI
      INCLUDE 'mpif.h'
      INTEGER :: ierr
#endif
      CHARACTER*(*)        :: message
      INTEGER, INTENT(IN)  :: irank

      IF (irank.EQ.0) CALL endXMLOutput()

      WRITE(0,*) "*****************************************"
      WRITE(0,*) "Run finished successfully"
      WRITE(0,*) "Stop message:"
      WRITE(0,*) "  ",message
      WRITE(0,*) "*****************************************"
      CALL writetimes()
      CALL priv_memory_info()
#ifdef CPP_MPI
      CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)
      CALL MPI_FINALIZE(ierr)
#endif
      STOP 'OK'
      END SUBROUTINE juDFT_END

      !this is a private subroutine that stops the calculations
      !different compilers might have to be added here
      SUBROUTINE juDFT_stop()


#ifdef __INTEL_COMPILER
      USE ifcore
#endif
#ifdef CPP_MPI
      INCLUDE 'mpif.h'
      INTEGER :: ierr
#endif
      LOGICAL :: calltrace
      !try to print times
      !call writelocation()
      !CALL writetimes(.true.)
      INQUIRE(FILE="JUDFT_TRACE",EXIST=calltrace)
      IF (calltrace) THEN
#ifdef __INTEL_COMPILER
         CALL tracebackqq(USER_EXIT_CODE=-1) !return after traceback
#elif (defined(CPP_AIX)&&!defined(__PGI))
         CALL xl__trbk()
#endif
      ENDIF

#if defined(CPP_MPI)
      CALL MPI_ABORT(MPI_COMM_WORLD,1,ierr)
#endif
      STOP 'juDFT-STOPPED'
      END SUBROUTINE juDFT_stop

      SUBROUTINE priv_memory_info()
      IMPLICIT NONE

      CHARACTER(LEN=1024):: line
      INTEGER            :: err
      OPEN(99,FILE="/proc/self/status",ERR=999)
      DO
         READ(99,"(a)",ERR=999,END=999) line
         WRITE(6,*) trim(line)
      ENDDO
 999  CLOSE(99,IOSTAT=err)
      END SUBROUTINE priv_memory_info


      END MODULE m_juDFT_stop


