!--------------------------------------------------------------------------------
! Copyright (c) 2016 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------
MODULE m_juDFT_time
!DEC$ NOOPTIMIZE
   !*****************************************************************
   !     DESC:Timer module for measuring the execution times of different
   !     parts of the code
   !     timestart and timestop should be
   !     called with suitable names for timers
   !     Daniel Wortmann, Fri Sep  6 11:53:08 2002
   !*****************************************************************
   USE m_judft_xmlOutput
   IMPLICIT NONE
   !     List of different timers
   PRIVATE
   INTEGER, PARAMETER        :: max_subtimer = 5 ! blocks of subtimers are allocated in this size
   REAL                      :: min_time = 0.02 ! minimal time to display in output (part of total)
   LOGICAL                   :: l_debug  !write out each start& stop of timer
   REAL                      :: debugtimestart = -1.0
   TYPE t_p
      TYPE(t_timer), POINTER:: p
   END TYPE t_p

   TYPE t_timer
      REAL                  :: starttime
      REAL                  :: time, mintime, maxtime
      INTEGER               :: no_calls
      CHARACTER(LEN=60)     :: name
      INTEGER               :: n_subtimers
      TYPE(t_p), ALLOCATABLE :: subtimer(:)
      TYPE(t_timer), POINTER :: parenttimer
   END TYPE t_timer

   TYPE(t_timer), POINTER, SAVE:: globaltimer => NULL()
   TYPE(t_timer), POINTER, SAVE:: current_timer => NULL()
   CHARACTER(LEN=256), SAVE   :: lastfile = ""
   INTEGER, SAVE   :: lastline = 0

   PUBLIC timestart, timestop, writetimes, writeTimesXML
   PUBLIC resetIterationDependentTimers, check_time_for_next_iteration
   PUBLIC juDFT_time_lastlocation !should be used for error handling only

CONTAINS

   SUBROUTINE priv_new_timer(name)
      IMPLICIT NONE
      CHARACTER(LEN=*), INTENT(IN)          ::name
      TYPE(t_timer), POINTER      ::t
      TYPE(t_p), ALLOCATABLE      :: tmp(:)

      ALLOCATE (t)
      t%starttime = cputime()
      t%time = 0.0
      t%mintime = 1E99
      t%maxtime = 0.0
      t%no_calls = 0
      t%name = name
      t%n_subtimers = 0
      t%parenttimer => null()
      ALLOCATE (t%subtimer(max_subtimer))

      IF (ASSOCIATED(current_timer)) THEN
         ALLOCATE (t%parenttimer)
         t%parenttimer => current_timer
         !now add timer to list of subtimers of parenttimer
         IF (current_timer%n_subtimers + 1 > SIZE(current_timer%subtimer)) THEN
            ALLOCATE (tmp(SIZE(current_timer%subtimer)))
            tmp = current_timer%subtimer
            DEALLOCATE (current_timer%subtimer)
            ALLOCATE (current_timer%subtimer(SIZE(tmp) + max_subtimer))
            current_timer%subtimer(:SIZE(tmp)) = tmp
            DEALLOCATE (tmp)
         ENDIF

         current_timer%n_subtimers = current_timer%n_subtimers + 1

         current_timer%subtimer(current_timer%n_subtimers)%p => t
      ELSE
         globaltimer => t
      ENDIF
      current_timer => t

   END SUBROUTINE priv_new_timer

   !<-- S: timestart(timer)

   SUBROUTINE timestart(ttimer, file, line)
      USE m_judft_args
      IMPLICIT NONE
      CHARACTER(LEN=*), INTENT(IN)          :: ttimer
      CHARACTER(LEN=*), INTENT(IN), OPTIONAL  :: file
      INTEGER, INTENT(IN), OPTIONAL           :: line

      INTEGER::n
      IF (PRESENT(file)) lastfile = file
      IF (PRESENT(line)) lastline = line
      IF (.NOT. ASSOCIATED(current_timer)) THEN
         CALL priv_new_timer("Total Run")
         l_debug = judft_was_Argument("-debugtime")
      ENDIF

      DO n = 1, current_timer%n_subtimers
         IF (TRIM(ttimer) == TRIM(current_timer%subtimer(n)%p%name)) THEN
            current_timer => current_timer%subtimer(n)%p
            IF (current_timer%starttime > 0) THEN
               WRITE (*, *) "Timer already running:", ttimer
               STOP "BUG:starttime"
            ENDIF
            current_timer%starttime = cputime()
            CALL priv_debug_output(" started ", current_timer%name)
            RETURN
         ENDIF
      ENDDO
      !new subtimer
      CALL priv_new_timer(ttimer)
      CALL priv_debug_output(" started ", current_timer%name)
   END SUBROUTINE timestart

   !>
   !<-- S:timestop(timer)

   SUBROUTINE timestop(ttimer)
      CHARACTER(LEN=*), INTENT(IN) :: ttimer

      REAL::time

      IF (.NOT. TRIM(ttimer) == TRIM(current_timer%name)) THEN
         WRITE (*, *) "Current timer:", current_timer%name, " could not stop:", ttimer
         STOP "BUG:timestop"
      ENDIF
      IF (current_timer%starttime < 0) THEN
         WRITE (*, *) "Timer not initialized:"//ttimer
         STOP "BUG:timestop"
      ENDIF
      time = cputime() - current_timer%starttime !This runtime
      current_timer%time = current_timer%time + time
      current_timer%no_calls = current_timer%no_calls + 1
      current_timer%mintime = MIN(current_timer%mintime, time)
      current_timer%maxtime = MAX(current_timer%maxtime, time)
      current_timer%starttime = -1

      CALL priv_debug_output(" stopped ", current_timer%name)

      current_timer => current_timer%parenttimer
   END SUBROUTINE timestop

   !>
   SUBROUTINE priv_debug_output(startstop, name)
      USE m_judft_sysinfo
      IMPLICIT NONE
      CHARACTER(LEN=*), INTENT(IN):: startstop, name
#ifdef CPP_MPI
      INTEGER::irank, ierr
      LOGICAL:: l_mpi
      INCLUDE 'mpif.h'
#endif
      IF (.NOT. l_debug) RETURN
      if (debugtimestart < 0) debugtimestart = cputime()
#ifdef CPP_MPI
      CALL MPI_INITIALIZED(l_mpi,ierr)
      IF (l_mpi) THEN
         CALL MPI_COMM_RANK(MPI_COMM_WORLD, irank, ierr)
         WRITE (*, "(i3,3a,f20.2,5x,a)") irank, startstop, name, " at:", cputime() - debugtimestart, memory_usage_string()
      ELSE
         WRITE (*, "(3a,f20.2,5x,a)") startstop, name, " at:", cputime() - debugtimestart, memory_usage_string()
      ENDIF
#else
      WRITE (*, "(3a,f20.2,5x,a)") startstop, name, " at:", cputime() - debugtimestart, memory_usage_string()
#endif
   END SUBROUTINE priv_debug_output

   RECURSIVE SUBROUTINE priv_writetimes_longest(timer, fid, timernames, timertimes)
      IMPLICIT NONE
      TYPE(t_timer), INTENT(IN)                 :: timer
      INTEGER, INTENT(IN), OPTIONAL              :: fid
      CHARACTER(LEN=60), INTENT(INOUT), OPTIONAL ::timernames(10)
      REAL, INTENT(INOUT), OPTIONAL              ::timertimes(10)

      REAL, ALLOCATABLE              ::times(:)
      CHARACTER(LEN=60), ALLOCATABLE ::names(:)
      REAL                          ::sum_time

      INTEGER :: n, i
      IF (.NOT. PRESENT(timernames)) THEN
         ALLOCATE (times(10), names(10))
         times = 0.0
         names = ""
         CALL priv_writetimes_longest(timer, timernames=names, timertimes=times)
         WRITE (fid, *)
         WRITE (fid, *) "-------------------------------------------------"
         WRITE (fid, *)
         WRITE (fid, *) "most relevant subroutines:"
         sum_time = 0.0
         DO n = 1, 10
            IF (MAXVAL(times) < 1E-4) EXIT
            i = MAXLOC(times, dim=1)
            WRITE (fid, "(a,T7,a,T46,a)") " ", TRIM(names(i)), timestring(times(i), timer%time, 2)
            sum_time = sum_time + times(i)
            times(i) = 0.0
         ENDDO
         WRITE (fid, "(t77,'Sum: ',f5.1,'%')") sum_time/timer%time*100.
         WRITE (fid, *)
         WRITE (fid, *) "-------------------------------------------------"
         WRITE (fid, *)
         RETURN
      ENDIF

      sum_time = timer%time
      DO n = 1, timer%n_subtimers
         CALL priv_writetimes_longest(timer%subtimer(n)%p, timernames=timernames, timertimes=timertimes)
         sum_time = sum_time - timer%subtimer(n)%p%time
      ENDDO
      IF (sum_time > MINVAL(timertimes)) THEN
         i = MINLOC(timertimes, dim=1)
         IF (ASSOCIATED(timer%parenttimer)) THEN
            WRITE (timernames(i), "(a18,'->',a40)") timer%parenttimer%name, timer%name
         ELSE
            timernames(i) = timer%name
         ENDIF
         timertimes(i) = sum_time
      ENDIF
   END SUBROUTINE priv_writetimes_longest

   RECURSIVE SUBROUTINE priv_writetimes(timer, level, fid, debug)
      IMPLICIT NONE
      TYPE(t_timer), INTENT(IN)    :: timer
      INTEGER, INTENT(IN)          :: level, fid
      LOGICAL, INTENT(IN), OPTIONAL ::debug

      INTEGER          :: n
      REAL             :: time
      CHARACTER(LEN=30):: timername, indentstring
      WRITE (timername, "(i0)") level - 1
      WRITE (indentstring, "(100a)") (" ", n=2, (MIN(level, 5))), TRIM(timername)
      IF (timer%starttime > 0) THEN
         time = timer%time + cputime() - timer%starttime
         timername = timer%name//" not term."
      ELSE
         time = timer%time
         timername = timer%name
      ENDIF
      IF (time < min_time*globaltimer%time) RETURN !do not print parts that take less than min_time of the totaltime
      IF (timer%no_calls > 1) THEN
         WRITE (fid, "(a,T7,a,T46,a,T85,a,i0,a,f9.3,a,f9.3,a)") &
            TRIM(indentstring), TRIM(timername), TRIM(timestring(time, globaltimer%time, level)), &
            "%  (", timer%no_calls, "calls:", timer%mintime, "sec -", timer%maxtime, "sec)"
      ELSE
         WRITE (fid, "(a,T7,a,T46,a)") TRIM(indentstring), TRIM(timername), timestring(time, globaltimer%time, level)
      END IF
      FLUSH(fid)
      IF (PRESENT(debug) .OR. timer%n_subtimers == 0) RETURN

      time = 0
      DO n = 1, timer%n_subtimers
         time = time + timer%subtimer(n)%p%time
      ENDDO
      WRITE (fid, "(a,a,T46,3a)") TRIM(indentstring), " measured in submodules:", TRIM(timestring(time, -.1, level))
      FLUSH(fid)
      DO n = 1, timer%n_subtimers
         CALL priv_writetimes(timer%subtimer(n)%p, level + 1, fid)
      ENDDO
   END SUBROUTINE priv_writetimes

   !<-- S:writetimes()
   
   RECURSIVE SUBROUTINE priv_genjson(timer, level, outstr, opt_idstr)
      use m_judft_string
      IMPLICIT NONE
      TYPE(t_timer), INTENT(IN)                    :: timer
      INTEGER, INTENT(IN)                          :: level
      CHARACTER(len=:), allocatable, INTENT(INOUT) :: outstr
      CHARACTER(len=:), allocatable, optional      :: opt_idstr

      CHARACTER(len=1), PARAMETER   :: nl=NEW_LINE("A")
      INTEGER, PARAMETER            :: indent_spaces=3

      INTEGER          :: n
      REAL             :: time
      CHARACTER(LEN=30):: timername 
      CHARACTER(LEN=:), allocatable :: idstr
   
      if(present(opt_idstr)) then
         idstr = opt_idstr
      else
         idstr = ""
      endif
     
      IF (timer%starttime > 0) THEN
         time = timer%time + cputime() - timer%starttime
         timername = timer%name//" not term."
      ELSE
         time = timer%time
         timername = timer%name
      ENDIF
      
      if(level > 1 ) outstr = outstr // nl 
      outstr = outstr // idstr // "{"
      idstr = idstr // repeat(" ", indent_spaces)

      outstr = outstr // nl // idstr // '"timername" : "' // trim(timername)         // '",'
      outstr = outstr // nl // idstr // '"totaltime" : '  // float2str(time)      
      if(level > 1) then
         outstr = outstr // ","
         outstr = outstr // nl // idstr // '"mintime"   : '  // float2str(timer%mintime)// ','
         outstr = outstr // nl // idstr // '"maxtime"   : '  // float2str(timer%maxtime)// ','
         outstr = outstr // nl // idstr // '"ncalls"    : '  // int2str(timer%no_calls) 
      endif

      time = 0
      DO n = 1, timer%n_subtimers
         time = time + timer%subtimer(n)%p%time
      ENDDO
      if(timer%n_subtimers > 0) then
         !add comma behind ncalls
         outstr = outstr // "," 
         outstr = outstr // nl // idstr // '"subtimers": '  // "[" 
         idstr = idstr // repeat(" ", indent_spaces)
         DO n = 1, timer%n_subtimers
            CALL priv_genjson(timer%subtimer(n)%p, level + 1, outstr, idstr)
            if(n /= timer%n_subtimers) outstr = outstr // ","
         ENDDO
         idstr  = idstr(:len(idstr)-indent_spaces) 
         outstr = outstr // nl // idstr // ']' 
      endif
      idstr  = idstr(:len(idstr)-indent_spaces) 
      outstr = outstr // nl // idstr // "}"
   END SUBROUTINE priv_genjson

   RECURSIVE SUBROUTINE writelocation(location)
      !writes the stack of current timers to std-out
      !usefull for debugging and error messages
      IMPLICIT NONE
      TYPE(t_timer), INTENT(IN), OPTIONAL::location

      IF (.NOT. PRESENT(location)) THEN
         IF (ASSOCIATED(current_timer)) CALL writelocation(current_timer)
      ELSE
         WRITE (0, *) "Timer:", location%name
         IF (ASSOCIATED(location%parenttimer)) CALL writelocation(location%parenttimer)
      ENDIF
   END SUBROUTINE writelocation

   ! writes all times to file
   SUBROUTINE writetimes(stdout)
     USE m_judft_usage
     USE m_judft_args
      IMPLICIT NONE
      LOGICAL, INTENT(IN), OPTIONAL::stdout
      INTEGER :: irank = 0
      CHARACTER(len=:), allocatable :: json_str
#ifdef CPP_MPI
      INCLUDE "mpif.h"
      INTEGER::err,isize
      LOGICAL:: l_mpi
      CHARACTER(len=30)::filename
      CALL mpi_initialized(l_mpi,err)
      if (l_mpi) CALL MPI_COMM_RANK(MPI_COMM_WORLD, irank, err)
#endif
      IF (.NOT. ASSOCIATED(globaltimer)) RETURN !write nothing if no timing recorded


      IF (irank == 0) THEN
         globaltimer%time = cputime() - globaltimer%starttime
         globaltimer%starttime = cputime()
         WRITE (6, "(//,'Total execution time: ',i0,'sec')") INT(globaltimer%time)
         CALL add_usage_data("Runtime", globaltimer%time)
         CALL priv_writetimes_longest(globaltimer, fid=6)

         WRITE (6, "('Total execution time: ',i0,'sec, minimal timing printed:',i0,'sec')") &
            INT(globaltimer%time), INT(min_time*globaltimer%time)

         CALL priv_writetimes(globaltimer, 1, 6)
#ifdef CPP_MPI
         IF (l_mpi) THEN
            CALL MPI_COMM_SIZE(MPI_COMM_WORLD, isize, err)
            WRITE (6, *) "Program used ", isize, " PE"
         ENDIF
#endif
      END IF
      IF (irank==0.OR.judft_was_argument("-all_times")) THEN
         json_str = ""
         CALL priv_genjson(globaltimer, 1, json_str)
         IF (irank==0) THEN
            OPEN(32, file="juDFT_times.json")
         ELSE
            WRITE(filename,"(a,i0,a)") "juDFT_times.",irank,".json"
            OPEN(32, file=trim(filename))
         END IF
         write (32,"(A)") json_str
         close(32)
      ENDIF
   END SUBROUTINE writetimes

   ! writes all times to out.xml file
   SUBROUTINE writeTimesXML()

      IMPLICIT NONE

      INTEGER                ::  irank = 0
      LOGICAL                :: l_out
      TYPE(t_timer), POINTER :: timer
#ifdef CPP_MPI
      INCLUDE "mpif.h"
      INTEGER::err, isize
      LOGICAL:: l_mpi
      CALL mpi_initialized(l_mpi,err)
      IF (l_mpi) CALL MPI_COMM_RANK(MPI_COMM_WORLD, irank, err)
#endif

      IF (irank .NE. 0) RETURN
      IF (.NOT. ASSOCIATED(globaltimer)) RETURN !write nothing if no timing recorded

      timer => NULL()
      CALL findIterationTimer(globaltimer, timer)
      IF (.NOT. ASSOCIATED(timer)) RETURN

      CALL openXMLElement('timing', (/'units'/), (/'sec'/))
      CALL privWriteTimesXML(timer, 1)
      CALL closeXMLElement('timing')

   END SUBROUTINE writeTimesXML

   RECURSIVE SUBROUTINE privWriteTimesXML(timer, level)

      IMPLICIT NONE

      TYPE(t_timer), INTENT(IN)    :: timer
      INTEGER, INTENT(IN)          :: level

      INTEGER            :: n, timerNameLength
      REAL               :: time
      CHARACTER(LEN=50)  :: timername
      CHARACTER(LEN=50)  :: attributes(2)

      IF (timer%starttime > 0) THEN
         time = timer%time + cputime() - timer%starttime
      ELSE
         time = timer%time
      END IF
      timername = ''
      DO n = 1, level - 2
         timername = timername(1:3*(n - 1))//'|  '
      END DO
      IF (level .GE. 2) THEN
         timername = timername(1:3*(level - 2))//'+--'
      END IF
      timername = TRIM(ADJUSTL(timername))//TRIM(ADJUSTL(timer%name))
      timerNameLength = LEN(TRIM(ADJUSTL(timername)))

      DO n = 1, timerNameLength
         IF (timername(n:n) .EQ. '&') THEN
            timername(n:n) = '+'
         END IF
      END DO

      WRITE (attributes(1), '(a)') TRIM(ADJUSTL(timername))
      WRITE (attributes(2), '(f12.3)') time
      IF (timer%n_subtimers .EQ. 0) THEN
         CALL writeXMLElementForm('timer', (/'name ', 'value'/), attributes, &
                                  RESHAPE((/14 + 20 - 3*level, 55 - timerNameLength, timerNameLength, 12/), (/2, 2/)))
      ELSE
         CALL openXMLElementForm('compositeTimer', (/'name ', 'value'/), attributes, &
                                 RESHAPE((/5 + 20 - 3*level, 55 - timerNameLength, timerNameLength, 12/), (/2, 2/)))
         DO n = 1, timer%n_subtimers
            CALL privWriteTimesXML(timer%subtimer(n)%p, level + 1)
         END DO
         CALL closeXMLElement('compositeTimer')
      END IF
   END SUBROUTINE privWriteTimesXML

   SUBROUTINE check_time_for_next_iteration(it, l_cont)
      USE m_judft_args
      IMPLICIT NONE
      INTEGER, INTENT(IN)     :: it
      LOGICAL, INTENT(INOUT)  :: l_cont
      CHARACTER(len=1000)::wtime_string
      INTEGER          :: wtime, time_used, time_per_iter
      INTEGER:: irank = 0
#ifdef CPP_MPI
      INCLUDE "mpif.h"
      INTEGER::err, isize
      LOGICAL:: l_mpi
      CALL mpi_initialized(l_mpi,err)
      if (l_mpi) CALL MPI_COMM_RANK(MPI_COMM_WORLD, irank, err)
#endif

      IF (.NOT. l_cont) RETURN !stop anyway
      IF (.NOT. ASSOCIATED(globaltimer)) RETURN !do nothing if no timing recorded
      IF (judft_was_argument("-wtime")) THEN
         wtime_string = judft_string_for_argument("-wtime")
         READ (wtime_string, *) wtime
         time_used = FLOOR((cputime() - globaltimer%starttime)/60.0) + 1
         time_per_iter = FLOOR((cputime() - globaltimer%starttime)/60.0/it) + 1
         IF (irank == 0) THEN
            WRITE (*, *) "Test for time of next iteration:"
            WRITE (*, *) "Time provided (min):", wtime
            WRITE (*, *) "Time used     (min):", time_used
            WRITE (*, *) "Time per iter (min):", time_per_iter
         ENDIF
         IF (time_used + time_per_iter > wtime) THEN
            l_cont = .FALSE.
            IF (irank == 0) WRITE (*, *) "No further iterations"
         ENDIF
      END IF
   END SUBROUTINE check_time_for_next_iteration

   SUBROUTINE resetIterationDependentTimers()

      IMPLICIT NONE

      INTEGER                ::  irank = 0
      LOGICAL                :: l_out
      TYPE(t_timer), POINTER :: timer, parenttimer
#ifdef CPP_MPI
      INCLUDE "mpif.h"
      INTEGER::err, isize
      LOGICAL:: l_mpi
      CALL mpi_initialized(l_mpi,err)
      if (l_mpi)  CALL MPI_COMM_RANK(MPI_COMM_WORLD, irank, err)
#endif
      !Check if not enough time for another iteration is left

      IF (irank .NE. 0) RETURN
      IF (.NOT. ASSOCIATED(globaltimer)) RETURN !write nothing if no timing recorded

      timer => NULL()
      CALL findIterationTimer(globaltimer, timer)
      IF (.NOT. ASSOCIATED(timer)) RETURN

      DO WHILE (TRIM(ADJUSTL(timer%name)) .NE. 'Iteration')
         IF (timer%n_subtimers .EQ. 0) RETURN
         timer => timer%subtimer(1)%p
      END DO

      parenttimer => timer%parenttimer
      CALL resetSubtimers(timer%parenttimer)
      ALLOCATE (parenttimer%subtimer(5))

   END SUBROUTINE resetIterationDependentTimers

   RECURSIVE SUBROUTINE findIterationTimer(parent, result)

      IMPLICIT NONE

      TYPE(t_timer), POINTER, INTENT(IN)    :: parent
      TYPE(t_timer), POINTER, INTENT(OUT)   :: result

      INTEGER numSubtimers, i

      IF (TRIM(ADJUSTL(parent%name)) .EQ. 'Iteration') THEN
         result => parent
      ELSE
         numSubtimers = parent%n_subtimers
         DO i = 1, numSubtimers
            CALL findIterationTimer(parent%subtimer(i)%p, result)
            IF (ASSOCIATED(result)) RETURN
         END DO
      END IF

   END SUBROUTINE findIterationTimer

   RECURSIVE SUBROUTINE resetSubtimers(timer)

      IMPLICIT NONE

      TYPE(t_timer), INTENT(INOUT)    :: timer

      INTEGER :: n

      DO n = 1, timer%n_subtimers
         CALL resetSubtimers(timer%subtimer(n)%p)
         DEALLOCATE (timer%subtimer(n)%p)
      END DO
      DEALLOCATE (timer%subtimer)
      timer%n_subtimers = 0
   END SUBROUTINE resetSubtimers

   !>
   !<-- private function timestring

   FUNCTION timestring(time, ttime, level)
      REAL, INTENT(IN)    :: time, ttime
      CHARACTER*(90)     :: timestring
      INTEGER, INTENT(IN) ::level
      !     .. Local Scalars ..
      REAL :: rest, seconds
      INTEGER :: ihours, iminutes

      rest = time
      ihours = INT(rest/3600.0)
      rest = rest - REAL(ihours)*3600
      iminutes = INT(rest/60.0)
      seconds = rest - REAL(iminutes)*60
      IF (ttime < 0) THEN
         WRITE (timestring, "(f9.2,'sec= ',i3,'h ',i2,'min ',i2,'sec')") time, ihours, iminutes, INT(seconds)
      ELSE
         WRITE (timestring, "(f9.2,'sec= ',i3,'h ',i2,'min ',i2,'sec ->',1x,f5.1,'%')") &
            time, ihours, iminutes, INT(seconds), time/ttime*100.0
      ENDIF
   END FUNCTION timestring

   !>
   !<-- F: cputime()
   !Private function to return the cpu_time in sec
   FUNCTION cputime()
!$    use omp_lib
#ifdef __INTEL_COMPILER
      USE ifport
#endif
      IMPLICIT NONE
#ifdef CPP_MPI
      INCLUDE 'mpif.h'
#endif
      REAL::cputime

      !TRY TO USE mpi OR openmp wall-clock functions
#ifdef CPP_MPI
      cputime = MPI_WTIME()
#elif __INTEL_COMPILER
      cputime = rtc()
#elif _OPENMP
      cputime = omp_get_wtime()
#else
      !use f95 intrinsic function
      CALL CPU_TIME(cputime)
#endif
   END FUNCTION cputime
   !>
   SUBROUTINE juDFT_time_lastlocation()
      IF (ASSOCIATED(current_timer)) THEN
         WRITE (0, *) "Last kown location:"
         WRITE (0, *) "Last timer:", current_timer%name
         IF (lastline > 0) THEN
            WRITE (0, *) "File:", TRIM(lastfile), ":", lastline
         ENDIF
         IF (ASSOCIATED(current_timer%parenttimer)) THEN
            WRITE (0, *) "Timerstack:"
            CALL writelocation(current_timer%parenttimer)
         ENDIF
         WRITE (0, *) "*****************************************"
      END IF
   END SUBROUTINE juDFT_time_lastlocation
END MODULE m_juDFT_time

