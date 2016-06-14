      MODULE m_juDFT_time
!*****************************************************************
!     DESC:Timer module for measuring the execution times of different
!     parts of the code
!     timestart and timestop should be
!     called with suitable names for timers
!     Daniel Wortmann, Fri Sep  6 11:53:08 2002
!*****************************************************************
      USE m_xmlOutput
      IMPLICIT NONE
!     List of different timers
      PRIVATE
      INTEGER,PARAMETER         :: max_subtimer=5 ! blocks of subtimers are allocated in this size
      REAL                      :: min_time=0.02 ! minimal time to display in output (part of total)
      TYPE t_p
         TYPE(t_timer),POINTER:: p
      END TYPE

      TYPE t_timer
        REAL                  :: starttime
        REAL                  :: time
        CHARACTER(LEN=60)     :: name
        INTEGER               :: n_subtimers
        TYPE(t_p),ALLOCATABLE :: subtimer(:)
        TYPE(t_timer),POINTER :: parenttimer
      END TYPE

      TYPE(t_timer),POINTER,SAVE:: globaltimer=>NULL()
      TYPE(t_timer),POINTER,SAVE:: current_timer=>NULL()
      CHARACTER(LEN=256),SAVE   :: lastfile=""
      INTEGER           ,SAVE   :: lastline=0

      PUBLIC timestart,timestop,writetimes,writelocation,writeTimesXML
      PUBLIC resetIterationDependentTimers
      PUBLIC juDFT_time_lastlocation !should not be used

      CONTAINS



      SUBROUTINE priv_new_timer(name)
      IMPLICIT NONE
      CHARACTER(LEN=*),INTENT(IN)          ::name



      TYPE(t_timer),POINTER      ::t
      TYPE(t_p),ALLOCATABLE      :: tmp(:)


      ALLOCATE(t)
      t%starttime=cputime()
      t%time=0.0
      t%name=name
      t%n_subtimers=0
      ALLOCATE(t%subtimer(max_subtimer))

      IF (associated(current_timer)) THEN
         ALLOCATE(t%parenttimer)
         t%parenttimer=>current_timer
          !now add timer to list of subtimers of parenttimer
         IF (current_timer%n_subtimers+1>size(current_timer%subtimer)) THEN
                ALLOCATE(tmp(size(current_timer%subtimer)))
                tmp=current_timer%subtimer
                DEALLOCATE(current_timer%subtimer)
                ALLOCATE(current_timer%subtimer(size(tmp)+max_subtimer))
                current_timer%subtimer(:size(tmp))=tmp
                DEALLOCATE(tmp)
         ENDIF

         current_timer%n_subtimers=current_timer%n_subtimers+1

         current_timer%subtimer(current_timer%n_subtimers)%p => t
      ELSE
         globaltimer=>t
      ENDIF
      current_timer=>t

      END SUBROUTINE priv_new_timer

      !<-- S: timestart(timer)

      SUBROUTINE timestart(ttimer,file,line)
      IMPLICIT NONE
      CHARACTER(LEN =*),INTENT(IN)          :: ttimer
      CHARACTER(LEN=*),INTENT(IN),OPTIONAL  :: file
      INTEGER,INTENT(IN),OPTIONAL           :: line

      INTEGER::n
#ifdef CPP_MPI
      INTEGER::irank,ierr
      include 'mpif.h'
#endif
      IF (present(file)) lastfile=file
      IF (present(line)) lastline=line
      IF (.NOT.associated(current_timer)) THEN
        CALL priv_new_timer("Total Run")
      ENDIF

      DO n=1,current_timer%n_subtimers
        IF (trim(ttimer)==trim(current_timer%subtimer(n)%p%name)) THEN
           current_timer=>current_timer%subtimer(n)%p
           IF (current_timer%starttime>0) THEN
               WRITE(*,*) "Timer already running:",ttimer
               STOP "BUG:starttime"
           ENDIF
           current_timer%starttime=cputime()
#ifdef CPP_DEBUG
#ifdef CPP_MPI
           CALL MPI_COMM_RANK(MPI_COMM_WORLD,irank,ierr)
#endif
           WRITE(*,"(i3,3a,i10)") irank," started: ",current_timer%name," at:",cputime()
#endif
           RETURN
        ENDIF
      ENDDO
      !new subtimer
      CALL priv_new_timer(ttimer)
#ifdef CPP_DEBUG
#ifdef CPP_MPI
           CALL MPI_COMM_RANK(MPI_COMM_WORLD,irank,ierr)
#endif
           WRITE(*,*) irank," started: ",current_timer%name
#endif
      END SUBROUTINE timestart

      !>
      !<-- S:timestop(timer)

      SUBROUTINE timestop(ttimer)
      CHARACTER(LEN =*),INTENT(IN) :: ttimer

      IF (.NOT.trim(ttimer)==trim(current_timer%name)) THEN
          WRITE(*,*)"Current timer:",current_timer%name," could not stop:",ttimer
          STOP "BUG:timestop"
      ENDIF
      IF (current_timer%starttime<0) THEN
             WRITE(*,*)   "Timer not initialized:"//ttimer
             STOP "BUG:timestop"
      ENDIF
      current_timer%time=current_timer%time+cputime()-current_timer%starttime
      current_timer%starttime=-1

#ifdef CPP_DEBUG
      CALL priv_writetimes(current_timer,5,6)
#endif
      current_timer=>current_timer%parenttimer
      END SUBROUTINE timestop

      !>

      RECURSIVE SUBROUTINE priv_writetimes_longest(timer,fid,timernames,timertimes)
      IMPLICIT NONE
      TYPE(t_timer),INTENT(IN) :: timer
      INTEGER,INTENT(IN),OPTIONAL              :: fid
      CHARACTER(LEN=60),INTENT(INOUT),OPTIONAL ::timernames(10)
      REAL,INTENT(INOUT),OPTIONAL              ::timertimes(10)

      REAL,ALLOCATABLE              ::times(:)
      CHARACTER(LEN=60),ALLOCATABLE ::names(:)
      real                          ::sum_time

      integer :: n,i
      IF (.NOT.present(timernames)) THEN
            ALLOCATE(times(10),names(10))
            times=0.0
            names=""
            CALL priv_writetimes_longest(timer,timernames=names,timertimes=times)
            WRITE(fid,*)
            WRITE(fid,*) "-------------------------------------------------"
            WRITE(fid,*)
            WRITE(fid,*) "most relevant subroutines:"
            sum_time=0.0
            DO n=1,10
                IF (maxval(times)<1E-4) exit
                i=maxloc(times,dim=1)
                WRITE(fid,"(a,T7,a,T46,a)") " ",trim(names(i)),timestring(times(i),timer%time,2)
                sum_time=sum_time+times(i)
                times(i)=0.0
            ENDDO
            WRITE(fid,"(t77,'Sum:  ',f4.1,'%')") sum_time/timer%time*100.
            WRITE(fid,*)
            WRITE(fid,*) "-------------------------------------------------"
            WRITE(fid,*)
            RETURN
      ENDIF

      sum_time=timer%time
      DO n=1,timer%n_subtimers
         CALL priv_writetimes_longest(timer%subtimer(n)%p,timernames=timernames,timertimes=timertimes)
         sum_time=sum_time-timer%subtimer(n)%p%time
      ENDDO
      IF (sum_time>minval(timertimes)) THEN
            i=minloc(timertimes,dim=1)
            if (associated(timer%parenttimer)) THEN
                  write(timernames(i),"(a18,'->',a40)") timer%parenttimer%name,timer%name
            else
                 timernames(i)=timer%name
            endif
            timertimes(i)=sum_time
      ENDIF
      END SUBROUTINE

      RECURSIVE SUBROUTINE priv_writetimes(timer,level,fid,debug)
      IMPLICIT NONE
      TYPE(t_timer),INTENT(IN)    :: timer
      INTEGER,INTENT(IN)          :: level,fid
      LOGICAL,INTENT(IN),OPTIONAL ::debug

      INTEGER          :: n
      REAL             :: time
      CHARACTER(LEN=30):: timername,indentstring
      WRITE(timername,"(i0)") level-1
      WRITE(indentstring,"(100a)") (" ",n=2,(min(level,5))),trim(timername)
      IF (timer%starttime>0) THEN
         time=timer%time+cputime()-timer%starttime
         timername=timer%name//" not term."
      ELSE
         time=timer%time
         timername=timer%name
      ENDIF
      IF (time<min_time*globaltimer%time) RETURN !do not print parts that take less than min_time
      WRITE(fid,"(a,T7,a,T46,a)") trim(indentstring),trim(timername),timestring(time,globaltimer%time,level)
      flush(fid)
      IF (present(debug).OR.timer%n_subtimers==0) RETURN

      time=0
      DO n=1,timer%n_subtimers
        time=time+timer%subtimer(n)%p%time
      ENDDO
      WRITE(fid,"(a,a,T46,a)") trim(indentstring)," measured in submodules:",timestring(time,-.1,level)
      flush(fid)
      DO n=1,timer%n_subtimers
         CALL priv_writetimes(timer%subtimer(n)%p,level+1,fid)
      ENDDO
      END SUBROUTINE priv_writetimes

      !<-- S:writetimes()

      RECURSIVE SUBROUTINE writelocation(location)
      !writes the stack of current timers to std-out
      !usefull for debugging and error messages
      IMPLICIT NONE
      TYPE(t_timer),INTENT(IN),OPTIONAL::location

      IF (.NOT.present(location)) THEN
         IF (associated(current_timer)) CALL writelocation(current_timer)
      ELSE
        WRITE(*,*) "Timer:",location%name
        IF (associated(location%parenttimer)) CALL writelocation(location%parenttimer)
      ENDIF
      END SUBROUTINE writelocation


      ! writes all times to file
      SUBROUTINE writetimes(stdout)
      IMPLICIT NONE
      LOGICAL,INTENT(IN),OPTIONAL::stdout
      INTEGER :: fn,irank=0
      LOGICAL :: l_out
#ifdef CPP_MPI
      include "mpif.h"
      INTEGER::err,isize

      CALL MPI_COMM_RANK(MPI_COMM_WORLD,irank,err)
#endif
      IF (.NOT.associated(globaltimer)) RETURN !write nothing if no timing recorded
      l_out=.FALSE.
      IF (present(stdout)) l_out=stdout
      IF (l_out) THEN
         fn=6
      ELSE
         IF (irank>0) RETURN
         fn=2
         OPEN(2,FILE ='juDFT_times',STATUS="replace")
      ENDIF
      !Stop globaltimer if still running
      IF (globaltimer%starttime>-1) THEN
        globaltimer%time=cputime()-globaltimer%starttime
        globaltimer%starttime=-1
      ENDIF
      write(fn,"('Total execution time: ',i0,'sec')") int(globaltimer%time)
      CALL priv_writetimes_longest(globaltimer,fid=fn)

      WRITE(fn,"('Total execution time: ',i0,'sec, minimal timing printed:',i0,'sec')") &
                  int(globaltimer%time),int(min_time*globaltimer%time)
#ifdef CPP_MPI
      CALL MPI_COMM_SIZE(MPI_COMM_WORLD,isize,err)
      WRITE(fn,*) "Program used ",isize," PE"
#endif
      CALL priv_writetimes(globaltimer,1,fn)



      WRITE(fn,*)
      WRITE(fn,*) "-------------------------------------------------"
      WRITE(fn,*)
      WRITE(fn,*) "Total timings:"
      min_time=0.0
      CALL priv_writetimes(globaltimer,1,fn)
      flush(fn)
      IF (.NOT.l_out) CLOSE(2)

      END SUBROUTINE writetimes

      ! writes all times to out.xml file
      SUBROUTINE writeTimesXML()

         IMPLICIT NONE

         INTEGER                :: fn,irank=0
         LOGICAL                :: l_out
         TYPE(t_timer), POINTER :: timer
#ifdef CPP_MPI
         include "mpif.h"
         INTEGER::err,isize

         CALL MPI_COMM_RANK(MPI_COMM_WORLD,irank,err)
#endif
         IF (irank.NE.0) RETURN
         IF (.NOT.associated(globaltimer)) RETURN !write nothing if no timing recorded

         timer => globaltimer
         DO WHILE (TRIM(ADJUSTL(timer%name)).NE.'Iteration')
            IF(timer%n_subtimers.EQ.0) RETURN
            timer => timer%subtimer(1)%p
         END DO

         CALL openXMLElement('timing',(/'units'/),(/'sec'/))
         CALL privWriteTimesXML(timer,1)
         CALL closeXMLElement('timing')

      END SUBROUTINE writeTimesXML

      RECURSIVE SUBROUTINE privWriteTimesXML(timer,level)

         IMPLICIT NONE

         TYPE(t_timer),INTENT(IN)    :: timer
         INTEGER,INTENT(IN)          :: level

         INTEGER            :: n, timerNameLength
         REAL               :: time
         CHARACTER(LEN=50)  :: timername
         CHARACTER(LEN=50)  :: attributes(2)

         IF (timer%starttime>0) THEN
            time=timer%time+cputime()-timer%starttime
         ELSE
            time=timer%time
         END IF
         timername=''
         DO n = 1, level-2
            timername = timername(1:3*(n-1))//'|  '
         END DO
         IF (level.GE.2) THEN
            timername = timername(1:3*(level-2))//'+--'
         END IF
         timername=TRIM(ADJUSTL(timername))//TRIM(ADJUSTL(timer%name))
         timerNameLength = LEN(TRIM(ADJUSTL(timername)))

         DO n = 1, timerNameLength
            IF (timername(n:n).EQ.'&') THEN
               timername(n:n) = '+'
            END IF
         END DO

         WRITE(attributes(1),'(a)') TRIM(ADJUSTL(timername))
         WRITE(attributes(2),'(f12.3)') time
         IF(timer%n_subtimers.EQ.0) THEN
            CALL writeXMLElementForm('timer',(/'name ','value'/),attributes,&
                                         reshape((/14+20-3*level,55-timerNameLength,timerNameLength,12/),(/2,2/)))
         ELSE
            CALL openXMLElementForm('compositeTimer',(/'name ','value'/),attributes,&
                                        reshape((/5+20-3*level,55-timerNameLength,timerNameLength,12/),(/2,2/)))
            DO n = 1, timer%n_subtimers
               CALL privWriteTimesXML(timer%subtimer(n)%p,level+1)
            END DO
            CALL closeXMLElement('compositeTimer')
         END IF
      END SUBROUTINE privWriteTimesXML

      SUBROUTINE resetIterationDependentTimers()

         IMPLICIT NONE

         INTEGER                :: fn,irank=0
         LOGICAL                :: l_out
         TYPE(t_timer), POINTER :: timer, parenttimer
#ifdef CPP_MPI
         include "mpif.h"
         INTEGER::err,isize

         CALL MPI_COMM_RANK(MPI_COMM_WORLD,irank,err)
#endif
         IF (irank.NE.0) RETURN
         IF (.NOT.associated(globaltimer)) RETURN !write nothing if no timing recorded

         timer => globaltimer
         DO WHILE (TRIM(ADJUSTL(timer%name)).NE.'Iteration')
            IF(timer%n_subtimers.EQ.0) RETURN
            timer => timer%subtimer(1)%p
         END DO

         parenttimer => timer%parenttimer
         CALL resetSubtimers(timer%parenttimer)
         ALLOCATE(parenttimer%subtimer(5))
         
      END SUBROUTINE resetIterationDependentTimers

      RECURSIVE SUBROUTINE resetSubtimers (timer)

         IMPLICIT NONE

         TYPE(t_timer),INTENT(INOUT)    :: timer

         INTEGER :: n

         DO n = 1, timer%n_subtimers
            CALL resetSubtimers (timer%subtimer(n)%p)
            DEALLOCATE(timer%subtimer(n)%p)
         END DO
         DEALLOCATE(timer%subtimer)
         timer%n_subtimers = 0
      END SUBROUTINE resetSubtimers

      !>
      !<-- private function timestring

      FUNCTION timestring(time,ttime,level)
      REAL,INTENT(IN)    :: time,ttime
      CHARACTER*(90)     :: timestring
      INTEGER,INTENT(IN) ::level
!     .. Local Scalars ..
      REAL :: rest,seconds
      INTEGER :: ihours,iminutes
      CHARACTER(LEN=90)::formatstring
      rest = time

      ihours = int(rest/3600.0)
      rest = rest - real(ihours)*3600

      iminutes = int(rest/60.0)
      seconds = rest - real(iminutes)*60
      IF (ttime<0) THEN
        formatstring="(f9.2,'sec = ',i3,'h ',i2,'min ',i2,'sec')"
        WRITE (timestring,FMT = formatstring) time,ihours,iminutes,int(seconds)
       ELSE
        WRITE(formatstring,"(a,i0,a)") "(f9.2,'sec= ',i3,'h ',i2,'min ',i2,'sec ->',",(2*level-1),"x,f5.1,'%')"
        WRITE (timestring,FMT = formatstring) time,ihours,iminutes,int(seconds),time/ttime*100.0
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
      cputime=MPI_WTIME()
#elif __INTEL_COMPILER
      cputime=rtc()
#elif _OPENMP
      cputime = omp_get_wtime ( )
#else
     !use f95 intrinsic function
      CALL cpu_time(cputime)
#endif
      END FUNCTION cputime
      !>
      SUBROUTINE juDFT_time_lastlocation(PE)
      CHARACTER(LEN=4),INTENT(IN):: PE
      IF (lastline>0) THEN
          WRITE(0,*) PE,"Last kown location:"
          WRITE(0,*) PE,"File:",trim(lastfile),":",lastline
          WRITE(0,*) PE,"*****************************************"
      ENDIF
      END SUBROUTINE juDFT_time_lastlocation
      END MODULE m_juDFT_time

