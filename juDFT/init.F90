!--------------------------------------------------------------------------------
! Copyright (c) 2016 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------

      MODULE m_juDFT_init
      USE m_judft_time
      IMPLICIT NONE
      PRIVATE
      PUBLIC juDFT_init,juDFT_was_argument
      CONTAINS

      FUNCTION juDFT_was_argument(arg) RESULT(OK)
      IMPLICIT NONE
      CHARACTER(len=*),INTENT(IN)::arg
      LOGICAL ok

      INTEGER:: i
      CHARACTER(LEN=30)::str
      ok=.false.
      DO i=1,command_argument_count()
         CALL get_command_argument(i,str)
         IF(adjustl(str)==adjustl(arg)) ok=.true.
      ENDDO
      IF (ok) return
      !Test for environment variable as well
      CALL get_environment_variable("juDFT",str,status=i)
      IF (i==0) ok=index(str,adjustl(arg))>0
 
      END FUNCTION

      SUBROUTINE juDFT_init()
      CALL signal_handler()
      END SUBROUTINE juDFT_init

      SUBROUTINE signal_handler()
      !Installs custom handlers for SIGTERM,SIGSEGV
#ifdef __INTEL_COMPILER
      USE ifport
      INTEGER :: result
      EXTERNAL intel_signal_handler
      result=signal(SIGTERM,intel_signal_handler,-1)
      result=signal(SIGSEGV,intel_signal_handler,-1)

#endif
      END SUBROUTINE signal_handler
      END MODULE m_juDFT_init

#ifdef __INTEL_COMPILER
      FUNCTION intel_signal_handler(signal)
      USE m_judft_time
      IMPLICIT NONE
      INTEGER :: signal
      INTEGER :: intel_signal_handler
#ifdef CPP_MPI
      include "mpif.h"
      INTEGER:: irank,ierr

      CALL MPI_COMM_RANK (MPI_COMM_WORLD,irank,ierr)
      WRITE(*,*) "Signal ",signal," detected on PE:",irank
      CALL writetimes()
#else
      WRITE(*,*) "Signal detected:",signal
      CALL writetimes()
#endif
      STOP "Signal"
      intel_signal_handler=0
      END FUNCTION intel_signal_handler
#endif
