!--------------------------------------------------------------------------------
! Copyright (c) 2016 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------

      MODULE m_juDFT_init

      USE m_judft_time
      USE m_judft_sysinfo
      USE m_judft_stop
      USE m_judft_args
      IMPLICIT NONE
      PRIVATE
      PUBLIC juDFT_init
      CONTAINS
     
      SUBROUTINE juDFT_init()
        IF (.NOT.judft_was_argument("-debugtime")) &
             CALL signal_handler()
      CALL checkstack()
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

      ! NOTE: The intel_signal_handler has to be outside the module
      !       as the OS has to have it under a certain name that
      !       would be changed if it would be defined in the module.
#ifdef __INTEL_COMPILER
      FUNCTION intel_signal_handler(signal)
      USE m_judft_time
      USE m_judft_sysinfo
      IMPLICIT NONE
      INTEGER :: signal
      INTEGER :: intel_signal_handler
#ifdef CPP_MPI
      include "mpif.h"
      INTEGER:: irank,ierr

      CALL MPI_COMM_RANK (MPI_COMM_WORLD,irank,ierr)
      WRITE(*,*) "Signal ",signal," detected on PE:",irank
#else
      WRITE(*,*) "Signal detected:",signal
#endif
      WRITE(*,*) "This might be due to either:"
      WRITE(*,*) " - A bug in FLEUR"
      WRITE(*,*) " - Your job running out of memory"
      WRITE(*,*) " - Your job got killed externally (e.g. no cpu-time left)"
      WRITE(*,*) " - ...." 
      WRITE(*,*) "Please check and report if you believe you found a bug"
      CALL writetimes()
      CALL PRINT_memory_info()
      STOP "Signal"
      intel_signal_handler=0
      END FUNCTION intel_signal_handler
#endif
