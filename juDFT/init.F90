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
      LOGICAL:: mpi_init
      CALL MPI_initialized(mpi_init,ierr)
      IF (mpi_init) THEN
         CALL MPI_COMM_RANK (MPI_COMM_WORLD,irank,ierr)
         WRITE(0,*) "Signal ",signal," detected on PE:",irank
      ELSE
         WRITE(0,*) "Signal detected:",signal
      END IF
#else
      WRITE(0,*) "Signal detected:",signal
#endif
      WRITE(0,*) "This might be due to either:"
      WRITE(0,*) " - A bug"
      WRITE(0,*) " - Your job running out of memory"
      WRITE(0,*) " - Your job got killed externally (e.g. no cpu-time left)"
      WRITE(0,*) " - ...." 
      WRITE(0,*) "Please check and report if you believe you found a bug"
      CALL writetimes()
      CALL PRINT_memory_info(0,.true.)
#ifdef CPP_MPI
      IF (mpi_init) CALL MPI_ABORT(MPI_COMM_WORLD,ierr)
#endif      
      STOP "Signal"
      intel_signal_handler=0
      END FUNCTION intel_signal_handler
#endif
