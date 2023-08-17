!--------------------------------------------------------------------------------
! Copyright (c) 2016 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------

MODULE m_juDFT_init
#ifdef CPP_MPI
      use mpi
#endif
      USE m_judft_time
      USE m_judft_sysinfo
      USE m_judft_stop
      USE m_judft_args
      USE m_juDFT_internalParams
      IMPLICIT NONE
      PRIVATE
      PUBLIC juDFT_init
      CONTAINS

      SUBROUTINE juDFT_init(outUnit,l_checkStack)

         INTEGER, INTENT(IN) :: outUnit
         LOGICAL, INTENT(IN) :: l_checkStack

         juDFT_outUnit = outUnit
         IF (.NOT.judft_was_argument("-debugtime")) CALL signal_handler()
         IF (l_checkStack) CALL checkstack()
#if defined(CPP_PATCH_INTEL)&&defined(__INTEL_COMPILER)
       call patch_intel()
#endif

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

#if defined(CPP_PATCH_INTEL)&&defined(__INTEL_COMPILER)
      subroutine patch_intel()
        !we try to patch the intel libraries to overwrite determination of 'INTEL' brand
        !otherwise performance on AMD CPUs is bad.
        INTERFACE
          subroutine mkl_patch() BIND(C, name="intel_mkl_patch")
          END subroutine
        END INTERFACE
        INTERFACE
          subroutine cpu_patch() BIND(C, name="intel_cpu_patch")
          END subroutine
        END INTERFACE
        print *,"INTEL PATCH applied"
        call cpu_patch()
        call mkl_patch()
      end subroutine


#endif



      END MODULE m_juDFT_init

      ! NOTE: The intel_signal_handler has to be outside the module
      !       as the OS has to have it under a certain name that
      !       would be changed if it would be defined in the module.
#ifdef __INTEL_COMPILER
      FUNCTION intel_signal_handler(signal)
#ifdef CPP_MPI
      use mpi
#endif
      USE iso_fortran_env
      USE m_judft_time
      USE m_judft_sysinfo
      IMPLICIT NONE
      INTEGER :: signal
      INTEGER :: intel_signal_handler
#ifdef CPP_MPI
      INTEGER:: irank,ierr
      LOGICAL:: l_mpi_init
      CALL MPI_initialized(l_mpi_init,ierr)
      IF (l_mpi_init) THEN
         CALL MPI_COMM_RANK (MPI_COMM_WORLD,irank,ierr)
         WRITE(*,*) "Signal ",signal," detected on PE:",irank
      ELSE
         WRITE(*,*) "Signal detected:",signal
      END IF
#else
      WRITE(*,*) "Signal detected:",signal
#endif
      WRITE(*,*) "This might be due to either:"
      WRITE(*,*) " - A bug"
      WRITE(*,*) " - Your job running out of memory"
      WRITE(*,*) " - Your job got killed externally (e.g. no cpu-time left)"
      WRITE(*,*) " - ...."
      WRITE(*,*) "Please check and report if you believe you found a bug"
      CALL writetimes()
      CALL PRINT_memory_info(output_unit,.true.)
#ifdef CPP_MPI
      IF (l_mpi_init) CALL MPI_ABORT(MPI_COMM_WORLD,0,ierr)
#endif
      STOP "Signal"
      intel_signal_handler=0
      END FUNCTION intel_signal_handler
#endif
