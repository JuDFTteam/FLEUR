!--------------------------------------------------------------------------------
! Copyright (c) 2026 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------
MODULE m_gfleur_driver
   !**********************************************************************
   !     Driver of the embedding Green-function code on the basis of
   !     the FLEUR code.
   !
   !     Daniel Wortmann 2000-2004
   !     Layer decomposition: Frank Freimuth, November 2007
   !
   !     CALL tree:
   !     gfleur
   !     |-gf_init             (per-layer init from inp.xml + gf_inp)
   !     |-gfleur_execute      (SCF / transport / CBS ... loops)
   !**********************************************************************
#ifdef CPP_MPI
   USE mpi
#endif
   USE m_juDFT
   IMPLICIT NONE
   PRIVATE
   PUBLIC :: gfleur_run
CONTAINS
   SUBROUTINE gfleur_run(l_mpi_multithreaded)
      USE m_types
      USE m_constants
      USE m_gf_types
      USE m_gf_init
      LOGICAL, INTENT(IN) :: l_mpi_multithreaded

      TYPE(t_gfmpi)   :: gmpi
      TYPE(t_embinp)  :: embinp
      TYPE(t_gfmix)   :: gfmix
      TYPE(t_layers)  :: layers
      TYPE(t_gflayer), ALLOCATABLE :: ld(:)
      TYPE(t_kpts)    :: gkpts

      gmpi%fmpi%l_mpi_multithreaded = l_mpi_multithreaded
#ifdef CPP_MPI
      gmpi%fmpi%mpi_comm = MPI_COMM_WORLD
#else
      gmpi%fmpi%mpi_comm = 1
#endif

      CALL gf_init(gmpi, embinp, gfmix, layers, ld, gkpts)

      CALL gfleur_execute(gmpi, embinp, gfmix, layers, ld, gkpts)
   END SUBROUTINE gfleur_run

   SUBROUTINE gfleur_execute(gmpi, embinp, gfmix, layers, ld, gkpts)
      !The self-consistency / transport / CBS loops of gfleur.
      !
      !Port status: the initialization chain is fully ported; the
      !computational kernels (H/S + T-matrix setup, propagation,
      !Green-function composition, density generation, layered Coulomb
      !solver, mixing) are re-enabled step by step on top of the modern
      !FLEUR kernels.
      USE m_types
      USE m_constants
      USE m_gf_types
      IMPLICIT NONE
      TYPE(t_gfmpi), INTENT(INOUT)  :: gmpi
      TYPE(t_embinp), INTENT(INOUT) :: embinp
      TYPE(t_gfmix), INTENT(INOUT)  :: gfmix
      TYPE(t_layers), INTENT(IN)    :: layers
      TYPE(t_gflayer), INTENT(INOUT) :: ld(:)
      TYPE(t_kpts), INTENT(IN)      :: gkpts

      INTEGER :: i

      IF (gmpi%pe0) THEN
         WRITE (oUnit, *)
         WRITE (oUnit, "(a)") "GFLEUR setup summary:"
         WRITE (oUnit, "(a,i0)") "  layers : ", layers%num_layers
         WRITE (oUnit, "(a,i0)") "  k-points (2D BZ): ", gkpts%nkpt
         WRITE (oUnit, "(a,i0)") "  spins  : ", ld(1)%fi%input%jspins
         DO i = 1, layers%num_layers
            WRITE (oUnit, "(a,i0,a,i0,a,i0,a,i0,a,i0)") "  layer ", i, &
               ": ntype=", ld(i)%fi%atoms%ntype, " nat=", ld(i)%fi%atoms%nat, &
               " ng3=", ld(i)%stars%ng3, " nvd=", ld(i)%nvd
         END DO
      END IF

      !The computational kernels are not yet wired up in the port -
      !stop cleanly after the verified initialization.
      CALL juDFT_end("GFLEUR init done - computational kernels not yet ported", &
                     gmpi%fmpi%irank)
   END SUBROUTINE gfleur_execute
END MODULE m_gfleur_driver

PROGRAM gfleur
   USE m_gfleur_driver
   USE m_fleur_jobs, ONLY: fleur_job_init
   USE m_juDFT
   IMPLICIT NONE
   LOGICAL :: l_mpi_multithreaded
   CALL fleur_job_init(l_mpi_multithreaded)
   CALL gfleur_run(l_mpi_multithreaded)
   CALL juDFT_end("GFLEUR done")
END PROGRAM gfleur
