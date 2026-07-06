!--------------------------------------------------------------------------------
! Copyright (c) 2026 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------
MODULE m_gfleur_driver
   !Driver of the Green's-function embedding code. Currently a stub that
   !verifies the gfleur executable can initialize from inp.xml using the
   !standard FLEUR machinery; it grows into gf_init/gfleur_execute during
   !the port.
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
      USE m_fleur_init
      LOGICAL, INTENT(IN) :: l_mpi_multithreaded

      TYPE(t_fleurinput) :: fi
      TYPE(t_mpi)        :: fmpi
      TYPE(t_sphhar)     :: sphhar
      TYPE(t_stars)      :: stars
      TYPE(t_enpara)     :: enpara
      TYPE(t_results)    :: results
      TYPE(t_nococonv)   :: nococonv
      TYPE(t_wann)       :: wann
      TYPE(t_hybdat)     :: hybdat
      TYPE(t_mpdata)     :: mpdata
      CLASS(t_forcetheo), ALLOCATABLE :: forcetheo
      CLASS(t_xcpot), ALLOCATABLE     :: xcpot

      CHARACTER(len=100) :: filename_add

      fmpi%l_mpi_multithreaded = l_mpi_multithreaded
#ifdef CPP_MPI
      fmpi%mpi_comm = MPI_COMM_WORLD
#else
      fmpi%mpi_comm = 1
#endif

      CALL timestart("GFLEUR Initialization")
      filename_add = ""
      CALL fleur_init(fmpi, fi, sphhar, stars, nococonv, forcetheo, enpara, xcpot, &
                      results, wann, hybdat, mpdata, filename_add)
      CALL timestop("GFLEUR Initialization")

      IF (fmpi%irank == 0) THEN
         WRITE (oUnit, *) "GFLEUR stub: FLEUR initialization completed"
         WRITE (oUnit, *) "   ntype=", fi%atoms%ntype, " nat=", fi%atoms%nat
         WRITE (oUnit, *) "   ng3  =", stars%ng3, " ng2=", stars%ng2
      END IF
   END SUBROUTINE gfleur_run
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
