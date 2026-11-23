!--------------------------------------------------------------------------------
! Copyright (c) 2026 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------
PROGRAM fleurjob
    USE m_fleur_jobs
    USE m_juDFT
    IMPLICIT NONE
    TYPE(t_job),ALLOCATABLE::jobs(:)
    logical :: l_mpi_multithreaded
    CALL fleur_job_init(l_mpi_multithreaded)
    CALL fleur_job_arguments(jobs)
    CALL fleur_job_distribute(jobs)
    CALL fleur_job_execute(jobs, l_mpi_multithreaded)
END PROGRAM fleurjob
