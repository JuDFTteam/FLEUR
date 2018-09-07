!--------------------------------------------------------------------------------
! Copyright (c) 2016 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------

MODULE m_mpi_dist_forcetheorem
CONTAINS
#ifndef CPP_OLDINTEL
  SUBROUTINE mpi_dist_forcetheorem(mpi,forcetheo)
    USE m_types_mpi
    USE m_types_forcetheo, ONLY: t_forcetheo
    USE m_types_forcetheo_extended
    IMPLICIT NONE
    TYPE(t_mpi),INTENT(in)::mpi
    CLASS(t_forcetheo),ALLOCATABLE,INTENT(INOUT)::forcetheo

    INTEGER::t,ierr
#ifdef CPP_MPI
    INCLUDE 'mpif.h'
    
    IF (mpi%irank==0) THEN
       SELECT TYPE(forcetheo)
       TYPE IS (t_forcetheo)
          t=1
       TYPE IS (t_forcetheo_mae)
          t=2
       END SELECT
    ENDIF   
    CALL MPI_BCAST(t,1,MPI_INTEGER,0,mpi%mpi_comm,ierr)
    IF (mpi%irank.NE.0) THEN
       IF (ALLOCATED(forcetheo)) DEALLOCATE(forcetheo)
       SELECT CASE (t)
       CASE (1)
          ALLOCATE(t_forcetheo::forcetheo)
       CASE(2)
          ALLOCATE(t_forcetheo_mae::forcetheo)
       END SELECT
    END IF

    !now we have the correct type, now we have to distribute the data
    SELECT TYPE(forcetheo)
    TYPE IS (t_forcetheo_mae)
       CALL forcetheo%dist(mpi)
    END SELECT
#endif
  END SUBROUTINE mpi_dist_forcetheorem
#else
#endif

END MODULE m_mpi_dist_forcetheorem
