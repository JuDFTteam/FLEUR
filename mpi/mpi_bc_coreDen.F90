!--------------------------------------------------------------------------------
! Copyright (c) 2016 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------

MODULE m_mpi_bc_coreden
CONTAINS
   SUBROUTINE mpi_bc_coreden(mpi,atoms,input,&
                             rhcs,tecs,qints)

   USE m_types
   IMPLICIT NONE
   INCLUDE 'mpif.h'

   TYPE(t_mpi),INTENT(IN)       :: mpi
   TYPE(t_atoms),INTENT(IN)     :: atoms
   TYPE(t_input),INTENT(IN)     :: input
   

   REAL, INTENT(INOUT) :: rhcs(atoms%jmtd,atoms%ntype,input%jspins)
   REAL, INTENT(INOUT) :: tecs(atoms%ntype,input%jspins)
   REAL, INTENT(INOUT) :: qints(atoms%ntype,input%jspins)

   INTEGER :: n, ierr(3)

    n = atoms%jmtd * atoms%ntype * input%jspins
    CALL MPI_BCAST(rhcs,n,MPI_DOUBLE,0,mpi%mpi_comm,ierr)
    n = atoms%ntype * input%jspins
    CALL MPI_BCAST(tecs,n,MPI_DOUBLE,0,mpi%mpi_comm,ierr)
    n = atoms%ntype * input%jspins
    CALL MPI_BCAST(qints,n,MPI_DOUBLE,0,mpi%mpi_comm,ierr)

   END SUBROUTINE mpi_bc_coreden
END MODULE m_mpi_bc_coreden
