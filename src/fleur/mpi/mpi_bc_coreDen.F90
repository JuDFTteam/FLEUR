!--------------------------------------------------------------------------------
! Copyright (c) 2026 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------

MODULE m_mpi_bc_coreden
#ifdef CPP_MPI
   use mpi
#endif
   use m_types_atoms
   use m_types_input
   use m_types_mpi
   implicit none
   private
   public :: mpi_bc_coreden
CONTAINS
   SUBROUTINE mpi_bc_coreden(fmpi,atoms,input,&
                             rhcs,tecs,qints)

   IMPLICIT NONE

   TYPE(t_mpi),INTENT(IN)       :: fmpi
   TYPE(t_atoms),INTENT(IN)     :: atoms
   TYPE(t_input),INTENT(IN)     :: input
   

   REAL, INTENT(INOUT) :: rhcs(atoms%jmtd,atoms%ntype,input%jspins)
   REAL, INTENT(INOUT) :: tecs(atoms%ntype,input%jspins)
   REAL, INTENT(INOUT) :: qints(atoms%ntype,input%jspins)
#ifdef CPP_MPI
   INTEGER :: n, ierr

    n = atoms%jmtd * atoms%ntype * input%jspins
    CALL MPI_BCAST(rhcs,n,MPI_DOUBLE,0,fmpi%mpi_comm,ierr)
    n = atoms%ntype * input%jspins
    CALL MPI_BCAST(tecs,n,MPI_DOUBLE,0,fmpi%mpi_comm,ierr)
    n = atoms%ntype * input%jspins
    CALL MPI_BCAST(qints,n,MPI_DOUBLE,0,fmpi%mpi_comm,ierr)
#endif
   END SUBROUTINE mpi_bc_coreden
END MODULE m_mpi_bc_coreden
