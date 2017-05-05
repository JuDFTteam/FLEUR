!--------------------------------------------------------------------------------
! Copyright (c) 2016 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------

MODULE m_mpi_bc_pot
CONTAINS
   SUBROUTINE mpi_bc_pot(mpi,stars,sphhar,atoms,input,vacuum,&
                         iter,fr,fpw,fz,fzxy)

   USE m_types
   IMPLICIT NONE
   INCLUDE 'mpif.h'

   TYPE(t_mpi),INTENT(IN)        :: mpi
   TYPE(t_input),INTENT(IN)      :: input
   TYPE(t_vacuum),INTENT(IN)     :: vacuum
   TYPE(t_stars),INTENT(IN)      :: stars
   TYPE(t_sphhar),INTENT(IN)     :: sphhar
   TYPE(t_atoms),INTENT(IN)      :: atoms

   INTEGER, INTENT (INOUT) :: iter

   COMPLEX, INTENT (INOUT) :: fpw(stars%ng3,input%jspins)
   COMPLEX, INTENT (INOUT) :: fzxy(vacuum%nmzxyd,stars%ng2-1,2,input%jspins)
   REAL,    INTENT (INOUT) :: fr(atoms%jmtd,0:sphhar%nlhd,atoms%ntype,input%jspins)
   REAL,    INTENT (INOUT) :: fz(vacuum%nmzd,2,input%jspins)

   INTEGER :: n, ierr(3)

   CALL MPI_BCAST(iter,1,MPI_INTEGER,0,mpi%mpi_comm,ierr)

    n = stars%ng3 * input%jspins
    CALL MPI_BCAST(fpw,n,MPI_DOUBLE_COMPLEX,0,mpi%mpi_comm,ierr)
    n = vacuum%nmzxyd * (stars%ng2-1) * 2 * input%jspins
    CALL MPI_BCAST(fzxy,n,MPI_DOUBLE_COMPLEX,0,mpi%mpi_comm,ierr)
    n = atoms%jmtd * (sphhar%nlhd+1) * atoms%ntype * input%jspins
    CALL MPI_BCAST(fr,n,MPI_DOUBLE,0,mpi%mpi_comm,ierr)
    n = vacuum%nmzd * 2 * input%jspins
    CALL MPI_BCAST(fz,n,MPI_DOUBLE,0,mpi%mpi_comm,ierr)

   END SUBROUTINE mpi_bc_pot
END MODULE m_mpi_bc_pot
