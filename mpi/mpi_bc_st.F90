!--------------------------------------------------------------------------------
! Copyright (c) 2016 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------

MODULE m_mpi_bc_st
  !**********************************************************************
  !     mpi_bc_st :  broadcast all information for qpw_to_nmt
  !     mpi_col_st:  collect the density from pe's 
  !**********************************************************************
CONTAINS
  SUBROUTINE mpi_bc_st(mpi,stars,qpwc)
    !
    USE m_types
    IMPLICIT NONE

    TYPE(t_mpi),INTENT(IN)     :: mpi
    TYPE(t_stars),INTENT(IN)   :: stars
    !     ..
    !     .. Array Arguments ..
    COMPLEX :: qpwc(stars%n3d)
    !     ..
    !     ..
    !     .. Local Arrays ..
    INTEGER ierr(3)
    !     ..
    !     .. External Subroutines.. 
    EXTERNAL MPI_BCAST
    !     ..
    INCLUDE 'mpif.h'
    !
    !
    ! -> Broadcast the arrays:

    CALL MPI_BCAST(qpwc,stars%n3d,MPI_DOUBLE_COMPLEX,0,mpi%mpi_comm,ierr)

  END SUBROUTINE mpi_bc_st
  !*********************************************************************
  SUBROUTINE mpi_col_st(mpi,atoms,sphhar,rho)
    !
#include"cpp_double.h"
    USE m_types
    IMPLICIT NONE

    TYPE(t_mpi),INTENT(IN)     :: mpi
    TYPE(t_sphhar),INTENT(IN)  :: sphhar
    TYPE(t_atoms),INTENT(IN)   :: atoms
    INCLUDE 'mpif.h'
    EXTERNAL MPI_REDUCE
    !     ..
    !     .. Scalar Arguments ..
    REAL, INTENT (INOUT) :: rho(atoms%jmtd,0:sphhar%nlhd,atoms%ntype)

    INTEGER n
    INTEGER ierr(3)
    REAL, ALLOCATABLE :: r_b(:)

    n = atoms%jmtd*(sphhar%nlhd+1)*atoms%ntype
    ALLOCATE(r_b(n))
    CALL MPI_REDUCE(rho,r_b,n,MPI_DOUBLE_PRECISION,MPI_SUM,0,&
         &                                       mpi%mpi_comm,ierr)
    IF (mpi%irank == 0) rho=reshape(r_b,(/atoms%jmtd,1+sphhar%nlhd,atoms%ntype/))

    DEALLOCATE(r_b) 

  END SUBROUTINE mpi_col_st
  !*********************************************************************
END MODULE m_mpi_bc_st
