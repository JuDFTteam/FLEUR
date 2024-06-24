!--------------------------------------------------------------------------------
! Copyright (c) 2016 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------

MODULE m_mpi_bc_st
#ifdef CPP_MPI
  use mpi
#endif
  !**********************************************************************
  !     mpi_bc_st :  broadcast all information for qpw_to_nmt
  !     mpi_col_st:  collect the density from pe's 
  !**********************************************************************
CONTAINS
  SUBROUTINE mpi_bc_st(fmpi,stars,qpwc)
    !
    USE m_types
    IMPLICIT NONE

    TYPE(t_mpi),INTENT(IN)     :: fmpi
    TYPE(t_stars),INTENT(IN)   :: stars
    !     ..
    !     .. Array Arguments ..
    COMPLEX :: qpwc(stars%ng3)
    !     ..
    !     ..
    !     .. Local Arrays ..
    INTEGER ierr
    !
    ! -> Broadcast the arrays:
#ifdef CPP_MPI
    CALL MPI_BCAST(qpwc,stars%ng3,MPI_DOUBLE_COMPLEX,0,fmpi%mpi_comm,ierr)
#endif

  END SUBROUTINE mpi_bc_st
  !*********************************************************************
  SUBROUTINE mpi_col_st(fmpi,atoms,sphhar,rho)
    !
    USE m_types
    IMPLICIT NONE

    TYPE(t_mpi),INTENT(IN)     :: fmpi
    TYPE(t_sphhar),INTENT(IN)  :: sphhar
    TYPE(t_atoms),INTENT(IN)   :: atoms
    !     ..
    !     .. Scalar Arguments ..
    REAL, INTENT (INOUT) :: rho(atoms%jmtd,0:sphhar%nlhd,atoms%ntype)

    INTEGER n
    INTEGER ierr
    REAL, ALLOCATABLE :: r_b(:)
#ifdef CPP_MPI
    n = atoms%jmtd*(sphhar%nlhd+1)*atoms%ntype
    ALLOCATE(r_b(n))
    CALL MPI_REDUCE(rho,r_b,n,MPI_DOUBLE_PRECISION,MPI_SUM,0,&
         &                                       fmpi%mpi_comm,ierr)
    IF (fmpi%irank == 0) rho=reshape(r_b,(/atoms%jmtd,1+sphhar%nlhd,atoms%ntype/))

    DEALLOCATE(r_b) 
#endif
  END SUBROUTINE mpi_col_st
  !*********************************************************************
END MODULE m_mpi_bc_st
