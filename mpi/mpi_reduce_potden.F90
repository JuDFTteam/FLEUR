!--------------------------------------------------------------------------------
! Copyright (c) 2016 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------

MODULE m_mpi_reduce_potden
#ifdef CPP_MPI 
  use mpi 
#endif
CONTAINS

  SUBROUTINE mpi_reduce_potden( fmpi, stars, sphhar, atoms, input, vacuum,   noco, potden )

    ! It is assumed that, if some quantity is allocated for some fmpi rank, that it is also allocated on fmpi rank 0. 

    USE m_types
    USE m_constants
    USE m_juDFT
    IMPLICIT NONE

    TYPE(t_mpi),     INTENT(IN)     :: fmpi
     
    TYPE(t_input),   INTENT(IN)     :: input
    TYPE(t_vacuum),  INTENT(IN)     :: vacuum
    TYPE(t_noco),    INTENT(IN)     :: noco
    TYPE(t_stars),   INTENT(IN)     :: stars
    TYPE(t_sphhar),  INTENT(IN)     :: sphhar
    TYPE(t_atoms),   INTENT(IN)     :: atoms
    TYPE(t_potden),  INTENT(INOUT)  :: potden
    
    INTEGER              :: n
    INTEGER              :: ierr
    REAL,    ALLOCATABLE :: r_b(:)

    ! reduce pw
    n = stars%ng3 * size( potden%pw, 2 )
    allocate( r_b(n) )
    call MPI_REDUCE( potden%pw, r_b, n, MPI_DOUBLE_COMPLEX, MPI_SUM, 0, fmpi%mpi_comm, ierr )
    if( fmpi%irank == 0 ) call zcopy( n, r_b, 1, potden%pw, 1 )
    deallocate( r_b )

    ! reduce mt
    n = atoms%jmtd * ( sphhar%nlhd + 1 ) * atoms%ntype * input%jspins
    allocate( r_b(n) )
    call MPI_REDUCE( potden%mt, r_b, n, MPI_DOUBLE_PRECISION, MPI_SUM, 0, fmpi%mpi_comm, ierr )
    if( fmpi%irank == 0 ) call dcopy( n, r_b, 1, potden%mt, 1 )
    deallocate( r_b )

    ! reduce pw_w
    if( allocated( potden%pw_w ) ) then
      n = stars%ng3 * size( potden%pw_w, 2 )
      allocate( r_b(n) )
      call MPI_REDUCE( potden%pw_w, r_b, n, MPI_DOUBLE_COMPLEX, MPI_SUM, 0, fmpi%mpi_comm, ierr )
      if( fmpi%irank == 0 ) call zcopy( n, r_b, 1, potden%pw_w, 1 )
      deallocate( r_b )
    end if

    ! reduce vac
    if( allocated( potden%vac ) ) then
      n = vacuum%nmzd * stars%ng2 * 2 * size( potden%vac, 4 )
      allocate( r_b(n) )
      call MPI_REDUCE( potden%vac, r_b, n, MPI_DOUBLE_COMPLEX, MPI_SUM, 0, fmpi%mpi_comm, ierr )
      if( fmpi%irank == 0 ) call zcopy( n, r_b, 1, potden%vac, 1 )
      deallocate( r_b )
    end if

    ! reduce mmpMat
    if( allocated( potden%mmpMat ) ) then
      n = size( potden%mmpMat, 1 ) * size( potden%mmpMat, 2 ) * size( potden%mmpMat, 3 ) * size( potden%mmpMat, 4 )
      allocate( r_b(n) )
      call MPI_REDUCE( potden%mmpMat, r_b, n, MPI_DOUBLE_COMPLEX, MPI_SUM, 0, fmpi%mpi_comm, ierr )
      if( fmpi%irank == 0 ) call zcopy( n, r_b, 1, potden%mmpMat, 1 )
      deallocate( r_b )
    end if

  END SUBROUTINE mpi_reduce_potden

END MODULE m_mpi_reduce_potden
