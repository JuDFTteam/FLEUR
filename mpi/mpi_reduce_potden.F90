!--------------------------------------------------------------------------------
! Copyright (c) 2016 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------

MODULE m_mpi_reduce_potden

CONTAINS

  SUBROUTINE mpi_reduce_potden( mpi, stars, sphhar, atoms, input, vacuum, oneD, noco, potden )

    ! It is assumed that, if some quantity is allocated for some mpi rank, that it is also allocated on mpi rank 0. 

#include"cpp_double.h"
    USE m_types
    USE m_constants
    USE m_juDFT
    IMPLICIT NONE

    TYPE(t_mpi),     INTENT(IN)     :: mpi
    TYPE(t_oneD),    INTENT(IN)     :: oneD
    TYPE(t_input),   INTENT(IN)     :: input
    TYPE(t_vacuum),  INTENT(IN)     :: vacuum
    TYPE(t_noco),    INTENT(IN)     :: noco
    TYPE(t_stars),   INTENT(IN)     :: stars
    TYPE(t_sphhar),  INTENT(IN)     :: sphhar
    TYPE(t_atoms),   INTENT(IN)     :: atoms
    TYPE(t_potden),  INTENT(INOUT)  :: potden
    INCLUDE 'mpif.h'
    
    INTEGER              :: n
    INTEGER              :: ierr(3)
    REAL,    ALLOCATABLE :: r_b(:)
    
    EXTERNAL CPP_BLAS_scopy,CPP_BLAS_ccopy,MPI_REDUCE

    ! reduce pw
    n = stars%ng3 * size( potden%pw, 2 )
    allocate( r_b(n) )
    call MPI_REDUCE( potden%pw, r_b, n, MPI_DOUBLE_COMPLEX, MPI_SUM, 0, mpi%mpi_comm, ierr )
    if( mpi%irank == 0 ) call CPP_BLAS_ccopy( n, r_b, 1, potden%pw, 1 )
    deallocate( r_b )

    ! reduce mt
    n = atoms%jmtd * ( sphhar%nlhd + 1 ) * atoms%ntype * input%jspins
    allocate( r_b(n) )
    call MPI_REDUCE( potden%mt, r_b, n, MPI_DOUBLE_PRECISION, MPI_SUM, 0, mpi%mpi_comm, ierr )
    if( mpi%irank == 0 ) call CPP_BLAS_scopy( n, r_b, 1, potden%mt, 1 )
    deallocate( r_b )

    ! reduce pw_w
    if( allocated( potden%pw_w ) ) then
      n = stars%ng3 * size( potden%pw_w, 2 )
      allocate( r_b(n) )
      call MPI_REDUCE( potden%pw_w, r_b, n, MPI_DOUBLE_COMPLEX, MPI_SUM, 0, mpi%mpi_comm, ierr )
      if( mpi%irank == 0 ) call CPP_BLAS_ccopy( n, r_b, 1, potden%pw_w, 1 )
      deallocate( r_b )
    end if

    ! reduce vacz
    if( allocated( potden%vacz ) ) then
      n = vacuum%nmz * 2 * size( potden%vacz, 3 )
      allocate( r_b(n) )
      call MPI_REDUCE( potden%vacz, r_b, n, MPI_DOUBLE_PRECISION, MPI_SUM, 0, mpi%mpi_comm, ierr )
      if( mpi%irank == 0 ) call CPP_BLAS_scopy( n, r_b, 1, potden%vacz, 1 )
      deallocate( r_b )
    end if

    ! reduce vacxy
    if( allocated( potden%vacxy ) ) then
      n = vacuum%nmzxy * ( stars%ng2 - 1 ) * 2 * size( potden%vacxy, 4 )
      allocate( r_b(n) )
      call MPI_REDUCE( potden%vacxy, r_b, n, MPI_DOUBLE_COMPLEX, MPI_SUM, 0, mpi%mpi_comm, ierr )
      if( mpi%irank == 0 ) call CPP_BLAS_ccopy( n, r_b, 1, potden%vacxy, 1 )
      deallocate( r_b )
    end if

    ! reduce mmpMat
    if( allocated( potden%mmpMat ) ) then
      n = size( potden%mmpMat, 1 ) * size( potden%mmpMat, 2 ) * size( potden%mmpMat, 3 ) * size( potden%mmpMat, 4 )
      allocate( r_b(n) )
      call MPI_REDUCE( potden%mmpMat, r_b, n, MPI_DOUBLE_COMPLEX, MPI_SUM, 0, mpi%mpi_comm, ierr )
      if( mpi%irank == 0 ) call CPP_BLAS_ccopy( n, r_b, 1, potden%mmpMat, 1 )
      deallocate( r_b )
    end if

  END SUBROUTINE mpi_reduce_potden

END MODULE m_mpi_reduce_potden
