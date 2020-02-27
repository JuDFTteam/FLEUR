!--------------------------------------------------------------------------------
! Copyright (c) 2016 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------
MODULE m_types_hybmpi
   TYPE t_hybmpi
      INTEGER :: comm
      INTEGER :: rank
      INTEGER :: size
      type(t_hybmpi), allocatable :: subcomm
   contains
      procedure :: copy_mpi => t_hybmpi_copy_mpi
      procedure :: barrier => t_hybmpi_barrier
      procedure :: split => t_hybmpi_split
   END TYPE t_hybmpi
contains
   subroutine t_hybmpi_copy_mpi(hybmpi, mpi)
      use m_types_mpi
      implicit none
      class(t_hybmpi), intent(inout) :: hybmpi
      type(t_mpi), intent(in)        :: mpi

      hybmpi%comm = mpi%mpi_comm
      hybmpi%size = mpi%isize
      hybmpi%rank = mpi%irank
   end subroutine

   subroutine t_hybmpi_split(hybmpi, color, key)
      use m_judft
      implicit NONE
      class(t_hybmpi), intent(inout) :: hybmpi
      integer, intent(in)            :: color
      integer, intent(in), optional  :: key
      integer :: ierr, use_key

      use_key = MERGE(key, hybmpi%rank, present(key))
      allocate(hybmpi%subcomm)
      call MPI_Comm_split(hybmpi%comm, color, use_key, hybmpi%subcomm, ierr)
      if(ierr /= 0) call judft_error("MPI splitting failed")

      call MPI_comm_rank(hybmpi%subcomm, hybmpi%subcomm%rank, ierr)
      if(ierr /= 0) call judft_error("MPI ranking failed")

      call MPI_comm_size(hybmpi%subcomm, hybmpi%subcomm%size, ierr)
      if(ierr /= 0) call judft_error("MPI sizing failed")
   end subroutine t_hybmpi_split

   subroutine t_hybmpi_barrier(hybmpi)
      use m_judft
      implicit none
      class(t_hybmpi), intent(inout) :: hybmpi
      integer :: ierr

      call MPI_Barrier(hybmpi%comm, ierr)
      if(ierr /= 0) call juDFT_error("barrier failed on process: " // &
                                      int2str(hybmpi%rank))
   end subroutine t_hybmpi_barrier
END MODULE m_types_hybmpi
