!--------------------------------------------------------------------------------
! Copyright (c) 2016 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------
MODULE m_types_hybmpi
#ifdef CPP_MPI
   use mpi
#endif
   TYPE t_hybmpi
      INTEGER :: comm
      INTEGER :: rank
      INTEGER :: size
   contains
      procedure :: copy_mpi => t_hybmpi_copy_mpi
      procedure :: barrier  => t_hybmpi_barrier
      procedure :: init     => t_hybmpi_init
      procedure :: root     => t_hybmpi_root
   END TYPE t_hybmpi
contains
   function t_hybmpi_root(mpi_var) result(l_root)
      implicit none
      class(t_hybmpi), intent(in) :: mpi_var 
      logical :: l_root 

      l_root = mpi_var%rank == 0
   end function t_hybmpi_root

   subroutine t_hybmpi_copy_mpi(glob_mpi, mpi_var)
      use m_types_mpi
      implicit none
      class(t_hybmpi), intent(inout) :: glob_mpi
      type(t_mpi), intent(in)        :: mpi_var

      glob_mpi%comm = mpi_var%mpi_comm
      glob_mpi%size = mpi_var%isize
      glob_mpi%rank = mpi_var%irank
   end subroutine

   subroutine t_hybmpi_barrier(glob_mpi)
      use m_judft
      implicit none
      class(t_hybmpi), intent(inout) :: glob_mpi
      integer :: ierr
#ifdef CPP_MPI
      call MPI_Barrier(glob_mpi%comm, ierr)
      if(ierr /= 0) call juDFT_error("barrier failed on process: " // &
                                      int2str(glob_mpi%rank))
#endif
   end subroutine t_hybmpi_barrier

   subroutine t_hybmpi_init(hybmpi, in_comm) 
      class(t_hybmpi), intent(inout) :: hybmpi 
      integer, intent(in)            :: in_comm 
      integer                        :: ierr

      hybmpi%comm = in_comm 

#ifdef CPP_MPI
      call MPI_Comm_size(hybmpi%comm, hybmpi%size, ierr)
      call MPI_Comm_rank(hybmpi%comm, hybmpi%rank, ierr)
#else
      hybmpi%size = 1
      hybmpi%rank = 0
#endif
   end subroutine t_hybmpi_init
END MODULE m_types_hybmpi
