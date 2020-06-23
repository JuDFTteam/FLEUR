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
   contains
      procedure :: copy_mpi => t_hybmpi_copy_mpi
      procedure :: barrier => t_hybmpi_barrier
   END TYPE t_hybmpi
contains
   subroutine t_hybmpi_copy_mpi(glob_mpi, mpi)
      use m_types_mpi
      implicit none
      class(t_hybmpi), intent(inout) :: glob_mpi
      type(t_mpi), intent(in)        :: mpi

      glob_mpi%comm = mpi%mpi_comm
      glob_mpi%size = mpi%isize
      glob_mpi%rank = mpi%irank
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
END MODULE m_types_hybmpi
