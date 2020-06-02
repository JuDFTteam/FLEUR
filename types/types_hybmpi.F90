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
   subroutine t_hybmpi_copy_mpi(hybmpi, mpi)
      use m_types_mpi
      implicit none
      class(t_hybmpi), intent(inout) :: hybmpi
      type(t_mpi), intent(in)        :: mpi

      hybmpi%comm = mpi%mpi_comm
      hybmpi%size = mpi%isize
      hybmpi%rank = mpi%irank
   end subroutine

   subroutine t_hybmpi_barrier(hybmpi)
      use m_judft
      implicit none
      class(t_hybmpi), intent(inout) :: hybmpi
      integer :: ierr
#ifdef CPP_MPI
      call MPI_Barrier(hybmpi%comm, ierr)
      if(ierr /= 0) call juDFT_error("barrier failed on process: " // &
                                      int2str(hybmpi%rank))
#endif
   end subroutine t_hybmpi_barrier
END MODULE m_types_hybmpi
