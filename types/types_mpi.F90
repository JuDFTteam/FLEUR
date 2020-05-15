!--------------------------------------------------------------------------------
! Copyright (c) 2016 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------
MODULE m_types_mpi
   TYPE t_mpi
      !k-point parallelism
      INTEGER :: mpi_comm !< replaces MPI_COMM_WORLD
      INTEGER :: irank    !< rank of task in mpi_comm
      INTEGER :: isize    !< no of tasks in mpi_comm
      INTEGER, ALLOCATABLE :: k_list(:)
      !Eigenvalue parallelism
      INTEGER :: sub_comm !< Sub-Communicator for eigenvalue parallelization (all PE working on same k-point)
      INTEGER :: n_rank   !< rank in sub_comm
      INTEGER :: n_size   !< PE per kpoint, i.e. "isize" for eigenvalue parallelization
      INTEGER, ALLOCATABLE :: ev_list(:)
      !Communicator for PE on same node
      INTEGER :: mpi_comm_same_node
   CONTAINS
      procedure :: set_errhandler    => t_mpi_set_errhandler
      procedure :: is_root => mpi_is_root
   END TYPE t_mpi
contains
   function mpi_is_root(mpi) result(is_root)
      implicit none 
      class(t_mpi), intent(in) :: mpi
      logical :: is_root 
      is_root = mpi%irank == 0
   end function mpi_is_root

   subroutine t_mpi_set_errhandler(self)
      use m_judft
#ifdef CPP_MPI
      use mpi
#endif
      implicit none 
      class(t_mpi), intent(in) :: self

#ifdef CPP_MPI
      integer                  :: err_handler, ierr 

      call MPI_Comm_create_errhandler(judft_mpi_error_handler, err_handler, ierr)
      if(ierr /= 0) call judft_error("Can't create Error handler")

      call MPI_Comm_Set_Errhandler(MPI_COMM_WORLD, err_handler, ierr)
      if(ierr /= 0) call judft_error("Can't assign Error handler to MPI_COMM_WORLD")

      call MPI_Comm_Set_Errhandler(self%mpi_comm, err_handler, ierr)
      if(ierr /= 0) call judft_error("Can't assign Error handler to self%mpi_comm")

      call MPI_Comm_Set_Errhandler(self%sub_comm, err_handler, ierr)
      if(ierr /= 0) call judft_error("Can't assign Error handler to self%sub_comm")
#endif
   end subroutine t_mpi_set_errhandler

   subroutine judft_mpi_error_handler(comm, error_code) 
      use m_judft
      implicit none 
      integer, intent(in) :: comm, error_code 

      call judft_error("MPI failed with Error_code = " // int2str(error_code))
   end subroutine judft_mpi_error_handler
END MODULE m_types_mpi
