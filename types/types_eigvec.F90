module m_types_eigvec
   use m_types_mat 
   USE m_types_mpi
   USE m_types_fleurinput
#ifdef CPP_MPI
   use mpi
#endif
   implicit none

   type t_eigvec 
      logical     :: l_participate = .False. 
      logical     :: l_recv        = .False.

      integer, allocatable :: root_pe(:)

      integer     :: nk = -1, jsp = -1
#ifdef CPP_MPI
      integer                :: comm = MPI_COMM_NULL ! communicator for this t_eigvec
#else
      integer                :: comm = -1
#endif
      type(t_mat) :: mat
   contains 
      procedure :: free => free_eigvec
   end type t_eigvec 
contains 
   subroutine free_eigvec(eigvec)
      implicit NONE
      class(t_eigvec) :: eigvec 

      integer :: ierr 

      eigvec%l_participate = .False.
      eigvec%l_recv = .false.

      if(allocated(eigvec%root_pe)) deallocate(eigvec%root_pe) 

      eigvec%nk  = -1 
      eigvec%jsp = -1 

#ifdef CPP_MPI
      if(eigvec%comm /= MPI_COMM_NULL) call MPI_Comm_free(eigvec%comm, ierr)
#endif
      call eigvec%mat%free()
   end subroutine free_eigvec
end module m_types_eigvec
