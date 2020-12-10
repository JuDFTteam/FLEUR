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

      integer     :: min_band = huge(1)
      integer     :: max_band = -1 
      integer     :: nk       = -1
#ifdef CPP_MPI
      integer                :: comm = MPI_COMM_NULL ! communicator for this t_eigvec
#else
      integer                :: comm = -1
#endif
      type(t_mat) :: mat
   contains 
      procedure :: create_comm => type_eigvec_create_comm
   end type t_eigvec 
contains 
   subroutine type_eigvec_create_comm(eigvec)
      use m_types_mpi
      implicit none
      class(t_eigvec), intent(inout)     :: eigvec

#ifdef CPP_MPI
      integer :: color
      
      if(eigvec%comm /= MPI_COMM_NULL) then
         color = merge(1,2,eigvec%l_participate)
         call judft_comm_split(MPI_COMM_WORLD, color, 1, eigvec%comm)
      endif
#endif
   end subroutine type_eigvec_create_comm
end module m_types_eigvec