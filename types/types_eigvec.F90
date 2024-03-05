module m_types_eigvec
   use m_types_mat 

#ifdef CPP_MPI
   use mpi
#endif
   implicit none
   private 
   type,public:: t_eigvec 
      logical     :: l_participate = .False. 
      logical     :: l_recv        = .False.
      integer     :: nk = -1, jsp = -1
      
#ifdef CPP_MPI
      integer                :: comm = MPI_COMM_NULL ! communicator for this t_eigvec
#else
      integer                :: comm = -1
#endif
      integer, allocatable :: root_pe(:)
      type(t_mat) :: mat
   end type t_eigvec 
   
end module m_types_eigvec