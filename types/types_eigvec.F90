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

      integer     :: min_band = huge(1)
      integer     :: max_band = -1 
      integer     :: nk = -1, jsp = -1
#ifdef CPP_MPI
      integer                :: comm = MPI_COMM_NULL ! communicator for this t_eigvec
#else
      integer                :: comm = -1
#endif
      type(t_mat) :: mat
   contains 
   end type t_eigvec 
contains 

end module m_types_eigvec