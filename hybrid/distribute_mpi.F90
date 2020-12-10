module m_distribute_mpi 
contains 
   subroutine distribute_mpi(weights, glob_mpi, group_mpi, group_rank)
      use m_types_hybmpi
      use m_types_mpi
      implicit none 
      integer, intent(in)           :: weights(:)
      type(t_hybmpi), intent(in)    :: glob_mpi
      type(t_hybmpi), intent(inout) :: group_mpi
      integer, intent(inout)        :: group_rank 

      integer, allocatable :: nprocs(:), color(:)
      integer              :: idx(1), i, j, cnt, n_grps, new_comm


      n_grps = size(weights)
      allocate(nprocs(n_grps), source=0)
      allocate(color(glob_mpi%size), source=0)
         
      do i = 1,glob_mpi%size 
         idx = minloc(1.0*nprocs/weights)
         nprocs(idx(1)) = nprocs(idx(1)) + 1
      enddo

      cnt = 1
      do i = 1,n_grps
         do j = 1,nprocs(i)
            color(cnt) = i - 1
            cnt = cnt + 1 
         enddo 
      enddo

      group_rank = color(glob_mpi%rank+1) 
      call judft_comm_split(glob_mpi%comm, group_rank, glob_mpi%rank, new_comm)

      call group_mpi%init(new_comm)
   end subroutine distribute_mpi
end module m_distribute_mpi