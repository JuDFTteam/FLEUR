module m_create_coul_comms
   use m_types
   implicit none
contains
   subroutine create_coul_comms(hybdat, fi, fmpi)
#ifdef CPP_MPI
      use mpi
#endif
      implicit none 
      type(t_hybdat), intent(inout)   :: hybdat 
      type(t_fleurinput), intent(in)  :: fi
      TYPE(t_mpi), INTENT(IN)         :: fmpi

      integer :: ik, color, key
      logical :: i_am_root

#ifdef CPP_MPI
      do ik = 1,fi%kpts%nkpt 
         if(hybdat%coul(ik)%comm == MPI_COMM_NULL) then
            i_am_root = (fmpi%n_rank == 0) .and. any(ik == fmpi%k_list)

            if(hybdat%coul(ik)%l_participate) then
               color = 1
            else
               color = 2
            endif

            ! put the root rank on 0, others don't care
            key = merge(0, fmpi%irank+1, i_am_root)

            call judft_comm_split(MPI_COMM_WORLD, color, key, hybdat%coul(ik)%comm)
         endif
      enddo
#endif
   end subroutine
end module m_create_coul_comms