MODULE m_balance_barriers
#ifdef CPP_MPI 
   use mpi 
#endif
   use m_types
   use m_work_package
   use m_judft

   type t_balance_wavef
      integer :: remaining_barries
   contains 
      procedure :: init => t_balance_wavef_init
      procedure :: balance => t_balance_wavef_balance
   end type t_balance_wavef


contains 
   subroutine balance_hsfock(work_pack)
      implicit none
      class(t_work_package), intent(in) :: work_pack
      integer :: i, ierr

#ifdef CPP_MPI 
      do i =1, (work_pack%max_kpacks - work_pack%n_kpacks)
         ! call timestart("Post read_z Barrier: hsfock(dangle)")
         ! call MPI_Barrier(MPI_COMM_WORLD, ierr)
         ! call timestop("Post read_z Barrier: hsfock(dangle)")
      enddo
#endif
   end subroutine balance_hsfock

   subroutine t_balance_wavef_init(bal, fi, work_pack)
      implicit none
      class(t_balance_wavef), intent(inout) :: bal 
      type(t_fleurinput), intent(in)    :: fi
      class(t_work_package), intent(in)     :: work_pack
      integer :: barrier_calls, max_barriers, ierr, i, ik, iq, jq, n_parts, start, stride, ipart
      integer :: rank
      


#ifdef CPP_MPI 
      barrier_calls = 0

      DO i = 1,work_pack%n_kpacks
         ik = work_pack%k_packs(i)%nk
         DO jq = 1, size(work_pack%k_packs(i)%q_packs)
            n_parts = size(work_pack%k_packs(i)%q_packs(jq)%band_packs)
            start  = work_pack%k_packs(i)%q_packs(jq)%submpi%rank + 1
            stride = work_pack%k_packs(i)%q_packs(jq)%submpi%size

            do ipart = start, n_parts, stride
               barrier_calls = barrier_calls + 1 
            enddo 
         enddo 
      enddo

      CALL MPI_Comm_rank(MPI_COMM_WORLD, rank, ierr)

      call MPI_Allreduce(barrier_calls, max_barriers, 1, MPI_INTEGER, MPI_MAX, MPI_COMM_WORLD, ierr)
      bal%remaining_barries = max_barriers - barrier_calls
#else 
      bal%remaining_barries = 0
#endif
   end subroutine t_balance_wavef_init

   subroutine t_balance_wavef_balance(bal)
      implicit none
      class(t_balance_wavef), intent(inout) :: bal 

      integer :: i, ierr

#ifdef CPP_MPI
      do i = 1,bal%remaining_barries
         ! call timestart("Post read_z Barrier: is_fft(dangle)")
         ! call MPI_Barrier(MPI_COMM_WORLD, ierr)
         ! call timestop("Post read_z Barrier: is_fft(dangle)")
      enddo
#endif
   end subroutine t_balance_wavef_balance
end module m_balance_barriers
