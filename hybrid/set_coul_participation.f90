module m_set_coul_participation
   use m_types
contains
   subroutine set_coul_participation(hybdat, fi, fmpi, work_pack)
      use m_work_package
      implicit none
      class(t_hybdat), intent(inout)   :: hybdat
      type(t_fleurinput), intent(in)   :: fi
      TYPE(t_mpi), INTENT(IN)          :: fmpi
      type(t_work_package), intent(in) :: work_pack(:)

      integer :: jsp, ik, jq, iq, iq_p

      do ik = 1,fi%kpts%nkpt
         hybdat%coul(ik)%l_participate = .false.
      enddo

      do jsp = 1,size(work_pack)
         do ik = 1,work_pack(jsp)%k_packs(1)%size
            do jq = 1, size(work_pack(jsp)%k_packs(ik)%q_packs)
               iq = work_pack(jsp)%k_packs(ik)%q_packs(jq)%ptr
               iq_p = fi%kpts%bkp(iq)

               hybdat%coul(iq_p)%l_participate = .True.
            enddo
         enddo
      enddo

      do ik = 1,fi%kpts%nkpt
         if(fmpi%n_rank == 0 .and. any(fmpi%k_list == ik)) then 
            hybdat%coul(ik)%l_participate = .true. 
         endif 
      enddo
   end subroutine set_coul_participation
end module m_set_coul_participation
