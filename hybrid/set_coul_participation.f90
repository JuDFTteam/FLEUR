module m_set_coul_participation
   use m_types
contains
   subroutine set_coul_participation(hybdat, fi, work_pack)
      use m_work_package
      implicit none
      class(t_hybdat), intent(inout) :: hybdat
      type(t_fleurinput), intent(in)   :: fi
      type(t_work_package), intent(in) :: work_pack(:)

      integer :: i, jsp, ik, jq, iq, iq_p

      ! do i = 1,fi%kpts%nkpt
      !    hybdat%coul(i)%l_participate = .false.
      ! enddo

      do jsp = 1,size(work_pack)
         do ik = 1,work_pack(jsp)%k_packs(1)%size
            do jq = 1, size(work_pack(jsp)%k_packs(ik)%q_packs)
               iq = work_pack(jsp)%k_packs(ik)%q_packs(jq)%ptr
               iq_p = fi%kpts%bkp(iq)

               hybdat%coul(iq_p)%l_participate = .True.
            enddo
         enddo
      enddo

   end subroutine set_coul_participation
end module m_set_coul_participation
