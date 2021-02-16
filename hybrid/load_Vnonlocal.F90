module m_load_Vnonlocal
   use m_types
   use m_io_matrix
contains
   subroutine load_Vnonlocal(fi, fmpi, mpdata, hybdat)
      implicit none
      type(t_fleurinput), intent(in) :: fi
      type(t_mpi), intent(in)        :: fmpi
      type(t_mpdata), intent(in)     :: mpdata
      TYPE(t_hybdat), INTENT(inout)  :: hybdat

      integer                       :: fid, mat_sz, nk, nk_i, jsp, no_records,  record
      character(len=:), allocatable :: filename
      logical                       :: files_present, l_tmp

#ifdef CPP_HDF
      ! IF (fi%hybinp%l_hybrid .and. (.not. hybdat%l_addhf)) then
      !    call timestart("load_Vnonlocal")

      !    INQUIRE (file="nbands.hdf", exist=files_present)
      !    INQUIRE (file="v_x.hdf", exist=l_tmp)
      !    files_present = files_present .and. l_tmp

      !    if (files_present) then
      !       call hybdat%read_nbands(fi)
      !       if (.not. allocated(hybdat%v_x)) allocate (hybdat%v_x(fi%kpts%nkpt, fi%input%jspins))

      !       no_records = fi%kpts%nkpt*fi%input%jspins

      !       fid = open_matrix(fi%sym%invs, -1, 2, no_records, "v_x")

      !       DO jsp = 1, fi%input%jspins
      !          DO nk_i = 1, size(fmpi%k_list)
      !             nk = fmpi%k_list(nk_i)
      !             record = (jsp -1) * fi%kpts%nkpt + nk
      !             call read_matrix(hybdat%v_x(nk, jsp), record, fid)
      !          end do
      !       end do
      !       call close_matrix(fid)

      !       ! prep for add_Vnonlocal and subvxc
      !       hybdat%l_addhf  = .True.
      !       hybdat%l_subvxc = .True. 
            
      !       call mpdata%set_num_radfun_per_l(fi%atoms)
      !       call hybdat%set_maxlmindx(fi%atoms, mpdata%num_radfun_per_l)

      !    end if
      !    call timestop("load_Vnonlocal")
      ! end if
#endif
   end subroutine load_Vnonlocal

end module m_load_Vnonlocal
