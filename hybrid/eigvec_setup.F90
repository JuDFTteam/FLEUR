module m_eigvec_setup
   use m_types
   use m_work_package
   implicit none

contains
   subroutine eigvec_setup(eigvec, fi, lapw, work_packs, fmpi, nbands, ik)
      implicit none
      class(t_eigvec), intent(inout)   :: eigvec
      type(t_fleurinput), intent(in)   :: fi
      TYPE(t_lapw), INTENT(IN)         :: lapw
      integer, intent(in)              :: ik
      type(t_work_package), intent(in) :: work_packs(:)
      type(t_mpi), intent(in)          :: fmpi
      integer, intent(in)              :: nbands! hybdat%nbands(ik,jsp) passed like this to avoid circular dependencies

      integer :: nbasfcn

      eigvec%nk = ik

      call eigvec_set_part_and_band(eigvec, fi, work_packs, fmpi, nbands)
      !communication only happen on reduced BZ
      if (ik <= fi%kpts%nkpt) call eigvec%create_comm()

      if (eigvec%l_recv) then
         nbasfcn = lapw%hyb_num_bas_fun(fi)
         call eigvec%mat%alloc(fi%sym%invs, nbasfcn, (eigvec%max_band - eigvec%min_band + 1))
      endif
   end subroutine eigvec_setup

   subroutine eigvec_set_part_and_band(eigvec, fi, work_packs, fmpi, nbands)
      implicit none
      class(t_eigvec), intent(inout)   :: eigvec
      type(t_fleurinput), intent(in)   :: fi
      type(t_mpi), intent(in)          :: fmpi
      type(t_work_package), intent(in) :: work_packs(:)
      integer, intent(in)              :: nbands ! hybdat%nbands(nk,jsp) passed like this to avoid circular dependencies

      integer :: jsp, i

      !set senders
      eigvec%l_participate = any(fmpi%k_list == eigvec%nk)

      !set recipients for k-side
      do jsp = 1, fi%input%jspins
         do i = 1, work_packs(jsp)%k_packs(1)%size
            if (eigvec%nk == work_packs(jsp)%k_packs(i)%nk) then
               eigvec%l_participate = .True.
               eigvec%l_recv = .True.

               eigvec%min_band = min(eigvec%min_band, 1)
               eigvec%max_band = max(eigvec%max_band, nbands)
            endif
         enddo
      enddo
   end subroutine eigvec_set_part_and_band

   subroutine bcast_eigvecs(hybdat, fi, nococonv, fmpi)
      use m_eig66_data
      USE m_eig66_io
      use m_eig66_mpi, only: priv_find_data
      implicit none
      type(t_hybdat), intent(inout)     :: hybdat
      type(t_fleurinput), intent(in)    :: fi
      TYPE(t_nococonv), INTENT(IN)      :: nococonv
      TYPE(t_mpi), INTENT(IN)           :: fmpi

      TYPE(t_data_MPI), POINTER, ASYNCHRONOUS :: d
      type(t_lapw)              :: lapw
      type(t_mat) :: tmp


      logical  :: l_zref
      integer  :: jsp, ik, nbasfcn, ieig, ierr, root

      l_zref = (fi%sym%zrfs .AND. (SUM(ABS(fi%kpts%bk(3, :fi%kpts%nkpt))) < 1e-9) .AND. .NOT. fi%noco%l_noco)
      select case (eig66_data_mode(hybdat%eig_id) )
      case( mpi_mode)
         CALL priv_find_data(hybdat%eig_id, d)

         do jsp = 1, fi%input%jspins
            do ik = 1, fi%kpts%nkpt
               if(hybdat%zmat(ik, jsp)%l_participate) then

                  CALL lapw%init(fi%input, fi%noco, nococonv, fi%kpts, fi%atoms, fi%sym, ik, fi%cell, l_zref)
                  !allocate tmp array
                  nbasfcn = lapw%hyb_num_bas_fun(fi)
                  call tmp%alloc(fi%sym%invs, nbasfcn, 1)
                  write (*,*) "nbands=", hybdat%nbands(ik,jsp)
                  do ieig = 1, hybdat%nbands(ik,jsp)
                     root = d%pe_ev(ik, jsp, ieig)
                     ! make sure read_eig is only run if I have it in mem
                     if (fmpi%irank == root) call read_eig(hybdat%eig_id, ik, jsp, zmat=tmp, list=[ieig])

#ifdef CPP_MPI
                     if (fi%sym%invs) then
                        call MPI_Bcast(tmp%data_r, nbasfcn, MPI_DOUBLE_PRECISION, root, hybdat%zmat(ik, jsp)%comm, ierr)
                     else
                        call MPI_Bcast(tmp%data_c, nbasfcn, MPI_DOUBLE_COMPLEX, root, hybdat%zmat(ik, jsp)%comm, ierr)
                     endif
#endif
                     ! deal with k-copies
                     if(hybdat%zmat(ik, jsp)%l_recv) then 
                        if(fi%sym%invs)then
                           hybdat%zmat(ik, jsp)%mat%data_r(:,ieig) = tmp%data_r(:,1)
                        else 
                           hybdat%zmat(ik, jsp)%mat%data_c(:,ieig) = tmp%data_c(:,1)
                        endif
                     endif
                  enddo
                  call tmp%free()

               endif
            enddo
         enddo
      case(mem_mode)
         call juDFT_error("need to implement this!")
      CASE DEFAULT
         CALL juDFT_error("The hybrid-code only supports eigvec comm via MEM or MPI")
      END select

   end subroutine bcast_eigvecs
end module m_eigvec_setup
