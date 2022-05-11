module m_eigvec_setup
   use m_types
   use m_work_package
   implicit none

contains
   subroutine eigvec_setup(eigvec, fi, lapw, work_packs, fmpi, nbands, ik, jsp, eig_id)
      implicit none
      class(t_eigvec), intent(inout)   :: eigvec
      type(t_fleurinput), intent(in)   :: fi
      TYPE(t_lapw), INTENT(IN)         :: lapw
      integer, intent(in)              :: ik, jsp, eig_id
      type(t_work_package), intent(in) :: work_packs(:)
      type(t_mpi), intent(in)          :: fmpi
      integer, intent(in)              :: nbands! hybdat%nbands(ik,jsp) passed like this to avoid circular dependencies

      integer :: nbasfcn

      eigvec%nk  = ik
      eigvec%jsp = jsp

      call eigvec_set_part_and_band(eigvec, fi, work_packs, fmpi, nbands, jsp)
      !communication only happen on reduced BZ
      if (ik <= fi%kpts%nkpt) call eigvec_create_comm(eigvec, fi, eig_id, ik, jsp, nbands)

      if (eigvec%l_recv) then
         nbasfcn = lapw%hyb_num_bas_fun(fi)
         call eigvec%mat%alloc(fi%sym%invs, nbasfcn, nbands)
      endif
   end subroutine eigvec_setup

   subroutine eigvec_set_part_and_band(eigvec, fi, work_packs, fmpi, nbands, jsp)
      implicit none
      class(t_eigvec), intent(inout)   :: eigvec
      type(t_fleurinput), intent(in)   :: fi
      type(t_mpi), intent(in)          :: fmpi
      type(t_work_package), intent(in) :: work_packs(:)
      integer, intent(in)              :: nbands ! hybdat%nbands(nk,jsp) passed like this to avoid circular dependencies
      integer, intent(in)              :: jsp

      integer :: i

      !set senders
      eigvec%l_participate = any(fmpi%k_list == eigvec%nk)

      !set recipients for k-side
      do i = 1, work_packs(jsp)%k_packs(1)%size
         if (eigvec%nk == work_packs(jsp)%k_packs(i)%nk) then
            eigvec%l_participate = .True.
            eigvec%l_recv = .True.            
         endif
      enddo
   end subroutine eigvec_set_part_and_band

   subroutine bcast_eigvecs(hybdat, fi, nococonv, fmpi)
      use m_eig66_data
      USE m_eig66_io
      use m_eig66_mpi, only: priv_find_data
      use m_judft
      use m_io_hybrid
      implicit none
      type(t_hybdat), intent(inout)     :: hybdat
      type(t_fleurinput), intent(in)    :: fi
      TYPE(t_nococonv), INTENT(IN)      :: nococonv
      TYPE(t_mpi), INTENT(IN)           :: fmpi

      type(t_lapw)              :: lapw
      type(t_mat) :: tmp


      logical  :: l_zref
      integer  :: jsp, ik, nbasfcn, ieig, ierr, root, me

      call timestart("bcast zmat")

      l_zref = (fi%sym%zrfs .AND. (SUM(ABS(fi%kpts%bk(3, :fi%kpts%nkpt))) < 1e-9) .AND. .NOT. fi%noco%l_noco)
      select case (eig66_data_mode(hybdat%eig_id) )
      case( mpi_mode)
#ifdef CPP_MPI
         do jsp = 1, fi%input%jspins
            do ik = 1, fi%kpts%nkpt
               if(hybdat%zmat(ik, jsp)%l_participate) then

                  CALL lapw%init(fi%input, fi%noco, nococonv, fi%kpts, fi%atoms, fi%sym, ik, fi%cell)
                  !allocate tmp array
                  nbasfcn = lapw%hyb_num_bas_fun(fi)
                  call tmp%alloc(fi%sym%invs, nbasfcn, 1)
                  do ieig = 1, hybdat%nbands(ik,jsp)
                     root = hybdat%zmat(ik, jsp)%root_pe(ieig)
                     call MPI_comm_rank(hybdat%zmat(ik, jsp)%comm, me, ierr)
                     ! make sure read_eig is only run if I have it in mem
                     if (me == root) call read_eig(hybdat%eig_id, ik, jsp, zmat=tmp, list=[ieig])

                     if (fi%sym%invs) then
                        call MPI_Bcast(tmp%data_r, nbasfcn, MPI_DOUBLE_PRECISION, root, hybdat%zmat(ik, jsp)%comm, ierr)
                     else
                        call MPI_Bcast(tmp%data_c, nbasfcn, MPI_DOUBLE_COMPLEX, root, hybdat%zmat(ik, jsp)%comm, ierr)
                     endif
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
#endif
      case(mem_mode)
         do jsp = 1, fi%input%jspins
            do ik = 1, fi%kpts%nkpt
               call read_z(fi%atoms, fi%cell, hybdat, fi%kpts, fi%sym, fi%noco, nococonv, fi%input, ik, jsp, hybdat%zmat(ik,jsp)%mat)
            enddo 
         enddo
      CASE DEFAULT
         CALL juDFT_error("The hybrid-code only supports eigvec comm via MEM or MPI")
      END select

      call timestop("bcast zmat")
   end subroutine bcast_eigvecs

   subroutine eigvec_create_comm(eigvec, fi, eig_id, ik, jsp, nbands)
      use m_types_mpi
      use m_types_lapw
      use m_eig66_data
      use m_eig66_io
      use m_eig66_mpi, only: priv_find_data
      implicit none
      class(t_eigvec), intent(inout)     :: eigvec
      type(t_fleurinput), intent(in)     :: fi
      integer, intent(in)                :: eig_id, ik, jsp, nbands
#ifdef CPP_MPI

      TYPE(t_data_MPI), POINTER, ASYNCHRONOUS :: d
      integer :: color, me_glob, me_loc,  ieig, ierr

      select case (eig66_data_mode(eig_id) )
      case( mpi_mode)
         CALL priv_find_data(eig_id, d)
         
         if(eigvec%comm == MPI_COMM_NULL) then
            color = merge(1,2,eigvec%l_participate)
            call judft_comm_split(MPI_COMM_WORLD, color, 1, eigvec%comm)
         endif


         if(eigvec%l_participate) then
            call mpi_comm_rank(MPI_COMM_WORLD, me_glob, ierr)
            call mpi_comm_rank(eigvec%comm, me_loc, ierr)
         
            if(allocated(eigvec%root_pe)) deallocate(eigvec%root_pe)
            allocate(eigvec%root_pe(nbands), source=-1)

            do ieig = 1,nbands
               if(me_glob == d%pe_ev(ik, jsp, ieig)) then 
                  eigvec%root_pe(ieig) = me_loc 
               endif 
            enddo
            call MPI_Allreduce(MPI_IN_PLACE, eigvec%root_pe, nbands, MPI_INTEGER, MPI_MAX, eigvec%comm, ierr)

            if(any(eigvec%root_pe < 0)) call judft_error("A vector can't be on a negative PE. Distrb failed.")
         endif
      case(mem_mode)
         
      CASE DEFAULT
         CALL juDFT_error("The hybrid-code only supports eigvec comm via MEM or MPI")
      END select
#endif
   end subroutine eigvec_create_comm
end module m_eigvec_setup
