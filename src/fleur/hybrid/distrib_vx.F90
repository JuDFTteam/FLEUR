module m_distrib_vx
   use m_types
#ifdef CPP_MPI
   use mpi
#endif
   use m_types_mpimat
   use m_glob_tofrom_loc
contains
   subroutine distrib_vx(fi, fmpi, nococonv, vx_loc, vx_tmp, hybdat)
      implicit none 
      type(t_fleurinput), intent(in)    :: fi
      type(t_nococonv), intent(in)      :: nococonv
      TYPE(t_mpi), INTENT(IN)           :: fmpi
      integer, intent(in)               :: vx_loc(:,:)
      type(t_mat), intent(inout)        :: vx_tmp(:,:)
      type(t_hybdat), intent(inout)     :: hybdat


      integer            :: jsp, nk

      call timestart("distrib_vx")

      if(.not. allocated(hybdat%v_x)) then
         IF (fmpi%n_size == 1) THEN
            ALLOCATE (t_mat::hybdat%v_x(fi%kpts%nkpt, fi%input%jspins))
         ELSE
            ALLOCATE (t_mpimat::hybdat%v_x(fi%kpts%nkpt, fi%input%jspins))
         END IF
      endif

      do jsp = 1,fi%input%jspins
         do nk = 1,fi%kpts%nkpt
            call distrib_single_vx(fi, fmpi, jsp, nk, vx_loc(nk,jsp), vx_tmp(nk,jsp), hybdat, nococonv=nococonv)
         enddo
      enddo

      call timestop("distrib_vx")
   end subroutine distrib_vx

   subroutine distrib_single_vx(fi, fmpi, jsp, nk, vx_loc, vx_tmp, hybdat, dims, nococonv)
      implicit none 
      type(t_fleurinput), intent(in)    :: fi
      TYPE(t_mpi), INTENT(IN)           :: fmpi
      integer, intent(in)               :: vx_loc, jsp, nk
      type(t_mat), intent(inout)        :: vx_tmp
      type(t_hybdat), intent(inout)     :: hybdat
      integer, optional, intent(in)     :: dims(2)
      type(t_nococonv), intent(in), optional      :: nococonv
   
      type(t_lapw)       :: lapw
      integer            :: mat_sz, recver, ierr

      call timestart("distrib_vx_single")

      if(present(nococonv) .and. (.not. present(dims))) then
         CALL lapw%init(fi, nococonv, nk) 
         mat_sz = lapw%nv(jsp) + fi%atoms%nlotot
      elseif((.not. present(nococonv)) .and. present(dims)) then
         mat_sz = dims(1)
      else
         call juDFT_error("I need either nococonv or dims")
      endif

      if(any(nk == fmpi%k_list)) then 
         CALL hybdat%v_x(nk, jsp)%init(fi%sym%invs, mat_sz, mat_sz, fmpi%sub_comm, .false.)
      else
         call hybdat%v_x(nk, jsp)%init(fi%sym%invs, 1,1, fmpi%sub_comm, .false.)
      endif

      call copy_vx_to_distr(fmpi, vx_loc, nk, mat_sz, vx_tmp, hybdat%v_x(nk, jsp))
      call vx_tmp%free()

      call timestop("distrib_vx_single")
   end subroutine distrib_single_vx

   subroutine copy_vx_to_distr(fmpi, sender, nk, mat_sz, vx_den, vx_distr)
      implicit none
      type(t_mpi), intent(in)       :: fmpi
      integer, intent(in)           :: sender, nk, mat_sz
      type(t_mat), intent(in)       :: vx_den 
      class(t_mat), intent(inout) :: vx_distr

      integer :: i_loc, pe_i, i, ierr, recver

      pe_i = 0
      i_loc = 1
      do i = 1, mat_sz
         call glob_to_loc(fmpi, i, pe_i, i_loc)
         if(any(fmpi%k_list == nk) .and. pe_i == fmpi%n_rank) then
            recver = fmpi%irank 
         else
            recver = -1
         endif
#ifdef CPP_MPI
         call MPI_Allreduce(MPI_IN_PLACE, recver, 1, MPI_INTEGER, MPI_MAX, fmpi%mpi_comm, ierr)
#endif

         if(fmpi%irank == sender .and. sender == recver) then
            if(vx_den%l_real) then
               vx_distr%data_r(:,i_loc) = vx_den%data_r(:,i)
            else
               vx_distr%data_c(:,i_loc) = vx_den%data_c(:,i)
            endif
#ifdef CPP_MPI
         elseif(fmpi%irank == sender) then
            if(vx_den%l_real) then
               call MPI_Send(vx_den%data_r(:,i), vx_den%matsize1, MPI_DOUBLE_PRECISION, recver, i, fmpi%mpi_comm, ierr)
            else
               call MPI_Send(vx_den%data_c(:,i), vx_den%matsize1, MPI_DOUBLE_COMPLEX, recver, i, fmpi%mpi_comm, ierr)
            endif 
         elseif(fmpi%irank == recver) then 
            if(vx_distr%l_real) then
               call MPI_Recv(vx_distr%data_r(:,i_loc), vx_distr%matsize1, MPI_DOUBLE_PRECISION, sender, i, fmpi%mpi_comm, MPI_STATUS_IGNORE, ierr)
            else
               call MPI_Recv(vx_distr%data_c(:,i_loc), vx_distr%matsize1, MPI_DOUBLE_COMPLEX, sender, i, fmpi%mpi_comm, MPI_STATUS_IGNORE, ierr)
            endif
#endif 
         endif
      enddo
   end subroutine copy_vx_to_distr
end module m_distrib_vx