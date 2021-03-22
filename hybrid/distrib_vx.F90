module m_distrib_vx
contains
   subroutine distrib_vx(fi, fmpi, nococonv, vx_loc, vx_tmp, hybdat)
      use m_types
      use mpi
      use m_types_mpimat
      implicit none 
      type(t_fleurinput), intent(in)    :: fi
      type(t_nococonv), intent(in)      :: nococonv
      TYPE(t_mpi), INTENT(IN)           :: fmpi
      integer, intent(in)               :: vx_loc(:,:)
      type(t_mat), intent(inout)        :: vx_tmp(:,:)
      type(t_hybdat), intent(inout)     :: hybdat


      type(t_lapw)       :: lapw
      integer            :: jsp, ik , recver, sender, ierr
      integer            :: matinfo(3)

      call timestart("distrib_vx")

      if(.not. allocated(hybdat%v_x)) then
         IF (fmpi%n_size == 1) THEN
            ALLOCATE (t_mat::hybdat%v_x(fi%kpts%nkpt, fi%input%jspins))
         ELSE
            ALLOCATE (t_mpimat::hybdat%v_x(fi%kpts%nkpt, fi%input%jspins))
         END IF
      endif

      do jsp = 1,fi%input%jspins
         do ik = 1,fi%kpts%nkpt
            CALL lapw%init(fi, nococonv, ik) 
            ! make sure vx_tmp is on n_rank = 0
            call timestart("move vx_tmp to n_rank 0")
            sender = vx_loc(ik, jsp)
            if(fmpi%n_rank == 0 .and. any(ik == fmpi%k_list)) then
               recver = fmpi%irank 
            else 
               recver = -1 
            endif 
            write (*,*)
#ifdef CPP_MPI
            call MPI_Allreduce(MPI_IN_PLACE, recver, 1, MPI_INTEGER, MPI_MAX, MPI_COMM_WORLD, ierr)

            if(sender /= recver) then
               if(fmpi%irank == sender) then
                  matinfo(1) = merge(1, 0, vx_tmp(ik,jsp)%l_real)
                  matinfo(2) = vx_tmp(ik,jsp)%matsize1
                  matinfo(3) = vx_tmp(ik,jsp)%matsize2
                  call MPI_Send(matinfo, 3, MPI_INTEGER, recver, 7, MPI_COMM_WORLD, ierr)
                  if(vx_tmp(ik,jsp)%l_real) then
                     call MPI_Send(vx_tmp(ik,jsp)%data_r, size(vx_tmp(ik,jsp)%data_r), MPI_DOUBLE_PRECISION, recver, 1000*recver+sender, MPI_COMM_WORLD, ierr)
                  else 
                     call MPI_Send(vx_tmp(ik,jsp)%data_c, size(vx_tmp(ik,jsp)%data_c), MPI_DOUBLE_COMPLEX, recver, 1000*recver+sender, MPI_COMM_WORLD, ierr)
                  endif
                  call vx_tmp(ik,jsp)%free()
               elseif(fmpi%irank == recver) then 
                  call MPI_Recv(matinfo, 3, MPI_INTEGER, sender, 7, MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr) 
                  call vx_tmp(ik,jsp)%init(matinfo(1) == 1, matinfo(2), matinfo(3))
                  if(vx_tmp(ik,jsp)%l_real) then 
                     call MPI_Recv(vx_tmp(ik,jsp)%data_r, size(vx_tmp(ik,jsp)%data_r), MPI_DOUBLE_PRECISION, sender, 1000*recver+sender, MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
                  else 
                     call MPI_Recv(vx_tmp(ik,jsp)%data_c, size(vx_tmp(ik,jsp)%data_c), MPI_DOUBLE_COMPLEX, sender, 1000*recver+sender, MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
                  endif 
               endif 
            endif
            call timestop("move vx_tmp to n_rank 0")
#endif

            if(any(ik == fmpi%k_list)) then 
               if(fmpi%n_size == 1) then 
                  select type(vx_mat => hybdat%v_x(ik, jsp))
                  class is(t_mat)
                     call hybdat%v_x(ik, jsp)%copy(vx_tmp(ik,jsp), 1, 1)
                  class default
                     call juDFT_error("hybdat%vx needs to be t_mat")
                  end select
               else 
                  select type(vx_mat => hybdat%v_x(ik, jsp))
                  class is(t_mpimat)
                     CALL vx_mat%init(fi%sym%invs, lapw%nv(jsp) + fi%atoms%nlotot, &
                                                             lapw%nv(jsp) + fi%atoms%nlotot, fmpi%sub_comm, .false.)

                     !prep vx_tmp for from_non_dist
                     call vx_tmp(ik,jsp)%init(fi%sym%invs, 1,1)
                     call MPI_Bcast(vx_tmp%matsize1, 1, MPI_INTEGER, 0, fmpi%sub_comm, ierr)
                     call MPI_Bcast(vx_tmp%matsize2, 1, MPI_INTEGER, 0, fmpi%sub_comm, ierr)

                     call vx_mat%from_non_dist(vx_tmp(ik,jsp))
                  class default
                     call juDFT_error("hybdat%vx needs to be t_mpimat")
                  end select
               endif
            else 
               select type(vx_mat => hybdat%v_x(ik, jsp))
               class is(t_mat)
                  call hybdat%v_x(ik, jsp)%init(fi%sym%invs, 1,1)
               class is(t_mpimat)
                  call hybdat%v_x(ik, jsp)%init(fi%sym%invs, 1,1, fmpi%sub_comm, .false.)
               end select
            endif

            call vx_tmp(ik,jsp)%free()

            call timestart("barrier")
            call MPI_Barrier(MPI_COMM_WORLD, ierr)
            call timestop("barrier")
         enddo
      enddo

      call timestop("distrib_vx")
   end subroutine distrib_vx
end module m_distrib_vx