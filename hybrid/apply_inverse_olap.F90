module m_apply_inverse_olap
   use m_glob_tofrom_loc
contains
   subroutine apply_inverse_olaps(mpdata, atoms, cell, hybdat, fmpi, sym, ikpt, coulomb)
      USE m_olap, ONLY: olap_pw
      USE m_types
      use m_judft
      implicit none
      type(t_mpdata), intent(in) :: mpdata
      type(t_atoms), intent(in)  :: atoms
      type(t_cell), intent(in)   :: cell
      type(t_hybdat), intent(in) :: hybdat
      type(t_mpi), intent(in)    :: fmpi
      type(t_sym), intent(in)    :: sym
      type(t_mat), intent(inout) :: coulomb
      integer, intent(in)        :: ikpt

      type(t_mat)     :: olap, coul_submtx
      integer         :: nbasm, loc_size, i, j, i_loc, ierr, pe_i, pe_j, pe_recv, pe_send, recv_loc, send_loc, j_loc
      complex         :: cdum

      call timestart("solve olap linear eq. sys")
      nbasm = hybdat%nbasp + mpdata%n_g(ikpt)
      CALL olap%alloc(sym%invs, mpdata%n_g(ikpt), mpdata%n_g(ikpt), 0.0)
      !calculate IR overlap-matrix
      CALL olap_pw(olap, mpdata%g(:, mpdata%gptm_ptr(:mpdata%n_g(ikpt), ikpt)), mpdata%n_g(ikpt), atoms, cell, fmpi)

      ! perform O^-1 * coulhlp%data_r(hybdat%nbasp + 1:, :) = x
      ! rewritten as O * x = C

      loc_size = 0
      do i = 1, nbasm
         call glob_to_loc(fmpi, i, pe_i, i_loc)
         if (fmpi%n_rank == pe_i) loc_size = loc_size + 1
      end do

      call timestart("copy in 1")
      call coul_submtx%alloc(sym%invs, mpdata%n_g(ikpt), loc_size)
      if (coul_submtx%l_real) then
         coul_submtx%data_r(:, :) = real(coulomb%data_c(hybdat%nbasp + 1:, :))
      else
         coul_submtx%data_c(:, :) = coulomb%data_c(hybdat%nbasp + 1:, :)
      end if
      call timestop("copy in 1")

      call olap%linear_problem(coul_submtx)
      call timestart("copy out 1")
      if (coul_submtx%l_real) then
         coulomb%data_c(hybdat%nbasp + 1:, :) = coul_submtx%data_r
      else
         coulomb%data_c(hybdat%nbasp + 1:, :) = coul_submtx%data_c
      end if
      call timestop("copy out 1")

      call timestart("copy in 2")
      do j = 1, mpdata%n_g(ikpt)
         call glob_to_loc(fmpi, hybdat%nbasp + j, pe_j, j_loc)
         do i = hybdat%nbasp + 1, nbasm
            call glob_to_loc(fmpi, i, pe_i, i_loc)
            if (pe_j == pe_i .and. fmpi%n_rank == pe_i) then
               if (coul_submtx%l_real) then
                  coul_submtx%data_r(j, i_loc) = real(coulomb%data_c(i, j_loc))
               else
                  coul_submtx%data_c(j, i_loc) = conjg(coulomb%data_c(i, j_loc))
               end if
#ifdef CPP_MPI
            elseif (pe_j == fmpi%n_rank) then
               call MPI_Send(coulomb%data_c(i, j_loc), 1, MPI_DOUBLE_COMPLEX, pe_i, j + 10000*i, fmpi%sub_comm, ierr)
            elseif (pe_i == fmpi%n_rank) then
               call MPI_Recv(cdum, 1, MPI_DOUBLE_COMPLEX, pe_j, j + 10000*i, fmpi%sub_comm, MPI_STATUS_IGNORE, ierr)
               if (coul_submtx%l_real) then
                  coul_submtx%data_r(j, i_loc) = real(cdum)
               else
                  coul_submtx%data_c(j, i_loc) = conjg(cdum)
               end if
#endif
            end if
         end do
      end do
      call timestop("copy in 2")

      ! perform  coulomb%data_r(hybdat%nbasp + 1:, :) * O^-1  = X
      ! rewritten as O^T * x^T = C^T

      ! reload O, since the solver destroys it.
      CALL olap_pw(olap, mpdata%g(:, mpdata%gptm_ptr(:mpdata%n_g(ikpt), ikpt)), mpdata%n_g(ikpt), atoms, cell, fmpi)
      ! Notice O = O^T since it's symmetric
      call olap%linear_problem(coul_submtx)

      call timestart("copy out 2")
      do j = 1, mpdata%n_g(ikpt)
         call glob_to_loc(fmpi, hybdat%nbasp + j, pe_recv, recv_loc)
         do i = hybdat%nbasp + 1, nbasm
            call glob_to_loc(fmpi, i, pe_send, send_loc)
            if (pe_send == pe_recv .and. fmpi%n_rank == pe_recv) then
               if (coul_submtx%l_real) then
                  coulomb%data_c(i, recv_loc) = coul_submtx%data_r(j, send_loc)
               else
                  coulomb%data_c(i, recv_loc) = conjg(coul_submtx%data_c(j, send_loc))
               end if
#ifdef CPP_MPI
            elseif (pe_send == fmpi%n_rank) then
               if (coul_submtx%l_real) then
                  cdum = coul_submtx%data_r(j, send_loc)
               else
                  cdum = conjg(coul_submtx%data_c(j, send_loc))
               end if
               call MPI_Send(cdum, 1, MPI_DOUBLE_COMPLEX, pe_recv, j + 10000*i, fmpi%sub_comm, ierr)
            elseif (pe_recv == fmpi%n_rank) then
               call MPI_Recv(coulomb%data_c(i, recv_loc), 1, MPI_DOUBLE_COMPLEX, pe_send, j + 10000*i, fmpi%sub_comm, MPI_STATUS_IGNORE, ierr)
#endif
            end if
         end do
      end do
      call timestop("copy out 2")

      call coul_submtx%free()
      call olap%free()
      call timestop("solve olap linear eq. sys")
   end subroutine apply_inverse_olaps
end module m_apply_inverse_olap
