module m_apply_inverse_olap
   use m_glob_tofrom_loc
   USE m_types_mpimat
contains
   subroutine apply_inverse_olaps(mpdata, atoms, cell, hybdat, fmpi, sym, ikpt, coulomb)
      USE m_olap, ONLY: olap_pw
      USE m_types
      use m_judft
      implicit none
      type(t_mpdata), intent(in)  :: mpdata
      type(t_atoms), intent(in)   :: atoms
      type(t_cell), intent(in)    :: cell
      type(t_hybdat), intent(in)  :: hybdat
      type(t_mpi), intent(in)     :: fmpi
      type(t_sym), intent(in)     :: sym
      class(t_mat), intent(inout) :: coulomb
      integer, intent(in)         :: ikpt

      type(t_mat)               :: olap
      class(t_mat), allocatable :: coul_submtx
      type(t_mpimat)            :: olap_mpi

      integer         :: nbasm, loc_size, i, j, i_loc, ierr, pe_i, pe_j, pe_recv, pe_send, recv_loc, send_loc, j_loc
      complex         :: cdum

      call timestart("solve olap linear eq. sys")
      nbasm = hybdat%nbasp + mpdata%n_g(ikpt)
      CALL olap%alloc(.false., mpdata%n_g(ikpt), mpdata%n_g(ikpt), 0.0)
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
      allocate(t_mat::coul_submtx)
      call coul_submtx%alloc(.false., mpdata%n_g(ikpt), loc_size)
      coul_submtx%data_c(:, :) = coulomb%data_c(hybdat%nbasp + 1:, :)
      call timestop("copy in 1")

      !$acc data copyin(olap, olap%data_r, olap%data_c, coul_submtx) copy(coul_submtx%data_r, coul_submtx%data_c)
         call olap%linear_problem(coul_submtx)
      !$acc end data
      call timestart("copy out 1")
      coulomb%data_c(hybdat%nbasp + 1:, :) = coul_submtx%data_c
      call coul_submtx%free()
      deallocate(coul_submtx)
      call timestop("copy out 1")


      ! perform  coulomb%data_r(hybdat%nbasp + 1:, :) * O^-1  = X
      ! rewritten as O^T * x^T = C^T
      call copy_in_2(fmpi, sym, mpdata, hybdat, coulomb, ikpt, coul_submtx)

      ! reload O, since the solver destroys it.
      CALL olap_pw(olap, mpdata%g(:, mpdata%gptm_ptr(:mpdata%n_g(ikpt), ikpt)), mpdata%n_g(ikpt), atoms, cell, fmpi)
      ! Notice O = O^T since it's symmetric

      SELECT TYPE(coul_submtx)
      CLASS is (t_mat)
         !$acc data copyin(olap, olap%data_r, olap%data_c, coul_submtx) copy(coul_submtx%data_r, coul_submtx%data_c)
            call olap%linear_problem(coul_submtx)
         !$acc end data
         call olap%free()
      class is (t_mpimat)
         call olap_mpi%init(coul_submtx,  olap%matsize1, olap%matsize2)
         call olap_mpi%from_non_dist(olap)
         call olap_mpi%linear_problem(coul_submtx)
         call olap_mpi%free()
      end select

      call copy_out_2(fmpi, sym, mpdata, hybdat, ikpt, coul_submtx, coulomb)
      deallocate(coul_submtx)
      call timestop("solve olap linear eq. sys")
   end subroutine apply_inverse_olaps

   subroutine copy_in_2(fmpi, sym, mpdata, hybdat, coulomb, ikpt, coul_submtx)
      USE m_types
      implicit none 
      type(t_mpi), intent(in)      :: fmpi 
      integer, intent(in)          :: ikpt
      type(t_sym), intent(in)      :: sym
      type(t_mpdata), intent(in)   :: mpdata 
      type(t_hybdat), intent(in)   :: hybdat
      class(t_mat), intent(in)     :: coulomb 
      class(t_mat), intent(inout), allocatable  :: coul_submtx

      integer :: i, j, ierr, i_loc, j_loc, pe_i, pe_j
      complex :: cdum

      call timestart("copy in 2")

      SELECT TYPE(coulomb)
      CLASS is (t_mat)
         allocate(t_mat::coul_submtx)
         call coul_submtx%alloc(.false., mpdata%n_g(ikpt), mpdata%n_g(ikpt))
         do j = 1, mpdata%n_g(ikpt)
            do i = 1, mpdata%n_g(ikpt)
               coul_submtx%data_c(j, i) = conjg(coulomb%data_c(hybdat%nbasp+i, hybdat%nbasp + j))
            enddo 
         enddo
      class is (t_mpimat)
#ifdef CPP_SCALAPACK
         allocate(t_mpimat::coul_submtx)
         call coul_submtx%init(.False., mpdata%n_g(ikpt), mpdata%n_g(ikpt), fmpi%sub_comm, .True.)
         select type(coul_submtx)
         class is (t_mpimat)
            ! copy bottom right corner of coulomb to coul_submtx
            !call pzgemr2d(m,              n,               a,                ia,           ja,             desca, 
            call pzgemr2d(mpdata%n_g(ikpt),mpdata%n_g(ikpt),coulomb%data_c, hybdat%nbasp+1, hybdat%nbasp+1, coulomb%blacsdata%blacs_desc,&
            !             b, ib, jb,             descb, ictxt)
                        coul_submtx%data_c, 1, 1, coul_submtx%blacsdata%blacs_desc, coulomb%blacsdata%blacs_desc(2))
            call coul_submtx%transpose()
         class default
            call judft_error("coul_submtx should also be mpimat")
         end select
#endif
      END SELECT
      call timestop("copy in 2")
   end subroutine copy_in_2

   subroutine copy_out_2(fmpi, sym, mpdata, hybdat, ikpt, coul_submtx, coulomb)
      USE m_types
      implicit none 
      type(t_mpi), intent(in)      :: fmpi 
      integer, intent(in)          :: ikpt
      type(t_sym), intent(in)      :: sym
      type(t_mpdata), intent(in)   :: mpdata 
      type(t_hybdat), intent(in)   :: hybdat
      class(t_mat), intent(inout)  :: coulomb 
      class(t_mat), intent(inout)  :: coul_submtx

      integer :: i, j

      call timestart("copy out 2")

      SELECT TYPE(coulomb)
      CLASS is (t_mat)
         do j = 1, mpdata%n_g(ikpt)
            do i = 1, mpdata%n_g(ikpt)
               coulomb%data_c(hybdat%nbasp+i, hybdat%nbasp + j) = conjg(coul_submtx%data_c(j, i))
            enddo 
         enddo
      class is (t_mpimat)
#ifdef CPP_SCALAPACK
         select type(coul_submtx)
         class is (t_mpimat)
            call coul_submtx%transpose()
            ! copy coul_submtx to bottom right corner of coulomb
            !call pzgemr2d(m,              n,               a,                  ia, ja,             desca, 
            call pzgemr2d(mpdata%n_g(ikpt),mpdata%n_g(ikpt),coul_submtx%data_c, 1, 1, coul_submtx%blacsdata%blacs_desc,&
            !             b,             ib,            jb,             descb, ictxt)
                        coulomb%data_c, hybdat%nbasp+1, hybdat%nbasp+1, coulomb%blacsdata%blacs_desc, coulomb%blacsdata%blacs_desc(2))
         class default
            call judft_error("coul_submtx should also be mpimat")
         end select
#endif
      end select
      call timestop("copy out 2")
   end subroutine copy_out_2
end module m_apply_inverse_olap
