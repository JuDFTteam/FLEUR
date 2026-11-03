!--------------------------------------------------------------------------------
! Copyright (c) 2026 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions 
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------
module m_apply_inverse_olap
   use m_glob_tofrom_loc
   USE m_types_mpimat
   USE m_olap, ONLY: olap_pw
   USE m_judft
   USE m_types_mat
   USE m_types_atoms
   USE m_types_cell
   USE m_types_hybdat
   USE m_types_mpdata
   USE m_types_mpi
   USE m_types_sym
   IMPLICIT NONE
   PRIVATE
   PUBLIC :: apply_inverse_olaps, copy_in_2, copy_out_2
contains
   subroutine apply_inverse_olaps(mpdata, atoms, cell, hybdat, fmpi, sym, ikpt, coulomb)
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
      integer         :: ir_first, ir_last
      complex         :: cdum

      call timestart("solve olap linear eq. sys")
      nbasm = hybdat%nbasm(ikpt)

      ir_first = hybdat%n_mt + 1
      ir_last = hybdat%n_mt + mpdata%n_g(ikpt)

      CALL olap%alloc(.false., mpdata%n_g(ikpt), mpdata%n_g(ikpt), 0.0)
      !calculate IR overlap-matrix
      CALL olap_pw(olap, mpdata%g(:, mpdata%gptm_ptr(:mpdata%n_g(ikpt), ikpt)), mpdata%n_g(ikpt), atoms, cell, fmpi)

      ! perform O^-1 * coulomb%data_c(ir_first:ir_last, :) = x
      ! rewritten as O * x = C

      loc_size = 0
      do i = 1, nbasm
         call glob_to_loc(fmpi, i, pe_i, i_loc)
         if (fmpi%n_rank == pe_i) loc_size = loc_size + 1
      end do

      call timestart("copy in 1")
      allocate(t_mat::coul_submtx)
      call coul_submtx%alloc(.false., mpdata%n_g(ikpt), loc_size)
      coul_submtx%data_c(:, :) = coulomb%data_c(ir_first:ir_last, :)
      call timestop("copy in 1")

      !$acc data copyin(olap, olap%data_r, olap%data_c, coul_submtx) copy(coul_submtx%data_r, coul_submtx%data_c)
         call olap%linear_problem(coul_submtx)
      !$acc end data
      call timestart("copy out 1")
      coulomb%data_c(ir_first:ir_last, :) = coul_submtx%data_c
      call coul_submtx%free()
      deallocate(coul_submtx)
      call timestop("copy out 1")


      ! perform coulomb%data_c(ir_first:, ir_first:ir_last) * O^-1 = X, rows IR and (films) VAC:
      ! the VAC x IR block is copied to the upper triangle by l2u later
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
      call coul_submtx%free()
      deallocate(coul_submtx)
      call timestop("solve olap linear eq. sys")
   end subroutine apply_inverse_olaps

   subroutine copy_in_2(fmpi, sym, mpdata, hybdat, coulomb, ikpt, coul_submtx)
      implicit none 
      type(t_mpi), intent(in)      :: fmpi 
      integer, intent(in)          :: ikpt
      type(t_sym), intent(in)      :: sym
      type(t_mpdata), intent(in)   :: mpdata 
      type(t_hybdat), intent(in)   :: hybdat
      class(t_mat), intent(in)     :: coulomb 
      class(t_mat), intent(inout), allocatable  :: coul_submtx

      integer :: i, j, ierr, i_loc, j_loc, pe_i, pe_j, n_g, n_row
      complex :: cdum
      type(t_mpimat) :: rect

      call timestart("copy in 2")

      n_g = mpdata%n_g(ikpt)
      n_row = hybdat%nbasm(ikpt) - hybdat%n_mt

      SELECT TYPE(coulomb)
      CLASS is (t_mat)
         allocate(t_mat::coul_submtx)
         call coul_submtx%alloc(.false., n_g, n_row)
         do j = 1, n_g
            do i = 1, n_row
               coul_submtx%data_c(j, i) = conjg(coulomb%data_c(hybdat%n_mt+i, hybdat%n_mt + j))
            enddo 
         enddo
      class is (t_mpimat)
#ifdef CPP_SCALAPACK
         allocate(t_mpimat::coul_submtx)
         select type(coul_submtx)
         class is (t_mpimat)
            if (n_row == n_g) then
               call coul_submtx%init(.False., n_g, n_g, fmpi%sub_comm, MPIMAT_2D_BLOCK_CYCLIC)
               call pzgemr2d(n_g, n_g, coulomb%data_c, hybdat%n_mt + 1, hybdat%n_mt + 1, coulomb%blacsdata%blacs_desc,&
                        coul_submtx%data_c, 1, 1, coul_submtx%blacsdata%blacs_desc, coulomb%blacsdata%blacs_desc(2))
               call coul_submtx%transpose()
            else
               ! films: rectangular, t_mpimat%transpose handles square matrices only
               call rect%init(.False., n_row, n_g, fmpi%sub_comm, MPIMAT_2D_BLOCK_CYCLIC)
               call pzgemr2d(n_row, n_g, coulomb%data_c, hybdat%n_mt + 1, hybdat%n_mt + 1, coulomb%blacsdata%blacs_desc,&
                           rect%data_c, 1, 1, rect%blacsdata%blacs_desc, coulomb%blacsdata%blacs_desc(2))
               call coul_submtx%init(rect, n_g, n_row)
               ! init_template keeps the leading dimension of the template
               coul_submtx%blacsdata%blacs_desc(9) = max(1, coul_submtx%matsize1)
               call pztranc(n_g, n_row, cmplx(1.0, 0.0), rect%data_c, 1, 1, rect%blacsdata%blacs_desc, &
                            cmplx(0.0, 0.0), coul_submtx%data_c, 1, 1, coul_submtx%blacsdata%blacs_desc)
               call rect%free()
            endif
         class default
            call judft_error("coul_submtx should also be mpimat")
         end select
#endif
      END SELECT
      call timestop("copy in 2")
   end subroutine copy_in_2

   subroutine copy_out_2(fmpi, sym, mpdata, hybdat, ikpt, coul_submtx, coulomb)
      implicit none 
      type(t_mpi), intent(in)      :: fmpi 
      integer, intent(in)          :: ikpt
      type(t_sym), intent(in)      :: sym
      type(t_mpdata), intent(in)   :: mpdata 
      type(t_hybdat), intent(in)   :: hybdat
      class(t_mat), intent(inout)  :: coulomb 
      class(t_mat), intent(inout)  :: coul_submtx

      integer :: i, j, n_g, n_row
      type(t_mpimat) :: rect

      call timestart("copy out 2")

      n_g = mpdata%n_g(ikpt)
      n_row = hybdat%nbasm(ikpt) - hybdat%n_mt

      SELECT TYPE(coulomb)
      CLASS is (t_mat)
         do j = 1, n_g
            do i = 1, n_row
               coulomb%data_c(hybdat%n_mt+i, hybdat%n_mt + j) = conjg(coul_submtx%data_c(j, i))
            enddo 
         enddo
      class is (t_mpimat)
#ifdef CPP_SCALAPACK
         select type(coul_submtx)
         class is (t_mpimat)
            if (n_row == n_g) then
               call coul_submtx%transpose()
               call pzgemr2d(n_g, n_g, coul_submtx%data_c, 1, 1, coul_submtx%blacsdata%blacs_desc,&
                        coulomb%data_c, hybdat%n_mt + 1, hybdat%n_mt + 1, coulomb%blacsdata%blacs_desc, coulomb%blacsdata%blacs_desc(2))
            else
               call rect%init(coul_submtx, n_row, n_g)
               ! init_template keeps the leading dimension of the template
               rect%blacsdata%blacs_desc(9) = max(1, rect%matsize1)
               call pztranc(n_row, n_g, cmplx(1.0, 0.0), coul_submtx%data_c, 1, 1, coul_submtx%blacsdata%blacs_desc, &
                            cmplx(0.0, 0.0), rect%data_c, 1, 1, rect%blacsdata%blacs_desc)
               call pzgemr2d(n_row, n_g, rect%data_c, 1, 1, rect%blacsdata%blacs_desc, &
                        coulomb%data_c, hybdat%n_mt + 1, hybdat%n_mt + 1, coulomb%blacsdata%blacs_desc, coulomb%blacsdata%blacs_desc(2))
               call rect%free()
            endif
         class default
            call judft_error("coul_submtx should also be mpimat")
         end select
#endif
      end select
      call timestop("copy out 2")
   end subroutine copy_out_2
end module m_apply_inverse_olap
