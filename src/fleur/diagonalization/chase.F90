!--------------------------------------------------------------------------------
! Copyright (c) 2025 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------
! Added MPI implementation, DW 2018
!--------------------------------------------------------------------------------
module m_chase
   use m_judft
   use m_constants
   use m_types_mat
   use m_types_mpimat
   use m_types_solver
   use , INTRINSIC :: iso_c_binding, ONLY: c_char
#ifdef CPP_MPI
   use mpi
#endif
#ifdef CPP_CHASE
   use chase_diag  !this is the interface to chase provided in the chase library
#endif
   implicit none

   private
   type, extends(t_solver)::t_solver_chase
   contains
      procedure        :: solve_std_sp => chase_diag_sp
      procedure        :: solve_std_dp => chase_diag_dp
   end type
   public :: get_solver_chase

contains

   function get_solver_chase() result(solver)
      class(t_solver), allocatable :: solver
      allocate (t_solver_chase :: solver)
      solver%name = "chase"
#ifdef CPP_CHASE
      solver%available = .true.
#else
      solver%available = .false.
#endif
      solver%parallel = .true.
      solver%serial = .true.
      solver%generalized = .false.
      solver%standard = .true.
      solver%single_precision = .true.
      solver%transform = .false.
      solver%GPU = .true.
      solver%use_sp = .false.
   end function

   subroutine chase_diag_dp(self, hmat, ne, eig, zmat)
      implicit none
      class(t_solver_chase)       :: self
      class(t_mat), intent(INOUT)  :: hmat
      integer, intent(INOUT)      :: ne
      class(t_mat), allocatable, intent(OUT)    :: zmat
      real, intent(OUT)           :: eig(:)

      select type (hmat)
      type is (t_mat)
         call chase_serial_dp(hmat, ne, eig, zmat)
      type is (t_mpimat)
         call chase_mpi_dp(hmat, ne, eig, zmat)
      end select
   end subroutine

   subroutine chase_diag_sp(self, hmat, ne, eig, zmat)
      implicit none
      class(t_solver_chase)       :: self
      class(t_mat), intent(INOUT)  :: hmat
      integer, intent(INOUT)      :: ne
      class(t_mat), allocatable, intent(OUT)    :: zmat
      real, intent(OUT)           :: eig(:)

      select type (hmat)
      type is (t_mat)
         call chase_serial_sp(hmat, ne, eig, zmat)
      type is (t_mpimat)
         call judft_bug("Chase in SP for parallel matrices not implemented")
         !call chase_mpi_diag_sp(hmat,ne,eig,zmat)
      end select
   end subroutine

   subroutine chase_serial_dp(hmat, ne, eig, zmat)
      !Simple driver to solve Standard Eigenvalue Problem using ChASE routine
      implicit none
      type(t_mat), intent(INOUT), VOLATILE  :: hmat
      integer, intent(INOUT)      :: ne
      class(t_mat), allocatable, intent(OUT)    :: zmat
      real, intent(OUT)           :: eig(:)

      !These chase parameters might to be adjusted
      real, parameter ::   tol = 1e-10
      character(kind=c_char), parameter ::  mode = 'R', opt = 'S', qr = 'C'
      integer, parameter  :: deg = 20
#ifdef CPP_CHASE
      integer :: nex !extra search space
      integer :: init  !status variable
      !chase will modify these variables in call to xchase even though these are not arguments!!
      real, allocatable, VOLATILE :: zr(:, :), eigval(:)
      complex, allocatable, VOLATILE :: zc(:, :)
      call timestart("CHASE STD")
      nex = 0.2*ne
      allocate (eigval(ne + nex))
      allocate (t_mat::zmat)
      call zmat%init(hmat%l_real, hmat%matsize1, ne)
      call timestart("CHASE STD U2L")
      call hmat%u2l() !chase needs full matrix not only upper part!
      call timestop("CHASE STD U2L")
      if (hmat%l_real) then
         allocate (zr(hmat%matsize1, ne + nex))
         ! Initialize of ChASE
         call timestart("CHASE STD INIT")
         call dchase_init(size(hmat%data_r, 1), ne, nex, hmat%data_r, size(hmat%data_r, 1), zr, eigval, init)
         call timestop("CHASE STD INIT")
         !Solve eigenvalue problem
         call timestart("CHASE STD SOLVE")
         call dchase(deg, tol, mode, opt, qr)
         call timestop("CHASE STD SOLVE")
         ! finalize and clean up
         call timestart("CHASE STD FINALIZE")
         call dchase_finalize(init)
         call timestop("CHASE STD FINALIZE")

         zmat%data_r = zr(:, :ne)
      else
         allocate (zc(hmat%matsize1, ne + nex))
         ! Initialize of ChASE
         call timestart("CHASE STD INIT")
         call zchase_init(hmat%matsize1, ne, nex, hmat%data_c, size(hmat%data_c, 1), zc, eigval, init)
         call timestop("CHASE STD INIT")
         !Solve eigenvalue problem
         call timestart("CHASE STD SOLVE")
         call zchase(deg, tol, mode, opt, qr)
         call timestop("CHASE STD SOLVE")
         ! finalize and clean up
         call timestart("CHASE STD FINALIZE")
         call zchase_finalize(init)
         call timestop("CHASE STD FINALIZE")
         zmat%data_c = zc(:, :ne)
      end if
      eig(:ne) = eigval(:ne)
#endif
      call timestop("CHASE STD")
   end subroutine chase_serial_dp

   subroutine chase_serial_sp(hmat, ne, eig, zmat)
      !Simple driver to solve Standard Eigenvalue Problem using ChASE routine in single precision
      implicit none
      type(t_mat), intent(INOUT)  :: hmat
      integer, intent(INOUT)      :: ne
      class(t_mat), allocatable, intent(OUT)    :: zmat
      real, intent(OUT)           :: eig(:)

      integer, parameter:: sp = selected_real_kind(6)
      !These chase parameters might to be adjusted
      real(sp), parameter ::   tol = 1e-6
      character, parameter ::  mode = 'R', opt = 'S', qr = 'C'
      integer, parameter  :: deg = 20

#ifdef CPP_CHASE
      integer :: nex !extra search space
      integer :: init  !status variable
      !chase will modify these variables in call to xchase even though these are not arguments!!
      real(sp), allocatable, volatile :: zr(:, :), eigval(:)
      complex(sp), allocatable, volatile :: zc(:, :)
      real(sp), allocatable :: hr(:, :)
      complex(sp), allocatable :: hc(:, :)
      call timestart("CHASE STD-SP")
      nex = 0.2*ne
      allocate (eigval(ne + nex))
      allocate (t_mat::zmat)
      call zmat%init(hmat%l_real, hmat%matsize1, ne)

      call timestart("CHASE STD-SP U2L")
      call hmat%u2l() !chase needs full matrix not only upper part!
      call timestop("CHASE STD-SP U2L")
      if (hmat%l_real) then
         allocate (zr(hmat%matsize1, ne + nex))
         allocate (hr(hmat%matsize1, hmat%matsize2))
         hr = hmat%data_r  !cast to sp
         ! Initialize of ChASE
         call timestart("CHASE STD-SP INIT")
         call schase_init(hmat%matsize1, ne, nex, hr, size(hr, 1), zr, eigval, init)
         call timestop("CHASE STD-SP INIT")
         !Solve eigenvalue problem
         call timestart("CHASE STD-SP SOLVE")
         call schase(deg, tol, mode, opt, qr)
         call timestop("CHASE STD-SP SOLVE")
         ! finalize and clean up
         call timestart("CHASE STD-SP FINALIZE")
         call schase_finalize(init)
         call timestop("CHASE STD-SP FINALIZE")
         zmat%data_r = zr(:, :ne)
      else
         allocate (zc(hmat%matsize1, ne + nex))
         allocate (hc(hmat%matsize1, hmat%matsize2))
         hc = hmat%data_c
         ! Initialize of ChASE
         call timestart("CHASE STD-SP INIT")
         call cchase_init(hmat%matsize1, ne, nex, hc, size(hc, 1), zc, eigval, init)
         call timestop("CHASE STD-SP INIT")
         !Solve eigenvalue problem
         call timestart("CHASE STD-SP SOLVE")
         call cchase(deg, tol, mode, opt, qr)
         call timestop("CHASE STD-SP SOLVE")
         ! finalize and clean up
         call timestart("CHASE STD-SP FINALIZE")
         call cchase_finalize(init)
         call timestop("CHASE STD-SP FINALIZE")
         zmat%data_c = zc(:, :ne)
      end if
      eig(:ne) = eigval(:ne)
#endif
      call timestop("CHASE STD-SP")
   end subroutine chase_serial_sp

   subroutine chase_mpi_dp(hmat, ne, eig, zmat)
      !Simple driver to solve Standard Eigenvalue Problem using ChASE routine
      implicit none
      type(t_mpimat), intent(INOUT)  :: hmat
      integer, intent(INOUT)         :: ne
      class(t_mat), allocatable, intent(OUT), volatile    :: zmat
      real, intent(OUT)           :: eig(:)

      !These chase parameters might to be adjusted
      real, parameter ::   tol = 1e-10
      character(kind=c_char), parameter ::  mode = 'R', opt = 'S', qr = 'C'
      character(kind=c_char), parameter ::  grid_major = 'C' !major of 2D MPI grid. Row major: grid_major=’R’, column major: grid_major=’C’
      integer, parameter  :: deg = 20
#ifdef CPP_CHASE
      integer:: mbsize, nbsize, irsrc, icsrc, dim0, dim1, myprow, mypcol
      integer :: comm_1d, comm_2d, ierr, myrank, ne_global
      integer :: nex !extra search space
      integer :: init  !status variable
      integer :: j, gcol
      integer, external :: indxl2g
      !chase will modify these variables in call to xchase even though these are not arguments!!
      real, allocatable, volatile :: eigval(:)
      real, allocatable :: zmat_chase_r(:, :)
      complex, allocatable :: zmat_chase_c(:, :)

      call timestart("CHASE MPI-STD")
      ne_global = ne
      nex = 0.2*ne
      allocate (eigval(ne + nex))

      !setup ChASE
      mbsize = hmat%blacsdata%blacs_desc(5) !block size for the block-cyclic distribution for the rows of global matrix
      nbsize = hmat%blacsdata%blacs_desc(6) !block size for the block-cyclic distribution for the cloumns of global matrix
      irsrc = hmat%blacsdata%blacs_desc(7) !process row over which the first row of the global matrix h is distributed
      icsrc = hmat%blacsdata%blacs_desc(8) !process column over which the first column of the global matrix h is distributed.

      !determine the processor grid
      !dim0 :row number of 2D MPI grid
      !dim1 :column number of 2D MPI grid
      call BLACS_GRIDINFO(hmat%blacsdata%blacs_desc(2), dim0, dim1, MYPROW, MYPCOL)
      call timestart("CHASE MPI COMMS")
      call create_mpi_comms(hmat%blacsdata%mpi_com, hmat%blacsdata%blacs_desc(2), comm_2d, comm_1d)
      call timestop("CHASE MPI COMMS")

      if (comm_1d == MPI_COMM_NULL .or. comm_2d == MPI_COMM_NULL) then
         call juDFT_error("ChASE communicator creation failed (MPI_COMM_NULL)", calledby="chase_mpi_dp")
      end if

      call timestart("CHASE MPI U2L")
      call hmat%u2l() !chase needs full matrix not only upper part!
      call timestop("CHASE MPI U2L")

      if (hmat%l_real) then
         allocate (zmat_chase_r(hmat%matsize1, ne + nex))
         ! Initialize of ChASE
         call timestart("CHASE MPI INIT")
         call pdchase_init_blockcyclic(hmat%global_size1, ne, nex, mbsize, nbsize, hmat%data_r, hmat%matsize1, &
                                       zmat_chase_r, eigval, dim0, dim1, grid_major, irsrc, icsrc, comm_2d, init)
         call timestop("CHASE MPI INIT")
!Solve eigenvalue problem
         call timestart("CHASE MPI SOLVE")
         call pdchase(deg, tol, mode, opt, qr)
         call timestop("CHASE MPI SOLVE")
         ! finalize and clean up
         call timestart("CHASE MPI FINALIZE")
         call pdchase_finalize(init)
         call timestop("CHASE MPI FINALIZE")
      else
         allocate (zmat_chase_c(hmat%matsize1, ne + nex))
         ! Initialize of ChASE
         call timestart("CHASE MPI INIT")
         call pzchase_init_blockcyclic(hmat%global_size1, ne, nex, mbsize, nbsize, hmat%data_c, hmat%matsize1, &
                                       zmat_chase_c, eigval, dim0, dim1, grid_major, irsrc, icsrc, comm_2d, init)
         call timestop("CHASE MPI INIT")
         !Solve eigenvalue problem
         call timestart("CHASE MPI SOLVE")
         call pzchase(deg, tol, mode, opt, qr)
         call timestop("CHASE MPI SOLVE")
         ! finalize and clean up
         call timestart("CHASE MPI FINALIZE")
         call pzchase_finalize(init)
         call timestop("CHASE MPI FINALIZE")
      end if

      ! Redistribute ChASE vectors into 2D block-cyclic layout expected by ELPA.
      allocate (t_mpimat::zmat)
      select type (zmat)
      type is (t_mpimat)
         call zmat%init(hmat)
         if (hmat%l_real) then
            do j = 1, ne
           gcol = indxl2g(j, zmat%blacsdata%blacs_desc(6), zmat%blacsdata%mycol, zmat%blacsdata%blacs_desc(8), zmat%blacsdata%npcol)
               zmat%data_r(:, j) = zmat_chase_r(:, gcol)
            end do
         else
            do j = 1, ne
           gcol = indxl2g(j, zmat%blacsdata%blacs_desc(6), zmat%blacsdata%mycol, zmat%blacsdata%blacs_desc(8), zmat%blacsdata%npcol)
               zmat%data_c(:, j) = zmat_chase_c(:, gcol)
            end do
         end if
      class default
         call juDFT_error("Wrong type for ChASE eigenvector output", calledby="chase_mpi_dp")
      end select

      call MPI_COMM_RANK(hmat%blacsdata%mpi_com, myrank, ierr)
      if (myrank == 0) eig(:ne_global) = eigval(:ne_global)
      call MPI_BCAST(ne_global, 1, MPI_INTEGER, 0, hmat%blacsdata%mpi_com, ierr)
      call MPI_BCAST(eig(:ne_global), ne_global, MPI_DOUBLE_PRECISION, 0, hmat%blacsdata%mpi_com, ierr)
      ne = ne_global

      if (allocated(zmat_chase_r)) deallocate (zmat_chase_r)
      if (allocated(zmat_chase_c)) deallocate (zmat_chase_c)

      call timestart("CHASE MPI COMM_FREE")
      call MPI_COMM_FREE(comm_1d, ierr)
      call MPI_COMM_FREE(comm_2d, ierr)
      call timestop("CHASE MPI COMM_FREE")
      call timestop("CHASE MPI-STD")
#endif
   end subroutine chase_mpi_dp

   subroutine create_mpi_comms(parent_comm, icontext, comm_2d, comm_1d)
      integer, intent(in) :: parent_comm, icontext
      integer, intent(out):: comm_2d, comm_1d
#ifdef CPP_CHASE

      integer:: isize, ierr, i, col, row, no_col, no_row, myprow, mypcol
      integer, allocatable:: map2d(:)
      integer :: group, group_2d, group_1d

      call timestart("CHASE MPI COMMS BUILD")
      call BLACS_GRIDINFO(icontext, no_row, no_col, MYPROW, MYPCOL)
      call MPI_COMM_SIZE(parent_comm, isize, ierr)
      allocate (map2d(0:isize - 1))

      do i = 0, isize - 1
         call blacs_pcoord(icontext, i, ROW, COL)
         map2d(row + no_row*col) = i
      end do

      !create 2d communicator
      call MPI_COMM_group(parent_comm, group, ierr)
      call MPI_Group_incl(group, isize, map2d, group_2d, ierr)
      call MPI_COMM_create_group(parent_comm, group_2d, 1, comm_2d, ierr)

      !create 1d communicator per process column (all ranks get a valid comm_1d)
      call MPI_COMM_split(parent_comm, MYPCOL, MYPROW, comm_1d, ierr)

      ! group_1d is unused with MPI_COMM_split, keep a valid null handle for cleanup
      group_1d = MPI_GROUP_NULL

      call MPI_group_free(group_2d, ierr)
      if (group_1d /= MPI_GROUP_NULL) call MPI_group_free(group_1d, ierr)
      call MPI_group_free(group, ierr)
      call timestop("CHASE MPI COMMS BUILD")
#endif

   end subroutine

end module m_chase
