!--------------------------------------------------------------------------------
! Copyright (c) 2025 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions 
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------
module m_writeout
   use m_types_solver
   private
   type, extends(t_solver):: t_solver_debugout
   contains
      procedure        :: solve_gev => diag_writeout  !solver for generalized eigenvalue problem
   end type
   public :: t_solver_debugout, get_solver_debugout
contains
   function get_solver_debugout() result(solver)
      class(t_solver), allocatable :: solver
      allocate(t_solver_debugout :: solver)
      solver%name = "debugout"
      solver%available = .true.
      solver%parallel = .true.
      solver%serial = .true.
      solver%generalized = .true.
      solver%standard = .false.
      solver%single_precision = .false.
      solver%transform = .false.
      solver%GPU = .false.
      solver%use_sp = .false.
   end function

   subroutine diag_writeout(self, hmat, smat, ne, eig, zmat, ikpt)
      !Dummy diver: does not solve actual eigenvalue problem but simply returns a set of orthogonal vectors.
      !Could be useful for performance testing workloads in which we do not want to look at the diagonalization.
      ! A Cholesky decomp is still done to be able to do a back transform so that the resulting vector are orthonormal
      ! with respect to overlapp matrix.

      use m_types_mat
      use m_judft
      use m_io_matrix
      use m_types_mpimat
#ifdef CPP_MPI
      use mpi
#endif

      implicit none
      class(t_solver_debugout)      :: self
      class(t_mat), intent(INOUT) :: hmat, smat
      integer, intent(INOUT) :: ne
      class(t_mat), allocatable, intent(OUT)   :: zmat
      real, intent(OUT)   :: eig(:)
      integer, intent(IN) :: ikpt

      !small subroutine that does only wite the matrix to a file
      integer:: i, ii, irank, ierr, matsize
      character(len=20)::filename
#ifdef CPP_MPI
      call MPI_COMM_RANK(MPI_COMM_WORLD, irank, ierr)
#else
      irank = 0
#endif

      select type (hmat)
      type is (t_mpimat)
         matsize = hmat%global_size1
      class default
         matsize = hmat%matsize1
      end select
      !First write binary file
#ifdef CPP_HDF
      i = open_matrix(hmat%l_real, matsize, 2, 2, "hs_mat")
#else
      i = open_matrix(hmat%l_real, hmat%matsize1, 1, 2, "hs_mat")
#endif
      call write_matrix(hmat, 1, i)
      call write_matrix(smat, 2, i)
      call close_matrix(i)

      !Now the formatted matrix
      write (filename, "(a,i0)") "hmat", irank
      open (999, file=trim(filename))
      write (filename, "(a,i0)") "smat", irank
      open (998, file=trim(filename))
      do i = 1, hmat%matsize2
         do ii = 1, hmat%matsize1
            if (hmat%l_real) then
               write (999, "(2i6,f15.6)") ii, i, hmat%data_r(ii, i)
               write (998, "(2i6,f15.6)") ii, i, smat%data_r(ii, i)
            else
               write (999, "(2i6,2f15.6)") ii, i, hmat%data_c(ii, i)
               write (998, "(2i6,2f15.6)") ii, i, smat%data_c(ii, i)
            end if
         end do
      end do
      close (999)
      close (998)
#ifdef CPP_MPI
      call MPI_BARRIER(MPI_COMM_WORLD, ierr)
#endif
      call judft_error("STOP in eigen_diag:debug_diag")
   end subroutine diag_writeout
end module m_writeout
