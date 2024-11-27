!--------------------------------------------------------------------------------
! Copyright (c) 2016 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------

module m_eigen_diag
   use m_juDFT
   use m_available_solvers
   use m_types_mpimat
   use m_types_mat
   implicit none
   private
   public :: eigen_diag

contains

   subroutine eigen_diag(hmat, smat, ne, eig, ev)
#ifdef CPP_MPI
      use mpi
#endif
      implicit none
      class(t_mat), intent(INOUT) :: smat, hmat
      class(t_mat), allocatable, intent(OUT)   :: ev         ! eigenvectors
      integer, intent(INOUT) :: ne         ! number of eigenpairs searched (and found) on this node
      !   on input, overall number of eigenpairs searched,
      !   on output, local number of eigenpairs found
      real, intent(OUT)   :: eig(:)     ! eigenvalues

      !Locals
      logical                    :: parallel
      class(t_solver), pointer    :: solver

      select type (smat)
      class IS (t_mpimat)
#ifdef CPP_MPI
         parallel = smat%blacsdata%mpi_com /= MPI_COMM_SELF
#endif
      class default
         parallel = .false.
      end select

      solver => select_solver(parallel)

      if (.not. associated(solver%transformer)) then
         ! We solve directly the generalized eigenvalue problem
         if (solver%generalized) then
            call timestart("Diagonalization")
            call solver%solve_gev(hmat, smat, ne, eig, ev)
            call timestop("Diagonalization")
         else
            call judft_bug("Generalized solver not available?")
         end if
      else
         ! We do a reduction, to a standard problem, then solve the standard problem and transform back
         call timestart("Reduction to S-EVP")
         call solver%transformer%to_std(smat, hmat)
         call timestop("Reduction to S-EVP")
         call timestart("Diagonalization")
         call solver%solve_std(hmat, ne, eig, ev)
         call timestop("Diagonalization")
         call timestart("Backtransform of eigenvectors")
         call solver%transformer%backtrans(smat, ev)
         call timestop("Backtransform of eigenvectors")
      end if

   end subroutine

end module m_eigen_diag
