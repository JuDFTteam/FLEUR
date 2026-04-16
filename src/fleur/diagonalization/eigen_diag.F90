!--------------------------------------------------------------------------------
! Copyright (c) 2026 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------

module m_eigen_diag
   !! Module provides the high level entry point for solving a generalized eigenvalue problem
   !! The solver actually used is determined by a call to [[select_solver]] from [[m_available_solvers]].
   use m_juDFT
   use m_available_solvers
   use m_types_mpimat
   use m_types_mat
   use m_types_solver
   use m_lapack
   implicit none
   private
   public :: eigen_diag,eigen_diag_std

contains

   subroutine eigen_diag(hmat, smat, ne, eig, ev, ikpt)
      !! Solve generalized eigenvalue problem
#ifdef CPP_MPI
      use mpi
#endif
      implicit none
      class(t_mat), intent(INOUT) :: smat, hmat !! overlapp matrix and Hamiltonian
      class(t_mat), allocatable, intent(OUT)   :: ev         !! eigenvectors
      integer, intent(INOUT) :: ne         !! number of eigenpairs searched (and found) on this node
      !!   on input, overall number of eigenpairs searched,
      !!   on output, local number of eigenpairs found
      real, intent(OUT)   :: eig(:)     !! eigenvalues (must be allocated to size ne before)
      integer, intent(IN) :: ikpt     !! The index of the k-point, just to have more information for debugging

      !Locals
      logical                       :: parallel
      class(t_solver),allocatable   :: solver,transform

      select type (hmat)
      class IS (t_mpimat)
#ifdef CPP_MPI
         parallel = hmat%blacsdata%mpi_com /= MPI_COMM_SELF
#endif
      class default
         parallel = .false.
      end select

      call select_solver(parallel,diag_solver=solver,diag_transform=transform)
      
      if (.not. allocated(transform)) then
         ! We solve directly the generalized eigenvalue problem
         if (solver%generalized) then
            call timestart("Diagonalization")
            call solver%solve_gev(hmat, smat, ne, eig, ev, ikpt)
            call timestop("Diagonalization")
         else
            call judft_bug("Generalized solver not available?")
         end if
      else
         ! We do a reduction, to a standard problem, then solve the standard problem and transform back
         call timestart("Reduction to S-EVP")
         call transform%to_std(hmat, smat)
         call timestop("Reduction to S-EVP")
         call timestart("Diagonalization")
         print *,"Solver:",solver%name
         call solver%solve_std(hmat, ne, eig, ev)
         call timestop("Diagonalization")
         call timestart("Backtransform of eigenvectors")
         call transform%backtrans(smat, ev)
         call timestop("Backtransform of eigenvectors")
      end if
   end subroutine

   subroutine eigen_diag_std(hmat, ne, eig, ev)
      !! Solve standard eigenvalue problem
#ifdef CPP_MPI
      use mpi
#endif
      implicit none
      class(t_mat), intent(INOUT) ::  hmat !!  Hamiltonian
      class(t_mat), allocatable, intent(OUT)   :: ev         !! eigenvectors
      integer, intent(INOUT) :: ne         !! number of eigenpairs searched (and found) on this node
      !!   on input, overall number of eigenpairs searched,
      !!   on output, local number of eigenpairs found
      real, intent(OUT)   :: eig(:)     !! eigenvalues (must be allocated to size ne before)
      
      !Locals
      logical                       :: parallel
      class(t_solver),allocatable   :: solver,transform

      select type (hmat)
      class IS (t_mpimat)
#ifdef CPP_MPI
         parallel = hmat%blacsdata%mpi_com /= MPI_COMM_SELF
#endif
      class default
         parallel = .false.
      end select

      call select_solver(parallel,diag_solver=solver,diag_transform=transform)
      if (solver%standard) THEN
         call timestart("Diagonalization")
         call solver%solve_std(hmat, ne, eig, ev)
         call timestop("Diagonalization")
      else  
         call judft_error("Solver does not support standard eigenvalue problems.")
      endif
      end subroutine

    
end module m_eigen_diag
