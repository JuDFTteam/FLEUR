!--------------------------------------------------------------------------------
! Copyright (c) 2024 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------

module m_types_solver
   use m_types_mat
   use m_judfT
   implicit none
   private
   ! This is the data-type all eigenvalue solvers supported by FLEUR should implement
   type, abstract:: t_solver
      character(len=10):: name = "none"
      logical :: available
      logical :: parallel
      logical :: serial
      logical :: generalized
      logical :: standard
      logical :: single_precision
      logical :: transform
      logical :: GPU
      logical :: use_sp = .false.
   contains
      procedure        :: solve_std  !solver for standard eigenvalue problem
      procedure        :: solve_std_dp  !solver for standard eigenvalue problem (double precision)
      procedure        :: solve_std_sp  !solver for standard eigenvalue problem (single precision)
      procedure        :: solve_gev  !solver for generalized eigenvalue problem
      procedure        :: to_std     !transform the H of the generalized problem to a std problem
      procedure        :: backtrans  !transform the Eigenvalue back to the generalized problem
   end type

   !Simple solver that only stops
   type, extends(t_solver):: t_solver_stop
   contains
      procedure        :: solve_gev => solve_stop  !solver for generalized eigenvalue problem
   end type

   public:: t_solver, t_solver_stop, get_solver_stop

contains
   function get_solver_stop() result(solver)
      type(t_solver_stop), pointer::solver
      allocate (solver)
      solver%name = "stop"
      solver%available = .true.
      solver%parallel = .true.
      solver%serial = .true.
      solver%generalized = .true.
      solver%standard = .false.
      solver%single_precision = .false.
      solver%transform = .false.
      solver%GPU = .true.
      solver%use_sp = .false.

   end function

   subroutine solve_stop(self, hmat, smat, ne, eig, zmat, ikpt)
      implicit none
      class(t_solver_stop)       :: self
      class(t_mat), intent(INOUT)  :: hmat, smat
      integer, intent(INOUT)       :: ne
      class(t_mat), allocatable, intent(OUT)    :: zmat
      real, intent(OUT)           :: eig(:)
      integer, intent(IN)         :: ikpt

      call judft_error("FLEUR stopped as -diag stop was choosen")
   end subroutine

   subroutine solve_std(self, hmat, ne, eig, zmat)
      implicit none
      class(t_solver)            :: self
      class(t_mat), intent(INOUT)  :: hmat
      integer, intent(INOUT)       :: ne
      class(t_mat), allocatable, intent(OUT)    :: zmat
      real, intent(OUT)           :: eig(:)

      if (self%use_sp) then
         call self%solve_std_sp(hmat, ne, eig, zmat)
      else
         call self%solve_std_dp(hmat, ne, eig, zmat)
      end if
   end subroutine

   subroutine solve_gev(self, hmat, smat, ne, eig, zmat, ikpt)
      implicit none
      class(t_solver)                    :: self
      class(t_mat), intent(INOUT)          :: hmat, smat
      integer, intent(INOUT)               :: ne
      class(t_mat), allocatable, intent(OUT):: zmat
      real, intent(OUT)                    :: eig(:)
      integer, intent(IN)                  :: ikpt

      call judft_bug("Not implemented", calledby="solve_std")
   end subroutine
   subroutine solve_std_sp(self, hmat, ne, eig, zmat)
      implicit none
      class(t_solver)                     :: self
      class(t_mat), intent(INOUT)          :: hmat
      integer, intent(INOUT)               :: ne
      class(t_mat), allocatable, intent(OUT):: zmat
      real, intent(OUT)                    :: eig(:)

      call judft_bug("Not implemented", calledby="solve_std_sp")
   end subroutine
   subroutine solve_std_dp(self, hmat, ne, eig, zmat)
      implicit none
      class(t_solver)                     :: self
      class(t_mat), intent(INOUT)          :: hmat
      integer, intent(INOUT)               :: ne
      class(t_mat), allocatable, intent(OUT):: zmat
      real, intent(OUT)                    :: eig(:)

      call judft_bug("Not implemented", calledby="solve_std_dp")
   end subroutine
   subroutine to_std(self, hmat, smat)
      implicit none
      class(t_solver)                     :: self
      class(t_mat), intent(INOUT)          :: hmat, smat

      call judft_bug("Not implemented", calledby="to_std")
   end subroutine
   subroutine backtrans(self, smat, zmat)
      implicit none
      class(t_solver)                     :: self
      class(t_mat), intent(INOUT)         :: smat
      class(t_mat), intent(INOUT)         :: zmat

      call judft_bug("Not implemented", calledby="backtrans")
   end subroutine

end module m_types_solver
