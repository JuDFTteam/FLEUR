! Copyright (c) 2018 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions
! of the MIT license as expressed in the LICENSE file in more detail.
!
! @authors: Miriam Hinzen, Gregor Michalicek
! Added MPI implementation, DW 2018
!--------------------------------------------------------------------------------
module m_dummy_diag
   use m_judft
   use m_constants
   use m_types_solver
   implicit none
   private
   type, extends(t_solver):: t_solver_dummy
   contains
      procedure        :: solve_gev => dummy_diag  !solver for generalized eigenvalue problem
   end type

   public t_solver_dummy, get_solver_dummy

contains

   function get_solver_dummy() result(solver)
      type(t_solver_dummy), pointer::solver
      allocate (solver)
      solver%name = "dummy"
      solver%available = .true.
      solver%parallel = .false.
      solver%serial = .true.
      solver%generalized = .true.
      solver%standard = .false.
      solver%single_precision = .false.
      solver%transform = .false.
      solver%GPU = .true.
      solver%use_sp = .false.
   end function

   subroutine dummy_diag(self, hmat, smat, ne, eig, zmat, ikpt)
      !Dummy diver: does not solve actual eigenvalue problem but simply returns a set of orthogonal vectors.
      !Could be useful for performance testing workloads in which we do not want to look at the diagonalization.
      ! A Cholesky decomp is still done to be able to do a back transform so that the resulting vector are orthonormal
      ! with respect to overlapp matrix.

      use m_types_mat
      use m_judft

      implicit none
      class(t_solver_dummy)                  :: self
      class(t_mat), intent(INOUT) :: hmat, smat
      integer, intent(INOUT) :: ne
      class(t_mat), allocatable, intent(OUT)   :: zmat
      real, intent(OUT)   :: eig(:)
      integer, intent(IN) :: ikpt

      integer            :: nev, lwork, liwork, n
      integer            :: info

      allocate (t_mat::zmat)
      call zmat%alloc(hmat%l_real, hmat%matsize1, ne)

      if (hmat%l_real) then
         ! --> start with Cholesky factorization of b ( so that b = l * l^t)
         ! --> b is overwritten by l
         call dpotrf('U', smat%matsize1, smat%data_r, size(smat%data_r, 1), info)
         if (info .ne. 0) then
            write (*, *) 'Error in dpotrf: info =', info
            call juDFT_error("Diagonalization failed", calledby="lapack_singlePrec_diag")
         end if

         ! --> solve a' * z' = eig * z' for eigenvalues eig between lb und ub
         zmat%data_r = 0.0
         do n = 1, ne
            eig(ne) = -0.1 + ne*1e-5
            zmat%data_r(ne, ne) = 1.0
         end do
         ! --> recover the generalized eigenvectors z by solving z' = l^t * z
         call dtrtrs('U', 'N', 'N', hmat%matsize1, nev, smat%data_r, smat%matsize1, zMat%data_r, zmat%matsize1, info)
         if (info .ne. 0) then
            write (oUnit, *) 'Error in dtrtrs: info =', info
            call juDFT_error("Diagonalization failed", calledby="lapack_singlePrec_diag")
         end if

      else

         ! --> start with Cholesky factorization of b ( so that b = l * l^t)
         ! --> b is overwritten by l
         call zpotrf('U', smat%matsize1, smat%data_c, size(smat%data_c, 1), info)
         if (info .ne. 0) then
            write (*, *) 'Error in zpotrf: info =', info
            call juDFT_error("Diagonalization failed", calledby="chase_diag")
         end if

         ! --> solve a' * z' = eig * z' for eigenvalues eig between lb und ub
         zmat%data_c = 0.0
         do n = 1, ne
            eig(ne) = -0.1 + ne*1e-5
            zmat%data_c(ne, ne) = 1.0
         end do

         ! --> recover the generalized eigenvectors z by solving z' = l^t * z
         call ztrtrs('U', 'N', 'N', hmat%matsize1, nev, smat%data_c, smat%matsize1, zMat%data_c, zmat%matsize1, info)
         if (info .ne. 0) then
            write (oUnit, *) 'Error in ztrtrs: info =', info
            call juDFT_error("Diagonalization failed", calledby="chase_diag")
         end if

      end if
   end subroutine dummy_diag

end module m_dummy_diag
