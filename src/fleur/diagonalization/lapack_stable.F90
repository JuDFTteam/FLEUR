!--------------------------------------------------------------------------------
! Copyright (c) 2026 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------
module m_lapack_stable
   ! LAPACK solver using the divide-and-conquer drivers (?sygvd/?hegvd, ?syevd/?heevd).
   ! All eigenpairs are computed and the lowest ne are returned. The reduction to the
   ! standard problem and the back-transformation are inherited from t_solver_lapack_fast.
   use m_types_solver
   use m_types_mat
   use m_lapack_fast
   use m_judft
   implicit none
   private
   type, extends(t_solver_lapack_fast)::t_solver_lapack_stable
   contains
      procedure        :: solve_std_dp => lapack_dc_diag  !solver for standard eigenvalue problem
      procedure        :: solve_std_sp => lapack_dc_diag_sp  !solver for standard eigenvalue problem
      procedure        :: solve_gev => lapack_dc_gev  !solver for generalized eigenvalue problem
   end type
   public :: t_solver_lapack_stable, get_solver_lapack_stable

contains

   function get_solver_lapack_stable() result(solver)
      class(t_solver), allocatable :: solver
      allocate (t_solver_lapack_stable :: solver)
      solver%name = "lapack_stable"
      solver%available = .true.
      solver%parallel = .false.
      solver%serial = .true.
      solver%generalized = .true.
      solver%standard = .true.
      solver%single_precision = .true.
      solver%transform = .true.
      solver%GPU = .false.
      solver%use_sp = .false.
   end function

   subroutine lapack_dc_gev(self, hmat, smat, ne, eig, zmat, ikpt)
      ! Generalized eigenvalue problem via dsygvd/zhegvd
      implicit none
      class(t_solver_lapack_stable)          :: self
      class(t_mat), intent(INOUT)            :: hmat, smat
      integer, intent(INOUT)                 :: ne
      class(t_mat), allocatable, intent(OUT) :: zmat
      real, intent(OUT)                      :: eig(:)
      integer, intent(IN)                    :: ikpt

      integer              :: n, info, lwork, lrwork, liwork(1)
      real                 :: rwork_dum(1)
      complex              :: work_dum(1)
      integer, allocatable :: iwork(:)
      real, allocatable    :: rwork(:)
      complex, allocatable :: work(:)
      real                 :: eigTemp(hmat%matsize1)
      character(len=6)     :: routine

      call timestart("LAPACK DC GEV")

      n = hmat%matsize1
      if (n /= hmat%matsize2 .or. n /= smat%matsize1 .or. n /= smat%matsize2) &
         call judft_error("Matrices not square in lapack_dc_gev")
      if (ne > n) call judft_error("More eigenvalues requested than matrix size in lapack_dc_gev")
      allocate (t_mat::zmat)
      call zmat%alloc(hmat%l_real, n, ne)
      if (hmat%l_real) then
         routine = "dsygvd"
         call dsygvd(1, 'V', 'U', n, hmat%data_r, size(hmat%data_r, 1), smat%data_r, size(smat%data_r, 1), &
                     eigTemp, rwork_dum, -1, liwork, -1, info)
         call check_info(info, routine//" (workspace query)", ikpt)
         lwork = rwork_dum(1)
         allocate (rwork(lwork), iwork(liwork(1)))
         call dsygvd(1, 'V', 'U', n, hmat%data_r, size(hmat%data_r, 1), smat%data_r, size(smat%data_r, 1), &
                     eigTemp, rwork, lwork, iwork, liwork(1), info)
      else
         routine = "zhegvd"
         call zhegvd(1, 'V', 'U', n, hmat%data_c, size(hmat%data_c, 1), smat%data_c, size(smat%data_c, 1), &
                     eigTemp, work_dum, -1, rwork_dum, -1, liwork, -1, info)
         call check_info(info, routine//" (workspace query)", ikpt)
         lwork = work_dum(1)
         lrwork = rwork_dum(1)
         allocate (work(lwork), rwork(lrwork), iwork(liwork(1)))
         call zhegvd(1, 'V', 'U', n, hmat%data_c, size(hmat%data_c, 1), smat%data_c, size(smat%data_c, 1), &
                     eigTemp, work, lwork, rwork, lrwork, iwork, liwork(1), info)
      end if
      if (info > n) then
         call check_info(info, routine, ikpt, hint="Overlap matrix is not positive definite (leading minor of order " &
                         //int2str(info - n)//"), the basis is probably linearly dependent")
      else
         call check_info(info, routine, ikpt, hint="Divide-and-conquer failed to compute an eigenvalue")
      end if
      ! eigenvectors overwrite hmat, keep the lowest ne
      if (hmat%l_real) then
         zmat%data_r = hmat%data_r(:n, :ne)
      else
         zmat%data_c = hmat%data_c(:n, :ne)
      end if
      eig(:min(size(eig), n)) = eigTemp(:min(size(eig), n))
      call timestop("LAPACK DC GEV")
   end subroutine lapack_dc_gev

   subroutine lapack_dc_diag(self, hmat, ne, eig, zmat)
      ! Standard eigenvalue problem via dsyevd/zheevd
      implicit none
      class(t_solver_lapack_stable)          :: self
      class(t_mat), intent(INOUT)            :: hmat
      integer, intent(INOUT)                 :: ne
      class(t_mat), allocatable, intent(OUT) :: zmat
      real, intent(OUT)                      :: eig(:)

      integer              :: n, info, lwork, lrwork, liwork(1)
      real                 :: rwork_dum(1)
      complex              :: work_dum(1)
      integer, allocatable :: iwork(:)
      real, allocatable    :: rwork(:)
      complex, allocatable :: work(:)
      real                 :: eigTemp(hmat%matsize1)

      call timestart("LAPACK DC STD")

      n = hmat%matsize1
      if (n /= hmat%matsize2) call judft_error("Non-square matrix in lapack_dc_diag")
      if (ne > n) call judft_error("More eigenvalues requested than matrix size in lapack_dc_diag")
      allocate (t_mat::zmat)
      call zmat%alloc(hmat%l_real, n, ne)
      if (hmat%l_real) then
         call dsyevd('V', 'U', n, hmat%data_r, size(hmat%data_r, 1), eigTemp, rwork_dum, -1, liwork, -1, info)
         call check_info(info, "dsyevd (workspace query)")
         lwork = rwork_dum(1)
         allocate (rwork(lwork), iwork(liwork(1)))
         call dsyevd('V', 'U', n, hmat%data_r, size(hmat%data_r, 1), eigTemp, rwork, lwork, iwork, liwork(1), info)
         call check_info(info, "dsyevd", hint="Divide-and-conquer failed to compute an eigenvalue")
         zmat%data_r = hmat%data_r(:n, :ne)
      else
         call zheevd('V', 'U', n, hmat%data_c, size(hmat%data_c, 1), eigTemp, work_dum, -1, rwork_dum, -1, &
                     liwork, -1, info)
         call check_info(info, "zheevd (workspace query)")
         lwork = work_dum(1)
         lrwork = rwork_dum(1)
         allocate (work(lwork), rwork(lrwork), iwork(liwork(1)))
         call zheevd('V', 'U', n, hmat%data_c, size(hmat%data_c, 1), eigTemp, work, lwork, rwork, lrwork, &
                     iwork, liwork(1), info)
         call check_info(info, "zheevd", hint="Divide-and-conquer failed to compute an eigenvalue")
         zmat%data_c = hmat%data_c(:n, :ne)
      end if
      eig(:min(size(eig), n)) = eigTemp(:min(size(eig), n))
      call timestop("LAPACK DC STD")
   end subroutine lapack_dc_diag

   subroutine lapack_dc_diag_sp(self, hmat, ne, eig, zmat)
      ! Standard eigenvalue problem in single precision via ssyevd/cheevd
      implicit none
      class(t_solver_lapack_stable)          :: self
      class(t_mat), intent(INOUT)            :: hmat
      integer, intent(INOUT)                 :: ne
      class(t_mat), allocatable, intent(OUT) :: zmat
      real, intent(OUT)                      :: eig(:)

      integer, parameter   :: sp = selected_real_kind(6)
      integer              :: n, info, lwork, lrwork, liwork(1)
      integer, allocatable :: iwork(:)

      call timestart("LAPACK DC STD-SP")

      n = hmat%matsize1
      if (n /= hmat%matsize2) call judft_error("Non-square matrix in lapack_dc_diag_sp")
      if (ne > n) call judft_error("More eigenvalues requested than matrix size in lapack_dc_diag_sp")
      allocate (t_mat::zmat)
      call zmat%alloc(hmat%l_real, n, ne)

      if (hmat%l_real) then
         block
            real(kind=sp), allocatable :: h(:, :), eigval_sp(:), work(:)
            real(kind=sp)              :: work_dum(1)
            allocate (h(n, n), eigval_sp(n))
            h = hmat%data_r(:n, :n)
            call ssyevd('V', 'U', n, h, n, eigval_sp, work_dum, -1, liwork, -1, info)
            call check_info(info, "ssyevd (workspace query)")
            lwork = work_dum(1)
            allocate (work(lwork), iwork(liwork(1)))
            call ssyevd('V', 'U', n, h, n, eigval_sp, work, lwork, iwork, liwork(1), info)
            call check_info(info, "ssyevd", hint="Divide-and-conquer failed to compute an eigenvalue")
            eig(:ne) = eigval_sp(:ne)
            zmat%data_r = h(:, :ne)
            deallocate (h, eigval_sp, work, iwork)
         end block
      else
         block
            complex(kind=sp), allocatable :: h(:, :), work(:)
            real(kind=sp), allocatable    :: eigval_sp(:), rwork(:)
            complex(kind=sp)              :: work_dum(1)
            real(kind=sp)                 :: rwork_dum(1)
            allocate (h(n, n), eigval_sp(n))
            h = hmat%data_c(:n, :n)
            call cheevd('V', 'U', n, h, n, eigval_sp, work_dum, -1, rwork_dum, -1, liwork, -1, info)
            call check_info(info, "cheevd (workspace query)")
            lwork = work_dum(1)
            lrwork = rwork_dum(1)
            allocate (work(lwork), rwork(lrwork), iwork(liwork(1)))
            call cheevd('V', 'U', n, h, n, eigval_sp, work, lwork, rwork, lrwork, iwork, liwork(1), info)
            call check_info(info, "cheevd", hint="Divide-and-conquer failed to compute an eigenvalue")
            eig(:ne) = eigval_sp(:ne)
            zmat%data_c = h(:, :ne)
            deallocate (h, eigval_sp, work, rwork, iwork)
         end block
      end if
      call timestop("LAPACK DC STD-SP")
   end subroutine lapack_dc_diag_sp

end module m_lapack_stable
