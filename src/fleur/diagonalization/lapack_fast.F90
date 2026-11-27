!--------------------------------------------------------------------------------
! Copyright (c) 2026 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------
module m_lapack_fast
   use m_types_solver
   use m_types_mat
   use m_judft
   implicit none
   private
   type, extends(t_solver)::t_solver_lapack_fast
   contains
      procedure        :: solve_std_dp => lapack_diag  !solver for standard eigenvalue problem
      procedure        :: solve_std_sp => lapack_diag_sp  !solver for standard eigenvalue problem
      procedure        :: solve_gev => lapack_gev  !solver for generalized eigenvalue problem
      procedure        :: to_std => lapack_reduction     !transform the H of the generalized problem to a std problem
      procedure        :: backtrans => lapack_recover  !transform the Eigenvalue back to the generalized problem
   end type
   public :: t_solver_lapack_fast, get_solver_lapack_fast, check_info

contains

   function get_solver_lapack_fast() result(solver)
      class(t_solver), allocatable :: solver
      allocate (t_solver_lapack_fast :: solver)
      solver%name = "lapack_fast"
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

   subroutine lapack_gev(self, hmat, smat, ne, eig, zmat, ikpt)
      !Simple driver to solve Generalized Eigenvalue Problem using LAPACK routine
      implicit none
      class(t_solver_lapack_fast)            :: self
      class(t_mat), intent(INOUT)  :: hmat, smat
      integer, intent(INOUT)      :: ne
      class(t_mat), allocatable, intent(OUT)    :: zmat
      real, intent(OUT)           :: eig(:)
      integer, intent(IN)         :: ikpt

      integer            :: lwork, info, m
      integer, allocatable:: ifail(:), iwork(:)
      complex, allocatable:: work(:)
      real, allocatable   :: rwork(:)
      real               :: dumrwork(1), abstol
      complex            :: dumwork(1)
      real, external      :: dlamch
      real               :: eigTemp(hmat%matsize1)
      character(len=6)   :: routine

      call timestart("LAPACK GEV")

      allocate (t_mat::zmat)
      call zmat%alloc(hmat%l_real, hmat%matsize1, ne)
      abstol = 2*dlamch('S')
      if (hmat%l_real) then
         routine = "dsygvx"
         allocate (iwork(5*hmat%matsize1), ifail(hmat%matsize1))
         call dsygvx(1, 'V', 'I', 'U', hmat%matsize1, hmat%data_r, size(hmat%data_r, 1), smat%data_r, size(smat%data_r, 1), &
                     0.0, 0.0, 1, ne, abstol, m, eigTemp, zmat%data_r, size(zmat%data_r, 1), dumrwork, -1, iwork, ifail, info)
         call check_info(info, routine//" (workspace query)", ikpt)
         lwork = dumrwork(1)
         allocate (rwork(lwork))
         call dsygvx(1, 'V', 'I', 'U', hmat%matsize1, hmat%data_r, size(hmat%data_r, 1), smat%data_r, size(smat%data_r, 1), &
                     0.0, 0.0, 1, ne, abstol, m, eigTemp, zmat%data_r, size(zmat%data_r, 1), rwork, lwork, iwork, ifail, info)
      else
         routine = "zhegvx"
         allocate (rwork(7*hmat%matsize1), iwork(5*hmat%matsize1), ifail(hmat%matsize1))
         !Do a workspace query
         call zhegvx(1, 'V', 'I', 'U', hmat%matsize1, hmat%data_c, size(hmat%data_c, 1), smat%data_c, size(smat%data_c, 1), &
                     0.0, 0.0, 1, ne, abstol, m, eigTemp, zmat%data_c, size(zmat%data_c, 1), dumwork, -1, rwork, iwork, ifail, info)
         call check_info(info, routine//" (workspace query)", ikpt)
         lwork = dumwork(1)
         allocate (work(lwork))
         !Perform diagonalization
         call zhegvx(1, 'V', 'I', 'U', hmat%matsize1, hmat%data_c, size(hmat%data_c, 1), smat%data_c, size(smat%data_c, 1), &
                     0.0, 0.0, 1, ne, abstol, m, eigTemp, zmat%data_c, size(zmat%data_c, 1), work, lwork, rwork, iwork, ifail, info)
      end if
      if (info > hmat%matsize1) then
         call check_info(info, routine, ikpt, hint="Overlap matrix is not positive definite (leading minor of order " &
                         //int2str(info - hmat%matsize1)//"), the basis is probably linearly dependent")
      else
         call check_info(info, routine, ikpt, hint=int2str(info)//" eigenvectors failed to converge")
      end if
      call check_count(m, ne, routine, ikpt)
      eig(:min(size(eig), size(eigTemp))) = eigTemp(:min(size(eig), size(eigTemp)))
      call timestop("LAPACK GEV")
   end subroutine lapack_gev

   subroutine lapack_diag(self, hmat, ne, eig, zmat)
      !Simple driver to solve Generalized Eigenvalue Problem using LAPACK routine
      implicit none
      class(t_solver_lapack_fast)            :: self
      class(t_mat), intent(INOUT)  :: hmat
      integer, intent(INOUT)      :: ne
      class(t_mat), allocatable, intent(OUT)    :: zmat
      real, intent(OUT)           :: eig(:)

      integer            :: info, m, n
      real               :: abstol
      real, external      :: dlamch
      real               :: eigTemp(hmat%matsize1)

      call timestart("LAPACK STD")

      n = hmat%matsize1
      if (n /= hmat%matsize2) call judft_error("Non-square matrix in lapack_diag")
      allocate (t_mat::zmat)
      call zmat%alloc(hmat%l_real, n, ne)
      abstol = 2*dlamch('S')
      if (hmat%l_real) then
         block  !workspace locally
            integer:: isuppz(2*n), lrwork, liwork(1)
            real   :: rwork_dum(1)
            real, allocatable     :: rwork(:)
            integer, allocatable  :: iwork(:)
            ! Workspace query
            call dsyevr('V', 'I', 'U', n, hmat%data_r, size(hmat%data_r, 1), &
                        0.0, 0.0, 1, ne, abstol, m, eigTemp, zmat%data_r, size(zmat%data_r, 1), &
                        isuppz, rwork_dum, -1, liwork, -1, info)
            call check_info(info, "dsyevr (workspace query)")
            lrwork = rwork_dum(1)
            allocate (rwork(lrwork), iwork(liwork(1)))
            call dsyevr('V', 'I', 'U', n, hmat%data_r, size(hmat%data_r, 1), &
                        0.0, 0.0, 1, ne, abstol, m, eigTemp, zmat%data_r, size(zmat%data_r, 1), &
                        isuppz, rwork, lrwork, iwork, liwork(1), info)
            call check_info(info, "dsyevr", hint="Internal error in the MRRR algorithm")
            call check_count(m, ne, "dsyevr")
         end block
      else
         block  !workspace locally
            integer:: isuppz(2*n), lwork, lrwork, liwork(1)
            complex:: work_dum(1)
            real   :: rwork_dum(1)
            complex, allocatable  :: work(:)
            real, allocatable     :: rwork(:)
            integer, allocatable  :: iwork(:)
            ! Workspace query
            call zheevr('V', 'I', 'U', n, hmat%data_c, size(hmat%data_c, 1), &
                        0.0, 0.0, 1, ne, abstol, m, eigTemp, zmat%data_c, size(zmat%data_c, 1), isuppz, work_dum, &
                        -1, rwork_dum, -1, liwork, -1, info)
            call check_info(info, "zheevr (workspace query)")
            lwork = work_dum(1)
            lrwork = rwork_dum(1)
            allocate (work(lwork), rwork(lrwork), iwork(liwork(1)))
            call zheevr('V', 'I', 'U', n, hmat%data_c, size(hmat%data_c, 1), &
                        0.0, 0.0, 1, ne, abstol, m, eigTemp, zmat%data_c, size(zmat%data_c, 1), isuppz, work, &
                        lwork, rwork, lrwork, iwork, liwork(1), info)
            call check_info(info, "zheevr", hint="Internal error in the MRRR algorithm")
            call check_count(m, ne, "zheevr")
         end block
      end if
      eig(:min(size(eig), size(eigTemp))) = eigTemp(:min(size(eig), size(eigTemp)))
      call timestop("LAPACK STD")
   end subroutine lapack_diag

   subroutine lapack_diag_sp(self, hmat, ne, eig, zmat)
      !Simple driver to solve Standard Eigenvalue Problem using LAPACK routine
      implicit none
      class(t_solver_lapack_fast)            :: self
      class(t_mat), intent(INOUT)  :: hmat
      integer, intent(INOUT)      :: ne
      class(t_mat), allocatable, intent(OUT)    :: zmat
      real, intent(OUT)           :: eig(:)

      integer, parameter:: sp = selected_real_kind(6)
      integer          :: info, m, n, lwork

      call timestart("LAPACK STD-SP")

      n = hmat%matsize1
      if (n /= hmat%matsize2) call judft_error("Non-square matrix in lapack_diag")
      allocate (t_mat::zmat)
      call zmat%alloc(hmat%l_real, n, ne)

      if (hmat%l_real) then
         BLOCK
            REAL(kind=sp), allocatable:: h(:, :), z(:, :), eigval_sp(:), work(:)
            integer, allocatable      :: iwork(:), ifail(:)
            Allocate (h(size(hmat%data_r, 1), size(hmat%data_r, 2)))
            Allocate (eigval_sp(size(hmat%data_r, 1)), ifail(size(hmat%data_r, 1)))
            Allocate (z(size(hmat%data_r, 1), ne))
            h = hmat%data_r

            allocate (work(1), iwork(5*size(h, 1)))
            call ssyevx('V','I','U',size(h,1),h,size(h,1),0.0,0.0,1,ne,1.0E-8_sp,m,eigval_sp,z,size(z,1),work,-1,iwork,ifail,info)
            call check_info(info, "ssyevx (workspace query)")
            lwork = work(1)
            deallocate (work)
            allocate (work(lwork))

           call ssyevx('V','I','U',size(h,1),h,size(h,1),0.0,0.0,1,ne,1.0E-8_sp,m,eigval_sp,z,size(z,1),work,lwork,iwork,ifail,info)
            call check_info(info, "ssyevx", hint=int2str(info)//" eigenvectors failed to converge")
            call check_count(m, ne, "ssyevx")

            eig(:ne) = eigval_sp(:ne)
            zmat%data_r = z(:, :ne)
            deallocate (h, z, eigval_sp, work, iwork)
         END BLOCK
      else
         BLOCK
            COMPLEX(kind=sp), allocatable:: h(:, :), z(:, :), work(:)
            REAL(kind=sp), allocatable:: eigval_sp(:), rwork(:)
            integer, allocatable      :: iwork(:), ifail(:)
            Allocate (h(size(hmat%data_c, 1), size(hmat%data_c, 2)))
            Allocate (eigval_sp(size(hmat%data_c, 1)), ifail(size(hmat%data_c, 1)))
            Allocate (z(size(hmat%data_c, 1), ne), rwork(7*size(hmat%data_c, 1)))
            h = hmat%data_c

            allocate (work(1), iwork(5*size(hmat%data_c, 1)))
            call cheevx('V','I','U',size(h,1),h,size(h,1),0.0,0.0,1,ne,0.0,m,eigval_sp,z,size(z,1),work,-1,rwork,iwork,ifail,info)
            call check_info(info, "cheevx (workspace query)")
            lwork = work(1)
            deallocate (work)
            allocate (work(lwork))

           call cheevx('V','I','U',size(h,1),h,size(h,1),0.0,0.0,1,ne,0.0,m,eigval_sp,z,size(z,1),work,lwork,rwork,iwork,ifail,info)
            call check_info(info, "cheevx", hint=int2str(info)//" eigenvectors failed to converge")
            call check_count(m, ne, "cheevx")
            eig(:ne) = eigval_sp(:ne)
            zmat%data_c = z(:, :ne)
            deallocate (h, z, eigval_sp, work, rwork, iwork)
         END BLOCK
      end if
      call timestop("LAPACK STD-SP")
   end subroutine lapack_diag_sp

   subroutine lapack_reduction(self, hmat, smat, ne)
      !Simple driver to solve Generalized Eigenvalue Problem using LAPACK routine
      class(t_solver_lapack_fast)            :: self
      class(t_mat), intent(INOUT)  :: hmat, smat
      integer, intent(IN)  :: ne

      integer            :: info, n

      call timestart("LAPACK REDUCTION")

      n = smat%matsize1 !Matrix size
      if (n /= smat%matsize2 .or. n /= hmat%matsize1 .or. n /= hmat%matsize2) &
         call judft_error("Matices not square in lapack_reduction")
      if (smat%l_real) then
         ! Perform Cholesky decomposition of B to obtain L (B = L * L^T)
         call dpotrf('U', n, smat%data_r, size(smat%data_r, 1), info)
         call check_info(info, "dpotrf", hint=posdef_hint(info))

         ! Transform A to A' = L^-1 * A * L^-T using chegst
         call dsygst(1, "U", n, hmat%data_r, size(hmat%data_r, 1), smat%data_r, size(smat%data_r, 1), info)
         call check_info(info, "dsygst")

      else
         ! Perform Cholesky decomposition of B to obtain L (B = L * L^T)
         call zpotrf('U', n, smat%data_c, size(smat%data_c, 1), info)
         call check_info(info, "zpotrf", hint=posdef_hint(info))

         ! Transform A to A' = L^-1 * A * L^-T using chegst
         call zhegst(1, "U", n, hmat%data_c, size(hmat%data_c, 1), smat%data_c, size(smat%data_c, 1), info)
         call check_info(info, "zhegst")
      end if

      call timestop("LAPACK REDUCTION")

   end subroutine lapack_reduction

   subroutine lapack_recover(self, smat, zmat)
      class(t_solver_lapack_fast)            :: self
      class(t_mat), intent(INOUT)  :: zmat, smat
      integer :: m, n, info

      call timestart("LAPACK BACKTRANSFORM")

      n = smat%matsize1
      m = zmat%matsize2
      if (n /= smat%matsize2 .or. n /= zmat%matsize1) call judft_error("Invalid matix sizes in reduction_lapack")
      if (smat%l_real) then
         ! recover the generalized eigenvectors z by solving z' = l^t * z
         call dtrtrs('U', 'N', 'N', n, m, smat%data_r, n, zMat%data_r, n, info)
         call check_info(info, "dtrtrs", hint=singular_hint(info))
      else
         ! --> recover the generalized eigenvectors z by solving z' = l^t * z
         call ztrtrs('U', 'N', 'N', n, m, smat%data_c, n, zMat%data_c, n, info)
         call check_info(info, "ztrtrs", hint=singular_hint(info))
      end if
      call timestop("LAPACK BACKTRANSFORM")
   end subroutine

   subroutine check_info(info, routine, ikpt, hint)
      ! Abort if a LAPACK routine returned info /= 0
      integer, intent(in)                    :: info
      character(len=*), intent(in)           :: routine
      integer, intent(in), optional          :: ikpt
      character(len=*), intent(in), optional :: hint

      character(len=:), allocatable :: msg

      if (info == 0) return
      if (info < 0) then
         msg = "LAPACK "//routine//": argument "//int2str(-info)//" had an illegal value"
      else
         msg = "LAPACK "//routine//" failed with info="//int2str(info)
      end if
      if (present(ikpt)) msg = msg//" for k-point "//int2str(ikpt)
      if (info < 0) then
         call judft_error(msg, calledby="m_lapack", bug=.true.)
      else
         call judft_error(msg, calledby="m_lapack", hint=hint)
      end if
   end subroutine check_info

   subroutine check_count(m, ne, routine, ikpt)
      ! Abort if fewer eigenpairs than requested were found
      integer, intent(in)           :: m, ne
      character(len=*), intent(in)  :: routine
      integer, intent(in), optional :: ikpt

      character(len=:), allocatable :: msg

      if (m == ne) return
      msg = "LAPACK "//routine//" found "//int2str(m)//" eigenvalues instead of "//int2str(ne)
      if (present(ikpt)) msg = msg//" for k-point "//int2str(ikpt)
      call judft_error(msg, calledby="m_lapack")
   end subroutine check_count

   function posdef_hint(info) result(hint)
      integer, intent(in)           :: info
      character(len=:), allocatable :: hint
      hint = "Overlap matrix is not positive definite (leading minor of order "//int2str(info)// &
             "), the basis is probably linearly dependent"
   end function posdef_hint

   function singular_hint(info) result(hint)
      integer, intent(in)           :: info
      character(len=:), allocatable :: hint
      hint = "Cholesky factor is singular (diagonal element "//int2str(info)//" is zero)"
   end function singular_hint

end module m_lapack_fast
