!--------------------------------------------------------------------------------
! Copyright (c) 2025 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------
module m_lapack
   use m_types_solver
   use m_types_mat
   use m_judft
   implicit none
   private
   type, extends(t_solver)::t_solver_lapack
   contains
   procedure        :: solve_std_dp => lapack_diag  !solver for standard eigenvalue problem
   procedure        :: solve_std_sp => lapack_diag_sp  !solver for standard eigenvalue problem
   procedure        :: solve_gev => lapack_gev  !solver for generalized eigenvalue problem
   procedure        :: to_std => lapack_reduction     !transform the H of the generalized problem to a std problem
   procedure        :: backtrans => lapack_recover  !transform the Eigenvalue back to the generalized problem
   end type
   public :: t_solver_lapack, get_solver_lapack

contains

   function get_solver_lapack() result(solver)
      class(t_solver), allocatable :: solver
      allocate(t_solver_lapack :: solver)
      solver%name = "lapack"
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
      class(t_solver_lapack)            :: self
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

      call timestart("LAPACK GEV")

      allocate (t_mat::zmat)
      call zmat%alloc(hmat%l_real, hmat%matsize1, ne)
      abstol = 2*dlamch('S')
      if (hmat%l_real) then
         allocate (iwork(5*hmat%matsize1), ifail(hmat%matsize1))
         call dsygvx(1, 'V', 'I', 'U', hmat%matsize1, hmat%data_r, size(hmat%data_r, 1), smat%data_r, size(smat%data_r, 1), &
                     0.0, 0.0, 1, ne, abstol, m, eigTemp, zmat%data_r, size(zmat%data_r, 1), dumrwork, -1, iwork, ifail, info)
         lwork = dumrwork(1)
         allocate (rwork(lwork))
         if (info .ne. 0) THEN
            WRITE(*,*) 'Error for k-point ', ikpt
            call judft_error("Diagonalization via LAPACK failed (Workspace)", no=info)
         END IF
         call dsygvx(1, 'V', 'I', 'U', hmat%matsize1, hmat%data_r, size(hmat%data_r, 1), smat%data_r, size(smat%data_r, 1), &
                     0.0, 0.0, 1, ne, abstol, m, eigTemp, zmat%data_r, size(zmat%data_r, 1), rwork, lwork, iwork, ifail, info)
      else
         allocate (rwork(7*hmat%matsize1), iwork(5*hmat%matsize1), ifail(hmat%matsize1))
         !Do a workspace query
         call zhegvx(1, 'V', 'I', 'U', hmat%matsize1, hmat%data_c, size(hmat%data_c, 1), smat%data_c, size(smat%data_c, 1), &
                     0.0, 0.0, 1, ne, abstol, m, eigTemp, zmat%data_c, size(zmat%data_c, 1), dumwork, -1, rwork, iwork, ifail, info)
         lwork = dumwork(1)
         allocate (work(lwork))
         if (info .ne. 0) THEN
            WRITE(*,*) 'Error for k-point ', ikpt
            call judft_error("Diagonalization via LAPACK failed (Workspace)", no=info)
         END IF
         !Perform diagonalization
         call zhegvx(1, 'V', 'I', 'U', hmat%matsize1, hmat%data_c, size(hmat%data_c, 1), smat%data_c, size(smat%data_c, 1), &
                     0.0, 0.0, 1, ne, abstol, m, eigTemp, zmat%data_c, size(zmat%data_c, 1), work, lwork, rwork, iwork, ifail, info)
      end if
      eig(:min(size(eig), size(eigTemp))) = eigTemp(:min(size(eig), size(eigTemp)))
      if (info .ne. 0) THEN
         WRITE(*,*) 'Error for k-point ', ikpt
         call judft_error("Diagonalization via LAPACK failed(zhegvx/dsygvx)", no=info)
      END IF
      if (m .ne. ne) THEN
         WRITE(*,*) 'Error for k-point ', ikpt
         call judft_error("Diagonalization via LAPACK failed failed without explicit errorcode.")
      END IF
      call timestop("LAPACK GEV")
   end subroutine lapack_gev

   subroutine lapack_diag(self, hmat, ne, eig, zmat)
      !Simple driver to solve Generalized Eigenvalue Problem using LAPACK routine
      implicit none
      class(t_solver_lapack)            :: self
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
            call dsyevr('V', 'I', 'U', n, hmat%data_r, size(hmat%data_r,1),&
             0.0, 0.0, 1, ne, abstol, m, eigTemp, zmat%data_r, size(zmat%data_r,1), &
                        isuppz, rwork_dum, -1, liwork, -1, info)
            lrwork = rwork_dum(1)
            allocate (rwork(lrwork), iwork(liwork(1)))
            call dsyevr('V', 'I', 'U', n, hmat%data_r, size(hmat%data_r,1), &
            0.0, 0.0, 1, ne, abstol, m, eigTemp, zmat%data_r, size(zmat%data_r,1), &
                        isuppz, rwork, lrwork, iwork, liwork(1), info)
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
            call zheevr('V', 'I', 'U', n, hmat%data_c, size(hmat%data_c,1), &
            0.0, 0.0, 1, ne, abstol, m, eigTemp, zmat%data_c, size(zmat%data_c,1), isuppz, work_dum, &
                        -1, rwork_dum, -1, liwork, -1, info)
            lwork = work_dum(1)
            lrwork = rwork_dum(1)
            allocate (work(lwork), rwork(lrwork), iwork(liwork(1)))
            call zheevr('V', 'I', 'U', n, hmat%data_c, size(hmat%data_c,1), &
            0.0, 0.0, 1, ne, abstol, m, eigTemp, zmat%data_c, size(zmat%data_c,1), isuppz, work, &
                        lwork, rwork, lrwork, iwork, liwork(1), info)
         end block
      end if
      eig(:min(size(eig), size(eigTemp))) = eigTemp(:min(size(eig), size(eigTemp)))
      call timestop("LAPACK STD")
   end subroutine lapack_diag

   subroutine lapack_diag_sp(self, hmat, ne, eig, zmat)
      !Simple driver to solve Standard Eigenvalue Problem using LAPACK routine
      implicit none
      class(t_solver_lapack)            :: self
      class(t_mat), intent(INOUT)  :: hmat
      integer, intent(INOUT)      :: ne
      class(t_mat), allocatable, intent(OUT)    :: zmat
      real, intent(OUT)           :: eig(:)

      integer, parameter:: sp = selected_real_kind(6)
      integer          :: info, m, n ,lwork
      real(sp)         :: eigval(hmat%matsize1)

      call timestart("LAPACK STD-SP")

      n = hmat%matsize1
      if (n /= hmat%matsize2) call judft_error("Non-square matrix in lapack_diag")
      allocate (t_mat::zmat)
      call zmat%alloc(hmat%l_real, n, ne)

      if (hmat%l_real) then
         BLOCK
            REAL(kind=sp),allocatable:: h(:,:),z(:,:),eigval(:),work(:)
            integer,allocatable      :: iwork(:),ifail(:)
            Allocate(h(size(hmat%data_r,1),size(hmat%data_r,2)))
            Allocate(eigval(size(hmat%data_r,1)),ifail(size(hmat%data_r,1)))
            Allocate(z(size(hmat%data_r,1),ne))
            h=hmat%data_r
    
            allocate(work(1),iwork(5*size(h,1)))
            call ssyevx('V','I','U',size(h,1),h,size(h,1),0.0,0.0,1,ne,1.0E-8_sp,m,eigval,z,size(z,1),work,-1,iwork,ifail,info)
            lwork=work(1)
            deallocate(work)
            allocate(work(lwork))
    
            call ssyevx('V','I','U',size(h,1),h,size(h,1),0.0,0.0,1,ne,1.0E-8_sp,m,eigval,z,size(z,1),work,lwork,iwork,ifail,info)
            
            eig(:ne)=eigval(:ne)
            zmat%data_r=z(:,:ne)
            deallocate(h,z,eigval,work,iwork)
           END BLOCK
      else
         BLOCK
            COMPLEX(kind=sp),allocatable:: h(:,:),z(:,:),work(:)
            REAL(kind=sp),allocatable:: eigval(:),rwork(:)
            integer,allocatable      :: iwork(:),ifail(:)
            Allocate(h(size(hmat%data_c,1),size(hmat%data_c,2)))
            Allocate(eigval(size(hmat%data_c,1)),ifail(size(hmat%data_c,1)))
            Allocate(z(size(hmat%data_c,1),ne),rwork(7*size(hmat%data_c,1)))
            h=hmat%data_c
    
            allocate(work(1),iwork(5*size(hmat%data_c,1)))
            call cheevx('V','I','U',size(h,1),h,size(h,1),0.0,0.0,1,ne,0.0,m,eigval,z,size(z,1),work,-1,rwork,iwork,ifail,info)
            lwork=work(1)
            deallocate(work)
            allocate(work(lwork))
    
            call cheevx('V','I','U',size(h,1),h,size(h,1),0.0,0.0,1,ne,0.0,m,eigval,z,size(z,1),work,lwork,rwork,iwork,ifail,info)
            eig=eigval(:ne)
            zmat%data_c=z(:,:ne)
            deallocate(h,z,eigval,work,rwork,iwork)
            END BLOCK   
      end if
      call timestop("LAPACK STD-SP")
   end subroutine lapack_diag_sp

   subroutine lapack_reduction(self, hmat, smat)
      !Simple driver to solve Generalized Eigenvalue Problem using LAPACK routine
      class(t_solver_lapack)            :: self
      class(t_mat), intent(INOUT)  :: hmat, smat

      integer            :: info, n

      call timestart("LAPACK REDUCTION")

      n = smat%matsize1 !Matrix size
      if (n /= smat%matsize2 .or. n /= hmat%matsize1 .or. n /= hmat%matsize2) &
         call judft_error("Matices not square in lapack_reduction")
      if (smat%l_real) then
         ! Perform Cholesky decomposition of B to obtain L (B = L * L^T)
         call dpotrf('U', n, smat%data_r, size(smat%data_r,1), info)
         if (info /= 0) call juDFT_error("Error in Cholesky decomposition of B")

         ! Transform A to A' = L^-1 * A * L^-T using chegst
         call dsygst(1, "U", n, hmat%data_r, size(hmat%data_r,1), smat%data_r, size(smat%data_r,1), info)
         if (info /= 0) call juDFT_error("Error in dsygst")

      else
         ! Perform Cholesky decomposition of B to obtain L (B = L * L^T)
         call zpotrf('U', n, smat%data_c, size(smat%data_c,1), info)
         if (info /= 0) call juDFT_error("Error in Cholesky decomposition of B")

         ! Transform A to A' = L^-1 * A * L^-T using chegst
         call zhegst(1, "U", n, hmat%data_c, size(hmat%data_c,1), smat%data_c, size(smat%data_c,1), info)
         if (info /= 0) call juDFT_error("Error in zhegst")
      end if

      call timestop("LAPACK REDUCTION")

   end subroutine lapack_reduction

   subroutine lapack_recover(self, smat, zmat)
      class(t_solver_lapack)            :: self
      class(t_mat), intent(INOUT)  :: zmat, smat
      integer :: m, n, info

      call timestart("LAPACK BACKTRANSFORM")

      n = smat%matsize1
      m = zmat%matsize2
      if (n /= smat%matsize2 .or. n /= zmat%matsize1) call judft_error("Invalid matix sizes in reduction_lapack")
      if (smat%l_real) then
         ! recover the generalized eigenvectors z by solving z' = l^t * z
         call dtrtrs('U', 'N', 'N', n, m, smat%data_r, n, zMat%data_r, n, info)
         if (info /= 0) call juDFT_error("Error in back transformation (dpotrs)")
      else
         ! --> recover the generalized eigenvectors z by solving z' = l^t * z
         call ztrtrs('U', 'N', 'N', n, m, smat%data_c, n, zMat%data_c, n, info)
         if (info /= 0) call juDFT_error("Error in back transformation (zpotrs)")
      end if
      call timestop("LAPACK BACKTRANSFORM")
   end subroutine

end module m_lapack
