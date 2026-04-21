!--------------------------------------------------------------------------------
! Copyright (c) 2025 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------
module m_nvlamath
   use m_types_solver
   use m_types_mat
   use m_judft
#ifdef CPP_GPU_NVLAMATH      
   use nvlamath
#endif
   implicit none
   private
   type, extends(t_solver)::t_solver_nvlamath
   contains
      !procedure        :: solve_std_dp => nvlamath_diag  !solver for standard eigenvalue problem
      !procedure        :: solve_std_sp => nvlamath_diag_sp  !solver for standard eigenvalue problem
      procedure        :: solve_gev => nvlamath_gev  !solver for generalized eigenvalue problem
      !procedure        :: to_std => nvlamath_reduction     !transform the H of the generalized problem to a std problem
      !procedure        :: backtrans => nvlamath_recover  !transform the Eigenvalue back to the generalized problem
   end type
   public :: t_solver_nvlamath, get_solver_nvlamath

contains

   function get_solver_nvlamath() result(solver)
      class(t_solver), allocatable :: solver
      allocate(t_solver_nvlamath :: solver)
      solver%name = "nvlamath"
#ifdef CPP_GPU_NVLAMATH      
      solver%available = .true.
#else
      solver%available = .false.
#endif      
      solver%parallel = .false.
      solver%serial = .true.
      solver%generalized = .true.
      solver%standard = .false.
      solver%single_precision = .false.
      solver%transform = .false.
      solver%GPU = .true.
      solver%use_sp = .false.
   end function

   subroutine nvlamath_gev(self, hmat, smat, ne, eig, zmat, ikpt)
      !Simple driver to solve Generalized Eigenvalue Problem using nvlamath routine
      implicit none
      class(t_solver_nvlamath)            :: self
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

      allocate (t_mat::zmat)
      call zmat%alloc(hmat%l_real, hmat%matsize1, ne)
      abstol = 2*dlamch('S')
      if (hmat%l_real) then
         allocate (iwork(5*hmat%matsize1), ifail(hmat%matsize1))
         !$acc data create(iwork,ifail,eigTemp,zmat%data_r)
         !$acc host_data use_device(hmat%data_r,smat%data_r,eigTemp,zmat%data_r,dumrwork,iwork,ifail)
         call dsygvx(1, 'V', 'I', 'U', hmat%matsize1, hmat%data_r, size(hmat%data_r, 1), smat%data_r, size(smat%data_r, 1), &
                     0.0, 0.0, 1, ne, abstol, m, eigTemp, zmat%data_r, size(zmat%data_r, 1), dumrwork, -1, iwork, ifail, info)
         !$acc end host_data
         !$acc update self(dumrwork)             
         lwork = dumrwork(1)
         allocate (rwork(lwork))
         !$acc data create(rwork)
         if (info .ne. 0) call judft_error("Diagonalization via nvlamath failed (Workspace)", no=info)
         !$acc host_data use_device(hmat%data_r,smat%data_r,eigTemp,zmat%data_r,rwork,iwork,ifail)
         call dsygvx(1, 'V', 'I', 'U', hmat%matsize1, hmat%data_r, size(hmat%data_r, 1), smat%data_r, size(smat%data_r, 1), &
                     0.0, 0.0, 1, ne, abstol, m, eigTemp, zmat%data_r, size(zmat%data_r, 1), rwork, lwork, iwork, ifail, info)
         !$acc end host_data
         !$acc end data
         !$acc update self(eigTemp,zmat%data_r)            
         !$acc end data            
      else
        
         allocate (rwork(7*hmat%matsize1), iwork(5*hmat%matsize1), ifail(hmat%matsize1))
         !$acc data create(iwork,ifail,eigTemp,zmat%data_c,rwork)
         !Do a workspace query
         !$acc host_data use_device(hmat%data_c,smat%data_c,eigTemp,zmat%data_c,dumwork,rwork,iwork,ifail)
         
         call zhegvx(1, 'V', 'I', 'U', hmat%matsize1, hmat%data_c, size(hmat%data_c, 1), smat%data_c, size(smat%data_c, 1), &
                     0.0, 0.0, 1, ne, abstol, m, eigTemp, zmat%data_c, size(zmat%data_c, 1), dumwork, -1, rwork, iwork, ifail, info)
         !$acc end host_data
         !$acc update self(dumwork)
         lwork = dumwork(1)
         allocate (work(lwork))
         !$acc data create(work)
         if (info .ne. 0) call judft_error("Diagonalization via nvlamath failed (Workspace)", no=info)
         !Perform diagonalization
         !$acc host_data use_device(hmat%data_c,smat%data_c,eigTemp,zmat%data_c,work,rwork,iwork,ifail)
         call zhegvx(1, 'V', 'I', 'U', hmat%matsize1, hmat%data_c, size(hmat%data_c, 1), smat%data_c, size(smat%data_c, 1), &
                      0.0, 0.0, 1, ne, abstol, m, eigTemp, zmat%data_c, size(zmat%data_c, 1), work, lwork, rwork, iwork, ifail, info)
         !$acc end host_data
         !$acc end data
         !$acc update self(eigTemp,zmat%data_c)            
         !$acc end data
      end if
      eig(:min(size(eig), size(eigTemp))) = eigTemp(:min(size(eig), size(eigTemp)))
      if (info .ne. 0) call judft_error("Diagonalization via nvlamath failed(zhegvx/dsygvx)", no=info)
      if (m .ne. ne) call judft_error("Diagonalization via nvlamath failed failed without explicit errorcode.")
   end subroutine nvlamath_gev

end module m_nvlamath
