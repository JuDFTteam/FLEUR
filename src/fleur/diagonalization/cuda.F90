!--------------------------------------------------------------------------------
! Copyright (c) 2025 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------
module m_cuda_diag
   use m_types_mat
   use m_types_mpimat
   use m_judft
#ifdef CPP_CUSOLVER
   use cusolverDn
#endif
   use m_types_solver
   implicit none
!**********************************************************
!     Solve the generalized eigenvalue problem
!     using the cusolver library
!**********************************************************
   type, extends(t_solver)::t_solver_cuda
   contains
      procedure        :: solve_gev => cuda_GEV
   end type
   public :: get_solver_cuda

#ifdef CPP_CUSOLVER
   type(cusolverDnHandle)  :: handle
#endif

contains

   function get_solver_cuda() result(solver)
      class(t_solver), allocatable :: solver
      allocate(t_solver_cuda :: solver)
      solver%name = "cuda"
#ifdef CPP_CUSOLVER
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

   subroutine cuda_gev(self, hmat, smat, ne, eig, zmat, ikpt)
    !!Simple driver to solve Generalized Eigenvalue Problem using CuSolverDN
      implicit none
      class(t_solver_cuda) ::self
      class(t_mat), intent(INOUT) :: hmat, smat
      integer, intent(INOUT)      :: ne
      class(t_mat), allocatable, intent(OUT)    :: zmat
      real, intent(OUT)           :: eig(:)
      integer, intent(IN)         :: ikpt

#ifdef CPP_CUSOLVER
      integer                 :: istat, ne_found, lwork_d, devinfo(1)
      real, allocatable        :: work_d(:), eig_tmp(:)
      complex, allocatable     :: work_c(:)

      call timestart("CUDA GEV")

      logical :: firstcall = .true.
      if (firstcall) then
         firstcall = .false.
         istat = cusolverDnCreate(handle)
         if (istat /= CUSOLVER_STATUS_SUCCESS) call judft_error('handle creation failed')
      end if

      allocate (t_mat::zmat)
      allocate (eig_tmp(hmat%matsize1))
      call zmat%alloc(hmat%l_real, hmat%matsize1, ne)
    !!$acc Data copyin(hmat,smat)
      if (hmat%l_real) then
         associate (h => hmat%data_r, s => smat%data_r)
            !$ACC DATA copyin(s)COPY(h)COPYOUT(eig_tmp)
            !$ACC HOST_DATA USE_DEVICE(s,h,eig_tmp)
            istat = cusolverDnDsygvdx_bufferSize(handle, CUSOLVER_EIG_TYPE_1, CUSOLVER_EIG_MODE_VECTOR, CUSOLVER_EIG_RANGE_I, &
                                                 CUBLAS_FILL_MODE_UPPER, hmat%matsize1, h, hmat%matsize1, &
                                                 s, smat%matsize1, 0.0, 0.0, 1, ne, ne_found, eig_tmp, lwork_d)
            !$acc end host_data
            if (istat /= CUSOLVER_STATUS_SUCCESS) call judft_error('cusolverDnZhegvdx_buffersize failed')
            allocate (work_d(lwork_d))
            !$ACC DATA create(work_d,devinfo)
            !$ACC HOST_DATA USE_DEVICE(s,h,eig_tmp,work_d,devinfo)
            istat = cusolverDnDsygvdx(handle, CUSOLVER_EIG_TYPE_1, CUSOLVER_EIG_MODE_VECTOR, CUSOLVER_EIG_RANGE_I, &
                                      CUBLAS_FILL_MODE_UPPER, hmat%matsize1, h, hmat%matsize1, &
                                      s, smat%matsize1, 0.0, 0.0, 1, ne, ne_found, eig_tmp, work_d, lwork_d, devinfo(1))
            !$ACC END HOST_DATA
            !$ACC END DATA
            !$ACC END DATA
            if (istat /= CUSOLVER_STATUS_SUCCESS) call judft_error('cusolverDnZhegvdx failed')
            ne = ne_found
            call zmat%alloc(hmat%l_real, hmat%matsize1, ne_found)
            zmat%data_r = h(:, :ne_found)
            eig = eig_tmp(:ne)
         end associate
      else
         associate (h => hmat%data_c, s => smat%data_c)
            !$ACC DATA copyin(s) COPY(h) COPYOUT(eig_tmp)
            !$ACC HOST_DATA USE_DEVICE(s,h,eig_tmp)
            istat = cusolverDnZhegvdx_bufferSize(handle, CUSOLVER_EIG_TYPE_1, CUSOLVER_EIG_MODE_VECTOR, CUSOLVER_EIG_RANGE_I, &
                                                 CUBLAS_FILL_MODE_UPPER, hmat%matsize1, h, hmat%matsize1, &
                                                 s, smat%matsize1, 0.0, 0.0, 1, ne, ne_found, eig_tmp, lwork_d)
            !$acc end host_data
            if (istat /= CUSOLVER_STATUS_SUCCESS) write (*, *) 'cusolverDnZhegvdx_buffersize failed'
            allocate (work_c(lwork_d))
            !$ACC DATA create(work_c,devinfo)
            !$ACC HOST_DATA USE_DEVICE(s,h,eig_tmp,work_c,devinfo)
            istat = cusolverDnZhegvdx(handle, CUSOLVER_EIG_TYPE_1, CUSOLVER_EIG_MODE_VECTOR, CUSOLVER_EIG_RANGE_I, &
                                      CUBLAS_FILL_MODE_UPPER, hmat%matsize1, h, hmat%matsize1, &
                                      s, smat%matsize1, 0.0, 0.0, 1, ne, ne_found, eig_tmp, work_c, lwork_d, devinfo(1))
            !$ACC END HOST_DATA
            !$acc update self(devinfo)
            if (istat /= CUSOLVER_STATUS_SUCCESS) then
               write (*, *) devinfo
               call judft_error('cusolverDnZhegvdx failed')
            end if
            !$ACC END DATA
            !$ACC END DATA
            ne = ne_found
            call zmat%alloc(hmat%l_real, hmat%matsize1, ne_found)
            zmat%data_c = h(:, :ne_found)
            eig = eig_tmp(:ne)

         end associate
      end if
#endif

      call timestop("CUDA GEV")

   end subroutine

end module m_cuda_diag
