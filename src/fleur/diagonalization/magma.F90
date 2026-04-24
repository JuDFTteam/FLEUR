!--------------------------------------------------------------------------------
! Copyright (c) 2025 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------

module m_magma
   !! Implementation of a t_solver using the MAGMA library from ICL
   use m_juDFT
   use m_types_solver
   use m_types_mat
#ifdef CPP_MAGMA
   use magma
   use openacc
#endif
   private
   logical ,save :: initialized=.false.
   !integer, save :: Magma_NumGPU = 1
   type, extends(t_solver)::t_solver_magma
   !! provides all solvers& transforms for a "serial" case on the GPU
   contains
      procedure        :: solve_gev => magma_GEV
      procedure        :: solve_std_dp => magma_diag  !solver for standard eigenvalue problem
      procedure        :: solve_std_sp => magma_diag_sp  !solver for standard eigenvalue problem
      procedure        :: to_std => magma_reduction     !transform the H of the generalized problem to a std problem
      procedure        :: backtrans => magma_recover  !transform the Eigenvalue back to the generalized problem
   end type
   public :: get_solver_magma

contains

   function get_solver_magma() result(solver)
      class(t_solver), allocatable :: solver
      allocate(t_solver_magma :: solver)
      solver%name = "magma"
#ifdef CPP_MAGMA
      solver%available = .true.
#else
      solver%available = .false.
#endif
      solver%parallel = .false.
      solver%serial = .true.
      solver%generalized = .true.
      solver%standard = .true.
      solver%single_precision = .true.
      solver%transform = .true.
      solver%GPU = .true.
      solver%use_sp = .false.
   end function

   subroutine init()
#ifdef CPP_MAGMA      
      if (.not. initialized) then
         initialized = .true.
         call magmaf_init()
         call magmaf_setdevice(acc_get_device_num(acc_device_nvidia))
         print *, acc_get_device_num(acc_device_nvidia)
      end if
#endif      
   end subroutine   

   subroutine magma_gev(self, hmat, smat, ne, eig, zmat, ikpt)

      use m_types_mat
      implicit none

      ! ... Arguments ...
      class(t_solver_magma) :: self
      class(t_mat), intent(INOUT)  :: hmat, smat
      integer, intent(INOUT)      :: ne
      class(t_mat), allocatable, intent(OUT)    :: zmat
      real, intent(OUT)           :: eig(:)
      integer, intent(IN)         :: ikpt

#ifdef CPP_MAGMA

      ! ... Local Variables ..
      integer :: lwork, liwork, lrwork, error, mout(1), i
      real    :: eigTemp(hmat%matsize1)
      logical :: initialized = .false.

      real, allocatable :: rwork(:)
      integer, allocatable :: iwork(:)
      complex, allocatable :: work(:)

      call timestart("MAGMA GEV")

      call init()

      if (hmat%l_real) then
         allocate (rwork(1), iwork(1))
         !CALL magmaf_dsygvdx_m(Magma_numGPU,1,'v','i','U',hmat%matsize1,hmat%data_r,SIZE(hmat%data_r,1),smat%data_r,&
         !                    SIZE(smat%data_r,1),0.0,0.0,1,ne,mout,eigTemp,rwork,-1,iwork,-1,error)
         call magmaf_dsygvdx(1, 'v', 'i', 'U', hmat%matsize1, hmat%data_r, size(hmat%data_r, 1), smat%data_r, &
                             size(smat%data_r, 1), 0.0, 0.0, 1, ne, mout, eigTemp, rwork, -1, iwork, -1, error)
         if (error /= 0) then
            write (*, *) 'magmaf_dsygvdx error code: ', error
            call juDFT_error("Failed to query workspaces (1)", calledby="magma.F90")
         end if
         print *, "Magma1"
         lrwork = rwork(1)
         liwork = iwork(1)
         deallocate (rwork, iwork)
         allocate (rwork(lrwork), iwork(liwork))
!       CALL magmaf_dsygvdx_m(Magma_numGPU,1,'v','i','U',hmat%matsize1,hmat%data_r,SIZE(hmat%data_r,1),smat%data_r,&
!                           SIZE(smat%data_r,1),0.0,0.0,1,ne,mout,eigTemp,rwork,lrwork,iwork,liwork,error)
         call magmaf_dsygvdx(1, 'v', 'i', 'U', hmat%matsize1, hmat%data_r, size(hmat%data_r, 1), smat%data_r, &
                             size(smat%data_r, 1), 0.0, 0.0, 1, ne, mout, eigTemp, rwork, lrwork, iwork, liwork, error)
         print *, "Magma2"
         if (error /= 0) then
            write (*, *) 'magmaf_dsygvdx error code: ', error
            call juDFT_error("Magma failed to diagonalize Hamiltonian (1)", calledby="magma.F90")
         end if
      else
         !Query the workspace size
         allocate (work(1), rwork(1), iwork(1))
         !CALL magmaf_zhegvdx_2stage_m(NGPU_CONST,&
         !CALL magmaf_zhegvdx_m(Magma_numGPU,1,'v','i','U',hmat%matsize1,hmat%data_c,SIZE(hmat%data_c,1),smat%data_c,&
         !                    SIZE(smat%data_c,1),0.0,0.0,1,ne,mout,eigTemp,work,-1,rwork,-1,iwork,-1,error)
         call magmaf_zhegvdx(1, 'v', 'i', 'U', hmat%matsize1, hmat%data_c, size(hmat%data_c, 1), smat%data_c, &
                             size(smat%data_c, 1), 0.0, 0.0, 1, ne, mout, eigTemp, work, -1, rwork, -1, iwork, -1, error)
         if (error /= 0) then
            write (*, *) 'magmaf_zhegvdx error code: ', error
            call juDFT_error("Failed to query workspaces (2)", calledby="magma.F90")
         end if
         print *, "Magma1"

         lwork = work(1)
         lrwork = rwork(1)
         liwork = iwork(1)
         deallocate (work, rwork, iwork)
         allocate (work(lwork), rwork(lrwork), iwork(liwork))
         !Now the diagonalization
         !CALL magmaf_zhegvdx_2stage_m(NGPU_CONST,&
         !CALL magmaf_zhegvdx_m(Magma_numGPU,1,'v','i','U',hmat%matsize1,hmat%data_c,SIZE(hmat%data_c,1),smat%data_c,&
         !                    SIZE(smat%data_c,1),0.0,0.0,1,ne,mout,eigTemp,work,lwork,rwork,lrwork,iwork,liwork,error)
         call magmaf_zhegvdx(1, 'v', 'i', 'U', hmat%matsize1, hmat%data_c, size(hmat%data_c, 1), smat%data_c, &
                             size(smat%data_c, 1), 0.0, 0.0, 1, ne, mout, eigTemp, work, lwork, rwork, lrwork, iwork, liwork, error)

         print *, "Magma2"

         if (error /= 0) then
            write (*, *) 'magmaf_zhegvdx error code: ', error
            call juDFT_error("Magma failed to diagonalize Hamiltonian (2)", calledby="magma.F90")
         end if
      end if

      allocate (t_mat::zmat)
      call zmat%alloc(hmat%l_real, hmat%matsize1, ne)
      do i = 1, ne
         eig(i) = eigTemp(i)
         if (hmat%l_real) then
            zmat%data_r(:, i) = hmat%data_r(:, i)
         else
            zmat%data_c(:, i) = hmat%data_c(:, i)
         end if
      end do
#endif
      call timestop("MAGMA GEV")
   end subroutine magma_gev


   subroutine magma_diag(self, hmat, ne, eig, zmat)
      !Simple driver to solve Generalized Eigenvalue Problem using magma routine
      implicit none
      class(t_solver_magma)            :: self
      class(t_mat), intent(INOUT)  :: hmat
      integer, intent(INOUT)      :: ne
      class(t_mat), allocatable, intent(OUT)    :: zmat
      real, intent(OUT)           :: eig(:)

      integer            :: info, m, n
      real               :: abstol
      real, external      :: dlamch
      real               :: eigTemp(hmat%matsize1)

      call timestart("MAGMA STD")
#ifdef CPP_MAGMA
      call init()
      n = hmat%matsize1
      if (n /= hmat%matsize2) call judft_error("Non-square matrix in magma_diag")
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
            call magmaf_dsyevr('V', 'I', 'U', n, hmat%data_r, size(hmat%data_r,1),&
             0.0, 0.0, 1, ne, abstol, m, eigTemp, zmat%data_r, size(zmat%data_r,1), &
                        isuppz, rwork_dum, -1, liwork, -1, info)
            lrwork = rwork_dum(1)
            allocate (rwork(lrwork), iwork(liwork(1)))
            call magmaf_dsyevr('V', 'I', 'U', n, hmat%data_r, size(hmat%data_r,1), &
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
            call magmaf_zheevr('V', 'I', 'U', n, hmat%data_c, size(hmat%data_c,1), &
            0.0, 0.0, 1, ne, abstol, m, eigTemp, zmat%data_c, size(zmat%data_c,1), isuppz, work_dum, &
                        -1, rwork_dum, -1, liwork, -1, info)
            lwork = work_dum(1)
            lrwork = rwork_dum(1)
            allocate (work(lwork), rwork(lrwork), iwork(liwork(1)))
            call magmaf_zheevr('V', 'I', 'U', n, hmat%data_c, size(hmat%data_c,1), &
            0.0, 0.0, 1, ne, abstol, m, eigTemp, zmat%data_c, size(zmat%data_c,1), isuppz, work, &
                        lwork, rwork, lrwork, iwork, liwork(1), info)
         end block
      end if
      eig(:min(size(eig), size(eigTemp))) = eigTemp(:min(size(eig), size(eigTemp)))
#endif      
      call timestop("MAGMA STD")
   end subroutine magma_diag
   subroutine magma_diag_sp(self, hmat, ne, eig, zmat)
      !Simple driver to solve Standard Eigenvalue Problem using magma routine
      implicit none
      class(t_solver_magma)            :: self
      class(t_mat), intent(INOUT)  :: hmat
      integer, intent(INOUT)      :: ne
      class(t_mat), allocatable, intent(OUT)    :: zmat
      real, intent(OUT)           :: eig(:)

      call timestart("MAGMA STD-SP")
#ifdef CPP_MAGMA
      integer, parameter:: sp = selected_real_kind(6)
      integer          :: info, m, n ,lwork
      real(sp)         :: eigval(hmat%matsize1)
      call init()
      n = hmat%matsize1
      if (n /= hmat%matsize2) call judft_error("Non-square matrix in magma_diag")
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
            call magmaf_ssyevx('V','I','U',size(h,1),h,size(h,1),0.0,0.0,1,ne,1.0E-8_sp,m,eigval,z,size(z,1),work,-1,iwork,ifail,info)
            lwork=work(1)
            deallocate(work)
            allocate(work(lwork))
    
            call magmaf_ssyevx('V','I','U',size(h,1),h,size(h,1),0.0,0.0,1,ne,1.0E-8_sp,m,eigval,z,size(z,1),work,lwork,iwork,ifail,info)
            
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
            call magmaf_cheevx('V','I','U',size(h,1),h,size(h,1),0.0,0.0,1,ne,0.0,m,eigval,z,size(z,1),work,-1,rwork,iwork,ifail,info)
            lwork=work(1)
            deallocate(work)
            allocate(work(lwork))
    
            call magmaf_cheevx('V','I','U',size(h,1),h,size(h,1),0.0,0.0,1,ne,0.0,m,eigval,z,size(z,1),work,lwork,rwork,iwork,ifail,info)
            eig=eigval(:ne)
            zmat%data_c=z(:,:ne)
            deallocate(h,z,eigval,work,rwork,iwork)
            END BLOCK   
      end if
#endif      
      call timestop("MAGMA STD-SP")
   end subroutine magma_diag_sp

   subroutine magma_reduction(self, hmat, smat)
      !Simple driver to solve Generalized Eigenvalue Problem using magma routine
      class(t_solver_magma)            :: self
      class(t_mat), intent(INOUT)  :: hmat, smat

      integer            :: info, n

      call timestart("MAGMA REDUCTION")
#ifdef CPP_MAGMA
      call init()
      n = smat%matsize1 !Matrix size
      if (n /= smat%matsize2 .or. n /= hmat%matsize1 .or. n /= hmat%matsize2) &
         call judft_error("Matices not square in magma_reduction")
      if (smat%l_real) then
         ! Perform Cholesky decomposition of B to obtain L (B = L * L^T)
         call magmaf_dpotrf('U', n, smat%data_r, size(smat%data_r,1), info)
         if (info /= 0) call juDFT_error("Error in Cholesky decomposition of B")

         ! Transform A to A' = L^-1 * A * L^-T using chegst
         call magmaf_dsygst(1, "U", n, hmat%data_r, size(hmat%data_r,1), smat%data_r, size(smat%data_r,1), info)
         if (info /= 0) call juDFT_error("Error in dsygst")

      else
         ! Perform Cholesky decomposition of B to obtain L (B = L * L^T)
         call magmaf_zpotrf('U', n, smat%data_c, size(smat%data_c,1), info)
         if (info /= 0) call juDFT_error("Error in Cholesky decomposition of B")

         ! Transform A to A' = L^-1 * A * L^-T using chegst
         call magmaf_zhegst(1, "U", n, hmat%data_c, size(hmat%data_c,1), smat%data_c, size(smat%data_c,1), info)
         if (info /= 0) call juDFT_error("Error in zhegst")
      end if
#endif
      call timestop("MAGMA REDUCTION")
   end subroutine magma_reduction

   subroutine magma_recover(self, smat, zmat)
      class(t_solver_magma)            :: self
      class(t_mat), intent(INOUT)  :: zmat, smat
      integer :: m, n, info

      call timestart("MAGMA BACKTRANSFORM")
#ifdef CPP_MAGMA
      call init()
      n = smat%matsize1
      m = zmat%matsize2
      if (n /= smat%matsize2 .or. n /= zmat%matsize1) call judft_error("Invalid matix sizes in reduction_magma")
      if (smat%l_real) then
         ! recover the generalized eigenvectors z by solving z' = l^t * z
         call magmaf_dtrtrs('U', 'N', 'N', n, m, smat%data_r, n, zMat%data_r, n, info)
         if (info /= 0) call juDFT_error("Error in back transformation (dpotrs)")
      else
         ! --> recover the generalized eigenvectors z by solving z' = l^t * z
         call magmaf_ztrtrs('U', 'N', 'N', n, m, smat%data_c, n, zMat%data_c, n, info)
         if (info /= 0) call juDFT_error("Error in back transformation (zpotrs)")
      end if
#endif      
      call timestop("MAGMA BACKTRANSFORM")
   end subroutine
end module m_magma
