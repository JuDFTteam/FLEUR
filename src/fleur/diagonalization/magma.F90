!--------------------------------------------------------------------------------
! Copyright (c) 2016 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------

module m_magma
   use m_juDFT
   use m_types_solver
   integer, save :: Magma_NumGPU = 1
   !**********************************************************
   !     Solve the generalized eigenvalue problem
   !     using the MAGMA library for multiple GPUs
   !**********************************************************
   type, extends(t_solver)::t_solver_magma
   contains
      procedure        :: solve_gev => magma_GEV
   end type
   public :: get_solver_magma

contains

   function get_solver_magma() result(solver)
      type(t_solver_magma), pointer::solver
      allocate (solver)
      solver%name = "magma"
#ifdef CPP_MAGMA
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
   end function

   subroutine magma_gev(self, hmat, smat, ne, eig, zmat)
#ifdef CPP_MAGMA
      use magma
      use openacc
#endif
      use m_types_mat
      implicit none

      ! ... Arguments ...
      class(t_solver_magma) :: self
      class(t_mat), intent(INOUT)  :: hmat, smat
      integer, intent(INOUT)      :: ne
      class(t_mat), allocatable, intent(OUT)    :: zmat
      real, intent(OUT)           :: eig(:)

#ifdef CPP_MAGMA

      ! ... Local Variables ..
      integer :: lwork, liwork, lrwork, error, mout(1), i
      real    :: eigTemp(hmat%matsize1)
      logical :: initialized = .false.

      real, allocatable :: rwork(:)
      integer, allocatable :: iwork(:)
      complex, allocatable :: work(:)

      if (.not. initialized) then
         initialized = .true.
         call magmaf_init()
         call magmaf_setdevice(acc_get_device_num(acc_device_nvidia))
         print *, acc_get_device_num(acc_device_nvidia)
      end if

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
   end subroutine magma_gev
end module m_magma
