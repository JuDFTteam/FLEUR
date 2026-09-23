!--------------------------------------------------------------------------------
! Copyright (c) 2026 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------
module m_dfpt_gga_kernel
   !! Region independent assembly of the GGA part of the xc potential response
   !! v_xc1 = driv1 - div(H) from the libxc kernels drivsigma, driv2rho2,
   !! driv2rhosigma and driv2sigma2. The divergence is expanded by the product
   !! rule: the terms containing a Laplacian are done in dfpt_gga_assemble, the
   !! remaining two grad.grad terms in dfpt_gga_grdotgr once the caller has
   !! supplied the gradients of drivsigma and drivsigma1.
   !!
   !! Everything is linear in the response, so the real and the imaginary
   !! channel of rho1 are two independent calls sharing the same kernels.

   implicit none
   private

   ! Symmetric packing of driv2sigma2, libxc order [11,12,13,22,23,33]
   integer, parameter :: sigmaIdx(3,3) = reshape([1,2,3, 2,4,5, 3,5,6],[3,3])

   public :: dfpt_gga_assemble, dfpt_gga_grdotgr

contains

   subroutine dfpt_gga_assemble(drivsigma, driv2rho2, driv2rhosigma, driv2sigma2, gradRho, gradRho1, rho1, v_xc1, drivsigma1)
      !! Local part of the GGA potential response for one channel. Returns v_xc1
      !! including the Laplacian terms of -div(H), and the response of drivsigma,
      !! whose gradient the caller has to feed back through dfpt_gga_grdotgr.
      use m_types

      implicit none

      real,              intent(in)  :: drivsigma(:,:)     ! (n_sigma, points)
      real,              intent(in)  :: driv2rho2(:,:)     ! (2*jspins-1, points)
      real,              intent(in)  :: driv2rhosigma(:,:) ! (jspins*n_sigma, points)
      real,              intent(in)  :: driv2sigma2(:,:)   ! (1 or 6, points)
      type(t_gradients), intent(in)  :: gradRho, gradRho1
      real,              intent(in)  :: rho1(:,:)          ! (points, jspins)
      real,              intent(out) :: v_xc1(:,:)         ! (points, jspins)
      real,              intent(out) :: drivsigma1(:,:)    ! (points, n_sigma)

      integer :: ipt, isig, jsig ! sig counts libxc convention of stored (up,up ; up,down ; down;down)
      real    :: sigma1(3)

      v_xc1 = 0.0
      drivsigma1 = 0.0

      if (size(rho1,2) == 1) then
         do ipt = 1, size(rho1,1)
            sigma1(1) = 2.0*dot_product(gradRho%gr(:,ipt,1),gradRho1%gr(:,ipt,1))

            v_xc1(ipt,1) = driv2rho2(1,ipt)*rho1(ipt,1) + driv2rhosigma(1,ipt)*sigma1(1)
            drivsigma1(ipt,1) = driv2rhosigma(1,ipt)*rho1(ipt,1) + driv2sigma2(1,ipt)*sigma1(1)

            v_xc1(ipt,1) = v_xc1(ipt,1) - 2.0*drivsigma1(ipt,1)*gradRho%laplace(ipt,1)
            v_xc1(ipt,1) = v_xc1(ipt,1) - 2.0*drivsigma(1,ipt)*gradRho1%laplace(ipt,1)
         end do
      else
         do ipt = 1, size(rho1,1)
            sigma1(1) = 2.0*dot_product(gradRho%gr(:,ipt,1),gradRho1%gr(:,ipt,1))
            sigma1(2) = dot_product(gradRho%gr(:,ipt,1),gradRho1%gr(:,ipt,2)) + dot_product(gradRho%gr(:,ipt,2),gradRho1%gr(:,ipt,1))
            sigma1(3) = 2.0*dot_product(gradRho%gr(:,ipt,2),gradRho1%gr(:,ipt,2))

            v_xc1(ipt,1) = driv2rho2(1,ipt)*rho1(ipt,1) + driv2rho2(2,ipt)*rho1(ipt,2)
            v_xc1(ipt,2) = driv2rho2(2,ipt)*rho1(ipt,1) + driv2rho2(3,ipt)*rho1(ipt,2)

            do isig = 1, 3
               v_xc1(ipt,1) = v_xc1(ipt,1) + driv2rhosigma(isig,ipt)*sigma1(isig)
               v_xc1(ipt,2) = v_xc1(ipt,2) + driv2rhosigma(3+isig,ipt)*sigma1(isig)
               drivsigma1(ipt,isig) = driv2rhosigma(isig,ipt)*rho1(ipt,1) + driv2rhosigma(3+isig,ipt)*rho1(ipt,2)
               do jsig = 1, 3
                  drivsigma1(ipt,isig) = drivsigma1(ipt,isig) + driv2sigma2(sigmaIdx(isig,jsig),ipt)*sigma1(jsig)
               end do
            end do

            v_xc1(ipt,1) = v_xc1(ipt,1) - 2.0*drivsigma1(ipt,1)*gradRho%laplace(ipt,1) - drivsigma1(ipt,2)*gradRho%laplace(ipt,2)
            v_xc1(ipt,1) = v_xc1(ipt,1) - 2.0*drivsigma(1,ipt)*gradRho1%laplace(ipt,1) - drivsigma(2,ipt)*gradRho1%laplace(ipt,2)
            v_xc1(ipt,2) = v_xc1(ipt,2) - 2.0*drivsigma1(ipt,3)*gradRho%laplace(ipt,2) - drivsigma1(ipt,2)*gradRho%laplace(ipt,1)
            v_xc1(ipt,2) = v_xc1(ipt,2) - 2.0*drivsigma(3,ipt)*gradRho1%laplace(ipt,2) - drivsigma(2,ipt)*gradRho1%laplace(ipt,1)
         end do
      end if
   end subroutine dfpt_gga_assemble

   subroutine dfpt_gga_grdotgr(gradRho, gradRho1, gradDrivsigma, gradDrivsigma1, v_xc1)
      !! Adds the two grad.grad terms of -div(H) to the potential response of one
      !! channel. gradDrivsigma and gradDrivsigma1 hold the n_sigma fields in the
      !! slot that a gradient usually reserves for the spin.
      use m_types

      implicit none

      type(t_gradients), intent(in)    :: gradRho, gradRho1
      type(t_gradients), intent(in)    :: gradDrivsigma, gradDrivsigma1
      real,              intent(inout) :: v_xc1(:,:)

      integer :: ipt

      if (size(v_xc1,2) == 1) then
         do ipt = 1, size(v_xc1,1)
            v_xc1(ipt,1) = v_xc1(ipt,1) - 2.0*dot_product(gradDrivsigma1%gr(:,ipt,1),gradRho%gr(:,ipt,1))
            v_xc1(ipt,1) = v_xc1(ipt,1) - 2.0*dot_product(gradDrivsigma%gr(:,ipt,1),gradRho1%gr(:,ipt,1))
         end do
      else
         do ipt = 1, size(v_xc1,1)
            v_xc1(ipt,1) = v_xc1(ipt,1) - 2.0*dot_product(gradDrivsigma1%gr(:,ipt,1),gradRho%gr(:,ipt,1))
            v_xc1(ipt,1) = v_xc1(ipt,1) - dot_product(gradDrivsigma1%gr(:,ipt,2),gradRho%gr(:,ipt,2))
            v_xc1(ipt,1) = v_xc1(ipt,1) - 2.0*dot_product(gradDrivsigma%gr(:,ipt,1),gradRho1%gr(:,ipt,1))
            v_xc1(ipt,1) = v_xc1(ipt,1) - dot_product(gradDrivsigma%gr(:,ipt,2),gradRho1%gr(:,ipt,2))

            v_xc1(ipt,2) = v_xc1(ipt,2) - 2.0*dot_product(gradDrivsigma1%gr(:,ipt,3),gradRho%gr(:,ipt,2))
            v_xc1(ipt,2) = v_xc1(ipt,2) - dot_product(gradDrivsigma1%gr(:,ipt,2),gradRho%gr(:,ipt,1))
            v_xc1(ipt,2) = v_xc1(ipt,2) - 2.0*dot_product(gradDrivsigma%gr(:,ipt,3),gradRho1%gr(:,ipt,2))
            v_xc1(ipt,2) = v_xc1(ipt,2) - dot_product(gradDrivsigma%gr(:,ipt,2),gradRho1%gr(:,ipt,1))
         end do
      end if
   end subroutine dfpt_gga_grdotgr

end module m_dfpt_gga_kernel
