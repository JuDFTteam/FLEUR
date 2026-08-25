!--------------------------------------------------------------------------------
! Copyright (c) 2026 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------
module m_dfpt_gga_kernel
   !! Region independent assembly of the GGA part of the xc potential response
   !! Dv_xc = Dvrho - div(H) from the libxc kernels vsigma, v2rho2, v2rhosigma
   !! and v2sigma2. The divergence is expanded by the product rule: the terms
   !! containing a Laplacian are done in dfpt_gga_local, the remaining two
   !! grad.grad terms in dfpt_gga_grterms once the caller has supplied the
   !! gradients of vsigma and of Dvsigma.
   !!
   !! Everything is linear in the response, so the real and the imaginary
   !! channel of Drho are two independent calls sharing the same kernels.

   implicit none
   private

   ! Symmetric packing of v2sigma2, libxc order [11,12,13,22,23,33]
   integer, parameter :: ss_idx(3,3) = reshape([1,2,3, 2,4,5, 3,5,6],[3,3])

   public :: dfpt_gga_local, dfpt_gga_grterms

contains

   subroutine dfpt_gga_local(vsigma, v2rho2, v2rhosigma, v2sigma2, grad0, grad1, rho1, v_xc1, dvsigma1)
      !! Local part of the GGA potential response for one channel. Returns v_xc1
      !! including the Laplacian terms of -div(H), and Dvsigma in dvsigma1, whose
      !! gradient the caller has to feed back through dfpt_gga_grterms.
      use m_types

      implicit none

      real,              intent(in)  :: vsigma(:,:)     ! (n_sigma, points)
      real,              intent(in)  :: v2rho2(:,:)     ! (2*jspins-1, points)
      real,              intent(in)  :: v2rhosigma(:,:) ! (jspins*n_sigma, points)
      real,              intent(in)  :: v2sigma2(:,:)   ! (1 or 6, points)
      type(t_gradients), intent(in)  :: grad0, grad1
      real,              intent(in)  :: rho1(:,:)       ! (points, jspins)
      real,              intent(out) :: v_xc1(:,:)      ! (points, jspins)
      real,              intent(out) :: dvsigma1(:,:)   ! (points, n_sigma)

      integer :: i, mu, nu
      real    :: ds(3)

      v_xc1 = 0.0
      dvsigma1 = 0.0

      if (size(rho1,2) == 1) then
         do i = 1, size(rho1,1)
            ds(1) = 2.0*dot_product(grad0%gr(:,i,1),grad1%gr(:,i,1))

            v_xc1(i,1) = v2rho2(1,i)*rho1(i,1) + v2rhosigma(1,i)*ds(1)
            dvsigma1(i,1) = v2rhosigma(1,i)*rho1(i,1) + v2sigma2(1,i)*ds(1)

            v_xc1(i,1) = v_xc1(i,1) - 2.0*dvsigma1(i,1)*grad0%laplace(i,1)
            v_xc1(i,1) = v_xc1(i,1) - 2.0*vsigma(1,i)*grad1%laplace(i,1)
         end do
      else
         do i = 1, size(rho1,1)
            ds(1) = 2.0*dot_product(grad0%gr(:,i,1),grad1%gr(:,i,1))
            ds(2) = dot_product(grad0%gr(:,i,1),grad1%gr(:,i,2)) + dot_product(grad0%gr(:,i,2),grad1%gr(:,i,1))
            ds(3) = 2.0*dot_product(grad0%gr(:,i,2),grad1%gr(:,i,2))

            v_xc1(i,1) = v2rho2(1,i)*rho1(i,1) + v2rho2(2,i)*rho1(i,2)
            v_xc1(i,2) = v2rho2(2,i)*rho1(i,1) + v2rho2(3,i)*rho1(i,2)

            do mu = 1, 3
               v_xc1(i,1) = v_xc1(i,1) + v2rhosigma(mu,i)*ds(mu)
               v_xc1(i,2) = v_xc1(i,2) + v2rhosigma(3+mu,i)*ds(mu)
               dvsigma1(i,mu) = v2rhosigma(mu,i)*rho1(i,1) + v2rhosigma(3+mu,i)*rho1(i,2)
               do nu = 1, 3
                  dvsigma1(i,mu) = dvsigma1(i,mu) + v2sigma2(ss_idx(mu,nu),i)*ds(nu)
               end do
            end do

            v_xc1(i,1) = v_xc1(i,1) - 2.0*dvsigma1(i,1)*grad0%laplace(i,1) - dvsigma1(i,2)*grad0%laplace(i,2)
            v_xc1(i,1) = v_xc1(i,1) - 2.0*vsigma(1,i)*grad1%laplace(i,1) - vsigma(2,i)*grad1%laplace(i,2)
            v_xc1(i,2) = v_xc1(i,2) - 2.0*dvsigma1(i,3)*grad0%laplace(i,2) - dvsigma1(i,2)*grad0%laplace(i,1)
            v_xc1(i,2) = v_xc1(i,2) - 2.0*vsigma(3,i)*grad1%laplace(i,2) - vsigma(2,i)*grad1%laplace(i,1)
         end do
      end if
   end subroutine dfpt_gga_local

   subroutine dfpt_gga_grterms(grad0, grad1, grad_vsigma, grad_dvsigma, v_xc1)
      !! Adds the two grad.grad terms of -div(H) to the potential response of one
      !! channel. grad_vsigma and grad_dvsigma hold the n_sigma fields in the slot
      !! that grad usually reserves for the spin.
      use m_types

      implicit none

      type(t_gradients), intent(in)    :: grad0, grad1
      type(t_gradients), intent(in)    :: grad_vsigma, grad_dvsigma
      real,              intent(inout) :: v_xc1(:,:)

      integer :: i

      if (size(v_xc1,2) == 1) then
         do i = 1, size(v_xc1,1)
            v_xc1(i,1) = v_xc1(i,1) - 2.0*dot_product(grad_dvsigma%gr(:,i,1),grad0%gr(:,i,1))
            v_xc1(i,1) = v_xc1(i,1) - 2.0*dot_product(grad_vsigma%gr(:,i,1),grad1%gr(:,i,1))
         end do
      else
         do i = 1, size(v_xc1,1)
            v_xc1(i,1) = v_xc1(i,1) - 2.0*dot_product(grad_dvsigma%gr(:,i,1),grad0%gr(:,i,1))
            v_xc1(i,1) = v_xc1(i,1) - dot_product(grad_dvsigma%gr(:,i,2),grad0%gr(:,i,2))
            v_xc1(i,1) = v_xc1(i,1) - 2.0*dot_product(grad_vsigma%gr(:,i,1),grad1%gr(:,i,1))
            v_xc1(i,1) = v_xc1(i,1) - dot_product(grad_vsigma%gr(:,i,2),grad1%gr(:,i,2))

            v_xc1(i,2) = v_xc1(i,2) - 2.0*dot_product(grad_dvsigma%gr(:,i,3),grad0%gr(:,i,2))
            v_xc1(i,2) = v_xc1(i,2) - dot_product(grad_dvsigma%gr(:,i,2),grad0%gr(:,i,1))
            v_xc1(i,2) = v_xc1(i,2) - 2.0*dot_product(grad_vsigma%gr(:,i,3),grad1%gr(:,i,2))
            v_xc1(i,2) = v_xc1(i,2) - dot_product(grad_vsigma%gr(:,i,2),grad1%gr(:,i,1))
         end do
      end if
   end subroutine dfpt_gga_grterms

end module m_dfpt_gga_kernel
