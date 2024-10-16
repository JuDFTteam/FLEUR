!-------------------------------------------------------------------------------
! Copyright (c) 2022 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions
! of the MIT license as expressed in the LICENSE file in more detail.
!-------------------------------------------------------------------------------

module m_vvac
   ! Legacy comment:
   ! ****************************************************************
   ! calculates the g(2-dim)=0 part of the vacuum coulomb potential *
   ! for general symmetry.          c.l.fu, r.podloucky             *
   ! ****************************************************************
contains
   subroutine vvac(vacuum, stars, cell, input, field, psq, rht, vnew, rhobar, sig1dh, vz1dh, vslope, l_bind, vmz1dh, sigma_disc, sigma_disc2)
      !! Calculates the \(\boldsymbol{G}_{||}=0\) part of the vacuum Coulomb potential.
      !! There are two possible cases for Dirichlet and von Neumann boundary conditions, respectively.
      !! von Neumann case:
      !! upper vacuum \( z\ge D/2)\):
      !! $$V^{0}(z)=-4\pi[\int_{z}^{\infty}\sigma_{+}^{0}(z')dz'+(z-z_{\sigma})\sigma_{+}]=: V_{+}^{0}(z)$$
      !! with
      !! $$\sigma_{+}^{0}(z):=\int_{z}^{\infty}n_{V}^{0}(z')dz'$$
      use m_constants
      use m_qsf
      use m_types
    
      implicit none

      type(t_vacuum), intent(in)  :: vacuum
      type(t_stars),  intent(in)  :: stars
      type(t_cell),   intent(in)  :: cell
      type(t_input),  intent(in)  :: input
      type(t_field),  intent(in)  :: field

      complex,        intent(in)  :: psq(stars%ng3)
      complex,        intent(in)  :: rht(vacuum%nmzd,2)
      
      complex,        intent(out) :: vnew(vacuum%nmzd,2)
      complex,        intent(out) :: rhobar
      complex,        intent(out) :: sig1dh, vz1dh
      complex,        intent(out) :: vslope
      logical,        intent(in)  :: l_bind
      complex,        intent(out) :: vmz1dh
      complex,        intent(in)  :: sigma_disc(2)

      complex, optional, intent(in) :: sigma_disc2(2)

      complex                     :: sumq, newdp, newdm, newdp2, newdm2
      real                        :: bj0, bj1, qzh, sigmaa(2)
      integer                     :: ig3, imz, ivac, ncsh
      real                        :: f(vacuum%nmzd), sig(vacuum%nmzd), vtemp(vacuum%nmzd)

      vmz1dh = cmplx(0.0,0.0)

      newdp =  cmplx(0.0,0.0)
      newdm =  cmplx(0.0,0.0)
      newdp2 =  cmplx(0.0,0.0)
      newdm2 =  cmplx(0.0,0.0)

      vnew(:,1:vacuum%nvac) = 0.0 ! initialize potential

      ! obtain mesh point (ncsh) of charge sheet for external electric field
      ncsh = field%efield%zsigma / vacuum%delz + 1.01
      sigmaa(1) = ( field%efield%sigma + field%efield%sig_b(1) ) / cell%area
      sigmaa(2) = ( field%efield%sigma + field%efield%sig_b(2) ) / cell%area

      ! g=0 vacuum potential due to neutral charge density
      ! inside slab and zero charge density outside

      rhobar = - psq(1)
      sumq = 0.0

      IF (l_bind) newdp = -psq(1) * cell%z1
      IF (l_bind) newdm = -psq(1) * cell%z1
      IF (l_bind) newdp2 = -psq(1) * cell%z1 * cell%z1
      IF (l_bind) newdm2 =  psq(1) * cell%z1 * cell%z1

      do ig3 = 2, stars%ng3
         if (stars%ig2(ig3) == 1) then           ! select g_|| = 0
            qzh = stars%kv3(3,ig3) * cell%bmat(3,3) * cell%z1
            bj0 = sin(qzh) / qzh
            rhobar = rhobar - psq(ig3) * bj0 * stars%nstr(ig3)
            if ( vacuum%nvac==2 ) then
               bj1 = ( sin(qzh) - qzh * cos(qzh) ) / ( qzh * qzh )
               sumq = sumq - 2. * fpi_const * ImagUnit * bj1 * psq(ig3) * cell%z1 *cell%z1

               ! New discontinuity correction
               IF (l_bind) newdp = newdp + ImagUnit * psq(ig3) * cmplx(cos(qzh), sin(qzh)) * cell%z1 / qzh
               IF (l_bind) newdm = newdm - ImagUnit * psq(ig3) * cmplx(cos(qzh),-sin(qzh)) * cell%z1 / qzh
               IF (l_bind) newdp2 = newdp2 + psq(ig3) * cmplx(cos(qzh), sin(qzh)) * cell%z1**2 / qzh**2
               IF (l_bind) newdm2 = newdm2 - psq(ig3) * cmplx(cos(qzh),-sin(qzh)) * cell%z1**2 / qzh**2

            end if
         end if
      end do ! --> qzh, bj0, bj1, psq finished; rhobar, sumq passed on and unchanged
       
      ! lower (nvac=2) vacuum
      if ( vacuum%nvac == 2 ) vnew(1:vacuum%nmz,2) =  sumq

      ! g=0 vacuum potential due to
      ! negative of rhobar + vacuum (g=0) charge ----> v2(z)

      ivac = 1 ! upper vacuum

      if ( field%efield%dirichlet ) then ! Dirichlet
         vnew(ncsh+1:vacuum%nmz,ivac) = field%efield%sig_b(1)
         call qsf( vacuum%delz, REAL(rht(:,ivac)), sig, ncsh, 1 )
         sig(1:ncsh) = sig(ncsh) - sig(1:ncsh)
         call qsf( vacuum%delz, sig, vtemp, ncsh, 1 )
         do imz = 1, ncsh
            vnew(imz,ivac) = - fpi_const * ( vtemp(ncsh) - vtemp(imz) ) + field%efield%sig_b(1)
         end do
         sig1dh = sig(1)
         vz1dh = vnew(1,ivac)   ! potential on vacuum boundary
         if ( vacuum%nvac == 1 ) return

         ivac = 2     ! lower vacuum
         call qsf( vacuum%delz, REAL(rht(:,ivac)), sig, ncsh, 1 )
         f(1:ncsh) = sig(1:ncsh) - rhobar*vacuum%dvac + sig1dh
         call qsf( vacuum%delz, f, vtemp, ncsh, 1 )
         do imz = 1,ncsh
            vnew(imz,ivac) = - fpi_const * ( vtemp(imz) + sig1dh * vacuum%dvac - rhobar * vacuum%dvac * vacuum%dvac / 2. ) + vnew(imz,ivac) + vz1dh
         end do

         ! force matching on the other side
         vslope = ( field%efield%sig_b(2) - vnew(ncsh,1) ) / ( 2 * vacuum%delz * ( ncsh + 1 ) + vacuum%dvac )
         ivac = 1
         do imz = 1, ncsh
            vnew(imz,ivac) = vnew(imz,ivac) + vslope * vacuum%delz * ( ncsh - imz + 1 )
         end do
         ivac = 2
         do imz = 1, ncsh
            vnew(imz,ivac) = vnew(imz,ivac) + vslope * ( vacuum%dvac + vacuum%delz * imz + vacuum%delz * ncsh )
         end do
         vnew(ncsh+1:vacuum%nmz,ivac) = field%efield%sig_b(2)
      else ! Neumann
         call qsf( vacuum%delz, REAL(rht(:,ivac)), sig, vacuum%nmz, 1 )
         sig1dh = sig(vacuum%nmz) - sigmaa(1)  ! need to include contribution from electric field
         sig(1:vacuum%nmz) = sig(vacuum%nmz) - sig(1:vacuum%nmz)
         call qsf( vacuum%delz, sig, vtemp, vacuum%nmz, 1 )
         ! external electric field contribution
         do imz = 1, ncsh
            vnew(imz,ivac) = - fpi_const * ( vtemp(vacuum%nmz) - vtemp(imz) ) + vnew(imz,ivac) - fpi_const * ( imz - ncsh ) * vacuum%delz * sigmaa(1)
         end do
         do imz = ncsh + 1, vacuum%nmz
            vnew(imz,ivac) = - fpi_const * ( vtemp(vacuum%nmz) - vtemp(imz) ) + vnew(imz,ivac)
         end do
         vz1dh = vnew(1,ivac)   ! potential on vacuum boundary
         if ( vacuum%nvac == 1 ) return

         ivac = 2 ! lower vacuum
         call qsf( vacuum%delz, REAL(rht(:,ivac)), sig, vacuum%nmz, 1 )
         f(1:vacuum%nmz) = sig(1:vacuum%nmz) - rhobar * vacuum%dvac + sig1dh
         call qsf( vacuum%delz, f, vtemp, vacuum%nmz, 1 )
         ! external electric field contribution
         do imz = 1, ncsh
            vnew(imz,ivac) = - fpi_const * ( vtemp(imz) + sig1dh * vacuum%dvac - rhobar * vacuum%dvac * vacuum%dvac / 2. ) + vz1dh + vnew(imz,ivac) &
                             - fpi_const * (sigma_disc(1) * ( vacuum%dvac / 2. - (imz-1) * vacuum%delz ) - sigma_disc(2) * ( vacuum%dvac / 2. + (imz-1) * vacuum%delz )) &! Discontinuity correction
                             - fpi_const * (-newdp2-newdm2+newdp * ( vacuum%dvac / 2. - (imz-1) * vacuum%delz ) - newdm * ( vacuum%dvac / 2. + (imz-1) * vacuum%delz ))! New discontinuity correction
            !if (present(sigma_disc2)) vnew(imz,ivac) = vnew(imz,ivac) - fpi_const * (sigma_disc2(1)+sigma_disc2(2))
         end do
         do imz = ncsh + 1, vacuum%nmz
            vnew(imz,ivac) = - fpi_const * ( vtemp(imz) + sig1dh * vacuum%dvac - rhobar * vacuum%dvac * vacuum%dvac / 2. ) + vz1dh + vnew(imz,ivac) &
                           + fpi_const * ( imz - ncsh ) * vacuum%delz * sigmaa(2) & ! Discontinuity correction
                           - fpi_const * (sigma_disc(1) * ( vacuum%dvac / 2. - (imz-1) * vacuum%delz ) - sigma_disc(2) * ( vacuum%dvac / 2. + (imz-1) * vacuum%delz )) &! Discontinuity correction
                           - fpi_const * (-newdp2-newdm2+newdp * ( vacuum%dvac / 2. - (imz-1) * vacuum%delz ) - newdm * ( vacuum%dvac / 2. + (imz-1) * vacuum%delz )) ! New discontinuity correction
            !if (present(sigma_disc2)) vnew(imz,ivac) = vnew(imz,ivac) - fpi_const * (sigma_disc2(1)+sigma_disc2(2))
         end do
         !if (l_bind) then
         !   ! Fix the potential to 0 at -infinity and save the resulting value at the vacuum border -D/2
         !   vnew(:,2) = vnew(:,2) - vnew(vacuum%nmz,2)
         !   vmz1dh = vnew(1,2)
         !end if
         ! Discontinuity correction
         !if (l_bind) then
            ! Fix the potential to 0 at -infinity and save the resulting value at the vacuum border -D/2
            !vnew(:,2) = cmplx(0.0,0.0) 
            !vmz1dh = vnew(1,2)
         !end if
      end if ! Dirichlet/Neumann
   end subroutine vvac
end module m_vvac
