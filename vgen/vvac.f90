module m_vvac
   ! ****************************************************************
   ! calculates the g(2-dim)=0 part of the vacuum coulomb potential *
   ! for general symmetry.          c.l.fu, r.podloucky             *
   ! ****************************************************************
contains
   subroutine vvac( vacuum, stars, cell,  input, field, psq, rht, vnew, rhobar, sig1dh, vz1dh ,vslope)
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
      real,           intent(in)  :: rht(vacuum%nmzd,2)
      
      complex,        intent(out) :: vnew(vacuum%nmzd,2)
      complex,        intent(out) :: rhobar
      real,           intent(out) :: sig1dh, vz1dh
      complex,        intent(out) :: vslope

      real                        :: sumq
      real                        :: bj0, bj1, qzh, sigmaa(2)
      integer                     :: ig3, imz, ivac, ncsh
      real                        :: f(vacuum%nmzd), sig(vacuum%nmzd), vtemp(vacuum%nmzd)

      vnew(:,1:vacuum%nvac) = 0.0 ! initialize potential

      ! obtain mesh point (ncsh) of charge sheet for external electric field
      ncsh = field%efield%zsigma / vacuum%delz + 1.01
      sigmaa(1) = ( field%efield%sigma + field%efield%sig_b(1) ) / cell%area
      sigmaa(2) = ( field%efield%sigma + field%efield%sig_b(2) ) / cell%area

      ! g=0 vacuum potential due to neutral charge density
      ! inside slab and zero charge density outside

      rhobar = - psq(1)
      sumq = 0.0

      do ig3 = 2, stars%ng3
         if (stars%ig2(ig3) == 1) then           ! select g_|| = 0
            qzh = stars%kv3(3,ig3) * cell%bmat(3,3) * cell%z1
            bj0 = sin(qzh) / qzh
            rhobar = rhobar - psq(ig3) * bj0 * stars%nstr(ig3)
            if ( vacuum%nvac==2 ) then
               bj1 = ( sin(qzh) - qzh * cos(qzh) ) / ( qzh * qzh )
               sumq = sumq - 2. * fpi_const * ImagUnit * bj1 * psq(ig3) * cell%z1 *cell%z1
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
         call qsf( vacuum%delz, rht(1,ivac), sig, ncsh, 1 )
         sig(1:ncsh) = sig(ncsh) - sig(1:ncsh)
         call qsf( vacuum%delz, sig, vtemp, ncsh, 1 )
         do imz = 1, ncsh
            vnew(imz,ivac) = - fpi_const * ( vtemp(ncsh) - vtemp(imz) ) + field%efield%sig_b(1)
         end do
         sig1dh = sig(1)
         vz1dh = vnew(1,ivac)   ! potential on vacuum boundary
         if ( vacuum%nvac == 1 ) return

         ivac = 2     ! lower vacuum
         call qsf( vacuum%delz, rht(1,ivac), sig, ncsh, 1 )
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
         call qsf( vacuum%delz, rht(1,ivac), sig, vacuum%nmz, 1 )
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
         call qsf( vacuum%delz, rht(1,ivac), sig, vacuum%nmz, 1 )
         f(1:vacuum%nmz) = sig(1:vacuum%nmz) - rhobar * vacuum%dvac + sig1dh
         call qsf( vacuum%delz, f, vtemp, vacuum%nmz, 1 )

         ! external electric field contribution
         do imz = 1, ncsh
            vnew(imz,ivac) = - fpi_const * ( vtemp(imz) + sig1dh * vacuum%dvac - rhobar * vacuum%dvac * vacuum%dvac / 2. ) + vz1dh + vnew(imz,ivac)
         end do
         do imz = ncsh + 1, vacuum%nmz
            vnew(imz,ivac) = - fpi_const * ( vtemp(imz) + sig1dh * vacuum%dvac - rhobar * vacuum%dvac * vacuum%dvac / 2. ) + vz1dh + vnew(imz,ivac) &
                           + fpi_const * ( imz - ncsh ) * vacuum%delz * sigmaa(2)
         end do
      end if ! Dirichlet/Neumann
   end subroutine vvac
end module m_vvac
