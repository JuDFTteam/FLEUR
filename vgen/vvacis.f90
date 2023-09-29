module m_vvacis
  !     **********************************************************
  !     g.ne.0 coefficients of vacuum coulomb potential          *
  !     due to the interstitial charge density inside slab       *
  !                                   c.l.fu, r.podloucky        *
  !     **********************************************************
  !     modified for thick films to avoid underflows gb`06
  !---------------------------------------------------------------
contains
   subroutine vvacis( stars, vacuum, cell, psq, input, field, vxy, l_dfptvgen )
      use m_constants
      use m_types

      implicit none

      type(t_input),  intent(in)  :: input
      type(t_field),  intent(in)  :: field
      type(t_vacuum), intent(in)  :: vacuum
      type(t_stars),  intent(in)  :: stars
      type(t_cell),   intent(in)  :: cell

      complex,        intent(in)  :: psq(stars%ng3)
      complex,        intent(out) :: vxy(vacuum%nmzxyd,stars%ng2,2)
      logical,        intent(in)  :: l_dfptvgen

      complex                     :: arg, c_ph, sumr(2)
      real                        :: g, qz, sign,  vcons, z, e_m, arg_r, arg_i
      integer                     :: ig3n, imz, ivac, k1, k2, kz, nrec2, start_star
  
      start_star = 2
      if (l_dfptvgen) then
         if (norm2(stars%center)<1e-8) start_star = 1
      end if

      vxy(:,start_star:,:) = cmplx( 0., 0. )

      do  nrec2 = start_star, stars%ng2
         k1 = stars%kv2(1,nrec2)
         k2 = stars%kv2(2,nrec2)
         g = stars%sk2(nrec2)
         if ( field%efield%dirichlet ) then
            vcons = 2.0 * tpi_const / ( g * sinh( g * 2.0 * ( field%efield%zsigma + cell%z1 ) ) )
            arg_r = g * ( cell%z1 + field%efield%zsigma + cell%z1 )
         else ! neumann
            vcons = tpi_const / g
            arg_r = exp_safe( - 2 * cell%z1 * g )
         end if
         do ivac = 1, vacuum%nvac
            sumr(ivac) = ( 0.0, 0.0 )
            sign = 3. - 2. * ivac
            do kz = -stars%mx3, stars%mx3
               ig3n = stars%ig(k1,k2,kz)
               c_ph = stars%rgphs(k1,k2,kz)
               ! use only stars within the g_max sphere
               if ( ig3n /= 0 ) then
                  qz = kz * cell%bmat(3,3)
                  ! sum over gz-stars
                  if ( field%efield%dirichlet ) then
                     ! prefactor
                     arg = exp( - ImagUnit *  qz * cell%z1 ) / ( 2 * ( g ** 2 + qz ** 2 ) ) * psq(ig3n)
                     if ( ivac == 1 ) then
                        sumr(ivac) = sumr(ivac) + c_ph * exp( - arg_r ) * arg * ( &                           ! c_ph not tested in this case
                            ( - exp( 2 * g * ( field%efield%zsigma + cell%z1 ) ) + exp( 2 * ( ImagUnit  * qz * cell%z1 + arg_r ) ) ) * ( g - ImagUnit *  qz ) &
                          + ( - exp( 2 * g * cell%z1 ) + exp( 2 * ImagUnit  * qz * cell%z1 ) )                                       * ( g + ImagUnit *  qz ) )
                     else
                        sumr(ivac) = sumr(ivac) + c_ph * arg * ( &
                             exp(   arg_r ) * ( g + ImagUnit *  qz ) &
                           + exp( - arg_r ) * ( g - ImagUnit *  qz ) &
                       + 2 * exp( 2 * ImagUnit * qz * cell%z1 ) &
                        * ( - g * cosh( g * ( - field%efield%zsigma ) ) &
                        + ImagUnit *  qz * sinh( - g * field%efield%zsigma ) ) )
                     end if
                  else
                     arg = g + sign * ImagUnit *  qz
                     arg_i = sign  * qz * cell%z1
                     sumr(ivac) = sumr(ivac) + c_ph * psq(ig3n) * ( cos( arg_i ) * ( 1 - arg_r ) + ImagUnit * sin( arg_i ) * ( 1 + arg_r ) ) / arg
                  end if
               end if
            end do
            z = 0 ! moved cell%z1 into above equations
            do imz = 1, vacuum%nmzxy
               if ( field%efield%dirichlet ) then
                  e_m = sinh( g * ( field%efield%zsigma - z ) )
               else ! neumann
                  e_m = exp_safe( - g * z )
               end if
               vxy(imz,nrec2,ivac) = vxy(imz,nrec2,ivac) + vcons * sumr(ivac) * e_m
               z = z + vacuum%delz
            end do
         end do
      end do
   end subroutine vvacis

   pure real function exp_safe(x)
      ! replace exp by a function that does not under/overflow

      implicit none

      real, intent(in) :: x
      real, parameter  :: maxexp = log( 2.0 ) * maxexponent( 2.0 )
      real, parameter  :: minexp = log( 2.0 ) * minexponent( 2.0 )

      if ( abs(x) > minexp .and. abs(x) < maxexp ) then
         exp_safe = exp( x )
      else
         if ( x > 0 ) then
            if ( x > minexp ) then
               exp_safe = exp( maxexp )
            else
               exp_safe = exp( minexp )
            end if
         else
            if ( -x > minexp ) then
               exp_safe = exp( -maxexp )
            else
               exp_safe = exp( -minexp )
            end if
         end if
      end if
  end function exp_safe
end module m_vvacis
