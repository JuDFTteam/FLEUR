module m_vvacis
  !     **********************************************************
  !     g.ne.0 coefficients of vacuum coulomb potential          *
  !     due to the interstitial charge density inside slab       *
  !                                   c.l.fu, r.podloucky        *
  !     **********************************************************
  !     modified for thick films to avoid underflows gb`06
  !---------------------------------------------------------------
  contains

  subroutine vvacis( stars, vacuum, sym, cell, psq, input, field, vxy )

    use m_constants
    use m_types
    implicit none

    type(t_input),  intent(in)  :: input
    type(t_field),  intent(in)  :: field
    type(t_vacuum), intent(in)  :: vacuum
    type(t_sym),    intent(in)  :: sym
    type(t_stars),  intent(in)  :: stars
    type(t_cell),   intent(in)  :: cell
    complex,        intent(in)  :: psq(stars%ng3)
    complex,        intent(out) :: vxy(vacuum%nmzxyd,stars%ng2-1,2)

    complex                     :: arg, c_ph, sumr(2)
    real                        :: dh, g, qz, sign, signz, vcons, z, e_m, arg_r, arg_i
    integer                     :: i2d, ig3n, imz, imzxy, ivac, k1, k2, kz, m0, nrec2, nrz, nz
    intrinsic exp

    vxy(:,:,:) = cmplx( 0., 0. )
    dh = cell%z1
    m0 = -stars%mx3
    if ( sym%zrfs ) m0 = 0
    do  nrec2 = 2, stars%ng2
      k1 = stars%kv2(1,nrec2)
      k2 = stars%kv2(2,nrec2)
      g = stars%sk2(nrec2)
      if ( field%efield%dirichlet ) then
        vcons = 2.0 * tpi_const / ( g * sinh( g * 2.0 * ( field%efield%zsigma + dh ) ) )
        arg_r = g * ( dh + field%efield%zsigma + dh )
      else ! neumann
        vcons = tpi_const / g
        arg_r = exp_save( - 2 * dh * g )
      end if
      do ivac = 1, vacuum%nvac
        sumr(ivac) = ( 0.0, 0.0 )
        sign = 3. - 2. * ivac
        do kz = m0, stars%mx3
          ig3n = stars%ig(k1,k2,kz)
          ! use only stars within the g_max sphere
          if ( ig3n /= 0 ) then
            c_ph = stars%rgphs(k1,k2,kz)
            nz = 1
            if (sym%zrfs) nz = stars%nstr(ig3n) / stars%nstr2(nrec2)
            qz = kz * cell%bmat(3,3)
            ! sum over gz-stars
            do  nrz = 1, nz
              signz = 3. - 2. * nrz
              if ( field%efield%dirichlet ) then
                ! prefactor
                arg = exp( - ImagUnit * signz * qz * dh ) / ( 2 * ( g ** 2 + qz ** 2 ) ) * psq(ig3n)
                if ( ivac == 1 ) then
                  sumr(ivac) = sumr(ivac) + c_ph * exp( - arg_r ) * arg * ( &                           ! c_ph not tested in this case
                      ( - exp( 2 * g * ( field%efield%zsigma + dh ) ) + exp( 2 * ( ImagUnit * signz * qz * dh + arg_r ) ) ) * ( g - ImagUnit * signz * qz ) &
                    + ( - exp( 2 * g * dh ) + exp( 2 * ImagUnit * signz * qz * dh ) )                                       * ( g + ImagUnit * signz * qz ) )
                else
                  sumr(ivac) = sumr(ivac) + c_ph * arg * ( &
                    exp(   arg_r ) * ( g + ImagUnit * signz * qz ) &
                    + exp( - arg_r ) * ( g - ImagUnit * signz * qz ) &
                    + 2 * exp( 2 * ImagUnit * signz * qz * dh ) &
                    * ( - g * cosh( g * ( - field%efield%zsigma ) ) &
                      + ImagUnit * signz * qz * sinh( - g * field%efield%zsigma ) ) )
                end if
              else
                arg = g + sign * ImagUnit * signz * qz
                arg_i = sign * signz * qz * dh
                sumr(ivac) = sumr(ivac) + c_ph * psq(ig3n) * ( cos( arg_i ) * ( 1 - arg_r ) + ImagUnit * sin( arg_i ) * ( 1 + arg_r ) ) / arg
              end if
            enddo
          endif
        enddo
        z = 0 ! moved cell%z1 into above equations
        do imz = 1, vacuum%nmzxy
          if ( field%efield%dirichlet ) then
            e_m = sinh( g * ( field%efield%zsigma - z ) )
          else ! neumann
            e_m = exp_save( - g * z )
          end if
          vxy(imz,nrec2-1,ivac) = vxy(imz,nrec2-1,ivac) + vcons * sumr(ivac) * e_m
          z = z + vacuum%delz
        enddo
      enddo
    enddo

  end subroutine vvacis


  pure real function exp_save(x)
    ! replace exp by a function that does not under/overflow

    implicit none
    real, intent(in) :: x
    real, parameter  :: maxexp = log( 2.0 ) * maxexponent( 2.0 )
    real, parameter  :: minexp = log( 2.0 ) * minexponent( 2.0 )

    if ( abs(x) > minexp .and. abs(x) < maxexp ) then
      exp_save = exp( x )
    else
      if ( x > 0 ) then
        if ( x > minexp ) then
          exp_save = exp( maxexp )
        else
          exp_save = exp( minexp )
        endif
      else
        if ( -x > minexp ) then
          exp_save = exp( -maxexp )
        else
          exp_save = exp( -minexp )
        endif
      endif
    endif
  end function exp_save

end module m_vvacis
