module m_vvacis
  !     **********************************************************
  !     g.ne.0 coefficients of vacuum coulomb potential          *
  !     due to the interstitial charge density inside slab       *
  !                                   c.l.fu, r.podloucky        *
  !     **********************************************************
  !     modified for thick films to avoid underflows gb`06
  !---------------------------------------------------------------
  contains

  subroutine vvacis( stars, vacuum, cell, psq, input, field, vxy )

    use m_constants
    use m_types
    use m_expsave
    implicit none

    type(t_input),  intent(in)  :: input
    type(t_field),  intent(in)  :: field
    type(t_vacuum), intent(in)  :: vacuum
    type(t_stars),  intent(in)  :: stars
    type(t_cell),   intent(in)  :: cell
    complex,        intent(in)  :: psq(stars%ng3)
    complex,        intent(out) :: vxy(vacuum%nmzxyd,stars%ng2-1,2)

    complex                     :: arg, c_ph, sumr(2)
    real                        :: g, qz, sign,  vcons, z, e_m, arg_r, arg_i
    integer                     :: i2d, ig3n, imz, imzxy, ivac, k1, k2, kz,  nrec2, nrz, nz
    intrinsic exp

    vxy(:,:,:) = cmplx( 0., 0. )
  
    do  nrec2 = 2, stars%ng2
      k1 = stars%kv2(1,nrec2)
      k2 = stars%kv2(2,nrec2)
      g = stars%sk2(nrec2)
      if ( field%efield%dirichlet ) then
        vcons = 2.0 * tpi_const / ( g * sinh( g * 2.0 * ( field%efield%zsigma + cell%z1 ) ) )
        arg_r = g * ( cell%z1 + field%efield%zsigma + cell%z1 )
      else ! neumann
        vcons = tpi_const / g
        arg_r = exp_save( - 2 * cell%z1 * g )
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


  

end module m_vvacis
