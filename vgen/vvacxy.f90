module m_vvacxy
  !     **********************************************************
  !     g/=0 coefficient of vacuum coulomb potential           *
  !     due to warped vacuum charged density                     *
  !                                 c.l.fu, r.podloucky          *
  !     **********************************************************
  !     modified for thick films to avoid underflows gb`06
  !---------------------------------------------------------------

  use m_judft

  contains

  subroutine vvacxy( stars, vacuum, cell, sym, input, field, rhtxy, vxy, alphm )

    use m_intgr, only: intgz1
    use m_constants
    use m_types
    use m_qsf
    implicit none

    type(t_input),  intent(in)     :: input
    type(t_field),  intent(in)     :: field
    type(t_vacuum), intent(in)     :: vacuum
    type(t_sym),    intent(in)     :: sym
    type(t_stars),  intent(in)     :: stars
    type(t_cell),   intent(in)     :: cell
    complex,        intent(in)     :: rhtxy(vacuum%nmzxy,stars%ng2-1,2)
    complex,        intent(inout)  :: vxy(vacuum%nmzxy,stars%ng2-1,2)
    complex,        intent(out)    :: alphm(stars%ng2,2)

    complex                        :: alph0, alph2, alph1, alphaz, betaz, test
    real                           :: g, vcons, z, e_m, e_p
    integer                        :: imz, ip, irec2, ivac, ncsh
    logical                        :: tail
    real                           :: fra(vacuum%nmzxy), frb(vacuum%nmzxy), fia(vacuum%nmzxy), fib(vacuum%nmzxy)
    real                           :: alpha(vacuum%nmzxy,2,2), beta(vacuum%nmzxy,2,2)
    real, allocatable              :: sig_top(:), sig_bot(:)
    intrinsic aimag, cmplx, conjg, exp, real

    if ( allocated( field%efield%rhoef ) .or. field%efield%dirichlet ) then
      ncsh = field%efield%zsigma / vacuum%delz + 1.01
      ! if nmzxy < ncsh, the inhomogenous field cannot be represented.
      ! if nmzxy > ncsh, the boundary condition is wrong - and the
      ! potential is very wavy - try it yourself, if you don't believe.
      if ( vacuum%nmzxy < ncsh ) then
        write (6,*) 'error vvacxy.f: vacuum%nmzxy =', vacuum%nmzxy, '< ncsh(zsigma) = ', ncsh
        call judft_error( "error: vacuum%nmzxy < ncsh", calledby="vvacxy" )
      else if ( vacuum%nmzxy > ncsh ) then
        write (6,*) 'warning vvacxy.f: vacuum%nmzxy =', vacuum%nmzxy, '> ncsh(zsigma) = ', ncsh
        write (0,*) 'warning vvacxy.f: vacuum%nmzxy =', vacuum%nmzxy, '> ncsh(zsigma) = ', ncsh
        call judft_warn( "nmzxyd > ncsh", calledby="vvacxy" )
      end if
    end if

    ! 2-dim star loop g/=0
    ip = vacuum%nmzxy + 1
    irec2_loop: do irec2 = 2, stars%ng2
      g = stars%sk2(irec2)
      vcons = tpi_const / g
      ! ********** dirichlet ************************************
      if ( field%efield%dirichlet ) then
        if ( allocated( field%efield%rhoef ) ) then
          vxy(ncsh:vacuum%nmzxy,irec2-1,1) = field%efield%rhoef(irec2-1, 1)
          if ( vacuum%nvac == 2 ) then
            vxy(ncsh:vacuum%nmzxy,irec2-1,2) = field%efield%rhoef(irec2-1, 2)
          end if
        else
          vxy(ncsh:vacuum%nmzxy,irec2-1,1:vacuum%nvac) = 0.0
        end if
        vcons = 2.0 * vcons / sinh( g * 2.0 * ( field%efield%zsigma + cell%z1 ) )
        ivac_loop1: do ivac = 1, vacuum%nvac
          z = cell%z1
          imz_loop1: do imz = 1, ncsh-1
            ! as "z" > 0 in this subroutine, the integrand is the same
            ! for both ivac -- but the integral bounds are reversed
            e_m = sinh( g * ( field%efield%zsigma + cell%z1 - z ) )
            e_p = sinh( g * ( field%efield%zsigma + cell%z1 + z ) )
            fra(ncsh-imz) = real(  rhtxy(imz,irec2-1,ivac) ) * e_m
            fia(ncsh-imz) = aimag( rhtxy(imz,irec2-1,ivac) ) * e_m
            frb(imz)      = real(  rhtxy(imz,irec2-1,ivac) ) * e_p
            fib(imz)      = aimag( rhtxy(imz,irec2-1,ivac) ) * e_p
            z = z + vacuum%delz
          end do imz_loop1
          call intgz1( fra, vacuum%delz, ncsh - 1, alpha(1,ivac,1), .false. )
          call intgz1( fia, vacuum%delz, ncsh - 1, alpha(1,ivac,2), .false. )
          call qsf( vacuum%delz, frb, beta(1,ivac,1), ncsh - 1, 1 )
          call qsf( vacuum%delz, fib, beta(1,ivac,2), ncsh - 1, 1 )
        end do ivac_loop1

        if ( ivac == 2 ) then
          ! honour reversed integral bounds
          alpha(:,ivac,1) = - alpha(:,ivac,1)
          alpha(:,ivac,2) = - alpha(:,ivac,2)
          beta(:,ivac,1)  = - beta(:,ivac,1)
          beta(:,ivac,2)  = - beta(:,ivac,2)
        end if

        alph1 = cmplx( alpha(ncsh-1,1,1), alpha(ncsh-1,1,2) )
        if ( vacuum%nvac == 1 ) then
          if ( sym%invs ) then
            alph2 = conjg( alph1 )
          else
            alph2 = alph1
          end if
        else
          alph2 = cmplx( alpha(ncsh-1,2,1), alpha(ncsh-1,2,2) )
        end if
        ivac_loop2: do ivac = 1, vacuum%nvac
          z = cell%z1
          if ( ivac == 1 ) alph0 = alph2
          if ( ivac == 2 ) alph0 = alph1
          imz_loop2: do imz = 1, ncsh-1
            betaz  = cmplx( beta(imz,ivac,1), beta(imz,ivac,2) )
            alphaz = cmplx( alpha(ncsh-imz,ivac,1), alpha(ncsh-imz,ivac,2) )
            e_m = sinh( g * ( field%efield%zsigma + cell%z1 - z ) )
            e_p = sinh( g * ( field%efield%zsigma + cell%z1 + z ) )
            test = e_m * ( alph0 + betaz ) + e_p * alphaz
            if ( 2.0 * test == test ) test = cmplx( 0.0, 0.0 )
            vxy(imz,irec2-1,ivac) = vxy(imz,irec2-1,ivac) + vcons * test
            if ( allocated( field%efield%c1 ) ) then
              e_m = exp_save( - g * z )
              e_p = exp_save(   g * z )
              if ( ivac == 1 ) then ! z > 0
                vxy(imz,irec2-1,ivac) = vxy(imz,irec2-1,ivac) &
                + field%efield%c1(irec2-1) * e_p &
                + field%efield%c2(irec2-1) * e_m
              else ! z < 0
                vxy(imz,irec2-1,ivac) = vxy(imz,irec2-1,ivac) &
                + field%efield%c1(irec2-1) * e_m &
                + field%efield%c2(irec2-1) * e_p
              end if
            end if
            z = z + vacuum%delz
          end do imz_loop2
        end do ivac_loop2
        alphm(irec2-1,1) = alph1
        alphm(irec2-1,2) = alph2

        ! ********** neumann ************************************
      else
        ivac_loop3: do ivac = 1, vacuum%nvac
          z = cell%z1
          imz_loop3: do imz = 1, vacuum%nmzxy
            e_m = exp_save( - g * z )
            e_p = exp_save(   g * z )
            fra(ip-imz) = real(  rhtxy(imz,irec2-1,ivac) ) * e_m
            fia(ip-imz) = aimag( rhtxy(imz,irec2-1,ivac) ) * e_m
            frb(imz)    = real(  rhtxy(imz,irec2-1,ivac) ) * e_p
            fib(imz)    = aimag( rhtxy(imz,irec2-1,ivac) ) * e_p
            z = z + vacuum%delz
          end do imz_loop3
          ! add external field, if segmented
          if ( allocated( field%efield%rhoef ) ) then
            z = cell%z1 + field%efield%zsigma
            e_m = exp_save( - g * z )
            e_p = exp_save(   g * z )
                ! the equation has a minus sign as "rhtxy" contains the electron density
            ! (a positive number representing a negative charge) while rhoef
            ! specifies the charges in terms of the (positive) elementary charge "e".
            fra(ip-ncsh) = fra(ip-ncsh) - real(  field%efield%rhoef(irec2-1,ivac) ) * e_m
            fia(ip-ncsh) = fia(ip-ncsh) - aimag( field%efield%rhoef(irec2-1,ivac) ) * e_m
            frb(ncsh)    = frb(ncsh)    - real(  field%efield%rhoef(irec2-1,ivac) ) * e_p
            fib(ncsh)    = fib(ncsh)    - aimag( field%efield%rhoef(irec2-1,ivac) ) * e_p
          end if
          call intgz1( fra, vacuum%delz, vacuum%nmzxy, alpha(1,ivac,1), .true. )
          call intgz1( fia, vacuum%delz, vacuum%nmzxy, alpha(1,ivac,2), .true. )
          call qsf( vacuum%delz, frb, beta(1,ivac,1), vacuum%nmzxy, 1 )
          call qsf( vacuum%delz, fib, beta(1,ivac,2), vacuum%nmzxy, 1 )
        end do ivac_loop3
        alph1 = cmplx( alpha(vacuum%nmzxy,1,1), alpha(vacuum%nmzxy,1,2) )
        if ( vacuum%nvac == 1 ) then
          if ( sym%invs ) then
            alph2 = conjg( alph1 )
          else
            alph2 = alph1
          end if
        else
          alph2 = cmplx( alpha(vacuum%nmzxy,2,1), alpha(vacuum%nmzxy,2,2) )
        end if
        ivac_loop4: do ivac = 1, vacuum%nvac
          z = cell%z1
          if ( ivac == 1 ) alph0 = alph2
          if ( ivac == 2 ) alph0 = alph1
          imz_loop4: do imz = 1, vacuum%nmzxy
            betaz  = cmplx( beta(imz,ivac,1), beta(imz,ivac,2) )
            alphaz = cmplx( alpha(ip-imz,ivac,1), alpha(ip-imz,ivac,2) )
            e_m = exp_save( - g * z )
            e_p = exp_save(   g * z )
            test = e_m * ( alph0 + betaz ) + e_p * alphaz
            if ( 2.0 * test == test ) test = cmplx( 0.0, 0.0 )
            vxy(imz,irec2-1,ivac) = vxy(imz,irec2-1,ivac) + vcons * test
            z = z + vacuum%delz
          end do imz_loop4
        end do ivac_loop4
        alphm(irec2-1,1) = alph1
        alphm(irec2-1,2) = alph2
      end if
    end do irec2_loop
  end subroutine vvacxy

  !------------------------------------------------------------------
  pure real function exp_save( x )
    ! replace exp by a function that does not under/overflow
    implicit none
    real, intent(in) :: x
    real, parameter  :: maxexp = log( 2.0 ) * maxexponent( 2.0 )
    real, parameter  :: minexp = log( 2.0 ) * minexponent( 2.0 )

    if ( abs( x ) > minexp .and. abs( x ) < maxexp ) then
      exp_save = exp( x )
    else
      if ( x > 0 ) then
        if ( x > minexp ) then
          exp_save = exp( maxexp )
        else
          exp_save = exp( minexp )
        endif
      else
        if ( - x > minexp ) then
          exp_save = exp( - maxexp )
        else
          exp_save = exp( - minexp )
        endif
      endif
    endif
  end function exp_save

end module m_vvacxy
