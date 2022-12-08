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

    use m_types
    use m_constants
    use m_intgr, only: intgz1
    use m_qsf

    implicit none

    type(t_input),  intent(in)     :: input
    type(t_field),  intent(in)     :: field
    type(t_vacuum), intent(in)     :: vacuum
    type(t_sym),    intent(in)     :: sym
    type(t_stars),  intent(in)     :: stars
    type(t_cell),   intent(in)     :: cell
    complex,        intent(in)     :: rhtxy(vacuum%nmzxyd,stars%ng2-1,2)
    complex,        intent(inout)  :: vxy(vacuum%nmzxyd,stars%ng2-1,2)
    complex,        intent(out)    :: alphm(stars%ng2,2)

    complex                        :: alph0, alphaz, betaz, test, phas
    real                           :: g, vcons, z, e_m, e_p,arg
    integer                        :: imz, ip, irec2, ivac, ncsh,irec2r
    logical                        :: tail
    real                           :: fra(vacuum%nmzxyd), frb(vacuum%nmzxyd), fia(vacuum%nmzxyd), fib(vacuum%nmzxyd)
    real                           :: alpha(vacuum%nmzxyd,2,2,stars%ng2), beta(vacuum%nmzxyd,2,2,stars%ng2)
    real, allocatable              :: sig_top(:), sig_bot(:)
    intrinsic aimag, cmplx, conjg, exp, real

    if ( allocated( field%efield%rhoef ) .or. field%efield%dirichlet ) then
      ncsh = field%efield%zsigma / vacuum%delz + 1.01
      ! if nmzxy < ncsh, the inhomogenous field cannot be represented.
      ! if nmzxy > ncsh, the boundary condition is wrong - and the
      ! potential is very wavy - try it yourself, if you don't believe.
      if ( vacuum%nmzxyd < ncsh ) then
        write (oUnit,*) 'error vvacxy.f: vacuum%nmzxyd =', vacuum%nmzxyd, '< ncsh(zsigma) = ', ncsh
        call judft_error( "error: vacuum%nmzxyd < ncsh", calledby="vvacxy" )
      else if ( vacuum%nmzxyd > ncsh ) then
        write (oUnit,*) 'warning vvacxy.f: vacuum%nmzxyd =', vacuum%nmzxyd, '> ncsh(zsigma) = ', ncsh
        write (0,*) 'warning vvacxy.f: vacuum%nmzxyd =', vacuum%nmzxyd, '> ncsh(zsigma) = ', ncsh
        call judft_warn( "nmzxyd > ncsh", calledby="vvacxy" )
      end if
    end if

    ! 2-dim star loop g/=0
    ip = vacuum%nmzxy + 1
!    irec2_loop: do irec2 = 2, stars%ng2
!      g = stars%sk2(irec2)
!      vcons = tpi_const / g
      ! ********** dirichlet ************************************
      if ( field%efield%dirichlet ) then
        do irec2 = 2, stars%ng2
          g = stars%sk2(irec2)
          vcons = tpi_const / g
          if ( allocated( field%efield%rhoef ) ) then
            vxy(ncsh:vacuum%nmzxy,irec2-1,1) = field%efield%rhoef(irec2-1, 1)
            if ( vacuum%nvac == 2 ) then
              vxy(ncsh:vacuum%nmzxy,irec2-1,2) = field%efield%rhoef(irec2-1, 2)
            end if
          else
            vxy(ncsh:vacuum%nmzxy,irec2-1,1:vacuum%nvac) = 0.0
          end if
          vcons = 2.0 * vcons / sinh( g * 2.0 * ( field%efield%zsigma + cell%z1 ) )
          do ivac = 1, vacuum%nvac
            z = cell%z1
            do imz = 1, ncsh-1
            ! as "z" > 0 in this subroutine, the integrand is the same
            ! for both ivac -- but the integral bounds are reversed
              e_m = sinh( g * ( field%efield%zsigma + cell%z1 - z ) )
              e_p = sinh( g * ( field%efield%zsigma + cell%z1 + z ) )
              fra(ncsh-imz) = real(  rhtxy(imz,irec2-1,ivac) ) * e_m
              fia(ncsh-imz) = aimag( rhtxy(imz,irec2-1,ivac) ) * e_m
              frb(imz)      = real(  rhtxy(imz,irec2-1,ivac) ) * e_p
              fib(imz)      = aimag( rhtxy(imz,irec2-1,ivac) ) * e_p
              z = z + vacuum%delz
            end do 
            call intgz1( fra, vacuum%delz, ncsh - 1, alpha(1,ivac,1,irec2), .false. )
            call intgz1( fia, vacuum%delz, ncsh - 1, alpha(1,ivac,2,irec2), .false. )
            call qsf( vacuum%delz, frb, beta(1,ivac,1,irec2), ncsh - 1, 1 )
            call qsf( vacuum%delz, fib, beta(1,ivac,2,irec2), ncsh - 1, 1 )
          end do 

          if ( ivac == 2 ) then
          ! honour reversed integral bounds
            alpha(:,ivac,:,irec2) = - alpha(:,ivac,:,irec2)
            beta(:,ivac,:,irec2)  = - beta(:,ivac,:,irec2)
          end if
        enddo
        do irec2 = 2, stars%ng2
          g = stars%sk2(irec2)
          vcons = tpi_const / g
          alphm(irec2-1,1) = cmplx( alpha(vacuum%nmzxy,1,1,irec2), alpha(vacuum%nmzxy,1,2,irec2) )
          if ( vacuum%nvac == 1 ) then
            call stars%map_2nd_vac(vacuum,irec2,irec2r,phas)
            alphm(irec2r-1,2) = phas * cmplx(alpha(vacuum%nmzxy,1,1,irec2),alpha(vacuum%nmzxy,1,2,irec2))
          else
            alphm(irec2-1,2) = cmplx( alpha(vacuum%nmzxy,2,1,irec2), alpha(vacuum%nmzxy,2,2,irec2) )
          end if
          do ivac = 1, vacuum%nvac
            z = cell%z1
            if ( ivac == 1 ) alph0 = alphm(irec2-1,2)
            if ( ivac == 2 ) alph0 = alphm(irec2-1,1)
            do imz = 1, ncsh-1
              betaz  = cmplx( beta(imz,ivac,1,irec2), beta(imz,ivac,2,irec2) )
              alphaz = cmplx( alpha(ncsh-imz,ivac,1,irec2), alpha(ncsh-imz,ivac,2,irec2) )
              e_m = sinh( g * ( field%efield%zsigma + cell%z1 - z ) )
              e_p = sinh( g * ( field%efield%zsigma + cell%z1 + z ) )
              test = e_m * ( alph0 + betaz ) + e_p * alphaz
              if ( 2.0 * test == test ) test = cmplx( 0.0, 0.0 )
              vxy(imz,irec2-1,ivac) = vxy(imz,irec2-1,ivac) + vcons * test
              if ( allocated( field%efield%c1 ) ) then
                e_m = exp_save( - g * z )
                e_p = exp_save(   g * z )
                if ( ivac == 1 ) then ! z > 0
                  vxy(imz,irec2-1,ivac) = vxy(imz,irec2-1,ivac) + field%efield%c1(irec2-1) * e_p + field%efield%c2(irec2-1) * e_m
                else ! z < 0
                  vxy(imz,irec2-1,ivac) = vxy(imz,irec2-1,ivac) + field%efield%c1(irec2-1) * e_m + field%efield%c2(irec2-1) * e_p
                end if
              end if
              z = z + vacuum%delz
            end do 
          end do 
        enddo 

        ! ********** neumann ************************************
      else
        do irec2 = 2, stars%ng2
          g = stars%sk2(irec2)
          vcons = tpi_const / g
        
          do ivac = 1, vacuum%nvac
            z = cell%z1
            do imz = 1, vacuum%nmzxy
              e_m = exp_save( - g * z )
              e_p = exp_save(   g * z )
              fra(ip-imz) = real(  rhtxy(imz,irec2-1,ivac) ) * e_m
              fia(ip-imz) = aimag( rhtxy(imz,irec2-1,ivac) ) * e_m
              frb(imz)    = real(  rhtxy(imz,irec2-1,ivac) ) * e_p
              fib(imz)    = aimag( rhtxy(imz,irec2-1,ivac) ) * e_p
              z = z + vacuum%delz
            end do 
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
            call intgz1( fra, vacuum%delz, vacuum%nmzxy, alpha(1,ivac,1,irec2), .true. )
            call intgz1( fia, vacuum%delz, vacuum%nmzxy, alpha(1,ivac,2,irec2), .true. )
            call qsf( vacuum%delz, frb, beta(1,ivac,1,irec2), vacuum%nmzxy, 1 )
            call qsf( vacuum%delz, fib, beta(1,ivac,2,irec2), vacuum%nmzxy, 1 )
          end do 
        enddo
        do irec2 = 2, stars%ng2
          g = stars%sk2(irec2)
          vcons = tpi_const / g
        
          alphm(irec2-1,1) = cmplx( alpha(vacuum%nmzxy,1,1,irec2), alpha(vacuum%nmzxy,1,2,irec2) )
          if ( vacuum%nvac == 1 ) then
            call stars%map_2nd_vac(vacuum,irec2,irec2r,phas)
            alphm(irec2r-1,2) = phas * cmplx(alpha(vacuum%nmzxy,1,1,irec2),alpha(vacuum%nmzxy,1,2,irec2))
          else
            alphm(irec2-1,2) = cmplx( alpha(vacuum%nmzxy,2,1,irec2), alpha(vacuum%nmzxy,2,2,irec2) )
          end if
          do ivac = 1, vacuum%nvac
            z = cell%z1
            if ( ivac == 1 ) alph0 = alphm(irec2-1,2)
            if ( ivac == 2 ) alph0 = alphm(irec2-1,1)
            do imz = 1, vacuum%nmzxy
              betaz  = cmplx( beta(imz,ivac,1,irec2), beta(imz,ivac,2,irec2) )
              alphaz = cmplx( alpha(ip-imz,ivac,1,irec2), alpha(ip-imz,ivac,2,irec2) )
              e_m = exp_save( - g * z )
              e_p = exp_save(   g * z )
              test = e_m * ( alph0 + betaz ) + e_p * alphaz
              if ( 2.0 * test == test ) test = cmplx( 0.0, 0.0 )
              vxy(imz,irec2-1,ivac) = vxy(imz,irec2-1,ivac) + vcons * test
              z = z + vacuum%delz
            end do 
          end do 
        enddo
      end if
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
