module m_VYukawaFilm

  ! computation of the film-case Yukawa potential for the preconditioning of the charge density residual

  contains 



  subroutine VYukawaFilm( stars, vacuum, cell, sym, input, mpi, atoms, sphhar, dimension, oneD, den, &
                          VYukawa )

    use m_constants
    use m_types
    use m_psqpw
    use m_vmts
    implicit none

    type(t_stars),      intent(in)    :: stars
    type(t_vacuum),     intent(in)    :: vacuum
    type(t_cell),       intent(in)    :: cell
    type(t_sym),        intent(in)    :: sym
    type(t_input),      intent(in)    :: input
    type(t_mpi),        intent(in)    :: mpi
    type(t_atoms),      intent(in)    :: atoms 
    type(t_sphhar),     intent(in)    :: sphhar
    type(t_dimension),  intent(in)    :: dimension
    type(t_oneD),       intent(in)    :: oneD
    type(t_potden),     intent(in)    :: den

    type(t_potden),     intent(inout) :: VYukawa

    complex                           :: psq(stars%ng3)
    complex                           :: alphm(stars%ng2,2)


    ! PSEUDO-CHARGE DENSITY

    call psqpw( mpi, atoms, sphhar, stars, vacuum, dimension, cell, input, sym, oneD, den%pw(:,1), den%mt(:,:,:,1), den%vacz(:,:,1), .false., VYukawa%potdenType, psq )


    ! VACUUM POTENTIAL

    call VYukawaFilmVacuum( stars, vacuum, cell, sym, input, &
                            psq, den%vacxy(:,:,:,1), den%vacz(:,:,1), &
                            VYukawa%vacxy, VYukawa%vacz, alphm )


    ! INTERSTITIAL POTENTIAL

    call VYukawaFilmInterstitial( stars, vacuum, cell, sym, input, &
                                  psq, VYukawa%vacxy, VYukawa%vacz, alphm, &
                                  VYukawa%pw(:,1) )


    ! MUFFIN-TIN POTENTIAL

    call Vmts( input, mpi, stars, sphhar, atoms, sym, cell, oneD, VYukawa%pw(:,1), den%mt(:,0:,:,1), VYukawa%potdenType, VYukawa%mt(:,0:,:,1) )


  end subroutine VYukawaFilm



  subroutine VYukawaFilmVacuum( stars, vacuum, cell, sym, input, &
                                psq, rhtxy, rht, &
                                VVxy, VVz, alphm )

    ! 1. part: compute the contribution from the interstitial charge density to the vacuum potential as a function of q_xy and z
    ! 2. part: compute the contribution from the vacuum charge density to the vacuum potential as a function of q_xy and z by numerical integration

    use m_ExpSave
    use m_constants
    use m_types
    use m_intgr, only: intgz1Reverse
    use m_qsf
    implicit none

    type(t_stars),  intent(in)  :: stars
    type(t_vacuum), intent(in)  :: vacuum
    type(t_cell),   intent(in)  :: cell
    type(t_sym),    intent(in)  :: sym
    type(t_input),  intent(in)  :: input
    complex,        intent(in)  :: psq(stars%ng3)
    complex,        intent(in)  :: rhtxy(vacuum%nmzxyd,stars%ng2-1,2)
    real,           intent(in)  :: rht(vacuum%nmzd,2) 

    complex,        intent(out) :: VVxy(vacuum%nmzxyd,2:stars%ng2,2) ! this is the qxy /= 0 part of the vacuum potential
    real,           intent(out) :: VVz(vacuum%nmzd,2)                ! this is the qxy  = 0 part of the vacuum potential
    complex,        intent(out) :: alphm(stars%ng2,2)                ! these are the integrals in upper and lower vacuum, now including the first star---integral for star ig2 is in alphm(ig2,ivac) 

    complex                     :: sum_qz(2,stars%ng2)
    complex                     :: c_ph(-stars%mx3:stars%mx3,stars%ng2)
    complex                     :: signIqz
    real                        :: g_damped(stars%ng2), qz, sign, vcons(stars%ng2) 
    real                        :: exp_m(vacuum%nmzd,stars%ng2), exp_p(vacuum%nmzd,stars%ng2)
    real                        :: expDhg(stars%ng2), expDg(stars%ng2)
    real                        :: z(vacuum%nmzd)
    integer                     :: iz, irec2, irec3, ivac, kz
    complex                     :: fa(vacuum%nmzxyd,2:stars%ng2), fb(vacuum%nmzxyd,2:stars%ng2)
    complex                     :: alpha(vacuum%nmzxyd,2:stars%ng2,2), beta(vacuum%nmzxyd,2:stars%ng2,2), gamma(vacuum%nmzxyd,2:stars%ng2)
    real                        :: ga(vacuum%nmzd), gb(vacuum%nmzd)
    real                        :: delta(vacuum%nmzd,2), epsilon(vacuum%nmzd,2), zeta(vacuum%nmzd)

    intrinsic aimag, cmplx, conjg, exp, real, mod


    ! definitions / initialisations:
    do iz = 1, vacuum%nmz
      z(iz) = ( iz - 1 ) * vacuum%delz ! z_i starting at 0
    end do
    do irec2 = 1, stars%ng2
      g_damped(irec2) = sqrt( stars%sk2(irec2) ** 2 + input%preconditioning_param ** 2 )
      vcons(irec2) = tpi_const / g_damped(irec2)
      do iz = 1, vacuum%nmz
        exp_m(iz,irec2) = exp_save( - g_damped(irec2) * z(iz) ) ! redefined in contribution of vacuum charge density part
        exp_p(iz,irec2) = exp_save(   g_damped(irec2) * z(iz) )
      end do
      expDhg(irec2) = exp_save( - cell%z1 * g_damped(irec2) )
      expDg(irec2)  = expDhg(irec2) ** 2
    end do
    sum_qz = (0.,0.)
    VVxy = (0.,0.)
    VVz = 0.
   

    ! CONTRIBUTION FROM THE INTERSTITIAL CHARGE DENSITY

    do irec2 = 1, stars%ng2
      do ivac = 1, vacuum%nvac
        sign = 3. - 2. * ivac
        do kz = -stars%mx3, stars%mx3
          irec3 = stars%ig(stars%kv2(1,irec2),stars%kv2(2,irec2),kz)
          ! use only stars within the g_max sphere -> stars outside the sphere have per definition index ig3n = 0
          if ( irec3 /= 0 ) then
            c_ph(kz,irec2) = stars%rgphs(stars%kv2(1,irec2),stars%kv2(2,irec2),kz)
            qz = kz * cell%bmat(3,3)
            signIqz = sign * ImagUnit * qz
            sum_qz(ivac,irec2) = sum_qz(ivac,irec2) + c_ph(kz,irec2) * psq(irec2) * ( exp( signIqz * cell%z1 ) - exp( - signIqz * cell%z1 ) * expDg(irec2) ) / ( signIqz + g_damped(irec2) )
          endif
        enddo
        if( irec2 /= 1 ) then
          VVxy(1:vacuum%nmzxy,irec2,ivac) = vcons(irec2) * sum_qz(ivac,irec2) * exp_m(1:vacuum%nmzxy,irec2)
        else
          VVz(1:vacuum%nmz,ivac)          = vcons(1)     * sum_qz(ivac,1)     * exp_m(1:vacuum%nmz,1)
        end if
      enddo
    enddo


    ! CONTRIBUTION FROM THE VACUUM CHARGE DENSITY

    ! case irec2 > 1:
    do irec2 = 2, stars%ng2
      exp_m(1:vacuum%nmzxy,irec2) = exp_m(1:vacuum%nmzxy,irec2) * expDhg(irec2) ! z_i can now be thought of as starting from D/2 -> expDhg is the correction
      exp_p(1:vacuum%nmzxy,irec2) = exp_p(1:vacuum%nmzxy,irec2) * expDhg(irec2)
      do ivac = 1, vacuum%nvac
        ! integrands:
        fa(1:vacuum%nmzxy,irec2) = rhtxy(1:vacuum%nmzxy,irec2-1,ivac) * exp_m(1:vacuum%nmzxy,irec2)
        fb(1:vacuum%nmzxy,irec2) = rhtxy(1:vacuum%nmzxy,irec2-1,ivac) * exp_p(1:vacuum%nmzxy,irec2)
        ! integrals:
        ! alpha(z,q_xy,ivac) = int_z^infty rho(z',q_xy,ivac) exp(-sqrt(q_xy**2+prec_param**2)*z') dz'
        ! beta (z,q_xy,ivac) = int_{D/2}^z rho(z',q_xy,ivac) exp(+sqrt(q_xy**2+prec_param**2)*z') dz'
        ! where for z < 0 the lower vacuum charge density (ivac=2) is defined by rho(q_xy,z,ivac=2) := rho(q_xy,-z,ivac=2)
        call intgz1Reverse( fa(:,irec2), vacuum%delz, vacuum%nmzxy, alpha(:,irec2,ivac), .true. ) 
        call qsfComplex( vacuum%delz, fb(:,irec2), beta(:,irec2,ivac), vacuum%nmzxy, 1 )
        ! alphm(q_xy,ivac) = alpha(D/2,q_xy,ivac) --- these integrals are also needed for the interstitial potential
        alphm(irec2,ivac) = alpha(1,irec2,ivac)
      end do
      if ( vacuum%nvac == 1 ) then
        if ( sym%invs ) then
          alphm(irec2,2) = conjg( alphm(irec2,1) )
        else
          alphm(irec2,2) = alphm(irec2,1)
        end if
      end if
      do ivac = 1, vacuum%nvac
        gamma(1:vacuum%nmzxy,irec2) = exp_m(1:vacuum%nmzxy,irec2) * ( alphm(irec2,mod(ivac,2)+1) + beta(1:vacuum%nmzxy,irec2,ivac) ) &
                                    + exp_p(1:vacuum%nmzxy,irec2) *                               alpha(1:vacuum%nmzxy,irec2,ivac) ! mod(ivac,2)+1 outputs the other ivac value 
        where ( 2. * gamma(:,irec2) == gamma(:,irec2) ) gamma(:,irec2) = cmplx( 0., 0. )
        VVxy(1:vacuum%nmzxy,irec2,ivac) = VVxy(1:vacuum%nmzxy,irec2,ivac) + vcons(irec2) * gamma(1:vacuum%nmzxy,irec2) 
      end do
    end do

    ! case irec2 = 1:
    exp_m(1:vacuum%nmz,1) = exp_m(1:vacuum%nmz,1) * expDhg(1)
    exp_p(1:vacuum%nmz,1) = exp_p(1:vacuum%nmz,1) * expDhg(1)
    do ivac = 1, vacuum%nvac
      ga(1:vacuum%nmz) = rht(1:vacuum%nmz,ivac) * exp_m(1:vacuum%nmz,1)
      gb(1:vacuum%nmz) = rht(1:vacuum%nmz,ivac) * exp_p(1:vacuum%nmz,1)
      call intgz1Reverse( ga(:), vacuum%delz, vacuum%nmz, delta(:,ivac), .true. ) ! integrals 
      call qsf( vacuum%delz, gb(:), epsilon(:,ivac), vacuum%nmz, 1 )
      alphm(1,ivac) = delta(1,ivac)
    end do
    if ( vacuum%nvac == 1 ) alphm(1,2) = alphm(1,1) ! is real, thus no conjg as in the irec2 > 1 case
    do ivac = 1, vacuum%nvac
      zeta(1:vacuum%nmz) = exp_m(1:vacuum%nmz,1) * ( alphm(1,mod(ivac,2)+1) + epsilon(1:vacuum%nmz,ivac) ) &
                         + exp_p(1:vacuum%nmz,1) *                              delta(1:vacuum%nmz,ivac) 
      where ( 2. * zeta == zeta ) zeta = 0.
      VVz(1:vacuum%nmz,ivac) = VVz(1:vacuum%nmz,ivac) + vcons(1) * zeta(1:vacuum%nmz) 
    end do


  end subroutine VYukawaFilmVacuum



  subroutine VYukawaFilmInterstitial( stars, vacuum, cell, sym, input, &
                                      psq, VVxy, VVz, alphm, &
                                      VIq )

    ! 1. part: compute the contribution from the interstitial charge density to the interstitial potential (largest part) as a function of q_xy and z
    ! 2. part: add the contribution from the vacuum charge density to the interstitial potential, which had already been computed earlier for the vacuum potential
    ! 3. part: compute the coefficients V^I(q_xy,q_z) from the function V^I(q_xy,z) by a 1D Fourier transform: V^I(q_xy,z) = sum_{q_z} V^I(q_xy,q_z) * exp( ImagUnit * q_z * z )

    use m_ExpSave
    use m_constants
    use m_types
    use m_cfft
    implicit none

    type(t_stars),  intent(in)  :: stars
    type(t_vacuum), intent(in)  :: vacuum
    type(t_cell),   intent(in)  :: cell
    type(t_sym),    intent(in)  :: sym
    type(t_input),  intent(in)  :: input
    complex,        intent(in)  :: psq(stars%ng3)
    complex,        intent(in)  :: VVxy(vacuum%nmzxyd,2:stars%ng2,2)
    real,           intent(in)  :: VVz(vacuum%nmzd,2)
    complex,        intent(in)  :: alphm(stars%ng2,2)

    complex,        intent(out) :: VIq(stars%ng3)

    real                        :: partitioning, rz, qz, fit, q
    integer                     :: irec2, irec3, iz, jz, ivac, kz, ifft
    complex                     :: VIz(3*stars%mx3,stars%ng2), sum_qz(3*stars%mx3,stars%ng2), eta(3*stars%mx3,stars%ng2)
    complex                     :: VIqz(-stars%mx3:stars%mx3,stars%ng2), c_ph(-stars%mx3:stars%mx3,stars%ng2)
    complex                     :: vcons1(stars%ng3)
    real                        :: VIzReal(3*stars%mx3,stars%ng2), VIzImag(3*stars%mx3,stars%ng2), exp_m(3*stars%mx3,stars%ng2), exp_p(3*stars%mx3,stars%ng2)
    real                        :: z(3*stars%mx3)
    real                        :: g_damped(stars%ng2), vcons2(stars%ng2)

    intrinsic abs, cmplx, conjg, cos, exp, sin, sqrt


    ! definitions / initialisations
    ifft = 3 * stars%mx3
    partitioning = 1. / real( ifft )
    do iz = 1, ifft
      z(iz) = cell%amat(3,3) * ( iz - 1 ) * partitioning               ! z_i are equidistantly distributed along ( 0, cell%amat(3,3) ]
      if( z(iz) > cell%amat(3,3) / 2. ) z(iz) = z(iz) - cell%amat(3,3) ! z_i are equidistantly distributed along ( -cell%amat(3,3) / 2, cell%amat(3,3) / 2 ]
    end do
    do irec2 = 1, stars%ng2
      g_damped(irec2) = sqrt( stars%sk2(irec2) ** 2 + input%preconditioning_param ** 2 )
      vcons2(irec2) = -1. / ( 2. * g_damped(irec2) )
      do iz = 1, ifft
        exp_m(iz,irec2) = exp_save( - g_damped(irec2) * ( cell%z1 + z(iz) ) )
        exp_p(iz,irec2) = exp_save( - g_damped(irec2) * ( cell%z1 - z(iz) ) )
      end do
    end do
    do irec3 = 1, stars%ng3
      vcons1(irec3) = fpi_const * psq(irec3) / ( stars%sk3(irec3) ** 2 + input%preconditioning_param ** 2 )
    end do
    sum_qz = (0.,0.)
    VIz = (0.,0.)
    VIq = (0.,0.)


    ! CONTRIBUTION FROM THE INTERSTITIAL CHARGE DENSITY

    ! compute V^I(q_xy,z) as a function of q_xy and z
    do irec2 = 1, stars%ng2
      do iz = 1, ifft

        ! in the transition region ( D/2, D~/2 ), smooth out the potential numerically
        if ( abs( z(iz) ) >= cell%z1 ) then
          ivac = 1 ! upper vacuum
          if ( z(iz) < 0.0 .and. .not. ( sym%invs .or. sym%zrfs ) ) ivac = 2 ! lower vacuum
          rz = ( abs( z(iz) ) - cell%z1 ) / vacuum%delz + 1.0 ! numbering the grid points in the vacuum region
          jz = rz
          q = rz - jz
          if ( irec2 == 1 ) then
            fit = 0.5     * ( q - 1. ) * ( q - 2. ) * VVz(jz,  ivac) &
                -       q *              ( q - 2. ) * VVz(jz+1,ivac) &
                + 0.5 * q * ( q - 1. )              * VVz(jz+2,ivac)
            VIz(iz,irec2) = cmplx( fit, 0.0 )
          else if ( jz + 2 <= vacuum%nmzxy ) then
            VIz(iz,irec2) = 0.5 *     ( q - 1. ) * ( q - 2. ) * VVxy(jz,  irec2-1,ivac) &
                          -       q              * ( q - 2. ) * VVxy(jz+1,irec2-1,ivac) &
                          + 0.5 * q * ( q - 1. )              * VVxy(jz+2,irec2-1,ivac)
            if ( ( sym%invs .and. .not. sym%zrfs ) .and. z(iz) < 0 ) VIz(iz,irec2) = conjg( VIz(iz,irec2) )
          end if

        ! z in ( -D/2, D/2 )
        else
          do kz = -stars%mx3, stars%mx3
            irec3 = stars%ig(stars%kv2(1,irec2),stars%kv2(2,irec2),kz)
            if ( irec3 /= 0 ) then ! use only stars within the g_max sphere
              c_ph(kz,irec2) = stars%rgphs(stars%kv2(1,irec2),stars%kv2(2,irec2),kz)
              qz = kz * cell%bmat(3,3)
              sum_qz(iz,irec2) = sum_qz(iz,irec2) + exp( ImagUnit * qz * z(iz) ) &
                               + vcons2(irec2) * ( ( g_damped(irec2) + ImagUnit * qz ) * exp_p(iz,irec2) * exp(   ImagUnit * qz * cell%z1 ) &
                                                 + ( g_damped(irec2) - ImagUnit * qz ) * exp_m(iz,irec2) * exp( - ImagUnit * qz * cell%z1 ) )
              sum_qz(iz,irec2) = vcons1(irec3) * c_ph(kz,irec2) * sum_qz(iz,irec2)
            end if
          enddo
          VIz(iz,irec2) = VIz(iz,irec2) + sum_qz(iz,irec2)
        end if
      end do


      ! CONTRIBUTION FROM THE VACUUM CHARGE DENSITY

      eta(:,irec2) = exp_m(:,irec2) * alphm(irec2,2) + exp_p(:,irec2) * alphm(irec2,1)
      where ( 2.0 * eta(:,irec2) == eta(:,irec2) ) eta(:,irec2) = cmplx( 0.0, 0.0 )
      VIz(:,irec2) = VIz(:,irec2) + tpi_const / g_damped(irec2) * eta(:,irec2)


      ! 1D FOURIER TRANSFORM TO FIND THE COEFFICIENTS V^I(q_xy,q_z)

      ! V^I(q_xy,z) = sum_{q_z} V^I(q_xy,q_z) * exp( ImagUnit * q_z * z )
      VIzReal =  real( VIz )
      VIzImag = aimag( VIz )
      call cfft( VIzReal(:,irec2), VIzImag(:,irec2), ifft, ifft, ifft, -1 )
      ! reorder:
      VIqz(0,irec2) = cmplx( VIzReal(1,irec2), VIzImag(1,irec2) )
      do kz = 1, stars%mx3
        VIqz( kz,irec2) = cmplx( VIzReal(kz+1,     irec2), VIzImag(kz+1,     irec2) )
        VIqz(-kz,irec2) = cmplx( VIzReal(ifft+1-kz,irec2), VIzImag(ifft+1-kz,irec2) )
      end do
      ! add the computed components to V^I(q_xy,q_z):
      do kz= -stars%mx3, stars%mx3
        irec3 = stars%ig(stars%kv2(1,irec2),stars%kv2(2,irec2),kz)
        VIq(irec3) = VIq(irec3) + VIqz(kz,irec2) / partitioning * stars%nstr(irec3) / stars%nstr2(irec2)
      end do
    end do


  end subroutine VYukawaFilmInterstitial



end module m_VYukawaFilm
