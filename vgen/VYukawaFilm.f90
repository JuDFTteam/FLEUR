module m_VYukawaFilm


  contains 


  subroutine VYukawaFilm( stars, vacuum, cell, sym, input, mpi, atoms, sphhar, dimension, oneD, den, &
                          VYukawa )

    use m_constants
    use m_types
    use m_psqpw
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

    !call VYukawaFilmInterstitial( stars, vacuum, cell, sym, input, &
    !                              psq, VYukawa%vacxy, VYukawa%vacz, alphm, &
    !                              VYukawa%pw(:,1) )


    ! MUFFIN-TIN POTENTIAL

    !call Vmts( input, mpi, stars, sphhar, atoms, sym, cell, oneD, VYukawa%pw(:,1), den%mt(:,0:,:,1), VYukawa%potdenType, VYukawa%mt(:,0:,:,1) )


  end subroutine VYukawaFilm


  subroutine VYukawaFilmVacuum( stars, vacuum, cell, sym, input, &
                                psq, rhtxy, rht, &
                                VVxy, VVz, alphm )

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
          ! use only stars within the g_max sphere -> stars outside the sphere
          ! have per definition index ig3n = 0
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
        fa(1:vacuum%nmzxy,irec2) = rhtxy(1:vacuum%nmzxy,irec2-1,ivac) * exp_m(1:vacuum%nmzxy,irec2)
        fb(1:vacuum%nmzxy,irec2) = rhtxy(1:vacuum%nmzxy,irec2-1,ivac) * exp_p(1:vacuum%nmzxy,irec2) ! betaz_i = beta_i = int_{z_1}^{z_i} fb(z') dz'
        call intgz1Reverse( fa(:,irec2), vacuum%delz, vacuum%nmzxy, alpha(:,irec2,ivac), .true. ) ! integrals 
        call qsfComplex( vacuum%delz, fb(:,irec2), beta(:,irec2,ivac), vacuum%nmzxy, 1 ) 
        ! alphm(:,1) = integral from D/2 to infinity in upper vacuum, over rho^+(g,z) * exp(-g*z) (rho^+ = charge density in upper vacuum)
        ! alphm(:,2) = integral from D/2 to infinity in lower vacuum, over rho^-(g,z) * exp(-g*z), where rho^-(g,z):=rho^-(g,-z) (rho^- = charge density in lower vacuum)
        alphm(irec2,ivac) = alpha(1,irec2,ivac)
      end do
      if ( vacuum%nvac == 1 ) then ! using symmetry relations, if there are any
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
        VVxy(:,irec2,ivac) = VVxy(:,irec2,ivac) + vcons(irec2) * gamma(:,irec2) 
      end do
    end do

    ! case irec2 = 1:
    exp_m(1:vacuum%nmz,1) = exp_m(1:vacuum%nmz,1) * expDhg(1)
    exp_p(1:vacuum%nmz,1) = exp_p(1:vacuum%nmz,1) * expDhg(1)
    do ivac = 1, vacuum%nvac
      ga(1:vacuum%nmz) = rht(1:vacuum%nmz,ivac) * exp_m(1:vacuum%nmz,1)
      gb(1:vacuum%nmz) = rht(1:vacuum%nmz,ivac) * exp_p(1:vacuum%nmz,1) ! betaz_i = beta_i = int_{z_1}^{z_i} fb(z') dz'
      call intgz1Reverse( ga, vacuum%delz, vacuum%nmz, delta(:,ivac), .true. ) ! integrals 
      call qsf( vacuum%delz, gb, epsilon(:,ivac), vacuum%nmz, 1 )
      alphm(1,ivac) = delta(1,ivac) ! integral from D/2 (z_1) to infty and conversion to complex
    end do
    if ( vacuum%nvac == 1 ) alphm(1,2) = alphm(1,1) ! is real, thus no conjg as in the irec2 > 1 case
    do ivac = 1, vacuum%nvac
      zeta(1:vacuum%nmz) = exp_m(1:vacuum%nmz,1) * ( alphm(1,mod(ivac,2)+1) + epsilon(1:vacuum%nmz,ivac) ) &
                         + exp_p(1:vacuum%nmz,1) *                              delta(1:vacuum%nmz,ivac) 
      where ( 2. * zeta == zeta ) zeta = 0.
      VVz(:,ivac) = VVz(:,ivac) + vcons(1) * zeta(:) 
    end do


  end subroutine VYukawaFilmVacuum



end module m_VYukawaFilm
