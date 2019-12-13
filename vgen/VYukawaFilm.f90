module m_VYukawaFilm

  ! Computation of the film-case Yukawa potential for the preconditioning of the
  ! residual charge density in 5 steps:
  ! 1. pseudo-charge density generation
  ! 2. vacuum potential generation
  ! 3. interstitial potential generation
  ! 4. muffin-tin potential generation
  ! 5. modification for charge neutrality

  ! The Yukawa potential is the solution to the modified Helmholtz equation
  ! ( Delta - lambda^2 ) V_lambda = -4 pi ( rho_out - rho_in )
  ! subject to some conditions. 
  ! The general scheme (steps 1 to 4) is the same as for the Poisson equation --
  ! we use Green function methods for the z-dependent vacuum and interstitial
  ! potentials as well as for the muffin-tin potential and apply Weinert's
  ! method.
  ! You can choose between two variants:
  ! 1. variant: 
  ! zero Dirichlet boundary conditions at +/- infinity; 
  ! multiplication with a decaying exponential in vacuum;
  ! modification in the film for charge neutrality (step 5)
  ! 2. variant:
  ! zero Dirichlet boundary conditions near the film boundary in vacuum (D/2+2R);
  ! modification in the film for charge neutrality (step 5) 
  ! In both cases charge neutrality is broken.
  ! To restore charge neutrality, we need the integral over the potential to be
  ! zero.
  ! In step 5 we therefore solve the modified Helmholtz equation again with
  ! constant right-hand side, for an additive correction to the potential.
  ! The constant is chosen such that the integral over the final potential is
  ! zero.

  contains 



  subroutine VYukawaFilm( stars, vacuum, cell, sym, input, mpi, atoms, sphhar, oneD, noco, den, &
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
    type(t_oneD),       intent(in)    :: oneD
    type(t_noco),       intent(in)    :: noco
    type(t_potden),     intent(inout) :: den

    type(t_potden),     intent(inout) :: VYukawa

    complex                           :: psq(stars%ng3)
    complex                           :: alphm(stars%ng2,2)
    real                              :: dh_prec
    real                              :: coshdh(stars%ng2)

 
    ! PSEUDO-CHARGE DENSITY

    call psqpw( mpi, atoms, sphhar, stars, vacuum, cell, input, sym, oneD, &
                den%pw(:,1), den%mt(:,:,:,1), den%vacz(:,:,1), .false., VYukawa%potdenType, &
                psq )


    ChooseVariant: if ( .true. ) then

      ! VACUUM POTENTIAL

      call VYukawaFilmVacuumVariant1( &
              stars, vacuum, cell, sym, input, atoms%rmt(1), &
              psq, den%vacxy(:,:,:,1), den%vacz(:,:,1), &
              VYukawa%vacxy, VYukawa%vacz, alphm )


      ! INTERSTITIAL POTENTIAL

      call VYukawaFilmInterstitialVariant1( &
              stars, vacuum, cell, sym, input, &
              psq, VYukawa%vacxy, VYukawa%vacz, alphm, &
              VYukawa%pw(:,1) )

    else ChooseVariant

      ! VACUUM POTENTIAL

      call VYukawaFilmVacuumVariant2( &
              stars, vacuum, cell, sym, input, 2*atoms%rmt(1), &
              psq, den%vacxy(:,:,:,1), den%vacz(:,:,1), &
              VYukawa%vacxy, VYukawa%vacz, alphm, dh_prec, coshdh )


      ! INTERSTITIAL POTENTIAL

      call VYukawaFilmInterstitialVariant2( &
              stars, vacuum, cell, sym, input, &
              psq, VYukawa%vacxy, VYukawa%vacz, alphm, dh_prec, coshdh, &
              VYukawa%pw(:,1) )

    end if ChooseVariant


    ! MUFFIN-TIN POTENTIAL

    call Vmts( input, mpi, stars, sphhar, atoms, sym, cell, oneD, .FALSE., &
               VYukawa%pw(:,1), den%mt(:,0:,:,1), VYukawa%potdenType, &
               VYukawa%mt(:,0:,:,1) )

 
    ! MODIFICATION FOR CHARGE NEUTRALITY

    call VYukawaModify( stars, vacuum, cell, sym, input, mpi, atoms, sphhar, oneD, noco, &
                        den, &
                        VYukawa )


  end subroutine VYukawaFilm



  subroutine VYukawaFilmVacuumVariant1( &
                stars, vacuum, cell, sym, input, rmt, &
                psq, rhtxy, rht, &
                VVxy, VVz, alphm )

    ! 1. part: Compute the contribution from the interstitial charge density to the vacuum potential as a function of q_xy and z (analytic expression for integral)
    ! 2. part: Compute the contribution from the vacuum charge density to the vacuum potential as a function of q_xy and z by numerical integration

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
    real,           intent(in)  :: rmt
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
    integer                     :: iz, irec2, irec3, ivac, iqz
    complex                     :: fa(vacuum%nmzxyd,2:stars%ng2), fb(vacuum%nmzxyd,2:stars%ng2)
    complex                     :: alpha(vacuum%nmzxyd,2:stars%ng2,2), beta(vacuum%nmzxyd,2:stars%ng2,2), gamma(vacuum%nmzxyd,2:stars%ng2)
    real                        :: ga(vacuum%nmzd), gb(vacuum%nmzd)
    real                        :: delta(vacuum%nmzd,2), epsilon(vacuum%nmzd,2), zeta(vacuum%nmzd)


    ! DEFINITIONS / ALLOCATIONS / INITIALISATIONS    
    
    do iz = 1, vacuum%nmz
      z(iz) = ( iz - 1 ) * vacuum%delz
    end do
    do irec2 = 1, stars%ng2
      g_damped(irec2) = sqrt( stars%sk2(irec2) ** 2 + input%preconditioning_param ** 2 )
      vcons(irec2) = tpi_const / g_damped(irec2)
      do iz = 1, vacuum%nmz
        exp_m(iz,irec2) = exp_save( - g_damped(irec2) * z(iz) )
        exp_p(iz,irec2) = exp_save(   g_damped(irec2) * z(iz) )
      end do
      expDhg(irec2) = exp_save( - cell%z1 * g_damped(irec2) )
      expDg(irec2)  = exp_save( -2 * cell%z1 * g_damped(irec2) )
    end do
    sum_qz = (0.,0.)
    VVxy = (0.,0.)
    VVz = 0.
   

    ! CONTRIBUTION FROM THE INTERSTITIAL CHARGE DENSITY

    do irec2 = 1, stars%ng2
      do ivac = 1, vacuum%nvac
        sign = 3. - 2. * ivac
        do iqz = -stars%mx3, stars%mx3
          irec3 = stars%ig(stars%kv2(1,irec2),stars%kv2(2,irec2),iqz)
          ! use only stars within the g_max sphere -> stars outside the sphere have per definition index ig3n = 0
          if ( irec3 /= 0 ) then
            c_ph(iqz,irec2) = stars%rgphs(stars%kv2(1,irec2),stars%kv2(2,irec2),iqz)
            qz = iqz * cell%bmat(3,3)
            signIqz = sign * ImagUnit * qz 
            sum_qz(ivac,irec2) = sum_qz(ivac,irec2) + c_ph(iqz,irec2) * psq(irec3) * ( exp( signIqz * cell%z1 ) - exp( - signIqz * cell%z1 ) * expDg(irec2) ) / ( signIqz + g_damped(irec2) )
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

    ! shifting z:
    do irec2 = 1, stars%ng2
      exp_m(1:vacuum%nmz,irec2) = exp_m(1:vacuum%nmz,irec2) * expDhg(irec2)
      exp_p(1:vacuum%nmz,irec2) = exp_p(1:vacuum%nmz,irec2) / expDhg(irec2)
    end do

    ! case irec2 > 1:
    do irec2 = 2, stars%ng2
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
    do ivac = 1, vacuum%nvac
      ga(1:vacuum%nmz) = rht(1:vacuum%nmz,ivac) * exp_m(1:vacuum%nmz,1)
      gb(1:vacuum%nmz) = rht(1:vacuum%nmz,ivac) * exp_p(1:vacuum%nmz,1)
      call intgz1Reverse( ga(:), vacuum%delz, vacuum%nmz, delta(:,ivac), .true. ) ! integrals 
      call qsf( vacuum%delz, gb(:), epsilon(:,ivac), vacuum%nmz, 1 )
      alphm(1,ivac) = delta(1,ivac)
    end do
    if ( vacuum%nvac == 1 ) alphm(1,2) = alphm(1,1)
    do ivac = 1, vacuum%nvac
      zeta(1:vacuum%nmz) = exp_m(1:vacuum%nmz,1) * ( alphm(1,mod(ivac,2)+1) + epsilon(1:vacuum%nmz,ivac) ) &
                         + exp_p(1:vacuum%nmz,1) *                              delta(1:vacuum%nmz,ivac) 
      where ( 2. * zeta == zeta ) zeta = 0.
      VVz(1:vacuum%nmz,ivac) = VVz(1:vacuum%nmz,ivac) + vcons(1) * zeta(1:vacuum%nmz)
    end do

    ! damping in vacuum:
    do ivac = 1, vacuum%nvac
      VVz(1:vacuum%nmz,ivac) = VVz(1:vacuum%nmz,ivac) * exp( -0.1 / rmt * z(1:vacuum%nmz) )
      do irec2 = 2, stars%ng2
        VVxy(1:vacuum%nmzxy,irec2,ivac) = VVxy(1:vacuum%nmzxy,irec2,ivac) * exp( -0.1 / rmt * z(1:vacuum%nmzxy) )
      end do
    end do


  end subroutine VYukawaFilmVacuumVariant1



  subroutine VYukawaFilmInterstitialVariant1( &
                stars, vacuum, cell, sym, input, &
                psq, VVxy, VVz, alphm, &
                VIq )

    ! main parts:
    ! 1. part: Compute the contribution from the interstitial charge density to the interstitial potential as a function of q_xy and z (analytic expression for integral)
    ! 2. part: Add the contribution from the vacuum charge density to the interstitial potential, which had already been computed earlier for the vacuum potential
 
    ! 4. part: Compute the coefficients V^I(q_xy,q_z) from the function V^I(q_xy,z) by a 1D Fourier transform:
    !          V^I(q_xy,z) = sum_{q_z} V^I(q_xy,q_z) * exp( ImagUnit * q_z * z )
    ! In order to be able to match the interstitial and vacuum potentials smoothly at the interface, the Fourier transform is done
    ! in a slightly larger region. -> 3. part
    ! 3. part: Interpolate the vacuum potential in a small region surrounding the slab

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

    real                        :: partitioning, rz, qz, q
    integer                     :: irec2, irec3, iz, jz, ivac, iqz, nfft, nzmax, nzmin, nzdh, nLower, nUpper, jvac
    complex, allocatable        :: VIz(:,:), eta(:,:)
    complex                     :: VIqz(-stars%mx3:stars%mx3,stars%ng2), c_ph(-stars%mx3:stars%mx3,stars%ng2)
    complex                     :: vcons1(stars%ng3)
    real                        :: VIzReal(3*stars%mx3,stars%ng2), VIzImag(3*stars%mx3,stars%ng2)
    real, allocatable           :: exp_m(:,:), exp_p(:,:)
    real, allocatable           :: z(:)
    real                        :: g_damped(stars%ng2), vcons2(stars%ng2), expDhg(stars%ng2)
    

    ! DEFINITIONS / ALLOCATIONS / INITIALISATIONS

    ! grid points z_i:
    nfft = 3 * stars%mx3                                  ! number of grid points for Fourier transform
    partitioning = 1. / real( nfft )
    nzmax = nfft / 2                                      ! index of maximal z below D~/2 
    nzmin = nzmax - nfft + 1                              ! index of minimal z above -D~/2
    nzdh = ceiling( cell%z1 / cell%amat(3,3) * nfft ) - 1 ! index of maximal z below D/2
    allocate( z(nzmin:nzmax) )
    ! definition of z_i:        ! indexing: z_0 = 0; positive indices for z > 0; negative indices for z < 0
    do iz = nzmin, nzmax
      z(iz) = cell%amat(3,3) * iz * partitioning
    end do
    ! other variables:
    allocate( VIz(nzmin:nzmax,stars%ng2), eta(nzmin:nzmax,stars%ng2) )
    allocate( exp_m(nzmin:nzmax,stars%ng2), exp_p(nzmin:nzmax,stars%ng2) )
    do irec2 = 1, stars%ng2
      g_damped(irec2) = sqrt( stars%sk2(irec2) ** 2 + input%preconditioning_param ** 2 )
      vcons2(irec2) = -1. / ( 2. * g_damped(irec2) )
      do iz = nzmin, nzmax
        exp_m(iz,irec2) = exp_save( - g_damped(irec2) * z(iz) )
        exp_p(iz,irec2) = exp_save(   g_damped(irec2) * z(iz) )
      end do
      expDhg(irec2) = exp_save( - cell%z1 * g_damped(irec2) )
    end do
    do irec3 = 1, stars%ng3
      vcons1(irec3) = fpi_const * psq(irec3) / ( stars%sk3(irec3) ** 2 + input%preconditioning_param ** 2 )
    end do
    VIz = (0.,0.)
    VIq = (0.,0.)


    ! CONTRIBUTION FROM THE INTERSTITIAL CHARGE DENSITY
    
    ! compute V^I(q_xy,z) as a function of q_xy and z
    do irec2 = 1, stars%ng2  
      do iz = -nzdh, nzdh
        do iqz = -stars%mx3, stars%mx3
          irec3 = stars%ig(stars%kv2(1,irec2),stars%kv2(2,irec2),iqz)
          if ( irec3 /= 0 ) then ! use only stars within the g_max sphere
            c_ph(iqz,irec2) = stars%rgphs(stars%kv2(1,irec2),stars%kv2(2,irec2),iqz)
            qz = iqz * cell%bmat(3,3)
            VIz(iz,irec2) = VIz(iz,irec2) &
                          + vcons1(irec3) * c_ph(iqz,irec2) * &
                            ( exp( ImagUnit * qz * z(iz) ) &
                            + vcons2(irec2) * expDhg(irec2) * &
                              ( ( g_damped(irec2) + ImagUnit * qz ) * exp_p(iz,irec2) * exp(   ImagUnit * qz * cell%z1 ) &
                              + ( g_damped(irec2) - ImagUnit * qz ) * exp_m(iz,irec2) * exp( - ImagUnit * qz * cell%z1 ) ) )
          end if
        enddo
      end do
    ! irec2 loop continues


    ! CONTRIBUTION FROM THE VACUUM CHARGE DENSITY

    ! irec2 loop continues
      eta(-nzdh:nzdh,irec2) = exp_m(-nzdh:nzdh,irec2) * alphm(irec2,2) + exp_p(-nzdh:nzdh,irec2) * alphm(irec2,1)
      where ( 2.0 * eta(:,irec2) == eta(:,irec2) ) eta(:,irec2) = cmplx( 0.0, 0.0 )
      VIz(-nzdh:nzdh,irec2) = VIz(-nzdh:nzdh,irec2) + tpi_const / g_damped(irec2) * eta(-nzdh:nzdh,irec2)
    ! irec2 loop continues
    
    
    ! INTERPOLATION IN VACUUM REGION

    ! use Lagrange polynomials of order 3 to interpolate the vacuum potential outside I
    ! q, q-1 and q-2 (scaled with +/-1 or +/-0.5) are the factors of the Lagrange basis polynomials
    ! irec2 loop continues
      do ivac = 1, 2
        select case( ivac )
          case( 1 )
            nUpper = nzmax; nLower =  nzdh + 1; jvac = 1
          case( 2 )
            nLower = nzmin; nUpper = -nzdh - 1; jvac = 2; if ( sym%invs .or. sym%zrfs ) jvac = 1
        end select
        do iz = nLower, nUpper
          rz = ( abs( z(iz) ) - cell%z1 ) / vacuum%delz + 1.0
          jz = rz      ! index of maximal vacuum grid point below z_i
          q = rz - jz  ! factor in Lagrange basis polynomials
          if ( irec2 == 1 ) then
            VIz(iz,irec2) = 0.5     * ( q - 1. ) * ( q - 2. ) * VVz(jz,  jvac) &
                          -       q *              ( q - 2. ) * VVz(jz+1,jvac) &
                          + 0.5 * q * ( q - 1. )              * VVz(jz+2,jvac)
          else if ( jz + 2 <= vacuum%nmzxy ) then
            VIz(iz,irec2) = 0.5 *     ( q - 1. ) * ( q - 2. ) * VVxy(jz,  irec2,jvac) &
                          -       q              * ( q - 2. ) * VVxy(jz+1,irec2,jvac) &
                          + 0.5 * q * ( q - 1. )              * VVxy(jz+2,irec2,jvac)
            if ( ( sym%invs .and. .not. sym%zrfs ) .and. ivac == 2 ) VIz(iz,irec2) = conjg( VIz(iz,irec2) )
          end if
        end do
      end do
    end do ! irec2


    ! 1D FOURIER TRANSFORM TO FIND THE COEFFICIENTS V^I(q_xy,q_z)

    ! change the indexing for the subroutine cfft, and split real and imaginary parts
    VIzReal(1:nzmax+1,:) =  real( VIz(0:nzmax,:) ); VIzReal(nzmax+2:nfft,:) =  real( VIz(nzmin:-1,:) )
    VIzImag(1:nzmax+1,:) = aimag( VIz(0:nzmax,:) ); VIzImag(nzmax+2:nfft,:) = aimag( VIz(nzmin:-1,:) )
 
    ! V^I(q_xy,z) = sum_{q_z} V^I(q_xy,q_z) * exp( ImagUnit * q_z * z )
    do irec2 = 1, stars%ng2
      call cfft( VIzReal(:,irec2), VIzImag(:,irec2), nfft, nfft, nfft, -1 )
    ! irec2 loop continues

    ! reorder
      VIqz(0,irec2) = cmplx( VIzReal(1,irec2), VIzImag(1,irec2) )
      do iqz = 1, stars%mx3
        VIqz( iqz,irec2) = cmplx( VIzReal(iqz+1,     irec2), VIzImag(iqz+1,     irec2) )
        VIqz(-iqz,irec2) = cmplx( VIzReal(nfft+1-iqz,irec2), VIzImag(nfft+1-iqz,irec2) )
      end do
   
    ! add the computed components to V^I(q_xy,q_z):
      do iqz= -stars%mx3, stars%mx3
        irec3 = stars%ig(stars%kv2(1,irec2),stars%kv2(2,irec2),iqz)
        if ( irec3 /= 0 ) VIq(irec3) = VIq(irec3) + VIqz(iqz,irec2) * partitioning / ( stars%nstr(irec3) / stars%nstr2(irec2) )
      end do
    end do


  end subroutine VYukawaFilmInterstitialVariant1



  subroutine VYukawaFilmVacuumVariant2( &
                stars, vacuum, cell, sym, input, rmt, &
                psq, rhtxy, rht, &
                VVxy, VVz, alphm, dh_prec, coshdh )

    ! 1. part: Compute the contribution from the interstitial charge density to the vacuum potential as a function of q_xy and z (analytic expression for integral)
    ! 2. part: Compute the contribution from the vacuum charge density to the vacuum potential as a function of q_xy and z by numerical integration

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
    real,           intent(in)  :: rmt
    complex,        intent(in)  :: psq(stars%ng3)
    complex,        intent(in)  :: rhtxy(vacuum%nmzxyd,stars%ng2-1,2)
    real,           intent(in)  :: rht(vacuum%nmzd,2) 

    complex,        intent(out) :: VVxy(vacuum%nmzxyd,2:stars%ng2,2) ! this is the qxy /= 0 part of the vacuum potential
    real,           intent(out) :: VVz(vacuum%nmzd,2)                ! this is the qxy  = 0 part of the vacuum potential
    complex,        intent(out) :: alphm(stars%ng2,2)                ! these are the integrals in upper and lower vacuum, now including the first star---integral for star ig2 is in alphm(ig2,ivac) 
    real,           intent(out) :: dh_prec
    real,           intent(out) :: coshdh(stars%ng2)

    integer                                             :: iz, irec2, irec3, ivac, iqz, nzdhprec, sign
    real                                                :: dh, qz, qxy_numerics 
    real, allocatable                                   :: z(:)
    complex, dimension(stars%ng2)                       :: sum_qz
    complex, dimension(stars%ng3)                       :: vcons3
    real, dimension(stars%ng2)                          :: g_damped, vcons2
    complex, dimension(-stars%mx3:stars%mx3,stars%ng2)  :: c_ph, quotq
    real, allocatable                                   :: quotz(:,:), sinhz(:,:)
    real, dimension(-1:1,stars%ng2)                     :: quotvardh
    complex, dimension(-stars%mx3:stars%mx3)            :: expdhqz
    complex, allocatable                                :: fa(:,:), fb(:,:)
    complex, allocatable                                :: alpha(:,:,:), beta(:,:,:), gamma(:,:)
    real, allocatable                                   :: ga(:), gb(:)
    real, allocatable                                   :: delta(:,:), epsilon(:,:), zeta(:)


    ! DEFINITIONS / ALLOCATIONS / INITIALISATIONS    
    
    dh = cell%z1 ! half the thickness of the film
    nzdhprec = ceiling( rmt / vacuum%delz ) ! index of dh_prec, see below
    dh_prec = dh + ( nzdhprec - 1 ) * vacuum%delz ! dh_prec is about dh + R; preconditioning boundary
    qxy_numerics = sqrt( ( 100 / dh_prec ) ** 2 - input%preconditioning_param ** 2 )
    allocate( z(nzdhprec) )
    allocate( quotz(-nzdhprec:nzdhprec,stars%ng2) )
    allocate( sinhz(nzdhprec,stars%ng2) )
    do iz = 1, nzdhprec ! new boundary
      z(iz) = ( iz - 1 ) * vacuum%delz + dh
    end do
    do irec2 = 1, stars%ng2
      g_damped(irec2) = sqrt( stars%sk2(irec2) ** 2 + input%preconditioning_param ** 2 )
      vcons2(irec2) = fpi_const / g_damped(irec2)
      coshdh(irec2) = cosh( g_damped(irec2) * ( dh_prec - dh ) )
      if( stars%sk2(irec2) < qxy_numerics ) then ! numerics ok
        do iz = 1, nzdhprec
          quotz( iz,irec2) = ( cosh( g_damped(irec2) * z(iz) ) / cosh( g_damped(irec2) * dh_prec ) &
                             + sinh( g_damped(irec2) * z(iz) ) / sinh( g_damped(irec2) * dh_prec ) ) / 2
          quotz(-iz,irec2) = ( cosh( g_damped(irec2) * z(iz) ) / cosh( g_damped(irec2) * dh_prec ) &
                             - sinh( g_damped(irec2) * z(iz) ) / sinh( g_damped(irec2) * dh_prec ) ) / 2
        end do
        quotvardh( 1,irec2) = ( cosh( g_damped(irec2) * dh ) / sinh( g_damped(irec2) * dh_prec ) &
                              + sinh( g_damped(irec2) * dh ) / cosh( g_damped(irec2) * dh_prec ) ) / 2
        quotvardh(-1,irec2) = ( cosh( g_damped(irec2) * dh ) / sinh( g_damped(irec2) * dh_prec ) &
                              - sinh( g_damped(irec2) * dh ) / cosh( g_damped(irec2) * dh_prec ) ) / 2
      else ! numerical treatment necessary
        do iz = 1, nzdhprec
          quotz( iz,irec2) = exp_save( g_damped(irec2) * (  z(iz) - dh_prec ) )
          quotz(-iz,irec2) = exp_save( g_damped(irec2) * ( -z(iz) - dh_prec ) )
        end do
        quotvardh( 1,irec2) = quotz( 1,irec2)
        quotvardh(-1,irec2) = quotz(-1,irec2)
      end if
      do iz = 1, nzdhprec
        sinhz(iz,irec2) = sinh( g_damped(irec2) * ( dh_prec - z(iz) ) )
      end do
      do iqz = -stars%mx3, stars%mx3
        quotq(iqz,irec2) = ImagUnit * iqz * cell%bmat(3,3) / g_damped(irec2)
      end do
    end do
    sum_qz = (0.,0.)
    VVxy = (0.,0.)
    VVz = 0.
    do irec3 = 1, stars%ng3
      vcons3(irec3) = fpi_const * psq(irec3) / ( stars%sk3(irec3) ** 2 + input%preconditioning_param ** 2 )
    end do
    do iqz = -stars%mx3, stars%mx3
      expdhqz(iqz) = exp( ImagUnit * iqz * cell%bmat(3,3) * dh )
    end do


    ! CONTRIBUTION FROM THE INTERSTITIAL CHARGE DENSITY

    do ivac = 1, vacuum%nvac
      sign = 3 - 2 * ivac
      do irec2 = 1, stars%ng2
        do iqz = -stars%mx3, stars%mx3
          irec3 = stars%ig(stars%kv2(1,irec2),stars%kv2(2,irec2),iqz)
          ! use only stars within the g_max sphere -> stars outside the sphere have per definition index ig3n = 0
          if ( irec3 /= 0 ) then
            c_ph(iqz,irec2) = stars%rgphs(stars%kv2(1,irec2),stars%kv2(2,irec2),iqz)
            sum_qz(irec2) = sum_qz(irec2) + &
                            c_ph(iqz,irec2) * vcons3(irec3) * &
                            ( ( quotvardh( 1,irec2) - sign * quotq(iqz,irec2) * quotz( 1,irec2) ) * expdhqz( sign*iqz) &
                            - ( quotvardh(-1,irec2) - sign * quotq(iqz,irec2) * quotz(-1,irec2) ) * expdhqz(-sign*iqz) )
          endif
        enddo
        if( irec2 /= 1 ) then
          VVxy(1:nzdhprec,irec2,ivac) = sum_qz(irec2) * sign * sinhz(1:nzdhprec,irec2)
        else
          VVz( 1:nzdhprec,      ivac) = sum_qz(1)     * sign * sinhz(1:nzdhprec,1)
        end if
      enddo
    enddo


    ! CONTRIBUTION FROM THE VACUUM CHARGE DENSITY

    allocate( fa(nzdhprec,2:stars%ng2), fb(nzdhprec,2:stars%ng2) )
    allocate( alpha(nzdhprec,2:stars%ng2,2), beta(nzdhprec,2:stars%ng2,2), gamma(nzdhprec,2:stars%ng2) )
    ! case irec2 > 1:
    do irec2 = 2, stars%ng2
      do ivac = 1, vacuum%nvac
        ! integrands:
        fa(1:nzdhprec,irec2) = rhtxy(1:nzdhprec,irec2-1,ivac) * sinhz(1:nzdhprec,irec2)
        fb(1:nzdhprec,irec2) = rhtxy(1:nzdhprec,irec2-1,ivac) * quotz(1:nzdhprec,irec2)
        ! integrals:
        ! alpha(z,q_xy,ivac) = int_z^infty rho(z',q_xy,ivac) sinhz dz'
        ! beta (z,q_xy,ivac) = int_{D/2}^z rho(z',q_xy,ivac) quotz dz'
        ! where for z < 0 the lower vacuum charge density (ivac=2) is defined by rho(q_xy,z,ivac=2) := rho(q_xy,-z,ivac=2)
        call intgz1Reverse( fa(:,irec2), vacuum%delz, nzdhprec, alpha(:,irec2,ivac), .false. ) 
        call qsfComplex( vacuum%delz, fb(:,irec2), beta(:,irec2,ivac), nzdhprec, 1 )
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
        gamma(1:nzdhprec,irec2) = quotz(-1:-nzdhprec:-1,irec2) * alphm(irec2,mod(ivac,2)+1) &
                                + quotz(1:nzdhprec,irec2) * alpha(1:nzdhprec,irec2,ivac) &
                                + sinhz(1:nzdhprec,irec2) * beta(1:nzdhprec,irec2,ivac)
        VVxy(1:nzdhprec,irec2,ivac) = VVxy(1:nzdhprec,irec2,ivac) + vcons2(irec2) * gamma(1:nzdhprec,irec2) 
      end do
    end do


    allocate( ga(nzdhprec), gb(nzdhprec) )
    allocate( delta(nzdhprec,2), epsilon(nzdhprec,2), zeta(nzdhprec) )
    ! case irec2 = 1:
    do ivac = 1, vacuum%nvac
      ga(1:nzdhprec) = rht(1:nzdhprec,ivac) * sinhz(1:nzdhprec,1)
      gb(1:nzdhprec) = rht(1:nzdhprec,ivac) * quotz(1:nzdhprec,1)
      call intgz1Reverse( ga(:), vacuum%delz, nzdhprec, delta(:,ivac), .false. ) 
      call qsf( vacuum%delz, gb(:), epsilon(:,ivac), nzdhprec, 1 )
      alphm(1,ivac) = delta(1,ivac)
    end do
    if ( vacuum%nvac == 1 ) alphm(1,2) = alphm(1,1)
    do ivac = 1, vacuum%nvac
      zeta(1:nzdhprec) = quotz(-1:-nzdhprec:-1,1) * alphm(1,mod(ivac,2)+1) &
                       + quotz(1:nzdhprec,1) * delta(1:nzdhprec,ivac) &
                       + sinhz(1:nzdhprec,1) * epsilon(1:nzdhprec,ivac) 
      VVz(1:nzdhprec,ivac) = VVz(1:nzdhprec,ivac) + vcons2(1) * zeta(1:nzdhprec)
    end do


  end subroutine VYukawaFilmVacuumVariant2



  subroutine VYukawaFilmInterstitialVariant2( &
                stars, vacuum, cell, sym, input, &
                psq, VVxy, VVz, alphm, dh_prec, coshdh, &
                VIq )

    ! main parts:
    ! 1. part: Compute the contribution from the interstitial charge density to the interstitial potential as a function of q_xy and z (analytic expression for integral)
    ! 2. part: Add the contribution from the vacuum charge density to the interstitial potential, which had already been computed earlier for the vacuum potential
 
    ! 4. part: Compute the coefficients V^I(q_xy,q_z) from the function V^I(q_xy,z) by a 1D Fourier transform:
    !          V^I(q_xy,z) = sum_{q_z} V^I(q_xy,q_z) * exp( ImagUnit * q_z * z )
    ! In order to be able to match the interstitial and vacuum potentials smoothly at the interface, the Fourier transform is done
    ! in a slightly larger region. -> 3. part
    ! 3. part: Interpolate the vacuum potential in a small region surrounding the slab

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
    real,           intent(in)  :: dh_prec
    real,           intent(in)  :: coshdh(stars%ng2)

    complex,        intent(out) :: VIq(stars%ng3)

    real                                               :: partitioning, rz, qz, q, qxy_numerics
    integer                                            :: irec2, irec3, iz, jz, ivac, iqz, jvac
    integer                                            :: nfft, nzmax, nzmin, nzdh, nLower, nUpper
    complex, allocatable                               :: VIz(:,:), eta(:,:)
    complex, allocatable                               :: expzqz(:,:)
    complex, dimension(-stars%mx3:stars%mx3,stars%ng2) :: VIqz, c_ph
    complex, dimension(-stars%mx3:stars%mx3,stars%ng2) :: qquot, qquottrigp, qquottrigm
    complex, dimension(stars%ng3)                      :: vcons3
    real, dimension(3*stars%mx3,stars%ng2)             :: VIzReal, VIzImag
    real, allocatable                                  :: quotz(:,:), sinhz(:,:)
    real, allocatable                                  :: z(:)
    real, dimension(stars%ng2)                         :: g_damped, vcons2
    

    ! DEFINITIONS / ALLOCATIONS / INITIALISATIONS

    ! grid points z_i:
    qxy_numerics = sqrt( ( 100 / dh_prec ) ** 2 - input%preconditioning_param ** 2 )
    nfft = 3 * stars%mx3                                  ! number of grid points for Fourier transform
    partitioning = 1. / real( nfft )
    nzmax = nfft / 2                                      ! index of maximal z below D~/2 
    nzmin = nzmax - nfft + 1                              ! index of minimal z above -D~/2
    nzdh = ceiling( cell%z1 / cell%amat(3,3) * nfft ) - 1 ! index of maximal z below D/2
    allocate( z(nzmin:nzmax) )
    ! definition of z_i:        ! indexing: z_0 = 0; positive indices for z > 0; negative indices for z < 0
    do iz = nzmin, nzmax
      z(iz) = cell%amat(3,3) * iz * partitioning
    end do
    ! other variables:
    allocate( VIz(nzmin:nzmax,stars%ng2) )
    allocate( eta(-nzdh:nzdh,stars%ng2) )
    allocate( sinhz(-nzdh:nzdh,stars%ng2) )
    allocate( quotz(-nzdh:nzdh,stars%ng2) )
    allocate( expzqz(-nzdh:nzdh,-stars%mx3:stars%mx3) )
    do irec2 = 1, stars%ng2
      g_damped(irec2) = sqrt( stars%sk2(irec2) ** 2 + input%preconditioning_param ** 2 )
      vcons2(irec2) = fpi_const / g_damped(irec2)
      if( stars%sk2(irec2) < qxy_numerics ) then ! numerics ok
        do iz = -nzdh, nzdh
          quotz( iz,irec2) = ( cosh( g_damped(irec2) * z(iz) ) / cosh( g_damped(irec2) * dh_prec ) &
                             + sinh( g_damped(irec2) * z(iz) ) / sinh( g_damped(irec2) * dh_prec ) ) / 2
        end do
      else ! numerical treatment necessary
        do iz = -nzdh, nzdh
          quotz( iz,irec2) = exp_save( g_damped(irec2) * (  z(iz) - dh_prec ) )
        end do
      end if
      do iz = -nzdh, nzdh
        sinhz(iz,irec2) = sinh( g_damped(irec2) * ( dh_prec - z(iz) ) )
        do iqz = -stars%mx3, stars%mx3
          expzqz(iz,iqz) = exp( ImagUnit * iqz * cell%bmat(3,3) * z(iz) )
        end do
      end do
      do iqz = -stars%mx3, stars%mx3
        qquot(iqz,irec2) = ImagUnit * iqz * cell%bmat(3,3) / g_damped(irec2)
        qquottrigp(iqz,irec2) = coshdh(irec2) + sinhz(nzdh,irec2) * qquot(iqz,irec2)
        qquottrigm(iqz,irec2) = coshdh(irec2) - sinhz(nzdh,irec2) * qquot(iqz,irec2)
      end do
    end do
    do irec3 = 1, stars%ng3
      vcons3(irec3) = fpi_const * psq(irec3) / ( stars%sk3(irec3) ** 2 + input%preconditioning_param ** 2 )
    end do
    VIz = (0.,0.)
    VIq = (0.,0.)


    ! CONTRIBUTION FROM THE INTERSTITIAL CHARGE DENSITY
    
    ! compute V^I(q_xy,z) as a function of q_xy and z
    do irec2 = 1, stars%ng2  
      do iz = -nzdh, nzdh
        do iqz = -stars%mx3, stars%mx3
          irec3 = stars%ig(stars%kv2(1,irec2),stars%kv2(2,irec2),iqz)
          if ( irec3 /= 0 ) then ! use only stars within the g_max sphere
            c_ph(iqz,irec2) = stars%rgphs(stars%kv2(1,irec2),stars%kv2(2,irec2),iqz)
            qz = iqz * cell%bmat(3,3)
            VIz(iz,irec2) = VIz(iz,irec2) &
                          + vcons3(irec3) * c_ph(iqz,irec2) * &
                            ( expzqz(iz,iqz) &
                            - qquottrigp(iqz,irec2) * expzqz( nzdh,iqz) * quotz( iz,irec2) &
                            - qquottrigm(iqz,irec2) * expzqz(-nzdh,iqz) * quotz(-iz,irec2) )
          end if
        enddo
      end do
    ! irec2 loop continues


    ! CONTRIBUTION FROM THE VACUUM CHARGE DENSITY

    ! irec2 loop continues
      eta(-nzdh:nzdh,irec2) = quotz(nzdh:-nzdh:-1,irec2) * alphm(irec2,2) &
                            + quotz(-nzdh:nzdh,   irec2) * alphm(irec2,1)
      VIz(-nzdh:nzdh,irec2) = VIz(-nzdh:nzdh,irec2) + vcons2(irec2) * eta(-nzdh:nzdh,irec2)
    ! irec2 loop continues
   
    
    ! INTERPOLATION IN VACUUM REGION

    ! use Lagrange polynomials of order 3 to interpolate the vacuum potential outside I
    ! q, q-1 and q-2 (scaled with +/-1 or +/-0.5) are the factors of the Lagrange basis polynomials
    ! irec2 loop continues
      do ivac = 1, 2
        select case( ivac )
          case( 1 )
            nUpper = nzmax; nLower =  nzdh + 1; jvac = 1
          case( 2 )
            nLower = nzmin; nUpper = -nzdh - 1; jvac = 2; if ( sym%invs .or. sym%zrfs ) jvac = 1
        end select
        do iz = nLower, nUpper
          rz = ( abs( z(iz) ) - cell%z1 ) / vacuum%delz + 1.0
          jz = rz      ! index of maximal vacuum grid point below z_i
          q = rz - jz  ! factor in Lagrange basis polynomials
          if ( irec2 == 1 ) then
            VIz(iz,irec2) = 0.5     * ( q - 1. ) * ( q - 2. ) * VVz(jz,  jvac) &
                          -       q *              ( q - 2. ) * VVz(jz+1,jvac) &
                          + 0.5 * q * ( q - 1. )              * VVz(jz+2,jvac)
          else if ( jz + 2 <= vacuum%nmzxy ) then
            VIz(iz,irec2) = 0.5 *     ( q - 1. ) * ( q - 2. ) * VVxy(jz,  irec2,jvac) &
                          -       q              * ( q - 2. ) * VVxy(jz+1,irec2,jvac) &
                          + 0.5 * q * ( q - 1. )              * VVxy(jz+2,irec2,jvac)
            if ( ( sym%invs .and. .not. sym%zrfs ) .and. ivac == 2 ) VIz(iz,irec2) = conjg( VIz(iz,irec2) )
          end if
        end do
      end do
    end do ! irec2


    ! 1D FOURIER TRANSFORM TO FIND THE COEFFICIENTS V^I(q_xy,q_z)

    ! change the indexing for the subroutine cfft, and split real and imaginary parts
    VIzReal(1:nzmax+1,:) =  real( VIz(0:nzmax,:) ); VIzReal(nzmax+2:nfft,:) =  real( VIz(nzmin:-1,:) )
    VIzImag(1:nzmax+1,:) = aimag( VIz(0:nzmax,:) ); VIzImag(nzmax+2:nfft,:) = aimag( VIz(nzmin:-1,:) )
 
    ! V^I(q_xy,z) = sum_{q_z} V^I(q_xy,q_z) * exp( ImagUnit * q_z * z )
    do irec2 = 1, stars%ng2
      call cfft( VIzReal(:,irec2), VIzImag(:,irec2), nfft, nfft, nfft, -1 )
    ! irec2 loop continues

    ! reorder
      VIqz(0,irec2) = cmplx( VIzReal(1,irec2), VIzImag(1,irec2) )
      do iqz = 1, stars%mx3
        VIqz( iqz,irec2) = cmplx( VIzReal(iqz+1,     irec2), VIzImag(iqz+1,     irec2) )
        VIqz(-iqz,irec2) = cmplx( VIzReal(nfft+1-iqz,irec2), VIzImag(nfft+1-iqz,irec2) )
      end do
   
    ! add the computed components to V^I(q_xy,q_z):
      do iqz= -stars%mx3, stars%mx3
        irec3 = stars%ig(stars%kv2(1,irec2),stars%kv2(2,irec2),iqz)
        if ( irec3 /= 0 ) VIq(irec3) = VIq(irec3) + VIqz(iqz,irec2) * partitioning / ( stars%nstr(irec3) / stars%nstr2(irec2) )
      end do
    end do


  end subroutine VYukawaFilmInterstitialVariant2



  subroutine VYukawaModify( stars, vacuum, cell, sym, input, mpi, atoms, sphhar, oneD, noco, den, &
                            VYukawa )

    ! This subroutine adds a potential to the previously computed Yukawa
    ! potential to ensure charge neutrality.
    ! The added potential itself is a solution to the modified Helmholtz
    ! equation with constant right-hand side, where the constant is chosen such
    ! that charge neutrality is obtained.
    ! The charge is distributed only over the film region, and we therefore 
    ! solve the differential equation subject to a boundary condition on the 
    ! film surface, in contrast to the basic Yukawa potential above. 

    use m_constants
    use m_types
    use m_vmts
    use m_constants
    use m_cdntot
    use m_cfft
    implicit none

    type(t_stars),      intent(in)    :: stars
    type(t_vacuum),     intent(in)    :: vacuum
    type(t_cell),       intent(in)    :: cell
    type(t_sym),        intent(in)    :: sym
    type(t_input),      intent(in)    :: input
    type(t_mpi),        intent(in)    :: mpi
    type(t_atoms),      intent(in)    :: atoms
    type(t_sphhar),     intent(in)    :: sphhar
    type(t_oneD),       intent(in)    :: oneD
    type(t_noco),       intent(in)    :: noco
    type(t_potden),     intent(inout) :: den
    type(t_potden),     intent(inout) :: VYukawa

    integer                           :: n, lh, irec3, iz, iqz, nfft, nzmax, nzmin, nzdh
    real                              :: q0, qhat, qbar, ldh, partitioning, dh
    real                              :: q(input%jspins), qis(input%jspins), qmt(atoms%ntype,input%jspins),qvac(2,input%jspins), qtot, qistot
    complex                           :: psq(stars%ng3)
    type(t_potden)                    :: VYukawaModification
    real, allocatable                 :: z(:)
    real                              :: VIzReal(3*stars%mx3), VIzImag(3*stars%mx3)
    complex, allocatable              :: VIz(:)
    complex                           :: VIqz(-stars%mx3:stars%mx3)


    ! DEFINITIONS / ALLOCATIONS / INITIALISATIONS

    ! constants:
    dh = cell%z1   ! half the width of the film
    ! indexing of grid points z_i:
    nfft = 3 * stars%mx3   ! number of grid points for Fourier transform
    partitioning = 1. / real( nfft )
    nzmax = nfft / 2           ! index of maximal z below D~/2 
    nzmin = nzmax - nfft + 1   ! index of minimal z above -D~/2
    nzdh = ceiling( dh / cell%amat(3,3) * nfft ) - 1   ! index of maximal z below D/2
    allocate( z(nzmin:nzmax) )
    ! definition of grid points z_i:
    ! indexing: z_0 = 0; positive indices for z > 0; negative indices for z < 0
    do iz = nzmin, nzmax
      z(iz) = cell%amat(3,3) * iz * partitioning
    end do


    ! INTEGRATION OF THE PREVIOUSLY COMPUTED YUKAWA POTENTIAL

    ! initialise VYukawaModification with in-going VYukawa and prepare for integration
    call VYukawaModification%init( stars, atoms, sphhar, vacuum, noco, input%jspins, 4 )
    call VYukawaModification%copyPotDen( VYukawa )
    do n = 1, atoms%ntype
      do lh = 0, sphhar%nlhd    
        VYukawaModification%mt(1:atoms%jri(n),lh,n,1) = VYukawaModification%mt(1:atoms%jri(n),lh,n,1) * atoms%rmsh(1:atoms%jri(n),n) ** 2
      end do
    end do

    ! integrate the potential over the film region
    call integrate_cdn( stars, atoms, sym, vacuum, input, cell, oneD, VYukawaModification, q, qis, qmt, qvac, qtot, qistot  )
    q0 = qtot / cell%area
    ldh = input%preconditioning_param * dh
    qhat = ( q0 / ( 2 * dh ) ) / ( sinh(ldh) / ( ldh * cosh( ldh ) ) - 1 )
    qbar = input%preconditioning_param ** 2 / fpi_const *  qhat


    ! SET UP CONSTANT CHARGE DENSITY

    ! instead of den%pw(1,1) = qbar we directly set the pseudo charge density
    den%mt = 0; den%pw = 0; den%vacxy = 0; den%vacz = 0    
    do n = 1, atoms%ntype
      den%mt(1:atoms%jri(n),0,n,1) = sfp_const * qbar * atoms%rmsh(1:atoms%jri(n),n) ** 2
    end do
    psq = cmplx(0.0,0.0); psq(1) = qbar


    ! CALCULATE THE INTERSTITIAL POTENTIAL AS A FUNCTION OF z

    ! initialise and calculate out-going modification potential; reuse VYukawaModification
    VYukawaModification%mt = 0; VYukawaModification%pw = 0; VYukawaModification%vacxy = 0; VYukawaModification%vacz = 0

    allocate( VIz(nzmin:nzmax) )
    VIz = (0.,0.)
    do iz = -nzdh+1, nzdh-1
      VIz(iz) = qhat * ( 1 - cosh( input%preconditioning_param * z(iz) ) / cosh( ldh ) ) 
    end do


    ! 1D FOURIER TRANSFORM TO FIND THE 3D-FOURIER COEFFICIENTS

    ! change the indexing for the subroutine cfft, and split real and imaginary parts
    VIzReal(1:nzmax+1) =  real( VIz(0:nzmax) ); VIzReal(nzmax+2:nfft) = real( VIz(nzmin:-1) )
    VIzImag(1:nzmax+1) = aimag( VIz(0:nzmax) ); VIzImag(nzmax+2:nfft) = aimag( VIz(nzmin:-1) )

    call cfft( VIzReal(:), VIzImag(:), nfft, nfft, nfft, -1 )

    ! reorder
    VIqz = 0
    VIqz(0) = cmplx( VIzReal(1), VIzImag(1) )
    do iqz = 1, stars%mx3
      VIqz( iqz) = cmplx( VIzReal(iqz+1     ), VIzImag(iqz+1     ) )
      VIqz(-iqz) = cmplx( VIzReal(nfft+1-iqz), VIzImag(nfft+1-iqz) )
    end do

    ! add the computed components
    do iqz= -stars%mx3, stars%mx3
      irec3 = stars%ig(stars%kv2(1,1),stars%kv2(2,1),iqz)
      if ( irec3 /= 0 ) VYukawaModification%pw(irec3,1) = VYukawaModification%pw(irec3,1) + VIqz(iqz) * partitioning / ( stars%nstr(irec3) / stars%nstr2(1) )
    end do


    ! MUFFIN-TIN POTENTIAL

    call Vmts( input, mpi, stars, sphhar, atoms, sym, cell, oneD, .FALSE., &
               VYukawaModification%pw(:,1), den%mt(:,0:,:,1), VYukawaModification%potdenType, &
               VYukawaModification%mt(:,0:,:,1) )


    ! APPLYING THE MODIFICATION TO THE YUKAWA POTENTIAL

    call VYukawa%AddPotDen( VYukawa, VYukawaModification )


  end subroutine VYukawaModify



end module m_VYukawaFilm
