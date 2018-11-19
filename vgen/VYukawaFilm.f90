module m_VYukawaFilm

  ! Computation of the film-case Yukawa potential for the preconditioning of the charge density residual

  ! 1. part: The pseudo-charge density is used for the vacuum and interstitial potentials.
  ! 4. part: The MT potential is calculated at the end.
  ! These two parts do not change compared to the bulk case.

  ! 2. and 3. part: The vacuum and interstitial potentials are calculated by integration of the product of the charge density and
  ! the 1D Green function with respect to z (from  minus infinity to infinity). The integrals for both potentials are split in 4:
  ! the seperating points are the boundaries of the slab and the point where the integration variable z' equals the point z where we
  ! want to evaluate the potential for. For the contribution from the vacuum charge density we do a numerical integration. For the
  ! contribution from the slab charge density we have analytical expressions of the integrals.

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
    complex,        intent(in)  :: psq(stars%ng3)
    complex,        intent(in)  :: rhtxy(vacuum%nmzxyd,stars%ng2-1,2)
    real,           intent(in)  :: rht(vacuum%nmz,2) 

    complex,        intent(out) :: VVxy(vacuum%nmzxyd,2:stars%ng2,2) ! this is the qxy /= 0 part of the vacuum potential
    real,           intent(out) :: VVz(vacuum%nmz,2)                ! this is the qxy  = 0 part of the vacuum potential
    complex,        intent(out) :: alphm(stars%ng2,2)                ! these are the integrals in upper and lower vacuum, now including the first star---integral for star ig2 is in alphm(ig2,ivac) 

    complex                     :: sum_qz(2,stars%ng2)
    complex                     :: c_ph(-stars%mx3:stars%mx3,stars%ng2)
    complex                     :: signIqz
    real                        :: g_damped(stars%ng2), qz, sign, vcons(stars%ng2) 
    real                        :: exp_m(vacuum%nmz,stars%ng2), exp_p(vacuum%nmz,stars%ng2)
    real                        :: expDhg(stars%ng2), expDg(stars%ng2)
    real                        :: z(vacuum%nmz)
    integer                     :: iz, irec2, irec3, ivac, iqz
    complex                     :: fa(vacuum%nmzxyd,2:stars%ng2), fb(vacuum%nmzxyd,2:stars%ng2)
    complex                     :: alpha(vacuum%nmzxyd,2:stars%ng2,2), beta(vacuum%nmzxyd,2:stars%ng2,2), gamma(vacuum%nmzxyd,2:stars%ng2)
    real                        :: ga(vacuum%nmz), gb(vacuum%nmz)
    real                        :: delta(vacuum%nmz,2), epsilon(vacuum%nmz,2), zeta(vacuum%nmz)


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

    ! case irec2 > 1:
    do irec2 = 2, stars%ng2
      do ivac = 1, vacuum%nvac
        ! integrands:
        fa(1:vacuum%nmzxy,irec2) = rhtxy(1:vacuum%nmzxy,irec2-1,ivac) * exp_m(1:vacuum%nmzxy,irec2) * expDhg(irec2)
        fb(1:vacuum%nmzxy,irec2) = rhtxy(1:vacuum%nmzxy,irec2-1,ivac) * exp_p(1:vacuum%nmzxy,irec2) * expDhg(irec2)
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
        VVxy(1:vacuum%nmzxy,irec2,ivac) = VVxy(1:vacuum%nmzxy,irec2,ivac) + vcons(irec2) * gamma(1:vacuum%nmzxy,irec2) * expDhg(irec2) 
      end do
    end do

    ! case irec2 = 1:
    do ivac = 1, vacuum%nvac
      ga(1:vacuum%nmz) = rht(1:vacuum%nmz,ivac) * exp_m(1:vacuum%nmz,1) * expDhg(1)
      gb(1:vacuum%nmz) = rht(1:vacuum%nmz,ivac) * exp_p(1:vacuum%nmz,1) * expDhg(1)
      call intgz1Reverse( ga(:), vacuum%delz, vacuum%nmz, delta(:,ivac), .true. ) ! integrals 
      call qsf( vacuum%delz, gb(:), epsilon(:,ivac), vacuum%nmz, 1 )
      alphm(1,ivac) = delta(1,ivac)
    end do
    if ( vacuum%nvac == 1 ) alphm(1,2) = alphm(1,1)
    do ivac = 1, vacuum%nvac
      zeta(1:vacuum%nmz) = exp_m(1:vacuum%nmz,1) * ( alphm(1,mod(ivac,2)+1) + epsilon(1:vacuum%nmz,ivac) ) &
                         + exp_p(1:vacuum%nmz,1) *                              delta(1:vacuum%nmz,ivac) 
      where ( 2. * zeta == zeta ) zeta = 0.
      VVz(1:vacuum%nmz,ivac) = VVz(1:vacuum%nmz,ivac) + vcons(1) * zeta(1:vacuum%nmz) * expDhg(1)
    end do


  end subroutine VYukawaFilmVacuum



  subroutine VYukawaFilmInterstitial( stars, vacuum, cell, sym, input, &
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
    real,           intent(in)  :: VVz(vacuum%nmz,2)
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
            nUpper = nzmin; nLower = -nzdh - 1; jvac = 2; if ( sym%invs .or. sym%zrfs ) jvac = 1
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

  end subroutine VYukawaFilmInterstitial



end module m_VYukawaFilm
