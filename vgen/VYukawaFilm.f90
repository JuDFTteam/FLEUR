module m_VYukawaFilm

    ! In a first step we solve the modified Helmholtz equation with right-hand
    ! side -4 pi ( rho_out - rho_in ) subject to the boundary condition that the
    ! potential equals zero on the surface of the film.
    ! In the interstitial we use a Green function method to solve the
    ! differential equation (separately for each q_xy) for a z-dependent function.
    ! The muffin-tin potential generation is the same as for bulk materials.

    ! Unfortunately, the potential cannot satisfy the condition that its
    ! integral over the film is zero at the same time, which is needed for
    ! charge neutrality of the preconditioned density.
    ! Therefore we add a solution to the differential equation with constant
    ! right-hand side sufficing the same boundary condition, where the constant
    ! is chosen such that the integral over the sum of these potentials is zero.
    ! So we basically added a constant to the charge density in order to make
    ! the resulting Yukawa potential charge neutral.

contains

  subroutine VYukawaFilm( stars, vacuum, cell, sym, input, mpi, atoms, sphhar, dimension, oneD, noco, den, &
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
    type(t_noco),       intent(in)    :: noco
    type(t_potden),     intent(inout) :: den

    type(t_potden),     intent(inout) :: VYukawa

    complex                           :: psq(stars%ng3)


    ! PSEUDO-CHARGE DENSITY

    call psqpw( mpi, atoms, sphhar, stars, vacuum, dimension, cell, input, sym, oneD, & 
                den%pw(:,1), den%mt(:,:,:,1), den%vacz(:,:,1), .false., VYukawa%potdenType, &
                psq )


    ! INTERSTITIAL POTENTIAL

    call VYukawaFilmInterstitial( stars, vacuum, cell, sym, input, &
                                  psq, &
                                  VYukawa%pw(:,1) )


    ! MUFFIN-TIN POTENTIAL

    call Vmts( input, mpi, stars, sphhar, atoms, sym, cell, oneD, &
               VYukawa%pw(:,1), den%mt(:,0:,:,1), VYukawa%potdenType, &
               VYukawa%mt(:,0:,:,1) )

 
    ! MODIFICATION FOR CHARGE NEUTRALITY

    call VYukawaModify( stars, vacuum, cell, sym, input, mpi, atoms, sphhar, dimension, oneD, noco, den, &
                        VYukawa )


  end subroutine VYukawaFilm



  subroutine VYukawaFilmInterstitial( stars, vacuum, cell, sym, input, &
                                      psq, &
                                      VIq )

    ! This subroutine computes the basic interstitial Yukawa potential that is
    ! zero on the film surface. Its integral is not zero.

    use m_ExpSave
    use m_constants
    use m_types
    use m_cfft
    implicit none

    type(t_stars),  intent(in)    :: stars
    type(t_vacuum), intent(in)    :: vacuum
    type(t_cell),   intent(in)    :: cell
    type(t_sym),    intent(in)    :: sym
    type(t_input),  intent(in)    :: input
    complex,        intent(inout) :: psq(stars%ng3)

    complex,        intent(out)   :: VIq(stars%ng3)

    real                          :: partitioning, rz, qz, q, q0, dh, qxy_numerics
    integer                       :: irec2, irec3, iz, ivac, iqz, nfft, nzmax, nzmin, nzdh
    complex, allocatable          :: VIz(:,:)
    complex                       :: VIqz(-stars%mx3:stars%mx3,stars%ng2), c_ph(-stars%mx3:stars%mx3,stars%ng2)
    complex                       :: vcons1(stars%ng3)
    real                          :: VIzReal(3*stars%mx3,stars%ng2), VIzImag(3*stars%mx3,stars%ng2)
    real, allocatable             :: cosh_quot(:,:), sinh_quot(:,:), exp_m(:,:), exp_p(:,:)
    real, allocatable             :: z(:)
    real                          :: g_damped(stars%ng2), cosh_dh(stars%ng2), sinh_dh(stars%ng2)
    

    ! DEFINITIONS / ALLOCATIONS / INITIALISATIONS

    ! constants:
    dh = cell%z1
    qxy_numerics = sqrt( ( 100 / dh ) ** 2 - input%preconditioning_param ** 2 )
    ! indexing of the grid points z_i:
    nfft = 3 * stars%mx3   ! number of grid points for Fourier transform
    partitioning = 1. / real( nfft )
    nzmax = nfft / 2           ! index of maximal z below D~/2 
    nzmin = nzmax - nfft + 1   ! index of minimal z above -D~/2
    nzdh = ceiling( dh / cell%amat(3,3) * nfft ) - 1 ! index of maximal z below D/2
    allocate( z(nzmin:nzmax) )
    ! definition of the grid points z_i:
    ! indexing: z_0 = 0; positive indices for z > 0; negative indices for z < 0
    do iz = nzmin, nzmax
      z(iz) = cell%amat(3,3) * iz * partitioning
    end do
    ! other variables:
    allocate( VIz(nzmin:nzmax,stars%ng2) )
    allocate( cosh_quot(-nzdh+1:nzdh-1,stars%ng2), sinh_quot(-nzdh+1:nzdh-1,stars%ng2) )
    allocate( exp_m(-nzdh+1:nzdh-1,stars%ng2), exp_p(-nzdh+1:nzdh-1,stars%ng2) )    
    do irec2 = 1, stars%ng2
      g_damped(irec2) = sqrt( stars%sk2(irec2) ** 2 + input%preconditioning_param ** 2 )
      if( stars%sk2(irec2) < qxy_numerics ) then ! numerics ok
        cosh_dh(irec2) = cosh( g_damped(irec2) * dh )
        sinh_dh(irec2) = sinh( g_damped(irec2) * dh )
        do iz = -nzdh+1, nzdh-1
          cosh_quot(iz,irec2) = cosh( g_damped(irec2) * z(iz) ) / cosh_dh(irec2)
          sinh_quot(iz,irec2) = sinh( g_damped(irec2) * z(iz) ) / sinh_dh(irec2)
        end do
      else ! numerical treatment necessary
        do iz = -nzdh+1, nzdh-1
          exp_m(iz,irec2) = exp_save( - g_damped(irec2) * ( dh - z(iz) ) )
          exp_p(iz,irec2) = exp_save( - g_damped(irec2) * ( dh + z(iz) ) )
        end do
      end if
    end do
    do irec3 = 1, stars%ng3
      vcons1(irec3) = fpi_const * psq(irec3) / ( stars%sk3(irec3) ** 2 + input%preconditioning_param ** 2 )
    end do
    VIz = (0.,0.)
    VIq = (0.,0.)


    ! CONTRIBUTION FROM THE INTERSTITIAL CHARGE DENSITY
    ! compute V^I(q_xy,z) as a function of q_xy and z
    do irec2 = 1, stars%ng2  
      if( stars%sk2(irec2) < qxy_numerics ) then ! numerics ok
        do iz = -nzdh+1, nzdh-1 ! for z=+-D/2 this is zero anyways
          do iqz = -stars%mx3, stars%mx3
            irec3 = stars%ig(stars%kv2(1,irec2),stars%kv2(2,irec2),iqz)
            if ( irec3 /= 0 ) then ! use only stars within the g_max sphere
              c_ph(iqz,irec2) = stars%rgphs(stars%kv2(1,irec2),stars%kv2(2,irec2),iqz)
              qz = iqz * cell%bmat(3,3)
              VIz(iz,irec2) = VIz(iz,irec2) &
                            + vcons1(irec3) * c_ph(iqz,irec2) * & 
                              ( exp( ImagUnit * qz * z(iz) ) &
                              - cosh_quot(iz,irec2) * cosh( ImagUnit * qz * dh ) &
                              - sinh_quot(iz,irec2) * sinh( ImagUnit * qz * dh ) )
            end if
          enddo
        end do
      else ! numerical treatment necessary
        do iz = -nzdh+1, nzdh-1 ! for z=+-D/2 this is zero anyways
          do iqz = -stars%mx3, stars%mx3
            irec3 = stars%ig(stars%kv2(1,irec2),stars%kv2(2,irec2),iqz)
            if ( irec3 /= 0 ) then ! use only stars within the g_max sphere
              c_ph(iqz,irec2) = stars%rgphs(stars%kv2(1,irec2),stars%kv2(2,irec2),iqz)
              qz = iqz * cell%bmat(3,3)
              VIz(iz,irec2) = VIz(iz,irec2) &
                            + vcons1(irec3) * c_ph(iqz,irec2) * &
                              ( exp( ImagUnit * qz * z(iz) ) &
                              - exp_m(iz,irec2) * exp(   ImagUnit * qz * dh ) &
                              - exp_p(iz,irec2) * exp( - ImagUnit * qz * dh ) )
            end if
          enddo
        end do
      end if
    end do


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



  subroutine VYukawaModify( stars, vacuum, cell, sym, input, mpi, atoms, sphhar, dimension, oneD, noco, den, &
                            VYukawa )

    ! This subroutine adds a potential to the previously computed Yukawa
    ! potential to ensure charge neutrality.
    ! The added potential itself is a solution to the modified Helmholtz
    ! equation with constant right-hand side, where the constant is chosen such
    ! that charge neutrality is obtained.

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
    type(t_dimension),  intent(in)    :: dimension
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
    call cdntot_integrate( stars, atoms, sym, vacuum, input, cell, oneD, VYukawaModification, q, qis, qmt, qvac, qtot, qistot  )
    q0 = qtot / cell%area
    ldh = input%preconditioning_param * cell%z1
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

    call Vmts( input, mpi, stars, sphhar, atoms, sym, cell, oneD, &
               VYukawaModification%pw(:,1), den%mt(:,0:,:,1), VYukawaModification%potdenType, &
               VYukawaModification%mt(:,0:,:,1) )


    ! APPLYING THE MODIFICATION TO THE YUKAWA POTENTIAL

    call VYukawa%AddPotDen( VYukawa, VYukawaModification )


  end subroutine VYukawaModify



end module m_VYukawaFilm
