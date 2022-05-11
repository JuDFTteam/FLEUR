!-----------------------------------------------------------------------------------------------------------------------------------
! Forschungszentrum Jülich, juPhon Plugin for the FLEUR program
!-----------------------------------------------------------------------------------------------------------------------------------
!
! MODULE: Calculation of first variation of charge density.
!
!> @author
!> Christian-Roman Gerhorst
!
!> @brief
!> This module contains almost all routines for calculating the linear variation of the density including core-tails.
!>
!> @note
!> Additional information and formulas pointing out the routines of this module can be found within this
!> <a href='jpDens1stVar.pdf'>document</a>.
!-----------------------------------------------------------------------------------------------------------------------------------
module m_jpDens1stVar

    use m_constants

  implicit none

  contains

  !---------------------------------------------------------------------------------------------------------------------------------
  !> @author
  !> Christian-Roman Gerhorst, Forschungszentrum Jülich: IAS1 / PGI1
  !>
  !> @brief
  !> Calculates the valence part of the interstitial first variation of the density with the double sum method
  !>
  !> @details
  !> See Equation 7.27 (dissertation CRG)
  !>
  !> @note
  !> k-loop ist outwards
  !>
  !> @param[in]    cell       : Unit cell type, see types.f90.
  !> @param[in]    results    : Results type, see types.f90.
  !> @param[in]    nobd       : Number of occupied bands per k-point and spin
  !> @param[in]    ikpt       : Index of k-point in k-point set
  !> @param[in]    iqpt       : Index of q in q-point set
  !> @param[in]    ikpq       : Index of k + q (k' backfolded) in k-point set
  !> @param[in]    idir       : Index of displacement direction
  !> @param[in]    nv         : Number of LAPW G-basis vectors for given k-point.
  !> @param[in]    GbasVec    : G-basis vectors
  !> @param[in]    z          : Kohn\--Sham eigenvectors.
  !> @param[in]    z1nG       : Linear variation of the Kohn\--Sham eigenvectors
  !> @param[in]    gdp2iLimi  : Stores the min and maxvals of gdp2Ind
  !> @param[in]    gdp2Ind    : Stores the index of a potential and density G-Vector. The dimensions are the G-vector components.
  !> @param[in]    mapGbas    : For various k-points G-basis vectors occur more than once, thus they are only stored once in juPhon.
  !> @param[in]    kpq2kPrVec : Backfolding vector from k + q in 2nd Brillouin zone to k' in 1st Brillouin zone
  !> @param[inout] rho1IR     : Plane-wave interstitial coefficients of the linear variation of the interstitial density
  !---------------------------------------------------------------------------------------------------------------------------------
  subroutine calcRho1IRValDS( cell, results, nobd, nv, ikpt, iqpt, ikpq, idir, GbasVec, z, z1nG, rho1IR, gdp2Ind, mapGbas,   &
      & gdp2iLimi, kpq2kPrVec )

    use m_types, only : t_cell, t_results

    implicit none

    ! Type parameter
    type(t_cell),    intent(in)    :: cell
    type(t_results), intent(in)    :: results

    ! Scalar parameter
    integer,         intent(in)    :: nobd
    integer,         intent(in)    :: ikpt
    integer,         intent(in)    :: iqpt
    integer,         intent(in)    :: ikpq
    integer,         intent(in)    :: idir

    ! Array parameter
    integer,         intent(in)    :: nv(:)
    integer,         intent(in)    :: GbasVec(:, :)
    complex,        intent(in)    :: z(:, :)
    complex,         intent(in)    :: z1nG(:, :, :)
    integer,         intent(in)    :: gdp2iLimi(2, 3)
    integer,         intent(in)    :: gdp2Ind(gdp2iLimi(1, 1):gdp2iLimi(2, 1), gdp2iLimi(1, 2):gdp2iLimi(2, 2),                    &
                                             &gdp2iLimi(1, 3):gdp2iLimi(2, 3))
    integer,         intent(in)    :: mapGbas(:, :, :)
    integer,         intent(in)    :: kpq2kPrVec(:, :, :)
    complex,         intent(inout) :: rho1IR(:, :)

    ! Scalar variables
    integer                        :: iGpp
    integer                        :: iGp
    integer                        :: iband
    integer                        :: Gind
    real                           :: prf

    ! Array variables
    integer                        :: Gvec(3)

    ! WARNING: If we have a q, which is not zero, then we can get G-vectors which have a greater norm than 2kmax. We have
    !  dimensioned gdp2Ind bigger, sothat they are found but they should be zero.
    ! WARNING: Time-reversal symmetry is broken for magnetism in combination with spin-orbit coupling!!!!! factor 2 instead of 4 we
    ! have to rethink about it!!!!

    ! factor 2 from product rule and time reversal symmetry and factor 2 from spin degeneracy
    prf = 4. / cell%omtil
    if (.FALSE.) then
      if (ikpt.eq.1.and.idir.eq.1) then
        open(209,file='000_rho1_pwsplit',form='FORMATTED',position='append',action='WRITE',status='replace')
      else
        open(209,file='000_rho1_pwsplit',form='FORMATTED',position='append',action='WRITE',status='unknown')
      end if
    end if
    do iGpp = 1, nv(ikpq)
      do iGp = 1, nv(ikpt)
        ! G = G" + G' + G_bf
        Gvec(:) = GbasVec(:, mapGbas(iGpp, ikpq, 1)) - GbasVec(:, mapGbas(iGp, ikpt, 1)) +  kpq2kPrVec(:, ikpt, iqpt)
        Gind = gdp2Ind(Gvec(1), Gvec(2), Gvec(3))
        if ( Gind /= 0 ) then
          do iband = 1, nobd

            rho1IR(Gind, idir) = rho1IR(Gind, idir) + prf * results%w_iks(iband, ikpt, 1) &
                                                        & * conjg( z(iGp, iband) ) &
                                                        & * z1nG(iGpp, iband, idir)
            if (.FALSE..and.iGpp.eq.1.and.iGp.eq.1) then
              write(209,*) ikpt, iband!iGpp, iGp, idir, iband
              write(209,*) prf * results%w_iks(iband, ikpt, 1)
              !write(209,*) conjg( z(iGp, iband) )
              !write(209,*) z1nG(iGpp, iband, idir)
            end if
          end do ! iband
        else
          write(*, *) 'Error in G-vector assignment.'
        end if
      end do ! iGp
    end do ! iGpp
    if (.FALSE.) then
      close(209)
    end if

  end subroutine calcRho1IRValDS

  !^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
  !     In this subroutine the star function expansion coefficients of
  !     the plane wave charge density is determined.
  !
  !     This subroutine is called for each k-point and each spin.
  !
  !
  !     Two methods are implemented to calculate the charge density
  !     1) which uses the FFT. The effort in calculating the charge
  !        density is proportional to M * N * log(N) , M being number of
  !        states and N being number of plane waves. This is the method
  !        which we use for production runs
  !     2) the traditional method for calculating the charge density
  !        using the double summation. In this case the effort scales as
  !        M * N * N. The method is only used for test purposes or for
  !        special cases.
  !
  !
  !     INPUT:    eigen vectors
  !               reciprocal lattice information
  !               Brillouine zone sampling
  !               FFT information
  !
  !     OUTPUT:   qpw(s)
  !               1) using FFT
  !
  !                2) traditional method
  !
  !                             -1             ef
  !                qpw  (g) = vol * sum{ sum{ sum{ sum{ w(k) * f(nu) *
  !                                  sp   k    nu   g'
  !                                     *
  !                                    c(g'-g,nu,k) * c(g',nu,k) } } } }
  !                or :
  !                             -1             ef
  !                qpw  (g) = vol * sum{ sum{ sum{ sum{ w(k) * f(nu) *
  !                                  sp   k    nu   g'
  !                                     *
  !                                    c(g',nu,k) * c(g'+g,nu,k) } } } }
  !
  !                qpw(g) are actuall
  !
  !                the weights w(k) are normalized: sum{w(k)} = 1
  !                                                  k                -6
  !                         a) 1                           for kT < 10
  !                f(nu) = {                           -1             -6
  !                         b){ 1+exp(e(k,nu) -ef)/kt) }   for kt >=10
  !
  !
  !                                      Stefan Bl"ugel, JRCAT, Feb. 1997
  !                                      Gustav Bihlmayer, UniWien
  !
  !     In non-collinear calculations the density becomes a hermitian 2x2
  !     matrix. This subroutine generates this density matrix in the
  !     interstitial region. The diagonal elements of this matrix
  !     (n_11 & n_22) are stored in qpw, while the real and imaginary part
  !     of the off-diagonal element are store in cdom.
  !
  !     Philipp Kurz 99/07
  !^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
  !
  !
  !this is fthe fft variant of the preivvous method
  ! todo review declarations of variables

  !------->          ABBREVIATIONS
  !
  !     rhon  : charge density in real space
  !     ne    : number of occupied states
  !     nv    : number of g-components in eigenstate
  !     cv=z  : wavefunction in g-space (reciprocal space)
  !     psir   : wavefunction in r-space (real-space)
  !     cwk   : complex work array: charge density in g-space (as stars)
  !     qpw   : charge density stored as stars
  !     trdchg: logical key, determines the mode of charge density
  !             calculation: false (default) : fft
  !                          true            : double sum over stars
  !     we    : weights for the BZ-integration for a particular k-point
  !     omtil : volume (slab) unit cell, between -.5*D_tilde and +.5*D_tilde
  !     k1   : reciprocal lattice vectors G=G(k1,k2,k3) for wavefunction
  !     k2   :                             =k1*a_1 + k2*a_2 + k3*a_3
  !     k3   : where a_i= Bravais lattice vectors in reciprocal space
  !             kwi, depend on k-point.
  !     kq1d  : dimension of the charge density FFT box in the pos. domain
  !     kq2d  : defined in dimens.f program (subroutine apws).1,2,3 indicate
  !     kq3d  ; a_1, a_2, a_3 directions.
  !     kq(i) : i=1,2,3 actual length of the fft-box for which FFT is done.
  !     nstr  : number of members (arms) of reciprocal lattice (g) vector
  !             of each star.
  !     ng3_fft: number of stars in the  charge density  FFT-box
  !     ng3   : number of 3 dim. stars in the charge density sphere defined
  !             by gmax
  !     kmxq_fft: number of g-vectors forming the ng3_fft stars in the
  !               charge density sphere
  !     kimax : number of g-vectors forming the ng3 stars in the gmax-sphere
  !     iv1d  : maps vector (k1,k2,k3) of wave function into one
  !             dimensional vector of cdn-fft box in positive domain.
  !     ifftq3d: elements (g-vectors) in the charge density  FFT-box
  !     igfft : pointer from the g-sphere (stored as stars) to fft-grid
  !             and     from fft-grid to g-sphere (stored as stars)
  !     pgfft : contains the phases of the g-vectors of sph.
  !     isn   : isn = +1, FFT transform for g-space to r-space
  !             isn = -1, vice versa
  ! ist array is a clever array defined in the declarations
  subroutine calcRho1IRValFFT( cell, stars, results, ikpt, ngdp2km, nobd, ne, nv, gpot, gbas, z0, z1, rho1IR )

    use m_juDFT_stop, only : juDFT_error, juDFT_warn
    use m_types
    use m_fft_interface

    implicit none

    type(t_cell),                          intent(in)    :: cell
    type(t_stars),                         intent(in)    :: stars
    type(t_results),                       intent(in)    :: results

    integer,                               intent(in)    :: ikpt
    integer,                               intent(in)    :: ngdp2km
    integer,                               intent(in)    :: nobd
    integer,                               intent(in)    :: ne
    integer,                               intent(in)    :: nv

    integer,                               intent(in)    :: gpot(:, :) ! gpot should be put only until ngdp2km
    integer,                               intent(in)    :: gbas(:, :) ! gbas should be put in with the ilst array!!!!
    complex,                              intent(in)    :: z0(:, :)
    complex,                               intent(in)    :: z1(:, :)
    complex,                               intent(inout) :: rho1IR(:) ! we have to see whether an inout or a out is better for performance

    integer                                              :: ifftq1
    integer                                              :: ifftq2
    integer                                              :: ifftq3
    integer                                              :: ifftq2d
    integer                                              :: ifftq3d
    integer                                              :: idir
    integer                                              :: iG
    integer                                              :: iG1
    integer                                              :: iG2
    integer                                              :: iG3
    integer                                              :: iband
    logical                                              :: forw
    real                                                 :: scaling
#ifdef DEBUG_MODE
    complex                                                  :: q0
    complex :: q0currentK
#endif

    integer,                    parameter                :: ist(-1:1) = [1, 0, 0]
    integer                                              :: length_zfft(3)
    real,          allocatable                           :: prf(:)
    integer,       allocatable                           :: iv1d(:)
    integer,       allocatable                           :: igfft(:)
    complex,       allocatable                           :: z0fft(:)
    complex,       allocatable                           :: z1fft(:)
    complex,       allocatable                           :: rho1g(:)
    integer                                                  :: nfft(3)
    integer                                                   :: ifftd
    integer                                                   :: gabs(3)
    integer :: idir1
    integer :: ii
    integer :: iG0


    INTEGER :: stars_kq1_fft,stars_kq2_fft,stars_kq3_fft


    ! setup FFT box dimensions
    ifftq1  = stars_kq1_fft
    ifftq2  = stars_kq1_fft * stars_kq2_fft
    ifftq3  = stars_kq1_fft * stars_kq2_fft * stars_kq3_fft
    ifftq3d = stars_kq1_fft * stars_kq2_fft * stars_kq3_fft
    ifftq2d = stars_kq1_fft * stars_kq2_fft

    ! Prefactor
    allocate( prf(nobd) )
    prf(:nobd) = 4 * results%w_iks(:nobd, ikpt, 1) / cell%omtil


    ! create mapping array from wavefunction expansion coefficients to FFT box
    allocate( iv1d(nv) )
    do iG = 1, nv
      ! -k1d <= iG1 <= k1d
      ! -k2d <= iG2 <= k2d
      ! -k3d <= iG3 <= k3d
      iG1 = gBas(1, iG)
      iG2 = gBas(2, iG)
      iG3 = gBas(3, iG)

      ! L,M,N LATTICE POINTS OF G-VECTOR IN POSITIVE DOMAIN
      ! (since charge density box = two times charge density box
      ! wrap arround error should not occur )
      ! 0 <= iG1 <= 2 * k1 - 1 = kq1_fft - 1
      ! 0 <= iG2 <= 2 * k2 - 1 = kq2_fft - 1
      ! 0 <= iG3 <= 2 * k3 - 1 = kq3_fft - 1
      iG1 = iG1 + stars_kq1_fft * ist( isign( 1, iG1) )
      iG2 = iG2 + stars_kq2_fft * ist( isign( 1, iG2) )
      iG3 = iG3 + stars_kq3_fft * ist( isign( 1, iG3) )

      iv1d(iG) = iG3 * ifftq2 + iG2 * ifftq1 + iG1
    end do

    ! create mapping array from FFT box to density plane-wave expansion coefficients
    !allocate( igfft(ngdp2km) )
    ! Probably it might be a good idea to define iggft outside because it is required several times!
    allocate( igfft(ngdp2km) )
    nfft = [stars_kq1_fft, stars_kq2_fft, stars_kq3_fft]
    gabs = 0
    igfft = 0
    do iG = 1, ngdp2km
#if DEBUG_MODE
      if ( all(gpot(:, iG) == 0 ) ) then
        iG0 = iG !for q0 test
      end if
#endif
      do idir1 = 1, 3
        if ( gpot(idir1, iG) >= 0 ) then
          gabs(idir1) = gpot(idir1, iG)
        else
          gabs(idir1) = gpot(idir1, iG) + nfft(idir1)
        end if
      end do
      igfft(iG) = gabs(1) + gabs(2) * nfft(1) + gabs(3) * nfft(1) * nfft(2)
    end do


    !kqi_fft are the dimensions of a fft box which are suitable for a 2 kmax density
    allocate( z0fft( 0 : stars_kq1_fft * stars_kq2_fft * stars_kq3_fft - 1 ), &
              z1fft( 0 : stars_kq1_fft * stars_kq2_fft * stars_kq3_fft - 1 ), &
              rho1g( 0 : stars_kq1_fft * stars_kq2_fft * stars_kq3_fft - 1 ) ) ! rho1 on fft Grid

    ! initialize first variation of MT charge density on FFT grid
    rho1g = cmplx( 0.0, 0.0 )

    DO iband = 1, nobd ! todo really nobd or ne?

      ! FFT transform c_iband,k(g) --> psi_iband,k(r), for each k-point and each iband-state
      z0fft = cmplx( 0.0, 0.0 )
      z1fft = cmplx( 0.0, 0.0 )

      z0fft( iv1d( :nv ) ) =  z0( :nv, iband)
      z1fft( iv1d( :nv ) ) =  z1( :nv, iband)


      ! do (real) inverse FFT; notice that the arrays z0/1fft are filled from 0 to ifftq3-1, but starts at -ifftq2 to give work space
      ! for rfft
      forw = .false. ! is equivalent to isn = 1
      ! FFT transform
      length_zfft(1) = stars_kq1_fft
      length_zfft(2) = stars_kq2_fft
      length_zfft(3) = stars_kq3_fft

      call fft_interface(3, length_zfft, z0fft, forw)
      call fft_interface(3, length_zfft, z1fft, forw)

      ! This is the comprehensive way of writing the product of two complex numbers
      !DO ir = 0,ifftq3d-1
      !  rho1r(ir) = rho1r(ir) + prf(iband) * ( psi0r(ir) * psi1pr(ir) + psi0i(ir) * psi1pi(ir) )
      !  rho1i(ir) = rho1i(ir) + prf(iband) * ( psi0r(ir) * psi1pi(ir) - psi0i(ir) * psi1pr(ir) )
      !ENDDO

      ! calculate first variation of density in real space on FFT grid for band iband
      rho1g( 0 : ifftq3d - 1 ) = rho1g( 0 : ifftq3d - 1 ) + prf(iband) * conjg(z0fft( 0 : ifftq3d - 1 )) * z1fft( 0 : ifftq3d - 1 )

!
!
    end do ! loop over occupied bands
    deallocate( iv1d )
    deallocate( prf, z0fft, z1fft )
    if ( any( abs( aimag( rho1g(:) ) )  > 1e-9 ) ) then
      call juDFT_warn( 'IR variation of density has imaginary components.', calledby = 'calcRho1IRValFFT', hint = 'Fix bug!' )
    end if

    ! set numerical noise to 0
    rho1g(:) = cmplx( real( rho1g(:) ), 0.0 )

    ! perform back FFT transform
    forw = .true. ! is equivalent to isn = -1
    call fft_interface( 3, length_zfft, rho1g, forw )

    scaling = 1.0 / real(ifftq3)
#ifdef DEBUG_MODE
    q0 = rho1IR(iG0)
#endif
    do iG = 1 , ngdp2km
      rho1IR(iG) = rho1IR(iG) + scaling * rho1g( igfft(iG) )
    end do

    deallocate(rho1g, igfft)

#ifdef DEBUG_MODE
    ! Consistency check for FFT, calculate the G = 0 component with double sum method and the FFT and compare
    q0currentK = 0.0
    do iband = 1, nobd ! only occupied?
      q0currentK = q0currentK + 4 * results%w_iks(iband, ikpt, 1) * dot_product( z0(:nv, iband), z1(:nv, iband) )
    end do
    q0currentK = q0currentK / cell%omtil
    q0 = q0 + q0currentK

    if ( abs( q0 ) > 1e-9) then
      if ( abs( q0 - rho1IR(iG0) ) > 1e-6 ) then
        call juDFT_warn( 'Consistency check of IR density first variation failed.', calledby = 'calcRho1IRValFFT', hint = 'Check FFT!' )
      end if
    end if
#endif

  end subroutine calcRho1IRValFFT

  !---------------------------------------------------------------------------------------------------------------------------------
  !> @author
  !> Christian-Roman Gerhorst, Forschungszentrum Jülich: IAS1 / PGI1
  !>
  !> @brief
  !> Non-LO spherical products of acof and bcof
  !>
  !> @note
  !> Based on adjusted from rhomt from new fleur new at density change
  !>
  !> @details
  !> See Equation 7.33b and Table 7.1 (dissertation CRG)
  !>
  !> @param[in]     atoms    : Atoms type, see types.f90
  !> @param[in]     nobd     : Number of occupied bands per k-point and spin
  !> @param[in]     acof     : (Any) large matching coefficient related to u
  !> @param[in]     bcof     : (Any) large matching coefficient related to udot
  !> @param[in]     acof2cjg : (Any) large matching coefficient related to u which is to be conjugated in the routine
  !> @param[in]     bcof2cjg : (Any) large matching coefficient related to udot which is to be conjugated in the routine
  !> @param[in]     we       : Unperturbed occupation number
  !> @param[inout]  uu       : Spherical (l = 0) product of acof and acof
  !> @param[inout]  dd       : Spherical (l = 0) product of bcof and acof
  !> @param[inout]  du       : Spherical (l = 0) product of acof and bcof
  !> @param[inout]  ud       : Spherical (l = 0) product of bcof and bcof
  !---------------------------------------------------------------------------------------------------------------------------------
  subroutine rhomt( atoms, we, nobd, acof2cjg, bcof2cjg, acof, bcof, uu, dd, du, ud )

    use m_types, only : t_atoms

    implicit none

    ! Type parameter
    type(t_atoms), intent(in)    :: atoms

    ! Scalar parameter
    integer,       intent(in)    :: nobd

    ! Array parameters (inout because of k-point sum in SternheimerSCC)
    complex,       intent(in)    :: acof(:,0:,:)
    complex,       intent(in)    :: bcof(:,0:,:)
    complex,       intent(in)    :: acof2cjg(:,0:,:)
    complex,       intent(in)    :: bcof2cjg(:,0:,:)
    real,          intent(in)    :: we(:)
    complex,       intent(inout) :: uu(0:, :)
    complex,       intent(inout) :: dd(0:, :)
    complex,       intent(inout) :: du(0:, :)
    complex,       intent(inout) :: ud(0:, :)

    ! Local scalar variables
    integer                      :: iband
    integer                      :: oqn_l
    integer                      :: mqn_m
    integer                      :: lm
    integer                      :: itype
    integer                      :: ieqat
    integer                      :: iatom


    ! the factor i^l or (i^l_p)*, respectively, which actually is part of the abcof cancels away, because we are diagonal in l !
    iatom = 0
    do itype = 1, atoms%ntype
      do ieqat = 1, atoms%neq(itype)
        iatom = iatom + 1
        do oqn_l = 0, atoms%lmax(itype)
          do mqn_m = -oqn_l, oqn_l
            lm = oqn_l * (oqn_l + 1) + mqn_m
            do iband = 1, nobd
              uu(oqn_l, iatom) = uu(oqn_l, iatom) + 2 * we(iband) * conjg( acof2cjg(iband, lm, iatom) ) * acof(iband, lm, iatom)
              dd(oqn_l, iatom) = dd(oqn_l, iatom) + 2 * we(iband) * conjg( bcof2cjg(iband, lm, iatom) ) * bcof(iband, lm, iatom)
              du(oqn_l, iatom) = du(oqn_l, iatom) + 2 * we(iband) * conjg( bcof2cjg(iband, lm, iatom) ) * acof(iband, lm, iatom)
              ud(oqn_l, iatom) = ud(oqn_l, iatom) + 2 * we(iband) * conjg( acof2cjg(iband, lm, iatom) ) * bcof(iband, lm, iatom)
            end do ! iband
          end do ! mqn_m
          if (.FALSE.) then
            write(109,*) oqn_l, oqn_l, 1, 1, 1
            write(109,*) uu(oqn_l,1)
            write(109,*) oqn_l, oqn_l, 1, 1, 2
            write(109,*) ud(oqn_l,1)
            write(109,*) oqn_l, oqn_l, 1, 2, 1
            write(109,*) du(oqn_l,1)
            write(109,*) oqn_l, oqn_l, 1, 2, 2
            write(109,*) dd(oqn_l,1)
          end if
        end do ! oqn_l
      end do ! ieqat
    end do ! itype

  end subroutine rhomt

  !---------------------------------------------------------------------------------------------------------------------------------
  !> @author
  !> Christian-Roman Gerhorst, Forschungszentrum Jülich: IAS1 / PGI1
  !>
  !> @brief
  !> This subroutine is the equivalent of rhomt for the local orbital
  !> contributions to the charge.
  !> aclo,bclo,cclo are the equivalents of uu,ud,dd in rhomt
  !> p.kurz sept. 1996
  !> @note
  !> adjusted from rhomt from new fleur new at density change
  !---------------------------------------------------------------------------------------------------------------------------------
  subroutine rhomtlo( atoms, ne, we, acof, bcof, ccof, acof2cjg, bcof2cjg, ccof2cjg, aclo, bclo, cclo )

    use m_types, only : t_atoms

    implicit none

    ! Type parameter
    type(t_atoms), intent(in)     :: atoms

    ! Scalar parameter
    integer,       intent(in)     :: ne

    ! Array parameter
    real,          intent (in)    :: we(:)!(nobd)
    complex,       intent (in)    :: acof(:,0:,:)!(nobd,0:dimension%lmd,atoms%nat)
    complex,       intent (in)    :: bcof(:,0:,:)!(nobd,0:dimension%lmd,atoms%nat)
    complex,       intent (in)    :: ccof(-atoms%llod:,:,:,:)!(-atoms%llod:llod,nobd,atoms%nlod,atoms%nat)
    complex,       intent (in)    :: acof2cjg(:,0:,:)!(nobd,0:dimension%lmd,atoms%nat)
    complex,       intent (in)    :: bcof2cjg(:,0:,:)!(nobd,0:dimension%lmd,atoms%nat)
    complex,       intent (in)    :: ccof2cjg(-atoms%llod:,:,:,:)!(-atoms%llod:llod,nobd,atoms%nlod,atoms%nat)
    complex,       intent (inout) :: aclo(:, :)
    complex,       intent (inout) :: bclo(:, :)
    complex,       intent (inout) :: cclo(:, :, :)

    ! Local scalar variable
    integer                       :: itype
    integer                       :: ieqat
    integer                       :: iatom
    integer                       :: oqn_l
    integer                       :: mqn_m
    integer                       :: lm
    integer                       :: iband
    integer                       :: ilo
    integer                       :: jlo

    iatom = 0
    !---> loop over atoms
    do itype = 1, atoms%ntype
      do ieqat = 1, atoms%neq(itype)
        iatom = iatom + 1
        !--->       loop over the local orbitals
        do ilo = 1, atoms%nlo(itype)
          oqn_l = atoms%llo(ilo, itype)
          !--->       contribution of cross terms flapw - local orbitals
          do mqn_m = -oqn_l, oqn_l
            lm = oqn_l * ( oqn_l + 1) + mqn_m
            do iband = 1, ne ! todo really ne or nobd?
              ! the factor i^l or (i^l_p)*, respectively, which actually is part of the abcof cancels away, because we are diagonal
              ! in l!.
              !todo why two factors of 2?
              aclo( ilo, iatom ) = aclo( ilo, iatom ) + 2 * we(iband) * 2 *  conjg( acof2cjg(iband, lm, iatom) ) &
                & * ccof(mqn_m, iband, ilo, iatom)
              bclo(ilo, iatom) = bclo(ilo, iatom) + 2 * we(iband) * 2 * conjg( bcof2cjg(iband, lm, iatom) ) &
                & * ccof(mqn_m, iband, ilo, iatom)
            end do
          end do
          !--->       contribution of local orbital - local orbital terms
          !--->       loop over lo'
          do jlo = 1, atoms%nlo(itype)
            if ( atoms%llo(jlo, itype) == oqn_l ) then
              do mqn_m = -oqn_l, oqn_l
                do iband = 1, ne
                  ! the factor i^l or (i^l_p)*, respectively, which actually is part of the abcof cancels away, because we are
                  ! diagonal in l!.
                  cclo(jlo, ilo, iatom) = cclo(jlo, ilo, iatom) + 2 * we(iband) * &
                       & conjg( ccof2cjg(mqn_m, iband, jlo, iatom) ) * ccof(mqn_m, iband, ilo ,iatom)
                end do
              end do
            end if
          end do
        end do
      end do
    end do

  end subroutine rhomtlo

  !---------------------------------------------------------------------------------------------------------------------------------
  !> @author
  !> Christian-Roman Gerhorst, Forschungszentrum Jülich: IAS1 / PGI1
  !>
  !> @brief
  !> Non-LO non-spherical product of acof bcof for linear density variation
  !>
  !> @details
  !> See Equation 7.33b and Table 7.1 (dissertation CRG)
  !>
  !> @param[in]     atoms    : Atoms type, see types.f90
  !> @param[in]     nobd     : Number of occupied bands per k-point and spin
  !> @param[in]     we       : Unperturbed occupation number
  !> @param[in]     acof     : (Any) large matching coefficient related to u
  !> @param[in]     bcof     : (Any) large matching coefficient related to udot
  !> @param[in]     acof2cjg : (Any) large matching coefficient related to u which is to be conjugated in the routine
  !> @param[in]     bcof2cjg : (Any) large matching coefficient related to udot which is to be conjugated in the routine
  !> @param[inout]  ddnmt    : Non-spherical (l> 0) product of bcof and bcof
  !> @param[inout]  dunmt    : Non-spherical (l> 0) product of acof and bcof
  !> @param[inout]  udnmt    : Non-spherical (l> 0) product of bcof and acof
  !> @param[inout]  uunmt    : Non-spherical (l> 0) product of acof and acof
  !>
  !---------------------------------------------------------------------------------------------------------------------------------
  subroutine rhonmt( atoms, we, nobd, acof, bcof, acof2cjg, bcof2cjg, uunmt, ddnmt, udnmt, dunmt )

    use m_types, only : t_atoms
    use m_gaunt, only : gaunt1

    implicit none

    ! Type parameters
    type(t_atoms),    intent(in)    :: atoms

    ! Scalar parameters
    integer,          intent(in)    :: nobd

    ! Array parameters
    real,             intent(in)    :: we(:)
    complex,          intent(in)    :: acof(:,0:,:)
    complex,          intent(in)    :: bcof(:,0:,:)
    complex,          intent(in)    :: acof2cjg(:,0:,:)
    complex,          intent(in)    :: bcof2cjg(:,0:,:)
    complex,          intent(inout) :: ddnmt(0:,0:,:,:)
    complex,          intent(inout) :: dunmt(0:,0:,:,:)
    complex,          intent(inout) :: udnmt(0:,0:,:,:)
    complex,          intent(inout) :: uunmt(0:,0:,:,:)

    ! Local scalar variables
    complex                         :: cconst
    complex                         :: cil
    real                            :: coef
    integer                         :: l
    integer                         :: lmv
    integer                         :: lm
    integer                         :: lmp
    integer                         :: lp
    integer                         :: lphi
    integer                         :: lplow
    integer                         :: lplow0
    integer                         :: lv
    integer                         :: mp
    integer                         :: mv
    integer                         :: ieqat
    integer                         :: iatom
    integer                         :: iband
    integer                         :: itype
    integer                         :: m
    integer                         :: lphi0

    iatom = 0
    do itype = 1, atoms%ntype
      do ieqat = 1, atoms%neq(itype)
        iatom = iatom + 1
        do lv = 1, atoms%lmax(itype)
          do mv = -lv, lv
            lmv = lv * ( lv + 1 ) + mv + 1
            do l = 0, atoms%lmaxd
              ! Gaunt selection rules |l - lv| <= l' <= |l + lv|
              lplow0 = abs( l - lv )
              lphi0 =  l + lv
              ! lp should not be larger than the numerical cutoff
              lphi = min(lphi0, atoms%lmax(itype))
              ! l + lv + lphi should be even and lphi the last loop iteration sothat it is reduced by 1 unlike this is the case
              lphi = lphi - mod(l + lv + lphi, 2)
              do m = -l, l
                lm = l * ( l + 1 ) + m
                ! Gaunt selection rule for the magnetic quantum numbers
                mp = m - mv
                ! lmin = max(|l' - lv|, |m' - mv|) to ensure that |m = m' - mv| <= l while fulfilling m = m' - mu.
                ! Imagine for example l' = 1, lamda = 1, m' = -1,mv = 1. Independent of the selection rule, l >= |m|
                lplow = max(lplow0, abs(mp))
                ! Serves only for ensuring l + l' + lambda = even, if either lmxx is odd or lmin odd, with the modulo, lphi gives
                ! its oddness or eveness to lplow, so as l' is interated in steps of size 2 l + l' is always even, because
                ! either l and l' are even or l and l' are odd, therefore l + l' + lambda is even because the eveness is ensured by
                ! the modulo operation
                lplow = lplow + mod( abs(lphi - lplow), 2)
                do lp = lplow, lphi, 2
                  ! this is the factor i^l or i^l_p, respectively, which actually is part of the abcof. For non-diagonal
                  ! contributions we have to calculate it explicetly. For diagonal contributions it cancels away because l = l_p!
                  cil = ImagUnit**( l - lp )
                  lmp = lp * ( lp + 1 ) + mp
                  !     -----> gaunt's coefficient
                  coef = gaunt1( l, lv, lp, m, mv, mp, atoms%lmaxd )
                  cconst = coef * cil
                  do iband = 1, nobd
                    uunmt(lp, l, lmv, iatom) = uunmt(lp, l, lmv, iatom) + 2 * we(iband) * cconst &
                                                                          & * conjg(acof2cjg(iband, lmp, iatom)) &
                                                                          & * acof(iband, lm, iatom)
                    ddnmt(lp, l, lmv, iatom) = ddnmt(lp, l, lmv, iatom) + 2 * we(iband) * cconst &
                                                                          & * conjg(bcof2cjg(iband, lmp, iatom)) &
                                                                          & * bcof(iband, lm, iatom)
                    udnmt(lp, l, lmv, iatom) = udnmt(lp, l, lmv, iatom) + 2 * we(iband) * cconst &
                                                                          & * conjg(acof2cjg(iband, lmp, iatom)) &
                                                                          & * bcof(iband, lm, iatom)
                    dunmt(lp, l, lmv, iatom) = dunmt(lp, l, lmv, iatom) + 2 * we(iband) * cconst &
                                                                          & * conjg(bcof2cjg(iband, lmp, iatom)) &
                                                                          & * acof(iband, lm, iatom)
                  end do ! iband
                  if (.FALSE..and.m.eq.l) then
                    write(109,*) lp, l, lmv, 1, 1
                    write(109,*) real(uunmt(lp,l,lmv,1))
                    write(109,*) lp, l, lmv, 1, 2
                    write(109,*) real(udnmt(lp,l,lmv,1))
                    write(109,*) lp, l, lmv, 2, 1
                    write(109,*) real(dunmt(lp,l,lmv,1))
                    write(109,*) lp, l, lmv, 2, 2
                    write(109,*) real(ddnmt(lp,l,lmv,1))
                  end if
                end do ! lp
              end do ! m
            end do ! l
          end do ! mv
        end do ! lv
      end do ! ieqat
    end do ! itype

  end subroutine rhonmt

  !---------------------------------------------------------------------------------------------------------------------------------
  !> @author
  !> Christian-Roman Gerhorst, Forschungszentrum Jülich: IAS1 / PGI1
  !>
  !> @brief
  !> This subroutine is the equivalent of rhomt for the local orbital
  !> contributions to the charge.
  !> acnmt, bcnmt, ccnmt are the equivalents of uunmt, ddnmt, udnmt dunmt
  !> in rhonmt
  !> p.kurz sept. 1996
  !>
  !> @details
  !> See Table 7.1
  !>
  !> @note
  !> adjusted from rhomt from new fleur new at density change
  !---------------------------------------------------------------------------------------------------------------------------------
  subroutine rhonmtlo( atoms, ne, we, acof, bcof, ccof, acof2cjg, bcof2cjg, ccof2cjg, acnmt, bcnmt, ccnmt )

    use m_types, only : t_atoms
    use m_gaunt, only : gaunt1

    implicit none

    type(t_atoms), intent(in)     :: atoms

    !     .. Scalar Arguments .   .
    integer,       intent (in)    :: ne

    !     .. Array Arguments ..
    real,          intent (in)    :: we(:)!(nobd)
    complex,       intent (in)    :: acof(:,0:,:)!(nobd,0:dimension%lmd,atoms%nat)
    complex,       intent (in)    :: bcof(:,0:,:)!(nobd,0:dimension%lmd,atoms%nat)
    complex,       intent (in)    :: ccof(-atoms%llod:,:,:,:)!(-llod:llod,nobd,atoms%nlod,atoms%nat)
    complex,       intent (in)    :: acof2cjg(:,0:,:)!(nobd,0:dimension%lmd,atoms%nat)
    complex,       intent (in)    :: bcof2cjg(:,0:,:)!(nobd,0:dimension%lmd,atoms%nat)
    complex,       intent (in)    :: ccof2cjg(-atoms%llod:,:,:,:)!(-llod:llod,nobd,atoms%nlod,atoms%nat)
    complex,       intent (inout) :: acnmt(0:,:,:,:)
    complex,       intent (inout) :: bcnmt(0:,:,:,:)
    complex,       intent (inout) :: ccnmt(:,:,:,:)

    !     .. Local Scalars ..
    complex                       :: ci
    complex                       :: fact
    complex                       :: cf1
    integer                       :: i
    integer                       :: l
    integer                       :: lmp
    integer                       :: lmpp
    integer                       :: lo
    integer                       :: lop
    integer                       :: lp
    integer                       :: lpmax
    integer                       :: lpmax0
    integer                       :: lpmin
    integer                       :: lpmin0
    integer                       :: m
    integer                       :: lpp
    integer                       :: mp
    integer                       :: mpp
    integer                       :: na
    integer                       :: neqat0
    integer                       :: nn
    integer                       :: ntyp

    ci = cmplx(0.0,1.0)

    !---> for optimal performance consider only
    !---> those combinations of l,l',l'',m,m',m'' that satisfy the three
    !---> conditions for non-zero gaunt-coeff. i.e.
    !---> |l - l''| <= l' <= l + l'' (triangular condition)
    !---> m' + m'' = m and l + l' + l'' even
    neqat0 = 0
    do ntyp = 1, atoms%ntype
      !--->    loop over the lattice harmonics
      do lpp = 1, atoms%lmax(ntyp) ! only non-spherical part! so l=0 is ignored
        do mpp = -lpp, lpp
          lmpp = lpp * ( lpp + 1 ) + mpp
          do lo = 1, atoms%nlo(ntyp)
            l = atoms%llo(lo, ntyp)
            lpmin0 = abs( l - lpp )
            lpmax0 = l + lpp
            !--->             check that lpmax is smaller than the max l of the
            !--->             wavefunction expansion at this atom
            lpmax = min( lpmax0, atoms%lmax(ntyp) )
            !--->             make sure that l + l'' + lpmax is even
            lpmax = lpmax - mod( l + lpp + lpmax, 2 )
            do m = -l,l
              !--->                add flapw - local orbital cross-terms

              !--->                add terms containing gaunt1(l,lp,lpp,m,mp,mpp)
              !--->                note that gaunt1(l,lp,lpp,m,mp,mpp) computes the
              !--->                integral of conjg(y(l,m))*y(lp,mp)*y(lpp,mpp),
              !--->                however, since the gaunt coef. are real, this is
              !--->                the same as int. y(l,m)*conjg(y(lp,mp)*y(lpp,mpp))
              mp = m - mpp
              lpmin = MAX( lpmin0, ABS( mp ) )
              !--->                make sure that l + l'' + lpmin is even
              lpmin = lpmin + MOD( ABS( lpmax - lpmin ), 2 )
              !--->                loop over l'
              do lp = lpmin, lpmax, 2
                lmp = lp * ( lp + 1 ) + mp
                ! this is the factor i^l or (i^l_p)*, respectively, which actually is part of the abcof. For non-diagonal
                ! contributions we have to calculate it explicetly. For diagonal contributions it cancels away because l = l_p!
                fact = ( ci** ( l - lp ) ) * gaunt1( l, lp, lpp, m, mp, mpp, atoms%lmaxd )
                na = neqat0
                do nn = 1, atoms%neq(ntyp)
                  na = na + 1
                  do i = 1, ne
                    cf1 = fact * ccof(m, i, lo, na)
                    acnmt( lp, lo, lmpp, na ) = acnmt(lp, lo, lmpp, na) + 2 * we(i) * cf1 * conjg(acof2cjg(i, lmp, na))
                    bcnmt( lp, lo, lmpp, na ) = bcnmt(lp, lo, lmpp, na) + 2 * we(i) * cf1 * conjg(bcof2cjg(i, lmp, na))
                  end do
                end do
              end do

              !--->                add terms containing gaunt1(lp,l,lpp,mp,m,mpp)
              mp = m + mpp
              lpmin = max(lpmin0, abs(mp) )
              !--->                make sure that l + l'' + lpmin is even
              lpmin = lpmin + mod( abs( lpmax-lpmin ), 2 )
              !--->                loop over l'
              do lp = lpmin, lpmax, 2
                lmp = lp * ( lp + 1 ) + mp
                fact = (ci** ( lp - l )) * gaunt1(lp, l, lpp, mp, m, mpp, atoms%lmaxd)
                na = neqat0 ! ensures that iatom starts at the first atom of the respective type
                do nn = 1, atoms%neq(ntyp)
                  na = na + 1
                  do i = 1, ne
                    cf1 = fact *  conjg(ccof2cjg(m, i, lo, na))
                    acnmt(lp,lo,lmpp,na) = acnmt(lp,lo,lmpp,na) + 2 * we(i) * cf1 * acof(i, lmp, na)
                    bcnmt(lp, lo, lmpp, na) = bcnmt(lp, lo, lmpp, na) + 2 * we(i) * cf1 * bcof(i, lmp, na)
                  end do
                end do
              end do

              !--->                add local orbital - local orbital terms
              do lop = 1, atoms%nlo(ntyp)
                lp = atoms%llo(lop, ntyp)

                !--->                   add terms containing gaunt1(l,lp,lpp,m,mp,mpp)
                mp = m - mpp
                if ( ( abs( l -lpp ) <= lp ) .and. ( lp <= ( l + lpp ) ) .and.( mod( l + lp + lpp, 2 ) == 0 ) .and.                 &
                  & ( abs( mp ) <= lp) ) then
                  fact = (ci** (l-lp))*gaunt1(l, lp, lpp, m, mp, mpp, atoms%lmaxd)
                  na = neqat0
                  do nn = 1,atoms%neq(ntyp)
                    na = na + 1
                    do i = 1,ne
                      ccnmt(lop, lo, lmpp, na) =&
                           ccnmt(lop,lo,lmpp,na) + 2 * we(i) * fact * conjg( ccof2cjg(mp, i, lop, na) ) * ccof(m, i, lo, na)
                    end do
                  end do
                end if
              end do
            end do
          end do
        end do
      end do
      neqat0 = neqat0 + atoms%neq(ntyp)
    end do

  end subroutine rhonmtlo

  !---------------------------------------------------------------------------------------------------------------------------------
  !> @author
  !> Christian-Roman Gerhorst, Forschungszentrum Jülich: IAS1 / PGI1
  !>
  !> @brief
  !> Calls methods that calculate the combinations of the matching coefficients for the spherical and the non-spherical
  !> contributions to the density variation, including LOs.
  !>
  !> @details
  !> See Table 7.1 and Equation 7.33b (dissertation CRG)
  !>
  !> @param[in]     atoms    : Atoms type, see types.f90
  !> @param[in]     nobd     : Number of occupied bands per k-point and spin
  !> @param[in]     rbas1    : Large components of radial solution, its energy derivative and u_LO
  !> @param[in]     rbas2    : Small components of radial solution, its energy derivative and u_LO
  !> @param[in]     we       : Unperturbed occupation number
  !> @param[in]     acof     : (Any) large matching coefficient related to u
  !> @param[in]     bcof     : (Any) large matching coefficient related to udot
  !> @param[in]     ccof     : (Any) large matching coefficient related to uLO
  !> @param[in]     acof2cjg : (Any) large matching coefficient related to u which is to be conjugated in the routine
  !> @param[in]     bcof2cjg : (Any) large matching coefficient related to udot which is to be conjugated in the routine
  !> @param[in]     ccof2cjg : (Any) large matching coefficient related to uLO which is to be conjugated in the routine
  !> @param[in]     ilo2p    : Mapping array giving the p value for given number of LO and itype
  !> @param[inout]  uu       : Spherical (l = 0) product of acof and acof
  !> @param[inout]  du       : Spherical (l = 0) product of bcof and acof
  !> @param[inout]  ud       : Spherical (l = 0) product of acof and bcof
  !> @param[inout]  dd       : Spherical (l = 0) product of bcof and bcof
  !> @param[inout]  aclo     : LO related, to be done!
  !> @param[inout]  bclo     : LO related, to be done!
  !> @param[inout]  cclo     : LO related, to be done!
  !> @param[inout]  acnmt    : LO related, to be done!
  !> @param[inout]  bcnmt    : LO related, to be done!
  !> @param[inout]  ccnmt    : LO related, to be done!
  !> @param[inout]  uunmt    : Non-spehrical (l> 0) product of acof and acof
  !> @param[inout]  udnmt    : Non-spehrical (l> 0) product of bcof and acof
  !> @param[inout]  dunmt    : Non-spehrical (l> 0) product of acof and bcof
  !> @param[inout]  ddnmt    : Non-spehrical (l> 0) product of bcof and bcof
  !---------------------------------------------------------------------------------------------------------------------------------
  subroutine calcVarZContrRho1MT( acof, bcof, ccof, acof2cjg, bcof2cjg, ccof2cjg, uu, du, dd, ud, aclo, bclo, cclo, uunmt, udnmt, &
      & dunmt, ddnmt, acnmt, bcnmt, ccnmt, atoms, nobd, rbas1, rbas2, we, ilo2p )

    use m_types, only : t_atoms

    implicit none

    ! Type parameter
    type(t_atoms), intent(in)    :: atoms

    ! Scalar parameter
    integer,       intent(in)    :: nobd

    ! Array parameter (inout because of the sum over k-points)
    real,          intent(in)    :: rbas1(:,:,0:,:,:)
    real,          intent(in)    :: rbas2(:,:,0:,:,:)
    real,          intent(in)    :: we(:)
    complex,       intent(in)    :: acof(:,0:,:)
    complex,       intent(in)    :: bcof(:,0:,:)
    complex,       intent(in)    :: ccof(-atoms%llod:,:,:,:)
    complex,       intent(in)    :: acof2cjg(:,0:,:)
    complex,       intent(in)    :: bcof2cjg(:,0:,:)
    complex,       intent(in)    :: ccof2cjg(-atoms%llod:,:,:,:)
    integer,       intent(in)    :: ilo2p(:, :)
    complex,       intent(inout) :: uu(:,:)
    complex,       intent(inout) :: du(:,:)
    complex,       intent(inout) :: ud(:,:)
    complex,       intent(inout) :: dd(:,:)
    complex,       intent(inout) :: aclo(:,:)
    complex,       intent(inout) :: bclo(:,:)
    complex,       intent(inout) :: cclo(:,:,:)
    complex,       intent(inout) :: acnmt(:,:,:,:)
    complex,       intent(inout) :: bcnmt(:,:,:,:)
    complex,       intent(inout) :: ccnmt(:,:,:,:)
    complex,       intent(inout) :: uunmt(:, :,:,:)
    complex,       intent(inout) :: udnmt(:, :,:,:)
    complex,       intent(inout) :: dunmt(:, :,:,:)
    complex,       intent(inout) :: ddnmt(:, :,:,:)

    if (.FALSE.) then
      open(109,file='000_coeffs',form='FORMATTED',position='append',action='WRITE',status='UNKNOWN')
    end if

    ! Spherical part of density without LOs
    call rhomt( atoms, we, nobd, acof2cjg, bcof2cjg, acof, bcof, uu, dd, du, ud )

    if (.false.) then
      ! Spherical part of the LOs
      call rhomtlo( atoms, nobd, we, acof, bcof, ccof, acof2cjg, bcof2cjg, ccof2cjg, aclo, bclo, cclo )
    end if

    ! Non spherical part of LAPWs
    call rhonmt( atoms, we, nobd, acof, bcof, acof2cjg, bcof2cjg, uunmt, ddnmt, udnmt, dunmt )

    if (.false.) then
      ! Non spherical part of LOs of 7.34d
      call rhonmtlo( atoms, nobd, we, acof, bcof, ccof, acof2cjg, bcof2cjg, ccof2cjg, acnmt, bcnmt, ccnmt )
    end if

    if (.FALSE.) then
      close(109)
    end if

  end subroutine calcVarZContrRho1MT

  !---------------------------------------------------------------------------------------------------------------------------------
  !> @author
  !> Christian-Roman Gerhorst, Forschungszentrum Jülich: IAS1 / PGI1
  !>
  !> @brief
  !> Calculates the valence contribution of the first variation of the density response in the muffin-tin which are k-dependent
  !>
  !> @details
  !> See k-dependent part in 7.31b (7.33b) including factor 4, Gaunt coefficient and sum over m' and m'' (dissertation CRG).
  !>
  !> @param[in]     atoms   : Atoms type, see types.f90
  !> @param[in]     dimens  : Dimension type, see types.f90.
  !> @param[in]     sym     : Symmetries type, see types.f90.
  !> @param[in]     cell    : Unit cell type, see types.f90.
  !> @param[in]     kpts    : K-points type, see types.f90.
  !> @param[in]     usdus   : Type containing quantities consisting of the radial solutions, see types.f90.
  !> @param[in]     results : Results type, see types.f90.
  !> @param[in]     ikpt    : Index of k-point in k-point set
  !> @param[in]     ikpq    : Index of k + q (k' backfolded) in k-point set
  !> @param[in]     idispA  : Index of displaced atom
  !> @param[in]     nv      : Number of LAPW G-basis vectors for given k-point.
  !> @param[in]     ilst    : For various k-points G-basis vectors occur more than once, thus they are only stored once in juPhon.
  !> @param[in]     GbasVec : G-basis vectors
  !> @param[in]     nobd    : Number of occupied bands per k-point and spin
  !> @param[in]     z0      : Kohn-Sham eigenvectors.
  !> @param[in]     z1      : Linear variation of the Kohn\--Sham eigenvectors
  !> @param[in]     kveclo  : Basis G-vectors of local orbitals.
  !> @param[in]     rbas1   : Large components of radial solution, its energy derivative and u_LO
  !> @param[in]     rbas2   : Small components of radial solution, its energy derivative and u_LO
  !> @param[in]     ilo2p   : mapping array giving the p value for given number of LO and itype
  !> @param[inout]  uu      : Spherical (l = 0) product of acof and acof
  !> @param[inout]  du      : Spherical (l = 0) product of bcof and acof
  !> @param[inout]  ud      : Spherical (l = 0) product of acof and bcof
  !> @param[inout]  dd      : Spherical (l = 0) product of bcof and bcof
  !> @param[inout]  aclo    : LO related, to be done!
  !> @param[inout]  bclo    : LO related, to be done!
  !> @param[inout]  cclo    : LO related, to be done!
  !> @param[inout]  uunmt   : Non-spehrical (l> 0) product of acof and acof
  !> @param[inout]  udnmt   : Non-spehrical (l> 0) product of bcof and acof
  !> @param[inout]  dunmt   : Non-spehrical (l> 0) product of acof and bcof
  !> @param[inout]  ddnmt   : Non-spehrical (l> 0) product of bcof and bcof
  !> @param[inout]  acnmt   : LO related, to be done!
  !> @param[inout]  bcnmt   : LO related, to be done!
  !> @param[inout]  ccnmt   : LO related, to be done!
  !---------------------------------------------------------------------------------------------------------------------------------
  subroutine calcKdepValRho1MT( fmpi, oneD, atoms, input, sym, cell, kpts, usdus, results, ikpt, ikpq, idispA, nv, ilst, GbasVec, nobd,&
      & uu, du, ud, dd, aclo, bclo, cclo, uunmt, udnmt, dunmt, ddnmt, acnmt, bcnmt, ccnmt, z0, z1, kveclo, rbas1, rbas2, ilo2p, &
      & uu2, du2, ud2, dd2, aclo2, bclo2, cclo2, uunmt2, udnmt2, dunmt2, ddnmt2, acnmt2, bcnmt2, ccnmt2 )

    use m_types
    use m_abcof, only : abcof

    implicit none

    ! Type parameter
    type(t_mpi),                  intent(in)    :: fmpi
    type(t_oneD),                  intent(in)    :: oneD
    type(t_atoms),                  intent(in)    :: atoms
    type(t_input),                  intent(in)    :: input
    type(t_sym),                    intent(in)    :: sym
    type(t_cell),                   intent(in)    :: cell
    type(t_kpts),                   intent(in)    :: kpts
    type(t_usdus),                  intent(in)    :: usdus
    type(t_results),                intent(in)    :: results

    ! Scalar parameter
    integer,                        intent(in)    :: ikpt
    integer,                        intent(in)    :: ikpq
    integer,                        intent(in)    :: idispA

    ! Array parameter
    integer,                        intent(in)    :: nv(:, :)
    integer,                        intent(in)    :: ilst(:, :, :)
    integer,                        intent(in)    :: GbasVec(:, :)
    integer,                        intent(in)    :: nobd(:)
    complex,                       intent(in)    :: z0(:, :)
    complex,                        intent(in)    :: z1(:, :, :)
    integer,                        intent(in)    :: kveclo(:,:)
    real,                           intent(in)    :: rbas1(:,:,:,:,:)
    real,                           intent(in)    :: rbas2(:,:,:,:,:)
    integer,                        intent(in)    :: ilo2p(:, :)
    complex,                        intent(inout) :: uu(0:,:,:)
    complex,                        intent(inout) :: du(0:,:,:)
    complex,                        intent(inout) :: ud(0:,:,:)
    complex,                        intent(inout) :: dd(0:,:,:)
    complex,                        intent(inout) :: aclo(:,:,:)
    complex,                        intent(inout) :: bclo(:,:,:)
    complex,                        intent(inout) :: cclo(:,:,:,:)
    complex,                        intent(inout) :: uunmt(0:,0:,:,:,:)
    complex,                        intent(inout) :: udnmt(0:,0:,:,:,:)
    complex,                        intent(inout) :: dunmt(0:,0:,:,:,:)
    complex,                        intent(inout) :: ddnmt(0:,0:,:,:,:)
    complex,                        intent(inout) :: acnmt(0:,:,:,:, :)
    complex,                        intent(inout) :: bcnmt(0:,:,:,:, :)
    complex,                        intent(inout) :: ccnmt(:,:,:,:, :)
    complex, optional,              intent(inout) :: uu2(0:,:,:)
    complex, optional,              intent(inout) :: du2(0:,:,:)
    complex, optional,              intent(inout) :: ud2(0:,:,:)
    complex, optional,              intent(inout) :: dd2(0:,:,:)
    complex, optional,              intent(inout) :: aclo2(:,:,:)
    complex, optional,              intent(inout) :: bclo2(:,:,:)
    complex, optional,              intent(inout) :: cclo2(:,:,:,:)
    complex, optional,              intent(inout) :: uunmt2(0:,0:,:,:,:)
    complex, optional,              intent(inout) :: udnmt2(0:,0:,:,:,:)
    complex, optional,              intent(inout) :: dunmt2(0:,0:,:,:,:)
    complex, optional,              intent(inout) :: ddnmt2(0:,0:,:,:,:)
    complex, optional,              intent(inout) :: acnmt2(0:,:,:,:, :)
    complex, optional,              intent(inout) :: bcnmt2(0:,:,:,:, :)
    complex, optional,              intent(inout) :: ccnmt2(:,:,:,:, :)

    ! Local type variables !todo beware maybe not take them from fleur_init might be dangerous
    type(t_noco)                                  :: noco
    type(t_nococonv)                                  :: nococonv
    type(od_inp)                                  :: odi
    type(od_sym)                                  :: ods
    TYPE(t_lapw) :: lapw
    TYPE (t_mat) :: zMat1, zMatikpG, zMat

    ! Local scalar variable
    integer                                       :: nmat
    integer                                       :: idir
    integer                                       :: iband, lm
    integer                                       :: iBas
    integer                                       :: iatom
    integer                                       :: itype
    integer                                       :: ieqat
    integer :: nk

    ! Local array variables
    real                                          :: Gpkext(3)
    complex,           allocatable                :: acofz1(:, :, :, :)
    complex,           allocatable                :: bcofz1(:, :, :, :)
    complex,           allocatable                :: ccofz1(:, :, :, :, :)
    complex,           allocatable                :: acofikpG(:, :, :, :)
    complex,           allocatable                :: bcofikpG(:, :, :, :)
    complex,           allocatable                :: ccofikpG(:, :, :, :, :)
    complex,           allocatable                :: acofSummed(:, :, :, :)
    complex,           allocatable                :: bcofSummed(:, :, :, :)
    complex,           allocatable                :: ccofSummed(:, :, :, :, :)
    complex,           allocatable                :: acof(:, :, :)
    complex,           allocatable                :: bcof(:, :, :)
    complex,           allocatable                :: ccof(:, :, :, :)
    complex,           allocatable                :: zBar(:, :)
    integer,           allocatable                :: ngoprI(:)



    ! Initialization
    allocate( acofz1(nobd(ikpt), 0:atoms%lmaxd*(atoms%lmaxd+2), atoms%nat, 3), bcofz1(nobd(ikpt), 0:atoms%lmaxd*(atoms%lmaxd+2), atoms%nat, 3), &
      &ccofz1(-atoms%llod:atoms%llod, nobd(ikpt), atoms%nlod, atoms%nat, 3) )
    allocate( acofikpG(nobd(ikpt), 0:atoms%lmaxd*(atoms%lmaxd+2), atoms%nat, 3), bcofikpG(nobd(ikpt), 0:atoms%lmaxd*(atoms%lmaxd+2), atoms%nat, 3), &
      &ccofikpG(-atoms%llod:atoms%llod, nobd(ikpt), atoms%nlod, atoms%nat, 3) )
    acofz1(:, :, :, :)      = cmplx(0., 0.)
    bcofz1(:, :, :, :)      = cmplx(0., 0.)
    ccofz1(:, :, :, :, :)   = cmplx(0., 0.)
    acofikpG(:, :, :, :)    = cmplx(0., 0.)
    bcofikpG(:, :, :, :)    = cmplx(0., 0.)
    ccofikpG(:, :, :, :, :) = cmplx(0., 0.)

    ! We have to avoid the rotation of the local coordinate systems
    allocate(ngoprI(atoms%nat))
    ngoprI(:) = 1

    ! We do not want to have noco
    allocate( nococonv%alph(atoms%ntype), nococonv%beta(atoms%ntype) )
    nococonv%alph(:) = 0
    nococonv%beta(:) = 0

    ! Set noco stuff to false.
    noco%l_soc  = .FALSE.
    noco%l_noco = .FALSE.
    noco%l_ss   = .FALSE.
    ALLOCATE(noco%l_constrained(atoms%ntype))
    ALLOCATE(noco%l_unrestrictMT(atoms%ntype))
    ALLOCATE(noco%l_spinoffd_ldau(atoms%ntype))
    noco%l_constrained(atoms%ntype)   = .FALSE.
    noco%l_unrestrictMT(atoms%ntype)  = .FALSE.
    noco%l_spinoffd_ldau(atoms%ntype) = .FALSE.

    allocate(zBar(SIZE(z0(:,1)), nobd(ikpt)))
    zBar(:, :) = cmplx(0., 0.)

    ! Calculate abcofs with z1 and with iG z0, see 7.30e
    do idir = 1, 3
      nmat = nv(1, ikpq) + atoms%nlotot
      !call abcof( atoms%lmaxd, atoms%ntype, input%neig, nobd(ikpt), atoms%nat, sym%nop, MAXVAL(nv), input%jspins, atoms%lmaxd*(atoms%lmaxd+2),     &
        !& SIZE(z0(:,1)), atoms%llod, atoms%nlod, atoms%nlotot, sym%invtab, atoms%ntype, sym%mrot, ngoprI, atoms%taual, atoms%neq, &
        !& atoms%lmax, atoms%rmt, cell%omtil, cell%bmat, cell%bbmat, kpts%bk(:, ikpq), GbasVec(1, ilst(:nv(1, ikpq), ikpq, 1)),     &
        !& GbasVec(2, ilst(:nv(1, ikpq), ikpq, 1)), GbasVec(3, ilst(:nv(1, ikpq), ikpq, 1)), nv(:, ikpq),  nmat, nobd(ikpt),        &
        !& z1(:, :, idir), usdus%us(:, :, 1), usdus%dus(:, :, 1), usdus%uds(:, :, 1), usdus%duds(:, :, 1), usdus%ddn(:, :, 1),      &
        !& sym%invsat, sym%invsatnr, usdus%ulos(:, :, 1), usdus%uulon(:, :, 1), usdus%dulon(:, :, 1),  usdus%dulos(:, :, 1),      &
        !& atoms%llo, atoms%nlo, atoms%l_dulo, atoms%lapw_l, noco%l_noco, noco%l_ss, 1, nococonv%alph, nococonv%beta, nococonv%qss,             &
        !& kveclo(:, ikpq), odi, ods, acofz1(:, :, :, idir), bcofz1(:, :, :, idir), ccofz1(:, :, :, :, idir) )

        nk=fmpi%k_list(ikpq)
        CALL lapw%init(input, noco, nococonv, kpts, atoms, sym, nk, cell, fmpi)
        CALL zMat1%init(.FALSE., nv(1, ikpq) + atoms%nlotot, nobd(ikpt))
        zMat1%data_c(:, :) = z1(:nv(1, ikpq) + atoms%nlotot, :nobd(ikpt), idir)
        CALL abcof(input, atoms, sym, cell, lapw, nobd(ikpt), usdus, noco, nococonv, 1, oneD, &
                 & acofz1(:, :, :, idir), bcofz1(:, :, :, idir), &
                 & ccofz1(-atoms%llod:, :, :, :, idir), zMat1)

      ! Calculate iG abcof, see 7.30a and 7.30f
      zBar(:, :) = cmplx(0., 0.)
      do iband = 1, nobd(ikpt)
        do iBas = 1, nv(1, ikpt)
          Gpkext(1:3) = matmul( cell%bmat(1:3, 1:3), GbasVec(1:3, ilst(iBas, ikpt, 1)) + kpts%bk(1:3, ikpt))
          zBar(iBas, iband) = ImagUnit * Gpkext(idir) * z0(iBas, iband)
        end do
        if (.false.) then
          do iBas = nv(1, ikpt) + 1, nv(1, ikpt) + atoms%nlotot
            Gpkext(1:3) = matmul( cell%bmat(1:3, 1:3), GbasVec(1:3, ilst(kveclo(iBas - nv(1, ikpt), ikpt), ikpt, 1)) &
              & + kpts%bk(1:3, ikpt) )
            zBar(iBas, iband) = ImagUnit * Gpkext(idir) * z0(iBas, iband)
          end do
        end if
      end do

      ! Calculating acofBar, bcofBar, ccofBar using zBar
      ! When rewriting abcof dimensions of acof, bcof, and ccof and be made even smaller
      ! the first neigd has to be neigd because the z0 read from fleur has to be neigd to flexibly use them within the Sternheimer
      ! equation
      nmat = nv(1, ikpt) + atoms%nlotot
      !call abcof ( atoms%lmaxd, atoms%ntype, nobd(ikpt), nobd(ikpt), atoms%nat, sym%nop, MAXVAL(nv), input%jspins, atoms%lmaxd*(atoms%lmaxd+2),      &
        !& SIZE(z0(:,1)), atoms%llod, atoms%nlod, atoms%nlotot, sym%invtab, atoms%ntype, sym%mrot, ngoprI, atoms%taual, atoms%neq, &
        !& atoms%lmax, atoms%rmt, cell%omtil, cell%bmat, cell%bbmat, kpts%bk(:, ikpt), GbasVec(1, ilst(:nv(1, ikpt), ikpt, 1)),     &
        !& GbasVec(2, ilst(:nv(1, ikpt), ikpt, 1)), GbasVec(3, ilst(:nv(1, ikpt), ikpt, 1)), nv(:, ikpt),  nmat, nobd(ikpt),        &
        !& zBar(:, :), usdus%us(:, :, 1), usdus%dus(:, :, 1), usdus%uds, usdus%duds(:, :, 1), usdus%ddn(:, :, 1), sym%invsat,     &
        !& sym%invsatnr, usdus%ulos(:, :, 1), usdus%uulon(:, :, 1), usdus%dulon(:, :, 1),  usdus%dulos(:, :, 1), atoms%llo,         &
        !& atoms%nlo, atoms%l_dulo, atoms%lapw_l, noco%l_noco, noco%l_ss, 1, nococonv%alph, nococonv%beta, nococonv%qss, kveclo(:, ikpt), odi,  &
        !& ods, acofikpG(:, :, :, idir), bcofikpG(:, :, :, idir), ccofikpG(:, :, :, :, idir) )
        nk=fmpi%k_list(ikpt)
        CALL lapw%init(input, noco, nococonv, kpts, atoms, sym, nk, cell, fmpi)
        CALL zMatikpG%init(.FALSE., nv(1, ikpt) + atoms%nlotot, nobd(ikpt))
        zMatikpG%data_c(:, :) = zBar(:nv(1, ikpt) + atoms%nlotot, :nobd(ikpt))
        CALL abcof(input, atoms, sym, cell, lapw, nobd(ikpt), usdus, noco, nococonv, 1, oneD, &
                 & acofikpG(:, :, :, idir), bcofikpG(:, :, :, idir), &
                 & ccofikpG(-atoms%llod:, :, :, :, idir), zMatikpG)
    end do

    ! Calculate classical abcof for multiplying it with k-vector, see 7.30f
    allocate( acof(nobd(ikpt), 0:atoms%lmaxd*(atoms%lmaxd+2), atoms%nat), bcof(nobd(ikpt), 0:atoms%lmaxd*(atoms%lmaxd+2), atoms%nat), &
            &ccof(-atoms%llod:atoms%llod, nobd(ikpt), atoms%nlod, atoms%nat) )
    acof(:, :, :) = cmplx(0., 0.)
    bcof(:, :, :) = cmplx(0., 0.)
    ccof(:, :, :, :) = cmplx(0., 0.)

    nmat = nv(1, ikpt) + atoms%nlotot
    !call abcof ( atoms%lmaxd, atoms%ntype, input%neig, nobd(ikpt), atoms%nat, sym%nop, MAXVAL(nv), input%jspins, atoms%lmaxd*(atoms%lmaxd+2),      &
    !     & SIZE(z0(:,1)), atoms%llod, atoms%nlod, atoms%nlotot, sym%invtab, atoms%ntype, sym%mrot, ngoprI, atoms%taual, atoms%neq,&
    !     & atoms%lmax, atoms%rmt, cell%omtil, cell%bmat, cell%bbmat, kpts%bk(:, ikpt), GbasVec(1, ilst(:nv(1, ikpt), ikpt, 1)),    &
    !     & GbasVec(2, ilst(:nv(1, ikpt), ikpt, 1)), GbasVec(3, ilst(:nv(1, ikpt), ikpt, 1)), nv(:, ikpt),  nmat, nobd(ikpt),       &
    !     & z0(:, :), usdus%us(:, :, 1), usdus%dus(:, :, 1), usdus%uds(:, :, 1), usdus%duds(:, :, 1), usdus%ddn(:, :, 1),           &
    !     & sym%invsat, sym%invsatnr, usdus%ulos(:, :, 1), usdus%uulon(:, :, 1), usdus%dulon(:, :, 1),  usdus%dulos(:, :, 1),     &
    !     & atoms%llo, atoms%nlo, atoms%l_dulo, atoms%lapw_l, noco%l_noco, noco%l_ss, 1, nococonv%alph, nococonv%beta, nococonv%qss,            &
    !     & kveclo(:, ikpt), odi, ods, acof, bcof, ccof )
    nk=fmpi%k_list(ikpt)
    CALL lapw%init(input, noco, nococonv, kpts, atoms, sym, nk, cell, fmpi)
    CALL zMat%init(.FALSE., nv(1, ikpt) + atoms%nlotot, nobd(ikpt))
    zMat%data_c(:, :) = z0(:nv(1, ikpt) + atoms%nlotot, :nobd(ikpt))
    CALL abcof(input, atoms, sym, cell, lapw, nobd(ikpt), usdus, noco, nococonv, 1, oneD, &
             & acof, bcof, &
             & ccof, zMat)

    allocate( acofSummed(nobd(ikpt), 0:atoms%lmaxd*(atoms%lmaxd+2), atoms%nat, 3), bcofSummed(nobd(ikpt), 0:atoms%lmaxd*(atoms%lmaxd+2), atoms%nat, 3), &
      &ccofSummed(-atoms%llod:atoms%llod, nobd(ikpt), atoms%nlod, atoms%nat, 3) )
    acofSummed(:, :, :, :) = cmplx(0., 0.)
    bcofSummed(:, :, :, :) = cmplx(0., 0.)
    ccofSummed(:, :, :, :, :) = cmplx(0., 0.)

    ! We care about this when not dealing with insulators any more.
    !todo what is if nobd(ikpt /= ikpq. we might run into problems at summing up the matching coefficients!
    if (nobd(ikpq) /= nobd(ikpt)) then
      write(*, *) 'Check following lines!!!'
    end if

    ! sum up matching coefficients from expansion coefficient variation and from basis variation, see 7.30g
    do idir = 1, 3
      iatom = 0
      do itype = 1, atoms%ntype
        do ieqat = 1, atoms%neq(itype)
          iatom = iatom + 1
          ! todo TEST: for the perturbuation test to work we have to comment out the if branch
          ! todo TEST: for the test where contributions cancel away for occupied occupied band combinations we have to switch out
          ! the basis set correction
          do iband = 1, nobd(ikpt)
            if ( iatom == idispA ) then

              acofSummed(iband, 0:atoms%lmaxd*(atoms%lmaxd+2), iatom, idir) = acofz1(iband, 0:atoms%lmaxd*(atoms%lmaxd+2), iatom, idir) &
                & + acofikpG(iband, 0:atoms%lmaxd*(atoms%lmaxd+2), iatom, idir)
              bcofSummed(iband, 0:atoms%lmaxd*(atoms%lmaxd+2), iatom, idir) = bcofz1(iband, 0:atoms%lmaxd*(atoms%lmaxd+2), iatom, idir) &
                & + bcofikpG(iband, 0:atoms%lmaxd*(atoms%lmaxd+2), iatom, idir)
              if (.false.) then
                if ( atoms%nlo(itype) > 0 ) then
                  ccofSummed(-atoms%llod:atoms%llod, iband, :atoms%nlod, iatom, idir) = &
                    & ccofz1(-atoms%llod:atoms%llod, iband, :atoms%nlod, iatom, idir) &
                    &+ ccofikpG(-atoms%llod:atoms%llod, iband, :atoms%nlod, iatom, idir)
                end if
              end if
            else
              acofSummed(iband, 0:atoms%lmaxd*(atoms%lmaxd+2), iatom, idir) = acofz1(iband, 0:atoms%lmaxd*(atoms%lmaxd+2), iatom, idir)
              bcofSummed(iband, 0:atoms%lmaxd*(atoms%lmaxd+2), iatom, idir) = bcofz1(iband, 0:atoms%lmaxd*(atoms%lmaxd+2), iatom, idir)
              if (.false.) then
                if ( atoms%nlo(itype) > 0 ) then
                  ccofSummed(-atoms%llod:atoms%llod, iband, :atoms%nlod, iatom, idir) = &
                    & ccofz1(-atoms%llod:atoms%llod, iband, :atoms%nlod, iatom, idir)
                end if
              end if
            end if
          end do
        end do
      end do
    end do

    if (.FALSE.) then
      if (ikpt.eq.1) then
        open(109,file='000_matchsum',form='FORMATTED',position='append',action='WRITE',status='REPLACE')
      else
         open(109,file='000_matchsum',form='FORMATTED',position='append',action='WRITE',status='UNKNOWN')
      end if
      do iband = 1, nobd(ikpt)
        do lm = 0, atoms%lmaxd*(atoms%lmaxd+2)
          write(109,*) ikpq, iband, 1, lm+1
          write(109,*) real(acofSummed(iband,lm,1,1)), aimag(acofSummed(iband,lm,1,1))
          write(109,*) ikpq, iband, 2, lm+1
          write(109,*) real(bcofSummed(iband,lm,1,1)), aimag(bcofSummed(iband,lm,1,1))
        end do
      end do
      close(109)
    end if

    ! Use adjusted Fleur routine for calculating k-dependent parts of the density variation
    ! Factor 2 because of spin degeneracy
    do idir = 1, 3
      call calcVarZContrRho1MT( acofSummed(:, :, :, idir), bcofSummed(:, :, :, idir), ccofSummed(:, :, :, :, idir), acof, bcof,    &
        &ccof, uu(:, :, idir), du(:, :, idir), dd(:, :, idir), ud(:, :, idir), aclo(:, :, idir), bclo(:, :, idir),                 &
        &cclo(:, :, :, idir), uunmt(:, :, :, :, idir), udnmt(:, :, :, :, idir), dunmt(:, :, :, :, idir), ddnmt(:, :, :, :, idir),  &
        &acnmt(:, :, :, :, idir), bcnmt(:, :, :, :, idir), ccnmt(:, :, :, :, idir), atoms, nobd(ikpt), rbas1, rbas2,               &
        &2 * results%w_iks(:, ikpt, 1),  ilo2p )
      if (present(uu2)) then
        call calcVarZContrRho1MT( acofikpG(:, :, :, idir), bcofikpG(:, :, :, idir), ccofikpG(:, :, :, :, idir), acof, bcof,              &
          &ccof, uu2(:, :, idir), du2(:, :, idir), dd2(:, :, idir), ud2(:, :, idir), aclo2(:, :, idir), bclo2(:, :, idir),               &
          &cclo2(:, :, :, idir), uunmt2(:, :, :, :, idir), udnmt2(:, :, :, :, idir), dunmt2(:, :, :, :, idir), ddnmt2(:, :, :, :, idir), &
          &acnmt2(:, :, :, :, idir), bcnmt2(:, :, :, :, idir), ccnmt2(:, :, :, :, idir), atoms, nobd(ikpt), rbas1, rbas2,                &
          &2 * results%w_iks(:, ikpt, 1),  ilo2p )
      end if
    end do

    deallocate( acofz1, bcofz1, ccofz1, acofikpG, bcofikpG, ccofikpG )

  end subroutine calcKdepValRho1MT

  !---------------------------------------------------------------------------------------------------------------------------------
  !> @author
  !> Christian-Roman Gerhorst, Forschungszentrum Jülich: IAS1 / PGI1
  !>
  !> @brief
  !> This subroutine multiplies the radial solutions to the MT density change contribution coming from the variation of the
  !> expansion coefficients
  !>
  !> @details
  !> See also 7.31b (7.33b) (dissertation CRG)
  !>
  !> @param[in]   atoms : Atoms type, see types.f90
  !> @param[in]   aclo  : LO related, to be done!
  !> @param[in]   bclo  : LO related, to be done!
  !> @param[in]   cclo  : LO related, to be done!
  !> @param[in]   acnmt : LO related, to be done!
  !> @param[in]   bcnmt : LO related, to be done!
  !> @param[in]   ccnmt : LO related, to be done!
  !> @param[in]   rbas1 : Large components of radial solution, its energy derivative and u_LO
  !> @param[in]   rbas2 : Small components of radial solution, its energy derivative and u_LO
  !> @param[in]   uu    : Spherical (l = 0) product of acof and acof
  !> @param[in]   du    : Spherical (l = 0) product of bcof and acof
  !> @param[in]   ud    : Spherical (l = 0) product of acof and bcof
  !> @param[in]   dd    : Spherical (l = 0) product of bcof and bcof
  !> @param[in]   uunmt : Non-spehrical (l> 0) product of acof and acof
  !> @param[in]   udnmt : Non-spehrical (l> 0) product of bcof and acof
  !> @param[in]   dunmt : Non-spehrical (l> 0) product of acof and bcof
  !> @param[in]   ddnmt : Non-spehrical (l> 0) product of bcof and bcof
  !> @param[in]   ilo2p : mapping array giving the p value for given number of LO and itype
  !> @param[out]  rho   : Spherical harmonic expansion coefficients of the full linear variation of the muffin-tin density
  !>
  !> @note
  !> adjusted from rhomt and rhosphnlo (from cdnmt actually!)from new fleur new at density change
  !> multiply radial solutions to the contribution coming from the variation of the expansion coefficients z to the density change
  !>
  !> @todo
  !> LO related
  !---------------------------------------------------------------------------------------------------------------------------------
  subroutine multRadSolVzcRho1MT( atoms, aclo, bclo, cclo, acnmt, bcnmt, ccnmt, rbas1, rbas2, uu, du, ud, dd, uunmt, udnmt,  &
      & dunmt, ddnmt, ilo2p, rho )

    use m_types, only : t_atoms
    use m_constants, only : sfp_const
    use m_juDFT_stop, only : juDFT_error

    implicit none

    ! Type parameter
    type(t_atoms), intent(in)    :: atoms

    ! Array parameter
    complex,       intent(in)    :: aclo(:,:,:)
    complex,       intent(in)    :: bclo(:,:,:)
    complex,       intent(in)    :: cclo(:,:,:,:)
    complex,       intent(in)    :: acnmt(0:,:,:,:,:)
    complex,       intent(in)    :: bcnmt(0:,:,:,:,:)
    complex,       intent(in)    :: ccnmt(:,:,:,:,:)
    real,          intent(in)    :: rbas1(:,:,0:,:,:)
    real,          intent(in)    :: rbas2(:,:,0:,:,:)
    complex,       intent(in)    :: uu(0:,:,:)
    complex,       intent(in)    :: du(0:,:,:)
    complex,       intent(in)    :: ud(0:,:,:)
    complex,       intent(in)    :: dd(0:,:,:)
    complex,       intent(in)    :: uunmt(0:,0:,:,:,:)
    complex,       intent(in)    :: udnmt(0:,0:,:,:,:)
    complex,       intent(in)    :: dunmt(0:,0:,:,:,:)
    complex,       intent(in)    :: ddnmt(0:,0:,:,:,:)
    integer,       intent(in)    :: ilo2p(:, :)
    complex,       intent(out)   :: rho(:,:,:,:)

    ! Local Scalars
    integer                      :: itype
    integer                      :: na
    integer                      :: l
    integer                      :: m
    integer                      :: lm
    integer                      :: lp
    integer                      :: imesh
    integer                      :: ieqat
    integer                      :: lo
    integer                      :: lop
    integer                      :: lv
    integer                      :: mv
    integer                      :: lmv
    integer                      :: idir
    real                         :: c_1
    real                         :: c_2
    complex                      :: s

    rho(:, :, :, :) = cmplx(0.0, 0.0)

    do idir = 1, 3
      na = 0
      do itype = 1, atoms%ntype
        do ieqat = 1, atoms%neq(itype)
          na = na + 1

          ! This is assigned to the l = 0, m = 0 spherical component of the density
          lm = 1
          do l = 0, atoms%lmax(itype)
            do imesh = 1, atoms%jri(itype) ! beware that these are the correct radial solutions here
              s = uu(l, itype, idir) * ( rbas1(imesh, 1, l, itype, 1) * rbas1(imesh, 1, l, itype, 1)   &
                                     & + rbas2(imesh, 1, l, itype, 1) * rbas2(imesh, 1, l, itype, 1) ) &
              & + dd(l, itype, idir) * ( rbas1(imesh, 2, l, itype, 1) * rbas1(imesh, 2, l, itype, 1)   &
                                     & + rbas2(imesh, 2, l, itype, 1) * rbas2(imesh, 2, l, itype, 1) ) &
              & + du(l, itype, idir) * ( rbas1(imesh, 2, l, itype, 1) * rbas1(imesh, 1, l, itype, 1)   &
                                     & + rbas2(imesh, 2, l, itype, 1) * rbas2(imesh, 1, l, itype, 1) ) &
              & + ud(l, itype, idir) * ( rbas1(imesh, 1, l, itype, 1) * rbas1(imesh, 2, l, itype, 1)   &
                                     & + rbas2(imesh, 1, l, itype, 1) * rbas2(imesh, 2, l, itype, 1) )
              ! The factor 1/ sqrt 4pi comes from the Gaunt factor if the middle l and m is zero. Therefore it is correct
              rho(imesh, lm, na, idir) = rho(imesh, lm, na, idir) &
                                       & + s / atoms%rmsh(imesh, itype) / atoms%rmsh(imesh, itype) / (sfp_const)
            end do ! imesh
          end do ! l

          if (.false.) then
            c_2 = 1.0 / (sfp_const)
            !--->       add the contribution of the local orbitals and flapw - lo
            !--->       cross-terms to rho, qmtl. the latter are stored in
            !--->       qmtllo. initialize qmtllo
            !---> add the contribution of the local orbitals and flapw - lo cross-
            !---> terms to the spherical chargedensity inside the muffin tins.
            ! This is assigned to the l = 0, m = 0 spherical component of the density
            lm = 1
            do lo = 1, atoms%nlo(itype)
              l = atoms%llo(lo, itype)
              do imesh = 1, atoms%jri(itype)
                rho(imesh, lm, na, idir) = rho(imesh, lm, na, idir) + c_2 * (aclo(lo, na, idir) * ( rbas1(imesh, 1, l, itype, 1) * &
                  & rbas1(imesh, ilo2p(lo, itype), l, itype, 1) + rbas2(imesh, 1, l, itype, 1) &
                  & * rbas2(imesh, ilo2p(lo, itype), l, itype, 1) ) + bclo(lo, na, idir) * ( rbas1(imesh, 2, l, itype, 1) &
                  & * rbas1(imesh, ilo2p(lo, itype), l, itype, 1) + rbas2(imesh, 2, l, itype, 1) &
                  & * rbas2(imesh, ilo2p(lo, itype), l, itype, 1) ) ) / atoms%rmsh(imesh, itype) / atoms%rmsh(imesh, itype)
              end do ! imesh
              do lop = 1,atoms%nlo(itype)
                if ( atoms%llo(lop, itype) == l ) then
                  do imesh = 1, atoms%jri(itype)
                    rho(imesh, lm, na, idir) = rho(imesh, lm, na, idir) + c_2 * cclo(lop, lo, na, idir) *                          &
                      & ( rbas1(imesh, ilo2p(lop, itype), l, itype, 1) * rbas1(imesh, ilo2p(lo, itype), l, itype, 1)               &
                      & + rbas2(imesh, ilo2p(lop, itype), l, itype, 1) * rbas2(imesh, ilo2p(lo, itype), l, itype, 1) )             &
                      & / atoms%rmsh(imesh, itype) / atoms%rmsh(imesh, itype)
                  end do ! imesh
                end if
              end do ! lop
            end do ! lo
          end if

          !---> add the contribution of the local orbitals and flapw - lo cross-
          !---> terms to the non-spherical chargedensity inside the muffin tins.
          c_1 = 1.0
          do l = 1, atoms%lmax(itype)
            do m = -l, l
              !TODO LO: lm should be incremented by one to start from 1, attention, this should also be done for the acnmt,...
              lm = l * (l + 1) + m
              do lp = 0, atoms%lmax(itype)
                do lo = 1, atoms%nlo(itype)
                  call juDFT_error('1. LO contribution not allowed.', calledby='multRadSolVzcRho1MT')
                  do imesh = 1, atoms%jri(itype)
                    rho(imesh, lm, na, idir) = rho(imesh, lm, na, idir) + c_1 * (acnmt(lp, lo, lm, na, idir)                       &
                      & * ( rbas1(imesh, 1, lp, itype, 1) * rbas1(imesh, ilo2p(lo, itype), atoms%llo(lo, itype), itype, 1)         &
                      & + rbas2(imesh, 1, lp, itype, 1) * rbas2(imesh, ilo2p(lo, itype), atoms%llo(lo, itype), itype, 1) )         &
                      & + bcnmt(lp, lo, lm, na, idir) * (rbas1(imesh, 2, lp, itype, 1)                                             &
                      & * rbas1(imesh, ilo2p(lo, itype), atoms%llo(lo, itype), itype, 1)                                           &
                      & + rbas2(imesh, 2, lp, itype, 1)  * rbas2(imesh, ilo2p(lo, itype), atoms%llo(lo, itype), itype, 1) ) )      &
                      & / atoms%rmsh(imesh, itype) / atoms%rmsh(imesh, itype)
                  end do ! imesh
                end do ! lo
              end do ! lp
              do lo = 1,atoms%nlo(itype)
                call juDFT_error('2. LO contribution not allowed.', calledby='multRadSolVzcRho1MT')
                do lop = 1,atoms%nlo(itype)
                  do imesh = 1,atoms%jri(itype)
                    rho(imesh,lm, na, idir) = rho(imesh,lm, na, idir) + c_1 * ccnmt(lop, lo, lm, na, idir)                         &
                      & * ( rbas1(imesh, ilo2p(lop, itype), atoms%llo(lop, itype), itype, 1)                                       &
                      & * rbas1(imesh, ilo2p(lo, itype), atoms%llo(lo, itype), itype, 1) +                                         &
                      & rbas2(imesh, ilo2p(lop, itype), atoms%llo(lop, itype), itype, 1) *                                         &
                      & rbas2(imesh, ilo2p(lo, itype), atoms%llo(lo, itype), itype, 1) ) / atoms%rmsh(imesh, itype)                &
                      & / atoms%rmsh(imesh, itype)
                  end do ! imesh
                end do ! lop
              end do ! lo
            end do ! m
          end do ! l

          !--->       non-spherical components
          do lv = 1, atoms%lmax(itype)
            do mv = -lv, lv
              lmv = lv * (lv + 1) + mv + 1
              do l = 0,atoms%lmax(itype)
                do lp = 0,atoms%lmax(itype)
                  do imesh = 1,atoms%jri(itype)
                    s = uunmt(lp, l, lmv, na, idir) * ( rbas1(imesh, 1, lp, itype, 1) * rbas1(imesh, 1, l, itype, 1)   &
                                                    & + rbas2(imesh, 1, lp, itype, 1) * rbas2(imesh, 1, l, itype, 1) ) &
                    & + ddnmt(lp, l, lmv, na, idir) * ( rbas1(imesh, 2, lp, itype, 1) * rbas1(imesh, 2, l, itype, 1)   &
                                                    & + rbas2(imesh, 2, lp, itype, 1) * rbas2(imesh, 2, l, itype, 1) ) &
                    & + udnmt(lp, l, lmv, na, idir) * ( rbas1(imesh, 1, lp, itype, 1) * rbas1(imesh, 2, l, itype, 1)   &
                                                    & + rbas2(imesh, 1, lp, itype, 1) * rbas2(imesh, 2, l, itype, 1) ) &
                    & + dunmt(lp, l, lmv, na, idir) * ( rbas1(imesh, 2, lp, itype, 1) * rbas1(imesh, 1, l, itype, 1)   &
                                                    & + rbas2(imesh, 2, lp, itype, 1) * rbas2(imesh, 1, l, itype, 1) )
                    rho(imesh, lmv, na, idir) = rho(imesh, lmv, na, idir) + s / atoms%rmsh(imesh, itype) / atoms%rmsh(imesh, itype)
                  end do ! imesh
                end do ! lp
              end do ! l
            end do ! mv
          end do ! lv
        end do ! ieqat
      end do ! end of loop over atom types
    end do ! idir

  end subroutine multRadSolVzcRho1MT

  !---------------------------------------------------------------------------------------------------------------------------------
  !> @author
  !> Christian-Roman Gerhorst, Forschungszentrum Jülich: IAS1 / PGI1
  !>
  !> @brief
  !> Calculates the parameters for the gradient of the Gauss function pseudo-core density and the resulting interstitial core-tail
  !> correction for q = 0 that still have to be modulated with the respective q.
  !>
  !> @details
  !> See also 7.26b (dissertation CRG)
  !>
  !> @note
  !> This routine is based on a routine from FLEUR
  !>
  !>  @param[in]  atoms   : Atoms type, see types.f90
  !>  @param[in]  cell    : Unit cell type, see types.f90.
  !>  @param[in]  sym     : Symmetries type, see types.f90.
  !>  @param[in]  stars   : Stars type, see types.f90.
  !>  @param[in]  dimens  : Dimension type, see types.f90.
  !>  @param[in]  input   : Input type, see types.f90.
  !>  @param[in]  logUnit : Unit number for juPhon.log.
  !>  @param[in]  ngdp2km : Number of G-vectors for potentials and densities which are smaller than 2 kmax.
  !>  @param[in]  gdp     : G-vectors of potentials and densities
  !>  @param[out] acoff   : Uppercase A of Gauss curve pseudo core-density, see 7.26b (dissertation CRG)
  !>  @param[out] alpha   : Lowercase a of Gauss curve pseudo core-density, see 7.26b (dissertation CRG)
  !>  @param[out] qpwcG   : Plane-wave expansion coefficients of the FFT of the gradient of the pseudo core-density ( Gauss curve )
  !>
  !---------------------------------------------------------------------------------------------------------------------------------
!!!! CALCPSDENSMT removed !!!
  !---------------------------------------------------------------------------------------------------------------------------------
  !> @author
  !> Christian-Roman Gerhorst, Forschungszentrum Jülich: IAS1 / PGI1
  !>
  !> @brief
  !> Calculates the IR part of the core-tail correction and its Rayleigh expansion for the muffin-tins
  !>
  !> @details
  !> See also 7.26a and 7.25c (dissertation CRG)
  !>
  !> @param[in]  atoms     : Atoms type, see types.f90
  !> @param[in]  cell      : Unit cell type, see types.f90.
  !> @param[in]  ngpqdp2km : Number of G-vectors for potentials and densities which are smaller than 2 kmax.
  !> @param[in]  gpqdp     : Number of G-vectors for shifted G-set for q-point with index iqpt
  !> @param[in]  qbk       : Qpoint in internal coordinates as in qpts%bk(1:3, iqpt)
  !> @param[in]  rhoPsC    : Plane-wave expansion coefficients of the FFT of the gradient of the pseudo core-density ( Gauss
  !>                         curve )
  !> @param[out] rho1IRctC : Plane-wave expansion coefficients of the Interstitial core-tail correction as in 7.25
  !>                         (dissertation CRG)
  !> @param[out] rho1MTctC : Muffin-tin expansion coefficients of the Rayleigh expansion of rho1IRctC as 7.26a (dissertation CRG)
  !---------------------------------------------------------------------------------------------------------------------------------
  subroutine calcVarCTcorr( atoms, cell, ngpqdp2km, gpqdp, qbk, rhoPsC, rho1IRctC, rho1MTctC )

    use m_types, only : t_atoms, t_cell
    use m_sphbes, only : sphbes
    use m_ylm_old, only : ylm4

    implicit none

    ! Type arguments
    type(t_atoms),              intent(in)  :: atoms
    type(t_cell),               intent(in)  :: cell

    ! Scalar arguments
    integer,                    intent(in)  :: ngpqdp2km

    ! Array arguments
    integer,                    intent(in)  :: gpqdp(:, :)
    real,                       intent(in)  :: qbk(:)
    complex,                    intent(in)  :: rhoPsC(:, :)
    complex,                    intent(out) :: rho1IRctC(:, :, :)
    complex,                    intent(out) :: rho1MTctC(:, :, :, :, :)

    ! Local variables
    integer                                 :: iG
    complex                                 :: temp1
    complex                                 :: temp2
    integer                                 :: itype
    integer                                 :: iatom
    integer                                 :: imesh
    integer                                 :: ieqat
    integer                                 :: oqn_l
    integer                                 :: idir
    integer                                 :: lm_temp
    integer                                 :: lm
    integer                                 :: mqn_m
    integer                                 :: iDtype
    integer                                 :: iDeqat
    integer                                 :: iDatom
    real                                    :: normedG
    real                                    :: sbesArg

    ! Local arrays
    complex,       allocatable             :: ylm(:, :)
    complex,       allocatable             :: phFac(:, :)
    real,          allocatable             :: sbes(:, :, :, :)
    complex,       allocatable             :: fpiul(:)
    real                                   :: Gqext(3)
    complex,       allocatable             :: ylmtemp(:)
    real,          allocatable             :: sbestemp(:)

    ! Init local variables
    allocate( ylm(ngpqdp2km, ( atoms%lmaxd + 1 )**2), ylmtemp(( atoms%lmaxd + 1 )**2), &
            & sbes( atoms%jmtd, 0: atoms%lmaxd, ngpqdp2km, atoms%ntype ), &
            & sbestemp(0: atoms%lmaxd), &
            & phFac(ngpqdp2km, atoms%nat), &
            & fpiul(0:atoms%lmaxd) )
    ylm(:, :)        = cmplx(0., 0.)
    phFac(:, :)      = cmplx(0., 0.)
    sbes(:, :, :, :) = 0.
    Gqext(:)         = 0.

    ! Init dummy variables for intent(out)
    rho1IRctC = cmplx(0., 0.)
    rho1MTctC = cmplx(0., 0.)

    ! to fill it like this (no sequential memory) is the better compromise for later
    do iG = 1, ngpqdp2km
      Gqext(1:3) = matmul( cell%bmat(1:3, 1:3), real(gpqdp(1:3, iG) + qbk(1:3)) )
      call ylm4(atoms%lmaxd, Gqext(1:3), ylmtemp)
      ylm(iG, :) = ylmtemp
      normedG = norm2(Gqext)
      iatom = 0
      do itype = 1, atoms%ntype
        do imesh = 1, atoms%jri(itype)
          sbesArg = atoms%rmsh(imesh, itype) * normedG
          call sphbes(atoms%lmax(itype), sbesArg , sbestemp)
          sbes(imesh, :, iG, itype) = sbestemp
        end do
        do ieqat = 1, atoms%neq(itype)
          iatom = iatom + 1
          do idir = 1, 3
            rho1IRctC(iG, iatom, idir) = -ImagUnit * Gqext(idir) * rhoPsC(iG, iatom)
          end do ! idir
          phFac(iG, iatom) = exp(ImagUnit * tpi_const * dot_product(gpqdp(:, iG) + qbk(:), atoms%taual(:, iatom)))
        end do
      end do
    end do

    do oqn_l = 0, atoms%lmaxd
       fpiul(oqn_l) = fpi_const * ImagUnit**(oqn_l)
    enddo

    iDatom = 0
    do iDtype = 1, atoms%ntype
      do iDeqat = 1, atoms%neq(iDtype)
        iDatom = iDatom + 1
        do idir = 1, 3
          iatom = 0
          do itype = 1, atoms%ntype
            do ieqat = 1, atoms%neq(itype)
              iatom = iatom + 1
              do iG = 1, ngpqdp2km
                temp1 = rho1IRctC(iG, iDatom, idir) * phFac(iG, iatom)
                do oqn_l = 0, atoms%lmax(itype)
                  lm_temp = oqn_l * (oqn_l + 1) + 1
                  do mqn_m = -oqn_l, oqn_l
                    lm = lm_temp + mqn_m
                    temp2 = fpiul(oqn_l) * conjg(ylm(iG, lm))
                    do imesh = 1, atoms%jri(itype)
                      rho1MTctC(imesh, lm, iatom, idir, iDatom) = rho1MTctC(imesh, lm, iatom, idir, iDatom) + temp1 * temp2 &
                        & * sbes(imesh, oqn_l, iG, itype)
                    end do
                  end do
                end do
              end do
            end do
          end do
        end do
      end do
    end do

  end subroutine calcVarCTcorr

  !---------------------------------------------------------------------------------------------------------------------------------
  !> @author
  !> Christian-Roman Gerhorst, Forschungszentrum Jülich: IAS1 / PGI1
  !>
  !> @brief
  !> Calculates the gradient of the pseudo-core Gauss function density in the displaced muffin-tin
  !>
  !> @details
  !> See also 7.26b (dissertation CRG)
  !>
  !> @param[in]   atoms             : Atoms type, see types.f90
  !> @param[in]   lcApsDens         : Lowercase A of Gauss curve pseudo core-density, see 7.26b (dissertation CRG)
  !> @param[in]   ucApsDens         : Uppercase A of Gauss curve pseudo core-density, see 7.26b (dissertation CRG)
  !> @param[out]  rho1MTAlphPSCcorr : Spherical-harmonic expansion coefficients of the gradient of the pseudo core-density ( Gauss
  !>                                  curve ) at the displaced atom.
  !---------------------------------------------------------------------------------------------------------------------------------
  subroutine calcPDCinAlph( atoms, lcApsDens, ucApsDens, rho1MTAlphPSCcorr )

    use m_types, only : t_atoms

    implicit none

    ! Type parameter
    type(t_atoms),              intent(in)  :: atoms

    ! Array parameter
    real,                       intent(in)  :: lcApsDens(:) ! lowercase A of pseudoDensity
    real,                       intent(in)  :: ucApsDens(:) ! uppercase A of pseudoDensity
    complex,                    intent(out) :: rho1MTAlphPSCcorr(:, : , :, :)

    ! Scalar variables
    integer                                 :: itype
    integer                                 :: imesh
    integer                                 :: lm
    integer                                 :: iatom
    integer                                 :: ieqat
    integer                                 :: idir
    integer                                 :: ikpt

    ! Array variables
    real,          allocatable              :: prFac(:)
    real,          allocatable              :: beforeSum(:, :)

    allocate( prFac(atoms%ntype), beforeSum(atoms%jmtd, atoms%ntype) )
    prFac(:) = 0.
    beforeSum(:, :) = 0.

    rho1MTAlphPSCcorr(:, :, :, :) = cmplx(0., 0.)

    ! todo is this really a itype story or iatom for the lc and uc Aps
    prFac(:) = 0.
    do itype = 1, atoms%ntype
      prFac(itype) = 2. * lcApsDens(itype) * ucApsDens(itype)
      beforeSum(:, :) = 0.
      do imesh = 1, atoms%jri(itype)
        beforeSum(imesh, itype) = prFac(itype) * atoms%rmsh(imesh, itype) * exp( -lcApsDens(itype) * atoms%rmsh(imesh, itype)**2 )
      end do
    end do

    ! todo with testing!!!! transpose c_im
    do idir = 1, 3
      iatom = 0
      do itype = 1, atoms%ntype
        do ieqat = 1, atoms%neq(itype)
          iatom = iatom + 1
          do lm = 2, 4
            do imesh = 1, atoms%jri(itype)
              rho1MTAlphPSCcorr(imesh, lm, iatom, idir) = rho1MTAlphPSCcorr(imesh, lm, iatom, idir) + beforeSum(imesh, itype) &
                                                                                                              & * c_im(idir, lm - 1)
            end do
          end do
        end do
      end do
    end do

  end subroutine calcPDCinAlph

end module m_jpDens1stVar
