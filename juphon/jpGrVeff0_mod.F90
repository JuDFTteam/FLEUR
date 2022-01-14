!----------------------------------------------------------------------------------------------------------------------------------------
! Forschungszentrum Jülich, juPhon Plugin for the FLEUR program
!----------------------------------------------------------------------------------------------------------------------------------------
!
! MODULE: Gradient of unperturbed effective potential
!
!> @author
!> Christian-Roman Gerhorst
!>
!> @note
!> Many routines are very similiar to the routines of FLEUR.
!>
!> @brief
!> This module contains routines related to the calculation of the gradient of the unperturbed effective potential.
!>
!> @todo
!> Complete documentation
!> Put in dimensions of array?
!>
!> @note
!> Additional information and formulas pointing out the routines of this module can be found within this
!> <a href='jpGrVeff0.pdf'>document</a>.
!----------------------------------------------------------------------------------------------------------------------------------------
module m_jpGrVeff0

    USE m_constants
  implicit none

  contains

  !--------------------------------------------------------------------------------------------------------------------------------------
  !> @author
  !> Christian-Roman Gerhorst, Forschungszentrum Jülich: IAS1 / PGI1
  !>
  !> @brief
  !> Main routine to calculate the gradient of the unperturbed effective potential.
  !>
  !> @details
  !> The first variation of the Coulomb potential is calculated using the Weinert method without using symmetry by not exploiting
  !> stars and lattice harmonics according to M.Weinert, J.Math.Phys. 22(11) (1981) p.2434 and the PhD thesis of Aaron Klüppelberg.
  !> For the xc-potential, the functional derivative of the kernel with respect to the density is connected with the gradient of the
  !> density. In the interstitial region, this is done by a convolution, whereas in the muffin-tin a projection onto a spherical
  !> harmonic is evaluted by using Gauss quadrature.
  !>
  !> @note
  !> At the moment, only the x-alpha kernel functional derivative is implemented, but its calculation is in seperate routines, that
  !> can be replaced by the libxc library later.
  !>
  !> @param[in]  atoms         : Atoms type, see types.f90
  !> @param[in]  cell          : Unit cell type, see types.f90
  !> @param[in]  dimensions    : Dimension type, see types.f90
  !> @param[in]  stars         : Stars type, see types.f90
  !> @param[in]  ngdp          : Number of G-vectors for potentials and densities
  !> @param[in]  harSw         : Switch to enable Hartree contribution
  !> @param[in]  extSw         : Switch to enable external contribution
  !> @param[in]  xcSw          : Switch to enable xc contribution
  !> @param[in]  testGoldstein : Switch needed by a test using the Goldstone condition
  !> @param[in]  grRhoTermSw   : Switch to activate volume term contribution of muffin-tin boundary problem
  !> @param[in]  gdp           : G-vectors of potentials and densities
  !> @param[in]  rho0IRpw      : Plane-wave coefficients of the unperturbed and converged interstitial density parsed from Fleur
  !> @param[in]  rho0MTsh      : Spherical harmonic coefficients of unperturbed and converged muffin-tin density parsed from Fleur
  !> @param[in]  grRho0IR      : Plane-wave coefficients of the gradient of the interstitial unperturbed density
  !> @param[in]  grRho0MT      : Spherical harmonic coefficients of the gradient of the muffin-tin unperturbed density
  !> @param[in]  grVxcIRKern   : Plane-wave interstitial coefficients of the functional derivative of the xc-kernel with respect to
  !>                             the density
  !> @param[in]  gausWts       : Gauss weights for Gauss quadrature created in m_jppotdenshelper::calcmtdvxckern
  !> @param[in]  ylm           : Set of spherical harmonics whose arguments are the unit vectors of the Gauss mesh points up to lmax
  !>                             created in m_jppotdenshelper::calcmtdvxckern
  !> @param[in]  dKernMTGpts   : Spherical harmonic muffin-tin coefficients of the functional derivative of the xc-kernel with
  !>                             respect to the density created in m_jppotdenshelper::calcmtdvxckern
  !> @param[out] grVeff0IR     : Plane-wave insterstitial coefficients of the gradient of the unperturbed effective potential
  !> @param[out] grVeff0MT     : Spherical harmonic muffin-tin coefficients of the gradient of the unperturbed effective potential
  !--------------------------------------------------------------------------------------------------------------------------------------
  subroutine GenGrVeff0( atoms, cell, stars, ngdp, harSw, extSw, xcSw, gdp, rho0IRpw, rho0MTsh, grRho0IR, grRho0MT, &
      & gausWts, ylm, dKernMTGPts, grVxcIRKern, testGoldstein, grRhoTermSw, grVeff0IR, grVeff0MT )

    USE m_constants
    use m_types, only : t_atoms, t_cell, t_stars

    implicit none

    ! Type parameters
    type(t_atoms),                  intent(in)  :: atoms
    type(t_cell),                   intent(in)  :: cell
    type(t_stars),                  intent(in)  :: stars

    ! Scalar parameter
    integer,                        intent(in)  :: ngdp
    logical,                        intent(in)  :: harSw
    logical,                        intent(in)  :: extSw
    logical,                        intent(in)  :: xcSw
    logical,                        intent(in)  :: testGoldstein
    logical,                        intent(in)  :: grRhoTermSw

    ! Array parameters
    integer,                        intent(in)  :: gdp(:, :)
    complex,                        intent(in)  :: rho0IRpw(:, :)
    complex,                        intent(in)  :: rho0MTsh(:, :, :, :)
    complex,                        intent(in)  :: grRho0IR(:, :)
    complex,                        intent(in)  :: grRho0MT(:, :, :, :)
    complex,                        intent(in)  :: grVxcIRKern(:)
    real,                           intent(in)  :: gausWts(:)
    complex,                        intent(in)  :: ylm(:, :)
    real,                           intent(in)  :: dKernMTGPts(:, :, :)
    complex,           allocatable, intent(out) :: grVeff0IR(:, :)
    complex,           allocatable, intent(out) :: grVeff0MT(:, :, :, :)

    ! Local scalar variables
    real                                        :: normGext
    integer                                     :: iG
    integer                                     :: idir
    integer                                     :: iatom
    integer                                     :: itype
    integer                                     :: ieqat
    integer                                     :: lm
    integer                                     :: oqn_l
    integer                                     :: mqn_m
    integer                                     :: nfftx
    integer                                     :: nffty
    integer                                     :: nfftz
    integer                                     :: nfftxy
    integer                                     :: GxFFT
    integer                                     :: GyFFT
    integer                                     :: GzFFT
    integer                                     :: imesh
    logical                                     :: didvext

    ! Local array variables
    complex,           allocatable              :: qlmGrVc0(:, :, :)
    complex,           allocatable              :: qlmGrVh0Vol(:, :, :)
    complex,           allocatable              :: psqGrVc0(:, :)
    complex,           allocatable              :: grVxc0IR(:, :)
    complex,           allocatable              :: grVxc0MT(:, :, :, :)
    integer,           allocatable              :: pdG2FouM(:)
    complex,           allocatable              :: qlmGrVh0Surf(:, :)
    real                                        :: Gext(3)


    if (harSw .or. extSw) then
      allocate( qlmGrVc0(( atoms%lmaxd + 1 )**2, atoms%nat, 3) )
      qlmGrVc0(:, :, :) = cmplx(0., 0.)
    end if ! harSw .or. extSw

    if (harSw) then
     ! Multipole moments to calculate the gradient of unperturbed Hartree potential using formulas 7.58, 7.59, 7.60, 4.17 and 3.30
     ! from PhDthesAK.
     ! NOTE: in Equation 4.17 within PhDthesAK, there are typos and errors, one should consult PhDthesCRG
      allocate( qlmGrVh0Vol(( atoms%lmaxd + 1 )**2, atoms%nat, 3) )
      qlmGrVh0Vol(:, :, :) = cmplx(0., 0.)
      call CalcQlmGrVh0Vol( atoms, cell, ngdp, gdp, grRho0IR, grRho0MT, qlmGrVh0Vol )

      do idir = 1, 3
        iatom = 0
        do itype = 1, atoms%ntype
          do ieqat = 1, atoms%neq(itype)
            iatom = iatom + 1
            do oqn_l = 0, atoms%lmax(itype)
              do mqn_m = -oqn_l, oqn_l
                lm = oqn_l * (oqn_l + 1) + 1 + mqn_m
                qlmGrVc0(lm, iatom, idir) = qlmGrVc0(lm, iatom, idir) + qlmGrVh0Vol(lm, iatom, idir)
              end do ! mqn_m
            end do ! oqn_l
          end do ! ieqat
        end do !itype
      end do ! idir

      ! Additional multipole moment dependent on discontinuity of density at the MT boundary that is shared for gradient of
      ! unperturbed Hartree potential and the linear variation of the Hartree potential using formulas 7.54, 3.30 and 7.56 from
      ! PhDthesAK
      iatom = 0
      do itype = 1, atoms%ntype
        do ieqat = 1, atoms%neq(itype)
          iatom = iatom + 1

          allocate( qlmGrVh0Surf( ( atoms%lmax(itype) + 1)**2, 3 ) )
          qlmGrVh0Surf(:, :) = cmplx(0., 0.)
          call CalcQlmHarSurf( atoms, cell, itype, iatom, ngdp, gdp, rho0IRpw, rho0MTsh, qlmGrVh0Surf )

          do idir = 1, 3
            do oqn_l = 0, atoms%lmax(itype)
              do mqn_m = -oqn_l, oqn_l
                lm = oqn_l * (oqn_l + 1) + 1 + mqn_m
                qlmgrvc0(lm, iatom, idir) = qlmgrvc0(lm, iatom, idir) - qlmGrVh0Surf(lm, idir)
              end do ! mqn_m
            end do ! oqn_l
          end do ! idir
          deallocate(qlmGrVh0Surf)

        end do ! ieqat
      end do ! itype

    end if ! harSw

    ! N_check: qlmgrvc0(lm, iatom, idir)-=\sum{\alpha} qlmGrVh0Surf(lm, idir)[\alpha]
    !
    ! Consistent with Aaron (7.62) and (2.41) in CRG 16.12.2020-19:00

    ! the atom loop is over all atoms alpha which are displaced this has to be considered later when adding to other quantities
    if ( extSw ) then
      do idir = 1, 3
        iatom = 0
        do itype = 1, atoms%ntype
          do ieqat = 1, atoms%neq(itype)
            iatom = iatom + 1
            do lm = 2, 4
              ! This is + instead of - because we have a - from 7.45 which can be drawn until here.
              qlmGrVc0(lm, iatom, idir) = qlmGrVc0(lm, iatom, idir) + atoms%zatom(itype) * 3 / 4 / pi_const * c_im(idir, lm - 1)
            end do ! lm
          end do ! ieqat
        end do ! itype
      end do ! idir
    end if ! extSw

    ! N_check: qlmGrVc0(1m, iatom, idir)+=\frac{3Z_{\alpha}}{4\pi}c_im(idir,m)
    !
    ! Consistent with Aaron (7.38) and (2.21) in CRG 16.12.2020-19:00 [- omitted; also in (2.27a)]

    ! Here, first the plane-wave contribution to the pseudocharge of the gradient of the unperturbed Hartree potential is
    ! evaluated. Then the pseudo charges are calculated for the Hartree case, as well as for the gradient of the unperturbed
    ! external potential. Finally, the pseudocharges (in the case of the external potential already summed over all atoms)
    ! are used to get the interstitial contribution of the potentials discussed here. According to 7.45b the external potential part
    ! needs a minus.

    ! Calculate the pseudo-charge for the Coulomb potential in general or the Hartree potential or the external potential
    if (harSw .or. extSw) then

      allocate(psqGrVc0(ngdp, 3))
      psqGrVc0(:, :) = cmplx(0., 0.)

      call psqpwVeclp1(atoms, cell, ngdp, grRho0IR, harSw, extSw, gdp, qlmGrVc0, psqGrVc0)

      deallocate(qlmGrVc0)

    end if ! harSw .or. extSw

    allocate(grVeff0IR(ngdp, 3))
    grVeff0IR(:, :) = cmplx(0., 0.)
    if ( harSW .or. extSw ) then
      do idir = 1, 3
        do iG = 1, ngdp
          Gext = matmul(cell%bmat, gdp(:, iG))
          normGext = norm2(Gext)
          if (normGext == 0) then
            cycle
          end if
          grVeff0IR(iG, idir) = fpi_const * psqGrVc0(iG, idir) / normGext**2
        end do ! iG
      end do ! idir
      deallocate(psqGrVc0)
    end if ! harSW .or. extSw

    allocate(grVeff0MT(atoms%jmtd, (atoms%lmaxd + 1)**2, 3, atoms%nat)) ! todo change order of direction and nat
    grVeff0MT(:, :, :, :) = cmplx(0., 0.)
    if ( harSW .or. extSw ) then

      iatom = 0
      do itype = 1, atoms%ntype
        do ieqat = 1, atoms%neq(itype)
          iatom = iatom + 1
          ! Note: We only have to subtract the external volume integral contribution in the displaced atom. However, the gradient of
          !       the external potential or the effective potential is only used in the muffin-tin that is displaced, except for the
          !       the case when it is plotted, but this does not justify another atom dimension. For plotting we have to create a
          !       new array in the plotting routine where grRhoTermSw is enabled in the displaced atom and disabled anywhere else.
          !
          !       Tests have shown that for kind=8 accuracy, it is not important whether one calculates the external and the
          !       Hartree potential seperately or not. By default and to have consistency with Fleur the interstitial boundary value
          !       is the Coulomb potential, i.e. the Hartree and the external potential.
          call vmtsCoul(atoms, cell, ngdp, iatom, itype, harSw, extSW, grVeff0IR, grRho0MT, gdp, &
            & grVeff0MT(:, :, :, iatom), testGoldstein, grRhoTermSw )
        end do ! ieqat
      end do ! itype

    end if ! harSW .or. extSw

    if ( xcSw ) then
      ! Prepare mapping array between a G-vector and its respective point on the FFT mesh, according to how the FFT routine requires it.
      allocate(pdG2FouM(ngdp))
      pdG2FouM(:) = 0
      nfftx = 3 * stars%mx1
      nffty = 3 * stars%mx2
      nfftz = 3 * stars%mx3
      nfftxy= 9 * stars%mx1 * stars%mx2

      do iG = 1, ngdp
        GxFFT = gdp(1, iG)
        GyFFT = gdp(2, iG)
        GzFFT = gdp(3, iG)
        if (GxFFT < 0) GxFFT = GxFFT + nfftx
        if (GyFFT < 0) GyFFT = GyFFT + nffty
        if (GzFFT < 0) GzFFT = GzFFT + nfftz
        pdG2FouM(iG) = GxFFT + GyFFT * nfftx + GzFFT * nfftxy
      end do
      allocate(grVxc0IR(ngdp, 3))
      grVxc0IR(:, :) = cmplx(0., 0.)
      do idir = 1, 3
        call convolGrRhoKern(stars, ngdp, ngdp, grRho0IR(:, idir), grVxcIRKern, pdG2FouM, pdG2FouM, grVxc0IR(:, idir), 1)
      end do
      deallocate(pdG2FouM)

      ! Add x-alpha xc-potential contribution for IR
      do idir = 1, 3
        do iG = 1, ngdp
          grVeff0IR(iG, idir) = grVeff0IR(iG, idir) + grVxc0IR(iG, idir)
        end do ! iG
      end do ! idir
      deallocate(grVxc0IR)

      ! During the normal Sternheimer iterations the density variation is passed to the routines without the gradient of the density
      ! coming from the basis set variation. Therefore, we do not calculate terms accounting for them an switch off terms canceling
      ! anyway.
      if ( grRhoTermSw ) then
      ! Add x-alpha xc-potential contribution for MT
        allocate( grVxc0MT(atoms%jmtd, ( atoms%lmaxd + 1)**2, atoms%nat, 3) )
        grVxc0MT(:, :, :, :) = cmplx(0., 0.)

        ! Note: we leave the lower block of dead code, because in near future, there will be an optimization, after which it will
        ! become required to divide out the r^2 again.
        !do idir = 1, 3
        !  iatom = 0
        !  do itype = 1, atoms%ntype
        !    do ieqat = 1, atoms%neq(itype)
        !      iatom = iatom + 1
        !      do oqn_l = 0, atoms%lmax(itype) + 1
        !        do mqn_m = -oqn_l, oqn_l
        !          lm = oqn_l * (oqn_l + 1) + 1 + mqn_m
        !          do imesh = 1, atoms%jri(itype)
        !             grRho0MT(imesh, lm, iatom, idir) = grRho0MT(imesh, lm, iatom, idir) / atoms%rmsh(imesh, itype)                    &
        !               & / atoms%rmsh(imesh, itype)
        !          end do
        !        end do
        !      end do
        !    end do
        !  end do
        !end do

        call convolMTgrVeff0dKern(atoms, grRho0MT, dKernMTGPts, gausWts, ylm, grVxc0MT)

        iatom = 0
        do itype = 1, atoms%ntype
          do ieqat = 1, atoms%neq(itype)
            iatom = iatom + 1
            do idir = 1, 3
            ! This is not really valid with the l + 1 and l+ 2 because the grVxc0MT is not calculated correctly with the gauß mesh
              do lm = 1, (atoms%lmax(itype) + 1)**2
                do imesh = 1, atoms%jri(itype)
                  grVeff0MT(imesh, lm, idir, iatom) = grVeff0MT(imesh, lm, idir, iatom) + grVxc0MT(imesh, lm, iatom, idir)
                end do
              end do
            end do
          end do
        end do
      end if
    end if

  end subroutine genGrVeff0

  !--------------------------------------------------------------------------------------------------------------------------------------
  !> @author
  !> Christian-Roman Gerhorst, Forschungszentrum Jülich: IAS1 / PGI1
  !>
  !> @brief
  !> Calculates the multipole moments of the gradient of the Hartree potential incorporating a gradient of the density as real
  !> charge.
  !>
  !> @param[in]  atoms    : Atoms type, see types.f90
  !> @param[in]  cell     : Unit cell type, see types.f90
  !> @param[in]  ngdp     : Number of G-vectors for potentials and densities
  !> @param[in]  gdp      : G-vectors of potentials and densities
  !> @param[in]  grRho0IR : Plane-wave coefficients of the gradient of the interstitial unperturbed density
  !> @param[in]  grRho0MT : Spherical harmonic coefficients of the gradient of the muffin-tin unperturbed density
  !> @param[out] qlmGrVc0 : Multipole coefficients of the gradient of the Hartree potential incorporating a volume integral over the
  !>                        gradient of the unperturbed density
  !--------------------------------------------------------------------------------------------------------------------------------------
  subroutine CalcQlmGrVh0Vol( atoms, cell, ngdp, gdp, grRho0IR, grRho0MT, qlmGrVc0 )

    use m_types, only : t_atoms, t_cell
    use m_intgr, only : Intgr3!LinIntp ! TODO: Is this ok?
    use m_sphbes, only : Sphbes

    implicit none

    ! Type Parameters
    type(t_atoms),               intent(in)  :: atoms
    type(t_cell),                intent(in)  :: cell

    ! Scalar Parameters
    integer,                     intent(in)  :: ngdp

    ! Array Parameters
    integer,                     intent(in)  :: gdp(:, :)
    complex,                     intent(in)  :: grRho0IR(:, :)
    complex,                     intent(in)  :: grRho0MT(:, :, :, :)
    complex,                     intent(out) :: qlmGrVc0(:, :, :)

    ! Local Variables:
    !
    ! iatom       : runs over all atoms
    ! lm          : encodes oqn_l and mqn_m
    ! itype       : runs over all atom types
    ! idir      : runs over 3 directions the atom can be displaced to
    ! tempGaunt1  : auxillary variable to store a Gaunt coefficient
    ! tempGaunt2  : auxillary variable to store a Gaunt coefficient
    ! imesh       : runs over mesh points of current grid
    ! rl2         : stores R^(l + 2)
    ! fint        : stores integral
    ! ll1         : auxillary variable to calculate lm
    ! mqn_m       : magnetic quantum number m
    ! sk3r        : stores argument of spherical Bessel function
    ! mqn_mpp     : magnetic quantum number m", also used for indexing 3 directions the atom can be displaced to
    ! oqn_l       : orbital quantum number l
    ! iGvec       : indexes current G-vector
    ! gradrho0PWi : stores the the second part of equation 7.58 using equation 7.60 from PhD thesis Aaron Klüppelberg
    ! pylm        : contains 4π i i^l G / |G| exp(i G τ)  Y*_lm(G / |G|)
    ! sbes        : stores the Bessel function
    ! cil         : stores everything except for pylm of the result of this routine
    ! f           : stores the integrand of the first term in (7.58, PhD thesis Aaron Klüppelberg)

    ! Local Scalar Variables
    integer                                  :: iatom
    integer                                  :: lm
    integer                                  :: itype
    integer                                  :: idir
    integer                                  :: imesh
    real                                     :: rl2
    real                                     :: intgrResR
    real                                     :: intgrResI
    real                                     :: ll1
    integer                                  :: mqn_m
    integer                                  :: ieqat
    real                                     :: sk3r
    integer                                  :: oqn_l
    integer                                  :: iG
    real                                     :: normGext
    complex                                  :: cil
#ifdef DEBUG_MODE
!todo review these variables
    real                                     :: analyticalInt
    logical                                  :: expensiveDebug
#endif DEBUG_MODE

    ! Local Array Variables
    complex,        allocatable              :: grRhoIRlm(:, :, :)
    complex,        allocatable              :: pylm(:, :)
    real,           allocatable              :: sbes(:)
    real,           allocatable              :: intgrR(:)
    real,           allocatable              :: intgrI(:)
    real                                     :: Gext(3)

#ifdef DEBUG_MODE
!todo review these variables
    complex,        allocatable              :: gradrho0PWiTest(:, :, :)
    real,           allocatable              :: integrandTest(:)
    real,           allocatable              :: sbesIntegr(:, :, :)
    complex,        allocatable              :: pylmOld(:, :, :)
#endif

    ! Initialize dummy assumed shape array
    allocate( intgrR(atoms%jmtd), intgrI(atoms%jmtd) )
    intgrR(:) = 0.
    intgrI(:) = 0.
    qlmGrVc0(:, :, :) = cmplx(0., 0.)

    ! Calculate the first term in 7.58 using 7.59 (in 7.59 there are mistakes, please refer to dissCRG).
    do idir = 1, 3
      iatom = 0
      do itype = 1, atoms%ntype
        do ieqat = 1, atoms%neq(itype)
          iatom = iatom + 1
          do oqn_l = 0, atoms%lmax(itype)
            do mqn_m = -oqn_l, oqn_l
              lm = oqn_l * ( oqn_l + 1 ) + 1 + mqn_m
              intgrR = 0
              intgrI = 0
              do imesh = 1, atoms%jri(itype)
                intgrR(imesh) = atoms%rmsh(imesh, itype)**(oqn_l + 2) * real( grRho0MT(imesh, lm, iatom, idir) )
                intgrI(imesh) = atoms%rmsh(imesh, itype)**(oqn_l + 2) * aimag( grRho0MT(imesh, lm, iatom, idir) )
              end do ! imesh
              call Intgr3( intgrR(1:atoms%jri(itype)), atoms%rmsh(1:atoms%jri(itype), itype), atoms%dx(itype), atoms%jri(itype), intgrResR )  ! TODO: Is this ok?
              call Intgr3( intgrI(1:atoms%jri(itype)), atoms%rmsh(1:atoms%jri(itype), itype), atoms%dx(itype), atoms%jri(itype), intgrResI )  ! TODO: Is this ok?
              qlmGrVc0(lm, iatom, idir) = cmplx( intgrResR, intgrResI )
            end do ! mqn_m
          end do ! oqn_l
        end do ! ieqat
      end do ! itype
    end do ! idir

    ! N_check: qlmGrVc0(lm,iatom,idir)=\int_{0}^{R_{MT^{\alpha}}}dr r^{l+2} grRho0MT(r, lm, iatom, idir)
    ! Consistent with Aaron (7.58) part 1 and (2.35) in CRG 16.12.2020-19:00

    ! This block calculates the second part of equation 7.58, where the integrand is rewritten using using equation 7.60 in PhDthesAK
    ! and then can be evaluated in similiar manner as in 7.57 PhDthesAK analytically.
    ! 4 π i^l i \sum_G  G / |G| ρ^(0)(G) Y*_lm(G / |G|) exp(i G τ_β) R_β^(l + 2) j_(l + 1)(|G|R_β) .

    allocate( grRhoIRlm(( atoms%lmaxd + 1 )**2, atoms%nat, 3) )
    allocate( pylm(( atoms%lmaxd + 1 )**2, atoms%nat ) )
    allocate( sbes(0:atoms%lmaxd + 1) )
    grRhoIRlm(:, :, :) = cmplx(0., 0.)
    pylm(:, :) = cmplx(0., 0.)
    sbes(:) = 0.

    do idir = 1, 3
      do iG = 1, ngdp
        !call Phasy1nSymVeclp1( atoms, cell, gdp(1:3, iG), pylm )
        ! calculates 4 π i^l exp(i Gvec * taual) Y*_lm(Gvec / |Gvec|)
        pylm(:, :) = cmplx(0., 0.)
        call Phasy1nSym( atoms, cell, gdp(1:3, iG), [0., 0., 0.], pylm)
        Gext(1:3) = matmul( cell%bmat(1:3, 1:3), gdp(1:3, iG) )
        normGext = norm2( Gext(1:3) )
        if (normGext < 1e-12) cycle
        iatom = 0
        do itype = 1, atoms%ntype
          do ieqat = 1, atoms%neq(itype)
            iatom = iatom + 1
            sk3r = normGext * atoms%rmt(itype)
            sbes(:) = 0.
            ! calculates spherical bessel function with argument sk3r
            call Sphbes( atoms%lmax(itype) + 1, sk3r, sbes )
            rl2 = atoms%rmt(itype)**2
            do oqn_l = 0, atoms%lmax(itype)
              cil = sbes(oqn_l + 1) * rl2 * grRho0IR(iG, idir) / normGext
              ll1 = oqn_l * (oqn_l + 1) + 1
              do mqn_m = -oqn_l, oqn_l
                lm = ll1 + mqn_m
                grRhoIRlm(lm, iatom, idir) = grRhoIRlm(lm, iatom, idir) + cil * pylm(lm, iatom)
              end do ! m
              rl2 = rl2 * atoms%rmt(itype)
            end do  ! l
          end do !ieqat
        end do ! itype
      end do ! iG
    end do ! idir

    ! N_check: grRhoIRlm(lm, iatom, idir)=\sum_{\bm{G}\neq\bm{0}} \frac{j_{l+1}(GR_{\alpha})R_{\alpha}^{l+2}}{G}grRho0IR (G, idir)*
    !                                     4\pi i^{l} e^{i2\pi \bm{G}_{int}\cdot\bm{tau}_{\alpha,int}}Y_{lm}^{*}(\har\bm{G})
    ! Consistent with Aaron (7.58) part 2 and (2.35) in CRG 16.12.2020-19:00

    ! Unit the MT and the IR contribution of 7.58 in PhDthesAK
    do idir = 1, 3
      iatom = 0
      do itype = 1, atoms%ntype
        do ieqat = 1, atoms%neq(itype)
          iatom = iatom + 1
          do oqn_l = 0, atoms%lmax(itype)
            do mqn_m = -oqn_l, oqn_l
              lm = oqn_l * (oqn_l + 1) + 1 + mqn_m
              qlmGrVc0(lm, iatom, idir) = qlmGrVc0(lm, iatom, idir) - grRhoIRlm(lm, iatom, idir)
            end do ! mqn_m
          end do ! oqn_l
        end do ! ieqat
      end do ! itype
    end do ! idir

    ! N_check:
    ! Consistent with (2.35) in CRG 16.12.2020-19:00

    if (.false.) then
      ! Debug output for maintenance
      do idir = 1, 3
        do oqn_l = 0, atoms%lmax(1)
          do mqn_m = -oqn_l, oqn_l
            lm = oqn_l * ( oqn_l + 1 ) + 1 + mqn_m
            write(2022, '(2(i8),2(f15.8))') lm, idir, qlmGrVc0(lm, 1, idir)
          end do ! mqn_m
        end do ! oqn_l
      end do ! idir
    end if

!#ifdef DEBUG_MODE
!
!  ! This tests whether the numerical and the analytical integral of r^(l+2) j_l(GR) is the same
!  !TODO not sure if this is working anymore, is this even important to test anymore, will not be formatted now?
!  !TODO format this and review
!  expensiveDebug = .false.
!  if (expensiveDebug) then
!    allocate( gradRho0PWiTest((atoms%lmaxd + 2)**2, atoms%nat, 3) )
!    allocate( integrandTest(atoms%jmtd) )
!    allocate( sbesIntegr(0:atoms%lmaxd + 1, ngdp, atoms%ntype) )
!    gradRho0PWiTest = cmplx(0., 0.)
!    integrandTest = 0.
!    sbesIntegr(:, :, :) = 0.
!
!    write (*, '(a)') 'l, numerical integral of spherical lth Bessel function, analytical version and differences'
!    do itype = 1, atoms%ntype
!    write (*, '(a, i1)') 'atom type ', itype
!      do iG = 1, ngdp
!        Gext = matmul(cell%bmat, gdp(:, iG))
!        normGext = norm2(Gext)
!        if (normGext == 0) cycle
!        do oqn_l = 0, atoms%lmax(itype) + 1
!          integrandTest = 0
!          do imesh = 1, atoms%jri(itype)
!            sk3r = normGext * atoms%rmsh(imesh, itype)
!            call sphbes(atoms%lmax(itype) + 2, sk3r, sbes)
!            integrandTest(imesh) = sbes(oqn_l) * atoms%rmsh(imesh, itype)**(oqn_l + 2)
!          end do ! imesh
!          !call intgr3(integrandTest, atoms%rmsh(:, itype), atoms%dx(itype), atoms%jri(itype), intgrResR)
!          call intgr3LinIntp(integrandTest, atoms%rmsh(:, itype), atoms%dx(itype), atoms%jri(itype), intgrResR, 1)
!          sbesIntegr(oqn_l, iG, itype) = intgrResR
!          analyticalInt = sbes(oqn_l + 1) * atoms%rmt(itype)**(oqn_l + 3) / atoms%rmt(itype)**(1) / normGext
!          if (iG == 1) then
!          if (abs(analyticalInt - intgrResR) > 1e-9) then
!            write (*, '(i2,2x,3(es15.8,2x))') oqn_l, intgrResR, analyticalInt, abs(analyticalInt - intgrResR)
!          end if
!        end if
!        end do
!      end do
!    end do
!
!
!    !todo optimize loop structure
!    allocate(pylmOld((atoms%lmaxd + 1)**2, atoms%nat, 3))
!    do idir = 1, 3
!      do iG = 1, ngdp
!        iatom = 0
!        call Phasy1nSymUVeclp1(atoms, cell, gdp(:, iG), pylmOld) ! todo is this method really necessary and not yet there?, review it! Can be replaced by normal phasy routine, by dividing through norm of G-vector
!        do itype = 1, atoms%ntype
!          do ieqat = 1, atoms%neq(itype)
!            iatom = iatom + 1
!            do oqn_l = 0, atoms%lmax(itype) + 1
!              cil = sbesIntegr(oqn_l, iG, itype) * rho0IRpw(iG)
!              ll1 = oqn_l * (oqn_l + 1) + 1
!              do mqn_m = -oqn_l, oqn_l
!                lm = ll1 + mqn_m
!                gradRho0PWiTest(lm, iatom, idir) = gradRho0PWiTest(lm, iatom, idir) + cil * pylmOld(lm, iatom, idir)
!              end do
!            end do
!          end do
!        end do
!      end do
!    end do
!
!    deallocate(pylmOld)
!
!    write (*, '(a)') 'Analytical version and numerical end result and its difference is displayed if difference not numerically zero'
!    do idir = 1, 3
!      iatom = 0
!      do itype = 1, atoms%ntype
!        do ieqat = 1, atoms%neq(itype)
!          iatom = iatom + 1
!          do oqn_l = 0, atoms%lmax(itype) + 1
!            do mqn_m = -oqn_l, oqn_l
!              lm = oqn_l * (oqn_l + 1) + 1 + mqn_m
!              if ( abs(real(grRhoIRlm(lm, iatom, idir) - gradRho0PWiTest(lm, iatom, idir)) > 1e-9) .or. &
!                &abs( aimag(grRhoIRlm(lm, iatom, idir) - gradRho0PWiTest(lm, iatom, idir)) > 1e-9) ) then
!              ! write (*, *) 'idir= ', idir, 'iatom= ', iatom, 'oqn_l= ', oqn_l, 'mqn_m= ', mqn_m
!                write (*, *) idir, oqn_l, mqn_m
!                write (*, '(3(es15.8,2x))') real(grRhoIRlm(lm, iatom, idir)),  real(gradRho0PWiTest(lm, iatom, idir)), &
!                  &abs(real(grRhoIRlm(lm, iatom, idir) - gradRho0PWiTest(lm, iatom, idir)))
!                write (*, '(3(es15.8,2x))') aimag(grRhoIRlm(lm, iatom, idir)),  aimag(gradRho0PWiTest(lm, iatom, idir)), &
!                  &aimag(grRhoIRlm(lm, iatom, idir) - gradRho0PWiTest(lm, iatom, idir))
!                write (*, *)
!              end if
!            end do
!          end do
!        end do
!      end do
!    end do
!  end if
!
!
!  do idir = 1, 1
!    iatom = 0
!    do itype = 1, atoms%ntype
!      do ieqat = 1, atoms%neq(itype)
!        iatom = iatom + 1
!        do oqn_l = 0, atoms%lmax(itype)
!          do mqn_m = -oqn_l, oqn_l
!            lm = oqn_l * (oqn_l + 1) + 1 + mqn_m
!            if ( ((abs(grRhoIRlm(lm, iatom, idir)) < 1e-9) .and. (abs(qlmGrVc0(lm, 1, idir)) > 1e-9))&
!              .or. ((abs(grRhoIRlm(lm, iatom, idir)) > 1e-9) .and. (abs(qlmGrVc0(lm, 1, idir)) < 1e-9))) then
!              ! lm channels where there is a MT part but no PW part
!              write (*, '(a)') 'No muffin-tin and plane-wave part for lm channels:'
!              write(*, '(i2, 2x, i3, 2x, 2(2(es12.5),2x))') oqn_l, mqn_m, grRhoIRlm(lm, iatom, idir), qlmGrVc0(lm, iatom, idir)
!            end if
!          end do
!        end do
!      end do
!    end do
!  end do
!
!#endif

  end subroutine CalcQlmGrVh0Vol

   !--------------------------------------------------------------------------------------------------------------------------------------
   !> @author
   !> Christian-Roman Gerhorst, Forschungszentrum Jülich: IAS1 / PGI1
   !>
   !> @brief
   !> Calculates the pseudocharges of the gradient of the Hartree, external or Coulomb unperturbed potential.
   !>
   !> @details
   !> This routine generates the fourier coefficients of the pseudo charge density. It
   !> is very similar to psqpw.F90 but lacks any consideration of symmetry. It is based
   !> on equation 28, 30 in M.Weinert J.Math.Phys. 22(11) (1981) p.2434, where a factor
   !> 1 / R^l has been forgotten which can be seen starting from the unsimplified equation
   !> above equation 28 in the Weinert paper. As this routine has been modified to deliver
   !> the multipole moments to calculate interstitial contributions of potentials for the
   !> determination of the first variation of the Coulomb potential due to a phonon
   !> pertubation the result psq of the method is 3-dimensional to mirror the 3 directions
   !> the atoms can be displaced due to a phonon mode and the phonon q-vector is considered.
   !> The respective formula can be found in equation 7.63 in the PhD thesis of Aaron Klüppelberg,
   !> where equations 28 and 30 have been united.
   !>
   !> @param[in]  atoms      : Atoms type, see types.f90
   !> @param[in]  cell       : Unit cell type, see types.f90
   !> @param[in]  dimensions : Dimension type, see types.f90
   !> @param[in]  ngdp       : Number of G-vectors for potentials and densities
   !> @param[in]  harSw      : Switch to enable Hartree contribution
   !> @param[in]  extSw      : Switch to enable external contribution
   !> @param[in]  grRho0IR   : Plane-wave coefficients of the gradient of the interstitial unperturbed density
   !> @param[in]  gdp        : G-vectors of potentials and densities
   !> @param[in]  qlmGrVc0   : Multipole coefficients of the gradient of the Hartree, external or Coulomb unperturbed potential
   !> @param[out] psqGrVc0   : Plane-wave coefficients of the pseudocharge to construct the gradient of the Hartree, external or
   !>                          Coulomb potential depending on what multipole moments are passed in.
   !--------------------------------------------------------------------------------------------------------------------------------------
   subroutine PsqpwVeclp1(atoms, cell, ngdp, grRho0IR, harSw, extSw, gdp, qlmGrVc0, psqGrVc0)

!#include "recycledRoutines/cpp_double.h"

     use m_types, only : t_atoms, t_cell
     use m_sphbes

     implicit none

     ! Variables:
     !
     ! atoms  : atoms type defined in m_types
     ! cell   : unit cell type defined in m_types
     ! dimensions : dimension type defined in m_types
     ! ngdp   : number of G-vectors used for the potential and the density
     ! grRho0IR  : planewave density in the interstitial region
     ! gdp    : array of G-vectors used for the potential and density
     ! qlmGrVc0    : multipole moments to construct the pseudo charge
     ! psq    : resulting 3-D pseudocharge of this routine
     ! itype  : runs over all types
     ! rmtl   : contains R^l for the prefactor
     ! oqn_l  : orbital quantum number l
     ! p      : auxiliary variable to calculate prefactor
     ! nc     : loop variable to calculate prefactor
     ! fpo    : auxiliary variable which contains 1 / Ω
     ! iGvec  : runs over G-vectors which are used for the potential and the density
     ! idir : runs over directions the atoms can be displaced to
     ! iatom  : runs over all atoms
     ! ncvn   : contains n + l for a certain atom type, while n is taken from a table in the Weinert paper
     ! ieqat   : runs over all equal atoms
     ! n1     : auxiliary variable which contains n - l + 1
     ! ll1    : auxiliary variable to calculate lm
     ! mqn_m  : magnetic quantum number
     ! lm     : variable encoding oqn_l and mqn_m
     ! pn     : variable which stores prefactor with double faculties times 1 / R^l
     ! pylm   : contains 4 π i^l exp((Gvec) * taual) Y*_lm(Gvec)
     ! Gext   : G-vector in external coordinates
     ! sbes   : contains spherical bessel function
     ! sl     : auxiliary variable to calculate psq
     ! sm     : auxiliary variable to calculate psq
     ! sa     : auxiliary variable to calculate psq

     !     .. Scalar Arguments ..
     type(t_atoms),                 intent(in)  :: atoms
     type(t_cell),                  intent(in)  :: cell

     integer,                       intent(in)  :: ngdp
     logical,                       intent(in)  :: harSw
     logical,                       intent(in)  :: extSw

     !     .. Array Arguments ..
     complex,                       intent(in)  :: grRho0IR(:, :)
     integer,                       intent(in)  :: gdp(:, :)
     complex,                       intent(in)  :: qlmGrVc0(:, :, :)
     complex,                       intent(out) :: psqGrVc0(:, :)

     !     .. Local Scalars ..
     integer                                    :: itype
     real                                       :: rmtl
     integer                                    :: oqn_l
     real                                       :: p
     integer                                    :: nc
     real                                       :: fpo
     integer                                    :: iG
     integer                                    :: idir
     integer                                    :: iatom
     integer                                    :: ncvn
     integer                                    :: ieqat
     integer                                    :: n1
     integer                                    :: ll1
     integer                                    :: mqn_m
     integer                                    :: lm
     complex                                    :: sl ! sum over l
     complex                                    :: sm ! sum over m
     complex                                    :: sa ! sum over atoms

     !     .. Local Arrays ..
     real,              allocatable             :: pn(:, :)
     complex,           allocatable             :: pylm(:, :)
     real,              allocatable             :: sbes(:)
     real                                       :: Gext(3)

     ! Init assumed shape array dummy
     psqGrVc0(:, :) = cmplx(0., 0.)

     ! This block calculates pn(l,itype) = (2l + 2nc(itype) + 3)!! / ((2l + 1)!! R^l) from equation 28 in the Weinert paper cited above or the
     ! second term of equation 7.63 (e.g.) in the PhD thesis of Aaron Klüppelberg. ncv(itype) is according to Weinert's paper still
     ! constant and defined as n + l and the formula to be calculated is simplified in the following
     ! way: Assume l + n = l', so the factorials cancel each other until l' = l - 1 since 3 is 1 + 2 and gives the next odd number.
     ! Since 2l + 2nc(itype) + 3 is odd we then multiply over all odd numbers from (2 l + 3) until (2ncv(itype) +3) which can be done
     ! by incrementing l until ncv(itype). That's why the loop starts at nc=l.
     ! The value of ncv per type is taken from Table 1 in the Weinert paper and given by the Fleur input .
     allocate( pn(0:atoms%lmaxd, atoms%ntype) )
     pn(:, :) = 0

     do itype = 1, atoms%ntype
        rmtl = 1.0
        do oqn_l = 0, atoms%lmax(itype)
           if (oqn_l < atoms%ncv(itype)) then ! N>=1
              p = 1.
              do nc = oqn_l, atoms%ncv(itype)
                 p = p * (2 * nc + 3)
              end do ! nc
              pn(oqn_l, itype) = p / rmtl
           end if ! oqn_l < atoms%ncv(itype)
           rmtl = rmtl * atoms%rmt(itype)
        end do ! oqn_l
     end do ! itype

     ! N_check: valid.

     ! This block calculates the G /= 0 term in equation 28 of the Weinert paper cited above:
     ! \tilde \rho_s (K) = 4 π / Ω \sum_{lmi} (-i)^l (2l + 2nc(n) + 3)!! / ((2l + 1)!! R^l) j_{l + n + 1}(|G + q|R_i) / (|G + q| R_i)^{n-l+1}
     ! qlmGrVc0 \exp{-i (G + q) τ} Y_{lm}(G + q), which is the second term in equation 7.63. It is similiar to equation 28 in the Weinert
     ! paper cited above.
     ! The G = 0 term is not needed for calculating the interstitial contributions of the potentials and should be zero anyway.
     allocate( pylm((atoms%lmaxd + 1)**2, atoms%nat) )
     allocate( sbes(0:MAXVAL(atoms%ncv) + 1))
     pylm(:, :) = cmplx(0., 0.)
     sbes(:) = 0.

     fpo = 1. / cell%omtil

     do iG = 1, ngdp
       Gext(1:3) = matmul( cell%bmat(1:3, 1:3), real(gdp(1:3, iG)) )
       if ( norm2( Gext(1:3) ) < 1e-12 ) cycle
       pylm(:, :) = cmplx(0., 0.)
       ! Return 4 pi i^l exp(iG tau) Y_lm^{*}(G)
       call Phasy1nSym( atoms, cell, gdp(1:3, iG), [0., 0., 0.], pylm )
       do idir = 1, 3
         ! variable accumulating sum over l resetted
         sa = cmplx(0., 0.)
         iatom = 0
         do itype = 1, atoms%ntype
           ncvn = atoms%ncv(itype)
           ! We need the spherical Bessel functions up to ncvn + 1 + 1 because, ncvn only assumes lmax not lmax + 1
           sbes(:) = 0.
           call sphbes(ncvn + 1, norm2(Gext(1:3)) * atoms%rmt(itype), sbes)
           do ieqat = 1, atoms%neq(itype)
             iatom  = iatom + 1
             ! variable accumulating sum over l resetted
             sl = cmplx(0., 0.)
             do oqn_l = 0, atoms%lmax(itype)
               if ( oqn_l >= ncvn ) cycle ! in earlier versions of FLEUR pn was set zero in this case, so zero was added to sl
               n1 = ncvn - oqn_l + 1
               ll1 = oqn_l * (oqn_l + 1) + 1
               ! variable accumulating sum over m resetted
               sm = cmplx(0., 0.)
               do mqn_m = -oqn_l, oqn_l
                 lm = ll1 + mqn_m
                 sm = sm + qlmGrVc0(lm, iatom, idir) * conjg( pylm(lm, iatom) )
               end do ! mqn_m
               sl = sl + sbes(ncvn + 1) * sm * (pn(oqn_l, itype) / ( ( norm2( Gext(1:3) ) * atoms%rmt(itype) )**n1 ))
             end do ! oqn_l
             sa = sa + sl
           end do ! ieqat
         end do ! itype
         if (harSw) then
           psqGrVc0(iG, idir) = grRho0IR(iG, idir) + fpo * sa
         else if (.not.harSw .and. extSw) then
           psqGrVc0(iG, idir) = fpo * sa
         end if ! special treatment if only external potential required
       end do ! idir
     end do ! iG

     ! N_check: psqGrVc0(iG, idir)=[grRho0IR(iG, idir) +] frac{1}{\Omega}\sum_{\alpha l m} qlmGrVc0(lm, iatom, idir)*
     !                             4\pi (-i)^{l} e^{-i2\pi \bm{G}_{int}\cdot\bm{tau}_{\alpha,int}}Y_{lm}(\hat{\bm{G})*
     !                             \frac{j_{l+N+1}(|GR_{\alpha})}{(GR_{\alpha})^{N+1}}*
     !                             \frac{(2l + 2N + 3)!!}{(2l + 1)!! R_{\alpha}^{l})}
     !
     ! Consistent with Aaron (7.64) and (2.40) in CRG 16.12.2020-19:00

   end subroutine psqpwVeclp1

   !--------------------------------------------------------------------------------------------------------------------------------------
   !> @author
   !> Christian-Roman Gerhorst, Forschungszentrum Jülich: IAS1 / PGI1
   !>
   !> @brief
   !> Solves a boundary problem according to Weinert to gain the muffin-tin contribution of the Hartree, external or
   !> Coulomb potential.
   !>
   !> @param[in]  atoms         : Atoms type, see types.f90
   !> @param[in]  cell          : Unit cell type, see types.f90
   !> @param[in]  ngdp          : Number of G-vectors for potentials and densities
   !> @param[in]  iatom         : Index of atom in whose muffin-tin the calculation is to be performed
   !> @param[in]  itype         : Index of atom type to which atom belongs in whose muffin-tin the calculation is to be performed
   !> @param[in]  harSw         : Switch to enable Hartree contribution
   !> @param[in]  extSw         : Switch to enable external contribution
   !> @param[in]  testGoldstein : Switch needed by a test using the Goldstone condition
   !> @param[in]  grRhoTermSw   : Switch to activate volume term contribution of muffin-tin boundary problem
   !> @param[in]  grVc0IR       : Plane-wave insterstitial coefficients of the gradient of the unperturbed effective potential
   !> @param[in]  grRho0MT      : Spherical harmonic coefficients of the gradient of the muffin-tin unperturbed density
   !> @param[in]  gdp           : G-vectors of potentials and densities
   !> @param[out] grVc0MT       : Spherical harmonic muffin-tin coefficients of the gradient of the unperturbed effective potential
   !--------------------------------------------------------------------------------------------------------------------------------------
  subroutine vmtsCoul( atoms, cell, ngdp, iatom, itype, harSw, extSw, grVc0IR, grRho0MT, gdp, grVc0MT, testGoldstein, grRhoTermSw)

!#include"recycledRoutines/cpp_double.h"
    use m_types, only : t_atoms, t_cell
    use m_intgr, only : intgr2!LinIntp ! TODO: Is this ok?
    use m_sphbes

    implicit none

    ! Type Parameter
    type(t_atoms),              intent(in)  :: atoms
    type(t_cell),               intent(in)  :: cell

    ! Scalar Parameter
    integer,                    intent(in)  :: ngdp
    integer,                    intent(in)  :: iatom
    integer,                    intent(in)  :: itype
    logical,                    intent(in)  :: harSw
    logical,                    intent(in)  :: extSw
    logical,                    intent(in)  :: testGoldstein
    logical,                    intent(in)  :: grRhoTermSw

    ! Array Parameter
    complex,                    intent(in)  :: grVc0IR(:, :)
    complex,                    intent(in)  :: grRho0MT(:,:,:,:)
    integer,                    intent(in)  :: gdp(:, :)
    complex,                    intent(out) :: grVc0MT(:, :, :)

    !----------------------------------------------------------------------------------------------------------------------------------
    ! Local Scalar Variables
    ! oqn_l    : orbital quantum number l
    ! mqn_m    : magnetic quantum number m
    ! l21      : contains 2 l + 1
    ! fpl21    : contains prefactor of volume integral
    ! rmtl     : auxiliary variable to calculate the volume integral
    ! rmt2l    : auxiliary variable to calculate the volume integral
    ! rrlr     : auxiliary variable to calculate the volume integral
    ! ror      : auxiliary variable to calculate the volume integral
    ! lm       : encodes orbital quantum number l and magnetic quantum number m
    ! iGvec    : runs over all G-vectors used here
    ! idir   : runs over all directions the atoms can be displaced to
    ! imesh    : runs over mesh points of local sphere mesh
    !----------------------------------------------------------------------------------------------------------------------------------
    integer                                 :: oqn_l
    integer                                 :: mqn_m
    integer                                 :: l21
    real                                    :: fpl21
    real                                    :: rmtl
    real                                    :: rmt2l
    integer                                 :: lm
    integer                                 :: iG
    integer                                 :: idir
    integer                                 :: imesh
    real                                    :: rrlr
    real                                    :: ror
    real                                    :: prefactor

    !----------------------------------------------------------------------------------------------------------------------------------
    ! Local Array Variables
    ! Gqext    : cartesian representation of G + q
    ! sbf      : contains spherical Bessel functions
    ! pylm     : contains 4 π i^l exp(i (Gvec) * taual) Y*_lm(Gvec )
    ! vtl      : contains the surface integral
    ! rrl      : auxiliary variable to calculate the volume integral
    ! rrl1     : auxiliary variable to calculate the volume integral
    ! f1r      : contains the integrand of the first type of integral occuring for the volume integral part
    ! f2r      : contains the integrand of the second type of integral occuring for the volume integral part
    ! x1r      : contains the result of the the first type of integral occuring for the volume integral part
    ! x2r      : contains the result of the the first type of integral occuring for the volume integral part
    !todo comment imag
    !----------------------------------------------------------------------------------------------------------------------------------
    complex,       allocatable              :: pylm(:, :)
    real,          allocatable              :: sbf(:)
    complex,       allocatable              :: vtl(:, :)
    real,          allocatable              :: rrl(:)
    real,          allocatable              :: rrl1(:)
    complex,       allocatable              :: f1r(:)
    complex,       allocatable              :: f2r(:)
    real,          allocatable              :: f1rReal(:)
    real,          allocatable              :: f2rReal(:)
    real,          allocatable              :: f1rImag(:)
    real,          allocatable              :: f2rImag(:)
    real,          allocatable              :: x1rReal(:)
    real,          allocatable              :: x2rReal(:)
    real,          allocatable              :: x1rImag(:)
    real,          allocatable              :: x2rImag(:)
    real                                    :: Gext(3)

    ! For a given iatom and itype the second term in equation (7.66 / 7.67, PhD thesis Aaron Klüppelberg) is evaluated, beginning from the sum
    allocate( vtl((atoms%lmaxd + 1)**2, 3) )
    allocate( sbf(0:atoms%lmaxd) )
    allocate( pylm(( atoms%lmaxd + 1 )**2, atoms%nat) )
    vtl(:, :) = cmplx(0., 0.)
    sbf(:) = 0.
    pylm(:, :) = cmplx(0., 0.)

    ! Init assumed shape dummy array
    grVc0MT(:, :, :) = cmplx(0., 0.)

    do iG = 1, ngdp
      Gext(1:3) = matmul( cell%bmat(1:3, 1:3), gdp(1:3, iG) )
      if ( norm2( Gext(1:3) ) < 1e-12 ) cycle
      ! Return 4 pi i^l exp(iG tau) Y_lm^{*}(G)
      pylm(:, :) = cmplx(0., 0.)
      call Phasy1nSym( atoms, cell, gdp(1:3, iG), [0., 0., 0.], pylm )
      sbf(:) = 0.
      call Sphbes( atoms%lmax(itype), norm2( Gext(1:3) ) * atoms%rmt(itype), sbf )
      do idir = 1, 3
        do oqn_l = 0, atoms%lmax(itype)
          do mqn_m = -oqn_l, oqn_l
            lm = oqn_l * (oqn_l + 1) + mqn_m + 1
            vtl(lm, idir) = vtl(lm, idir) + grVc0IR(iG, idir) * sbf(oqn_l) * pylm(lm, iatom)
          end do ! mqn_m
        end do ! oqn_l
      end do ! idir
    end do ! iG
    deallocate(pylm)
    deallocate(sbf)

    ! This is not similar to the code beneath because we return.
    if ( .not.grRhoTermSw ) then
      do idir = 1, 3
        do oqn_l = 0, atoms%lmax(itype)
          rmtl = 1. / atoms%rmt(itype)**oqn_l ! 1 / R^l
          do mqn_m = -oqn_l, oqn_l
            lm = oqn_l * ( oqn_l + 1 ) + 1 + mqn_m
            do imesh = 1, atoms%jri(itype)
              ror = atoms%rmsh(imesh, itype)**oqn_l * rmtl ! r^l / R^l
              grVc0MT(imesh, lm, idir) =  ror * vtl(lm, idir)
            end do
          end do
        end do
      end do

      ! N_check: vtl(lm, idir)=\sum_{\bm{G}\new\bm{0}} grVc0IR(G, idir)*
      !                             4\pi i^{l} e^{i2\pi (\bm{G}_{int}+\bm{q}_{int})\cdot\bm{tau}_{\alpha,int}}Y_{lm}^{*}(\hat{\bm{G}})*
      !                             j_{l}(G|R_{\alpha})+(\frac{r_{\alpha}}{\bm{R}_{\alpha}})^{l}
      !
      ! Consistent with Aaron (7.67) line 2 and (2.44) line 2/3 without lm-sum in CRG 16.12.2020-19:00

      ! Maintenance
      if (.false.) then
        do idir = 1, 3
          do oqn_l = 0, atoms%lmax(itype)
            rmtl = 1. / atoms%rmt(itype)**oqn_l ! 1 / R^l
            do mqn_m = -oqn_l, oqn_l
              lm = oqn_l * ( oqn_l + 1 ) + 1 + mqn_m
              do imesh = 1, atoms%jri(itype)
                write(2130, '(3(i8),2(f15.8))') imesh, lm, idir, grVc0MT(imesh, lm, idir)
              end do
            end do
          end do
        end do
      end if
      return
    end if

    ! For a given iatom and itype the first term in equation (7.66 / 7.67, PhD thesis Aaron Klüppelberg) is calculated first and then
    ! the second surface integral contribution inclucing the correct prefactor is added so that the complete Dirichlet boundary value
    ! problem is solved.
    ! Concerning the first term, basically, the integral is splitted into four integrals with respect to the definitions of r_> and r_<.
    ! This gives two types of integrands x1r = r^(l + 2) * ρ(r) dr and x2r = r^(1 - l) * ρ(r) dr. Considering the bounds of the
    ! integrals they can be rearranged to:
    ! 4π / (2 l + 1) * (1 / r^(l + 1) \int_0^r s^(l + 2) ρ ds - r^l / R^(2 * l + 1) \int_0^R ρ s^(l + 2)
    ! + r^l (\int_0^R 1 / s^(l + 1) ρ ds - \int_0^r 1 / s^(l + 1) ρ))

    ! Prefactor is for external contribution
    prefactor = atoms%zatom(itype) * 4 * pi_const / 3
    grVc0MT = cmplx(0., 0.)
    if (harSw) then
      allocate( rrl1(atoms%jmtd), rrl(atoms%jmtd), x1rReal(atoms%jmtd), x2rReal(atoms%jmtd), x1rImag(atoms%jmtd), &
        & x2rImag(atoms%jmtd), f1rReal(atoms%jmtd), f1rImag(atoms%jmtd), f1r(atoms%jmtd), f2rReal(atoms%jmtd), &
        & f2rImag(atoms%jmtd), f2r(atoms%jmtd) )
    end if
    do idir = 1, 3
      ! grRhoTermSw = true implicitely
      if (harSw) then
        rrl1 = 0.
        rrl = 0.
        x1rReal = 0.
        x2rReal = 0.
        x1rImag = 0.
        x2rImag = 0.
        f1rReal = 0.
        f1rImag = 0.
        f1r = cmplx(0., 0.)
        f2rReal = 0.
        f2rImag = 0.
        f2r = cmplx(0., 0.)
        do oqn_l = 0, atoms%lmax(itype)
          l21 = 2 * oqn_l + 1
          fpl21 = fpi_const / l21 ! prefactor 1st term
          rmtl = 1. / atoms%rmt(itype)**oqn_l ! 1 / R^l
          rmt2l = 1. / atoms%rmt(itype)**l21 ! 1 / R^(2 l + 1)
          do imesh = 1, atoms%jri(itype)
            rrl(imesh) = atoms%rmsh(imesh, itype)**oqn_l ! r^l
            rrl1(imesh) = 1. / (rrl(imesh) * atoms%rmsh(imesh, itype)) ! 1 / r^(l + 1)
          end do ! imesh
          do mqn_m = -oqn_l, oqn_l
            lm = oqn_l * (oqn_l + 1) + mqn_m + 1
            do imesh = 1, atoms%jri(itype)
              x1rReal(imesh) = rrl(imesh)  * real( grRho0MT(imesh, lm, iatom, idir) ) * atoms%rmsh(imesh, itype)**2
              x2rReal(imesh) = rrl1(imesh) * real( grRho0MT(imesh, lm, iatom, idir) ) * atoms%rmsh(imesh, itype)**2
              x1rImag(imesh) = rrl(imesh)  * aimag( grRho0MT(imesh, lm, iatom, idir) )* atoms%rmsh(imesh, itype)**2
              x2rImag(imesh) = rrl1(imesh) * aimag( grRho0MT(imesh, lm, iatom, idir) )* atoms%rmsh(imesh, itype)**2
            end do ! mqn_m
            call intgr2(x1rReal, atoms%rmsh(1, itype), atoms%dx(itype), atoms%jri(itype), f1rReal) ! TODO: Is this ok?
            call intgr2(x2rReal, atoms%rmsh(1, itype), atoms%dx(itype), atoms%jri(itype), f2rReal) ! TODO: Is this ok?
            call intgr2(x1rImag, atoms%rmsh(1, itype), atoms%dx(itype), atoms%jri(itype), f1rImag) ! TODO: Is this ok?
            call intgr2(x2rImag, atoms%rmsh(1, itype), atoms%dx(itype), atoms%jri(itype), f2rImag) ! TODO: Is this ok?
            x1rReal(:) = 0.
            x2rReal(:) = 0.
            x1rImag(:) = 0.
            x2rImag(:) = 0.
            f1r(1:atoms%jri(itype)) = cmplx( f1rReal(1:atoms%jri(itype)), f1rImag(1:atoms%jri(itype)) )
            f2r(1:atoms%jri(itype)) = cmplx( f2rReal(1:atoms%jri(itype)), f2rImag(1:atoms%jri(itype)) )
            f1rReal(:) = 0.
            f2rReal(:) = 0.
            f1rImag(:) = 0.
            f2rImag(:) = 0.
            do imesh = 1, atoms%jri(itype)
              rrlr = rrl(imesh) * rmt2l ! r^l / R^(2 l +  1)
              ror = rrl(imesh) * rmtl ! r^l / R^l
              grVc0MT(imesh, lm, idir) = fpl21 * (rrl1(imesh) * f1r(imesh) - rrlr * f1r(atoms%jri(itype))     &
                & + rrl(imesh) * (f2r(atoms%jri(itype)) - f2r(imesh))) + ror * vtl(lm, idir)
            end do ! imesh
            f1r = cmplx(0., 0.)
            f2r = cmplx(0., 0.)
          end do ! mqn_m
          rrl = 0.
          rrl1 = 0.
        end do ! oqn_l
      ! .and. grRhoTermSw implicitely
      else if (.not.harSw .and. extSw) then
        ! This part can be outsourced to a seperate routine
        do oqn_l = 0, atoms%lmax(itype)
          rmtl = 1. / atoms%rmt(itype)**oqn_l ! 1 / R^l
          do mqn_m = -oqn_l, oqn_l
            lm = oqn_l * ( oqn_l + 1 ) + 1 + mqn_m
            do imesh = 1, atoms%jri(itype)
              ror = atoms%rmsh(imesh, itype)**oqn_l * rmtl ! r^l / R^l
              grVc0MT(imesh, lm, idir) =  ror * vtl(lm, idir)
            end do ! imesh
          end do ! mqn_m
        end do ! oqn_l

        ! N_check: grVc0MT(r, lm, idir)=\frac{4\pi}{2l+1}\int_{0}^{R_{\alpha}}ds_{\alpha} s_{\alpha}^{2}\frac{r_{<}^{l}}{r_{>}^{l+1}}
        !                                   grRho0MT(r, lm, iatom, idir)*[1-(\frac{r_{>}}{\bm{R}_{\alpha}})^{2l+1}]+boundary term
        !
        ! Consistent with Aaron (7.67) line 1 and (2.44) line 1 in CRG 16.12.2020-19:00

        ! Maintenance
        if (.false.) then
          do oqn_l = 0, atoms%lmax(itype)
            rmtl = 1. / atoms%rmt(itype)**oqn_l ! 1 / R^l
            do mqn_m = -oqn_l, oqn_l
              lm = oqn_l * ( oqn_l + 1 ) + 1 + mqn_m
              do imesh = 1, atoms%jri(itype)
                write(2020, '(3(i8),2(f15.8))') imesh, lm, idir, grVc0MT(imesh, lm, idir)
              end do
            end do
          end do
        end if ! maintenance

      end if
      ! .and. grRhoTermSw implicitely
      ! This has to stand here because it is valid for both ifs
      if ( extSw ) then
        do lm = 2, 4
          do imesh = 1, atoms%jri(itype)
            grVc0MT(imesh, lm, idir) = grVc0MT(imesh, lm, idir) + prefactor / atoms%rmsh(imesh, itype)**2 &
              & * ( 1 - (atoms%rmsh(imesh, itype) / atoms%rmt(itype))**3) * 3 / 4 / pi_const * c_im(idir, lm - 1)
          end do
        end do
      end if
    end do

    ! N_check: grVc0MT(r, lm, idir)+=\frac{Z_{\alpha}4\pi}{3}r_{\alpha}^{-2}*
    !                               [1-(\frac{r_{\alpha}}{\bm{R}_{\alpha}})^{3}]\frac{3}{4\pi}c_im(idir, m)\delta_{l1}+boundary term
    !
    ! Consistent with Aaron (7.48) and (2.32)/(2.33) without lm-sum in CRG 16.12.2020-19:00

    ! Maintenance
    if (.false.) then
      do idir = 1, 3
        do oqn_l = 0, atoms%lmax(itype)
          do mqn_m = -oqn_l, oqn_l
            lm = oqn_l * (oqn_l + 1) + 1 + mqn_m
            do imesh = 1, atoms%jri(itype)
              write(2027, '(3(i8),2(f15.8))') imesh, lm, idir, grVc0MT(imesh, lm, idir)
            end do ! imesh
          end do ! mqn_m
        end do ! oqn_l
      end do ! idir
    end if

  end subroutine vmtsCoul

   !--------------------------------------------------------------------------------------------------------------------------------------
   !> @author
   !> Christian-Roman Gerhorst, Forschungszentrum Jülich: IAS1 / PGI1
   !>
   !> @brief
   !> Alternative way to calculate the external potential in the interstitial according to Equation 7.43 in PhdAK
   !>
   !> @note
   !> This routine is out of order at the moment and has to be reviewed before being reactivated. But it has run once in the past.
   !--------------------------------------------------------------------------------------------------------------------------------------
subroutine psqpwVecExt(atoms, cell, gdp, ngdp, psq)
! *************************************************************************************
!  This routine calculates the Fourier coefficients of the pseudo-charge to gain the
!  gradient of the unperturbed external potential in the interstitial region. For reasons
!  of convenience, the sum over all atoms is already evaluated.
! *************************************************************************************

!#include "recycledRoutines/cpp_double.h"
use m_sphbes
use m_types

implicit none

! Variables:
!
! atoms  : atoms type defined in m_types
! cell   : unit cell type defined in m_types
! dimens : dimension type defined in m_types
! ngdp   : number of G-vectors used for the potential and the density
! gdp    : array of G-vectors used for the potential and density
! psq    : resulting 3-D pseudocharge of this routine
! itype  : runs over all types
! oqn_l  : orbital quantum number l
! df     : loop variable to calculate prefactor
! pfac   : auxiliary variable which contains 1 / Ω
! iGvec  : runs over G-vectors which are used for the potential and the density
! idirec : runs over directions the atoms can be displaced to
! iatom  : runs over all atoms
! ncvn   : contains n + l for a certain atom type, while n is taken from a table in the Weinert paper
! ineq   : runs over all equal atoms
! pn     : variable which stores prefactor with double faculties times 1 / R^l
! Gext   : G-vector in external coordinates
! sbes   : contains spherical bessel function

!     .. Scalar Arguments ..
type(t_atoms),     intent(in)    :: atoms
type(t_cell),      intent(in)    :: cell
integer,           intent(in)    :: ngdp

!     .. Array Arguments ..
integer,           intent(in)    :: gdp(:, :)
complex,           intent(out)   :: psq(:, :)

!     .. Local Scalars ..
integer                          :: itype
integer                          :: oqn_l
integer                          :: df
complex                          :: pfac
integer                          :: iGvec
integer                          :: idirec
integer                          :: iatom
integer                          :: ncvn
integer                          :: ineq
integer                          :: ufacultb
integer                          :: nc

!     .. Local Arrays ..
real                             :: pn(atoms%ntype)
real                             :: Gext(3)
real                             :: nGext
real                             :: sbes(0:atoms%lmaxd+MAXVAL(atoms%ncv)+1)

! This block calculates pn(itype) = (2 n + 5)!! in equation 7.43a PhD thesis of Aaron Klüppelberg.

!      pn = 0
!      do itype = 1, atoms%ntype
!        oqn_l = 1
!        if (oqn_l < atoms%ncv(itype)) then
!          pn(itype) = 1.
!          do nc = oqn_l, atoms%ncv(itype)
!            pn(itype) = pn(itype) * (2 * nc + 3)
!          enddo
!        end if
!      enddo
!      write (*, *) pn * 3

do itype = 1, atoms%ntype
  oqn_l = 1
  if (oqn_l >= atoms%ncv(itype)) then
    pn(itype) = 0.0
  else
    pn(itype) = 1.
    ufacultb = 2 * atoms%ncv(itype) +  3 ! + 2l is hidden in atoms%ncv, as this is n + l
    do df = ufacultb, 1, -2
      pn(itype) = pn(itype) * df
    enddo
  endif
enddo


! This block calculates equation 7.43a in the PhD thesis of Aaron Klüppelberg:
! There, the N is the bare n without the l in ncvn, l is already 1 and added to the 1 which was already there in l + n + 1
! Thus, since ncvn is constant for all l per type we only add 1 instead of two to ncvn.
  psq = cmplx(0., 0.)
  do iGvec = 1, ngdp
    Gext = matmul(cell%bmat, gdp(:, iGvec))
    nGext = norm2(Gext)
    if ( nGext == 0 ) cycle
    do idirec = 1, 3
      iatom = 1
      do itype = 1, atoms%ntype
        pfac = ImagUnit * cmplx(atoms%zatom(itype), 0) / cmplx(cell%omtil, 0)
        ncvn = atoms%ncv(itype)
        call sphbes(ncvn + 1, norm2(Gext) * atoms%rmt(itype), sbes) ! again one +1 is hidden in ncvn
        do ineq = 1, atoms%neq(itype)
          psq(iGvec, idirec) = psq(iGvec, idirec) + pfac * pn(itype)* sbes(ncvn + 1) / (nGext * atoms%rmt(itype))**(ncvn + 1)&
            &* exp(-ImagUnit * tpi_const * dot_product(gdp(:, iGvec), atoms%taual(:, iatom))) * Gext(idirec)
          iatom  = iatom + 1
        enddo
      enddo
    enddo
  enddo

end subroutine psqpwVecExt

subroutine CalcQlmHarSurf( atoms, cell, iDtype, iDatom, ngdp, gdp, rho0IRpw, rho0MTsh, qlmHartSurf )

  use m_types, only : t_atoms, t_cell

  implicit none


  ! Type parameters
  type(t_atoms),             intent(in)  :: atoms
  type(t_cell),              intent(in)  :: cell

  ! Scalar parameter
  integer,                   intent(in)  :: iDtype
  integer,                   intent(in)  :: iDatom
  integer,                   intent(in)  :: ngdp

  ! Array parameters
  integer,                   intent(in)  :: gdp(:, :)
  complex,                   intent(in)  :: rho0IRpw(:, :)
  complex,                   intent(in)  :: rho0MTsh(:, :, :, :)
  complex,                   intent(out) :: qlmHartSurf(:, :)

  ! Scalar variables
  integer                                :: oqn_l
  integer                                :: mqn_m
  integer                                :: lm
  integer                                :: idir

  ! Array variables
  complex,       allocatable             :: qlmHartSurfIR(:, :)
  complex,       allocatable             :: qlmHartSurfMT(:, :)

  qlmHartSurf(:, :) = cmplx(0., 0.)

  allocate( qlmHartSurfIR(3, (atoms%lmax(iDtype) + 1)**2) )
  allocate( qlmHartSurfMT(3, (atoms%lmax(iDtype) + 1)**2) )
  qlmHartSurfIR(:, :) = cmplx(0., 0.)
  qlmHartSurfMT(:, :) = cmplx(0., 0.)

  call CalcQlmHarSurfIR( atoms, cell, ngdp, iDtype, iDatom, gdp, rho0IRpw(:, 1), qlmHartSurfIR )
  call CalcQlmHarSurMT(atoms, iDtype, iDatom, rho0MTsh, qlmHartSurfMT)

  do oqn_l = 0, atoms%lmax(iDtype)
    do mqn_m = -oqn_l, oqn_l
      lm = oqn_l * (oqn_l + 1) + 1 + mqn_m
      do idir = 1, 3
        qlmHartSurf(lm, idir) = qlmHartSurfMT(idir, lm)
        qlmHartSurf(lm, idir) = qlmHartSurf(lm, idir) - qlmHartSurfIR(idir, lm)
        if (.false.) then
          write(2401, '(3i8,2f15.8)') oqn_l, mqn_m, idir, qlmHartSurfIR(idir, lm)
          write(2402, '(3i8,2f15.8)') oqn_l, mqn_m, idir, qlmHartSurfMT(idir, lm)
          write(2403, '(3i8,2f15.8)') oqn_l, mqn_m, idir, qlmHartSurf(lm, idir)
        end if
      end do ! idir
    end do ! mqn_m
  end do ! oqn_l

end subroutine CalcQlmHarSurf

subroutine CalcQlmHarSurfIR( atoms, cell, ngdp, iDtype, iDatom, gdp, rho0IRpw, qlmHartSurfIR )

  use m_types, only : t_atoms, t_cell
  use m_ylm, only : ylm4
  use m_sphbes, only : sphbes
  use m_gaunt, only : Gaunt1

  implicit none

  ! Type parameter
  type(t_atoms),              intent(in)  :: atoms
  type(t_cell),               intent(in)  :: cell

  ! Scalar parameter
  integer,                    intent(in)  :: ngdp
  integer,                    intent(in)  :: iDtype
  integer,                    intent(in)  :: iDatom

  ! Array parameter
  integer,                    intent(in)  :: gdp(:, :)
  complex,                    intent(in)  :: rho0IRpw(:)
  complex,                    intent(out) :: qlmHartSurfIR(:, :)

  ! Scalar variables
  integer                                 :: idir
  integer                                 :: iG
  integer                                 :: oqn_l
  integer                                 :: mqn_m
  integer                                 :: oqn_l1p
  integer                                 :: mqn_m1p
  integer                                 :: mqn_m2p
  integer                                 :: lm
  integer                                 :: lm1p
  complex                                 :: phaseFac
  complex                                 :: pref
  complex                                 :: temp1
  complex                                 :: temp2
  complex                                 :: temp3
  real                                    :: gauntFactor

  ! Array variables
  complex,       allocatable              :: ylm(:)
  real,          allocatable              :: sbes(:)
  real                                    :: gExt(3)

  ! Init assumed-shape dummy array
  qlmHartSurfIR(:, :) = cmplx(0., 0.)

  allocate( ylm( (atoms%lmax(iDtype) + 1)**2 ) )
  allocate( sbes(0:atoms%lmax(iDtype)) )
  ylm(:)   = cmplx(0., 0.)
  sbes(:) = 0.

  pref = fpi_const * atoms%rmt(iDtype) * atoms%rmt(iDtype)
  do iG = 1, ngdp
    gExt(1:3) = matmul( cell%bmat(1:3, 1:3), real(gdp(1:3, iG)) )

    ylm(:) = cmplx(0., 0.)
    call ylm4( atoms%lmax(iDtype), gExt(1:3), ylm )

    sbes(:) = 0
    call sphbes(atoms%lmax(iDtype), norm2(gExt(1:3)) * atoms%rmt(iDtype), sbes)

    phaseFac = exp( ImagUnit * tpi_const * dot_product(gdp(1:3, iG), atoms%taual(1:3, iDatom)))

    do oqn_l = 0, atoms%lmax(iDtype)
      temp1 = pref * phaseFac * atoms%rmt(iDtype)**oqn_l * rho0IRpw(iG)
      do mqn_m = -oqn_l, oqn_l
        lm = oqn_l * (oqn_l + 1) + 1 + mqn_m
        do oqn_l1p = 0, atoms%lmax(iDtype)
          temp2 = temp1 * sbes(oqn_l1p) * ImagUnit**oqn_l1p
          do mqn_m1p = -oqn_l1p, oqn_l1p
            lm1p = oqn_l1p * (oqn_l1p + 1) + 1 + mqn_m1p
            temp3 = temp2 * conjg(ylm(lm1p))
            do mqn_m2p = -1, 1
              gauntFactor = Gaunt1( oqn_l, oqn_l1p, 1, mqn_m, mqn_m1p, mqn_m2p, atoms%lmax(iDtype))
              do idir = 1, 3
                qlmHartSurfIR(idir, lm) = qlmHartSurfIR(idir, lm) + c_im(idir, mqn_m2p + 2) * gauntFactor * temp3
              end do ! idir
            end do ! mqn_mpp
          end do ! mqn_mp
        end do ! oqn_lp
      end do ! mqn_m
    end do ! oqn_l
  end do ! iG

end subroutine CalcQlmHarSurfIR

subroutine CalcQlmHarSurMT(atoms, iDtype, iDatom, rho0MTsh, qlmHartSurfMT)

  use m_types, only : t_atoms
  use m_gaunt, only : Gaunt1

  implicit none

  ! Type parameters
  type(t_atoms), intent(in)  :: atoms

  ! Scalar parameter
  integer,       intent(in)  :: iDtype
  integer,       intent(in)  :: iDatom

  ! Array parameters
  complex,       intent(in)  :: rho0MTsh(:, :, :, :)
  complex,       intent(out) :: qlmHartSurfMT(:, :)

  ! Scalar variables
  integer                    :: oqn_l
  integer                    :: oqn_l1p
  integer                    :: mqn_m1p
  integer                    :: mqn_m2p
  integer                    :: mqn_m
  integer                    :: lm
  integer                    :: lm1p
  integer                    :: idir
  real                       :: gauntFactor

  qlmHartSurfMT(:, :) = cmplx(0., 0.)

  do oqn_l = 0, atoms%lmax(iDtype)
    do mqn_m = -oqn_l, oqn_l
      lm = oqn_l * (oqn_l + 1) + 1 + mqn_m
      do oqn_l1p = 0, atoms%lmax(iDtype)
        do mqn_m1p = -oqn_l1p, oqn_l1p
          lm1p = oqn_l1p * (oqn_l1p + 1) + 1 + mqn_m1p
          do mqn_m2p = -1, 1
            gauntFactor = Gaunt1( oqn_l, oqn_l1p, 1, mqn_m, mqn_m1p, mqn_m2p, atoms%lmax(iDtype))
            do idir = 1, 3
              qlmHartSurfMT(idir, lm) = qlmHartSurfMT(idir, lm) + c_im(idir, mqn_m2p + 2) * atoms%rmt(iDtype)**(oqn_l + 2)     &
                & * rho0MTsh(atoms%jri(iDtype), lm1p, iDatom, 1) * gauntFactor
            end do ! idir
          end do ! mqn_m2p
        end do ! mqn_m1p
      end do ! oqn_l1p
    end do ! mqn_m
  end do ! oqn_l

end subroutine CalcQlmHarSurMT

subroutine phasy1nSym(atoms, cell, Gvec, qptn, pylm)

  use m_ylm
  use m_types

  implicit none

  ! Scalar Type Arguments
  type(t_atoms),  intent(in)  ::  atoms
  type(t_cell),   intent(in)  ::  cell

  ! Array Arguments
  integer,        intent(in)  ::  Gvec(:)
  real,           intent(in)  ::  qptn(:)
  complex,        intent(out) ::  pylm(:, :)

  !-------------------------------------------------------------------------------------------------------------------------------
  ! Local Scalar Variables
  ! iatom : runs over all atoms
  ! itype : runs over all types
  ! ineq  : runs over all equivalent atoms of one atom type
  ! lm    : encodes oqn_l and mqn_m
  ! sf    : stores exponential function
  ! csf   : stores exponential function times 4 π i^l
  ! x     : stores argument of exponential function
  ! mqn_m : magnetic quantum number m
  ! oqn_l : orbital quantum number l
  ! ll1   : auxiliary variable to calculate lm
  !-------------------------------------------------------------------------------------------------------------------------------
  integer                     ::  iatom
  integer                     ::  itype
  integer                     ::  ineq
  integer                     ::  lm
  complex                     ::  sf
  complex                     ::  csf
  real                        ::  x
  integer                     ::  mqn_m
  integer                     ::  oqn_l
  integer                     ::  ll1

  !-------------------------------------------------------------------------------------------------------------------------------
  ! Local Array Variables
  ! fpiul: stores 4 π i^l
  ! Gqext: stores G + q in external coordinates
  ! ylm  : stores Y_lm
  !-------------------------------------------------------------------------------------------------------------------------------
  complex, allocatable        ::  fpiul(:)
  complex, allocatable        ::  ylm(:)
  real                        ::  Gqext(3)

  allocate(fpiul(0:atoms%lmaxd))
  fpiul(:) = cmplx(0., 0.)
  ! calculates 4 π i^l resolved for every l, not divided by nop because no loop over symmetry operations
  do oqn_l = 0, atoms%lmaxd
     fpiul(oqn_l) = fpi_const * ImagUnit**oqn_l
  enddo


  ! calculates Y*_lm(\vec{G} + \vec{q}) for every l and m. The argument Gqext must be in external coordinates.
  allocate(ylm((atoms%lmaxd + 1)**2))
  ylm(:) = cmplx(0., 0.)
  Gqext(1:3) = matmul(cell%bmat(1:3, 1:3), real(Gvec(1:3) + qptn(1:3)))
  call ylm4(atoms%lmaxd, Gqext(1:3), ylm)
  ylm = conjg(ylm)


  ! calculates first exp(i (G + q) tau)  and multiplies recent factors before storing the final result to pylm
  iatom = 1
  pylm(:, :) = cmplx(0.,0.)
  do itype = 1, atoms%ntype
     do ineq = 1, atoms%neq(itype)
        x = tpi_const * dot_product(real(Gvec(1:3) + qptn(1:3)), atoms%taual(1:3, iatom))
        sf = exp(ImagUnit *  x)
        do oqn_l = 0, atoms%lmax(itype)
           ll1 = oqn_l * (oqn_l + 1) + 1
           csf = fpiul(oqn_l) * sf
           do mqn_m = -oqn_l, oqn_l
              lm = ll1 + mqn_m
              pylm(lm, iatom) = csf * ylm(lm)
           enddo ! mqn_m
        enddo ! oqn_l
        iatom = iatom + 1
     enddo ! ineq
  enddo ! itype

end subroutine phasy1nSym

subroutine convolMTgrVeff0dKern(atoms, grRho0MT, dKernMTGPts, gWghts, ylm, grVxc0MT)

  use m_types, only : t_atoms

  implicit none

  ! Type parameter
  type(t_atoms),     intent(in)  :: atoms

  ! Array parameter
  complex,           intent(in)  :: grRho0MT(:, :, :, :)
  complex,           intent(in)  :: ylm(:, :)
  real,              intent(in)  :: gWghts(:) ! gaussian weights belonging to gausPts
  real,              intent(in)  :: dKernMTGPts(:, :, :)
  complex,           intent(out) :: grVxc0MT(:, :, :, :)

  ! Local scalar variables
  integer                        :: idir
  integer                        :: iatom
  integer                        :: itype
  integer                        :: ieqat
  integer                        :: igmesh ! Loop variable over sampling points of spherical Gauss mesh
  integer                        :: irmesh ! Loop variable over radial MT mesh
  integer                        :: oqn_l
  integer                        :: lm_lonly !reduce multiplication when calculating lm
  integer                        :: mqn_m
  integer                        :: lm
  complex                        :: grVxcMTKernAdd

  ! Local allocatable variables
  complex, allocatable           :: grRhoMTGpts(:, :) !grRhoMT on Gauss mesh
  real, allocatable              :: grVxcMTKernGPts(:, :)

  grVxc0MT(:, :, :, :) = cmplx(0., 0.)

  allocate( grRhoMTGpts(atoms%nsp(), atoms%jmtd), grVxcMTKernGPts(atoms%nsp(), atoms%jmtd) )
  grRhoMTGpts(:, :) = cmplx(0., 0.)
  grVxcMTKernGPts(:, :) = 0.

  do idir = 1, 3
    iatom = 0
    do itype = 1, atoms%ntype
      do ieqat = 1, atoms%neq(itype)
        iatom = iatom + 1
        grRhoMTGpts(:, :) = cmplx(0., 0.)
        grVxcMTKernGPts(:, :) = 0.
        ! Density's gradient has l + 1 entries
        do oqn_l = 0, atoms%lmax(itype)
          lm_lonly = oqn_l * ( oqn_l + 1 ) + 1
          do mqn_m = -oqn_l, oqn_l
            lm = lm_lOnly + mqn_m
            ! Evaluate grRho on spherical Gauss mesh in order to apply Gauss quadrature.
            do irmesh = 1, atoms%jri(itype)
              do igmesh = 1, atoms%nsp()
                grRhoMTGpts(igmesh, irmesh) = grRhoMTGpts(igmesh, irmesh) + grRho0MT(irmesh, lm, iatom, idir) * ylm(igmesh, lm)
              end do ! igmesh
            end do ! irmesh
          end do ! mqn_m
        end do ! oqn_l
#if DEBUG_MODE
        if ( any( aimag( grRhoMTGpts ) > 1e-7 ) ) then
          write(*, *) 'Warning rhoMTGpts has imaginary components.'
        end if
#endif
        ! On the spherical Gauss mesh the integral reduces to a weighted (gWghts) sum (over all sampling points on the Gauss mesh)
        ! of the MT exchange-correlation kernel and either the density's gradient or the first variation of the gradient.
        do irmesh = 1, atoms%jri(itype)
          do igmesh = 1, atoms%nsp()
            grVxcMTKernGPts(igmesh, irmesh) = grVxcMTKernGpts(igmesh, irmesh) + real(grRhoMTGpts(igmesh, irmesh)) &
                                            & * dKernMTGPts(igmesh, irmesh, iatom) * gWghts(igmesh)
          end do
        end do
        do oqn_l = 0, atoms%lmax(itype)
          lm_lonly = oqn_l * ( oqn_l + 1 ) + 1
          do mqn_m = -oqn_l, oqn_l
            lm = lm_lOnly + mqn_m
            do irmesh = 1, atoms%jri(itype)
            ! Back-transformation of the MT coefficients. Now they are expansion coefficients of the MT grid.
              grVxcMTKernAdd = dot_product( grVxcMTKernGPts(:atoms%nsp(), irmesh), conjg(ylm(: atoms%nsp(), lm)) )
            ! Add this contribution to MT exchange-correlation contribution to the potential
              grVxc0MT(irmesh, lm, iatom, idir) = grVxc0MT(irmesh, lm, iatom, idir) + grVxcMTKernAdd
            end do ! irmesh
          end do ! mqn_m
        end do ! oqn_l
      end do ! ieqat
    end do ! itype
  end do ! itype

end subroutine convolMTgrVeff0dKern

subroutine convolGrRhoKern(stars, ngdp, ngpqdp, f1IR, f2IR, pdG2FouM, pdG2FouMv, f3IR, iqpt)

  use m_cfft
  use m_types

  implicit none

  ! Type parameter
  ! -------------
  type(t_stars),        intent(in)  :: stars

  ! Scalar parameter
  integer,              intent(in)  :: ngdp
  integer,              intent(in)  :: ngpqdp
  integer,              intent(in)  :: iqpt

  ! Array parameter
  !----------------
  ! 2 quantities to convolute
  complex,              intent(in)  :: f1IR(:)
  complex,              intent(in)  :: f2IR(:)
  ! maps from a G to respective mesh entry of FFT mesh
  integer,              intent(in)  :: pdG2FouM(:)
  integer,              intent(in)  :: pdG2FouMv(:)
  ! Convoluted quantity
  complex,              intent(out) :: f3IR(:)

  ! Maximal length of FFT mesh
  integer                           :: ifftd
  ! Loop index
  integer                           :: iG
  ! Scaling factor for canceling FFT artefacts.
  real                              :: scaling
  integer                           :: ii

  ! Real and imaginary parts of quantities to be convoluted (1, 2) and convoluted quantity 3
  real,    allocatable              :: rf1(:)
  real,    allocatable              :: if1(:)
  real,    allocatable              :: rf2(:)
  real,    allocatable              :: if2(:)
  real,    allocatable              :: rf3(:)
  real,    allocatable              :: if3(:)


  ! Length of FFT mesh. The cartesian dimensions are stored sequentially. The G's expand maximally from -k_i to k_i (from this
  ! cube a ball of radius gmax is cut off) so 3 * k_i should be in principle large enough for all given G's components. And we
  ! have some free space to avoid aliasing. However, this factor of 3 is only working well for symmetrized quantities which is
  ! why we have to expand the mx1, mx2, mx3 if we expand our quantities in plane waves.
  ifftd = 27 * stars%mx1 * stars%mx2 * stars%mx3

  ! All quantities have the size of the FFT mesh.
  allocate( rf1(0:ifftd - 1), if1(0:ifftd - 1), rf2(0:ifftd - 1), if2(0:ifftd - 1), rf3(0:ifftd - 1), if3(0:ifftd - 1) )

  rf1(:) = 0
  rf2(:) = 0
  rf3(:) = 0
  if1(:) = 0
  if2(:) = 0
  if3(:) = 0

  ! Extract the real and imaginary part of the expansion coefficients of the quantities to be convoluted for every G-vector and
  ! map them with pdG2FouM to its respective mesh point.
  do iG = 1, ngpqdp
    rf1(pdG2FouMv(iG)) = real(f1IR(iG))
    if1(pdG2FouMv(iG)) = aimag(f1IR(iG))
  end do

  do iG = 1, ngdp ! is it kimax - 1 or kimax?
    rf2(pdG2FouM(iG)) = real(f2IR(iG))
    if2(pdG2FouM(iG)) = aimag(f2IR(iG))
  end do

    ! Complex FFT of 1st quantity into real space, this is done as it is done in FLEUR fft3d
    call Cfft(rf1, if1, ifftd, 3 * stars%mx1, 3 * stars%mx1, 1)
    call Cfft(rf1, if1, ifftd, 3 * stars%mx2, 9 * stars%mx1 * stars%mx2, 1)
    call Cfft(rf1, if1, ifftd, 3 * stars%mx3, ifftd, 1)

    ! Complex FFT of 2nd quantity into real space
    call Cfft(rf2, if2, ifftd, 3 * stars%mx1, 3 * stars%mx1, 1)
    call Cfft(rf2, if2, ifftd, 3 * stars%mx2, 9 * stars%mx1 * stars%mx2, 1)
    call Cfft(rf2, if2, ifftd, 3 * stars%mx3, ifftd, 1)

    ! Exploiting the convolution theorem
    do ii = 0, ifftd - 1
      rf3(ii) = real( cmplx(rf1(ii),if1(ii)) * cmplx(rf2(ii), if2(ii)))
      if3(ii) = aimag( cmplx(rf1(ii),if1(ii)) * cmplx(rf2(ii), if2(ii)))
    end do

# ifdef DEBUG_MODE
    ! It is sufficient to test it only for iqpt == 1 because the gradient of rho is calculated outside the q-loop and the density
    ! variation has only to be real for q = 0
    if(iqpt == 1) then
      if ( any(abs(if3) > 1e-7 )) then
        write(*, *) 'Convolution FFT has imaginary components'
        !NOstopNO'Convolution FFT has imaginary components'
      end if
      if3 = 0
    end if
#endif

    ! Complex FFT of convoluted quantity into reciprocal space
    call Cfft(rf3, if3, ifftd, 3 * stars%mx1, 3 * stars%mx1, -1)
    call Cfft(rf3, if3, ifftd, 3 * stars%mx2, 9 * stars%mx1 * stars%mx2, -1)
    call Cfft(rf3, if3, ifftd, 3 * stars%mx3, ifftd, -1)

    ! We have to care for the artefacts of this FFT
    scaling = 1. / ifftd
    f3IR = cmplx( 0.0, 0.0 )
    ! Map convoluted quantity from FFT mesh to plane-wave expansion coefficient representation.
    do iG = 1, ngpqdp !kimax is max G-vector
      f3IR(iG) = f3IR(iG) + cmplx(rf3(pdG2FouMv(iG)), if3(pdG2FouMv(iG))) * scaling
    end do

end subroutine convolGrRhoKern

end module m_jpGrVeff0
