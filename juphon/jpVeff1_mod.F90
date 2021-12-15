!-----------------------------------------------------------------------------------------------------------------------------------
! Forschungszentrum Jülich, juPhon Plugin for the FLEUR program
!-----------------------------------------------------------------------------------------------------------------------------------
!
! MODULE: Calculation of first variation of effective potential
!
!> @author
!> Christian-Roman Gerhorst
!
!> @brief
!> This module contains routines related to the calculation of the first variation of the effective potential.
!>
!> @todo
!> Complete documentation
!>
!> @note
!> Additional information and formulas pointing out the routines of this module can be found within this
!> <a href='jpVcoul1.pdf'>document</a>.
!-----------------------------------------------------------------------------------------------------------------------------------
module m_jpVeff1

  implicit none

  contains

  !---------------------------------------------------------------------------------------------------------------------------------
  !> @author
  !> Christian-Roman Gerhorst, Forschungszentrum Jülich: IAS1 / PGI1
  !>
  !> @brief
  !> Main routine which calculates the first variation of the effective potential calling other helping subroutines.
  !>
  !> @details
  !> The first varation of the Coulomb potential is calculated using the Weinert method without using symmetry by not exploiting
  !> stars and lattice harmonics according to M.Weinert, J.Math.Phys. 22(11) (1981) p.2434 and the PhD thesis of Aaron Klüppelberg.
  !> For the xc-potential, the functional derivative of the kernel with respect to the density is connected with the linear
  !> variation of the density. In the interstitial region, this is done by a convolution, whereas in the muffin-tin a projection
  !> onto a spherical harmonic is evaluted by using Gauss quadrature.
  !>
  !> @note
  !> At the moment, only the x-alpha kernel functional derivative is implemented, but its calcualtion is in seperate routines, that
  !> can be replaced by the libxc library later.
  !>
  !> @param[in]   stars       : Stars type, see types.f90
  !> @param[in]   cell        : Unit cell type, see types.f90
  !> @param[in]   atoms       : Atoms type, see types.f90
  !> @param[in]   dimens      : Dimension type, see types.f90
  !> @param[in]   harSw       : Switch to enable Hartree contribution
  !> @param[in]   extSw       : Switch to enable external contribution
  !> @param[in]   xcSw        : Switch to enable xc contribution
  !> @param[in]   VextFull    : Switch to activate volume term contribution of muffin-tin boundary problem
  !> @param[in]   ngdp        : Number of G-vectors for potentials and densities
  !> @param[in]   iDatom      : Index of displaced atom
  !> @param[in]   iDtype      : Index of atom type to which displaced atom belongs
  !> @param[in]   iqpt        : Index of q-point in kpts file that is evaluated
  !> @param[in]   vHarNum     : Switch to enable special terms that are needed to optimize the integrals rho1 Veff1 and rho1 Vext1
  !> @param[in]   ngpqdp      : Number of G-vectors for shifted G-set for q-point with index iqpt
  !> @param[in]   qpoint      : Qpoint in internal coordinates as in qpts%bk(1:3, iqpt)
  !> @param[in]   rho0IRpw    : Plane-wave coefficients of the unperturbed and converged interstitial density parsed from Fleur
  !> @param[in]   rho0MTsh    : Spherical harmonic coefficients of unperturbed and converged muffin-tin density parsed from Fleur
  !> @param[in]   rho1IR      : Plane-wave interstitial expansion coefficients of the linear density variation
  !> @param[in]   rho1MT      : Spherical harmonic muffin-tin expansion coefficients of the linear density variation without the
  !>                            muffin-tin gradient of the density
  !> @param[in]   grRho0MT    : Spherical harmonic coefficients of the gradient of the muffin-tin unperturbed density
  !> @param[in]   gdp         : G-vectors of potentials and densities
  !> @param[in]   ylm         : Set of spherical harmonics whose arguments are the unit vectors of the Gauss mesh points up to lmax
  !>                            created in m_jppotdenshelper::calcmtdvxckern
  !> @param[in]   gWghts      : Gauss weights for Gauss quadrature created in m_jppotdenshelper::calcmtdvxckern
  !> @param[in]   dKernMTGPts : Spherical harmonic muffin-tin coefficients of the functional derivative of the xc-kernel with
  !>                            respect to the density created in m_jppotdenshelper::calcmtdvxckern
  !> @param[in]   vxc1IRKern  : Plane-wave interstitial coefficients of the functional derivative of the xc-kernel with respect to
  !>                            the density
  !> @param[in]   gpqdp       : G-vectors for shifted G-set for q-point with index iqpt
  !> @param[out]  vEff1IR     : Plane-wave insterstitial coefficients of the gradient of the unperturbed effective potential
  !> @param[out]  vEff1MT     : Spherical harmonic muffin-tin coefficients of the gradient of the unperturbed effective potential
  !---------------------------------------------------------------------------------------------------------------------------------
  subroutine GenVeff1(stars, cell, atoms, dimens, harSw, extSw, xcSW, VextFull, ngdp, qpoint, rho0IRpw, rho0MTsh, rho1IR, rho1MT, &
      & grRho0MT, gdp, vEff1IR, vEff1MT, vxc1IRKern, ylm, dKernMTGPts, gWghts, iDatom, iDtype, iqpt, ngpqdp, gpqdp, vHarNum)

    use m_types, only : t_stars, t_cell, t_atoms, t_dimension
    use m_JPConstants, only : fpi, iu, compPhon
    use m_jpPotDensHelper, only : convolGrRhoKern, convolMTVeff1dKern
    use m_juDFT_NOstopNO, only : juDFT_warn

    implicit none

    ! Variables:
    !
    ! stars        : star type, see m_types for documentation
    ! cell         : unit cell type, see m_types for documentation
    ! atoms        : atoms type, see m_types for documentation
    ! dimens       : dimension type, see m_types for documentation
    ! rho0IR       : density in interstitial region
    ! rho0MT       : density in MT region
    ! qpoint       : current q-point
    ! rho1PW       : first variation of the density in the interstitial region
    ! vEff1IR     : Coulomb potential in the interstitial region
    ! vEff1MT     : Coulomb potential in the muffin-tin region
    ! qlmPhon1Coul : multipole moment for calculation of the IR part of the first variation of the Coulomb potential
    ! qlmGrad0Hart : multipole moment for calculation of the IR part of the gradient part of the unperturbed Hartree potential
    ! psqPhon1Coul : pseudo charge constructed from multipole moments qlmPhon1Coul
    ! psqGrad0Hart : pseudo charge constructed from qlmGrad0Hartr
    ! psqGrad0Ext  : pseudo charge to calculate the interstital contribution of the gradient of the unperturbed external potential
    ! Gint         : stores the G-vector in internal coordinates
    ! Gext         : stores the G-vector in cartesian coordinates
    ! Gqext        : stores G + q in cartesian coordinates
    ! gdp          : stores all G-vectors used for the potential and the density
    ! gdptemp      : auxiliary variable to determine G-vectors used for the potential and the density
    ! rho0HPW      : plane wave contribution to the hartree pseudocharge
    ! vGradHart0IR : gradient of the unperturbed Hartree potential in the interstitial region
    ! vGradExt0IR  : gradient of the unperturbed external potential in the interstitial region
    ! Gx           : auxiliary variable to determine required G-vectors
    ! Gy           : auxiliary variable to determine required G-vectors
    ! Gz           : auxiliary variable to determine required G-vectors
    ! ngdp         : number of G-vectors needed for the potential and the density
    ! iGvec        : runs over the G-vectors
    ! idirec       : runs over the directions the atoms can be displaced to
    ! iatom        : runs over all atoms
    ! itype        : runs over all types
    ! ineq         : runs over equivalent atoms of an atom type
    ! symType      : stores point symmetry of a type
    ! lh           : auxiliary variable to store lattice harmonic index
    ! mem          : runs over members of a lattice harmonic
    ! mems         : stores total number of members for a lattice harmonic
    ! oqn_l        : orbital quantum number l
    ! lm           : encodes l and m
    ! mqn_m        : magnetic quantum number m

    ! Type Parameter
    type(t_stars),                intent(in)    :: stars
    type(t_cell),                 intent(in)    :: cell
    type(t_atoms),                intent(in)    :: atoms
    type(t_dimension),            intent(in)    :: dimens

    ! Scalar Parameter
    logical,                      intent(in)    :: harSw
    logical,                      intent(in)    :: extSw
    logical,                      intent(in)    :: xcSw
    logical,                      intent(in)    :: VextFull
    integer,                      intent(in)    :: ngdp
    integer,                      intent(in)    :: iDatom
    integer,                      intent(in)    :: iDtype
    integer,                      intent(in)    :: iqpt
    logical,                      intent(in)    :: vHarNum
    integer,                      intent(in)    :: ngpqdp

    ! Array Parameter
    real,                         intent(in)    :: qpoint(:)
    complex,                      intent(in)    :: rho0IRpw(:, :)
    complex,                      intent(in)    :: rho0MTsh(:, :, :, :)
    complex,                      intent(in)    :: rho1IR(:, :)
    complex,                      intent(in)    :: rho1MT(:, :, :, :)
    complex,                      intent(in)    :: grRho0MT(:, :, :, :)
    integer,                      intent(in)    :: gdp(:,:)
    complex,                      intent(in)    :: ylm(:, :)
    real,                         intent(in)    :: gWghts(:) ! gaussian weights belonging to gausPts
    real,                         intent(in)    :: dKernMTGPts(:, :, :)
    complex,                      intent(in)    :: vxc1IRKern(:)
    integer,                      intent(in)    :: gpqdp(:, :)
    complex,        allocatable,  intent(out)   :: vEff1IR(:, :)
    complex,        allocatable,  intent(out)   :: vEff1MT(:, :, :, :)

    ! Local scalars
    integer                                     :: iG
    integer                                     :: iatom
    integer                                     :: itype
    integer                                     :: ineq
    integer                                     :: oqn_l
    integer                                     :: lm
    integer                                     :: mqn_m
    real                                        :: noGqext
    integer                                     :: nfftx
    integer                                     :: nffty
    integer                                     :: nfftz
    integer                                     :: nfftxy
    integer                                     :: GxFFT
    integer                                     :: GyFFT
    integer                                     :: GzFFT
    integer                                     :: imesh
    integer                                     :: idir
    integer                                     :: ieqat
    logical                                     :: linIntp

    ! Local arrays
    complex,        allocatable                 :: qlmPhon1Coul(:,:,:)
    complex,        allocatable                 :: psqPhon1Coul(:, :)
    complex,       allocatable                  :: psqPhon1CoulVar(:, :)
    integer,        allocatable                 :: pdG2FouM(:)
    integer,        allocatable                 :: pdG2FouMv(:)
    complex,        allocatable                 :: Vxc1MT(:, :, :, :)
    complex,        allocatable                 :: Vxc1IR(:, :)
    complex,        allocatable                 :: rho1IRzero(:, :)
    complex,        allocatable                 :: rho1MTzero(:, :, :, :)
    complex,        allocatable                 :: rho1MTfull(:, :, :, :)
    real                                        :: Gqext(3)

    allocate(vEff1MT(atoms%jmtd, (atoms%lmaxd + 1)**2, atoms%nat, 3))

    ! This has to be initialized here because if only xc potential calculated. It is not initiliazed within vmts
    vEff1MT(:, :, :, :) = cmplx(0., 0.)

    ! create multipole moments of the original charge to calculate in the interstitial region the first variation of the Coulomb potential
    ! and the gradients of the unperturbed Hartree and external potentials.
    ! the switches harSw, and extSw are passed to mpmomVc1 to manipulate the multipole moments. If only external potential is needed
    ! in the first iteration it is sufficient that the multipole moments are build such that only the displaced atom has a contribution and
    ! only l = 1 is relevant.
    call mpmomVc1(atoms, cell, qpoint, gpqdp, harSw, extSw, ngpqdp, ngdp, gdp, rho0IRpw, rho0MTsh, grRho0MT, rho1IR, rho1MT, qlmPhon1Coul, iDatom, iDtype)


    ! This block calculates the interstitial contribution of the first variation of the Coulomb potential. First the multipole moments
    ! qlmPhon1Coul are used to generate the respective pseudo charges. Then, the interstitial contribution of the 1st variation Coulomb
    ! potential is evaluated using the solution of the Poisson equation.
    ! Only exteranl potential

    allocate(psqPhon1Coul(3, ngpqdp))
    psqPhon1Coul(:, :) = cmplx(0., 0.)
    if(.not.harSw .and. extSw) then
      ! In case of only the external potential we still use the Weinert way of calculating the pseudodensity, however for the
      ! external potential there is no contribution by the interstitial density as all the sources lie within the muffin-tins.
      allocate(rho1IRzero(ngpqdp, 3))
      rho1IRzero(:, :) = cmplx(0.0, 0.0)
      call psqpwVec(atoms, cell, dimens, gpqdp, qpoint, ngpqdp, rho1IRzero, qlmPhon1Coul, psqPhon1Coul)
    else if (harSw) then
      call psqpwVec(atoms, cell, dimens, gpqdp, qpoint, ngpqdp, rho1IR, qlmPhon1Coul, psqPhon1Coul)
    end if

    if (.false.) then
      allocate(psqPhon1CoulVar(3, ngpqdp))
      psqPhon1CoulVar = 0
      ! Note: rho1IR migh also be rho1IRzero, which then has to be initialized
      call psqpwVecVar(atoms, cell, dimens, gpqdp, qpoint, ngpqdp, rho1IR, qlmPhon1Coul, psqPhon1CoulVar, iDatom)

      do iG = 1, ngpqdp
        do idir = 1, 3
          write(2185, '(5(i8),1x,2(f30.15))') iG, idir, gpqdp(1:3, iG), psqPhon1Coul(idir, iG)
          write(1020, '(5(i8),1x,2(f30.15))') iG, idir, gpqdp(1:3, iG), psqPhon1CoulVar(idir, iG)
        end do ! idir
      end do ! iG
    end if

    allocate(vEff1IR(ngpqdp, 3))
    vEff1IR(:, :) = cmplx(0., 0.)
    do iG = 1, ngpqdp
      Gqext(1:3) = matmul(cell%bmat(1:3, 1:3), gpqdp(1:3, iG) + qpoint(1:3))
      noGqext = norm2(Gqext(1:3))
      if (noGqext < 1E-9) then
        !if (compPhon) then
          write(110,*) iG
          write(110,*) 0.0,0.0
          write(113,*) iG
          write(113,*) 0.0,0.0
          write(114,*) iG
          write(114,*) 0.0,0.0
          write(111,*) iG
          write(111,*) gpqdp(1,iG), gpqdp(2, iG),gpqdp(3, iG)
        !end if

        cycle
      end if
      do idir = 1, 3
        vEff1IR(iG, idir) = fpi * psqPhon1Coul(idir, iG) / noGqext**2

        !if (compPhon) then
          if (idir.eq.3) then
            write(110,*) iG
            write(110,*) real(vEff1IR(iG, 1)), aimag(vEff1IR(iG, 1))
            write(113,*) iG
            write(113,*) real(vEff1IR(iG, 2)), aimag(vEff1IR(iG, 2))
            write(114,*) iG
            write(114,*) real(vEff1IR(iG, 3)), aimag(vEff1IR(iG, 3))
            write(111,*) iG
            write(111,*) gpqdp(1,iG), gpqdp(2, iG),gpqdp(3, iG)
          end if
        !end if

      end do
    end do


    ! This block solves the Dirichlet boundary value problem in order to get the MT contribution of the first variation of the Coulomb
    ! potential plus the gradient of the unperturbed Coulomb potential for displaced atoms. In this context, it is considered, whether
    ! the current atom is displaced or not.
    ! This has to be made here because in mpmom also the gradient is summed to the density sothat it would be doubled if we would
    ! give a decorated rho from outside.
    if (vExtFull .and. (harSw .or. xcSw)) then
      allocate( rho1MTfull(atoms%jmtd, (atoms%lmaxd + 1)**2, atoms%nat, 3) )
      rho1MTfull(:, :, :, :) = cmplx(0.0, 0.0)
      rho1MTfull(1:atoms%jmtd, 1:(atoms%lmaxd + 1)**2, 1:atoms%nat, 1:3) = rho1MT(1:atoms%jmtd, 1:(atoms%lmaxd + 1)**2, 1:atoms%nat, 1:3)
      do idir = 1, 3
        do oqn_l = 0, atoms%lmax(iDtype)
          do mqn_m = -oqn_l, oqn_l
            lm = oqn_l * (oqn_l + 1) + mqn_m + 1
            do imesh = 1, atoms%jri(iDtype)
              rho1MTfull(imesh, lm, iDatom, idir) = rho1MTfull(imesh, lm, iDatom, idir) - grRho0MT(imesh, lm, iDatom, idir)
            end do
          end do
        end do
      end do
    end if

    if (.not.harSw .and. extSw) then
      allocate(rho1MTzero(atoms%jmtd, atoms%lmaxd * (atoms%lmaxd + 1)**2, atoms%nat, 3))
      rho1MTzero = cmplx(0.0, 0.0)
    end if
    linIntp = .true.
    do idir = 1, 3
      iatom = 0
      do itype = 1, atoms%ntype
        do ineq = 1, atoms%neq(itype)
          iatom = iatom + 1
          if (.not.harSw .and. extSw) then
            call vmtsNSym(atoms, cell, ngpqdp, idir, iatom, itype, iDatom, vExtFull, gpqdp, qpoint, vEff1IR, rho1MTzero, vEff1MT(:, :, iatom, idir), harSw, extSw, vHarNum, linIntp )
          else if (harSw) then
            if ( vExtFull ) then
              call vmtsNSym(atoms, cell, ngpqdp, idir, iatom, itype, iDatom, vExtFull, gpqdp, qpoint, vEff1IR, rho1MTfull, vEff1MT(:, :, iatom, idir), harSw, extSw, vHarNum, linIntp)
            else
              call vmtsNSym(atoms, cell, ngpqdp, idir, iatom, itype, iDatom, vExtFull, gpqdp, qpoint, vEff1IR, rho1MT, vEff1MT(:, :, iatom, idir), harSw, extSw, vHarNum, linIntp)
            end if
          end if
        end do
      end do
    end do

    if (compPhon) then
      if (.not.harSw .and. extSw) then
        do idir = 1, 3
          iatom = 0
          do itype = 1, atoms%ntype
            do ieqat = 1, atoms%neq(itype)
              iatom = iatom + 1
              do lm = 1, (atoms%lmax(itype) + 1)**2
                do imesh = 1, atoms%jri(itype)
                  if (idir.eq.3) then
                    write(109,*) imesh, lm, iatom, 1
                    write(109,*) real(vEff1MT(imesh, lm, iatom, 1)), aimag(vEff1MT(imesh, lm, iatom, 1))
                    write(109,*) imesh, lm, iatom, 2
                    write(109,*) real(vEff1MT(imesh, lm, iatom, 2)), aimag(vEff1MT(imesh, lm, iatom, 2))
                    write(109,*) imesh, lm, iatom, 3
                    write(109,*) real(vEff1MT(imesh, lm, iatom, 3)), aimag(vEff1MT(imesh, lm, iatom, 3))
                  end if
                end do
              end do
            end do
          end do
        end do
      end if
    end if

    if ( xcSw ) then
      ! Construct mapping array for mapping a G-vector to its respective mesh point on the FFT mesh.
      nfftx = 3 * stars%k1d
      nffty = 3 * stars%k2d
      nfftz = 3 * stars%k3d
      nfftxy= 9 * stars%k1d * stars%k2d

      ! Mapping array for kernel
      allocate(pdG2FouM(ngdp))
      pdG2FouM(:) = 0
      do iG = 1, ngdp
        GxFFT = gdp(1, iG)
        GyFFT = gdp(2, iG)
        GzFFT = gdp(3, iG)
        if (GxFFT < 0) GxFFT = GxFFT + nfftx
        if (GyFFT < 0) GyFFT = GyFFT + nffty
        if (GzFFT < 0) GzFFT = GzFFT + nfftz
        pdG2FouM(iG) = GxFFT + GyFFT * nfftx + GzFFT * nfftxy
      end do

      ! Mapping array vor kernel
      allocate(pdG2FouMv(ngpqdp))
      pdG2FouMv = 0
      do iG = 1, ngpqdp
        GxFFT = gpqdp(1, iG)
        GyFFT = gpqdp(2, iG)
        GzFFT = gpqdp(3, iG)
        if (GxFFT < 0) GxFFT = GxFFT + nfftx
        if (GyFFT < 0) GyFFT = GyFFT + nffty
        if (GzFFT < 0) GzFFT = GzFFT + nfftz
        pdG2FouMv(iG) = GxFFT + GyFFT * nfftx + GzFFT * nfftxy
      end do

      ! Convolute rho1IR with functional derivative of kernel to get Vxc within the IR.
      allocate(vxc1IR(ngpqdp, 3))
      vxc1IR(:, :) = 0
      do idir = 1, 3
        call convolGrRhoKern(stars, ngdp, ngpqdp, rho1IR(:, idir), vxc1IRKern, pdG2FouM, pdG2FouMv, vxc1IR(:, idir), iqpt)
      end do

      ! Add x-alpha xc-potential contribution for IR
      do idir = 1, 3
        do iG = 1, ngpqdp
          vEff1IR(iG, idir) = vEff1IR(iG, idir) + vxc1IR(iG, idir)

          if (compPhon) then
            if (idir.eq.3) then
              write(112,*) iG, 1
              write(112,*) real(vEff1IR(iG, 1)), aimag(vEff1IR(iG, 1))
              write(112,*) iG, 2
              write(112,*) real(vEff1IR(iG, 2)), aimag(vEff1IR(iG, 2))
              write(112,*) iG, 3
              write(112,*) real(vEff1IR(iG, 3)), aimag(vEff1IR(iG, 3))
            end if
          end if

        end do ! iG
      end do ! idir
      deallocate(vxc1IR)

      ! Add x-alpha xc-potential contribution for MT
      allocate( Vxc1MT(atoms%jmtd, ( atoms%lmaxd + 1)**2, atoms%nat, 3) )
      Vxc1MT(:, :, :, :) = 0
      iatom = 0
      do itype = 1, atoms%ntype
        do ieqat = 1, atoms%neq(itype)
          iatom = iatom + 1
          if (vExtFull) then
            call convolMTVeff1dKern(atoms, dimens, iqpt, rho1MTfull, gWghts, dKernMTGPts, ylm, Vxc1MT)
          else
            call convolMTVeff1dKern(atoms, dimens, iqpt, rho1MT, gWghts, dKernMTGPts, ylm, Vxc1MT)
          end if
        end do
      end do

      do idir = 1, 3
        iatom = 0
        do itype = 1, atoms%ntype
          do ieqat = 1, atoms%neq(itype)
            iatom = iatom + 1
            do lm = 1, (atoms%lmax(itype) + 1)**2
              do imesh = 1, atoms%jri(itype)
                vEff1MT(imesh, lm, iatom, idir) = vEff1MT(imesh, lm, iatom, idir) + Vxc1MT(imesh, lm, iatom, idir)

                if (compPhon) then
                  if (idir.eq.3) then
                    write(109,*) imesh, lm, iatom, 1
                    write(109,*) real(vEff1MT(imesh, lm, iatom, 1)), aimag(vEff1MT(imesh, lm, iatom, 1))
                    write(109,*) imesh, lm, iatom, 2
                    write(109,*) real(vEff1MT(imesh, lm, iatom, 2)), aimag(vEff1MT(imesh, lm, iatom, 2))
                    write(109,*) imesh, lm, iatom, 3
                    write(109,*) real(vEff1MT(imesh, lm, iatom, 3)), aimag(vEff1MT(imesh, lm, iatom, 3))
                  end if
                end if

              end do
            end do
          end do
        end do
      end do
    end if


  end subroutine GenVeff1

  !--------------------------------------------------------------------------------------------------------------------------------------
  !> @author
  !> Christian-Roman Gerhorst, Forschungszentrum Jülich: IAS1 / PGI1
  !>
  !> @brief
  !> This method determines the multipole moments needed for constructing the first variation of the Coulomb potential in the
  !> interstitial region due to a phonon perturbation.
  !>
  !>  @param[in]  atoms        : Atoms type, see types.f90
  !>  @param[in]  cell         : Unit cell type, see types.f90
  !>  @param[in]  ngpqdp       : Number of G-vectors for shifted G-set for q-point with index iqpt
  !>  @param[in]  ngdp         : Number of G-vectors for potentials and densities
  !>  @param[in]  harSw        : Switch to enable Hartree contribution
  !>  @param[in]  extSW        : Switch to enable external contribution
  !>  @param[in]  iDatom       : Index of displaced atom
  !>  @param[in]  iDtype       : Index of atom type to which displaced atom belongs
  !>  @param[in]  gdp          : G-vectors of potentials and densities
  !>  @param[in]  gpqdp        : G-vectors for shifted G-set for q-point with index iqpt
  !>  @param[in]  rho0IRpw     : Plane-wave coefficients of the unperturbed and converged interstitial density parsed from Fleur
  !>  @param[in]  rho0MTsh     : Spherical harmonic coefficients of unperturbed and converged muffin-tin density parsed from Fleur
  !>  @param[in]  qpoint       : Qpoint in internal coordinates as in qpts%bk(1:3, iqpt)
  !>  @param[in]  rho1MT       : Spherical harmonic muffin-tin expansion coefficients of the linear density variation without the
  !>                             muffin-tin gradient of the density
  !>  @param[in]  rho1PW       : Plane-wave interstitial expansion coefficients of the linear density variation
  !>  @param[in]  grRho0MT     : Spherical harmonic coefficients of the gradient of the muffin-tin unperturbed density
  !>  @param[out] qlmPhon1Coul : Multipole moments for the Hartree, external or Coulomb potential depending on which switches
  !>
  !--------------------------------------------------------------------------------------------------------------------------------------
  subroutine mpmomVc1( atoms, cell, qpoint, gpqdp, harSw, extSw, ngpqdp, ngdp, gdp, rho0IRpw, rho0MTsh, grRho0MT, rho1PW, rho1MT, qlmPhon1Coul, iDatom, iDtype )

    use m_types, only : t_atoms, t_cell
    use m_JPConstants, only : fpi, c_im
    use m_jpPotDensHelper, only : CalcQlmHarSurf

    implicit none

    ! Type Parameters
    type(t_atoms),              intent(in)  :: atoms
    type(t_cell),               intent(in)  :: cell

    ! Scalar Paramters
    integer,                    intent(in)  :: ngpqdp
    integer,                    intent(in)  :: ngdp
    logical,                    intent(in)  :: harSw
    logical,                    intent(in)  :: extSW
    integer,                    intent(in)  :: iDatom
    integer,                    intent(in)  :: iDtype

    ! Array Parameters
    integer,                    intent(in)  :: gdp(:, :)
    integer,                    intent(in)  :: gpqdp(:, :)
    complex,                    intent(in)  :: rho0IRpw(:, :)
    complex,                    intent(in)  :: rho0MTsh(:, :, :, :)
    real,                       intent(in)  :: qpoint(:)
    complex,                    intent(in)  :: rho1MT(:, :, :, :)
    complex,                    intent(in)  :: rho1PW(:, :)
    complex,                    intent(in)  :: grRho0MT(:, :, :, :)
    complex,       allocatable, intent(out) :: qlmPhon1Coul(:,:,:)

    !------------------------------------------------------------------------------------------------------------------------------------
    ! Local Scalar Variabels
    ! itype        : runs over all types of atoms
    ! mqn_m        : magnetic quantum number m
    ! idirec       : runs over all directions of displacement
    !------------------------------------------------------------------------------------------------------------------------------------
    integer                                 :: itype
    integer                                 :: iatom
    integer                                 :: ieqat
    integer                                 :: mqn_m
    integer                                 :: idirec
    integer                                 :: oqn_l
    integer                                 :: lm

    !------------------------------------------------------------------------------------------------------------------------------------
    ! Local Array Variabels
    ! qlmo1Hartr   : MT charge contribution for multipole moments of first variation of the Hartree potential
    ! qlmp1Hartr   : MT plane wave contribution for multipole moments of first variation of the Hartree potential
    ! qlmGrad0Ext  : multipole moments needed to construct the gradient of the unperturbed external potential
    !------------------------------------------------------------------------------------------------------------------------------------
    complex,       allocatable              :: qlmo1Hartr(:,:, :)
    complex,       allocatable              :: qlmp1Hartr(:,:, :)
    complex,       allocatable              :: qlmGrad0Ext( :, :, :)
    complex,       allocatable              :: qlmHartSurf(:, :)


    allocate( qlmPhon1Coul((atoms%lmaxd + 1)**2, atoms%nat, 3) )
    qlmPhon1Coul = cmplx(0.0, 0.0)
    if ( harSW ) then
      allocate( qlmo1Hartr( (atoms%lmaxd + 1)**2, atoms%nat, 3) )
      allocate( qlmp1Hartr( (atoms%lmaxd + 1)**2, atoms%nat, 3) )
      qlmo1Hartr = cmplx(0.0, 0.0)
      qlmp1Hartr = cmplx(0.0, 0.0)
      ! calculate the multipole moments (equation 7.57, PhD thesis Aaron Klüppelberg) for IR part of first variation of Hartree potential
      call mpMom1Hart(atoms, iDatom, iDtype, rho1MT, grRho0MT, qlmo1Hartr)
      call pwMom1Hart(atoms, cell, gpqdp, ngpqdp, qpoint, rho1PW, qlmp1Hartr)

      do idirec = 1, 3
        iatom = 0
        do itype = 1, atoms%ntype
          do ieqat = 1, atoms%neq(itype)
            iatom = iatom + 1
            do oqn_l = 0, atoms%lmax(itype)
              do mqn_m = -oqn_l, oqn_l
                lm = oqn_l * (oqn_l + 1) + 1 + mqn_m
                qlmPhon1Coul(lm, iatom, idirec) = qlmPhon1Coul(lm, iatom, idirec) + qlmo1Hartr(lm, iatom, idirec)
                qlmPhon1Coul(lm, iatom, idirec) = qlmPhon1Coul(lm, iatom, idirec) - qlmp1Hartr(lm, iatom, idirec)
              end do ! mqn_m
            end do ! oqn_l
          end do ! ieqat
        end do !itype
      end do ! idirec

      ! Consistent with Aaron (7.57) and (2.34) in CRG 16.12.2020-19:00

      ! Additional multipole moment dependent on discontinuity of density at the MT boundary that is shared for gradient of
      ! unperturbed Hartree potential and the linear variation of the Hartree potential using formulas 7.54, 3.30 and 7.56 from
      ! PhDthesAK
      allocate( qlmHartSurf( (atoms%lmaxd + 1 )**2, 3) )
      qlmHartSurf(:, :) = cmplx(0., 0.)
      call CalcQlmHarSurf( atoms, cell, iDtype, iDatom, ngdp, gdp, rho0IRpw, rho0MTsh, qlmHartSurf )
      do idirec = 1, 3
        do oqn_l = 0, atoms%lmax(iDtype)
          do mqn_m = -oqn_l, oqn_l
            lm = oqn_l * (oqn_l + 1) + 1 + mqn_m
            qlmPhon1Coul(lm, iDatom, idirec) = qlmPhon1Coul(lm, iDatom, idirec) + qlmHartSurf(lm, idirec)
          end do ! mqn_m
        end do ! oqn_l
      end do ! idirec
    end if

    ! N_check: qlmPhon1Coul(lm, iDatom, idir)+=qlmHartSurf(lm, idir)[\alpha_{0}]
    !
    ! Consistent with Aaron (7.61) and (2.39) in CRG 16.12.2020-19:00

    ! Maintenence
    if (.false.) then
      iatom = 0
      do itype = 1, atoms%ntype
        do ieqat = 1, atoms%neq(itype)
          iatom = iatom + 1
          do idirec = 1, 3
            do oqn_l = 0, atoms%lmax(1)
              do mqn_m = -oqn_l, oqn_l
                lm = oqn_l * ( oqn_l + 1 ) + 1 + mqn_m
                write(2023, '(3(i8),2(f15.8))') lm, idirec, iatom, -qlmPhon1Coul(lm, iatom, idirec)
              end do ! mqn_m
            end do ! oqn_l
          end do ! idir
        end do ! ieqat
      end do ! itype
    end if

    ! calculate the multipole moments (equation 7.38, PhD thesis Aaron Klüppelberg) needed for the first variation of the external
    ! potential and the gradient of the unperturbed external potential. The calculation can be restricted to l = 1 (see PhD thesis Aaron
    ! Klüppelberg for details). We need to evaluate for different types only due to the gradient of the unperturbed external potential.
    ! Considering the first variation only the value for the displaced muffin tin type is required. The gradient of r^l Y*_lm(\hat{r})
    ! is calculated using (7.49, Phd thesis Aaron). (Markus Betzinger): In the end, we can start from (7.42) and can expand the unit
    ! vector in spherical harmonics; the expansions coefficients are the the pseuocharges.The c_(i, m) are defined in
    ! (4.28, PhD thesis Aaron Klüppelberg).
    allocate(qlmGrad0Ext(-1 : 1, atoms%nat, 3))
    qlmGrad0Ext = cmplx(0.0, 0.0)
    if ( extSw ) then
      do idirec = 1, 3
        iatom = 0
        do itype = 1, atoms%ntype
          do ieqat = 1, atoms%neq(itype)
            iatom = iatom + 1
            if (iatom == iDatom) then
              do mqn_m = -1, 1
                qlmGrad0Ext(mqn_m, iatom, idirec) = -3. / fpi * atoms%zatom(iDtype) * c_im(idirec, mqn_m + 2)
              end do ! mqn_m
            end if
          end do ! ieqat
        end do ! itype
      end do ! idirec
    end if

    ! N_check: qlmGrad0Ext(1m, iatom, idir)=-\frac{3Z_{\alpha_{0}}}{4\pi}c_im(idir,m)
    !
    ! Consistent with Aaron (7.38) and (2.21) in CRG 16.12.2020-19:00

    ! Get multipole moments for interstitial part of first variation of the Coulomb potential by summing those for the first variation
    ! of the Hartree and the external potential for the displaced atom type (which is effectively a single atom if no symmetry at all).
    do idirec = 1, 3
      iatom = 0
      do itype = 1, atoms%ntype
        do ieqat = 1, atoms%neq(itype)
          iatom = iatom + 1
          oqn_l = 1
          do mqn_m = -oqn_l, oqn_l
            lm = oqn_l * ( oqn_l + 1 ) + 1 + mqn_m
            qlmPhon1Coul(lm, iatom, idirec) = qlmPhon1Coul(lm, iatom, idirec) + qlmGrad0Ext(mqn_m, iatom, idirec)
          end do ! mqn_m
        end do !ieqat
      end do ! itype
    end do ! idirec

  end subroutine mpmomVc1

  !--------------------------------------------------------------------------------------------------------------------------------------
  !> @author
  !> Christian-Roman Gerhorst, Forschungszentrum Jülich: IAS1 / PGI1
  !>
  !> @brief
  !> Calculates the multipole moments of the muffin-tin charge which are needed to calculate the first variation of the Hartree potential.
  !>
  !> @details
  !> Calculates the multipole moments of the original muffin-tin charge for every l and m of every atom which can be displaced in three
  !> directions which are needed to calculate the first variation of the Hartree potential in the interstitial region. This routine
  !> evaluates the first term in (7.57, PhD thesis Aaron Klueppelberg) for the unit cell at R' = 0. If β = α or not, i.e. whether the
  !> current MT is displaced or not, is considered during setup of rho1MT.
  !>
  !> @param[in]   atoms      : Atoms type, see types.f90
  !> @param[in]   iDatom     : Index of displaced atom
  !> @param[in]   iDtype     : Index of atom type to which displaced atom belongs
  !> @param[in]   rho1MT     : Spherical harmonic muffin-tin expansion coefficients of the linear density variation without the
  !>                           muffin-tin gradient of the density
  !> @param[in]   grRho0MT   : Spherical harmonic coefficients of the gradient of the muffin-tin unperturbed density
  !> @param[out]  qlmo1Hartr : Multipole moments of the muffin-tin charge which are needed to calculate the first variation of the
  !>                           Hartree potential
  !--------------------------------------------------------------------------------------------------------------------------------------
  subroutine mpMom1Hart(atoms, iDatom, iDtype, rho1MT, grRho0MT, qlmo1Hartr)

    use m_intgr, only : intgr3LinIntp
    use m_types, only : t_atoms

    implicit none

    ! Type Parameter
    type(t_atoms),  intent(in)   :: atoms

    ! Scalar Parameter
    integer,        intent(in)   :: iDatom
    integer,        intent(in)   :: iDtype

    ! Array Parameter
    complex,        intent(in)   :: rho1MT( :, :, :, :)
    complex,        intent(in)   :: grRho0MT(:, :, :, :)
!    complex,        intent(out)  :: qlmo1Hartr(-atoms%lmaxd - 1:atoms%lmaxd + 1,0:atoms%lmaxd + 1,atoms%nat, 3)
    complex,        intent(out)  :: qlmo1Hartr(:, :, :)

    !------------------------------------------------------------------------------------------------------------------------------------
    ! Local Scalar Variables
    ! fintReal : contains real result of integral
    ! fintImag : contains imaginary result of integral
    ! iatom    : runs over all atoms
    ! itype    : runs over all types
    ! oqn_l    : orbital quantum number
    ! idir   : runs over all displacement directions
    ! mqn_m    : magnetic_quantum number
    ! imesh    : runs over the grid points of the MT mesh
    ! ieq      : runs over all equivalent atoms of an atom type
    ! lm       : encodes oqn_l and mqn_m
    !------------------------------------------------------------------------------------------------------------------------------------
    real                         :: fintReal
    real                         :: fintImag
    integer                      :: iatom
    integer                      :: itype
    integer                      :: oqn_l
    integer                      :: idir
    integer                      :: mqn_m
    integer                      :: imesh
    integer                      :: ieqat
    integer                      :: lm

    !------------------------------------------------------------------------------------------------------------------------------------
    ! Local Array Variables
    ! fReal : Real part of integrand
    ! fImag : Imaginary part of integrand
    !------------------------------------------------------------------------------------------------------------------------------------
    real, allocatable            :: fReal(:)
    real, allocatable            :: fImag(:)
    complex, allocatable         :: rho1MTfull(:, :, :, :)


    ! Init dummy assumed-shape array
    qlmo1Hartr(:, :, :) = 0.0

    !todo This can be moved to an outer scope!!!
    ! Within the MT contribution of the Coulomb potential the gradient of the density cancels in the Dirichelet boundary problem. That's
    ! why this contribution is not part of rho1MT during the convergence of rho1MT. For the integrand (fReal, fImag), however, the full
    ! density variation is required. It is better to add it once here than to subtract it as often as the amount of k-points used for the
    ! calculation. Note, that the gradient is only to be subtracted for the displaced atom
    allocate( rho1MTfull(atoms%jmtd, (atoms%lmaxd + 1)**2, atoms%nat, 3) )
    rho1MTfull(:, :, :, :) = cmplx(0.0, 0.0)
    rho1MTfull(1:atoms%jmtd, 1:(atoms%lmaxd + 1)**2, 1:atoms%nat, 1:3) = &
                                                                    & rho1MT(1:atoms%jmtd, 1:(atoms%lmaxd + 1)**2, 1:atoms%nat, 1:3)

    do idir = 1, 3
      do oqn_l = 0, atoms%lmax(iDtype)
        do mqn_m = -oqn_l, oqn_l
          lm = oqn_l * (oqn_l + 1) + mqn_m + 1
          do imesh = 1, atoms%jri(iDtype)
            rho1MTfull(imesh, lm, iDatom, idir) = rho1MTfull(imesh, lm, iDatom, idir) - grRho0MT(imesh, lm, iDatom, idir)
          end do
        end do
      end do
    end do

    ! First term of 7.57 is stored in qlmo1Hartr
    allocate( fReal(atoms%jmtd), fImag(atoms%jmtd) )
    fReal(:) = 0.0
    fImag(:) = 0.0
    do idir = 1, 3
      iatom = 1
      do itype = 1, atoms%ntype
        do ieqat =  1, atoms%neq(itype)
          do oqn_l = 0, atoms%lmax(itype)
            do mqn_m = -oqn_l, oqn_l
              lm = oqn_l * (oqn_l + 1) + 1 + mqn_m
              fReal(:) = 0.0
              fImag(:) = 0.0
              do imesh = 1, atoms%jri(itype)
                ! NOTE: We have to multiply the Jacobi determinant as rho1MT is not decorated with any factor by default
                fReal(imesh) = (atoms%rmsh(imesh, itype)**(oqn_l + 2)) * real(rho1MTfull(imesh, lm, iatom, idir))
                fImag(imesh) = (atoms%rmsh(imesh, itype)**(oqn_l + 2)) * aimag(rho1MTfull(imesh, lm, iatom, idir))
              enddo ! imesh
              call intgr3LinIntp(fReal, atoms%rmsh(:, itype), atoms%dx(itype), atoms%jri(itype), fintReal, 1)
              call intgr3LinIntp(fImag, atoms%rmsh(:, itype), atoms%dx(itype), atoms%jri(itype), fintImag, 1)
              qlmo1Hartr(lm, iatom, idir) = cmplx(fintReal, fintImag)
            enddo ! mqn_m
          enddo ! oqn_l
          iatom = iatom + 1
        enddo ! ieqat
      enddo ! itype
    enddo ! idir

    ! N_check: qlmoHartr(lm,iatom,idir)=\int_{0}^{R_{MT^{\alpha}}}dr r^{l+2} rho1MTfull(r, lm, iatom, idir)
    ! Consistent with Aaron (7.57) part 1 and (2.34) in CRG 16.12.2020-19:00

  end subroutine mpMom1Hart

  !--------------------------------------------------------------------------------------------------------------------------------------
  !> @author
  !> Christian-Roman Gerhorst, Forschungszentrum Jülich: IAS1 / PGI1
  !>
  !> @brief
  !> Calculates the multipole moments of the real interstitial charge which are needed for the first variation of the Coulomb
  !> potential.
  !>
  !> @details
  !> This method calculates the real charge multipole moments contribution stemming from the plane-waves within the MTs for every
  !> atom type, every l and m and all 3 directions the atoms can be displaced to.
  !> They are needed to calculate the interstitial part of the first variation of the Hartree potential. The formula for this routine
  !> is the second term of equation 7.57 in the thesis of Aaron Klüppelberg, where we choose the unit cell specified by R' = 0.
  !>
  !> @param[in]   atoms     : Atoms type, see types.f90
  !> @param[in]   cell      : Unit cell type, see types.f90
  !> @param[in]   gpqdp     : G-vectors for shifted G-set for q-point with index iqpt
  !> @param[in]   ngpqdp    : Number of G-vectors for shifted G-set for q-point with index iqpt
  !> @param[in]   qpoint    : Qpoint in internal coordinates as in qpts%bk(1:3, iqpt)
  !> @param[in]   rho1IR    : Plane-wave interstitial expansion coefficients of the linear density variation
  !> @param[out] qlmp1Hartr : Multipole moments of real interstitial charge.
  !>
  !> @todo
  !> Markus: rein mathematisch sieht diese Routine korrekt aus; je nach Zeitaufwand der Routine solle man die Loopstruktur noch einmal überdenken; es ist moeglicherweise ratsam den Loop über iG soweit nach innen zu ziehen wie möglich, da ngpd schnell mehrere tausend gross werden kann
  !--------------------------------------------------------------------------------------------------------------------------------------
  subroutine pwMom1Hart(atoms, cell, gpqdp, npqgdp, qpoint, rho1PW, qlmp1Hartr)

    use m_constants, only : sfp_const
    use m_types
    use m_sphbes
    use m_jpPotDensHelper, only : Phasy1nSym

    implicit none

    ! Type and Scalar Arguments, documentation in mpmomVc1
    type(t_atoms),  intent(in)  :: atoms
    type(t_cell),   intent(in)  :: cell
    integer,        intent(in)  :: npqgdp

    ! Array Arguemts, documentation in mpmomVc1
    integer,        intent(in)  :: gpqdp(:, :)
    real,           intent(in)  :: qpoint(:)
    complex,        intent(in)  :: rho1PW(:, :)
!    complex,        intent(out) :: qlmp1Hartr(-atoms%lmaxd - 1:atoms%lmaxd + 1,0:atoms%lmaxd + 1,atoms%nat, 3)
    complex,        intent(out) :: qlmp1Hartr(:,:, :)

    ! Local Scalar Variables
    integer                     :: idirec ! current direction
    integer                     :: itype ! current type
    real                        :: sk3r ! argument of bessel function which the whole equation is also divided by
    real                        :: rl3 !stores R^(l + 3)
    integer                     :: oqn_l ! orbital quantum number l
    integer                     :: ll1 ! auxillary variable to store proper l for calculating lm
    integer                     :: mqn_m ! magnetic quantum number
    integer                     :: lm ! combination of l and m
    integer                     :: iG ! runs over all G vectors
    integer                     :: iatom ! runs over all atoms
    integer                     :: ineq
    complex                     :: sk3i ! stores ρ(G) * j_(l + 1)(|G + q| R_MT) / |G + q| R_MT
    real                        :: normedG
    complex                     :: cil ! stores everything which does not come from pylm

    ! Local Array Variables
    complex, allocatable        :: pylm(: , :) !stores result of phasy1nSym
    real,    allocatable        :: sbes(:) ! bessel function
    real                        :: Gqext(3)

    qlmp1Hartr(:, :, :) = cmplx(0.0, 0.0)

    allocate( pylm((atoms%lmaxd + 1)**2 , atoms%nat) )
    allocate( sbes(0 : atoms%lmaxd + 1) )
    pylm(:, :) = cmplx(0., 0.)
    sbes(:) = 0.

    ! we don't need the nstr factor because we cover all Gs
    do iG = 1, npqgdp
      Gqext(1:3) = matmul(cell%bmat(1:3, 1:3), gpqdp(1:3, iG) + qpoint(1:3))
      normedG = norm2(Gqext(1:3))

      ! Should be left here, because q can also be chosen zero for Veff1
      if ( normedG < 1e-12 ) cycle

      ! stores 4 π i^l exp(i(G + q) τ) Y*_lm(G + q) into pylm for every atom, where G and q are vectors
      pylm(:, :) = cmplx(0., 0.)
      call Phasy1nSym(atoms, cell, gpqdp(1:3, iG), qpoint, pylm)


      ! calculates ρ(G) * j_(l + 1)(|G + q| R_MT) / |G + q| R_MT and stores it to sk3i
      do idirec = 1, 3
        iatom  = 1
        do itype = 1, atoms%ntype
          sk3r = normedG * atoms%rmt(itype)
          sbes(:) = 0.
          call sphbes(atoms%lmax(itype) + 1, sk3r , sbes)
          sk3i = rho1PW(iG, idirec) / sk3r
          do ineq = 1, atoms%neq(itype)
          ! todo one R can be cvanceld with sk3i
            rl3 = atoms%rmt(itype)**3
            do oqn_l = 0, atoms%lmax(itype)
              cil = sk3i * sbes(oqn_l + 1)  * rl3
              ll1 = oqn_l * (oqn_l + 1) + 1
              do mqn_m = -oqn_l, oqn_l
                lm = ll1 + mqn_m
                qlmp1Hartr(lm, iatom, idirec) = qlmp1Hartr(lm, iatom, idirec) + cil * pylm(lm, iatom)
              enddo ! mqn_m
              rl3 = rl3 * atoms%rmt(itype)
            enddo ! oqn_l
            iatom = iatom + 1
          enddo ! ineq
        enddo ! itype
      enddo !idirec
    enddo !iG

    ! N_check: qlmp1Hartr(lm, iatom, idir)=\sum_{\bm{G}\neq\bm{q}} \frac{j_{l+1}(|\bm{G}+\bm{q}|R_{\alpha})R_{\alpha}^{l+2}}{|\bm{G}+\bm{q}|}rho1PW(iG, idir)*
    !                                     4\pi i^{l} e^{i2\pi (\bm{G}_{int}+\bm{q}_{int})\cdot\bm{tau}_{\alpha,int}}Y_{lm}^{*}(\hat{\bm{G}+\bm{q}})
    ! Consistent with Aaron (7.57) part 2 and (2.34) in CRG 16.12.2020-19:00

  end subroutine pwMom1Hart

  !--------------------------------------------------------------------------------------------------------------------------------------
  !> @author
  !> Christian-Roman Gerhorst, Forschungszentrum Jülich: IAS1 / PGI1
  !>
  !> @brief
  !> Calculates the pseudo charge density from the multipole moments.
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
  !> @param[out] psq     : pseudo charges, result of the routine
  !>
  !> @note
  !> Calculates the pseudo charge by directly calculating 7.43! it has been quickly adapted from psqpwVec so a lot of arguments are not required
  !> It is out of order at the moment and has to be reviewed.
  !>
  !--------------------------------------------------------------------------------------------------------------------------------------
  subroutine psqpwVecVar(atomsT, cellT, dimensT, gdp, qptn, ngdp, rhoPW, qlm, psq, iDatom)

#include "recycledRoutines/cpp_double.h"
     use m_sphbes
     use m_types
     use m_JPConstants, only : iu, tpi, fpi

     implicit none

     ! Scalar Parameters
     type(t_atoms),     intent(in)    :: atomsT
     type(t_cell),      intent(in)    :: cellT
     type(t_dimension), intent(in)    :: dimensT
     integer,           intent(in)    :: ngdp
     integer,           intent(in)    :: iDatom

     ! Array Parameters
     complex,           intent(in)    :: rhoPW(:, :)
     integer,           intent(in)    :: gdp(:, :)
     complex,           intent(in)    :: qlm(:,:, :)
     real,              intent(in)    :: qptn(:)
     complex,           intent(out)   :: psq(3, ngdp)

     !-----------------------------------------------------------------------------------------------------------------------------------
     ! Local Scalar Variables
     ! itype  : runs over all types
     ! rmtl   : contains R^l for the prefactor
     ! oqn_l  : orbital quantum number l
     ! p      : auxiliary variable to calculate prefactor
     ! nc     : loop variable to calculate prefactor
     ! fpo    : auxiliary variable which contains 1 / Ω
     ! iGvec  : runs over G-vectors which are used for the potential and the density
     ! idirec : runs over directions the atoms can be displaced to
     ! iatom  : runs over all atoms
     ! ncvn   : contains n + l for a certain atom type, while n is taken from a table in the Weinert paper
     ! ieq   : runs over all equal atoms
     ! n1     : auxiliary variable which contains n - l + 1
     ! ll1    : auxiliary variable to calculate lm
     ! mqn_m  : magnetic quantum number
     ! lm     : variable encoding oqn_l and mqn_m
     !-----------------------------------------------------------------------------------------------------------------------------------
     integer                          :: itype
     real                             :: rmtl
     integer                          :: oqn_l
     real                             :: p
     integer                          :: nc
     real                             :: fpo
     integer                          :: iGvec
     integer                          :: idirec
     integer                          :: iatom
     integer                          :: ncvn
     integer                          :: ieq
     integer                          :: n1
     integer                          :: ll1
     integer                          :: mqn_m
     integer                          :: lm

     !-----------------------------------------------------------------------------------------------------------------------------------
     ! Local Array Variables
     ! pn     : variable which stores prefactor with double faculties times 1 / R^l
     ! pylm   : contains 4 π i^l exp((Gvec + qpoint) * taual) Y*_lm(Gvec + qpoint)
     ! Gext   : G-vector in external coordinates
     ! sbes   : contains spherical bessel function
     ! sl     : auxiliary variable to calculate psq
     ! sm     : auxiliary variable to calculate psq
     ! sa     : auxiliary variable to calculate psq
     !-----------------------------------------------------------------------------------------------------------------------------------
     real                             :: pn(0:atomsT%lmaxd,atomsT%ntypd)
     complex                          :: pylm((atomsT%lmaxd + 1)**2, atomsT%nat)
     real                             :: Gqext(3)
     real                             :: sbes(0:atomsT%lmaxd + dimensT%ncvd + 1)
     complex                          :: sl(3)
     complex                          :: sm(3)
     complex                          :: sa(3)

     real :: x
     complex :: sf


     psq = cmplx(0.0, 0.0)
     ! This block calculates pn(l,itype) = (2l + 2nc(itype) + 3)!! / ((2l + 1)!! R^l) from equation 28 in the Weinert paper cited above or the
     ! second term of equation 7.63 (e.g.) in the PhD thesis of Aaron Klüppelberg. However, the formula is simplified in the following
     ! way: Assume l + n = l', so the factorials cancel each other until l' = l since 3 is 1 + 2 and gives the next odd number.
     ! That's why the loop starts at nc=l. Delivering only odd numbers is ensured by the term
     ! itself which gives all odd numbers until 2ncv(n) +  3 !).
     ! atoms%ncv(itype)=n+l is constant per atom type and thus not dependent on l, the value of ncv per type is taken from Table 1 in the
     ! Weinert paper.
     ! Current: calculate (2N + 5) !!
     !write(*, *) 'ncv', atomsT%ncv(1)
     pn = 0
     do itype = 1, atomsT%ntype
   !    rmtl = 1.0
       !do oqn_l = 0, atomsT%lmax(itype)
        oqn_l = 1
   !     if (oqn_l < atomsT%ncv(itype)) then
          p = 1.
          do nc = 1, (2 * atomsT%ncv(itype) + 3), 2 !ncv = N + l!!!!!!!
   !       do nc = oqn_l, atomsT%ncv(itype)
           ! p = p * (2 * nc + 3)
            p = p * nc
          end do ! nc
         ! pn(oqn_l, itype) = p / rmtl
          pn(1, itype) = p!!/ rmtl
   !     end if
      !   rmtl = rmtl * atomsT%rmt(itype)
       !end do ! oqn_l
     end do ! itype

     !write(1022, *) pn(1, 1)

     ! This block calculates the G /= 0 term in equation 28 of the Weinert paper cited above:
     ! \tilde \rho_s (K) = 4 π / Ω \sum_{lmi} (-i)^l (2l + 2nc(n) + 3)!! / ((2l + 1)!! R^l) j_{l + n + 1}(|G + q|R_i) / (|G + q| R_i)^{n-l+1}
     ! qlm \exp{-i (G + q) τ} Y_{lm}(G + q), which is the second term in equation 7.63. It is similiar to equation 28 in the Weinert
     ! paper cited above.
     ! The G = 0 term is not needed for calculating the interstitial contributions of the potentials and should be zero anyway.
     fpo = 1. / cellT%omtil

     ! todo loop over atoms correct?, in FLEUR there is only loop over types
     do iGvec = 1, ngdp
    !  Gqext = matmul(cellT%bmat, Gvec + qptn)


      ! calculates first exp(i (G + q) tau)  and multiplies recent factors before storing the final result to pylm
     !  call phasy1nSym(atomsT, cellT, gdp(:, iGvec), qptn, pylm)
       Gqext = matmul(cellT%bmat, gdp(:, iGvec) + qptn)
       if ( norm2(Gqext) == 0 ) cycle
       sa = 0.
       do idirec = 1, 3
         iatom = 0
         do itype = 1, atomsT%ntype
           ncvn = atomsT%ncv(itype)
           if (1 >= ncvn) cycle ! in earlier versions of FLEUR pn was set zero in this case, so zero was added to sl
           call sphbes(ncvn + 1, norm2(Gqext) * atomsT%rmt(itype), sbes)
           do ieq = 1, atomsT%neq(itype)
             iatom = iatom + 1
             if (iatom == iDatom) then
               x = tpi * dot_product(gdp(:, iGvec) + qptn, atomsT%taual(:, iatom))
               sf = exp(-iu *  x)
               sl(idirec) = 0.
               n1 = ncvn + 1
               !sl(idirec) = iu * atomsT%zatom(itype) * pn(1, itype) * sbes(ncvn + 2) / ((norm2(Gqext) * atomsT%rmt(itype))**n1) * sf * Gqext(idirec)
               sl(idirec) =  pn(1, itype) * iu * atomsT%zatom(itype) * sf * sbes(ncvn + 1) * Gqext(idirec) / ((norm2(Gqext) * atomsT%rmt(itype))**n1)
          !     write(1022, *) real(sl(1))
          !     NOstopNO
               sa(idirec) = sl(idirec)
             end if
           enddo ! ieq
         enddo ! itype
         psq(idirec, iGvec) = fpo * sa(idirec)
       enddo ! idirec
     enddo ! iGvec
   end subroutine psqpwVecVar

  !--------------------------------------------------------------------------------------------------------------------------------------
  !> @author
  !> Christian-Roman Gerhorst, Forschungszentrum Jülich: IAS1 / PGI1
  !>
  !> @brief
  !> Calculates the pseudo charge density from the multipole moments.
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
  !> @param[in]  atoms  : Atoms type, see types.f90
  !> @param[in]  cell   : Unit cell type, see types.f90
  !> @param[in]  dimens : Dimension type, see types.f90
  !> @param[in]  gpqdp  : G-vectors for shifted G-set for q-point with index iqpt
  !> @param[in]  qptn   : Qpoint in internal coordinates as in qpts%bk(1:3, iqpt)
  !> @param[in]  ngpqdp : Number of G-vectors for shifted G-set for q-point with index iqpt
  !> @param[in]  rho1IR : Plane-wave interstitial expansion coefficients of the linear density variation
  !> @param[in]  qlm    : Multipole moments for the Hartree, external or Coulomb potential depending on which switches
  !> @param[out] psq    : Pseudo charges to construct the interstitial Hartree, external or Coulomb potential depending on the
  !>                      the switches.
  !>
  !> @todo
  !> unite with same routine psqpwVeclp1 from grVeff0
  !--------------------------------------------------------------------------------------------------------------------------------------
  subroutine psqpwVec(atoms, cell, dimens, gpqdp, qptn, ngpqdp, rho1PW, qlm, psq)

#include "recycledRoutines/cpp_double.h"
     use m_sphbes
     use m_types
     use m_jpPotDensHelper, only : phasy1nSym

     implicit none

     ! Type Parameters
     type(t_atoms),                  intent(in)  :: atoms
     type(t_cell),                   intent(in)  :: cell
     type(t_dimension),              intent(in)  :: dimens

     ! Scalar Parameters
     integer,                        intent(in)  :: ngpqdp

     ! Array Parameters
     complex,                        intent(in)  :: rho1PW(:, :)
     integer,                        intent(in)  :: gpqdp(:, :)
     complex,                        intent(in)  :: qlm(:,:, :)
     real,                           intent(in)  :: qptn(:)
     complex,                        intent(out) :: psq(:, :)

     !-----------------------------------------------------------------------------------------------------------------------------------
     ! Local Scalar Variables
     ! itype  : runs over all types
     ! rmtl   : contains R^l for the prefactor
     ! oqn_l  : orbital quantum number l
     ! p      : auxiliary variable to calculate prefactor
     ! nc     : loop variable to calculate prefactor
     ! fpo    : auxiliary variable which contains 1 / Ω
     ! iGvec  : runs over G-vectors which are used for the potential and the density
     ! idirec : runs over directions the atoms can be displaced to
     ! iatom  : runs over all atoms
     ! ncvn   : contains n + l for a certain atom type, while n is taken from a table in the Weinert paper
     ! ieq   : runs over all equal atoms
     ! n1     : auxiliary variable which contains n - l + 1
     ! ll1    : auxiliary variable to calculate lm
     ! mqn_m  : magnetic quantum number
     ! lm     : variable encoding oqn_l and mqn_m
     !-----------------------------------------------------------------------------------------------------------------------------------
     integer                                     :: itype
     real                                        :: rmtl
     integer                                     :: oqn_l
     real                                        :: p
     integer                                     :: nc
     real                                        :: fpo
     integer                                     :: iGvec
     integer                                     :: idirec
     integer                                     :: iatom
     integer                                     :: ncvn
     integer                                     :: ieq
     integer                                     :: n1
     integer                                     :: ll1
     integer                                     :: mqn_m
     integer                                     :: lm
     complex                                     :: sl
     complex                                     :: sm
     complex                                     :: sa

     !-----------------------------------------------------------------------------------------------------------------------------------
     ! Local Array Variables
     ! pn     : variable which stores prefactor with double faculties times 1 / R^l
     ! pylm   : contains 4 π i^l exp((Gvec + qpoint) * taual) Y*_lm(Gvec + qpoint)
     ! Gext   : G-vector in external coordinates
     ! sbes   : contains spherical bessel function
     ! sl     : auxiliary variable to calculate psq
     ! sm     : auxiliary variable to calculate psq
     ! sa     : auxiliary variable to calculate psq
     !-----------------------------------------------------------------------------------------------------------------------------------
     real,              allocatable              :: pn(:, :)
     complex,           allocatable              :: pylm(:, :)
     real,              allocatable              :: sbes(:)
     real                                        :: Gqext(3)


     ! Init assumed shape array
     psq(:, :) = cmplx(0.0, 0.0)

     ! This block calculates pn(l,itype) = (2l + 2nc(itype) + 3)!! / ((2l + 1)!! R^l) from equation 28 in the Weinert paper cited above or the
     ! second term of equation 7.63 (e.g.) in the PhD thesis of Aaron Klüppelberg. ncv(itype) is according to Weinert's paper still
     ! constant and defined as n + l and the formula to be calculated is simplified in the following
     ! way: Assume l + n = l', so the factorials cancel each other until l' = l - 1 since 3 is 1 + 2 and gives the next odd number.
     ! Since 2l + 2nc(itype) + 3 is odd we then multiply over all odd numbers from (2 l + 3) until (2ncv(itype) +3) which can be done
     ! by incrementing l until ncv(itype). That's why the loop starts at nc=l.
     ! The value of ncv per type is taken from Table 1 in the Weinert paper and given by the Fleur input .
     allocate(pn(0:atoms%lmaxd + 1,atoms%ntype))
     pn(:, :) = 0.
     do itype = 1, atoms%ntype
       rmtl = 1.0
       do oqn_l = 0, atoms%lmax(itype)
         if (oqn_l < atoms%ncv(itype)) then
           p = 1.
           do nc = oqn_l, atoms%ncv(itype)
             p = p * (2 * nc + 3)
           end do ! nc
           pn(oqn_l, itype) = p / rmtl
         end if
         rmtl = rmtl * atoms%rmt(itype)
       end do ! oqn_l
     end do ! itype


     ! This block calculates the G /= 0 term in equation 28 of the Weinert paper cited above:
     ! \tilde \rho_s (K) = 4 π / Ω \sum_{lmi} (-i)^l (2l + 2nc(n) + 3)!! / ((2l + 1)!! R^l) j_{l + n + 1}(|G + q|R_i) / (|G + q| R_i)^{n-l+1}
     ! qlm \exp{-i (G + q) τ} Y_{lm}(G + q), which is the second term in equation 7.63. It is similiar to equation 28 in the Weinert
     ! paper cited above.
     ! The G = 0 term is not needed for calculating the interstitial contributions of the potentials and should be zero anyway.
     allocate(pylm((atoms%lmaxd + 1)**2, atoms%nat))
     allocate(sbes(0:dimens%ncvd + 1 ))
     pylm(:, :) = cmplx(0., 0.)
     sbes(:) = 0.

     fpo = 1. / cell%omtil

     do iGvec = 1, ngpqdp
       Gqext = matmul( cell%bmat(1:3, 1:3), gpqdp(1:3, iGvec) + qptn(1:3) )
       if ( norm2(Gqext(1:3)) < 1e-11 ) cycle
       pylm(:, :) = cmplx(0., 0.)
       call Phasy1nSym(atoms, cell, gpqdp(:, iGvec), qptn, pylm)
       do idirec = 1, 3
         sa = cmplx(0., 0.)
         iatom = 1
         do itype = 1, atoms%ntype
           ncvn = atoms%ncv(itype)
           sbes(:) = 0.
           call sphbes(ncvn + 1, norm2(Gqext) * atoms%rmt(itype), sbes)
           do ieq = 1, atoms%neq(itype)
             sl = cmplx(0., 0.)
             do oqn_l = 0, atoms%lmax(itype)
               if (oqn_l >= ncvn) cycle ! in earlier versions of FLEUR pn was set zero in this case, so zero was added to sl
               n1 = ncvn - oqn_l + 1
               ll1 = oqn_l * (oqn_l + 1) + 1
               sm = cmplx(0., 0.)
               do mqn_m = -oqn_l, oqn_l
                  lm = ll1 + mqn_m
                  sm = sm + qlm(lm, iatom, idirec) * conjg(pylm(lm, iatom))
               enddo ! mqn_m
               ! if only external potential the fact that we have only contribution is accounted for by the fact that sm has only contributions for l = 1
               sl = sl + sbes(ncvn + 1) * sm * (pn(oqn_l, itype) / ((norm2(Gqext) * atoms%rmt(itype))**n1))
             enddo ! oqn_l
             sa = sa + sl
             iatom  = iatom + 1
           enddo ! ieq
         enddo ! itype
         ! In the first iteration rhoPW is set to zero so that we have now rhoPW contribution to the external pseudo density. For effectiveness we avoid an if clause this way
         psq(idirec, iGvec) = rho1PW(iGvec, idirec) + fpo * sa
       enddo ! idirec
     enddo ! iGvec

     ! N_check: psq(iG, idir)=[rho1PW(iG, idir) +] frac{1}{\Omega}\sum_{\alpha l m} qlm(lm, iatom, idir)*
     !                             4\pi (-i)^{l} e^{-i2\pi (\bm{G}_{int}+\bm{q}_{int})\cdot\bm{tau}_{\alpha,int}}Y_{lm}(\hat{\bm{G}+\bm{q}})*
     !                             \frac{j_{l+N+1}(|\bm{G}+\bm{q}|R_{\alpha})}{(|\bm{G}+\bm{q}|R_{\alpha})^{N+1}}*
     !                             \frac{(2l + 2N + 3)!!}{(2l + 1)!! R_{\alpha}^{l})}
     !
     ! Consistent with Aaron (7.63) and (2.38) in CRG 16.12.2020-19:00
   end subroutine psqpwVec


  !--------------------------------------------------------------------------------------------------------------------------------------
  !> @author
  !> Christian-Roman Gerhorst, Forschungszentrum Jülich: IAS1 / PGI1
  !>
  !> @brief
  !> Solves Dirichelet boundary value problem without symmetry
  !>
  !> @details
  !> This routine solves the Dirichlet boundary value problem. It is based on the paper
  !> M.Weinert, J.Math.Phys. 22(11) (1981) p.2434 using equations 7-9. The resulting formula
  !> can be extracted from equation 7.66 or 7.67 in the Phd thesis of Aaron Klueppelberg which
  !> is finally implemented in this routine. The difference to the existing routine vmts in FLEUR
  !> is that here no stars and lattice harmonics are exploited.
  !>
  !>
  !> @param[in]   atoms    : Atoms type, see types.f90
  !> @param[in]   cell     : Unit cell type, see types.f90
  !> @param[in]   ngpqdp   : Number of G-vectors for shifted G-set for q-point with index iqpt
  !> @param[in]   idirec   : Displacement direction
  !> @param[in]   iatom    : Index of atom in which the boundary value problem is solved
  !> @param[in]   itype    : Index of atom type to which atom belongs in which the boundary value problem is solved
  !> @param[in]   iDatom   : Index of displaced atom
  !> @param[in]   VextFull : Switch to activate volume term contribution of muffin-tin boundary problem
  !> @param[in]   gpqdp    : G-vectors for shifted G-set for q-point with index iqpt
  !> @param[in]   qpoint   : Qpoint in internal coordinates as in qpts%bk(1:3, iqpt)
  !> @param[in]   rho1IR   : Plane-wave interstitial expansion coefficients of the linear density variation
  !> @param[in]   rho1MT   : Spherical harmonic muffin-tin expansion coefficients of the linear density variation without the
  !>                         muffin-tin gradient of the density
  !> @param[in]   harSw    : Switch to enable Hartree contribution
  !> @param[in]   extSw    : Switch to enable external contribution
  !> @param[in]   vHarNum  : Switch to enable special terms that are needed to optimize the integrals rho1 Veff1 and rho1 Vext1
  !> @param[in]   linIntp  : Switch to enable linear interpolation at the integration of intgr2
  !> @param[out]  vr       : Result of the Dirichelet boundary value problem, depending on the switches it incorporates the Hartree
  !>                         external or Coulomb potential
  !>
  !> @note              : This routine cannot calculate the full first variation of the external potential, but always the sum of Hartree
  !> and external potential. Furthermore, the first term of 7.47 cancels away with the first term in 7.48 which is a reason of their sum
  !> in the Sternheimer equation but disables the user to calculate the full first variation of the external potential autonomously.
  !--------------------------------------------------------------------------------------------------------------------------------------
   subroutine vmtsNsym(atoms, cell,  ngpqdp, idirec, iatom, itype, iDatom, vExtFull, gpqdp, qpoint,  vpw, rho, vr, harSw, extSw, vHarNum, linIntp)

#include"recycledRoutines/cpp_double.h"
      use m_JPConstants
      use m_types
      use m_intgr, only : intgr2, intgr2LinIntp
      use m_sphbes
      use mod_juPhonUtils
      use m_jpPotDensHelper, only : phasy1nSym

      implicit none

      ! Type paramter
      type(t_atoms),    intent(in)  :: atoms
      type(t_cell),     intent(in)  :: cell

      ! Scalar parameter
      integer,          intent(in)  :: ngpqdp
      integer,          intent(in)  :: idirec
      integer,          intent(in)  :: iatom
      integer,          intent(in)  :: itype
      integer,          intent(in)  :: iDatom
      logical,          intent(in)  :: vExtFull
      logical,          intent(in)  :: harSw
      logical,          intent(in)  :: extSw
      logical,          intent(in)  :: vHarNum
      logical,          intent(in)  :: linIntp

      ! Array Parameter
      complex,          intent(in)  :: vpw(:, :)
      complex,          intent(in)  :: rho(:,:,:,:)
      integer,          intent(in)  :: gpqdp(:, :)
      real,             intent(in)  :: qpoint(:)
      complex,          intent(out) :: vr(:, :)

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
      ! idirec   : runs over all directions the atoms can be displaced to
      ! imesh    : runs over mesh points of local sphere mesh
      !----------------------------------------------------------------------------------------------------------------------------------
      integer                       :: oqn_l
      integer                       :: mqn_m
      integer                       :: l21
      real                          :: fpl21
      real                          :: rmtl
      real                          :: rmt2l
      integer                       :: lm
      integer                       :: iGvec
      integer                       :: imesh
      real                          :: rrlr
      real                          :: ror

      !----------------------------------------------------------------------------------------------------------------------------------
      ! Local Array Variables
      ! Gqext    : cartesian representation of G + q
      ! sbf      : contains spherical Bessel functions
      ! pylm     : contains 4 π i^l exp((Gvec + qpoint) * taual) Y*_lm(Gvec + qpoint)
      ! vtl      : contains the surface integral
      ! rrl      : auxiliary variable to calculate the volume integral
      ! rrl1     : auxiliary variable to calculate the volume integral
      ! f1r      : contains the integrand of the first type of integral occuring for the volume integral part
      ! f2r      : contains the integrand of the second type of integral occuring for the volume integral part
      ! x1r      : contains the result of the the first type of integral occuring for the volume integral part
      ! x2r      : contains the result of the the first type of integral occuring for the volume integral part
      !todo comment imag
      !----------------------------------------------------------------------------------------------------------------------------------
      real      ,allocatable                    :: sbf(:)
      complex   ,allocatable                    :: pylm(:, :)
      complex   ,allocatable                    :: vtl(:)
      real      ,allocatable                    :: rrl(:)
      real      ,allocatable                    :: rrl1(:)
      complex   ,allocatable                    :: f1r(:)
      complex   ,allocatable                    :: f2r(:)
      real      ,allocatable                    :: f1rReal(:)
      real      ,allocatable                    :: f2rReal(:)
      real      ,allocatable                    :: f1rImag(:)
      real      ,allocatable                    :: f2rImag(:)
      real      ,allocatable                    :: x1rReal(:)
      real      ,allocatable                    :: x2rReal(:)
      real      ,allocatable                    :: x1rImag(:)
      real      ,allocatable                    :: x2rImag(:)
      real                                      :: Gqext(3)


      allocate(vtl((atoms%lmaxd + 1)**2))
      allocate(sbf(0:atoms%lmaxd))
      allocate(pylm((atoms%lmaxd + 1)**2, atoms%nat))

      ! For a given iatom and itype the second term in equation (7.66 / 7.67, PhD thesis Aaron Klüppelberg) is evaluated, beginning from the sum
      ! todo optimize this routine by only taking this contribution into account for the first iteration. NOTE: Take also the (r / R)^l
      ! contribution into account!
      vtl(:) = cmplx(0., 0.)

      do iGvec = 1, ngpqdp
         Gqext(1:3) = matmul(cell%bmat(1:3, 1:3), gpqdp(1:3, iGvec) + qpoint(1:3))
         if (norm2(Gqext(1:3)) < 1e-9 ) cycle
         pylm(:, :) = cmplx(0., 0.)
         call phasy1nSym(atoms, cell, gpqdp(1:3, iGvec), qpoint, pylm)
         sbf(:) = 0.
         call sphbes(atoms%lmax(itype), norm2(Gqext(1:3)) * atoms%rmt(itype), sbf)
         do oqn_l = 0, atoms%lmax(itype)
            do mqn_m = -oqn_l, oqn_l
               lm = oqn_l * (oqn_l + 1) + mqn_m + 1
               ! todo we have to decide whether idirec should stay outside or not, then idirec dimension of vpw can also stay there
               vtl(lm) = vtl(lm) + vpw(iGvec, idirec) * sbf(oqn_l) * pylm(lm, iatom)
            end do ! iGvec
         end do !oqn_l
      end do ! iGvec
      deallocate(pylm)
      deallocate(sbf)

      allocate(rrl(atoms%jmtd))
      rrl = 0.
      ! For a given iatom and itype the first term in equation (7.66 / 7.67, PhD thesis Aaron Klüppelberg) is calculated first and then
      ! the second surface integral contribution inclucing the correct prefactor is added so that the complete Dirichelet boundary value
      ! problem is solved.
      ! Concerning the first term, basically, the integral is splitted into four integrals with respect to the definitions of r_> and r_<.
      ! This gives two types of integrands x1r = r^(l + 2) * ρ(r) dr and x2r = r^(1 - l) * ρ(r) dr. Considering the bounds of the
      ! integrals they can be rearranged to:
      ! 4π / (2 l + 1) * (1 / r^(l + 1) \int_0^r s^(l + 2) ρ ds - r^l / R^(2 * l + 1) \int_0^R ρ s^(l + 2)
      ! + r^l (\int_0^R 1 / s^(l + 1) ρ ds - \int_0^r 1 / s^(l + 1) ρ))
      vr(:, :) = cmplx(0., 0.)
      if (harSw) then
        allocate(rrl1(atoms%jmtd))
        allocate(f1r(atoms%jmtd))
        allocate(f2r(atoms%jmtd))
        allocate(f1rReal(atoms%jmtd))
        allocate(f2rReal(atoms%jmtd))
        allocate(f1rImag(atoms%jmtd))
        allocate(f2rImag(atoms%jmtd))
        allocate(x1rReal(atoms%jmtd), x1rImag(atoms%jmtd))
        allocate(x2rReal(atoms%jmtd), x2rImag(atoms%jmtd))
        rrl1 = 0.
        f1rReal = 0.
        f1rImag = 0.
        f1r = cmplx(0., 0.)
        f2rReal = 0.
        f2rImag = 0.
        f2r = cmplx(0., 0.)
        x1rReal = 0.
        x2rReal = 0.
        x1rImag = 0.
        x2rImag = 0.
        do oqn_l = 0, atoms%lmax(itype)
          l21 = 2 * oqn_l + 1
          fpl21 = fpi / l21 ! prefactor 1st term
          rmtl = 1. / atoms%rmt(itype)**oqn_l ! 1 / R^l
          rmt2l = 1. / atoms%rmt(itype)**l21 ! 1 / R^(2 l + 1)
          rrl(:) = 0.
          rrl1(:) = 0.
          do imesh = 1, atoms%jri(itype)
            rrl(imesh) = atoms%rmsh(imesh, itype)**oqn_l ! r^l
            rrl1(imesh) = 1. / (rrl(imesh) * atoms%rmsh(imesh, itype)) ! 1 / r^(l + 1)
          end do
          do mqn_m = -oqn_l, oqn_l
            lm = oqn_l * (oqn_l + 1) + mqn_m + 1
            do imesh = 1, atoms%jri(itype)
              ! In contrast to Fleur we have to multiply the Jacobi determinant in this place!
              x1rReal(imesh) = rrl(imesh) *  real(rho(imesh, lm, iatom, idirec)) * atoms%rmsh(imesh, itype) * atoms%rmsh(imesh, itype)
              x2rReal(imesh) = rrl1(imesh) * real(rho(imesh, lm, iatom, idirec)) * atoms%rmsh(imesh, itype) * atoms%rmsh(imesh, itype)
              x1rImag(imesh) = rrl(imesh) *  aimag(rho(imesh, lm, iatom, idirec))* atoms%rmsh(imesh, itype) * atoms%rmsh(imesh, itype)
              x2rImag(imesh) = rrl1(imesh) * aimag(rho(imesh, lm, iatom, idirec)) * atoms%rmsh(imesh, itype) * atoms%rmsh(imesh, itype)
             enddo
             if (linIntp) then
               call intgr2LinIntp(x1rReal, atoms%rmsh(1, itype), atoms%dx(itype), atoms%jri(itype), f1rReal)
               call intgr2LinIntp(x2rReal, atoms%rmsh(1, itype), atoms%dx(itype), atoms%jri(itype), f2rReal)
               call intgr2LinIntp(x1rImag, atoms%rmsh(1, itype), atoms%dx(itype), atoms%jri(itype), f1rImag)
               call intgr2LinIntp(x2rImag, atoms%rmsh(1, itype), atoms%dx(itype), atoms%jri(itype), f2rImag)
             else
               call intgr2(x1rReal, atoms%rmsh(1, itype), atoms%dx(itype), atoms%jri(itype), f1rReal)
               call intgr2(x2rReal, atoms%rmsh(1, itype), atoms%dx(itype), atoms%jri(itype), f2rReal)
               call intgr2(x1rImag, atoms%rmsh(1, itype), atoms%dx(itype), atoms%jri(itype), f1rImag)
               call intgr2(x2rImag, atoms%rmsh(1, itype), atoms%dx(itype), atoms%jri(itype), f2rImag)
             end if
             x1rReal(:) = 0.
             x2rReal(:) = 0.
             x1rImag(:) = 0.
             x2rImag(:) = 0.
             f1r(1:atoms%jri(itype)) = cmplx(f1rReal(1:atoms%jri(itype)), f1rImag(1:atoms%jri(itype)))
             f2r(1:atoms%jri(itype)) = cmplx(f2rReal(1:atoms%jri(itype)), f2rImag(1:atoms%jri(itype)))
             f1rReal(:) = 0.
             f2rReal(:) = 0.
             f1rImag(:) = 0.
             f2rImag(:) = 0.
             do imesh = 1, atoms%jri(itype)
               rrlr = rrl(imesh) * rmt2l ! r^l / R^(2 l +  1)
               ror = rrl(imesh) * rmtl ! r^l / R^l
               ! We neglect here the first term in the equation 7.48 and equation 7.47 which cancel away for the external potential because the muffin tin potentials are added in the Sternheimer equation and cannot considered to be independent. Everything before the last term only becomes relevant when rho is not set to zero as it is after the first iteration so where we have a Hartree contribution.
               ! The second term in 7.48 is relevant in all atoms but only has the boundary value of alpha.
               ! The first term cancels away within the Sternheimer equation so we don't have to care for the Kronecker delta.
               if (vHarNum) then
                 vr(imesh, lm) = ror * vtl(lm)
               else
                 vr(imesh, lm) = fpl21 * (rrl1(imesh) * f1r(imesh) - rrlr * f1r(atoms%jri(itype))      &
                                &+ rrl(imesh) * (f2r(atoms%jri(itype)) - f2r(imesh))) + ror * vtl(lm)
               end if
             end do ! imesh
             f1r(:) = cmplx(0., 0.)
             f2r(:) = cmplx(0., 0.)
          end do ! mqn_m
        end do ! oqn_l
      else if (.not.harSw .and. extSw) then
        do oqn_l = 0, atoms%lmax(itype)
          rmtl = 1. / atoms%rmt(itype)**oqn_l ! 1 / R^l
          rrl(:) = 0.
          do imesh = 1, atoms%jri(itype)
            rrl(imesh) = atoms%rmsh(imesh, itype)**oqn_l ! r^l
          end do
          do mqn_m = -oqn_l, oqn_l
            lm = oqn_l * ( oqn_l + 1 ) + 1 + mqn_m
            do imesh = 1, atoms%jri(itype)
              ror = rrl(imesh) * rmtl ! r^l / R^l
              vr(imesh, lm) = ror * vtl(lm)
            end do ! imesh
          end do ! mqn_m
        end do ! oqn_l
      end if

      ! N_check: vr(r, lm)=\frac{4\pi}{2l+1}\int_{0}^{R_{\alpha}}ds_{\alpha} s_{\alpha}^{2}\frac{r_{<}^{l}}{r_{>}^{l+1}}
      !                                   rho(r, lm, iatom)*[1-(\frac{r_{>}}{\bm{R}_{\alpha}})^{2l+1}]+boundary term
      !
      ! Consistent with Aaron (7.66) line 1 and (2.43) line 1 in CRG 16.12.2020-19:00

      !if ((extSw.and.(compPhon.or.anfix)).or.(vExtFull.and.extSw)) then
      if (extSw.and.vExtFull) then
        if ( iatom == iDatom ) then
          do lm = 2, 4
            do imesh = 1, atoms%jri(itype)
              vr(imesh, lm) = vr(imesh, lm) - atoms%zatom(itype) / atoms%rmsh(imesh, itype)**2 * ( 1 - (atoms%rmsh(imesh, itype) / atoms%rmt(itype))**3 ) * c_im(idirec, lm - 1)
            end do
          end do
        end if
      end if

      ! N_check: vr(r, lm)-=\frac{Z_{\alpha}4\pi}{3}r_{\alpha}^{-2}*
      !                               [1-(\frac{r_{\alpha}}{\bm{R}_{\alpha}})^{3}]\frac{3}{4\pi}c_im(idir, m)\delta_{l1}
      !
      ! Consistent with Aaron (7.47) and (2.30) in CRG 16.12.2020-19:00

   end subroutine vmtsNsym

end module m_jpVeff1
