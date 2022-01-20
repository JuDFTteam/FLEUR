!----------------------------------------------------------------------------------------------------------------------------------------
! Forschungszentrum Jülich, juPhon Plugin for the FLEUR program
!------------------------------------------------------------------------------
!
! MODULE: Test Routines for Potential Calculations
!
!> @author
!> Christian-Roman Gerhorst
!>
!> @brief
!> Routines within this module test the routines needed for the calculation of the first variation of the Coulomb potential and the
!> gradient of the unperturbed Coulomb potential.
!>
!> @note
!> Additional information and formulas pointing out the routines of this module can be found within this
!> <a href='jpTestPotential.pdf'>document</a>.
!----------------------------------------------------------------------------------------------------------------------------------------
module m_jpTestPotential

  use mod_juPhonUtils
  implicit none

  contains

  subroutine Test1stOrdPotentials( atoms, stars, lathar, sym, cell, dimens, input, V0Fleur, qpts, ngdp, paPoX, paPoY, paPoZ, harSw,      &
      & extSw, xcSw, logUnit, testCompareGrVeff0FleurSw, testVeff1Sw, testUnfoldStarsSw, testRadDerivativeSw, testGauntCoeffSw,    &
      & testGradLhExpandFuncSw, testContGrVeff0Sw, gdp, mlh_atom, nmem_atom, clnu_atom, rho0IR, rho0MT, memd_atom, noPtsCon, testVeff1IRqLatPeriod )

  use m_types
  use m_jpGrVeff0
  use m_jpVeff1
    use m_jpPotDensHelper
    use m_jPConstants, only : iu
    !use m_jpTestPotential, only : checkjuPhPots

  implicit none

  ! Type parameters
  type(t_atoms),     intent(in) :: atoms
  type(t_stars),     intent(in) :: stars
  type(t_sphhar),    intent(in) :: lathar
  type(t_sym),       intent(in) :: sym
  type(t_cell),      intent(in) :: cell
  type(t_dimension), intent(in) :: dimens
  type(t_input),     intent(in) :: input
  type(t_potential), intent(in) :: V0Fleur
  type(t_kpts),      intent(in) :: qpts

  ! Scalar parameters
  integer,           intent(in) :: ngdp
  real,              intent(in) :: paPoX
  real,              intent(in) :: paPoY
  real,              intent(in) :: paPoZ
  logical,           intent(in) :: harSw
  logical,           intent(in) :: extSw
  logical,           intent(in) :: xcSw
  integer,           intent(in) :: memd_atom
  integer,           intent(in) :: noPtsCon
  integer,           intent(in) :: logUnit
  logical,           intent(in) :: testCompareGrVeff0FleurSw
  logical,           intent(in) :: testVeff1Sw
  logical,           intent(in) :: testUnfoldStarsSw
  logical,           intent(in) :: testRadDerivativeSw
  logical,           intent(in) :: testGradLhExpandFuncSw
  logical,           intent(in) :: testContGrVeff0Sw
  logical,           intent(in) :: testGauntCoeffSw
  logical,           intent(in) :: testVeff1IRqLatPeriod

  ! Array parameters
  integer,           intent(in) :: gdp(:, :)
  integer,           intent(in) :: mlh_atom(:,0:,:)
  integer,           intent(in) :: nmem_atom(0:, :)
  complex,           intent(in) :: clnu_atom(:,0:,:)
  complex,           intent(in) :: rho0IR(:,:) !Interstitial density
  real,              intent(in) :: rho0MT(:,0:,:,:)  ! MT sphere density

  ! Scalar variables
  integer                       :: ispin
  integer                       :: itype
  integer                       :: ilh
  integer                       :: imesh
  logical                       :: testGoldstein
  ! Array variables
  complex, allocatable          :: grRho0MT(:, :, :, :)
    complex,           allocatable              :: rho0MTsh(:, :, :, :)
  complex, allocatable          :: grVeff0IR(:, :)
  complex, allocatable          :: grVeff0MT(:, :, :, :)
  real,    allocatable          :: myrho0MT(:,:,:,:)  ! MT sphere density
    real,              allocatable :: gausWts(:) ! gaussian weights belonging to gausPts
    complex,           allocatable :: ylm(:, : )
    real,              allocatable :: dKernMTGPts(:, :, :)
    complex,           allocatable              :: grVxcIRKern(:)

    complex,        allocatable                 :: psqPhon1CoulVar(:, :)
    complex,           allocatable              :: rho0IRpw(:, :)
    real,       allocatable :: r2Rho0MT(:, :, :, :)
    complex,       allocatable :: r2GrRho0MT(:, :, :, :)
    integer :: idir
    integer :: iatom
    integer :: ieqat
    integer :: oqn_l
    integer :: lm_pre
    integer :: imem
    integer :: mqn_m
    integer :: lm
    integer :: ptsym
    logical :: grRhoTermSw
    complex, allocatable :: w_grVeff0IR(:, :)
!    complex, allocatable :: grRhoFoo(:, :, :, :)
    complex, allocatable :: grRho0IR(:, :)
    integer :: iG
    real :: Gext(3)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!    allocate(rho0IRPW(ngdp, 3), rho0MTSH(atoms%jmtd, (atoms%lmaxd + 2)**2, atoms%nat, 3))
!    rho0IRPW = cmplx(0.0, 0.0)
!    rho0MTSH = cmplx(0.0, 0.0)
!    call convertStar2G(rho0IR(:, 1), rho0IRPW(:, 1), stars, ngdp, gdp)
!      rho0MTsh = 0
!      iatom = 0
!      do itype = 1, atoms%ntype
!        do ieqat = 1, atoms%neq(itype)
!          iatom = iatom + 1
!          ptsym = atoms%ntypsy(iatom)
!          do ilh = 0, lathar%nlh(ptsym)
!            oqn_l = lathar%llh(ilh, ptsym)
!            do imem = 1, nmem_atom(ilh, iatom)
!              mqn_m = mlh_atom(imem, ilh, iatom)
!              lm = oqn_l * (oqn_l + 1) + 1 + mqn_m
!              do imesh = 1, atoms%jri(itype)
!                rho0MTsh(imesh, lm, iatom, 1) = rho0MTsh(imesh, lm, iatom, 1) + rho0MT(imesh, ilh, itype, 1) * clnu_atom(imem, ilh, iatom)
!              end do ! imesh
!            end do ! imem
!          end do ! ilh
!        end do
!      end do
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !todo move to input tests
  if ( testUnfoldStarsSw ) then
    ! Calculate every g component of the ir potential from the original unperturbed FLEUR potential expanded in stars. These can be
    ! compared to a similiar system where all symmetries have been manually deleted in sym.out except for the idendity.
    call checkIRSta2gconv (V0Fleur%vpw_uw, stars, input, ngdp, gdp, logUnit)
  end if

!  if ( testCompareGrVeff0FleurSw .or. testContGrVeff0Sw .or. testVeff1Sw .or. testGradLhExpandFuncSw .or. testRadDerivativeSw .or. testGauntCoeffSw .or. testVeff1IRqLatPeriod ) then
!    write(*, *)
!    write(*, '(a)') 'Initiating potential test(s)...'
!    write(*, '(a)') '--------------------------------'
!    write(logUnit, *)
!  write ( logUnit, * )
!    write(logUnit, '(a)') 'Potential test(s)'
!    write(logUnit, '(a)') '*****************'
!    write(logUnit, *)
!  else
!    write(*, '(a)') '----------------------------'
!    write(*, '(a)') 'DISABLED potential test(s)!'
!    write(*, '(a)') '----------------------------'
!    write ( logUnit, * )
!    write(logUnit, '(a)') 'DISABLED potential tests!'
!    write(logUnit, '(a)') '*************************'
!    return
!  end if

!  if ( testCompareGrVeff0FleurSw .or. testContGrVeff0Sw ) then
    allocate( grRho0IR(ngdp, 3), rho0IRpw(ngdp, 1) )
    rho0IRpw(:, :) = cmplx(0., 0.)
    call convertStar2G( rho0IR(:, 1), rho0IRpw(:, 1), stars, ngdp, gdp )
    grRho0IR(:, :) = cmplx(0., 0.)
    do idir = 1, 3
      do iG = 1, ngdp
        Gext(1:3) = matmul(cell%bmat(1:3, 1:3), gdp(1:3, iG))
        grRho0IR(iG, idir)  = iu * Gext(idir) * rho0IRpw(iG, 1)
      end do ! iG
    end do ! idir
        allocate( rho0MTsh( atoms%jmtd, (atoms%lmaxd + 1)**2, atoms%nat, 1) )
        rho0MTsh(:, :, :, :) = cmplx(0., 0.)
        iatom = 0
        do itype = 1, atoms%ntype
          do ieqat = 1, atoms%neq(itype)
            iatom = iatom + 1
            ! Expand the coefficients of the lattice-harmonic given potential into spherical-harmonic coefficients for the given atom.
            ptsym = atoms%ntypsy(iatom)
            do ilh = 0, lathar%nlh(ptsym)
              oqn_l = lathar%llh(ilh, ptsym)
              lm_pre = oqn_l * (oqn_l + 1) + 1
              do imem = 1, nmem_atom(ilh, iatom)
                mqn_m = mlh_atom(imem, ilh, iatom)
                lm = lm_pre + mqn_m
                !todo one could only evaluate the vEff0MtSpH which have a contribution, i.e. the oqn_l and mqn_m which are in llh and mlh_atom
                ! maybe construct a pointer and run only over them to make it memory efficient.
                do imesh = 1, atoms%jri(itype)
                  rho0MTsh(imesh, lm, iatom, 1) = rho0MTsh(imesh, lm, iatom, 1) + rho0MT(imesh, ilh, itype, 1) &
                                                                                                        & * clnu_atom(imem, ilh, iatom)
                end do ! imesh
              end do ! imem
            end do ! ilh
          end do ! ieqat
        end do ! itype
        ! The factor r^2 has beeen divided out so that rho0MT is only rho0MT and not r^2 rho0MT as it is done in Fleur. The factor
        ! sqrt(4pi) for the zeroth component was already dividied out when constructing rho0MT in cdnmt routine in Fleur. Here to
        ! improve stability of the gradient routine we derive r2Rho0MT and divide out the r^2 again later. Doing so avoids the
        ! subtraction of small numbers close to the core.
        allocate( r2Rho0MT( atoms%jmtd, 0:lathar%nlhd, atoms%ntype, 1) )
        allocate( grRho0MT( atoms%jmtd, ( atoms%lmaxd + 1 )**2, atoms%nat, 3) )
        r2Rho0MT(:, :, :, :) = 0.
        grRho0MT(:, :, :, :) = cmplx(0., 0.)

        do itype = 1, atoms%ntype
          do ilh = 0, lathar%nlhd
            do imesh = 1, atoms%jri(itype)
              r2Rho0MT(imesh, ilh, itype, 1) = rho0MT(imesh, ilh, itype, 1) * atoms%rmsh(imesh, itype) * atoms%rmsh(imesh, itype)
            end do
          end do
        end do

        call calcGrR2FinLH( atoms, lathar, clnu_atom, nmem_atom, mlh_atom, r2Rho0MT(:, :, :, 1), r2GrRho0MT )

!        allocate(grRho0IR(ngdp, 3))
!        allocate(grRhoFoo(atoms%jmtd, (atoms%lmaxd + 2)**2, 3, atoms%nat))
!        grRho0IR = cmplx(0., 0.)
!        grRhoFoo = cmplx(0., 0.)
      !  do iG = 1, ngdp
      !      write(4000, '(2i8,2f15.8)') iG, idir, rho0IR(stars%ig(gdp(1, iG), gdp(2, iG), gdp(3, iG)), 1)
      !  end do ! iG
!        do idir = 1, 3
!          do iG = 1, ngdp
!            Gext = matmul( cell%bmat, gdp(:, iG) )
!            grRho0IR(iG, idir) = iu * Gext(idir) * rho0IR(stars%ig(gdp(1, iG), gdp(2, iG), gdp(3, iG)), 1) &
!              & * stars%rgphs(gdp(1, iG), gdp(2, iG), gdp(3, iG))
!            write(4000, '(2i8,2f15.8)') iG, idir, grRho0IR(iG, idir)
!          end do
!        end do
        do idir = 1, 3
          iatom = 0
          do itype = 1, atoms%ntype
            do ieqat = 1, atoms%neq(itype)
              iatom = iatom + 1
              do oqn_l = 0, atoms%lmax(itype)
                lm_pre = oqn_l * (oqn_l + 1) + 1
                do mqn_m = -oqn_l, oqn_l
                  lm = lm_pre + mqn_m
                  do imesh = 1, atoms%jri(itype)
                    grRho0MT(imesh, lm, iatom, idir) = r2GrRho0MT(imesh, lm, iatom, idir) / atoms%rmsh(imesh, itype)**2
!                    grRhoFoo(imesh, lm, idir, iatom) = grRho0MT(imesh, lm, iatom, idir)
!                    write(4001, '(4i8,2f15.8)') imesh, lm, idir, iatom, grRhoFoo(imesh, lm, idir, iatom)
                  end do ! imesh
                end do ! mqn_m
              end do ! oqn_l
            end do ! ieqat
          end do ! itype
        end do ! idir
        !write(100, *) 'foo'
           ! call checkjuPhPots(30, atoms, ngdp, gdp, cell, lathar, sym, grRho0IR, grRhoFoo, 100)
    call calcIRdVxcKern(stars, gdp, ngdp, rho0IR(:, 1), grVxcIRKern) ! weg
    call calcMTdVxcKern(atoms, dimens, lathar, rho0MT(:, :, :, 1), nmem_atom, clnu_atom, mlh_atom, gausWts, ylm, dKernMTGPts)
    testGoldstein = .false.
    grRhoTermSw = .true.
    call GenGrVeff0( atoms, cell, dimens, stars, ngdp, harSw, extSw, xcSw, gdp, rho0IRpw, rho0MTsh, grRho0IR, grRho0MT, gausWts, ylm,     &
      & dKernMTGPts, grVxcIRKern, testGoldstein, grRhoTermSw, grVeff0IR, grVeff0MT )
!  end if


    if (.false.) then
      call checkGradPot0sXSF(sym, cell, input, lathar, stars, atoms, grVeff0IR, ngdp, gdp, mlh_atom, nmem_atom, clnu_atom, rho0MT, grVeff0MT, paPoX, paPoY, paPoZ, harSw, extSw, xcSw, dimens)
    end if
  if ( testCompareGrVeff0FleurSw ) then
    write(*, '(2x,a)') 'Performing CheckGradPot0s'
    ! Compare gradient of the coulomb potential calculated analytically with the Weinert method
    call checkGradPot0s( sym, cell, input, lathar, stars, atoms, grVeff0IR, ngdp, gdp, mlh_atom, nmem_atom, clnu_atom, rho0MT,   &
      & grVeff0MT, paPoX, paPoY, paPoZ, harSw, extSw, xcSw )
    write ( logUnit, '(a)' ) '---------------------------------------------------------------------------------------------------&
      &-------'
    write ( logUnit, '(a)' ) 'Comparison between analytical gradient and Weinert method gradient of the unperturbed effective pot&
      &ential.'
    write ( logUnit, '(a)' ) '  --> Check by plotting the just created data within pathPot, pathPotCont using plotPotential.py'
    write ( logUnit, '(a)' ) '---------------------------------------------------------------------------------------------------&
      &-------'
    write ( logUnit, * )
  else
    write(*, '(2x,a)') 'DISABLED CheckGradPot0s!'
    write(logUnit, *)
    write ( logUnit, '(a)' ) 'Comparison between analytical gradient and Weinert method gradient of the unperturbed effective pot&
      &ential.'
    write ( logUnit, '(a)' ) '---------------------------------------------------------------------------------------------------&
      &-------'
    write ( logUnit, '(a)' ) '                                                                                                   &
      &      |_DISABLED!'
  end if

  if ( testContGrVeff0Sw ) then
    write(*, '(2x,a)') 'Performing CheckjuPhPots...'
    write ( logUnit, '(a)' ) "Testing the continuity of Weinert generated potential's gradient (warped and unwarped)."
    write ( logUnit, '(a)' ) "---------------------------------------------------------------------------------------"
    call checkjuPhPots(noPtsCon, atoms, ngdp, gdp, cell, lathar, sym, grVeff0IR, grVeff0MT, logUnit)
    ! Warp interstitial second-order external potential
    allocate(w_grVeff0IR(ngdp, 3))
    w_grVeff0IR = cmplx(0.0, 0.0)
    do idir = 1, 3
      call warpIRPot(stars, ngdp, idir, gdp, grVeff0IR, w_grVeff0IR(:, idir))
    end do ! idirC
    write ( logUnit, * )
    write ( logUnit, '(a)' ) 'WARPED VERSION IN IR, IS KNOWN TO BE NOT VERY CONTINIOUS'
    call checkjuPhPots(noPtsCon, atoms, ngdp, gdp, cell, lathar, sym, w_grVeff0IR, grVeff0MT, logUnit)
    write ( logUnit, * )
  else
    write(*, '(2x,a)') 'DISABLED CheckjuPhPots!'
      write(logUnit, *)
    write ( logUnit, '(a)' ) "Testing the continuity of Weinert generated potential's gradient (warped and unwarped)."
    write ( logUnit, '(a)' ) "---------------------------------------------------------------------------------------"
    write ( logUnit, '(a)' ) "                                                                                      |_DISABLED"
  end if

  !todo we might check also the continuity of veff1

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!todo make this not needed
!  allocate( myrho0MT(atoms%jmtd, 0:lathar%nlhd, atoms%ntype, input%jspins))
!    do ispin = 1, input%jspins
!      do itype = 1, atoms%ntype
!        do ilh = 0, lathar%nlhd
!          do imesh = 1, atoms%jri(itype)
!            myrho0MT(imesh, ilh, itype, ispin) = rho0MT(imesh, ilh, itype, ispin) * atoms%rmsh(imesh, itype) * atoms%rmsh(imesh, itype)
!        end do
!      end do
!    end do
!  end do
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  if ( testVeff1Sw ) then
    write(*, '(2x,a)') 'Performing CheckVCoul1...'
    ! Reproduce the Fleur Coulomb potential with the juPhon routines for calculating the linear variation of the Coulomb potential
    call checkVCoul1( atoms, lathar, cell, stars, dimens, input, rho0MT, clnu_atom, nmem_atom, mlh_atom, gdp, ngdp, rho0IR, logUnit )
  else
    write(*, '(2x,a)') 'DISABLED CheckVCoul1!'
    write(logUnit, *)
    write(logUnit, '(a)') 'Calculate FLEUR Coulomb potential with V_coulomb1 routines test'
    write(logUnit, '(a)') '---------------------------------------------------------------'
    write (logUnit, '(a)')'                                                              |__ DISABLED!'
  end if

  if ( testGradLhExpandFuncSw ) then
    write(*, '(2x,a)') 'Performing CheckMTDerivative...'
    ! test mt derivative
    call checkMTDerivative(atoms, lathar, nmem_atom, mlh_atom, memd_atom, logUnit)
  else
    write(*, '(2x,a)') 'DISABLED CheckMTDerivative!'
    write(logUnit, *)
    write(logUnit, '(a)') 'Testing routine to calculate the muffin-tin gradient of a lattice-harmonic expanded quantity'
    write(logUnit, '(a)') '--------------------------------------------------------------------------------------------'
    write (logUnit, '(a)')'                                                                                           |__ DISABLED!'
  end if

  if ( testRadDerivativeSw ) then
    write(*, '(2x,a)') 'Performing CheckDerivative...'
  ! checks whether derivative routine works correctly
    call checkDerivative(atoms, logUnit)
  else
    write(*, '(2x,a)') 'DISABLED CheckDerivative!'
    write(logUnit, *)
    write(logUnit, '(a)') 'Test of radial derivative'
    write(logUnit, '(a)') '-------------------------'
    write(logUnit, '(a)') '                        |_DISABLED '
  end if

  if ( testGauntCoeffSw ) then
    write(*, '(2x,a)') 'Performing TestGauntCoeff...'
     call testGauntCoeff(atoms, logUnit)
  else
    write(*, '(2x,a)') 'DISABLED TestGauntCoeff!'
    write(logUnit, *)
    write(logUnit, '(a)') 'Test for Gaunt coefficients routines'
    write(logUnit, '(a)') '------------------------------------'
    write(logUnit, '(a)') '                                   |__ DISABLED!'
  end if

  if (testVeff1IRqLatPeriod) then
    write(*, '(2x,a)') 'Performing testVeff1IRqLatPeriod...'
    call TestQLatticePeriodicity( atoms, stars, cell, lathar, dimens, sym, input, qpts, ngdp, memd_atom, logUnit, rho0IR, rho0MT, mlh_atom,     &
                                                                                                       & nmem_atom, clnu_atom, gdp )
  else
    write(*, '(2x,a)') 'DISABLED testVeff1IRqLatPeriod!'
  end if

end subroutine Test1stOrdPotentials

  !todo unite with checkDensnPot from jpTestInput
  ! this only checks the continuity of the potential and writes it to the log file
  ! This routine also supports lmax functions instead of lmax + 1 and complex functions
subroutine checkjuPhPots1(noPts2chCont, atomsT, ngdp, gdp, cellT, sphharT, symT, fIR, fMufT, logFileUnit, qpt)

    use m_types
    use m_cotra
    use m_ylm_old
    use mod_juPhonUtils
    use m_jPConstants
    implicit none
    type(t_atoms),        intent(in) :: atomsT
    type(t_cell),         intent(in) :: cellT
    type(t_sphhar),       intent(in) :: sphharT
    type(t_sym),          intent(in) :: symT
    integer,              intent(in) :: noPts2chCont !number of points to check continuity
    integer,              intent(in) :: ngdp
    integer,              intent(in) :: gdp(:, :)
    integer,              intent(in)            :: logFileUnit ! why not as intent(in)?
    complex,             intent(in) :: fIR(:, :)
    complex,             intent(in) :: fMufT(:, :, :, :)
    real, intent(in) :: qpt(:)
!    real,                 intent(in)  :: rmt(:) !MT-radii of atom type itype
!    integer,              intent(in)  :: ntype ! number of atom types
!    integer,              intent(in)  :: nop !number of symmetry operations
!    integer,              intent(in)  :: nq3 ! number of stars
!    logical,              intent(in)  :: symor !are symmetry operations symorphic or noct
!    integer,              intent(in)  :: kv3(:, :) !reciprocal g-vector of star
!    integer,              intent(in)  :: mrot(:, :, :) !symmetry operations, rotation matrices
!    real,              intent(in)  :: tau(:, :) !symmetry operations, translation vectors
!    integer,              intent(in)  :: invtab(:) !lists which symmetry operations are inverse to symmetry operations with indesx isym
!    real,                 intent(in)  :: bmat(:, :) !reciprocal Bravais matrix transposed
!    real,                 intent(in)  :: amat(:, :) !Bravais matrix
!    real,                 intent(in)  :: pos(:, :) !Postions of atoms in Brillouin zone
!    integer,              intent(in)  :: nmem(:, :) !number of lattice harmonics
!    integer,              intent(in)  :: llh(:, :) !ls of  lattice harmonics
!    complex,              intent(in)  :: clnu(:, :, :) !phasefactors in linear combination of lattice harmonics
!    logical,              intent(in)  :: density ! is function density?
!    complex,              intent(in)  :: fpw(:, :)
!    real,                 intent(in)  :: fr(:, :, :, :)
!    integer,              intent(in)  :: nstr(:)
!    integer,              intent(in)  :: nat
!    integer,              intent(in)  :: ntypsy(:)
!    integer,              intent(in)  :: ngopr(:)
!    integer,              intent(in)  :: lmax(:)
!    integer,              intent(in)  :: nlh(:)
!    integer,              intent(in)  :: mlh(:, :, :)
!    integer,              intent(in)  :: jri(:)
!    integer,              intent(in)  :: lmaxd
!    integer,              intent(in)  :: neq(:)
!
!
    integer                           :: oqn_l
    integer                           :: mqn_m
    integer                           :: lm_lcontr
    integer                           :: lm
    integer                           :: irandPt
    integer                           :: ieq !loop variable
    real                              :: euclidNorm !variable to store eculidian L2 norm of vector
    real,   allocatable               :: randPtsCart(:, :, :) !cartesian dimension, number of Pts, itype, random points
    real,   allocatable               :: randPtsGrid(:, :, :) !random points on grid, coordinates, number of points
    real,   allocatable               :: randPtsMTLoc(:, :, :) !random points in local atomic coordinate (coordinates, number of points, ntype)
    integer                           :: itype !certain atom type
!    complex                           :: starAtRanPt(nq3) !value of star at certain point ! ng3 = nq3
    complex,   allocatable               :: sumfuncvalI(:, :, :) !sum of function values coming from interstitial
    integer                         :: iGvec
    integer                          :: idirec
!    real                              :: scaling !scaling factor which makes difference between density and potential or wavefunction
!    integer                           :: icoord !loop variable
    integer                           :: imatmulc, imatmulr !loop variable
    real                              :: rvec(3) !randoom point
    real                              :: rvec_int(3)
    real  :: rotatedVector(3) !rotated vector after symmetry operation
!    integer                           :: istar !loop variable
    integer                           :: iatom
!    integer                           :: symAt_temp
    integer                           :: symOpr_temp
    complex                           :: ylm((atomsT%lmaxd + 1)**2) !TODO eigentlich lmax(itype)
    complex,    allocatable              :: sumOfMTFuncVal(:, :,:)
!    integer                           :: ilath
!    real                              :: linCombBas
!    integer                           :: l_temp
!    integer                           :: imem
!    integer                           :: lm_temp
    real                              :: av(3), dms(3), rms(3)
!    integer                           :: ispin

    allocate(randPtsCart(3, noPts2chCont, atomsT%nat), randPtsGrid(3, noPts2chCont, atomsT%ntype), randPtsMTLoc(3, noPts2chCont, atomsT%ntype) )

    allocate(sumfuncvalI(3, noPts2chCont,atomsT%nat), sumOfMTFuncVal(3, noPts2chCont,atomsT%nat))
   !Generate noPts2chCont random points on the MT sphere
    randPtsCart(:, :, :) = cmplx(0., 0.)
    iatom = 0
    do itype = 1, atomsT%ntype
      do ieq = 1, atomsT%neq(itype)
        iatom = iatom + 1
        call sphpts(randPtsCart(:, :, iatom), noPts2chCont, atomsT%rmt(itype), atomsT%pos(:, iatom))
!      do irandPt = 1, noPts2chCont
!        call random_number(randPtsCart(:, irandPt, itype))
!        euclidNorm = norm2(randPtsCart(:, irandPt, itype))
!        randPtsCart(:, irandPt, itype) = randPtsCart(:, irandPt, itype) / euclidNorm * atomsT%rmt(itype)
      end do
    end do


    !Evaluate Stars for this random points at each MT sphere boundary in the unit cell
    sumfuncValI = cmplx(0.0) !todo this is probably unneccessary, because it is already set to 0 in the loop
    iatom       = 0
    do itype = 1, atomsT%ntype
      do ieq = 1, atomsT%neq(itype)
        iatom = iatom + 1
        do irandPt = 1, noPts2chCont
!          rvec = atomsT%pos(:,iatom) + randPtsCart(:,irandPt,itype)
          rvec(:) = randPtsCart(:, irandPt, iatom)
          ! transform rvec in internal, real-space coordinates --> rvec_int
          call cotra1(rvec, rvec_int, cellT%bmat) !todo do this inline from now on!!!
          ! Evaluate stars at random point
          sumfuncValI(:, irandPt,iatom) = cmplx(0., 0.) !todo somehow this does not make sense
          do iGvec = 1, ngdp !nq3 is number of stars
            do idirec = 1, 3
              sumfuncValI(idirec, irandPt,iatom) = sumfuncvalI(idirec, irandPt,iatom) + fIR(iGvec, idirec) * exp(iu * tpi * dot_product(gdp(:, iGvec) + qpt(:), rvec_int(:))) !jsp noch einfügen, warum Multiplikation with G-vectors
            end do ! idirec
          end do  !iGvec
        end do  !irandPt
      end do  !ieq
    end do  !itype

    if (all(qpt(:) < 1e-6)) then
      if (any(aimag(sumfuncValI(:, :, :)) > 1e-8)) then
        write(*, *) "Imaginary part of observable potential should be 0"
        NOstopNO
      end if
    end if


    write(*, *) 'Attention: the MT quantitiy goes until lmax not lmax + 1'
    iatom = 0
    sumOfMTFuncVal = cmplx(0., 0.)
    do itype = 1, atomsT%ntype
      !Evaluate MT sphere function at same random points
      randPtsGrid = 0.0
      do ieq = 1, atomsT%neq(itype)
        iatom       = iatom + 1
        !symOpr_temp = 1!atomsT%ngopr(iatom) ! symmetry operation mapping local to global coordinate system
        do irandPt = 1, noPts2chCont
          !rvec = randPtsCart(:,irandPt,itype)
          ! We have to subtract the atom position to get into a local coordinate system of atom iatom
          rvec = randPtsCart(:,irandPt,iatom) - atomsT%pos(:, iatom)
         ! if (symOpr_temp /= 1) then
         !   !Transform into internal real space coordinates
         !   call cotra1(rvec, rvec_int, cellT%bmat) !allocate and deallocate randPtsGrid
         !   do imatmulr = 1, 3 !todo optimize this
         !     rotatedVector(imatmulr) = 0.0
         !     do imatmulc = 1, 3
         !       rotatedVector(imatmulr) = rotatedVector(imatmulr) + symT%mrot(imatmulr, imatmulc, symOpr_temp) *rvec_int(imatmulc)
         !     end do
         !   end do
         !   call cotra0(rotatedVector, rvec, cellT%amat)
         ! end if
    !      call ylmnorm_init(atomsT%lmaxd + 1)
          call ylm4(atomsT%lmax(itype), rvec, ylm)
    !      call ylmnorm_init(atomsT%lmaxd)

          do oqn_l = 0, atomsT%lmax(itype) !+ 1
            lm_lcontr = oqn_l * (oqn_l + 1) + 1
            do mqn_m = -oqn_l, oqn_l
              lm = lm_lcontr + mqn_m
              do idirec = 1, 3
                sumOfMTFuncVal(idirec, irandPt, iatom) = sumOfMTFuncVal(idirec, irandPt, iatom) + fMufT(atomsT%jri(itype), lm, idirec, iatom) * ylm(lm)
              end do ! idir
            end do
          end do
        end do
      end do
    end do

    if (all(qpt(:) < 1e-6)) then
      if (any(aimag(sumOfMTFuncVal(:, :, :)) > 1e-8)) then
        write(*, *) "Imaginary part of observable potential should be 0"
        NOstopNO
      end if
    end if

    write ( logFileUnit, '(a)' ) "Continuity of the unperturbed potential's gradient:"
    write(logFileUnit,*)         '---------------------------------------------------'
    iatom = 0
    do itype = 1,atomsT%ntype
      do ieq = 1,atomsT%neq(itype)
        iatom = iatom + 1
        write(logFileUnit,*) '  Atom:',iatom
        write(logFileUnit,'(a)') '  coordinate(cartesian)      IR x-comp.  MT x-comp.     IR y-comp.  MT y-comp.     IR z-comp.  MT z-comp.'
        do irandPt = 1,noPts2chCont
            write(logFileUnit,'(3f8.4,3x,3(4f12.8, 3x))') randPtsCart(:,irandPt,iatom),sumfuncValI(1, irandPt,iatom),sumOfMTFuncVal(1, irandPt,iatom), sumfuncValI(2, irandPt,iatom),sumOfMTFuncVal(2, irandPt,iatom), sumfuncValI(3, irandPt,iatom),sumOfMTFuncVal(3, irandPt,iatom)
        end do
        do idirec = 1, 3
          call fitchkNeg(real(sumfuncvalI(idirec, :,iatom)), real(sumOfMTFuncVal(idirec, :,iatom)), noPts2chCont, av(idirec), rms(idirec), dms(idirec))
        end do
        write(logFileUnit, *)
      !  write(logFileUnit,'(a)') '  average deviation          IR x-comp.  MT x-comp.     IR y-comp.  MT y-comp.     IR z-comp.  MT z-comp.'
        write(logFileUnit,'(2x, A,11x, 3(f15.8, 15x))'   )  'real average absolute value interstitial :',av(1), av(2), av(3)
        write(logFileUnit,'(2x,A,6f15.8,A/)')               'real rms(x, y, z), real dms(x, y, z):          ',rms(1), rms(2), rms(3), dms(1), dms(2), dms(3),' in %'
        av(:) = 0.
        rms(:) = 0.
        dms(:) = 0.
        do idirec = 1, 3
          call fitchkNeg(aimag(sumfuncvalI(idirec, :,iatom)), aimag(sumOfMTFuncVal(idirec, :,iatom)), noPts2chCont, av(idirec), rms(idirec), dms(idirec))
        end do
        write(logFileUnit, *)
      !  write(logFileUnit,'(a)') '  average deviation          IR x-comp.  MT x-comp.     IR y-comp.  MT y-comp.     IR z-comp.  MT z-comp.'
        write(logFileUnit,'(2x, A,11x, 3(f15.8, 15x))'   )  'aimag average absolute value interstitial :',av(1), av(2), av(3)
        write(logFileUnit,'(2x,A,6f15.8,A/)')               'aimag rms(x, y, z), dms(x, y, z):          ',rms(1), rms(2), rms(3), dms(1), dms(2), dms(3),' in %'
      end do
    end do
  end subroutine checkjuPhPots1
  !todo unite with checkDensnPot from jpTestInput
  ! this only checks the continuity of the potential and writes it to the log file
subroutine checkjuPhPots(noPts2chCont, atomsT, ngdp, gdp, cellT, sphharT, symT, fIR, fMufT, logFileUnit)

    use m_types
    use m_cotra
    use m_ylm_old
    use mod_juPhonUtils
    use m_jPConstants
    implicit none
    type(t_atoms),        intent(in) :: atomsT
    type(t_cell),         intent(in) :: cellT
    type(t_sphhar),       intent(in) :: sphharT
    type(t_sym),          intent(in) :: symT
    integer,              intent(in) :: noPts2chCont !number of points to check continuity
    integer,              intent(in) :: ngdp
    integer,              intent(in) :: gdp(:, :)
    integer,              intent(in)            :: logFileUnit ! why not as intent(in)?
    complex,             intent(in) :: fIR(:, :)
    complex,             intent(in) :: fMufT(:, :, :, :)
!    real,                 intent(in)  :: rmt(:) !MT-radii of atom type itype
!    integer,              intent(in)  :: ntype ! number of atom types
!    integer,              intent(in)  :: nop !number of symmetry operations
!    integer,              intent(in)  :: nq3 ! number of stars
!    logical,              intent(in)  :: symor !are symmetry operations symorphic or noct
!    integer,              intent(in)  :: kv3(:, :) !reciprocal g-vector of star
!    integer,              intent(in)  :: mrot(:, :, :) !symmetry operations, rotation matrices
!    real,              intent(in)  :: tau(:, :) !symmetry operations, translation vectors
!    integer,              intent(in)  :: invtab(:) !lists which symmetry operations are inverse to symmetry operations with indesx isym
!    real,                 intent(in)  :: bmat(:, :) !reciprocal Bravais matrix transposed
!    real,                 intent(in)  :: amat(:, :) !Bravais matrix
!    real,                 intent(in)  :: pos(:, :) !Postions of atoms in Brillouin zone
!    integer,              intent(in)  :: nmem(:, :) !number of lattice harmonics
!    integer,              intent(in)  :: llh(:, :) !ls of  lattice harmonics
!    complex,              intent(in)  :: clnu(:, :, :) !phasefactors in linear combination of lattice harmonics
!    logical,              intent(in)  :: density ! is function density?
!    complex,              intent(in)  :: fpw(:, :)
!    real,                 intent(in)  :: fr(:, :, :, :)
!    integer,              intent(in)  :: nstr(:)
!    integer,              intent(in)  :: nat
!    integer,              intent(in)  :: ntypsy(:)
!    integer,              intent(in)  :: ngopr(:)
!    integer,              intent(in)  :: lmax(:)
!    integer,              intent(in)  :: nlh(:)
!    integer,              intent(in)  :: mlh(:, :, :)
!    integer,              intent(in)  :: jri(:)
!    integer,              intent(in)  :: lmaxd
!    integer,              intent(in)  :: neq(:)
!
!
    integer                           :: oqn_l
    integer                           :: mqn_m
    integer                           :: lm_lcontr
    integer                           :: lm
    integer                           :: irandPt
    integer                           :: ieq !loop variable
    real                              :: euclidNorm !variable to store eculidian L2 norm of vector
    real,   allocatable               :: randPtsCart(:, :, :) !cartesian dimension, number of Pts, itype, random points
    real,   allocatable               :: randPtsGrid(:, :, :) !random points on grid, coordinates, number of points
    real,   allocatable               :: randPtsMTLoc(:, :, :) !random points in local atomic coordinate (coordinates, number of points, ntype)
    integer                           :: itype !certain atom type
!    complex                           :: starAtRanPt(nq3) !value of star at certain point ! ng3 = nq3
    complex,   allocatable               :: sumfuncvalI(:, :, :) !sum of function values coming from interstitial
    integer                         :: iGvec
    integer                          :: idirec
!    real                              :: scaling !scaling factor which makes difference between density and potential or wavefunction
!    integer                           :: icoord !loop variable
    integer                           :: imatmulc, imatmulr !loop variable
    real                              :: rvec(3) !randoom point
    real                              :: rvec_int(3)
    real  :: rotatedVector(3) !rotated vector after symmetry operation
!    integer                           :: istar !loop variable
    integer                           :: iatom
!    integer                           :: symAt_temp
    integer                           :: symOpr_temp
    complex                           :: ylm((atomsT%lmaxd + 1)**2) !TODO eigentlich lmax(itype)
    complex,    allocatable              :: sumOfMTFuncVal(:, :,:)
!    integer                           :: ilath
!    real                              :: linCombBas
!    integer                           :: l_temp
!    integer                           :: imem
!    integer                           :: lm_temp
    real                              :: av(3), dms(3), rms(3)
!    integer                           :: ispin

    allocate(randPtsCart(3, noPts2chCont, atomsT%nat), randPtsGrid(3, noPts2chCont, atomsT%ntype), randPtsMTLoc(3, noPts2chCont, atomsT%ntype) )

    allocate(sumfuncvalI(3, noPts2chCont,atomsT%nat), sumOfMTFuncVal(3, noPts2chCont,atomsT%nat))
   !Generate noPts2chCont random points on the MT sphere
    randPtsCart(:, :, :) = cmplx(0., 0.)
    iatom = 0
    do itype = 1, atomsT%ntype
      do ieq = 1, atomsT%neq(itype)
        iatom = iatom + 1
        call sphpts(randPtsCart(:, :, iatom), noPts2chCont, atomsT%rmt(itype), atomsT%pos(:, iatom))
!      do irandPt = 1, noPts2chCont
!        call random_number(randPtsCart(:, irandPt, itype))
!        euclidNorm = norm2(randPtsCart(:, irandPt, itype))
!        randPtsCart(:, irandPt, itype) = randPtsCart(:, irandPt, itype) / euclidNorm * atomsT%rmt(itype)
      end do
    end do


    !Evaluate Stars for this random points at each MT sphere boundary in the unit cell
    sumfuncValI = cmplx(0.0) !todo this is probably unneccessary, because it is already set to 0 in the loop
    iatom       = 0
    do itype = 1, atomsT%ntype
      do ieq = 1, atomsT%neq(itype)
        iatom = iatom + 1
        do irandPt = 1, noPts2chCont
!          rvec = atomsT%pos(:,iatom) + randPtsCart(:,irandPt,itype)
          rvec(:) = randPtsCart(:, irandPt, iatom)
          ! transform rvec in internal, real-space coordinates --> rvec_int
          call cotra1(rvec, rvec_int, cellT%bmat) !todo do this inline from now on!!!
          ! Evaluate stars at random point
          sumfuncValI(:, irandPt,iatom) = cmplx(0., 0.) !todo somehow this does not make sense
          do iGvec = 1, ngdp !nq3 is number of stars
            do idirec = 1, 3
              sumfuncValI(idirec, irandPt,iatom) = sumfuncvalI(idirec, irandPt,iatom) + fIR(iGvec, idirec) * exp(iu * tpi * dot_product(gdp(:, iGvec), rvec_int(:))) !jsp noch einfügen, warum Multiplikation with G-vectors
            end do ! idirec
          end do  !iGvec
        end do  !irandPt
      end do  !ieq
    end do  !itype

    if (any(aimag(sumfuncValI(:, :, :)) > 1e-8)) then
      write(*, *) "Imaginary part of observable IR potential should be 0"
      NOstopNO
    end if


    iatom = 0
    sumOfMTFuncVal = cmplx(0., 0.)
    do itype = 1, atomsT%ntype
      !Evaluate MT sphere function at same random points
      randPtsGrid = 0.0
      do ieq = 1, atomsT%neq(itype)
        iatom       = iatom + 1
        !symOpr_temp = 1!atomsT%ngopr(iatom) ! symmetry operation mapping local to global coordinate system
        do irandPt = 1, noPts2chCont
          !rvec = randPtsCart(:,irandPt,itype)
          ! We have to subtract the atom position to get into a local coordinate system of atom iatom
          rvec = randPtsCart(:,irandPt,iatom) - atomsT%pos(:, iatom)
         ! if (symOpr_temp /= 1) then
         !   !Transform into internal real space coordinates
         !   call cotra1(rvec, rvec_int, cellT%bmat) !allocate and deallocate randPtsGrid
         !   do imatmulr = 1, 3 !todo optimize this
         !     rotatedVector(imatmulr) = 0.0
         !     do imatmulc = 1, 3
         !       rotatedVector(imatmulr) = rotatedVector(imatmulr) + symT%mrot(imatmulr, imatmulc, symOpr_temp) *rvec_int(imatmulc)
         !     end do
         !   end do
         !   call cotra0(rotatedVector, rvec, cellT%amat)
         ! end if
        !  call ylmnorm_init(atomsT%lmaxd + 1)
          call ylm4(atomsT%lmax(itype), rvec, ylm)
        !  call ylmnorm_init(atomsT%lmaxd)

          do oqn_l = 0, atomsT%lmax(itype)
            lm_lcontr = oqn_l * (oqn_l + 1) + 1
            do mqn_m = -oqn_l, oqn_l
              lm = lm_lcontr + mqn_m
              do idirec = 1, 3
                sumOfMTFuncVal(idirec, irandPt, iatom) = sumOfMTFuncVal(idirec, irandPt, iatom) + fMufT(atomsT%jri(itype), lm, idirec, iatom) * ylm(lm)
              end do ! idir
            end do
          end do
        end do
      end do
    end do

    if (any(aimag(sumOfMTFuncVal(:, :, :)) > 1e-8)) then
      write(*, *) "Imaginary part of observable MT potential should be 0"
      NOstopNO
    end if

    write ( logFileUnit, '(a)' ) "Continuity of the unperturbed potential's gradient:"
    write(logFileUnit,*)         '---------------------------------------------------'
    iatom = 0
    do itype = 1,atomsT%ntype
      do ieq = 1,atomsT%neq(itype)
        iatom = iatom + 1
        write(logFileUnit,*) '  Atom:',iatom
        write(logFileUnit,'(a)') '  coordinate(cartesian)      IR x-comp.  MT x-comp.     IR y-comp.  MT y-comp.     IR z-comp.  MT z-comp.'
        do irandPt = 1,noPts2chCont
            write(logFileUnit,'(3f8.4,3x,3(2f12.8, 3x))') randPtsCart(:,irandPt,iatom),real(sumfuncValI(1, irandPt,iatom)),real(sumOfMTFuncVal(1, irandPt,iatom)), real(sumfuncValI(2, irandPt,iatom)),real(sumOfMTFuncVal(2, irandPt,iatom)), real(sumfuncValI(3, irandPt,iatom)),real(sumOfMTFuncVal(3, irandPt,iatom))
        end do
        do idirec = 1, 3
          call fitchkNeg(real(sumfuncvalI(idirec, :,iatom)), real(sumOfMTFuncVal(idirec, :,iatom)), noPts2chCont, av(idirec), rms(idirec), dms(idirec))
        end do
        write(logFileUnit, *)
      !  write(logFileUnit,'(a)') '  average deviation          IR x-comp.  MT x-comp.     IR y-comp.  MT y-comp.     IR z-comp.  MT z-comp.'
        write(logFileUnit,'(2x, A,11x, 3(f15.8, 15x))'   )  'average absolute value interstitial :',av(1), av(2), av(3)
        write(logFileUnit,'(2x,A,6f15.8,A/)')               'rms(x, y, z), dms(x, y, z):          ',rms(1), rms(2), rms(3), dms(1), dms(2), dms(3),' in %'
      end do
    end do
  end subroutine checkjuPhPots

  ! This is a fork of fitchk to handle quantities where the interstitial varies a lot so the average is almost zero
  subroutine fitchkNeg(f1, f2, n, av, rms, dmx)

    implicit none

    integer, intent(in) :: n

    real,    intent(in) :: f1(n)
    real,    intent(in) :: f2(n)

    real,    intent(out) :: av
    real,    intent(out) :: dmx
    real,    intent(out) :: rms


    real                 :: d
    integer              :: i

    av = 0.
    rms = 0.
    dmx = 0.
    do i = 1,n
       av = av + abs(f1(i))
       d = (f1(i)-f2(i))**2
       dmx = max(d,dmx)
       rms = rms + d
    end do
    av = av/n
    if (abs(av).LT.1.e-30) then
       rms = 0.
       dmx = 0.
       return
    end if
    rms = sqrt(rms/n)/av*100.
    dmx = sqrt(dmx)/av*100.

  end subroutine fitchkNeg


  ! Writes plot files to compare grVeff which has been generated by the Weinert method and by the gradient
  subroutine checkGradPot0s(symT, cellT, inputT, sphharT, starsT, atomsT, vGradCoul0IR, ngdp, gdp, mlh_atom, nmem_atom, clnu_atom, rho0MT, gradVrjuPhon, paPoXCo, paPoYCo, paPoZCo, harSw, extSw, xcSw)

    use m_types
    use m_jPConstants, only : iu, pi, tpi, fpi, Tmatrix
    use m_gaunt
    use m_ylm_old
    use m_cotra
    use m_jpPotDensHelper
    use m_juDFT_NOstopNO

    implicit none

    type(t_sym),        intent(in) :: symT
    type(t_cell),       intent(in) :: cellT
    type(t_sphhar),     intent(in) :: sphharT
    type(t_stars),      intent(in) :: starsT
    type(t_atoms),      intent(in) :: atomsT
    type(t_input),      intent(in) :: inputT
    complex,            intent(in) :: vGradCoul0IR(:, :)
    integer,            intent(in) :: ngdp
    integer,            intent(in) :: gdp(:, :)
    integer,            intent(in) :: mlh_atom(:, 0:, :)
    integer,            intent(in) :: nmem_atom(0:, :)
    complex,            intent(in) :: clnu_atom(:, 0:, :)
    real,               intent(in) :: rho0MT(:, 0:, :, :)
    real,               intent(in) :: paPoXCo
    real,               intent(in) :: paPoYCo
    real,               intent(in) :: paPoZCo
    logical,            intent(in) :: harSw
    logical,            intent(in) :: extSw
    logical,            intent(in) :: xcSw

    integer                        :: iGvec
    real                           :: Gext(3)
    complex                        :: vCoulBenchmark(ngdp, 3)
    complex                        :: gradVcMT(atomsT%jmtd, (atomsT%lmaxd + 1)**2, atomsT%nat, 3)
    real                           :: sqr4pi3
    real                           :: Vc0nSym(atomsT%jmtd, (atomsT%lmaxd + 1 +  1)**2, atomsT%nat) ! todo is this the correct dimension
    integer                        :: iatom
    integer                        :: itype
    integer                        :: ineq
    integer                        :: symType
    integer                        :: lh
    integer                        :: oqn_l
    integer                        :: lm_temp
    integer                        :: mem
    integer                        :: mems
    integer                        :: mqn_m
    integer                        :: mqn_mpp
    integer                        :: lm
    integer                        :: lmMinus
    integer                        :: lmPlus
    real                           :: Vc0MTDerived(atomsT%jmtd)
    real                           :: tempGaunt1, tempGaunt2
    integer                        :: imesh
    complex                        :: vpw_G(ngdp)
    complex                        :: vpw_G_coul(ngdp)
    complex                        :: vpw_G_hart(ngdp)
    complex                        :: vpw_G_xc(ngdp)
    complex              :: vpwCoulUwStars1(starsT%n3d, 1)
    complex              :: vpwCoulUwStars2(starsT%n3d, 1)
    complex              :: vpwStar(starsT%n3d, 1)
    complex              :: vpwStarDiff(starsT%n3d, 1)
    complex                        :: vpwRfleur(3)
    complex                        :: vpwRjuPhon(3)
    complex                        :: exponentialjuPhon
    complex                        :: exponentialFLEUR
    complex                        :: exponential
    integer           :: isIR
    integer           :: wasIR
    integer          :: c1, c2, c3
    real                           :: ucpath(3)
    real                              :: ucpathA(3), ucpathB(3), ucpathC(3)
    real                           :: ucpathExt(3)
    real                      :: ucpathExtc(3), ucpathExta(3)
    real                           :: dx
    integer                        :: ii
    real                           :: x
    complex                        :: newvpwFLEUR(starsT%n3d, 1) ! make this declaration more general!
    complex                        :: newvpw_G(ngdp)
    complex                        :: onlyPotFLEUR
    complex                        :: onlypsqFLEUR
    integer                        :: ieq
    integer                        :: ext2brav(3)
    integer                        :: ext2bravLock(3)
    integer                        :: atomLock
    integer                        :: itypeLock
    logical                        :: pseudoPot = .true.
    logical                        :: realPot = .true.
    integer                   :: iterations
  real   :: vrFleur(atomsT%jmtd, 0:sphharT%nlhd, atomsT%ntypd, inputT%jspins)
  real   :: vrFleurDiff(atomsT%jmtd, 0:sphharT%nlhd, atomsT%ntypd, inputT%jspins)
  real   :: vrFleurxc(atomsT%jmtd, 0:sphharT%nlhd, atomsT%ntypd, inputT%jspins)
  real   :: vrTemp(atomsT%jmtd, 0:sphharT%nlhd, atomsT%ntypd, inputT%jspins)
  complex,             intent(in) :: gradVrjuPhon(:, :, :, :)
  complex,      allocatable :: gradVrFleur(:, :, :, :)
  real    :: direc(3)
  real :: direcExt(3)
  complex                                     :: ylm((atomsT%lmaxd + 1)**2)
  complex                                     :: vMTfleur(3, atomsT%jmtd)
  complex                                     :: vMTjuPhon(3, atomsT%jmtd)
  integer                          :: idirec
  real                              :: av, dmx, rms
  real                              :: IRvalue(1), MTvalue(1)
  integer                             :: jspin
  integer                               :: ilh
  real        :: potFleurMT(atomsT%jmtd)
  integer    :: symAt_temp, symOpr_temp
  real      :: rvec(3), rvec_int(3)
  integer    :: imatmulr, imatmulc
  real :: rotatedVector(3)
  integer :: ilath
  real :: linCombBas
  integer :: l_temp
  integer :: imem

logical :: fBenchAvailSw

    ! in this test the gradient of the Coulomb potential in the interstitial region calculated with Aaron's method must be equal to iG qpw
    fBenchAvailSw = .true.
    vpwStar = 0
    vpwStarDiff = 0
    vrFLEUR = 0
    vrFLEURDiff = 0

    !todo do this with pottot and potcoul
    if ( harSw .and. .not.extSw .and. .not.xcSw ) then ! Hartree potential
      call fopen( 1000, name='vpw_hart', status='old', action='read', form='unformatted' )
        read( 1000 ) vpwStar
      call fclose( 1000 )

      call fopen(1000, name='v0MTFLEUR_hart', status='old', action='read', form='unformatted')
        read(1000) vrFLEUR
      call fclose(1000)

    else if ( .not.harSw .and. extSw .and. .not.xcSw ) then ! external potential
      call fopen( 1000, name='vpw_coul', status='old', action='read', form='unformatted' )
        read( 1000 ) vpwStar
      call fclose( 1000 )

      call fopen( 1000, name='vpw_hart', status='old', action='read', form='unformatted' )
        read( 1000 ) vpwStarDiff
      call fclose( 1000 )

      vpwStar = vpwStar - vpwStarDiff


      call fopen(1000, name='v0MTFLEUR_coul', status='old', action='read', form='unformatted')
        read(1000) vrFLEUR
      call fclose(1000)

      call fopen(1000, name='v0MTFLEUR_hart', status='old', action='read', form='unformatted')
        read(1000) vrFLEURDiff
      call fclose(1000)

      vrFLEUR = vrFLEUR - vrFLEURDiff

    else if ( .not.harSw .and. .not.extSw .and. xcSw ) then ! exchange-correlation potential
      call fopen( 1000, name='vpw_xc', status='old', action='read', form='unformatted' )
        read( 1000 ) vpwStar
      call fclose( 1000 )

      call fopen(1000, name='v0MTFLEUR_xc', status='old', action='read', form='unformatted')
        read(1000) vrFLEUR
      call fclose(1000)

    else if ( harSw .and. extSw .and. .not.xcSw ) then ! Coulomb potential
      call fopen( 1000, name='vpw_coul', status='old', action='read', form='unformatted' )
        read( 1000 ) vpwStar
      call fclose( 1000 )

      call fopen(1000, name='v0MTFLEUR_coul', status='old', action='read', form='unformatted')
        read(1000) vrFLEUR
      call fclose(1000)

    else if ( harSw .and. extSw .and. xcSw ) then ! effective potential
      call fopen( 1000, name='vpw_eff', status='old', action='read', form='unformatted' )
        read( 1000 ) vpwStar
      call fclose( 1000 )
      call fopen(1000, name='v0MTFLEUR_eff', status='old', action='read', form='unformatted')
        read(1000) vrFLEUR
      call fclose(1000)

    else
      call juDFT_warn('FLEUR benchmark for this potential combination not available. Setting benchmark curve equals zero.', &
        & calledby='checkGradPot0s', hint='Choose Hartree, external, exchange-correlation, Coulomb or effective potential', &
        & file='jpTestPotential_mod.f90')
      fBenchAvailSw = .false.

    end if

    call ConvertStar2G( vpwStar(:, 1), vpw_G, starsT, ngdp, gdp )

    do iGvec = 1, ngdp
      Gext(:) = matmul(cellT%bmat, gdp(:, iGvec))
      vCoulBenchmark(iGvec, :) = iu * Gext(:) * vpw_G(iGvec)
    enddo
    call calcGrFinLH(atomsT, sphharT, clnu_atom, nmem_atom, mlh_atom, vrFleur(:, :, :, 1), gradVrFleur)


    !call fopen(1111, name='pathPseudoPot', status='replace', action='write', form='formatted')
    !call fopen(1222, name='potFLEUR', status='replace', action='write', form='formatted')
    call fopen(1333, name='pathPot', status='replace', action='write', form='formatted')
    call fopen(1444, name='pathPotCont', status='replace', action='write', form='formatted')

    dx = 1. / 500.
    x = 0 - dx
    wasIR = -1
    do ii = 0, 500
      x = x + dx
      ucpath = [paPoXCo * x, paPoYCo * x, paPoZCo * x] !CHANGE
      iatom = 1
      isIR = 0
      do itype = 1, atomsT%ntype
        do ieq = 1, atomsT%neq(itype)
          do c1 = -1, 1
            do c2 = -1, 1
              do c3 = -1, 1
                ext2brav = [c1, c2, c3]
                ucpathExt= matmul(cellT%amat, ucpath - atomsT%taual(:, iatom) - ext2brav)
                if ( norm2 ( ucpathExt ) < atomsT%rmt(itype) ) then
                  isIR = iatom ! todo atomLock is redundant!
                  atomLock = iatom ! maybe also possible with two exit statements if MT found then lock variables are unneccssary
                  ext2bravLock = ext2brav
                  itypeLock = itype
                end if
              end do
            end do
          end do
          iatom = iatom + 1
        end do
      end do

      if ( pseudoPot ) then

        vpwRfleur = cmplx(0.,0.)
        vpwRjuPhon = cmplx(0.,0.)

        do iGvec = 1, ngdp
          exponential = exp(iu * tpi * dot_product(gdp(:, iGvec), ucpath))
          vpwRfleur(:)  = vpwRfleur(:)  + vCoulBenchmark(iGvec, :) * exponential
          vpwRjuPhon(:) = vpwRjuPhon(:) + vGradCoul0IR(iGvec, :)   * exponential
        end do

      !  write (1111, '(13(es16.8E3, 2x),i2, 3(2x, i2))') x, real(vpwRfleur(1)), aimag(vpwRfleur(1)), real(vpwRfleur(2)), &
      !    & aimag(vpwRfleur(2)), real(vpwRfleur(3)), aimag(vpwRfleur(3)), real(vpwRjuPhon(1)), aimag(vpwRjuPhon(1)),     &
      !    & real(vpwRjuPhon(2)), aimag(vpwRjuPhon(2)), real(vpwRjuPhon(3)), aimag(vpwRjuPhon(3)), isIR, ext2bravLock(1), &
      !    & ext2bravLock(2), ext2bravLock(3)

      end if

      if ( realPot ) then ! if no pseudoPotMode
        if ( wasIR == -1 ) then ! initial point of path
          if ( isIR == 0) then ! is in interstitial, so point can be simply plotted

            vpwRfleur = cmplx(0.,0.)
            vpwRjuPhon = cmplx(0.,0.)
            onlyPotFLEUR = cmplx(0.,0.)

            do iGvec = 1, ngdp ! calculate Interstitial Potential
              exponential = exp(iu * tpi * dot_product(gdp(:, iGvec), ucpath))
              vpwRfleur(:) = vpwRfleur(:) + vCoulBenchmark(iGvec, :) * exponential
              vpwRjuPhon(:) = vpwRjuPhon(:) + vGradCoul0IR(iGvec, :) * exponential ! todo has been changed, seems to be added so that starting point ist correct!
              !todo floating is inexact here for some reasons....
              onlyPotFLEUR = onlyPotFLEUR +1!+ !vpw_G(iGvec) * exponentialFLEUR ! pure potential
            end do

            !write (1222, '(13(es16.8E3, 2x))') x, real(onlyPotFLEUR), aimag(onlyPotFLEUR)
            write (1333, '(13(es16.8E3, 2x),i2, 3(2x, i2))') & ! write out calculated point
              &x, real(vpwRfleur(1)), aimag(vpwRfleur(1)), real(vpwRfleur(2)), aimag(vpwRfleur(2)),  real(vpwRfleur(3)), &
              &aimag(vpwRfleur(3)), real(vpwRjuPhon(1)), aimag(vpwRjuPhon(1)), real(vpwRjuPhon(2)),  aimag(vpwRjuPhon(2)), &
              &real(vpwRjuPhon(3)), aimag(vpwRjuPhon(3)), isIR, ext2bravLock(1), ext2bravLock(2), ext2bravLock(3)

          end if ! if not in interstitial we wait till reaching the MT surface enabling to solve the Dirichelet boundary problem
        else ! wasIR is available
          if ( ( wasIR == 0 .and. isIR == 0 )) then ! is between atoms

            vpwRfleur = cmplx(0.,0.)
            vpwRjuPhon = cmplx(0.,0.)
            onlyPotFLEUR = cmplx(0.,0.)

            do iGvec = 1, ngdp ! calculate Interstitial Potential
              exponential = exp(iu * tpi * dot_product(gdp(:, iGvec), ucpath))
              vpwRfleur(:) = vpwRfleur(:) + vCoulBenchmark(iGvec, :) * exponential
              vpwRjuPhon(:) = vpwRjuPhon(:) + vGradCoul0IR(iGvec, :) * exponential
              onlyPotFLEUR = onlyPotFLEUR + vpw_G(iGvec) * exponential ! pure potential
            end do

            !write (1222, '(13(es16.8E3, 2x))') x, real(onlyPotFLEUR), aimag(onlyPotFLEUR)
            write (1333, '(13(es16.8E3, 2x),i2, 3(2x, i2))') & ! write out calculated point
              &x, real(vpwRfleur(1)), aimag(vpwRfleur(1)), real(vpwRfleur(2)), aimag(vpwRfleur(2)),  real(vpwRfleur(3)), &
              &aimag(vpwRfleur(3)), real(vpwRjuPhon(1)), aimag(vpwRjuPhon(1)), real(vpwRjuPhon(2)),  aimag(vpwRjuPhon(2)), &
              &real(vpwRjuPhon(3)), aimag(vpwRjuPhon(3)), isIR, ext2bravLock(1), ext2bravLock(2), ext2bravLock(3)

          else if (wasIR == 0 .and. isIR >= 0 ) then ! crosses MT surface into the atom

            ucpathA = ucpath ! ucpathA is in atom
            ucpathB = ucpath - [paPoXCo * dx, paPoYCo * dx, paPoZCo * dx] !todo make these dx vector general! ucpathB is outside atom CHANGE
            do iterations = 1, 1000
              ucpathC = (ucpathA + ucpathB) / 2
              ucpathExtc = matmul(cellT%amat, ucpathC - atomsT%taual(:, atomLock) - ext2bravLock)
              if ( ( abs( norm2 ( ucpathExtc ) - atomsT%rmt(itypeLock) ) < 1e-7 ) .or. norm2( ucpathB - ucpathA ) < 1e-7 ) then
              !  ! mtsphere found
                exit
              else
                if ( ( norm2(ucpathExtc) - atomsT%rmt( itypeLock ) ) < 0) then
                  ucpathA = ucpathC
                else
                  ucpathB = ucpathC
                end if
              end if
            end do
            if (iterations > 1000) then
              write(*, *) 'Warning: Iteration not converged' !todo juDFT warning
            end if

            vpwRfleur = cmplx(0.,0.)
            vpwRjuPhon = cmplx(0.,0.)
            onlyPotFLEUR = cmplx(0., 0.)

            do iGvec = 1, ngdp ! calculate Interstitial Potential
              exponential = exp(iu * tpi * dot_product(gdp(:, iGvec), ucpathC))
              vpwRfleur(:) = vpwRfleur(:) + vCoulBenchmark(iGvec, :) * exponential
              vpwRjuPhon(:) = vpwRjuPhon(:) + vGradCoul0IR(iGvec, :) * exponential
              !onlyPotFLEUR = onlyPotFLEUR + vpw_G(iGvec) * exponential ! pure potential
            end do
            !write (1222, '(13(es16.8E3, 2x))') x, real(onlyPotFLEUR), aimag(onlyPotFLEUR)
            write (1333, '(13(es16.8E3, 2x),i2, 3(2x, i2))') & ! write out calculated point
              &x, real(vpwRfleur(1)), aimag(vpwRfleur(1)), real(vpwRfleur(2)), aimag(vpwRfleur(2)),  real(vpwRfleur(3)), &
              &aimag(vpwRfleur(3)), real(vpwRjuPhon(1)), aimag(vpwRjuPhon(1)), real(vpwRjuPhon(2)),  aimag(vpwRjuPhon(2)), &
              &real(vpwRjuPhon(3)), aimag(vpwRjuPhon(3)), wasIR, ext2bravLock(1), ext2bravLock(2), ext2bravLock(3)

            direc = ucpathC - atomsT%taual(:, atomLock) - ext2bravLock
            direcExt = matmul(cellT%amat, direc)
            direcExt = direcExt / norm2(direcExt)
!            call ylmnorm_init(atomsT%lmaxd + 1)
            !call ylm4(atomsT%lmaxd + 1, direcExt, ylm) ! todo really until l + 1
            call ylm4(atomsT%lmaxd, direcExt, ylm) ! todo really until l + 1
!            call ylmnorm_init(atomsT%lmaxd)

            vMTfleur = 0
            vMTjuPhon = 0
            do idirec = 1, 3
              do oqn_l = 0, atomsT%lmax(itypeLock)
                do mqn_m = -oqn_l, oqn_l
                  lm = oqn_l * (oqn_l + 1) + 1 + mqn_m
                  do imesh = 1, atomsT%jri(itypeLock)
                    vMTfLEUR(idirec, imesh) = vMTFLEUR(idirec, imesh) + gradVrFleur(imesh, lm, atomLock, idirec) * ylm(lm)
                    vMTjuPhon(idirec, imesh) = vMTjuPhon(idirec, imesh) + gradVrjuPhon(imesh, lm, idirec, atomLock) * ylm(lm)
                  end do
                end do
              end do
            end do

            do imesh = atomsT%jri(itypeLock), 1, -1
              write(1333, '(13(es20.10E3, 2x), i2, 3(2x, i2))')&
                &-atomsT%rmsh(imesh, itypeLock), real(vMTfleur(1, imesh)), aimag(vMTfleur(1, imesh)), &
                &real(vMTfleur(2, imesh)), aimag(vMTfleur(2, imesh)), real(vMTfleur(3, imesh)), &
                &aimag(vMTfleur(3, imesh)), real(vMTjuPhon(1, imesh)) ,aimag(vMTjuPhon(1, imesh)),&
                &real(vMTjuPhon(2, imesh)), aimag(vMTjuPhon(2, imesh)), real(vMTjuPhon(3, imesh)), &
                aimag(vMTjuPhon(3, imesh)), atomLock, ext2bravLock(1), ext2bravLock(2), ext2bravLock(3)
            end do

            do idirec = 1, 3
              if ( ( vpwRfleur(idirec) /= 0 ) .and. ( vpwRjuPhon(idirec) /= 0 ) ) then
                write(1444, '(8(es16.8E3, 2x))') &
                  &real(vMTfleur(idirec, atomsT%jri(itypeLock))), real(vpwRfleur(idirec)), &
                  &abs(real(vMTfleur(idirec, atomsT%jri(itypeLock))) - real(vpwRfleur(idirec))), &
                  &abs(real(vMTfleur(idirec, atomsT%jri(itypeLock))) - real(vpwRfleur(idirec))) / abs(real(vpwRfleur(idirec))),&
                  &real(vMTjuPhon(idirec, atomsT%jri(itypeLock))), real(vpwRjuPhon(idirec)), &
                  &abs(real(vMTjuPhon(idirec, atomsT%jri(itypeLock))) - real(vpwRjuPhon(idirec))), &
                  &abs(real(vMTjuPhon(idirec, atomsT%jri(itypeLock))) - real(vpwRjuPhon(idirec))) / abs(real(vpwRjuPhon(idirec)))
              else if ( ( vpwRfleur(idirec) == 0 ) .and. ( vpwRjuPhon(idirec) /= 0 ) ) then
                write(1444, '(8(es16.8E3, 2x))') &
                  &real(vMTfleur(idirec, atomsT%jri(itypeLock))), real(vpwRfleur(idirec)), &
                  &abs(real(vMTfleur(idirec, atomsT%jri(itypeLock))) - real(vpwRfleur(idirec))), -9e9,&
                  &real(vMTjuPhon(idirec, atomsT%jri(itypeLock))), real(vpwRjuPhon(idirec)), &
                  &abs(real(vMTjuPhon(idirec, atomsT%jri(itypeLock))) - real(vpwRjuPhon(idirec))), &
                  &abs(real(vMTjuPhon(idirec, atomsT%jri(itypeLock))) - real(vpwRjuPhon(idirec))) / abs(real(vpwRjuPhon(idirec)))
              else if ( ( vpwRfleur(idirec) /= 0 ) .and. ( vpwRjuPhon(idirec) == 0 ) ) then
                write(1444, '(8(es16.8E3, 2x))') &
                  &real(vMTfleur(idirec, atomsT%jri(itypeLock))), real(vpwRfleur(idirec)), &
                  &abs(real(vMTfleur(idirec, atomsT%jri(itypeLock))) - real(vpwRfleur(idirec))), &
                  &abs(real(vMTfleur(idirec, atomsT%jri(itypeLock))) - real(vpwRfleur(idirec))) / abs(real(vpwRfleur(idirec))),&
                  &real(vMTjuPhon(idirec, atomsT%jri(itypeLock))), real(vpwRjuPhon(idirec)), &
                  &abs(real(vMTjuPhon(idirec, atomsT%jri(itypeLock))) - real(vpwRjuPhon(idirec))), -9e9
              else
                write(1444, '(8(es16.8E3, 2x))') &
                  &real(vMTfleur(idirec, atomsT%jri(itypeLock))), real(vpwRfleur(idirec)), &
                  &abs(real(vMTfleur(idirec, atomsT%jri(itypeLock))) - real(vpwRfleur(idirec))), -9e9,&
                  &real(vMTjuPhon(idirec, atomsT%jri(itypeLock))), real(vpwRjuPhon(idirec)), &
                  &abs(real(vMTjuPhon(idirec, atomsT%jri(itypeLock))) - real(vpwRjuPhon(idirec))), -9e9
              end if
            end do

          else if (wasIR >= 0 .and. isIR == 0 ) then ! crosses MT surface out from atom

            ucpathA = ucpath ! is out of atom
            ucpathB = ucpath - [paPoXCo * dx, paPoYCo * dx, paPoZCo * dx] ! is in atom, todo make these dx vector general!
            do iterations = 1, 1000 ! really so much iterations needed?
              ucpathC = (ucpathA + ucpathB) / 2
              ucpathExtc = matmul(cellT%amat, ucpathC - atomsT%taual(:, atomLock) - ext2bravLock)
              if ( ( abs( norm2 ( ucpathExtc ) - atomsT%rmt(itypeLock) ) <= 1e-7 ) .or. norm2( ucpathB - ucpathA ) < 1e-7 ) then
                ! mtsphere found
                exit
              else
                if ( ( norm2( ucpathExtc ) - atomsT%rmt( itypeLock ) ) > 0) then
                  ucpathA = ucpathC
                else
                  ucpathB = ucpathC
                end if
              end if
            end do
            if (iterations >= 1000) then
              write(*, *) 'Warning: Iteration not converged'
            end if

            direc = ucpathC - atomsT%taual(:, atomLock) - ext2bravLock
            direcExt = matmul(cellT%amat, direc)
            direcExt = direcExt / norm2(direcExt)
  !          call ylmnorm_init(atomsT%lmaxd + 1)
            !call ylm4(atomsT%lmaxd + 1, direcExt, ylm) ! todo really until l + 1
            call ylm4(atomsT%lmaxd, direcExt, ylm) ! todo really until l + 1
  !          call ylmnorm_init(atomsT%lmaxd)

            vMTfleur = 0
            vMTjuPhon = 0
            do idirec = 1, 3
              do oqn_l = 0, atomsT%lmax(itypeLock)
                do mqn_m = -oqn_l, oqn_l
                  lm = oqn_l * (oqn_l + 1) + 1 + mqn_m
                  do imesh = 1, atomsT%jri(itypeLock)
                    vMTfLEUR(idirec, imesh) = vMTFLEUR(idirec, imesh) + gradVrFleur(imesh, lm, atomLock, idirec) * ylm(lm)
                    vMTjuPhon(idirec, imesh) = vMTjuPhon(idirec, imesh) + gradVrjuPhon(imesh, lm, idirec, atomLock) * ylm(lm)
                  end do
                end do
              end do
            end do
            do imesh = 1, atomsT%jri(itypeLock)
              write(1333, '(13(es20.10E3, 2x),i2, 3(2x, i2))')&
                &atomsT%rmsh(imesh, itypeLock), real(vMTfleur(1, imesh)), aimag(vMTfleur(1, imesh)), &
                &real(vMTfleur(2, imesh)), aimag(vMTfleur(2, imesh)), real(vMTfleur(3, imesh)), &
                &aimag(vMTfleur(3, imesh)), real(vMTjuPhon(1, imesh)) ,aimag(vMTjuPhon(1, imesh)),&
                &real(vMTjuPhon(2, imesh)), aimag(vMTjuPhon(2, imesh)), real(vMTjuPhon(3, imesh)), &
                &aimag(vMTjuPhon(3, imesh)), atomLock, ext2bravLock(1), ext2bravLock(2), ext2bravLock(3)
            end do

            !todo seems to be the real ptoential here
            !iatom = 0
            potFleurMT = 0
              !Evaluate MT sphere function at same random points
            !symAt_temp  = atomsT%ntypsy(atomLock) ! symmetry of atom iatom
            !symOpr_temp = atomsT%ngopr (atomLock) ! symmetry operation mapping local to global coordinate system
            !do imesh = 1, atomsT%jri(itypeLock)
            !  rvec = direcExt
            !  if (symOpr_temp /= 1) then
            !  !Transform into internal real space coordinates
            !  !todo might be wrong bmat
            !    call cotra1(rvec, rvec_int, cellT%bmat) !allocate and deallocate randPtsGrid
            !    do imatmulr = 1, 3
            !      rotatedVector(imatmulr) = 0.0
            !      do imatmulc = 1, 3
            !        rotatedVector(imatmulr) = rotatedVector(imatmulr) + symT%mrot(imatmulr, imatmulc, symOpr_temp) * rvec_int(imatmulc)
            !      end do
            !    end do
            !    call cotra0(rotatedVector, rvec, cellT%amat)
            !  end if
            !  call ylm4(atomsT%lmax(itypeLock), rvec, ylm)

            !  do ilath = 0, sphharT%nlh(symAt_temp)
            !    linCombBas = 0.0
            !    l_temp = sphharT%llh(ilath, symAt_temp) * (sphharT%llh(ilath, symAt_temp) + 1) + 1 !understand this
            !    do imem = 1, sphharT%nmem(ilath, symAt_temp)
            !      lm_temp = l_temp + sphharT%mlh(imem, ilath, symAt_temp)
            !      linCombBas = linCombBas + real(clnu_atom(imem, ilath, symAt_temp) * ylm(lm_temp))
            !    end do
            !      potFleurMT(imesh) = potFleurMT(imesh) + vrFleur(imesh, ilath, itypeLock, 1) * linCombBas
            !  end do
            !end do
            !call fopen(1020, name='potMTfleur', status='replace', action='write', form='formatted')
            !do imesh = 1, atomsT%jri(itypeLock)
            !  write(1020, '(3(es16.8E3, 2x))') atomsT%rmsh(imesh, itypeLock), potFleurMT(imesh), real(vMTfLEUR(1, imesh))
            !end do
            !call fclose(1020)


            vpwRfleur = cmplx(0.,0.)
            vpwRjuPhon = cmplx(0.,0.)
            onlyPotFLEUR = cmplx(0., 0.)
            do iGvec = 1, ngdp ! calculate Interstitial Potential
              exponential = exp(iu * tpi * dot_product(gdp(:, iGvec), ucpathC))
              vpwRfleur(:) = vpwRfleur(:) + vCoulBenchmark(iGvec, :) * exponential
              vpwRjuPhon(:) = vpwRjuPhon(:) + vGradCoul0IR(iGvec, :) * exponential
              onlyPotFLEUR = onlyPotFLEUR + vpw_G(iGvec) * exponential ! pure potential
            end do
            !write (1222, '(13(es16.8E3, 2x))') x, real(onlyPotFLEUR), aimag(onlyPotFLEUR) ! todo the lines down really needed?
            write (1333, '(13(es16.8E3, 2x),i2, 3(2x, i2))') & ! write out calculated point
              &x, real(vpwRfleur(1)), aimag(vpwRfleur(1)), real(vpwRfleur(2)), aimag(vpwRfleur(2)),  real(vpwRfleur(3)), &
              &aimag(vpwRfleur(3)), real(vpwRjuPhon(1)), aimag(vpwRjuPhon(1)), real(vpwRjuPhon(2)),  aimag(vpwRjuPhon(2)), &
              &real(vpwRjuPhon(3)), aimag(vpwRjuPhon(3)), isIR, ext2bravLock(1), ext2bravLock(2), ext2bravLock(3)

            do idirec = 1, 3
              if ( ( vpwRfleur(idirec) /= 0 ) .and. ( vpwRjuPhon(idirec) /= 0 ) ) then
                write(1444, '(8(es16.8E3, 2x))') &
                  &real(vMTfleur(idirec, atomsT%jri(itypeLock))), real(vpwRfleur(idirec)), &
                  &abs(real(vMTfleur(idirec, atomsT%jri(itypeLock))) - real(vpwRfleur(idirec))), &
                  &abs(real(vMTfleur(idirec, atomsT%jri(itypeLock))) - real(vpwRfleur(idirec))) / abs(real(vpwRfleur(idirec))),&
                  &real(vMTjuPhon(idirec, atomsT%jri(itypeLock))), real(vpwRjuPhon(idirec)), &
                  &abs(real(vMTjuPhon(idirec, atomsT%jri(itypeLock))) - real(vpwRjuPhon(idirec))), &
                  &abs(real(vMTjuPhon(idirec, atomsT%jri(itypeLock))) - real(vpwRjuPhon(idirec))) / abs(real(vpwRjuPhon(idirec)))
              else if ( ( vpwRfleur(idirec) == 0 ) .and. ( vpwRjuPhon(idirec) /= 0 ) ) then
                write(1444, '(8(es16.8E3, 2x))') &
                  &real(vMTfleur(idirec, atomsT%jri(itypeLock))), real(vpwRfleur(idirec)), &
                  &abs(real(vMTfleur(idirec, atomsT%jri(itypeLock))) - real(vpwRfleur(idirec))), -9e9,&
                  &real(vMTjuPhon(idirec, atomsT%jri(itypeLock))), real(vpwRjuPhon(idirec)), &
                  &abs(real(vMTjuPhon(idirec, atomsT%jri(itypeLock))) - real(vpwRjuPhon(idirec))), &
                  &abs(real(vMTjuPhon(idirec, atomsT%jri(itypeLock))) - real(vpwRjuPhon(idirec))) / abs(real(vpwRjuPhon(idirec)))
              else if ( ( vpwRfleur(idirec) /= 0 ) .and. ( vpwRjuPhon(idirec) == 0 ) ) then
                write(1444, '(8(es16.8E3, 2x))') &
                  &real(vMTfleur(idirec, atomsT%jri(itypeLock))), real(vpwRfleur(idirec)), &
                  &abs(real(vMTfleur(idirec, atomsT%jri(itypeLock))) - real(vpwRfleur(idirec))), &
                  &abs(real(vMTfleur(idirec, atomsT%jri(itypeLock))) - real(vpwRfleur(idirec))) / abs(real(vpwRfleur(idirec))),&
                  &real(vMTjuPhon(idirec, atomsT%jri(itypeLock))), real(vpwRjuPhon(idirec)), &
                  &abs(real(vMTjuPhon(idirec, atomsT%jri(itypeLock))) - real(vpwRjuPhon(idirec))), -9e9
              else
                write(1444, '(8(es16.8E3, 2x))') &
                  &real(vMTfleur(idirec, atomsT%jri(itypeLock))), real(vpwRfleur(idirec)), &
                  &abs(real(vMTfleur(idirec, atomsT%jri(itypeLock))) - real(vpwRfleur(idirec))), -9e9,&
                  &real(vMTjuPhon(idirec, atomsT%jri(itypeLock))), real(vpwRjuPhon(idirec)), &
                  &abs(real(vMTjuPhon(idirec, atomsT%jri(itypeLock))) - real(vpwRjuPhon(idirec))), -9e9
              end if
            end do

          end if! if still in MT then wait for crossing

        end if !is not first point on path
      end if ! is not in pseudopot mode
      wasIR = isIR
    end do
    !call fclose(1111)
    !call fclose(1222)
    call fclose(1333)
    call fclose(1444)
  end subroutine checkGradPot0s

  ! Writes plot files to compare grVeff which has been generated by the Weinert method and by the gradient
  subroutine checkGradPot0sXSF(symT, cellT, inputT, sphharT, starsT, atomsT, vGradCoul0IR, ngdp, gdp, mlh_atom, nmem_atom, clnu_atom, rho0MT, gradVrjuPhon, paPoXCo, paPoYCo, paPoZCo, harSw, extSw, xcSw, dimensT)

    use m_types
    use m_jPConstants, only : iu, pi, tpi, fpi, Tmatrix
    use m_gaunt
    use m_ylm_old
    use m_cotra
    use m_jpPotDensHelper
    use m_juDFT_NOstopNO
    use m_jpPlotObservables, only : plotVestPotDens1

    implicit none

    type(t_sym),        intent(in) :: symT
    type(t_dimension),        intent(in) :: dimensT
    type(t_cell),       intent(in) :: cellT
    type(t_sphhar),     intent(in) :: sphharT
    type(t_stars),      intent(in) :: starsT
    type(t_atoms),      intent(in) :: atomsT
    type(t_input),      intent(in) :: inputT
    complex,            intent(in) :: vGradCoul0IR(:, :)
    integer,            intent(in) :: ngdp
    integer,            intent(in) :: gdp(:, :)
    integer,            intent(in) :: mlh_atom(:, 0:, :)
    integer,            intent(in) :: nmem_atom(0:, :)
    complex,            intent(in) :: clnu_atom(:, 0:, :)
    real,               intent(in) :: rho0MT(:, 0:, :, :)
    real,               intent(in) :: paPoXCo
    real,               intent(in) :: paPoYCo
    real,               intent(in) :: paPoZCo
    logical,            intent(in) :: harSw
    logical,            intent(in) :: extSw
    logical,            intent(in) :: xcSw

    integer                        :: iGvec
    real                           :: Gext(3)
    complex                        :: vCoulBenchmark(ngdp, 3)
    complex                        :: gradVcMT(atomsT%jmtd, (atomsT%lmaxd + 1)**2, atomsT%nat, 3)
    real                           :: sqr4pi3
    real                           :: Vc0nSym(atomsT%jmtd, (atomsT%lmaxd + 1 +  1)**2, atomsT%nat) ! todo is this the correct dimension
    integer                        :: iatom
    integer                        :: itype
    integer                        :: ineq
    integer                        :: symType
    integer                        :: lh
    integer                        :: oqn_l
    integer                        :: lm_temp
    integer                        :: mem
    integer                        :: mems
    integer                        :: mqn_m
    integer                        :: mqn_mpp
    integer                        :: lm
    integer                        :: lmMinus
    integer                        :: lmPlus
    real                           :: Vc0MTDerived(atomsT%jmtd)
    real                           :: tempGaunt1, tempGaunt2
    integer                        :: imesh
    complex                        :: vpw_G(ngdp)
    complex                        :: vpw_G_coul(ngdp)
    complex                        :: vpw_G_hart(ngdp)
    complex                        :: vpw_G_xc(ngdp)
    complex              :: vpwCoulUwStars1(starsT%n3d, 1)
    complex              :: vpwCoulUwStars2(starsT%n3d, 1)
    complex              :: vpwStar(starsT%n3d, 1)
    complex              :: vpwStarDiff(starsT%n3d, 1)
    complex                        :: vpwRfleur(3)
    complex                        :: vpwRjuPhon(3)
    complex                        :: exponentialjuPhon
    complex                        :: exponentialFLEUR
    complex                        :: exponential
    integer           :: isIR
    integer           :: wasIR
    integer          :: c1, c2, c3
    real                           :: ucpath(3)
    real                              :: ucpathA(3), ucpathB(3), ucpathC(3)
    real                           :: ucpathExt(3)
    real                      :: ucpathExtc(3), ucpathExta(3)
    real                           :: dx
    integer                        :: ii
    real                           :: x
    complex                        :: newvpwFLEUR(starsT%n3d, 1) ! make this declaration more general!
    complex                        :: newvpw_G(ngdp)
    complex                        :: onlyPotFLEUR
    complex                        :: onlypsqFLEUR
    integer                        :: ieq
    integer                        :: ext2brav(3)
    integer                        :: ext2bravLock(3)
    integer                        :: atomLock
    integer                        :: itypeLock
    logical                        :: pseudoPot = .true.
    logical                        :: realPot = .true.
    integer                   :: iterations
  real   :: vrFleur(atomsT%jmtd, 0:sphharT%nlhd, atomsT%ntypd, inputT%jspins)
  real   :: r2vrFleur(atomsT%jmtd, 0:sphharT%nlhd, atomsT%ntypd, inputT%jspins)
  real   :: vrFleurDiff(atomsT%jmtd, 0:sphharT%nlhd, atomsT%ntypd, inputT%jspins)
  real   :: vrFleurxc(atomsT%jmtd, 0:sphharT%nlhd, atomsT%ntypd, inputT%jspins)
  real   :: vrTemp(atomsT%jmtd, 0:sphharT%nlhd, atomsT%ntypd, inputT%jspins)
  complex,             intent(in) :: gradVrjuPhon(:, :, :, :)
  complex,      allocatable :: gradVrFleur(:, :, :, :)
  complex,      allocatable :: r2GrVrFleur(:, :, :, :)
  real    :: direc(3)
  real :: direcExt(3)
  complex                                     :: ylm((atomsT%lmaxd + 1)**2)
  complex                                     :: vMTfleur(3, atomsT%jmtd)
  complex                                     :: vMTjuPhon(3, atomsT%jmtd)
  integer                          :: idirec
  real                              :: av, dmx, rms
  real                              :: IRvalue(1), MTvalue(1)
  integer                             :: jspin
  integer                               :: ilh
  real        :: potFleurMT(atomsT%jmtd)
  integer    :: symAt_temp, symOpr_temp
  real      :: rvec(3), rvec_int(3)
  integer    :: imatmulr, imatmulc
  real :: rotatedVector(3)
  integer :: ilath
  real :: linCombBas
  integer :: l_temp
  integer :: imem
  integer :: idir
  integer :: ieqat
  integer :: lm_pre
    complex,        allocatable :: r2AbsError(:, :, :, :)

logical :: fBenchAvailSw
      character(len=30)                          :: plotname

    allocate( r2AbsError(atomsT%jmtd, ( atomsT%lmaxd + 1)**2, atomsT%nat, 3))
     r2AbsError(:, :, :, :) = cmplx(0., 0.)
    ! in this test the gradient of the Coulomb potential in the interstitial region calculated with Aaron's method must be equal to iG qpw
    fBenchAvailSw = .true.
    vpwStar = 0
    vpwStarDiff = 0
    vrFLEUR = 0
    vrFLEURDiff = 0

    !todo do this with pottot and potcoul
    if ( harSw .and. .not.extSw .and. .not.xcSw ) then ! Hartree potential
      call fopen( 1000, name='vpw_hart', status='old', action='read', form='unformatted' )
        read( 1000 ) vpwStar
      call fclose( 1000 )

      call fopen(1000, name='v0MTFLEUR_hart', status='old', action='read', form='unformatted')
        read(1000) vrFLEUR
      call fclose(1000)

    else if ( .not.harSw .and. extSw .and. .not.xcSw ) then ! external potential
      call fopen( 1000, name='vpw_coul', status='old', action='read', form='unformatted' )
        read( 1000 ) vpwStar
      call fclose( 1000 )

      call fopen( 1000, name='vpw_hart', status='old', action='read', form='unformatted' )
        read( 1000 ) vpwStarDiff
      call fclose( 1000 )

      vpwStar = vpwStar - vpwStarDiff


      call fopen(1000, name='v0MTFLEUR_coul', status='old', action='read', form='unformatted')
        read(1000) vrFLEUR
      call fclose(1000)

      call fopen(1000, name='v0MTFLEUR_hart', status='old', action='read', form='unformatted')
        read(1000) vrFLEURDiff
      call fclose(1000)

      vrFLEUR = vrFLEUR - vrFLEURDiff

    else if ( .not.harSw .and. .not.extSw .and. xcSw ) then ! exchange-correlation potential
      call fopen( 1000, name='vpw_xc', status='old', action='read', form='unformatted' )
        read( 1000 ) vpwStar
      call fclose( 1000 )

      call fopen(1000, name='v0MTFLEUR_xc', status='old', action='read', form='unformatted')
        read(1000) vrFLEUR
      call fclose(1000)

    else if ( harSw .and. extSw .and. .not.xcSw ) then ! Coulomb potential
      call fopen( 1000, name='vpw_coul', status='old', action='read', form='unformatted' )
        read( 1000 ) vpwStar
      call fclose( 1000 )

      call fopen(1000, name='v0MTFLEUR_coul', status='old', action='read', form='unformatted')
        read(1000) vrFLEUR
      call fclose(1000)

    else if ( harSw .and. extSw .and. xcSw ) then ! effective potential
      call fopen( 1000, name='vpw_eff', status='old', action='read', form='unformatted' )
        read( 1000 ) vpwStar
      call fclose( 1000 )
      call fopen(1000, name='v0MTFLEUR_eff', status='old', action='read', form='unformatted')
        read(1000) vrFLEUR
      call fclose(1000)

    else
      call juDFT_warn('FLEUR benchmark for this potential combination not available. Setting benchmark curve equals zero.', &
        & calledby='checkGradPot0s', hint='Choose Hartree, external, exchange-correlation, Coulomb or effective potential', &
        & file='jpTestPotential_mod.f90')
      fBenchAvailSw = .false.

    end if

    call ConvertStar2G( vpwStar(:, 1), vpw_G, starsT, ngdp, gdp )

    do iGvec = 1, ngdp
      Gext(:) = matmul(cellT%bmat, gdp(:, iGvec))
      vCoulBenchmark(iGvec, :) = iu * Gext(:) * vpw_G(iGvec)
    enddo

    do itype = 1, atomsT%ntype
      do ilh = 0, sphharT%nlhd
        do imesh = 1, atomsT%jri(itype)
          r2VrFleur(imesh, ilh, itype, 1) = vrFLEUR(imesh, ilh, itype, 1) * atomsT%rmsh(imesh, itype) * atomsT%rmsh(imesh, itype)
        end do
      end do
    end do

    call calcGrR2FinLH( atomsT, sphharT, clnu_atom, nmem_atom, mlh_atom, r2VrFleur(:, :, :, 1), r2GrVrFleur )

    do idir = 1, 3
      iatom = 0
      do itype = 1, atomsT%ntype
        do ieqat = 1, atomsT%neq(itype)
          iatom = iatom + 1
          do oqn_l = 0, atomsT%lmax(itype)
            lm_pre = oqn_l * (oqn_l + 1) + 1
            do mqn_m = -oqn_l, oqn_l
              lm = lm_pre + mqn_m
              do imesh = 1, atomsT%jri(itype)
                r2AbsError(imesh, lm, iatom, idir) = abs(gradVrjuPhon(imesh, lm, idir, iatom) * atomsT%rmsh(imesh, itype)**2 - r2GrVrFleur(imesh, lm, iatom, idir))
              end do ! imesh
            end do ! mqn_m
          end do ! oqn_l
        end do ! ieqat
      end do ! itype
    end do ! idir



    write(*, *) 'Plot comparison'
    do idir = 1, 3
      plotname = ''
      write(plotname, '(a,i1,a,i1,a)') 'plotCompMTGradNWeinert', idir, 'at',1, '.xsf'
      ! Create file of density variation for Vesta
      call plotVestPotDens1(dimensT, starsT, sphharT, atomsT, inputT, symT, cellT, gdp, ngdp, [0., 0., 0.], &
                              & vCoulBenchmark(:, idir) - vCoulBenchmark(:, idir), r2AbsError(:, :, :, idir), plotname, 1, .FALSE., "plot_inpMTGradComp")
    end do ! idir

!    !call fopen(1111, name='pathPseudoPot', status='replace', action='write', form='formatted')
!    !call fopen(1222, name='potFLEUR', status='replace', action='write', form='formatted')
!    call fopen(1333, name='pathPot', status='replace', action='write', form='formatted')
!    call fopen(1444, name='pathPotCont', status='replace', action='write', form='formatted')
!
!    dx = 1. / 500.
!    x = 0 - dx
!    wasIR = -1
!    do ii = 0, 500
!      x = x + dx
!      ucpath = [paPoXCo * x, paPoYCo * x, paPoZCo * x] !CHANGE
!      iatom = 1
!      isIR = 0
!      do itype = 1, atomsT%ntype
!        do ieq = 1, atomsT%neq(itype)
!          do c1 = -1, 1
!            do c2 = -1, 1
!              do c3 = -1, 1
!                ext2brav = [c1, c2, c3]
!                ucpathExt= matmul(cellT%amat, ucpath - atomsT%taual(:, iatom) - ext2brav)
!                if ( norm2 ( ucpathExt ) < atomsT%rmt(itype) ) then
!                  isIR = iatom ! todo atomLock is redundant!
!                  atomLock = iatom ! maybe also possible with two exit statements if MT found then lock variables are unneccssary
!                  ext2bravLock = ext2brav
!                  itypeLock = itype
!                end if
!              end do
!            end do
!          end do
!          iatom = iatom + 1
!        end do
!      end do
!
!      if ( pseudoPot ) then
!
!        vpwRfleur = cmplx(0.,0.)
!        vpwRjuPhon = cmplx(0.,0.)
!
!        do iGvec = 1, ngdp
!          exponential = exp(iu * tpi * dot_product(gdp(:, iGvec), ucpath))
!          vpwRfleur(:)  = vpwRfleur(:)  + vCoulBenchmark(iGvec, :) * exponential
!          vpwRjuPhon(:) = vpwRjuPhon(:) + vGradCoul0IR(iGvec, :)   * exponential
!        end do
!
!      !  write (1111, '(13(es16.8E3, 2x),i2, 3(2x, i2))') x, real(vpwRfleur(1)), aimag(vpwRfleur(1)), real(vpwRfleur(2)), &
!      !    & aimag(vpwRfleur(2)), real(vpwRfleur(3)), aimag(vpwRfleur(3)), real(vpwRjuPhon(1)), aimag(vpwRjuPhon(1)),     &
!      !    & real(vpwRjuPhon(2)), aimag(vpwRjuPhon(2)), real(vpwRjuPhon(3)), aimag(vpwRjuPhon(3)), isIR, ext2bravLock(1), &
!      !    & ext2bravLock(2), ext2bravLock(3)
!
!      end if
!
!      if ( realPot ) then ! if no pseudoPotMode
!        if ( wasIR == -1 ) then ! initial point of path
!          if ( isIR == 0) then ! is in interstitial, so point can be simply plotted
!
!            vpwRfleur = cmplx(0.,0.)
!            vpwRjuPhon = cmplx(0.,0.)
!            onlyPotFLEUR = cmplx(0.,0.)
!
!            do iGvec = 1, ngdp ! calculate Interstitial Potential
!              exponential = exp(iu * tpi * dot_product(gdp(:, iGvec), ucpath))
!              vpwRfleur(:) = vpwRfleur(:) + vCoulBenchmark(iGvec, :) * exponential
!              vpwRjuPhon(:) = vpwRjuPhon(:) + vGradCoul0IR(iGvec, :) * exponential ! todo has been changed, seems to be added so that starting point ist correct!
!              !todo floating is inexact here for some reasons....
!              onlyPotFLEUR = onlyPotFLEUR +1!+ !vpw_G(iGvec) * exponentialFLEUR ! pure potential
!            end do
!
!            !write (1222, '(13(es16.8E3, 2x))') x, real(onlyPotFLEUR), aimag(onlyPotFLEUR)
!            write (1333, '(13(es16.8E3, 2x),i2, 3(2x, i2))') & ! write out calculated point
!              &x, real(vpwRfleur(1)), aimag(vpwRfleur(1)), real(vpwRfleur(2)), aimag(vpwRfleur(2)),  real(vpwRfleur(3)), &
!              &aimag(vpwRfleur(3)), real(vpwRjuPhon(1)), aimag(vpwRjuPhon(1)), real(vpwRjuPhon(2)),  aimag(vpwRjuPhon(2)), &
!              &real(vpwRjuPhon(3)), aimag(vpwRjuPhon(3)), isIR, ext2bravLock(1), ext2bravLock(2), ext2bravLock(3)
!
!          end if ! if not in interstitial we wait till reaching the MT surface enabling to solve the Dirichelet boundary problem
!        else ! wasIR is available
!          if ( ( wasIR == 0 .and. isIR == 0 )) then ! is between atoms
!
!            vpwRfleur = cmplx(0.,0.)
!            vpwRjuPhon = cmplx(0.,0.)
!            onlyPotFLEUR = cmplx(0.,0.)
!
!            do iGvec = 1, ngdp ! calculate Interstitial Potential
!              exponential = exp(iu * tpi * dot_product(gdp(:, iGvec), ucpath))
!              vpwRfleur(:) = vpwRfleur(:) + vCoulBenchmark(iGvec, :) * exponential
!              vpwRjuPhon(:) = vpwRjuPhon(:) + vGradCoul0IR(iGvec, :) * exponential
!              onlyPotFLEUR = onlyPotFLEUR + vpw_G(iGvec) * exponential ! pure potential
!            end do
!
!            !write (1222, '(13(es16.8E3, 2x))') x, real(onlyPotFLEUR), aimag(onlyPotFLEUR)
!            write (1333, '(13(es16.8E3, 2x),i2, 3(2x, i2))') & ! write out calculated point
!              &x, real(vpwRfleur(1)), aimag(vpwRfleur(1)), real(vpwRfleur(2)), aimag(vpwRfleur(2)),  real(vpwRfleur(3)), &
!              &aimag(vpwRfleur(3)), real(vpwRjuPhon(1)), aimag(vpwRjuPhon(1)), real(vpwRjuPhon(2)),  aimag(vpwRjuPhon(2)), &
!              &real(vpwRjuPhon(3)), aimag(vpwRjuPhon(3)), isIR, ext2bravLock(1), ext2bravLock(2), ext2bravLock(3)
!
!          else if (wasIR == 0 .and. isIR >= 0 ) then ! crosses MT surface into the atom
!
!            ucpathA = ucpath ! ucpathA is in atom
!            ucpathB = ucpath - [paPoXCo * dx, paPoYCo * dx, paPoZCo * dx] !todo make these dx vector general! ucpathB is outside atom CHANGE
!            do iterations = 1, 1000
!              ucpathC = (ucpathA + ucpathB) / 2
!              ucpathExtc = matmul(cellT%amat, ucpathC - atomsT%taual(:, atomLock) - ext2bravLock)
!              if ( ( abs( norm2 ( ucpathExtc ) - atomsT%rmt(itypeLock) ) < 1e-7 ) .or. norm2( ucpathB - ucpathA ) < 1e-7 ) then
!              !  ! mtsphere found
!                exit
!              else
!                if ( ( norm2(ucpathExtc) - atomsT%rmt( itypeLock ) ) < 0) then
!                  ucpathA = ucpathC
!                else
!                  ucpathB = ucpathC
!                end if
!              end if
!            end do
!            if (iterations > 1000) then
!              write(*, *) 'Warning: Iteration not converged' !todo juDFT warning
!            end if
!
!            vpwRfleur = cmplx(0.,0.)
!            vpwRjuPhon = cmplx(0.,0.)
!            onlyPotFLEUR = cmplx(0., 0.)
!
!            do iGvec = 1, ngdp ! calculate Interstitial Potential
!              exponential = exp(iu * tpi * dot_product(gdp(:, iGvec), ucpathC))
!              vpwRfleur(:) = vpwRfleur(:) + vCoulBenchmark(iGvec, :) * exponential
!              vpwRjuPhon(:) = vpwRjuPhon(:) + vGradCoul0IR(iGvec, :) * exponential
!              !onlyPotFLEUR = onlyPotFLEUR + vpw_G(iGvec) * exponential ! pure potential
!            end do
!            !write (1222, '(13(es16.8E3, 2x))') x, real(onlyPotFLEUR), aimag(onlyPotFLEUR)
!            write (1333, '(13(es16.8E3, 2x),i2, 3(2x, i2))') & ! write out calculated point
!              &x, real(vpwRfleur(1)), aimag(vpwRfleur(1)), real(vpwRfleur(2)), aimag(vpwRfleur(2)),  real(vpwRfleur(3)), &
!              &aimag(vpwRfleur(3)), real(vpwRjuPhon(1)), aimag(vpwRjuPhon(1)), real(vpwRjuPhon(2)),  aimag(vpwRjuPhon(2)), &
!              &real(vpwRjuPhon(3)), aimag(vpwRjuPhon(3)), wasIR, ext2bravLock(1), ext2bravLock(2), ext2bravLock(3)
!
!            direc = ucpathC - atomsT%taual(:, atomLock) - ext2bravLock
!            direcExt = matmul(cellT%amat, direc)
!            direcExt = direcExt / norm2(direcExt)
!!            call ylmnorm_init(atomsT%lmaxd + 1)
!            !call ylm4(atomsT%lmaxd + 1, direcExt, ylm) ! todo really until l + 1
!            call ylm4(atomsT%lmaxd, direcExt, ylm) ! todo really until l + 1
!!            call ylmnorm_init(atomsT%lmaxd)
!
!            vMTfleur = 0
!            vMTjuPhon = 0
!            do idirec = 1, 3
!              do oqn_l = 0, atomsT%lmax(itypeLock)
!                do mqn_m = -oqn_l, oqn_l
!                  lm = oqn_l * (oqn_l + 1) + 1 + mqn_m
!                  do imesh = 1, atomsT%jri(itypeLock)
!                    vMTfLEUR(idirec, imesh) = vMTFLEUR(idirec, imesh) + gradVrFleur(imesh, lm, atomLock, idirec) * ylm(lm)
!                    vMTjuPhon(idirec, imesh) = vMTjuPhon(idirec, imesh) + gradVrjuPhon(imesh, lm, idirec, atomLock) * ylm(lm)
!                  end do
!                end do
!              end do
!            end do
!
!            do imesh = atomsT%jri(itypeLock), 1, -1
!              write(1333, '(13(es20.10E3, 2x), i2, 3(2x, i2))')&
!                &-atomsT%rmsh(imesh, itypeLock), real(vMTfleur(1, imesh)), aimag(vMTfleur(1, imesh)), &
!                &real(vMTfleur(2, imesh)), aimag(vMTfleur(2, imesh)), real(vMTfleur(3, imesh)), &
!                &aimag(vMTfleur(3, imesh)), real(vMTjuPhon(1, imesh)) ,aimag(vMTjuPhon(1, imesh)),&
!                &real(vMTjuPhon(2, imesh)), aimag(vMTjuPhon(2, imesh)), real(vMTjuPhon(3, imesh)), &
!                aimag(vMTjuPhon(3, imesh)), atomLock, ext2bravLock(1), ext2bravLock(2), ext2bravLock(3)
!            end do
!
!            do idirec = 1, 3
!              if ( ( vpwRfleur(idirec) /= 0 ) .and. ( vpwRjuPhon(idirec) /= 0 ) ) then
!                write(1444, '(8(es16.8E3, 2x))') &
!                  &real(vMTfleur(idirec, atomsT%jri(itypeLock))), real(vpwRfleur(idirec)), &
!                  &abs(real(vMTfleur(idirec, atomsT%jri(itypeLock))) - real(vpwRfleur(idirec))), &
!                  &abs(real(vMTfleur(idirec, atomsT%jri(itypeLock))) - real(vpwRfleur(idirec))) / abs(real(vpwRfleur(idirec))),&
!                  &real(vMTjuPhon(idirec, atomsT%jri(itypeLock))), real(vpwRjuPhon(idirec)), &
!                  &abs(real(vMTjuPhon(idirec, atomsT%jri(itypeLock))) - real(vpwRjuPhon(idirec))), &
!                  &abs(real(vMTjuPhon(idirec, atomsT%jri(itypeLock))) - real(vpwRjuPhon(idirec))) / abs(real(vpwRjuPhon(idirec)))
!              else if ( ( vpwRfleur(idirec) == 0 ) .and. ( vpwRjuPhon(idirec) /= 0 ) ) then
!                write(1444, '(8(es16.8E3, 2x))') &
!                  &real(vMTfleur(idirec, atomsT%jri(itypeLock))), real(vpwRfleur(idirec)), &
!                  &abs(real(vMTfleur(idirec, atomsT%jri(itypeLock))) - real(vpwRfleur(idirec))), -9e9,&
!                  &real(vMTjuPhon(idirec, atomsT%jri(itypeLock))), real(vpwRjuPhon(idirec)), &
!                  &abs(real(vMTjuPhon(idirec, atomsT%jri(itypeLock))) - real(vpwRjuPhon(idirec))), &
!                  &abs(real(vMTjuPhon(idirec, atomsT%jri(itypeLock))) - real(vpwRjuPhon(idirec))) / abs(real(vpwRjuPhon(idirec)))
!              else if ( ( vpwRfleur(idirec) /= 0 ) .and. ( vpwRjuPhon(idirec) == 0 ) ) then
!                write(1444, '(8(es16.8E3, 2x))') &
!                  &real(vMTfleur(idirec, atomsT%jri(itypeLock))), real(vpwRfleur(idirec)), &
!                  &abs(real(vMTfleur(idirec, atomsT%jri(itypeLock))) - real(vpwRfleur(idirec))), &
!                  &abs(real(vMTfleur(idirec, atomsT%jri(itypeLock))) - real(vpwRfleur(idirec))) / abs(real(vpwRfleur(idirec))),&
!                  &real(vMTjuPhon(idirec, atomsT%jri(itypeLock))), real(vpwRjuPhon(idirec)), &
!                  &abs(real(vMTjuPhon(idirec, atomsT%jri(itypeLock))) - real(vpwRjuPhon(idirec))), -9e9
!              else
!                write(1444, '(8(es16.8E3, 2x))') &
!                  &real(vMTfleur(idirec, atomsT%jri(itypeLock))), real(vpwRfleur(idirec)), &
!                  &abs(real(vMTfleur(idirec, atomsT%jri(itypeLock))) - real(vpwRfleur(idirec))), -9e9,&
!                  &real(vMTjuPhon(idirec, atomsT%jri(itypeLock))), real(vpwRjuPhon(idirec)), &
!                  &abs(real(vMTjuPhon(idirec, atomsT%jri(itypeLock))) - real(vpwRjuPhon(idirec))), -9e9
!              end if
!            end do
!
!          else if (wasIR >= 0 .and. isIR == 0 ) then ! crosses MT surface out from atom
!
!            ucpathA = ucpath ! is out of atom
!            ucpathB = ucpath - [paPoXCo * dx, paPoYCo * dx, paPoZCo * dx] ! is in atom, todo make these dx vector general!
!            do iterations = 1, 1000 ! really so much iterations needed?
!              ucpathC = (ucpathA + ucpathB) / 2
!              ucpathExtc = matmul(cellT%amat, ucpathC - atomsT%taual(:, atomLock) - ext2bravLock)
!              if ( ( abs( norm2 ( ucpathExtc ) - atomsT%rmt(itypeLock) ) <= 1e-7 ) .or. norm2( ucpathB - ucpathA ) < 1e-7 ) then
!                ! mtsphere found
!                exit
!              else
!                if ( ( norm2( ucpathExtc ) - atomsT%rmt( itypeLock ) ) > 0) then
!                  ucpathA = ucpathC
!                else
!                  ucpathB = ucpathC
!                end if
!              end if
!            end do
!            if (iterations >= 1000) then
!              write(*, *) 'Warning: Iteration not converged'
!            end if
!
!            direc = ucpathC - atomsT%taual(:, atomLock) - ext2bravLock
!            direcExt = matmul(cellT%amat, direc)
!            direcExt = direcExt / norm2(direcExt)
!  !          call ylmnorm_init(atomsT%lmaxd + 1)
!            !call ylm4(atomsT%lmaxd + 1, direcExt, ylm) ! todo really until l + 1
!            call ylm4(atomsT%lmaxd, direcExt, ylm) ! todo really until l + 1
!  !          call ylmnorm_init(atomsT%lmaxd)
!
!            vMTfleur = 0
!            vMTjuPhon = 0
!            do idirec = 1, 3
!              do oqn_l = 0, atomsT%lmax(itypeLock)
!                do mqn_m = -oqn_l, oqn_l
!                  lm = oqn_l * (oqn_l + 1) + 1 + mqn_m
!                  do imesh = 1, atomsT%jri(itypeLock)
!                    vMTfLEUR(idirec, imesh) = vMTFLEUR(idirec, imesh) + gradVrFleur(imesh, lm, atomLock, idirec) * ylm(lm)
!                    vMTjuPhon(idirec, imesh) = vMTjuPhon(idirec, imesh) + gradVrjuPhon(imesh, lm, idirec, atomLock) * ylm(lm)
!                  end do
!                end do
!              end do
!            end do
!            do imesh = 1, atomsT%jri(itypeLock)
!              write(1333, '(13(es20.10E3, 2x),i2, 3(2x, i2))')&
!                &atomsT%rmsh(imesh, itypeLock), real(vMTfleur(1, imesh)), aimag(vMTfleur(1, imesh)), &
!                &real(vMTfleur(2, imesh)), aimag(vMTfleur(2, imesh)), real(vMTfleur(3, imesh)), &
!                &aimag(vMTfleur(3, imesh)), real(vMTjuPhon(1, imesh)) ,aimag(vMTjuPhon(1, imesh)),&
!                &real(vMTjuPhon(2, imesh)), aimag(vMTjuPhon(2, imesh)), real(vMTjuPhon(3, imesh)), &
!                &aimag(vMTjuPhon(3, imesh)), atomLock, ext2bravLock(1), ext2bravLock(2), ext2bravLock(3)
!            end do
!
!            !todo seems to be the real ptoential here
!            !iatom = 0
!            potFleurMT = 0
!              !Evaluate MT sphere function at same random points
!            !symAt_temp  = atomsT%ntypsy(atomLock) ! symmetry of atom iatom
!            !symOpr_temp = atomsT%ngopr (atomLock) ! symmetry operation mapping local to global coordinate system
!            !do imesh = 1, atomsT%jri(itypeLock)
!            !  rvec = direcExt
!            !  if (symOpr_temp /= 1) then
!            !  !Transform into internal real space coordinates
!            !  !todo might be wrong bmat
!            !    call cotra1(rvec, rvec_int, cellT%bmat) !allocate and deallocate randPtsGrid
!            !    do imatmulr = 1, 3
!            !      rotatedVector(imatmulr) = 0.0
!            !      do imatmulc = 1, 3
!            !        rotatedVector(imatmulr) = rotatedVector(imatmulr) + symT%mrot(imatmulr, imatmulc, symOpr_temp) * rvec_int(imatmulc)
!            !      end do
!            !    end do
!            !    call cotra0(rotatedVector, rvec, cellT%amat)
!            !  end if
!            !  call ylm4(atomsT%lmax(itypeLock), rvec, ylm)
!
!            !  do ilath = 0, sphharT%nlh(symAt_temp)
!            !    linCombBas = 0.0
!            !    l_temp = sphharT%llh(ilath, symAt_temp) * (sphharT%llh(ilath, symAt_temp) + 1) + 1 !understand this
!            !    do imem = 1, sphharT%nmem(ilath, symAt_temp)
!            !      lm_temp = l_temp + sphharT%mlh(imem, ilath, symAt_temp)
!            !      linCombBas = linCombBas + real(clnu_atom(imem, ilath, symAt_temp) * ylm(lm_temp))
!            !    end do
!            !      potFleurMT(imesh) = potFleurMT(imesh) + vrFleur(imesh, ilath, itypeLock, 1) * linCombBas
!            !  end do
!            !end do
!            !call fopen(1020, name='potMTfleur', status='replace', action='write', form='formatted')
!            !do imesh = 1, atomsT%jri(itypeLock)
!            !  write(1020, '(3(es16.8E3, 2x))') atomsT%rmsh(imesh, itypeLock), potFleurMT(imesh), real(vMTfLEUR(1, imesh))
!            !end do
!            !call fclose(1020)
!
!
!            vpwRfleur = cmplx(0.,0.)
!            vpwRjuPhon = cmplx(0.,0.)
!            onlyPotFLEUR = cmplx(0., 0.)
!            do iGvec = 1, ngdp ! calculate Interstitial Potential
!              exponential = exp(iu * tpi * dot_product(gdp(:, iGvec), ucpathC))
!              vpwRfleur(:) = vpwRfleur(:) + vCoulBenchmark(iGvec, :) * exponential
!              vpwRjuPhon(:) = vpwRjuPhon(:) + vGradCoul0IR(iGvec, :) * exponential
!              onlyPotFLEUR = onlyPotFLEUR + vpw_G(iGvec) * exponential ! pure potential
!            end do
!            !write (1222, '(13(es16.8E3, 2x))') x, real(onlyPotFLEUR), aimag(onlyPotFLEUR) ! todo the lines down really needed?
!            write (1333, '(13(es16.8E3, 2x),i2, 3(2x, i2))') & ! write out calculated point
!              &x, real(vpwRfleur(1)), aimag(vpwRfleur(1)), real(vpwRfleur(2)), aimag(vpwRfleur(2)),  real(vpwRfleur(3)), &
!              &aimag(vpwRfleur(3)), real(vpwRjuPhon(1)), aimag(vpwRjuPhon(1)), real(vpwRjuPhon(2)),  aimag(vpwRjuPhon(2)), &
!              &real(vpwRjuPhon(3)), aimag(vpwRjuPhon(3)), isIR, ext2bravLock(1), ext2bravLock(2), ext2bravLock(3)
!
!            do idirec = 1, 3
!              if ( ( vpwRfleur(idirec) /= 0 ) .and. ( vpwRjuPhon(idirec) /= 0 ) ) then
!                write(1444, '(8(es16.8E3, 2x))') &
!                  &real(vMTfleur(idirec, atomsT%jri(itypeLock))), real(vpwRfleur(idirec)), &
!                  &abs(real(vMTfleur(idirec, atomsT%jri(itypeLock))) - real(vpwRfleur(idirec))), &
!                  &abs(real(vMTfleur(idirec, atomsT%jri(itypeLock))) - real(vpwRfleur(idirec))) / abs(real(vpwRfleur(idirec))),&
!                  &real(vMTjuPhon(idirec, atomsT%jri(itypeLock))), real(vpwRjuPhon(idirec)), &
!                  &abs(real(vMTjuPhon(idirec, atomsT%jri(itypeLock))) - real(vpwRjuPhon(idirec))), &
!                  &abs(real(vMTjuPhon(idirec, atomsT%jri(itypeLock))) - real(vpwRjuPhon(idirec))) / abs(real(vpwRjuPhon(idirec)))
!              else if ( ( vpwRfleur(idirec) == 0 ) .and. ( vpwRjuPhon(idirec) /= 0 ) ) then
!                write(1444, '(8(es16.8E3, 2x))') &
!                  &real(vMTfleur(idirec, atomsT%jri(itypeLock))), real(vpwRfleur(idirec)), &
!                  &abs(real(vMTfleur(idirec, atomsT%jri(itypeLock))) - real(vpwRfleur(idirec))), -9e9,&
!                  &real(vMTjuPhon(idirec, atomsT%jri(itypeLock))), real(vpwRjuPhon(idirec)), &
!                  &abs(real(vMTjuPhon(idirec, atomsT%jri(itypeLock))) - real(vpwRjuPhon(idirec))), &
!                  &abs(real(vMTjuPhon(idirec, atomsT%jri(itypeLock))) - real(vpwRjuPhon(idirec))) / abs(real(vpwRjuPhon(idirec)))
!              else if ( ( vpwRfleur(idirec) /= 0 ) .and. ( vpwRjuPhon(idirec) == 0 ) ) then
!                write(1444, '(8(es16.8E3, 2x))') &
!                  &real(vMTfleur(idirec, atomsT%jri(itypeLock))), real(vpwRfleur(idirec)), &
!                  &abs(real(vMTfleur(idirec, atomsT%jri(itypeLock))) - real(vpwRfleur(idirec))), &
!                  &abs(real(vMTfleur(idirec, atomsT%jri(itypeLock))) - real(vpwRfleur(idirec))) / abs(real(vpwRfleur(idirec))),&
!                  &real(vMTjuPhon(idirec, atomsT%jri(itypeLock))), real(vpwRjuPhon(idirec)), &
!                  &abs(real(vMTjuPhon(idirec, atomsT%jri(itypeLock))) - real(vpwRjuPhon(idirec))), -9e9
!              else
!                write(1444, '(8(es16.8E3, 2x))') &
!                  &real(vMTfleur(idirec, atomsT%jri(itypeLock))), real(vpwRfleur(idirec)), &
!                  &abs(real(vMTfleur(idirec, atomsT%jri(itypeLock))) - real(vpwRfleur(idirec))), -9e9,&
!                  &real(vMTjuPhon(idirec, atomsT%jri(itypeLock))), real(vpwRjuPhon(idirec)), &
!                  &abs(real(vMTjuPhon(idirec, atomsT%jri(itypeLock))) - real(vpwRjuPhon(idirec))), -9e9
!              end if
!            end do
!
!          end if! if still in MT then wait for crossing
!
!        end if !is not first point on path
!      end if ! is not in pseudopot mode
!      wasIR = isIR
!    end do
!    !call fclose(1111)
!    !call fclose(1222)
!    call fclose(1333)
!    call fclose(1444)
  end subroutine checkGradPot0sXSF


  !--------------------------------------------------------------------------------------------------------------------------------------
  !> @author
  !> Christian-Roman Gerhorst, Forschungszentrum Jülich: IAS1 / PGI1
  !>
  !> @brief
  !> Test of the routines involved in calculating the first variation of the
  !> Coulomb potential.
  !>
  !> @details
  !> The idea is to use the routines for calculating the first variation of the
  !> Coulomb potential as much as possible in order to acquire the unperturbed
  !> Coulomb potential both in the interstitial and the muffin-tin region. It
  !> then can be compared to the Coulomb potentials calculated with FLEUR. This
  !> routine tests the routines m_jppotdens::mpmom1hart,
  !> m_jppotdens::pwmom1hart(), m_jppotdens::psqpwvec() and
  !> m_jppotdens::vmtsnsym().
  !>
  !> @note This test requires additional output of FLEUR which can be compared
  !> to the results of the juPhon routines.
  !>
  !> @param[in] atomsT     : Contains atoms-related quantities; definition of its members in types.F90 file.
  !> @param[in] latharT    : Contains quantities describing the lattice harmonics; definition of its members in types.F90 file.
  !> @param[in] cellT      : Contains unit-cell related quantities; definition of its members in types.F90 file.
  !> @param[in] starsT     : Contains stars-related quantities; definition of its members in types.F90 file.
  !> @param[in] dimensT    : Contains dimension-related quantities; definition of its members in types.F90 file.
  !> @param[in] inputT     : Contains input-related quantities; definition of its members in types.F90 file.
  !> @param[in] rho0MT     : Electronic density in muffin-tin regions; see for dimensions.
  !> @param[in] clnu_atom  : Expansion coefficients of lattice harmonics, resolved for every atom, see for dimensions.
  !> @param[in] nmem_atom  : Number of members per lattice harmonics, resolved for every atom, see for dimensions.
  !> @param[in] mlh_atom   : Magnetic quantum numbers of lattice harmonic members, resolved for every atom see for dimensions.
  !> @param[in] gdp        : Contains G-vectors used for density and potential, see for dimensions.
  !> @param[in] ngpd       : Number of G-vectors used for density and potential.
  !> @param[in] rho0PW     : Electronic density in the interstitial region, see for dimensions.
  !--------------------------------------------------------------------------------------------------------------------------------------
  subroutine checkVCoul1(atomsT, latharT, cellT, starsT, dimensT, inputT, rho0MT, clnu_atom, nmem_atom, mlh_atom, gdp, ngdp, rho0PW, logUnit)

    use m_jPConstants,   only : fpi
    use m_types
    use mod_juPhonUtils, only : fopen, fclose
    use m_jpVeff1,      only : mpMom1Hart, pwMom1Hart, psqpwVec, vmtsNSym
    use m_juDFT_NOstopNO, only : juDFT_warn
    use m_jpPotDensHelper, only : convertstar2g

    implicit none

    ! Type parameters
    type(t_atoms),      intent(in) :: atomsT
    type(t_cell),       intent(in) :: cellT
    type(t_dimension),  intent(in) :: dimensT
    type(t_input),      intent(in) :: inputT
    type(t_sphhar),     intent(in) :: latharT
    type(t_stars),      intent(in) :: starsT

    ! Scalar parameters
    integer,            intent(in) :: ngdp
    integer,            intent(in) :: logUnit

    ! Array parameters
    integer,            intent(in) :: gdp(:, :)
    integer,            intent(in) :: mlh_atom(:, 0:, :)
    integer,            intent(in) :: nmem_atom(0:, :)
    real,               intent(in) :: rho0MT(:, 0:, :, :)
    complex,            intent(in) :: clnu_atom(:, 0:, :)
    complex,            intent(in) :: rho0PW(:, :)

    !------------------------------------------------------------------------------------------------------------------------------------
    ! Scalar local variabels
    !
    ! iatom  : runs over all atoms in unit cell
    ! idirec : runs over all directions the atoms are displaced to
    ! iGvec  : runs over all G-vectors used for the potentials and the densities
    ! ilh    : runs over all lattice harmonics
    ! imem   : runs over all members of a lattice harmonics
    ! imesh  : runs over muffin-tin mesh
    ! ieq    : runs over all equivalent atoms of an atom type
    ! itype  : runs over all atom types
    ! lm     : index describing the orbital quantum number and the magnetic quantum number
    ! mqn_m  : magnetic quantum number
    ! oqn_l  : orbital quantum number
    ! pre_lm : auxiliary index for the lms
    ! ptsym  : contains index of point symmetry
    ! noGext : norm of a external G-vector used for the potentials and the density
    !------------------------------------------------------------------------------------------------------------------------------------
    integer                        :: iatom
    integer                        :: idirec
    integer                        :: iGvec
    integer                        :: ilh
    integer                        :: imem
    integer                        :: imesh
    integer                        :: ieq
    integer                        :: itype
    integer                        :: lm
    integer                        :: mqn_m
    integer                        :: oqn_l
    integer                        :: pre_lm
    integer                        :: ptsym
    real                           :: noGext
    logical                        :: passed = .true.
    logical                        :: vHarNum

    !------------------------------------------------------------------------------------------------------------------------------------
    ! Local array variabels
    !
    ! filename        : contains filenames of output files.
    ! Gext            : contains a external G-vector.
    ! vCoulMTFLEURLH  : Lattice harmonic expansion coefficients of the Coulomb potential in the muffin-tin region calculated with FLEUR.
    ! vCoulMTFLEURSpH : Spher. harmonic expansion coefficients of the Coulomb potential in the muffin-tin region calculated with juPhon.
    ! vCoulMTJP       : Expansion coefficients of Coulomb potential in muffin tins calculated with juPhon.
    ! psqpwJP         : Pseudo charges of the Coulomb potential calculated with juPhon.
    ! psqpwFLEUR      : Pseudo charges of the Coulomb potential calculated with FLEUR expanded in stars.
    ! psqpwFLEURGs    : Pseudo charges of the Coulomb potential calculated with FLEUR expanded in plane waves.
    ! qlmCoul         : Multipole moments of the Coulomb potential.
    ! qlmoCoul        : Muffin-tin part of the Coulomb potential multipole moments.
    ! qlmpCoul        : Plane-wave part of the Coulomb potential multipole moments.
    ! rho0MTSpH       : Spherical harmonic expansion coeficients of unperturbed density.
    ! rho0PWG         : Plane wave expansion coefficients of interstitial unperturbed density.
    ! rho0PWG3D       : Same as rho0PW but copied to other 2 dimensions.
    ! vpwFLEUR        : Expansion coefficients of interstitial Coulomb potential calculated with FLEUR.
    ! vpwJP           : Expansion coefficients of interstitial Coulomb potential calculated with juPhon.
    !------------------------------------------------------------------------------------------------------------------------------------
     ! todo put this to heap with allocate!
    character(len=15)              :: filename
    real                           :: Gext(3)
    real                           :: vCoulMTFLEURLH(atomsT%jmtd, 0: latharT%nlhd, atomsT%ntype, inputT%jspins)
    complex                        :: vCoulMTFLEURSpH(atomsT%jmtd, (atomsT%lmaxd + 1)**2, atomsT%ntype)
    complex                        :: psqpwJP(3, ngdp)
    complex                        :: psqpwFLEUR(starsT%n3d)
    complex                        :: psqpwFLEURGs(ngdp)
!    complex                        :: qlmCoul(-atomsT%lmaxd-1:atomsT%lmaxd+1,0:atomsT%lmaxd+1,atomsT%nat, 3)
!    complex                        :: qlmoCoul(-atomsT%lmaxd-1:atomsT%lmaxd+1,0:atomsT%lmaxd+1,atomsT%nat, 3)
!    complex                        :: qlmpCoul(-atomsT%lmaxd-1:atomsT%lmaxd+1,0:atomsT%lmaxd+1,atomsT%nat, 3)
    complex, allocatable           :: qlmCoul(:,:, :)
    complex, allocatable           :: qlmoCoul(:,:, :)
    complex, allocatable           :: qlmpCoul(:,:, :)
    complex                        :: rho0MTSpH(atomsT%jmtd, (atomsT%lmaxd + 2)**2, atomsT%nat, 3)
    complex                        :: rho0PWG(ngdp)
    complex                        :: rho0PWG3D(ngdp, 3)
    complex                        :: vCoulMTJP(atomsT%jmtd, (atomsT%lmaxd + 1)**2, 3, atomsT%ntype)
    complex                        :: vpwFLEUR(ngdp)
    complex                        :: vpwJP(ngdp, 3)
    complex,                      allocatable  :: grRho0MT(:, :, :, :)
    integer                         :: idir
    logical                         :: linIntp

    allocate( qlmCoul((atomsT%lmaxd + 1)**2,  atomsT%nat, 3))
    allocate( qlmoCoul((atomsT%lmaxd + 1)**2, atomsT%nat, 3))
    allocate( qlmpCoul((atomsT%lmaxd + 1)**2, atomsT%nat, 3))
    qlmCoul(:, :, :) = cmplx(0., 0.)
    qlmoCoul(:, :, :) = cmplx(0., 0.)
    qlmpCoul(:, :, :) = cmplx(0., 0.)

    write(logUnit, '(a)') 'Calculate FLEUR Coulomb potential with V_coulomb1 routines test'
    write(logUnit, '(a)')    '---------------------------------------------------------------'
    allocate( grRho0MT(atomsT%jmtd, ( atomsT%lmaxd + 1 )**2, atomsT%nat, 3) )
    grRho0MT(:, :, :, :) = cmplx(0.0, 0.0)

    ! Unfold symmetry-compressed muffin-tin density because juPhon routines do not implement symmetry
    rho0MTSpH = 0
    do idirec = 1, 3
      iatom = 1
      do itype = 1, atomsT%ntype
        do ieq = 1, atomsT%neq(itype)
          ptsym = atomsT%ntypsy(iatom)
          do ilh = 0, latharT%nlh(ptsym)
            oqn_l = latharT%llh(ilh, ptsym)
            pre_lm = oqn_l * (oqn_l + 1) + 1
            do imem = 1, nmem_atom(ilh, iatom)
              mqn_m = mlh_atom(imem, ilh, iatom)
              lm = pre_lm + mqn_m
              do imesh = 1, atomsT%jri(itype)
                rho0MTSpH(imesh, lm, iatom, idirec) = rho0MTSpH(imesh, lm, iatom, idirec) + rho0MT(imesh, ilh, itype, 1)&
                  &                                                                                         * clnu_atom(imem, ilh, iatom)
              end do ! imesh
            end do ! imem
          end do ! ilh
          iatom = iatom + 1
        end do ! ieq
      end do ! itype
    end do ! idirec

    ! Calculates muffin-tin multipole moments for the Hartree potential
    ! the two 1 are for the displaced atom type and displaced atom where grRho0MT is subtracted which is zero in this test in order that it works, so setting to 1 and 1 leads to a well-defined indexing of the array, because one atom and one atom type is always defined.
    call mpMom1Hart(atomsT, 1, 1, rho0MTSpH, grRho0MT, qlmoCoul)

    ! todo discuss with gustav
    ! Adding muffin-tin multipole moments for the external potential to gain the Coulomb potential multipole moments
    do idirec = 1, 3
      iatom = 1
      do itype = 1, atomsT%ntype
        do ieq = 1, atomsT%neq(itype)
          qlmoCoul(1, iatom, idirec) = qlmoCoul(1, iatom, idirec) - atomsT%zatom(itype) / sqrt(fpi) ! subtract contribution of ions
          iatom = iatom + 1
        end do ! ieq
      end do ! itype
    end do ! idirec


    ! Unfold symmetry-compressed plane-wave density since juPhon routines do not implement symmetry, then expand it to three dimensions.
    call convertStar2G(rho0PW(:, 1), rho0PWG, starsT, ngdp, gdp)
    do idirec = 1, 3
      rho0PWG3D(:, idirec) = rho0PWG(:)
    end do ! idirec

    ! Calculates plane-wave multipole moments for the Hartree potential.
    call pwMom1Hart(atomsT, cellT, gdp, ngdp, [0., 0., 0.], rho0PWG3D, qlmpCoul)
!    do idirec = 1, 3
!      iatom = 1
!      do itype = 1, atomsT%ntype
!        do ieq = 1, atomsT%neq(itype)
!          oqn_l = atomsT%lmax(itype) + 1
!          do mqn_m = -oqn_l, oqn_l
!          lm = oqn_l * (oqn_l + 1) + 1 + mqn_m
!          ! We don't need the l + 1 contribution here as this is due to the Sternheimer equation and would deliver wrong results
!          ! if we benchmark this array to the FLEUR results
!          qlmpCoul(lm, iatom, idirec) = cmplx(0.0, 0.0)
!          end do ! imem
!          iatom = iatom + 1
!        end do ! ieq
!      end do ! itype
!    end do ! idirec

    ! todo discuss with Gustav
    ! Add G = 0 contribution
    do idirec = 1, 3
      iatom = 1
      do itype = 1, atomsT%ntype
        do ieq = 1, atomsT%neq(itype)
          qlmpCoul(1, iatom, idirec) = qlmpCoul(1, iatom, idirec) + rho0PW(1, 1) * starsT%nstr(1) * atomsT%volmts(itype) / sqrt(fpi)
          iatom = iatom + 1
        end do ! ieq
      end do ! itype
    end do ! idirec

    ! Subtract plane-wave contribution in muffin-tins to get correct Coulomb multipole moments
    qlmCoul = qlmoCoul - qlmpCoul

!    ! Output of multipole moments
!    call fopen(1000, name='mpmomjP', status='replace', action='write', form='formatted')
!    iatom = 1
!    do itype = 1, atomsT%ntype
!      write (1000, *) itype
!      ptsym = atomsT%ntypsy(iatom)
!      do ilh = 0, latharT%nlh(ptsym)
!        oqn_l = latharT%llh(ilh,ptsym)
!        do imem = 1, latharT%nmem(ilh,ptsym)
!          mqn_m = latharT%mlh(imem, ilh, ptsym)
!          write (1000, '(1x,i2,2x,i2,2x,2 (5x,2e15.5))') oqn_l, mqn_m, qlmoCoul(mqn_m, oqn_l, iatom, 1), &
!            &                                                                                          qlmpCoul(mqn_m, oqn_l, iatom, 1)
!        end do ! imem
!      end do ! ilh
!      iatom = iatom + atomsT%neq(itype)
!    end do ! itype
!    call fclose(1000)


    ! Calculates pseudo charges from the multipole moments.
    call psqpwVec(atomsT, cellT, dimensT, gdp, [0., 0., 0.], ngdp, rho0PWG3D, qlmCoul, psqpwJP)

    ! Read in FLEUR-calculated pseudo-charge
    call fopen(1000, name='rho0pwFLEUR', status='old', action='read', form='unformatted')
    read(1000) psqpwFLEUR
    call fclose(1000)
    call convertStar2G(psqpwFLEUR, psqpwFLEURGs, starsT, ngdp, gdp)

    ! Calculate interstitial potential from pseudocharge
    vpwJP = 0
    vpwFLEUR = 0
    do iGvec = 1, ngdp
      Gext = matmul(cellT%bmat, gdp(:, iGvec))
      noGext = norm2(Gext)
      if ( noGext == 0 ) cycle
      do idirec = 1, 3
        vpwJP(iGvec, idirec) = fpi * psqpwJP(idirec, iGvec) / noGext**2
      end do ! idirec
        vpwFLEUR(iGvec) = fpi * psqpwFLEURGs(iGvec) / noGext**2
    end do ! iGvec

    ! Compare interstitial Coulomb potential from FLEUR with juPhon-calculated one.
    do iGvec = 1, ngdp
      if ( abs(vpwJP(iGvec, 1) - vpwFLEUR(iGvec)) >= 1e-9 ) then
        passed = .false.
        write (logUnit, *) 'Unperturbed potential calculated with juPhon and FLEUR do not coincide at ', iGvec, 'th G-vector.'
      end if
    end do ! iGvec


    !! Calculate MT potential solving the Dirichelet boundary problem
    idir = 1
    iatom = 1
    vHarNum = .false.
    linIntp = .false. ! Linear interpolation does not lead to the same results as in FLEUR.
    do itype = 1, atomsT%ntype
    !TODO check this test before commiting!!!!!!!!!!!!!!!!!!!!!!!
      call vmtsNSym(atomsT, cellT, ngdp, idir, iatom, itype, 1, .false., gdp, [0., 0., 0.], vpwJP, rho0MTSpH, vCoulMTJP(:, :, itype, 1), .true., .true., vHarNum, linIntp)
      iatom = iatom + atomsT%neq(itype)
    end do ! itype

    ! todo gustav
    ! Add G = 0 component
    iatom = 1
    do iGvec = 1, ngdp
      if ( starsT%ig(gdp(1, iGvec), gdp(2, iGvec), gdp(3, iGvec)) == 1) then
        do itype = 1, atomsT%ntype
          oqn_l = latharT%llh(0, atomsT%ntypsy(iatom))
          do imem = 1, nmem_atom(0, iatom)
            mqn_m = mlh_atom(imem, 0, iatom)
            lm = oqn_l * (oqn_l + 1) + 1 + mqn_m
            do imesh = 1, atomsT%jri(itype)
              vCoulMTJP(imesh, lm, 1, itype) = vCoulMTJP(imesh, lm, 1, itype) + sqrt(fpi) * vpwJP(iGvec, 1)&
                &                                                           * atomsT%rmsh(imesh, itype)**oqn_l / atomsT%rmt(itype)**oqn_l
            end do !imesh
          end do !imem
          iatom = iatom + atomsT%neq(itype)
        end do !itype
      end if
    end do !iGvec

    !todo discuss with gustav: g = 0 lacks in psqpw, not implemented here at all

    ! Adding additional external potential contribution to gain the correct Coulomb potential. These are only calculated for every atom
    ! type because FLEUR only calculates it per atom type
    iatom = 1
    do itype = 1, atomsT%ntype
      oqn_l = latharT%llh(0, atomsT%ntypsy(iatom))
      do imem = 1, nmem_atom(0, iatom)
        mqn_m = mlh_atom(imem, 0, iatom)
        lm = oqn_l * (oqn_l + 1) + 1 + mqn_m
        do imesh = 1, atomsT%jri(itype)
          vCoulMTJP(imesh, lm, 1, itype) = vCoulMTJP(imesh, lm, 1, itype) - sqrt(fpi) * (1 - atomsT%rmsh(imesh, itype)&
            &                                                      / atomsT%rmt(itype)) / atomsT%rmsh(imesh, iatom) * atomsT%zatom(itype)
        end do ! imesh
      end do ! imem
      iatom = iatom + atomsT%neq(itype)
    end do !itype

    ! read in FLEUR Coulomb muffin-tin potential
    call fopen(1000, name='v0MTFLEUR_coul', status='old', action='read', form='unformatted')
    read(1000) vCoulMTFLEURLH
    call fclose(1000)

    ! Unfold symmetry-compressed FLEUR-calculated Coulomb potential into spherical harmonics to compare with juPhon result
    vCoulMTFLEURSpH = 0
    iatom = 1
    do itype = 1, atomsT%ntype
      ptsym = atomsT%ntypsy(iatom)
      do ilh = 0, latharT%nlh(ptsym)
        oqn_l = latharT%llh(ilh, ptsym)
        pre_lm = oqn_l * (oqn_l + 1) + 1
        do imem = 1, nmem_atom(ilh, iatom)
          mqn_m = mlh_atom(imem, ilh, iatom)
          lm = pre_lm + mqn_m
          do imesh = 1, atomsT%jri(itype)
            vCoulMTFLEURSpH(imesh, lm, itype) = vCoulMTFLEURSpH(imesh, lm, itype) + vCoulMTFLEURLH(imesh, ilh, itype, 1)&
              &                                                                                             * clnu_atom(imem, ilh, iatom)
          end do !imesh
        end do ! imem
      end do ! ilh
      iatom = iatom + atomsT%neq(itype)
    end do ! itype

    do itype = 1, atomsT%ntype
      do lm = 1, (atomsT%lmax(itype) + 1)**2
        do imesh = 1, atomsT%jri(itype)
          if (abs(vCoulMTFLEURSpH(imesh, lm, itype) - vCoulMTJP(imesh, lm, 1, itype)) > 3e-6) then
            write(*, '(2f15.8)') abs(vCoulMTFLEURSpH(imesh, lm, itype) - vCoulMTJP(imesh, lm, 1, itype))
            passed = .false.
          end if
        end do
      end do
    end do
    ! Output of FLEUR-calcualted and juPhon-calculated muffin-tin Coulomb potential to compare it in plots.
!    do itype = 1, atomsT%ntype
!      write(*, *) 'lm cols: ', (atomsT%lmax(itype) + 1)**2
!      write (filename, "(A10,I1,A4)") 'vCoul0JPFL', itype, '.ssv'
!      call fopen(1000, name=filename, status='replace', action='write', form='formatted')
!      write(1000, '(A67,I1,A25,I1,A11,I3)') '# Unperturbed Coulomb potential from FLEUR and juphon for atom type', itype, &
!        & ' and its equivalent atom ', ieq, ', lm cols: ', (atomsT%lmax(itype) + 1)**2
!      do imesh = 1, atomsT%jri(itype)
!        write(1000,'(2x,*(es15.7))') atomsT%rmsh(imesh, itype), &
!          & (real(vCoulMTFLEURSpH(imesh, lm, itype)), lm = 1, (atomsT%lmax(itype) + 1)**2), &
!          & (real(vCoulMTJP(imesh, lm, 1, itype)), lm = 1, (atomsT%lmax(itype) + 1)**2), &
!          & (aimag(vCoulMTFLEURSpH(imesh, lm, itype)), lm = 1, (atomsT%lmax(itype) + 1)**2), &
!          & (aimag(vCoulMTJP(imesh, lm, 1, itype)), lm = 1, (atomsT%lmax(itype) + 1)**2)
!      end do ! imesh
!      call fclose(1000)
!    end do ! itype
    if (passed) then
      write (logUnit, '(a)')   '                                                              |__ passed!'
    else
      write (logUnit, '(a)')   '                                                              |__ failed!'
      call JuDFT_warn('Calculate FLEUR Coulomb potential with V_coulomb1 routines test.', calledby='CheckVCoul1', hint='Debug!')
    end if

  end subroutine checkVCoul1




  !--------------------------------------------------------------------------------------------------------------------------------------
  !> @author
  !> Christian-Roman Gerhorst, Forschungszentrum Jülich: IAS1 / PGI1
  !>
  !> @brief
  !> Tests the method mod_juphonutils::convertstar2g
  !>
  !> @details
  !> This method tests the method mod_juphonutils::convertstar2g by calling it with the interstitial potential of FLEUR and then writing
  !> out the transformed expansion coefficients to pottotGvec.ssv. At the same time the star expansion coefficients are written out. To
  !> ensure a comparability the star expansion coefficients should stem from the same system, which has been calculated without symmetry,
  !> so that every star contains only one component, i.e. is equivalent to a plane-wave approach.
  !>
  !> @note
  !> This routine needs 2 FLEUR calculations one with and one without symmetry.
  !>
  !> @param[in] vpw        : Unperturbed plane-wave potential
  !> @param[in] starsT     : Contains stars-related quantities; definition of its members in types.F90 file.
  !> @param[in] inputT     : Contains input-related quantities; definition of its members in types.F90 file.
  !> @param[in] ngpd       : Number of G-vectors used for density and potential.
  !> @param[in] gdp        : Contains G-vectors used for density and potential, see for dimensions.
  !--------------------------------------------------------------------------------------------------------------------------------------
  subroutine checkIRSta2GConv(vpw, starsT, inputT, ngdp, gdp, logUnit)

    use m_jpPotDensHelper, only : ConvertStar2G
    use m_types

    implicit none

    ! Scalar parameters
    type(t_stars), intent(in) :: starsT
    type(t_input), intent(in) :: inputT
    integer,       intent(in) :: ngdp
    integer,       intent(in) :: logUnit

    ! Array parameters
    complex,       intent(in) :: vpw(:, :)
    integer,       intent(in) :: gdp(:, :)

    ! Local scalar parameter
    integer                   :: jspin
    integer                   :: iGvec

    ! Local array parameter
    complex                   :: vpwG(ngdp, inputT%jspins)

    write(logUnit, * )
    write(logUnit, '(a)') '--------------------------------------------------------------------------------------------------------&
                                                                                                                               &---'
    write(logUnit, '(a)') 'Test unfolding of the stars-expanded into G-vector expanded quantities.'

    ! Convert interstitial potential
    do jspin = 1, inputT%jspins
      call convertStar2G(vpw(:, jspin), vpwG(:, jspin), starsT, ngdp, gdp)
    end do

    ! Write out plane-wave expansion coefficients
    call fopen(1000, name='pottotGvec.ssv', status='replace', action='write', form='formatted')
    do iGvec = 1, ngdp
      write(1000, '(*(es16.8,2x,es16.8))') (real(vpwG(iGvec, jspin)), jspin = 1, inputT%jspins),&
        &                                  (aimag(vpwG(iGvec, jspin)), jspin = 1, inputT%jspins)
    end do
    call fclose(1000)


    !Write out star expansion coefficients
    call fopen(1000, name='pottotStars.ssv', status='replace', action='write', form='formatted')
    do iGvec = 1, ngdp
      write(1000, '(*(es16.8,2x,es16.8))') &
        &                         (real(vpw(starsT%ig(gdp(1, iGvec), gdp(2, iGvec), gdp(3, iGvec)), jspin)), jspin = 1, inputT%jspins),&
        &                         (aimag(vpw(starsT%ig(gdp(1, iGvec), gdp(2, iGvec), gdp(3, iGvec)), jspin)), jspin = 1, inputT%jspins)
    end do
    call fclose(1000)

    write(logUnit, '(a)') '  --> Converge the same system without symmetry, i.e. only using the identity, in a subfolder called noS&
                                                                                                                               &ym.'
    write(logUnit, '(a)') '  --> Execute juPhon with testUnfoldStarsSw activated within the noSym folder.'
    write(logUnit, '(a)') '  --> Execute starsTest.py in the parent folder to check correctnes of unfolding.'
    write(logUnit, '(a)') '--------------------------------------------------------------------------------------------------------&
                                                                                                                               &---'
  end subroutine checkIRSta2GConv


  !--------------------------------------------------------------------------------------------------------------------------------------
  !> @author
  !> Christian-Roman Gerhorst, Forschungszentrum Jülich: IAS1 / PGI1
  !>
  !> @brief
  !> Checks mod_juphonutils::derivative
  !>
  !> @details
  !> The routine mod_juphonutils::derivative is tested by deriving \f$\sin(x)\f$, \f$\cos(x)\f$, \f$x\f$ and \f$1\f$. The functions and
  !> its derivatives are plotted to testDeriv_sin.ssv, testDeriv_cos.ssv testDeriv_x.ssv and testDeriv_1.ssv.
  !>
  !> @param[in] atomsT     : Contains atoms-related quantities; definition of its members in types.F90 file.
  !--------------------------------------------------------------------------------------------------------------------------------------
  subroutine checkDerivative(atomsT, logUnit)

    use mod_juPhonUtils

    implicit none

    ! Scalar parameter
    type(t_atoms), intent(in) :: atomsT

    ! Local scalar variables ! todo documentation
    integer                   :: imesh
    integer                   :: logUnit

    !Local array variables
    real                      :: f1(atomsT%jri(1))
    real                      :: f2(atomsT%jri(1))
    real                      :: f3(atomsT%jri(1))
    real                      :: f4(atomsT%jri(1))
    real                      :: f5(atomsT%jri(1))
    real                      :: f6(atomsT%jri(1))
    real                      :: df1(atomsT%jri(1))
    real                      :: df2(atomsT%jri(1))
    real                      :: df3(atomsT%jri(1))
    real                      :: df4(atomsT%jri(1))
    real                      :: df5(atomsT%jri(1))
    real                      :: df6(atomsT%jri(1))


    f1(:) = 0.
    f2(:) = 0.
    f3(:) = 0.
    f4(:) = 0.
    f5(:) = 0.
    f6(:) = 0.

    ! Create test functions
    do imesh = 1, atomsT%jri(1)
      f1(imesh) = sin(atomsT%rmsh(imesh, 1))
      f2(imesh) = cos(atomsT%rmsh(imesh, 1))
      f3(imesh) = atomsT%rmsh(imesh, 1)
      f4(imesh) = 1
      f5(imesh) = -1. / atomsT%rmsh(imesh, 1)
      f6(imesh) = exp(-atomsT%rmsh(imesh, 1))
    end do ! imesh


    ! Derive test functions
    call derivative(f1, 1, atomsT, df1)
    call derivative(f2, 1, atomsT, df2)
    call derivative(f3, 1, atomsT, df3)
    call derivative(f4, 1, atomsT, df4)
    call derivative(f5, 1, atomsT, df5)
    call derivative(f6, 1, atomsT, df6)


    ! Write out test functions and their derivatives
    call fopen (1000, name='testDeriv_sin.ssv', status='replace', action='write', form='formatted')
    do imesh = 1, atomsT%jri(1)
      write(1000, '(es28.14,2x,es28.14,2x,es28.14)') atomsT%rmsh(imesh, 1), f1(imesh), df1(imesh)
    end do ! imesh
    call fclose(1000)

    call fopen (1000, name='testDeriv_cos.ssv', status='replace', action='write', form='formatted')
    do imesh = 1, atomsT%jri(1)
      write(1000, '(es28.14,2x,es28.14,2x,es28.14)') atomsT%rmsh(imesh, 1), f2(imesh), df2(imesh)
    end do ! imesh
    call fclose(1000)

    call fopen (1000, name='testDeriv_x.ssv', status='replace', action='write', form='formatted')
    do imesh = 1, atomsT%jri(1)
      write(1000, '(es28.14,2x,es28.14,2x,es28.14)') atomsT%rmsh(imesh, 1), f3(imesh), df3(imesh) - 1.
    end do ! imesh
    call fclose(1000)

    call fopen (1000, name='testDeriv_1.ssv', status='replace', action='write', form='formatted')
    do imesh = 1, atomsT%jri(1)
      write(1000, '(es28.14,2x,es28.14,2x,es28.14)') atomsT%rmsh(imesh, 1), f4(imesh), df4(imesh)
    end do !imesh
    call fclose(1000)

    call fopen (1000, name='testDeriv_r-1.ssv', status='replace', action='write', form='formatted')
    do imesh = 1, atomsT%jri(1)
      write(1000, '(es28.14,2x,es28.14,2x,es28.14)') atomsT%rmsh(imesh, 1), f5(imesh), (abs(df5(imesh) - (1 / atomsT%rmsh(imesh, 1)**2)) * atomsT%rmsh(imesh, 1)**2)
    end do !imesh
    call fclose(1000)

    call fopen (1000, name='testDeriv_exp-r.ssv', status='replace', action='write', form='formatted')
    do imesh = 1, atomsT%jri(1)
      write(1000, '(es28.14,2x,es28.14,2x,es28.14)') atomsT%rmsh(imesh, 1), f6(imesh), (abs(df6(imesh) + exp(-atomsT%rmsh(imesh, 1))) / exp(-atomsT%rmsh(imesh, 1)))
    end do !imesh
    call fclose(1000)

    write(logUnit, *)
    write(logUnit, '(a)') '--------------------------------------------------------------------------------------------------------&
                                                                                                                          &--------'
    write(logUnit, '(a)') 'Test of radial derivative'
    write(logUnit, '(a)') '  --> Check created plot data of test functions and their numerical derivatives plotting with testDeriva&
                                                                                                                          &tives.py'
    write(logUnit, '(a)') '  --> Note: Test functions are only calculated until muffin-tin radius.'
    write(logUnit, '(a)') '--------------------------------------------------------------------------------------------------------&
                                                                                                                          &--------'
    write(logUnit, *)

  end subroutine checkDerivative

  ! Compares the output of the routine with manually set up lm channels of the MT gradient
  subroutine checkMTDerivative(atomsT, latharT, nmem_atom, mlh_atom, memd_atom, logUnit)

    use m_gaunt, only : gaunt1
    use m_types
    use m_jpPotDensHelper
    use m_juDFT_NOstopNO, only : juDFT_warn

    implicit none

    type(t_atoms),  intent(in)  :: atomsT
    type(t_sphhar), intent(in)  :: latharT
    integer,        intent(in)  :: logUnit
    integer,        intent(in)  :: nmem_atom(0 :, :)
    integer,        intent(in)  :: mlh_atom(:, 0 :, :)

    integer,        intent(in)  :: memd_atom

    complex                     :: clnu_atom(memd_atom, 0:latharT%nlhd, atomsT%nat)
    real                        :: testFunc(atomsT%jmtd, 0:latharT%nlhd, atomsT%ntype)
    complex,        allocatable :: gradrho0MT(:, :, :, :)

    integer                     :: imesh
    integer                     :: iatom
    integer                     :: idirec
    integer                     :: itype
    integer                     :: ieq
    integer                     :: lm
    integer                     :: oqn_l
    integer                     :: mqn_m
    logical                     :: realL
    logical                     :: imagL
    integer                     :: lh
    integer                     :: mem
    logical                     :: passed = .true.

    write(logUnit, *)
    write(logUnit, '(a)') 'Testing routine to calculate the muffin-tin gradient of a lattice-harmonic expanded quantity'
    write(logUnit, '(a)') '--------------------------------------------------------------------------------------------'
    clnu_atom = cmplx(1, 0) ! we need no transformation here, this was tested elsewhere
    ! generate testfunction r Y_00
    testFunc = 0
    do itype = 1, atomsT%ntype
      do imesh = 1, atomsT%jri(itype)
        testFunc(imesh, 0, itype) = atomsT%rmsh(imesh, itype) ! only l = 0, m = 0 component is not equals zero
      end do
    end do
    call calcGrFinLH(atomsT, latharT, clnu_atom, nmem_atom, mlh_atom, testFunc, gradrho0MT)

    ! remove Jacobi determinant which complicates the comparison with the analytical result
!    do itype = 1, atomsT%ntype
!      do imesh = 1, atomsT%jri(itype)
!        gradrho0MT(imesh, :, :, :) = gradrho0MT(imesh, :, :, :) / atomsT%rmsh(imesh, itype)**2
!      end do
!    end do

    ! compare to analytical result
    do idirec = 1, 3
      iatom = 1
      do itype = 1, atomsT%ntype
        do ieq = 1, atomsT%neq(itype)
          do oqn_l = 0, atomsT%lmax(itype)
            do mqn_m = -oqn_l, oqn_l
              lm = oqn_l * (oqn_l + 1) + 1 + mqn_m
              if ( (idirec == 1) .and. ( lm == 2 )) then
                if ( ( any ( ( real( gradRho0MT(:, lm, iatom, idirec)) - (1 / sqrt(6.)) ) > 1e-6 ) ) .or. ( any (aimag( gradRho0MT(:, lm, iatom, idirec)) > 1e-6 ) ) ) then
                  write (logUnit, *) 'Inconsistency at MT Derivative'
                  write (logUnit, *) 'at direction, atom index, l, m'
                  write (logUnit, *) idirec, iatom, oqn_l, mqn_m
                  passed = .false.
                end if
              else if ( (idirec == 1) .and. ( lm == 4 ) ) then
                if ( ( any (real( gradRho0MT(:, lm, iatom, idirec)) - (-1 / sqrt(6.)) > 1e-6 ) ) .or. ( any (aimag( gradRho0MT(:, lm, iatom, idirec)) > 1e-6 ) ) ) then
                  write (logUnit, *) 'Inconsistency at MT Derivative'
                  write (logUnit, *) 'at direction, atom index, l, m'
                  write (logUnit, *) idirec, iatom, oqn_l, mqn_m
                  passed = .false.
                end if
              !todo square root is made positive, probably the T matrix was corrected in the mean time commented out
              else if ( (idirec == 2) .and. ( lm == 2 ) ) then
                if ( (any (real( gradRho0MT(:, lm, iatom, idirec)) > 1e-6  ) ) .or. ( any (aimag( gradRho0MT(:, lm, iatom, idirec)) - (+1 / sqrt(6.)) > 1e-6 ) ) ) then
                  write (logUnit, *) 'Inconsistency at MT Derivative'
                  write (logUnit, *) 'at direction, atom index, l, m'
                  write (logUnit, *) idirec, iatom, oqn_l, mqn_m
                  write (logUnit, *) (abs(aimag( gradRho0MT(:, lm, iatom, idirec)) - (-1 / sqrt(6.))))
                  passed = .false.
                end if
              !todo square is made positive, probably the T matrix was corrected in the mean time, while this test was commented out
              else if ( (idirec == 2) .and. ( lm == 4 ) ) then
                if ( ( any (real( gradRho0MT(:, lm, iatom, idirec)) > 1e-6  ) ) .or. ( any (aimag( gradRho0MT(:, lm, iatom, idirec)) - (+1 / sqrt(6.)) > 1e-6 ) ) ) then
                  write (logUnit, *) 'Inconsistency at MT Derivative'
                  write (logUnit, *) 'at direction, atom index, l, m'
                  write (logUnit, *) idirec, iatom, oqn_l, mqn_m
                  passed = .false.
                end if
              else if ( (idirec == 3) .and. ( lm == 3 ) ) then
                if ( ( any (real( gradRho0MT(:, lm, iatom, idirec)) - 1 / sqrt(3.)  > 1e-6) ) .or. ( any (aimag( gradRho0MT(:, lm, iatom, idirec)) > 1e-6) ) ) then
                  write (logUnit, *) 'Inconsistency at MT Derivative'
                  write (logUnit, *) 'at direction, atom index, l, m'
                  write (logUnit, *) idirec, iatom, oqn_l, mqn_m
                  passed = .false.
                end if
              else if ( any( abs( gradRho0MT( :, lm, iatom, idirec) ) > 1e-6 ) ) then
                write (logUnit, *) 'Inconsistency at MT Derivative'
                write (logUnit, *) 'at direction, atom index, l, m'
                write (logUnit, *) idirec, iatom, oqn_l, mqn_m
                  passed = .false.
              end if
            end do
          end do
          iatom = iatom + 1
        end do
      end do
    end do

   ! write (*, *) 'l and m in lattice harmonics'
   ! iatom = 1
   ! do itype = 1, atomsT%ntype
   !   do ieq = 1, atomsT%neq(itype)
   !     write (*, *) 'Atom', iatom
   !     do lh = 0, latharT%nlh(atomsT%ntypsy(1))
   !       do mem = 1, nmem_atom(lh, 1)
   !         write (*, '(a,i2,a,i2,a,i3,a,i3)') 'lh=', lh, ' mem=', mem, ' l=', latharT%llh(lh, atomsT%ntypsy(itype)), ' m=', mlh_atom(mem, lh, iatom)
   !       end do
   !     end do
   !     iatom = iatom + 1
   !     write (*, *)
   !   end do
   ! end do

    ! This test is only working if we don't have any prefactors like Gaunt coefficients or other. For this test potential the terms with
    ! the density are neutralized as they are 1 all the time!
    !testFunc = 0
    !do itype = 1, atomsT%ntype
    !  do imesh = 1, atomsT%jri(itype)
    !    testFunc(imesh, :, itype) = atomsT%rmsh(imesh, itype) ! all components are equals r
    !  end do
    !end do

    !call calcGrFinLH(atomsT, latharT, clnu_atom, nmem_atom, mlh_atom, testFunc, gradrho0MT)

    !do idirec = 1, 3
    !  iatom = 1
    !  do itype = 1, atomsT%ntype
    !    do ieq = 1, atomsT%neq(itype)
    !      do oqn_l = 0, atomsT%lmax(itype)
    !        do mqn_m = -oqn_l, oqn_l
    !          lm = oqn_l * (oqn_l + 1) + 1 + mqn_m
    !          if ( any( real( gradRho0MT(:, lm, iatom, idirec) ) > 1e-6 ) .or. any( aimag ( gradRho0MT( :, lm, iatom, idirec) ) > 1e-6 )) then
    !              write (*, '(a,1x,i1,a,1x,i1, a,i3,a,i3)') 'Finite value at direction', idirec, ', atom index', iatom, ', l=', oqn_l, ', m=', mqn_m
    !              write (*, *) idirec, oqn_l, mqn_m, gradRho0MT(:, lm, iatom, idirec)
    !              read(*, *)
    !          end if
    !        end do
    !      end do
    !      iatom = iatom + 1
    !    end do
    !  end do
    !end do
    if (passed) then
      write (logUnit, '(a)')'                                                                                           |__ passed!'
    else
      write (logUnit, '(a)')'                                                                                           |__ failed!'
      call JuDFT_warn('Testing routine to calculate the muffin-tin gradient of a lattice-harmonic expanded quantity.', calledby='checkMTDerivative', &
        & hint='Check log output')
    end if

  end subroutine checkMTDerivative

  ! Tests the gaunt1 routine with manually calculated Gaunt coefficients.
  subroutine testGauntCoeff(atomsT, logUnit)

    use m_gaunt
    use m_jPConstants
    use m_juDFT_NOstopNO, only : juDFT_warn

    implicit none

    type(t_atoms), intent(in) :: atomsT
    integer, intent(in) :: logUnit
    real          :: tempGaunt
    real          :: tempFG
    integer       :: oqn_l
    integer       :: mqn_m
    logical       :: passed = .true.

    write(logUnit, '(a)') 'Test for Gaunt coefficients routines'
    write(logUnit, '(a)') '------------------------------------'

    do oqn_l = 0, 30
      do mqn_m = -oqn_l, oqn_l

        tempGaunt = - sqrt(4 * pi / 3) * gaunt1(oqn_l + 1, oqn_l, 1, mqn_m + 1, mqn_m, 1, 31)
        tempFG = - sqrt( (float( ( oqn_l + mqn_m + 1 ) * (oqn_l + mqn_m + 2 ) )) / (2 * (2 * oqn_l + 1 ) * ( 2 * oqn_l + 3 ) ) )
        if ((tempGaunt - tempFG) > 1e-9) then
          write(logUnit, '(a)')'Inconsistency in Gaunt function routine test 1'
          write (logUnit, *) oqn_l, mqn_m, tempFG, tempGaunt
          passed = .false.
        end if

        if ( abs(mqn_m + 1) <= oqn_l - 1) then
          tempGaunt = - sqrt(4 * pi / 3) * gaunt1(oqn_l - 1, 1, oqn_l , mqn_m + 1, 1, mqn_m, 31)
          tempFG = sqrt( (float( ( oqn_l - mqn_m ) * (oqn_l - mqn_m - 1) ) ) / (2 * (2 * oqn_l - 1 ) * (2 * oqn_l + 1) ) )
          if ((tempGaunt - tempFG) > 1e-9) then
            write(logUnit, '(a)')'Inconsistency in Gaunt function routine test 2'
            write (logUnit, *) oqn_l, mqn_m, tempFG, tempGaunt
            passed = .false.
          end if
        end if

        tempGaunt = sqrt(4 * pi / 3) * gaunt1(oqn_l + 1, oqn_l, 1, mqn_m, mqn_m, 0, 31)
        tempFG = sqrt( float(( oqn_l - mqn_m + 1 ) * ( oqn_l + mqn_m + 1 )) / ((2 * oqn_l + 1) * (2 * oqn_l + 3)))
        if ((tempGaunt - tempFG) > 1e-9) then
          write(logUnit, '(a)')'Inconsistency in Gaunt function routine test 3'
          write (logUnit, *) oqn_l, mqn_m, tempFG, tempGaunt
          passed = .false.
        end if

        if ( abs(mqn_m) <= oqn_l - 1 ) then
          tempGaunt = sqrt(4 * pi / 3) * gaunt1(oqn_l - 1, oqn_l, 1, mqn_m, mqn_m, 0, 31)
          tempFG = sqrt(float(( oqn_l - mqn_m ) * ( oqn_l + mqn_m )) / (( 2 * oqn_l - 1) * (2 * oqn_l + 1)))
          if ((tempGaunt - tempFG) > 1e-9) then
            write(logUnit, '(a)')'Inconsistency in Gaunt function routine test 4'
            write (logUnit, *) oqn_l, mqn_m, tempFG, tempGaunt
            passed = .false.
          end if
        end if

        tempGaunt = -sqrt(4 * pi / 3) * gaunt1(oqn_l + 1, oqn_l, 1, mqn_m - 1, mqn_m, -1, 31)
        tempFG = - sqrt( (float( ( oqn_l - mqn_m + 1 ) * (oqn_l - mqn_m + 2 ) )) / (2 * (2 * oqn_l + 1 ) * ( 2 * oqn_l + 3 ) ) )
        if ((tempGaunt - tempFG) > 1e-9) then
          write(logUnit, '(a)')'Inconsistency in Gaunt function routine test 5'
          write (logUnit, *) oqn_l, mqn_m, tempFG, tempGaunt
          passed = .false.
        end if

        if ( abs( mqn_m - 1 ) <= oqn_l - 1 ) then
          tempGaunt = -sqrt(4 * pi / 3) * gaunt1(oqn_l - 1, oqn_l, 1, mqn_m - 1, mqn_m, -1, 31)
          tempFG = sqrt( (float( ( oqn_l + mqn_m ) * (oqn_l + mqn_m - 1) ) ) / (2 * (2 * oqn_l - 1 ) * (2 * oqn_l + 1) ) )
          if ((tempGaunt - tempFG) > 1e-9) then
            write(logUnit, '(a)')'Inconsistency in Gaunt function routine test 6'
            write (logUnit, *) oqn_l, mqn_m, tempFG, tempGaunt
            passed = .false.
          end if
        end if
      end do
    end do

    if (passed) then
      write(logUnit, '(a)') '                                   |__ passed!'
    else
      write(logUnit, '(a)') '                                   |__ failed!'
      call JuDFT_warn('Test for Gaunt coefficients routines', calledby='testGauntCoeff', &
        & hint='Check log output and debug')
    end if

  end subroutine testGauntCoeff

  ! This test checks whether the Vext1 for a q is the complex conjugate of the Vext1 at the partner q ( 0 0 0.25 ) <-> ( 0 0 0.75 )
  ! Note: The G-sets at +q and -q are related, so are the Gset at -q + Gf
!24.01.2021
!If I want to access 0.75 from 0.25, I have to go to -0.25 which gives a -G and then go + G_t to 0.75 so that the Gs have to be subtracted by -Gt. If I want to make the new G-vectors equal, I first go back to -0.25 and then to 0.25 so that it gives me a -(Gz + 1) = -Gz - 1 in analogy to the previous test
  subroutine TestQLatticePeriodicity( atoms, stars, cell, lathar, dimens, sym, input, qpts, ngdp, memd_atom, logUnit, rho0IR, &
                                                                                     & rho0MT, mlh_atom, nmem_atom, clnu_atom, gdp )

    use m_types
    use m_jpVeff1, only : GenVeff1
    use m_jpPotDensHelper, only : WarpIRPot, genPertPotDensGvecs, genPotDensGvecs
    use m_juDFT_NOstopNO, only : juDFT_warn

    implicit none

    ! Type parameters
    type(t_atoms),                  intent(in)  :: atoms
    type(t_stars),                  intent(in)  :: stars
    type(t_cell),                   intent(in)  :: cell
    type(t_sphhar),                 intent(in)  :: lathar
    type(t_dimension),              intent(in)  :: dimens
    type(t_sym),                    intent(in)  :: sym
    type(t_input),                  intent(in)  :: input
    type(t_kpts),                   intent(in)  :: qpts

    ! Scalar parameters
    integer,                        intent(in)  :: ngdp
    integer,                        intent(in)  :: memd_atom
    integer,                        intent(in)  :: logUnit

    ! Array parameters
    complex,                        intent(in)  :: rho0IR(:, :)
    real,                           intent(in)  :: rho0MT(:, 0:, :, :)
    integer,                        intent(in)  :: mlh_atom(:, 0:, :)
    integer,                        intent(in)  :: nmem_atom(0:, :)
    complex,                        intent(in)  :: clnu_atom(:, 0:, :)
    integer,                        intent(in)  :: gdp(:, :)

    ! Scalar variables
    logical                                     :: harSw
    logical                                     :: extSw
    logical                                     :: xcSw
    logical                                     :: vExtFull
    integer                                     :: iDatom
    integer                                     :: iDtype
    integer                                     :: iqpt
    integer                                     :: iGvar
    integer                                     :: iG
    integer                                     :: idir
    integer                                     :: ngdp2
    integer                                     :: ngdp2km
    integer                                     :: ngpqdp
    integer                                     :: ngpqdp3
    integer                                     :: ngpqdp2km
    logical                                     :: passed
    integer                                     :: ii
    integer                                     :: jj
    integer                                     :: kk
    logical                                     :: vHarNum


    ! Array variables
    real                                        :: qpoint(3)
    integer                                     :: gdp2iLim(2, 3)
    integer                                     :: gShift(3)
    integer                                     :: gpqdp2iLim(2, 3)
    complex,           allocatable              :: rho1PW(:, :)
    complex,           allocatable              :: rho1MT(:, :, :, :)
    complex,           allocatable              :: grRho0MT(:, :, :, :)
    complex,           allocatable              :: vEff1IR(:, :)
    complex,           allocatable              :: vEff1IRRef(:, :)
    complex,           allocatable              :: vEff1MT(:, :, :, :)
    complex,           allocatable              :: vxc1IRKern(:)
    complex,           allocatable              :: ylm(:, :)
    real,              allocatable              :: dKernMTGPts(:, :, :)
    real,              allocatable              :: gWghts(:) ! gaussian weights belonging to gausPts
    integer,           allocatable              :: gdpIndex(:, :, :)
    complex,           allocatable              :: w_vEff1IR(:, :)
    complex,           allocatable              :: w_vEff1IRRef(:, :)
    integer,           allocatable              :: gdp2(:, :)
    integer,           allocatable              :: gdp2Ind(:, :, :)
    integer,           allocatable              :: gpqdp(:, :)
    integer,           allocatable              :: gpqdp2(:, :)
    integer,           allocatable              :: gpqdp3(:, :)
    integer,           allocatable              :: gpqdp2Ind(:, :, :)
    complex,           allocatable              :: rho0IRDummy(:, :)
    complex,           allocatable              :: rho0MTDummy(:, :, :, :)

    write(logUnit, '(a)')
    write(logUnit, '(a)') 'Performing testVeff1IRqLatPeriod'
    write(logUnit, '(a)') '--------------------------------'

    if (.false.) then

      ! Understand how the G-set are shifted due to a q-vector. Deativated by default because no benefit for the test. Left here
      ! for debugging. This was not reviewed again when test was added to test suite.
      !-------------------------------------------------------------------------------------------------------------------------

      ! Generate conventional G-set
      call genPotDensGvecs(stars, cell, input, ngdp2, ngdp2km, gdp2, gdp2Ind, gdp2iLim, .true. )

      write(2169, '(i8)') ngdp2
      write(2170, '(i8)') ngdp2km
      write(2171, '(6(i8))') gdp2iLim(1, 1), gdp2iLim(2, 1), gdp2iLim(1, 2), gdp2iLim(2, 2), gdp2iLim(1, 3), gdp2iLim(2, 3)

      do iG = 1, ngdp2
        write(2172, '(3(i8))') gdp2(1, iG), gdp2(2, iG), gdp2(3, iG)
      end do

      do ii = gdp2iLim(1, 1), gdp2iLim(2, 1)
        do jj = gdp2iLim(1, 2), gdp2iLim(2, 2)
          do kk = gdp2iLim(1, 3), gdp2iLim(2, 3)
            write(2173, '(5(i8))') ii, jj, kk, gdp2Ind(ii, jj, kk)
          end do ! kk
        end do ! jj
      end do ! ii

      ! Shift G-set with q -> -q
      call genPertPotDensGvecs( stars, cell, input, ngpqdp, ngpqdp2km, [0., 0., -0.25], gpqdp, gpqdp2Ind, gpqdp2iLim )
      write(2174, '(i8)') ngpqdp
      write(2175, '(i8)') ngpqdp2km
      write(2176, '(6(i8))') gpqdp2iLim(1, 1), gpqdp2iLim(2, 1), gpqdp2iLim(1, 2), gpqdp2iLim(2, 2), gpqdp2iLim(1, 3), gpqdp2iLim(2, 3)

      do iG = 1, ngpqdp
        write(2177, '(3(i8))') gpqdp(1, iG), gpqdp(2, iG), gpqdp(3, iG)
      end do

      do ii = gpqdp2iLim(1, 1), gpqdp2iLim(2, 1)
        do jj = gpqdp2iLim(1, 2), gpqdp2iLim(2, 2)
          do kk = gpqdp2iLim(1, 3), gpqdp2iLim(2, 3)
            write(2178, '(5(i8))') ii, jj, kk, gpqdp2Ind(ii, jj, kk)
          end do ! kk
        end do ! jj
      end do ! ii

      deallocate(gpqdp, gpqdp2Ind)
      gpqdp2iLim = cmplx(0., 0.)
      ngpqdp = 0
      ngpqdp2km = 0
      call genPertPotDensGvecs( stars, cell, input, ngpqdp, ngpqdp2km, [0., 0., 0.75], gpqdp, gpqdp2Ind, gpqdp2iLim )
      write(2179, '(i8)') ngpqdp
      write(2180, '(i8)') ngpqdp2km
      write(2181, '(6(i8))') gpqdp2iLim(1, 1), gpqdp2iLim(2, 1), gpqdp2iLim(1, 2), gpqdp2iLim(2, 2), gpqdp2iLim(1, 3), gpqdp2iLim(2, 3)

      do iG = 1, ngpqdp
        write(2182, '(3(i8))') gpqdp(1, iG), gpqdp(2, iG), gpqdp(3, iG)
      end do

      do ii = gpqdp2iLim(1, 1), gpqdp2iLim(2, 1)
        do jj = gpqdp2iLim(1, 2), gpqdp2iLim(2, 2)
          do kk = gpqdp2iLim(1, 3), gpqdp2iLim(2, 3)
            write(2183, '(5(i8))') ii, jj, kk, gpqdp2Ind(ii, jj, kk)
          end do ! kk
        end do ! jj
      end do ! ii

    end if

    ! Calculate G-set for q and Veff1 for same q-point
    harSw = .false.
    extSw = .true.
    xcSw = .false.
    vExtFull = .false.
    vHarNum = .false.

    iDatom = 1
    iDtype = 1
    iqpt = 1

    passed = .true.

    allocate( rho1PW( ngdp, 3 ),                                                   &
             &rho1MT( atoms%jmtd, (atoms%lmaxd + 1)**2, atoms%nat, 3 ), &
             &grRho0MT( atoms%jmtd, ( atoms%lmaxd + 1 )**2, atoms%nat, 3 ),        &
             &vxc1IRKern(ngdp),                                                    &
             &ylm(dimens%nspd, ( atoms%lmaxd + 1)**2 ),                            &
             &dKernMTGpts(dimens%nspd, atoms%jmtd, atoms%nat),                     &
             &gWghts(dimens%nspd) )

     allocate(  rho0IRDummy(ngdp, 1) )
     rho0IRDummy(:, :) = cmplx(0., 0.)

     allocate( rho0MTDummy( atoms%jmtd, (atoms%lmaxd + 1)**2, atoms%nat, 1) )
     rho0MTDummy(:, :, :, :) = cmplx(0., 0.)

    rho1PW(:, :)          = cmplx(0., 0.)
    rho1MT( :, :, :, :)   = cmplx(0., 0.)
    grRho0MT(:, :, :, :)  = cmplx(0., 0.)
    vxc1IRKern(:)         = cmplx(0., 0.)
    ylm(:, :)             = cmplx(0., 0.)
    dKernMTGpts(:, :, :)  = 0.
    gWghts(:)             = 0.
    do iqpt = 1, qpts%nkpt
    do idir = 1, 3
      if (qpts%bk(idir, iqpt) > 1e-12) then
        gShift(idir) = 1
      else
        gshift(idir) = 0
      end if
    end do ! idir
    qpoint(1:3) = qpts%bk(1:3, iqpt)
    call genPertPotDensGvecs( stars, cell, input, ngpqdp, ngpqdp2km, qpoint, gpqdp, gpqdp2Ind, gpqdp2iLim )
    call GenVeff1( stars, cell, atoms, dimens, harSw, extSw, xcSw, vExtFull, ngdp, qpoint, rho0IRDummy, rho0MTDummy, rho1PW, rho1MT(:, :, :, :), grRho0MT, &
      & gdp, vEff1IR, vEff1MT, vxc1IRKern, ylm, dKernMTGPts, gWghts, iDatom, iDtype, iqpt, ngpqdp, gpqdp, vHarNum ) ! add spin

    allocate(w_vEff1IR(ngpqdp, 3))
    do idir = 1, 3
      call warpIRPot(stars, ngpqdp, idir, gpqdp, vEff1IR, w_vEff1IR(:, idir))
    end do

    if (.false.) then
      do idir = 1, 3
        do iG = 1, ngpqdp
          write(2184, '(5i8,2f15.8)') idir, iG, gpqdp(1, iG), gpqdp(2, iG), gpqdp(3, iG), vEff1IR(iG, idir)
          write(2185, '(5i8,2f15.8)') idir, iG, gpqdp(1, iG), gpqdp(2, iG), gpqdp(3, iG), w_vEff1IR(iG, idir)
        end do
      end do
    end if
    deallocate(vEff1MT, gpqdp2Ind)

    ! Gset for -q and Veff1
    qpoint(1:3) = -qpts%bk(1:3, iqpt)
    ngpqdp = 0

    call genPertPotDensGvecs( stars, cell, input, ngpqdp, ngpqdp2km, qpoint, gpqdp2, gpqdp2Ind, gpqdp2iLim )
    allocate(gdpIndex(minval(gpqdp2(1, :)) : maxval(gpqdp2(1, :)), minval(gpqdp2(2, :)) : maxval(gpqdp2(2, :)), minval(gpqdp2(3, :))-1 : maxval(gpqdp2(3, :))))
    gdpIndex(:, :, :) = 0

    do iG = 1, ngpqdp
      gdpIndex(gpqdp2(1, iG), gpqdp2(2, iG), gpqdp2(3, iG)) = iG
    end do ! iG
    call GenVeff1( stars, cell, atoms, dimens, harSw, extSw, xcSw, vExtFull, ngdp, qpoint, rho0IRDummy, rho0MTDummy, rho1PW, rho1MT(:, :, :, :), grRho0MT, &
      & gdp, vEff1IRRef, vEff1MT, vxc1IRKern, ylm, dKernMTGPts, gWghts, iDatom, iDtype, iqpt, ngpqdp, gpqdp2, vHarNum ) ! add spin

    allocate(w_vEff1IRREf(ngpqdp, 3))
    do idir = 1, 3
      call warpIRPot(stars, ngpqdp, idir, gpqdp2, vEff1IRRef, w_vEff1IRRef(:, idir))
    end do

    ! Note: We find the values of +q for a certain G-vector at -q at the negative G-vector. Especially in Vext we see that a minus
    !       is generated in front of the (G + q) vector. This leads to the complex conjugate which is also should be because changing
    !       from +q to minus q gives a complex conjugate.
    do idir = 1, 3
      do iG = 1, ngpqdp
        iGvar = gdpIndex(-gpqdp(1, iG), -gpqdp(2, iG), -gpqdp(3, iG))
        if ( (abs(conjg(vEff1IRRef(iGvar, idir)) - vEff1IR(iG, idir)) > 1e-8) .or. (abs(conjg(w_vEff1IRRef(iGvar, idir)) - w_vEff1IR(iG, idir)) > 1e-8) ) then
          passed = .false.
        end if

        if (.false.) then
          write(2186, '(5i8,2f15.8)') idir, iG, gpqdp2(1, iGvar), gpqdp2(2, iGvar), gpqdp2(3, iGvar), vEff1IRRef(iGvar, idir)
          write(2187, '(5i8,2f15.8)') idir, iG, -gpqdp2(1, iGvar), -gpqdp2(2, iGvar), -gpqdp2(3, iGvar), conjg(vEff1IRRef(iGvar, idir))
          write(2188, '(5i8,2f15.8)') idir, iG, -gpqdp2(1, iGvar), -gpqdp2(2, iGvar), -gpqdp2(3, iGvar), conjg(w_vEff1IRRef(iGvar, idir))
        end if
      end do
    end do
    deallocate(vEff1IRRef, vEff1MT, w_vEff1IRRef, gpqdp2Ind, gdpIndex)

    ! Gset for -q projected back to the shifted first Brillouin zone (0 to 1)
    qpoint(1:3) = -qpts%bk(1:3, iqpt) + gShift(1:3)

    call genPertPotDensGvecs( stars, cell, input, ngpqdp3, ngpqdp2km, qpoint, gpqdp3, gpqdp2Ind, gpqdp2iLim )
    allocate(gdpIndex(minval(gpqdp3(1, :)) : maxval(gpqdp3(1, :)), minval(gpqdp3(2, :)) : maxval(gpqdp3(2, :)), minval(gpqdp3(3, :))-1 : maxval(gpqdp3(3, :))))
    gdpIndex(:, :, :) = -1

    do iG = 1, ngpqdp3
      gdpIndex(gpqdp3(1, iG), gpqdp3(2, iG), gpqdp3(3, iG)) = iG
    end do ! iG

    call GenVeff1( stars, cell, atoms, dimens, harSw, extSw, xcSw, vExtFull, ngdp, qpoint, rho0IRDummy, rho0MTDummy, rho1PW, rho1MT(:, :, :, :), grRho0MT, &
      & gdp, vEff1IRRef, vEff1MT, vxc1IRKern, ylm, dKernMTGPts, gWghts, iDatom, iDtype, iqpt, ngpqdp3, gpqdp3, vHarNum ) ! add spin

    allocate(w_vEff1IRRef(ngpqdp3, 3))
    do idir = 1, 3
      call warpIRPot(stars, ngpqdp3, idir, gpqdp3, vEff1IRRef, w_vEff1IRRef(:, idir))
    end do

    if (ngpqdp3 /= ngpqdp) NOstopNO'Gsets should have same size at symmetry related q-points'
    do idir = 1, 3
      do iG = 1, ngpqdp
        iGvar = gdpIndex(-(gpqdp(1, iG) + gShift(1)), -(gpqdp(2, iG) + gShift(2)), -(gpqdp(3, iG) + gShift(3)))
        if (iGvar == -1) then
          write(*, *) 'Error in Gset at -q - Gf'
          write(*, *) iG
          NOstopNO'index out of scope'
        end if
        ! Note: if we start from q, we find the partner in the first Brillouin zone if we additionaly shift -q by the shiftvector Gf
        if ( (abs(conjg(vEff1IRRef(iGvar, idir)) - vEff1IR(iG, idir)) > 1e-8) .or. (abs(conjg(w_vEff1IRRef(iGvar, idir)) - w_vEff1IR(iG, idir)) > 1e-8) ) then
          passed = .false.
        end if
        if (.false.) then
          ! If I want to access 0.75 from 0.25, I have to go to -0.25 which gives a -G and then go + G_T to 0.75 so that the Gs
          ! have to be subtracted by - G_T. If I want to make the new G-Vectors equal, I first go back to -0.25 and than to 0.25
          ! so that it gives me a -(Gz + 1) = - Gz - 1 in analogy to the previous test.
          write(2189, '(5i8,2f15.8)') idir, iG, gpqdp3(1, iGvar), gpqdp3(2, iGvar), gpqdp3(3, iGvar), vEff1IRRef(iGvar, idir)
          write(2190, '(5i8,2f15.8)') idir, iG, -(gpqdp3(1, iGvar) + gShift(1)), -(gpqdp3(2, iGvar) + gShift(2)), -(gpqdp3(3, iGvar) +gShift(3)), conjg(vEff1IRRef(iGvar, idir))
          write(2191, '(5i8,2f15.8)') idir, iG, -(gpqdp3(1, iGvar) + gShift(1)), -(gpqdp3(2, iGvar) + gShift(2)), -(gpqdp3(3, iGvar) +gShift(3)), conjg(w_vEff1IRRef(iGvar, idir))
        end if
      end do ! iG
    end do ! idir
    deallocate(vEff1IRRef, vEff1MT, w_vEff1IRRef, w_vEff1IR, vEff1IR, gdpIndex, gpqdp2Ind)
    end do !iqpt


    if (passed) then
      write(logUnit, '(a)') '                               |__ passed!'
    else
      write(logUnit, '(a)') '                               |__ failed!'
      call JuDFT_warn('Test for q-lattice integrity of Vext1 failed.', calledby='testGauntCoeff', &
        & hint='Check log output and debug')
    end if

    write(logUnit, '(a)')


  end subroutine TestQLatticePeriodicity

end module m_jpTestPotential
