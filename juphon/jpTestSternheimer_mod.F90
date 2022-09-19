module m_jpTestSternheimer

  implicit none

  contains

  !--------------------------------------------------------------------------------------------------------------------------------------
  !> @author
  !> Christian-Roman Gerhorst, Gregor Michalicek
  !>
  !> @brief
  !> Main routine operating the tests for the Sternheimer equation terms and their ingredients.
  !>
  !> @details
  !> This routine calls testings routines from the potential test module as well as from the Sternheimer test module.
  !>
  !>
  !> @note
  !> Additional information and formulas pointing out the routines of this module can be found within this
  !> <a href='jpTestSternheimer.pdf'>document</a>.
  !--------------------------------------------------------------------------------------------------------------------------------------
  subroutine TestSternheimerSetup( atoms, cell, sym, stars, lathar, enpara, kpts, qpts, input, dimens, ud, V0Fleur, results, iqpt, &
      & ngdp, ngdp2km, paPoX, paPoY, paPoZ, harSw, extSw, xcSw, testCompareGrVeff0FleurSw, testVeff1Sw, testUnfoldStarsSw,         &
      & testRadDerivativeSw, testGauntCoeffSw, testGradLhExpandFuncSw, testContGrVeff0Sw, testWarpingSw, testSternhHSMEtestSw,     &
      & testSternhSchroedPertTheoSw, testz1Phi0ContSw, testRho1BasCorrSw, testPlotRho03Dsw, testGradRho0PathSw, testRho1IRsw,      &
      & testRho1MTsw, logUnit, gdp, mlh_atom, nmem_atom, clnu_atom, rho0IR, rho0MT, ne, nobd, eig, rbas1, rbas2, uuilon, duilon,   &
      & ulouilopn, GbasVec, ilst, nv, z, kveclo, iloTable, nRadFun, El, mapKpq2K, ilo2p, memd_atom, noPtsCon, gdp2Ind, gdp2iLim,   &
      & kpq2kPrVec, testrho0ContSw, testCompareSurfInt, testSplitMTSurfIntSterh,                &
      & testVeff1IRMESternh, testEps1q0, testVeff1IRMatqBackFold, testVeff1IRqLatPeriod, vEff0MTsh, vEff0IRpwUw )

#include "cppmacro.h"
    use m_jpTestPotential
    use m_types
    use jpTest1stVarDens
    use m_jpPotDensHelper, only : convertStar2G

    implicit none

    ! Type parameters
    type(t_atoms),     intent(in) :: atoms
    type(t_cell),      intent(in) :: cell
    type(t_sym),       intent(in) :: sym
    type(t_stars),     intent(in) :: stars
    type(t_sphhar),    intent(in) :: lathar
    type(t_enpara),    intent(in) :: enpara
    type(t_kpts),      intent(in) :: kpts
    type(t_kpts),      intent(in) :: qpts
    type(t_input),     intent(in) :: input
    type(t_dimension), intent(in) :: dimens
    type(t_usdus),     intent(in) :: ud
    type(t_potential), intent(in) :: V0Fleur
    type(t_results),   intent(in) :: results


    ! Scalar parameters
    integer,           intent(in) :: iqpt
    integer,           intent(in) :: ngdp
    integer,           intent(in) :: ngdp2km
    real,              intent(in) :: paPoX
    real,              intent(in) :: paPoY
    real,              intent(in) :: paPoZ
    logical,           intent(in) :: harSw
    logical,           intent(in) :: extSw
    logical,           intent(in) :: xcSw
    logical,           intent(in) :: testCompareGrVeff0FleurSw
    logical,           intent(in) :: testVeff1Sw
    logical,           intent(in) :: testUnfoldStarsSw
    logical,           intent(in) :: testrho0ContSw
    logical,           intent(in) :: testRadDerivativeSw
    logical,           intent(in) :: testGauntCoeffSw
    logical,           intent(in) :: testGradLhExpandFuncSw
    logical,           intent(in) :: testContGrVeff0Sw
    logical,           intent(in) :: testWarpingSw 
    logical,           intent(in) :: testSternhHSMEtestSw 
    logical,           intent(in) :: testSternhSchroedPertTheoSw
    logical,           intent(in) :: testz1Phi0ContSw
    logical,           intent(in) :: testRho1BasCorrSw
    logical,           intent(in) :: testPlotRho03Dsw
    logical,           intent(in) :: testGradRho0PathSw
    logical,           intent(in) :: testRho1IRSw
    logical,           intent(in) :: testRho1MTSw
    integer,           intent(in) :: logUnit
    integer,           intent(in) :: memd_atom
    integer,           intent(in) :: noPtsCon
    logical,           intent(in) :: testCompareSurfInt
    logical,           intent(in) :: testSplitMTSurfIntSterh
    logical,           intent(in) :: testVeff1IRMESternh
    logical,           intent(in) :: testEps1q0
    logical,           intent(in) :: testVeff1IRMatqBackFold
    logical,           intent(in) :: testVeff1IRqLatPeriod

    ! Array parameters
    integer,           intent(in) :: gdp(:, :)
    integer,           intent(in) :: mlh_atom(:,0:,:)
    integer,           intent(in) :: nmem_atom(0:, :)
    complex,           intent(in) :: clnu_atom(:,0:,:)
    complex,           intent(in) :: rho0IR(:,:) !Interstitial density
    real,              intent(in) :: rho0MT(:,0:,:,:)  ! MT sphere density
    integer,           intent(in) :: ne(:)
    integer,           intent(in) :: nobd(:, :)
    real,              intent(in) :: eig(:, :, :)
    real,              intent(in) :: rbas1(:, :, 0:, :, :)
    real,              intent(in) :: rbas2(:, :, 0:, :, :)
    real,              intent(in) :: uuilon(:, :)
    real,              intent(in) :: duilon(:, :)
    real,              intent(in) :: ulouilopn(:, :, :)
    integer,           intent(in) :: GbasVec(:, :)
    integer,           intent(in) :: ilst(:, :, :)
    integer,           intent(in) :: nv(:, :)
    MCOMPLEX,          intent(in) :: z(:,:,:,:)
    integer,           intent(in) :: kveclo(:,:)
    integer,           intent(in) :: iloTable(:, 0:, :)
    integer,           intent(in) :: nRadFun(0:, :)
    real,              intent(in) :: El(:, 0:, :, :)
    integer,           intent(in) :: mapKpq2K(:, :) ! todo possible error here
    integer,           intent(in) :: ilo2p(:, :)
    integer,           intent(in) :: gdp2Ind( :, :, : )
    integer,           intent(in) :: gdp2iLim(2, 3)
    integer,           intent(in) :: kpq2kPrVec(:, :, :)
    complex,           intent(in) :: vEff0MTsh(:, :, :, :)
    complex,           intent(in) :: vEff0IRpwUw(:, :)

    complex, allocatable :: rho0IRPW(:, :)
    complex, allocatable :: rho0MTSH(:, :, :, :)
    complex, allocatable :: surfIntMT(:, :, :)
    complex, allocatable :: surfIntMTSafe(:, :, :)
    integer :: iatom, itype, ieqat, ptsym, ilh, oqn_l, imem, mqn_m, lm, imesh
    ! todo put them into Juphon input file
    !todo maybe type for the testvariables or for the whole input respectively


    call Test1stOrdPotentials( atoms, stars, lathar, sym, cell, dimens, input, V0Fleur, qpts, ngdp, paPoX, paPoY, paPoZ, harSw, extSw,   &
      & xcSw, logUnit, testCompareGrVeff0FleurSw, testVeff1Sw, testUnfoldStarsSw, testRadDerivativeSw, testGauntCoeffSw,           &
      & testGradLhExpandFuncSw, testContGrVeff0Sw, gdp, mlh_atom, nmem_atom, clnu_atom, rho0IR, rho0MT, memd_atom, noPtsCon, testVeff1IRqLatPeriod )

    call Test1stVarDens( atoms, sym, lathar, cell, input, results, stars, kpts, dimens, qpts, ud, ngdp2km, ngdp, logUnit,          &
      & testz1Phi0ContSw, testRho1IRsw, testRho1MTSw, testRho1BasCorrSw, testPlotRho03Dsw, testGradRho0PathSw, noPtsCon, iqpt,     &
      & paPoX, paPoY, paPoZ, harSw, extSw, xcSw, gdp, GbasVec, z, gdp2Ind, rho0IR, ne, ilst, nv, kveclo, rbas1, rbas2, ilo2p, nobd,&
      & mlh_atom, nmem_atom, clnu_atom, rho0MT, nRadFun, iloTable, kpq2kPrVec, mapKpq2K, gdp2iLim )

    ! Tests whether the Sternheimer converges immediatelly with the solution z1 = -i(k + G) z0
    if (.false.) then
      call testInstq0Conv( atoms, stars, cell, kpts, qpts, dimens, input, results, lathar, sym, enpara, ud, V0Fleur, ngdp, nobd, kveclo,&
        & memd_atom, rbas1, rbas2, uuilon, duilon, ilo2p, ulouilopn, z, nv, GbasVec, ilst, gdp2iLim, kpq2kPrVec, mapKpq2K,  &
        & gdp2Ind, clnu_atom, mlh_atom, nmem_atom, gdp, rho0IR, rho0MT, ne, eig, El, nRadFun, iloTable, noPtsCon, testrho0ContSw, logUnit, vEff0MTsh, vEff0IRpwUw)
    end if

    if (.false.) then
      ! Compares z1 in band G-vector represenattion to analytical solution for q = 0
      call Testz1ikGz0q0Diff( atoms, dimens, kpts, cell, nv, nobd, GbasVec, ilst, kveclo, z)
    end if

    if (.false.) then
      ! Compares z1 in band band represenattion to analytical solution for q = 0
      call TestSternhBandBandAnalyt( atoms, dimens, cell, kpts, nv, ne, Gbasvec, nobd, ilst, eig, z, kveclo )
    end if

    ! Tests whether the backfolding within Psi Veff1 Psi does function correctly
    if (.false.) then
      write(*, '(2x,a)') 'Performing TestIRbackFolding...'
      call testIRbackFolding(atoms, cell, lathar, input, stars, dimens, V0Fleur, kpts, sym, ngdp, ne, nv, nobd, GbasVec, ilst, gdp, kpq2kPrVec, mapKpq2K, z, rho0IR, rho0MT, nmem_atom, clnu_atom, mlh_atom, memd_atom)
    end if

    if ( testWarpingSw .or. testSternhHSMEtestSw  .or. testSternhSchroedPertTheoSw .or. testCompareSurfInt .or. testSplitMTSurfIntSterh .or. testVeff1IRMESternh .or. testEps1q0 ) then
      write(*, *)
      write(*, '(a)') 'Initiating Sternheimer test(s)...'
      write(*, '(a)') '---------------------------------'
      write(logUnit, *)
    write ( logUnit, * )
      write(logUnit, '(a)') 'Sternheimer test(s)'
      write(logUnit, '(a)') '*******************'
      write(logUnit, *)
    else
      write ( logUnit, * )
      write ( logUnit, * )
      write(*, '(a)') '----------------------------'
      write(*, '(a)') 'DISABLED Sternheimer test(s)!'
      write(*, '(a)') '----------------------------'
      write(logUnit, '(a)') 'DISABLED Sternheimer tests!'
      write(logUnit, '(a)') '**************************'
      return
    end if

    if ( testWarpingSw ) then
      write(*, '(2x,a)') 'Performing TestWarping...'
      ! Tests the warping routine and external parameter settings for them by warping a potential already warped by FLEUR
      call TestWarping(stars, ngdp, logUnit, gdp)
    else
      write(*, '(2x,a)') 'DISABLED TestWarping!'
      write ( logUnit, '(a)' ) 'Test of Sternheimer warping routine'
      write ( logUnit, '(a)' ) '-----------------------------------'
      write ( logUnit, '(a)' ) '                                  |_DISABLED'
    end if

    if ( testSternhHSMEtestSw ) then
      write(*, '(2x,a)') 'Performing TestHSME...'
      ! compares the setuped matrix element <Psi0|H|Psi0>_MT and <Psi0|Psi0>_MT with the Fleur result for them
      call TestHSME( atoms, sym, cell, lathar, dimens, enpara, ud, input, kpts, V0Fleur, logUnit, rbas1, rbas2, z, uuilon, duilon, &
        & ulouilopn, GbasVec, ilst, nv, ne, nobd, El, iloTable, nRadFun, kveclo, ilo2p, clnu_atom, nmem_atom, mlh_atom, vEff0MTsh )
    else
      write(*, '(2x,a)') 'DISABLED TestHSME!'
      write (logUnit, '(a)') 'Test comparison of muffin-tin Hamilton matrix in Sternheimer equation and Fleur with 1e-7 accuracy!'
      write (logUnit, '(a)') '---------------------------------------------------------------------------------------------------'
      write (logUnit, '(a)') '                                                                                                  |__&
                                                                                                                            &DISABLED'
      write (logUnit, '(a)') 'Test comparison of muffin-tin overlap matrix in Sternheimer equation and Fleur with 1e-7 accuracy!'
      write (logUnit, '(a)') '--------------------------------------------------------------------------------------------------'
      write (logUnit, '(a)') '                                                                                                 |__D&
                                                                                                                             &ISABLED'
    end if

    if ( testSternhSchroedPertTheoSw ) then
      write(*, '(2x,a)') 'Performing TestPotMEnEnergDiff...'
!!!!!!    !todo use clnu_atom not clnu within routines!!!!!
      ! sets a constant potential only taking the Hellmann-Feynman terms within Sternheimer equation into account. For a constant
      ! potential, the wavefunction should not be changed
      call TestPotMEnEnergDiff( atoms, cell, sym, stars, lathar, enpara, kpts, qpts, ud, input, dimens, ngdp, logUnit, ne, nobd,   &
        & eig, rbas1, rbas2, uuilon, duilon, ulouilopn, GbasVec, ilst, nv, z, kveclo, gdp, iloTable, nRadFun, El, mapKpq2K,        &
        & V0Fleur, nmem_atom, clnu_atom, mlh_atom, ilo2p, kpq2kPrVec )
    else
      write(logUnit, '(a)') 'Test of potential matrix element left side in Sternheimer equation'
      write(logUnit, '(a)') '------------------------------------------------------------------'
      write (logUnit, '(a)') '                                                                |__DISABLED'
      write(*, '(2x,a)') 'DISABLED TestPotMEnEnergDiff!'
    end if

    if (testVeff1IRMESternh) then
     write(*, '(2x,a)') 'Performing TestIRV1ME...'
      call TestIRV1ME( atoms, kpts, cell, lathar, input, stars, dimens, V0Fleur, results, memd_atom, logUnit, rho0MT, clnu_atom, nmem_atom, mlh_atom, &
        & GbasVec, ilst, nv, ne, nobd, z, rho0IR, kpq2kPrVec)
    else
      write(*, '(2x,a)') 'DISABLED TestIRV1ME!'
      write(logUnit, '(a)') 'Testing different ways to calculate <Psi|Veff|Psi>_IR'
      write(logUnit, '(a)') '-----------------------------------------------------'
      write(logUnit, '(a)') '                                                    |_DISABLED'
    end if

    if (testEps1q0) then
      write(*, '(2x,a)') 'Performing TestKSENerg1stVar...'
      call TestKSEnerg1stVar( atoms, cell, lathar, enpara, stars, dimens, sym, input, ud, V0Fleur, kpts, qpts, results,       &
      & memd_atom, ngdp, logUnit, gdp, rho0IR, rho0MT, mlh_atom, nmem_atom, clnu_atom, rbas1, rbas2, uuilon, duilon, ulouilopn, ilo2p,&
      & ne, El, nv, nobd, GbasVec, ilst, kveclo, iloTable, z, kpq2kPrVec, eig, nRadFun, vEff0MTsh, vEff0IRpwUw )
    else
      write(*, '(2x,a)') 'DISABLED TestKSEnerg1stVar!'
      write(logUnit, '(a)') 'Test vanishing of epsilon1 for q = 0 for monoatomic system.'
      write(logUnit, '(a)') '-----------------------------------------------------------'
      write(logUnit, '(a)') '                                                          |_DISABLED'
    end if

    ! Compare surface integral for Sternheimer in IR and MT
    if (testCompareSurfInt) then
     write(*, '(2x,a)') 'Performing CompareSurfInt...'
      call CompareSurfInt( atoms, lathar, kpts, qpts, dimens, sym, V0Fleur, cell, input, ud, results, stars, ngdp, gdp, El, rbas1, &
        & rbas2, nRadFun, mlh_atom, nmem_atom, clnu_atom, ilst, nv, GbasVec, ne, z, kveclo, eig, iloTable, nobd, kpq2kPrVec, vEff0IRpwUw )
    else
      write(*, '(2x,a)') 'DISABLED CompareSurfInt!'
    end if

    if (testSplitMTSurfIntSterh) then
      ! Check splitted routine to calculate MT
     write(*, '(2x,a)') 'Performing CheckSplittedMTSurfIntRoutine...'
      call CheckSplittedMTSurfIntRoutine( atoms, dimens, kpts, lathar, V0Fleur, sym, cell, input, ud, kveclo, nobd, ne, iloTable, &
        & eig, nv, GbasVec, ilst, nRadFun, El, rbas1(:, :, :, :, 1), rbas2(:, :, :, :, 1), nmem_atom, mlh_atom, clnu_atom, z )
    else
      write(*, '(2x,a)') 'DISABLED CheckSplittedMTSurfIntRoutine!'
    end if

    !todo test the consistency of the kinetic energy treatment in the MT (it should be the same than in the IR region.

  end subroutine TestSternheimerSetup

  !--------------------------------------------------------------------------------------------------------------------------------------
  !> @author
  !> Christian-Roman Gerhorst, Gregor Michaliczek
  !>
  !> @brief
  !> Tests for Sternheimer matrix element containing potential's gradient and first variation of potential, as well as for the left side
  !> of the Sternheimer equation.
  !>
  !> @details
  !> A constant potential which is set to one within the interstitial region and the muffin-tins is is used as parameter for calculating
  !> the matrix element for the potentials within the Sternheimer equation.
  !> The matrices for IR and MT should be a unit matrix if q = 0, if q /= 0 they should be zero. In any case the resulting z1 should be
  !> zero.
  !--------------------------------------------------------------------------------------------------------------------------------------
  subroutine TestPotMEnEnergDiff( atoms, cell, sym, stars, lathar, enpara, kpts, qpts, ud, input, dimens, ngdp, logUnit, ne, nobd, eig, &
      & rbas1, rbas2, uuilon, duilon, ulouilopn, GbasVec, ilst, nv, z, kveclo, gdp, iloTable, nRadFun, El, mapKpq2K, V0Fleur, nmem_atom, clnu_atom, mlh_atom, ilo2p, kpq2kPrVec )

    use m_types, only : t_atoms, t_cell, t_sym, t_stars, t_sphhar, t_enpara, t_kpts, t_usdus, t_input, t_dimension, t_tlmplm, t_noco, t_potential
    use m_abcof
     
    use m_jpConstants, only : fpi
    use m_juDFT_NOstopNO, only : juDFT_warn
    use m_jpSternhPulaySurface, only : tlmplm4H0
    use m_jpSternhHF, only : tlmplm4V
    use m_jpSternheimer, only : GenRecEdifME
    use m_jpPlotObservables, only : plotPotPathUC
    use m_jpPotDensHelper, only : warpIRPot
    implicit none

    ! Type parameter
    type(t_atoms),                  intent(in) :: atoms
    type(t_cell),                   intent(in) :: cell
    type(t_sym),                    intent(in) :: sym
    type(t_stars),                  intent(in) :: stars
    type(t_sphhar),                 intent(in) :: lathar
    type(t_enpara),                 intent(in) :: enpara
    type(t_kpts),                   intent(in) :: kpts
    type(t_kpts),                   intent(in) :: qpts
    type(t_usdus),                  intent(in) :: ud
    type(t_input),                  intent(in) :: input
    type(t_dimension),              intent(in) :: dimens
    type(t_potential),              intent(in) :: V0Fleur

    ! Scalar parameter
    integer,                        intent(in) :: ngdp
    integer,                        intent(in) :: logUnit

    ! Array parameter
    integer,                        intent(in) :: ne(:)
    integer,                        intent(in) :: nobd(:, :)
    real,                           intent(in) :: eig(:, :, :)
    real,                           intent(in) :: rbas1(:, :, 0:, :, :)
    real,                           intent(in) :: rbas2(:, :, 0:, :, :)
    real,                           intent(in) :: uuilon(:, :)
    real,                           intent(in) :: duilon(:, :)
    real,                           intent(in) :: ulouilopn(:, :, :)
    integer,                        intent(in) :: GbasVec(:, :)
    integer,                        intent(in) :: ilst(:, :, :)
    integer,                        intent(in) :: nv(:, :)
    MCOMPLEX,                       intent(in) :: z(:,:,:,:)
    integer,                        intent(in) :: kveclo(:,:)
    integer,                        intent(in) :: gdp(:, :)
    integer,                        intent(in) :: iloTable(:, 0:, :)
    integer,                        intent(in) :: nRadFun(0:, :)
    real,                           intent(in) :: El(:, 0:, :, :)
    integer,                        intent(in) :: mapKpq2K(:, :) ! todo possible error here
    integer, intent(in) :: ilo2p(:, :)
    integer, intent(in) :: kpq2kPrVec(:, :, :)

    ! Local type variable
    type(t_tlmplm)                             :: tdVx
    type(t_tlmplm)                             :: tdVy
    type(t_tlmplm)                             :: tdVz
    type(t_tlmplm)                             :: tdHS0
    type(t_noco)                               :: noco
    type(od_inp)                               :: odi
    type(od_sym)                               :: ods

    ! Local scalar variables
    integer                                    :: iqpt
    integer                                    :: ikpt
    integer                                    :: dispAtInd
    integer                                    :: dispType
    integer                                    :: idir
    integer                                    :: ieqat
    integer                                    :: nmat
    integer                                    :: itype
    integer                                    :: iatom
    integer                                    :: lmp
    integer                                    :: lm
    integer                                    :: oqn_l
    integer                                    :: mqn_m
    integer                                    :: p
    integer                                    :: ilo
    integer                                    :: nband
    integer                                    :: pband
    real                                       :: testPotMTval
    integer                                    :: maxlmp
    integer                                    :: ieqatp
    complex, allocatable                                   :: testPotME(:, :)

    ! Local array variables
    real,              allocatable             :: recEdiffME(:, :)
    complex,           allocatable             :: testPotIR(:, :)
    complex,           allocatable             :: w_testPotIR(:, :)
    complex,           allocatable             :: testIRPotME(:, :)
    complex,           allocatable             :: testPotMT(:, :, :, :)
    complex,           allocatable             :: testMTPotME(:, :)



    integer :: ii
    complex                        :: vr0SpH(atoms%jmtd, (atoms%lmaxd + 1)**2, atoms%nat, 3)
    integer                        :: ptsym
    integer :: ilh,  imem, imesh
    integer,            intent(in) :: nmem_atom(0:, :)
    integer,            intent(in) :: mlh_atom(:, 0:, :)
    complex,            intent(in) :: clnu_atom(:, 0:, :)
    complex,           allocatable              :: grRho0IR(:, :)
    integer :: iG
    integer, allocatable :: nlo_atom(:)
    real,              allocatable   :: cutContr(:, :)

    ! Set boundary for lm with LOs
    maxlmp = maxval( (/ (sum( (/ ((2*oqn_l+1)* nRadFun(oqn_l,dispType), oqn_l = 0,atoms%lmax(dispType)) /) ),dispType=1,atoms%ntype) /) )

    ! Construct constant potential and plot it
    ! Create constant potential (set value to 1) within whole unit cell, i.e. interstitial region and muffin-tin.
    allocate ( testPotIR(ngdp, 3), testPotMT(atoms%jmtd, (atoms%lmaxd + 1)**2, 3, atoms%nat) ) !todo l should start at 0!!!
    testPotIR = 0
    testPotMT = 0
    allocate( recEdiffME(dimens%neigd, maxval(nobd)) )
    allocate( cutContr(dimens%neigd, maxval(nobd)) )

    do ii = 1, ngdp
      if (all( gdp(:, ii) == 0 )) then
        testPotIR(ii, :) = 1
        exit
      end if
    end do

    testPotMTval = sqrt(fpi)
    testPotMT(:, 1, :, :) = testPotMTval

    !todo if activated we overwrite the results of the gradient of the potential test!!!!!
 !   call plotPotPathUC(cell, atoms, testPotIR, ngdp, gdp, testPotMT, 1., 1., 1.)

    ! Warp interstitial potential
!    allocate(w_testPotIR(ngdp, 3))
!    do idir = 1, 3
!    call WarpIRPot(stars, ngdp, idir, gdp, testPotIR, w_testPotIR(:, idir))
!    end do

    !todo change order of indices in vr so nat and idir
   iqpt = 1 !todo beware that the right loop over the types is performed
         ! This is written 3 times because of the type array which is not necessarily be stored linearily, so due to performance reasons
         ! Unfold symmetry-compressed muffin-tin density because juPhon routines do not implement symmetry
!         iatom = 0
!         !todo look at the indices!!!!!!
!          vr0SpH = 0
!        do itype = 1, atoms%ntype
!          do ieqat = 1, atoms%neq(itype)
!          iatom = iatom + 1
!          ptsym = atoms%ntypsy(iatom)
!          do ilh = 0, lathar%nlh(ptsym)
!            oqn_l = lathar%llh(ilh, ptsym)
!            do imem = 1, nmem_atom(ilh, iatom)
!              mqn_m = mlh_atom(imem, ilh, iatom)
!              lm = oqn_l * (oqn_l + 1) + 1 + mqn_m
!              do imesh = 1, atoms%jri(itype)
!                vr0SpH(imesh, lm, 1, iatom) = vr0SpH(imesh, lm, 1, iatom) + V0Fleur%vr(imesh, ilh, itype, 1) * clnu_atom(imem, ilh, iatom)
!              end do ! imesh
!            end do ! imem
!          end do ! ilh
!        end do
!      end do
          call tlmplm4V( atoms, lathar, dimens, enpara, ud, input, tdVx, input%jspins, input%jspins, 1, testPotMT, rbas1, rbas2, &
          &uuilon, duilon, ulouilopn, ilo2p, nlo_atom )

        !if (all(atoms%nlo(:) == 0)) allocate(nlo_atom(atoms%nat))

!    write ( 7550, '(2(f24.15))' ) tdVx%tuu!(:, 1, 1) !all lmplm, frist atom, first spin
!    write ( 7551, '(2(f24.15))' ) tdVx%tud!(:, 1, 1) ! all lmplm, first atom, first spin
!    write ( 7552, '(2(f24.15))' ) tdVx%tdu!(:, 1, 1) ! all lmplm, first atom, first spin
!    write ( 7553, '(2(f24.15))' ) tdVx%tdd!(:, 1, 1) ! all lmplm, first atom, first spin
!    write ( 7554, '(i8)' ) tdVx%ind!(:, :, 1, 1)
!    write ( 7555,  '(2(f24.15))') tdVx%tuulo
!    write ( 7556,  '(2(f24.15))') tdVx%tdulo
!    write ( 7557,  '(2(f24.15))') tdVx%tuloulo

   !       call tlmplm4V( atoms, lathar, dimens, enpara, ud, input, tdVy, input%jspins, input%jspins, 2, testPotMT, rbas1, rbas2, &
   !       &uuilon, duilon, ulouilopn )
   !       call tlmplm4V( atoms, lathar, dimens, enpara, ud, input, tdVz, input%jspins, input%jspins, 3, testPotMT, rbas1, rbas2, &
   !       &uuilon, duilon, ulouilopn )
          do ikpt = 1, kpts%nkpt
            recEdiffME(:, :)             = 0.
            cutContr(:, :)               = 0.
            call GenRecEdifME( ne(mapKpq2K(ikpt, iqpt)), nobd(ikpt, 1), eig(:, mapKpq2K(ikpt, iqpt), 1), eig(:, ikpt, 1), recEdiffME, cutContr )
            call TestPotMEnEnergDiffHelp( atoms, stars, input, dimens, sym, cell, kpts, logUnit, ud, tdVx, ne, El, nobd, kveclo, &
              & iloTable, nRadFun, recEdiffME, ikpt, mapKpq2K(ikpt, iqpt), iqpt, maxlmp, 1, ngdp, nv, gdp, GbasVec, ilst, z, &
        !      & grRho0IR, testPotME)
              & testPotIR, testPotME, nlo_atom, kpq2kPrVec)
   !         call TestPotMEnEnergDiffHelp( atoms, stars, input, dimens, sym, cell, kpts, logUnit, ud, tdVy, ne, El, nobd, kveclo, &
   !           & iloTable, nRadFun, recEdiffME,ikpt, mapKpq2K(ikpt, iqpt), iqpt, maxlmp, 2, ngdp, nv, gdp, GbasVec, ilst, z, &
   !           & w_testPotIR, testPotME)
   !         call TestPotMEnEnergDiffHelp( atoms, stars, input, dimens, sym, cell, kpts, logUnit, ud, tdVz, ne, El, nobd, kveclo, &
   !           & iloTable, nRadFun, recEdiffME, ikpt,mapKpq2K(ikpt, iqpt), iqpt, maxlmp, 3, ngdp, nv, gdp, GbasVec, ilst, z, &
   !           & w_testPotIR,  testPotME)

          end do ! ikpt
        !end do ! ieqt
      !end do ! dispType
   ! end do ! iqpt

   ! deallocate(testMTPotME)

   ! deallocate(testIRPotME)
    write(logUnit, '(a)') 'Test of potential matrix element left side in Sternheimer equation'
    write(logUnit, '(a)') '------------------------------------------------------------------'
    write (logUnit, '(a)') '                                                                |__passed'

  end subroutine TestPotMEnEnergDiff
! Helper of TEstPotMenEnergyDiff
  subroutine TestPotMEnEnergDiffHelp(atoms, stars, input, dimens, sym, cell, kpts, logUnit, ud, tdV, ne, El, nobd, kveclo, iloTable, &
      & nRadFun, recEdiffME, ikpt, ikpq, iqpt, maxlmp, idir, ngdp, nv,  gdp, GbasVec, ilst, z, w_testPotIR, testPotME, nlo_atom, kpq2kPrVec)

    use m_types, only : t_atoms, t_cell, t_sym, t_tlmplm, t_stars, t_kpts, t_usdus, t_input, t_dimension, t_noco
    use m_abcof
     
    use m_jpSternhHF, only : calcMEPotIR, calcVsumMT

    use m_juDFT_NOstopNO, only : juDFT_warn

    implicit none

    type(t_atoms),        intent(in) :: atoms
    type(t_stars),        intent(in) :: stars
    type(t_input),        intent(in) :: input
    type(t_dimension),    intent(in) :: dimens
    type(t_sym),          intent(in) :: sym
    type(t_cell),         intent(in) :: cell
    type(t_kpts),         intent(in) :: kpts
    type(t_usdus),        intent(in) :: ud
    type(t_tlmplm),       intent(in) :: tdV

    integer,              intent(in) :: ikpt
    integer,              intent(in) :: ikpq
    integer,              intent(in) :: iqpt
    integer,              intent(in) :: ngdp
    integer,              intent(in) :: idir
    integer,              intent(in) :: maxlmp

    integer,              intent(in) :: nv(:, :)
    integer, intent(in) :: kpq2kPrVec(:, :, :)
    real,                           intent(in) :: El(:, 0:, :, :)
    complex,           intent(in)             :: w_testPotIR(:, :)
    integer,                        intent(in) :: gdp(:, :)
    integer,                        intent(in) :: GbasVec(:, :) ! can be reduced in ket and bra
    integer,                        intent(in) :: ilst(:, :, :)
    MCOMPLEX,                       intent(in) :: z(:,:,:,:)
    integer,                        intent(in) :: ne(:)
    integer,                        intent(in) :: nobd(:, :)
    integer,                        intent(in) :: kveclo(:,:)
    integer,                        intent(in) :: iloTable(:, 0:, :)
    integer,                        intent(in) :: logUnit
    integer,                        intent(in) :: nRadFun(0:, :)
    real,                        intent(in) :: recEdiffME(:, :)
    complex,      allocatable,               intent(out) :: testPotME(:, :)

    complex,      allocatable                            :: testMTPotME(:, :)
    integer, intent(in) :: nlo_atom(:)

    type(t_noco)                               :: noco
    type(od_inp)                               :: odi
    type(od_sym)                               :: ods

    integer                          :: nmat
    integer :: iatom
    integer :: itype
    integer :: ieqat
    integer :: lmp
    integer :: lm
    integer :: oqn_l
    integer :: mqn_m
    integer :: p
    integer :: ilo

    complex, allocatable             :: testIRPotME(:, :) ! reduce to 1 dir!
    complex, allocatable             :: testIRPotME2(:, :) ! reduce to 1 dir!
    complex,           allocatable             :: acofKet(:, :, :)
    complex,           allocatable             :: bcofKet(:, :, :)
    complex,           allocatable             :: ccofKet(:, :, :, :)
    complex,           allocatable             :: acofBra(:, :, :)
    complex,           allocatable             :: bcofBra(:, :, :)
    complex,           allocatable             :: ccofBra(:, :, :, :)
    complex,           allocatable             :: mCoefK(:, :, :)
    complex,           allocatable             :: mCoefB(:, :, :)
    complex,           allocatable             :: z1Band(:, :)
    integer :: nband
    integer :: pband
    integer :: braBands
    integer :: ketBands
    integer :: iBasf
    complex, allocatable :: z1G(:, :)
    integer                          :: nobdHelp(64, 1)

    real :: foo, bar
    integer:: k1(dimens%nvd), k2(dimens%nvd), k3(dimens%nvd)
    integer:: k1Help(dimens%nvd,1), k2Help(dimens%nvd,1), k3Help(dimens%nvd,1)
   ! nobdHelp(:, 1) = ne(:)

   integer :: ii
   integer :: gdpP(ngdp)
   integer :: gdpR(0:stars%k1d, 0:stars%k2d + 1, 0:stars%k3d + 1)
   logical :: passed = .true.
   integer, allocatable :: ngoprI(:)

    braBands = ne(ikpq)
    ketBands = nobd(ikpt, 1)
    ! Generate inner part of braket
    ! Generate IR part of potential matrix element
    allocate( testIRPotME(braBands, ketBands) )
    ! nmat has to be from bra
    nmat = nv(1, ikpq) + atoms%nlotot

    allocate(ngoprI(atoms%nat))

    ngoprI(:) = 1

    !the z which is used in calcMEPotIR is the one taken from the eig file in fleur in aline before conjugating it!
    call calcMEPotIR( stars, dimens, GbasVec(:, ilst(:nv(input%jspins, ikpq), ikpq, input%jspins)), &
      & GbasVec(:, ilst(:nv(input%jspins, ikpt), ikpt, input%jspins)), nv, w_testPotIR(:,idir), z(:, :, ikpq, 1), &
      & z(:, :, ikpt, 1), gdp, nmat, ne(ikpq), nobd(ikpt, 1), ikpt, 1, ikpq, ngdp, testIRPotME, kpq2kPrVec ) !todo spin-relevant



    ! todo This is the complex conjugate of the potential we want to
   ! do nband = 1, ketBands
   !   do pband = 1, braBands
   !     write (167, '(i4,i4,2(f24.12))') pband, nband, real(testIRPotME(pband, nband)), -aimag(testIRPotME(pband, nband))
   !   end do
   ! end do



    ! Generate matching coefficients for the bra and ket
    ! Needed otherwise it does not run
    allocate(acofKet(dimens%neigd, 0:dimens%lmd, atoms%nat), bcofKet(dimens%neigd, 0:dimens%lmd, atoms%nat), &
      &ccofKet(-atoms%llod:atoms%llod, dimens%neigd, atoms%nlod, atoms%nat))
    nmat = nv(1, ikpt) + atoms%nlotot
    allocate(noco%alph(atoms%ntype), noco%beta(atoms%ntype)) !Up to now those variables are only of dummy character
    noco%alph = 0
    noco%beta = 0
    call abcof ( atoms%lmaxd, atoms%ntype, dimens%neigd, dimens%neigd, atoms%nat, sym%nop, dimens%nvd, dimens%jspd, dimens%lmd, &
      & dimens%nbasfcn, atoms%llod, atoms%nlod, atoms%nlotot, sym%invtab, atoms%ntype, sym%mrot, ngoprI, atoms%taual, atoms%neq, &
      & atoms%lmax, atoms%rmt, cell%omtil, cell%bmat, cell%bbmat, kpts%bk(:, ikpt), &
      & GbasVec(1, ilst(:nv(input%jspins, ikpt), ikpt, input%jspins)), GbasVec(2, ilst(:nv(input%jspins, ikpt), ikpt, input%jspins)), &
      & GbasVec(3, ilst(:nv(input%jspins, ikpt), ikpt, input%jspins)), nv(:, ikpt),  nmat, ne(ikpt), z(:, :, ikpt, 1), ud%us(:, :, 1), &
      & ud%dus(:, :, 1), ud%uds, ud%duds(:, :, 1), ud%ddn(:, :, 1), atoms%invsat, sym%invsatnr, ud%ulos(:, :, 1), &
      & ud%uulon(:, :, 1), ud%dulon(:, :, 1),  ud%dulos(:, :, 1), atoms%llo, atoms%nlo, atoms%l_dulo, &
      & atoms%lapw_l, noco%l_noco, noco%l_ss, 1, noco%alph, noco%beta, noco%qss, kveclo(:, ikpt), odi, ods, &
      & acofKet, bcofKet, ccofKet)

    allocate(acofBra(dimens%neigd, 0:dimens%lmd, atoms%nat), bcofBra(dimens%neigd, 0:dimens%lmd, atoms%nat),&
      &ccofBra(-atoms%llod:atoms%llod, dimens%neigd, atoms%nlod, atoms%nat))
    nmat = nv(1, ikpq) + atoms%nlotot
    call abcof ( atoms%lmaxd, atoms%ntype, dimens%neigd, dimens%neigd, atoms%nat, sym%nop, dimens%nvd, dimens%jspd, dimens%lmd, &
      & dimens%nbasfcn, atoms%llod, atoms%nlod, atoms%nlotot, sym%invtab, atoms%ntype, sym%mrot, ngoprI, atoms%taual, atoms%neq, &
      & atoms%lmax, atoms%rmt, cell%omtil, cell%bmat, cell%bbmat, kpts%bk(:, ikpq), &
      & GbasVec(1, ilst(:nv(input%jspins, ikpq), ikpq, input%jspins)), GbasVec(2, ilst(:nv(input%jspins, ikpq), ikpq, input%jspins)), &
      & GbasVec(3, ilst(:nv(input%jspins, ikpq), ikpq, input%jspins)), nv(:, ikpq),  nmat, braBands, z(:, :, ikpq, 1), ud%us(:, :, 1), &
      & ud%dus(:, :, 1), ud%uds, ud%duds(:, :, 1), ud%ddn(:, :, 1), atoms%invsat, sym%invsatnr, ud%ulos(:, :, 1), &
      & ud%uulon(:, :, 1), ud%dulon(:, :, 1), ud%dulos(:, :, 1), atoms%llo, atoms%nlo, atoms%l_dulo, &
      & atoms%lapw_l, noco%l_noco, noco%l_ss, 1, noco%alph, noco%beta, noco%qss, kveclo(:, ikpq), odi, ods, &
      & acofBra, bcofBra, ccofBra)

    allocate(mCoefK(ketBands, maxlmp, atoms%nat))
    allocate(mCoefB(braBands, maxlmp, atoms%nat))
    iatom  = 0
    mCoefK =  0
    mCoefB = 0
    do itype = 1, atoms%ntype
      do ieqat = 1, atoms%neq(itype)
        iatom = iatom + 1
        lmp   =   0
        lm    =  -1
        do oqn_l = 0, atoms%lmax(itype)
          do mqn_m = -oqn_l, oqn_l
            lm = lm + 1
            !p = 1
            lmp = lmp + 1
            mCoefK(:ketBands, lmp, iatom) = acofKet(:ketBands, lm, iatom)
            mCoefB(:braBands, lmp, iatom) = acofBra(:braBands, lm, iatom)
            !p = 2
            lmp = lmp + 1
            mCoefK(:ketBands, lmp, iatom) = bcofKet(:ketBands, lm, iatom)
            mCoefB(:braBands, lmp, iatom) = bcofBra(:braBands, lm, iatom)
            !LOs
            do p = 3, nRadFun(oqn_l, itype)
              ilo = iloTable(p, oqn_l, itype)
              lmp = lmp + 1
              mCoefK(:ketBands, lmp, iatom) = ccofKet(mqn_m, :ketBands, ilo, iatom)
              mCoefB(:braBands, lmp, iatom) = ccofBra(mqn_m, :braBands, ilo, iatom)
            end do
          end do
        end do
      end do
    end do
    deallocate(acofKet, bcofKet, ccofKet, acofBra, bcofBra, ccofBra)
    ! Rearrange matching coefficients for a nice treatment of LOs
    ! Create MT part of potentials matrix element

    allocate( testMTPotME(braBands, ketBands ) )
    !todo why ist tdV equal here????
    call calcVsumMT( atoms, tdV, tdV, ikpt, ikpq, ne, nobd, mCoefB, mCoefK, nRadFun, iloTable, nlo_atom, testMTPotME )
    deallocate( mCoefK, mCoefB )


    allocate( testPotME(braBands, ketBands) )
    testPotME = testIRPotME + testMTPotME !todo this conjg is wrong here!!!!
!    deallocate( testIRPotME, testMTPotME )

    ! Create summand of perturbation theory for first-order state
    allocate( z1Band(braBands, ketBands) )
    z1Band = 0
    z1Band = - recEdiffME * testPotME ! todo here might be an error this should be multiplied elementwise

    do nBand = 1, ketBands
      do pBand = 1, braBands
        !write (586, '(i4, i4, 8(f15.8))') pBand, nBand, z1Band(pBand, nBand), testIRPotME(pBand, nBand), testMTPotME(pBand, nBand), testPotME(pBand, nBand)
        !todo this conjg is wrong here we have to trace back it!
        !write (586, '(i4, i4, 2(f24.12))') pBand, nBand, testPotME(pBand, nBand)
      end do
    end do
    do nBand = 1, ketBands
      do pBand = 1, braBands
        if ( pBand == nBand ) then
          !10^-9 is still very accurate due to the usage of the integration routines providing only accuracy of 1e-9
          if ( real( abs(testPotME(pBand, nBand) ) - 1.0 )  > 9e-8 .or. aimag( testPotME(pBand, nBand) ) > 9e-8 )  then
            write (logUnit, '(a)') 'Test of potential matrix element and left side: Diagonal element of right side is not equals 1!'
            write (logUnit, '("ikpt = ",i3,", pBand = ",i5,",&
               & nBand = ",i5,", z1 = ", 2(es15.8))') ikpt, pBand, nBand, testPotME(pBand, nBand)
            write(logUnit , *)
            passed = .false.
          end if
          if ( abs( z1Band(pBand, nBand) ) > 9e-8 ) then
            write (logUnit, '(a)') 'Test of potential matrix element and left side: z1 is not (0.0, 0.0)'
            write (logUnit, '("ikpt = ",i3,", pBand = ",i5,",&
               & nBand = ",i5,", z1 = ", 2(es15.8))') ikpt, pBand, nBand, testPotME(pBand, nBand)
            write(logUnit, *)
            passed = .false.
          end if
        else if ( pBand /= nBand ) then
          if ( abs(testPotME(pBand, nBand))  >= 9e-8 ) then
            write (logUnit, '(a)') 'Test of potential matrix element and left side: Non-diagonal element of right side is not equals 1!'
            write (logUnit, '("ikpt = ",i3,", pBand = ",i5,",&
               & nBand = ",i5,", z1 = ", 2(es15.8))') ikpt, pBand, nBand, testPotME(pBand, nBand)
            write(logUnit, *)
            passed = .false.
          end if
          if ( abs( z1Band(pBand, nBand) ) >= 9e-8 ) then
            write (logUnit, '(a)') 'Test of potential matrix element and left side: z1 is not (0.0, 0.0)'
            write (logUnit, '("ikpt = ",i3,", pBand = ",i5,",&
               & nBand = ",i5,", z1 = ", 2(es15.8))') ikpt, pBand, nBand, z1Band(pBand, nBand)
            write(logUnit, *)
            passed = .false.
          end if
        end if
      end do ! pBand
    end do ! nBand
    if ( .not.passed ) then
      write(logUnit, '(a)') 'Test of potential matrix element left side in Sternheimer equation'
      write(logUnit, '(a)') '------------------------------------------------------------------'
      write (logUnit, '(a)') '                                                                |__failed'
      call juDFT_warn('Test of potential matrix element and left side failed', hint='See juPhon.log for details.', &
        &calledby='TestPotMEnEnergDiffHelp')
    end if
  end subroutine TestPotMEnEnergDiffHelp
  !--------------------------------------------------------------------------------------------------------------------------------------
  !> @author
  !> Christian-Roman Gerhorst, Gregor Michaliczek
  !>
  !> @brief
  !> Test for Sternheimer matrix element containing unperturbed Hamiltonian and unperturbed overlap matrix.
  !>
  !> @details
  !> The Hamiltonian (using read-in Fleur potential) and the overlap matrix is calculated in the regular way using regular and same abcof
  !> coefficients for bra and ket instead of the special abcofs used in the Sternheimer equation. This is compared to a read-in
  !> Hamiltonian from a Fleur calculation as well as a read in overlap matrix.
  !>
  !> @todo review documentation when method is known
  !--------------------------------------------------------------------------------------------------------------------------------------
  subroutine TestHSME( atoms, sym, cell, lathar, dimens, enpara, ud, input, kpts, V0Fleur, logUnit, rbas1, rbas2, z, uuilon, duilon,&
    & ulouilopn, GbasVec, ilst, nv, ne, nobd, El, iloTable, nRadFun, kveclo, ilo2p, clnu_atom, nmem_atom, mlh_atom, vEff0MTsh )

#include "cppmacro.h"
    use m_types, only : t_atoms, t_sym, t_cell, t_sphhar, t_dimension, t_enpara, t_usdus, t_input, t_kpts, t_potential, t_tlmplm, t_noco
     
    use m_abcof
    use mod_juPhonUtils, only : fopen, fclose
    use m_juDFT_NOstopNO, only : juDFT_warn
    use m_jpSternhPulaySurface, only : calcHS0MT, tlmplm4H0

    implicit none

    ! Type parameters
    type(t_atoms),                  intent(in) :: atoms
    type(t_sym),                    intent(in) :: sym
    type(t_cell),                   intent(in) :: cell
    type(t_sphhar),                 intent(in) :: lathar
    type(t_dimension),              intent(in) :: dimens
    type(t_enpara),                 intent(in) :: enpara
    type(t_usdus),                  intent(in) :: ud
    type(t_input),                  intent(in) :: input
    type(t_kpts),                   intent(in) :: kpts
    type(t_potential),              intent(in) :: V0Fleur

    ! Scalar parameters
    integer,                        intent(in) :: logUnit

    ! Array parameters
    real,                           intent(in) :: rbas1(:, :, 0:, :, :)
    real,                           intent(in) :: rbas2(:, :, 0:, :, :)
    MCOMPLEX,                       intent(in) :: z(:,:,:,:)
    real,                           intent(in) :: uuilon(:, :)
    real,                           intent(in) :: duilon(:, :)
    real,                           intent(in) :: ulouilopn(:, :, :)
    integer,                        intent(in) :: GbasVec(:, :)
    integer,                        intent(in) :: ilst(:, :, :)
    integer,                        intent(in) :: nv(:, :)
    integer,                        intent(in) :: ne(:)
    integer,                        intent(in) :: nobd(:, :)
    real,                           intent(in) :: El(:, 0:, :, :)
    integer,                        intent(in) :: iloTable(:, 0:, :)
    integer,                        intent(in) :: nRadFun(0:, :)
    integer,                        intent(in) :: kveclo(:,:)
    integer,                        intent(in) :: ilo2p(:, :)
    complex,                           intent(in)  :: clnu_atom(:, 0:, :)
    integer,                           intent(in)  :: nmem_atom(0:, :)
    integer,                           intent(in)  :: mlh_atom(:, 0:, :)
    complex,                         intent(in) :: vEff0MTsh(:, :, :, :)

    ! Local type variables
    type(t_tlmplm)                             :: tdHS0
!    type(t_tlmplm)                             :: tdHS0Dum
    type(t_noco)                               :: noco
    type(od_inp)                               :: odi
    type(od_sym)                               :: ods

    ! Local scalar variables
    integer                                    :: maxlmp
    integer                                    :: oqn_l
    integer                                    :: dispType
    integer                                    :: ikpt
    integer                                    :: itype
    integer                                    :: ieqat
    integer                                    :: iatom
    integer                                    :: lmp
    integer                                    :: p
    integer                                    :: ilo
    integer                                    :: nmat
    integer                                    :: lm
    integer                                    :: mqn_m
    integer                                    :: nband
    integer                                    :: pband
    integer                                    :: matsize
    integer                                    :: ieig

    ! Local array variables
    complex,           allocatable             :: acof(:, :, :)
    complex,           allocatable             :: bcof(:, :, :)
    complex,           allocatable             :: ccof(:, :, :, :)
!    complex,           allocatable             :: acof2(:, :, :)
!    complex,           allocatable             :: bcof2(:, :, :)
!    complex,           allocatable             :: ccof2(:, :, :, :)
    complex,           allocatable             :: acofFleur(:, :, :)
    complex,           allocatable             :: bcofFleur(:, :, :)
    complex,           allocatable             :: ccofFleur(:, :, :, :)
    complex,           allocatable             :: mCoef(:, :, :)
    MCOMPLEX,          allocatable             :: hJuPhon(:, :)
    MCOMPLEX,          allocatable             :: sJuPhon(:, :)
    MCOMPLEX,          allocatable             :: hJuPhonTemp(:, :)
    MCOMPLEX,          allocatable             :: sJuPhonTemp(:, :)
    MCOMPLEX,          allocatable             :: hFleur(:, :)
    MCOMPLEX,          allocatable             :: sFleur(:, :)
    MCOMPLEX,          allocatable             :: a(:)
    MCOMPLEX,          allocatable             :: b(:)
    MCOMPLEX,          allocatable             :: help(:)

    complex                                    :: CPP_BLAS_cdotc
    external                                   :: CPP_BLAS_cdotc, CPP_BLAS_chpmv

    complex                                    :: myzs(dimens%nbasfcn, dimens%neigd)
    integer :: ii, jj
    integer, allocatable :: nobdNe(:, :)
    logical :: passedS, passedH
!    integer ,allocatable                :: k1(:)                   ! temporary x coordinates of basis vectors for a k-point
!    integer ,allocatable                :: k2(:)                   ! temporary y coordinates of basis vectors for a k-point
!    integer ,allocatable                :: k3(:)                   ! temporary z coordinates of basis vectors for a k-point
!    integer:: k1Help(dimens%nvd,1), k2Help(dimens%nvd,1), k3Help(dimens%nvd,1)
!    integer :: ll

    integer, allocatable :: ngoprI(:)

    integer, allocatable              :: nlo_atom(:)
    integer :: mytype
    integer :: myatom
    integer :: myeqat
    !shift for LOs
    allocate(nlo_atom(atoms%nat))
    myatom = 0
    do mytype = 1, atoms%ntype
      do myeqat = 1, atoms%neq(mytype)
        myatom = myatom + 1
        nlo_atom(myatom) = atoms%nlo(mytype)
      end do
    end do
    ! Create inner part for Hamilton and overlap braket
    call tlmplm4H0( atoms, dimens, enpara, ud, input, tdHS0, 1, logUnit, rbas1, rbas2, uuilon, duilon, ulouilopn, ilo2p, &
                                                                                                           & vEff0MTsh(:, :, :, 1) )



    ! Generate bra and ket
    !todo what is ikpt, don't choose a ikpt which is the Gamma point
    ikpt = 28
    itype = 1
    iatom = 1

    allocate( nobdNe(size(nobd(:, 1)), 1) )
    nobdNe = 0
    nobdNe(:, 1) = ne(:)

    allocate(ngoprI(atoms%nat))
    ngoprI(:) = 1

    allocate(noco%alph(atoms%ntype), noco%beta(atoms%ntype)) !Up to now those variables are only of dummy character
    noco%alph = 0
    noco%beta = 0
    allocate(acof(dimens%neigd, 0:dimens%lmd, atoms%nat), bcof(dimens%neigd, 0:dimens%lmd, atoms%nat), ccof(-atoms%llod:atoms%llod, dimens%neigd, atoms%nlod, atoms%nat))
    nmat = nv(1, ikpt) + atoms%nlotot
    call abcof ( atoms%lmaxd, atoms%ntype, dimens%neigd, dimens%neigd, atoms%nat, sym%nop, dimens%nvd, dimens%jspd, dimens%lmd, &
      & dimens%nbasfcn, atoms%llod, atoms%nlod, atoms%nlotot, sym%invtab, atoms%ntype, sym%mrot, ngoprI, atoms%taual, atoms%neq, &
      & atoms%lmax, atoms%rmt, cell%omtil, cell%bmat, cell%bbmat, kpts%bk(:, ikpt), &
      & GbasVec(1, ilst(:nv(input%jspins, ikpt), ikpt, input%jspins)), GbasVec(2, ilst(:nv(input%jspins, ikpt), ikpt, input%jspins)), &
      & GbasVec(3, ilst(:nv(input%jspins, ikpt), ikpt, input%jspins)), nv(:, ikpt),  nmat, ne(ikpt), z(:, :, ikpt, 1), ud%us(:, :, 1), &
      & ud%dus(:, :, 1), ud%uds, ud%duds(:, :, 1), ud%ddn(:, :, 1), atoms%invsat, sym%invsatnr, ud%ulos(:, :, 1), &
      & ud%uulon(:, :, 1), ud%dulon(:, :, 1),  ud%dulos(:, :, 1), atoms%llo, atoms%nlo, atoms%l_dulo, &
      & atoms%lapw_l, noco%l_noco, noco%l_ss, 1, noco%alph, noco%beta, noco%qss, kveclo(:, ikpt), odi, ods, &
      & acof, bcof, ccof)

    ! Rearrange abcof to treat LOs more nicely
      maxlmp = maxval( (/ (sum( (/ ((2*oqn_l+1)* nRadFun(oqn_l,dispType), oqn_l = 0,atoms%lmax(dispType)) /) ),dispType=1,atoms%ntype) /) )
      allocate(mCoef(ne(ikpt), maxlmp, atoms%nat))
      iatom  = 0
      mCoef = 0
      do itype = 1, atoms%ntype
        do ieqat = 1, atoms%neq(itype)
          iatom = iatom + 1
          lmp   =   0
          lm    =  -1
          do oqn_l = 0, atoms%lmax(itype)
            do mqn_m = -oqn_l, oqn_l
              lm = lm + 1
              !p = 1
              lmp = lmp + 1
              mCoef(:ne(ikpt), lmp, iatom) = acof(:ne(ikpt), lm, iatom)
              !p = 2
              lmp = lmp + 1
              mCoef(:ne(ikpt), lmp, iatom) = bcof(:ne(ikpt), lm, iatom)
              !LOs
              do p = 3, nRadFun(oqn_l, itype)
                ilo = iloTable(p, oqn_l, itype)
                lmp = lmp + 1
                mCoef(:ne(ikpt), lmp, iatom) = ccof(mqn_m, :ne(ikpt), ilo, iatom)
              end do
            end do
          end do
        end do
      end do



    itype = 1
    iatom = 1

    ! Calculate Hamiltonian and overlap matrix
    allocate( hJuPhon(ne(ikpt), ne(ikpt)), sJuPhon(ne(ikpt), ne(ikpt)) )
    allocate( hJuPhonTemp(ne(ikpt), ne(ikpt)), sJuPhonTemp(ne(ikpt), ne(ikpt)) )
    hJuPhon = 0
    hJuPhonTemp = 0
    sJuPhon = 0
    sJuPhonTemp = 0
    iatom = 0
    do itype = 1, atoms%ntype
      do ieqat = 1, atoms%neq(itype)
        iatom = iatom + 1
        call calcHS0MT( atoms, ud, tdHS0, ikpt, ikpt, itype, iatom, ne, nobdNe, El, mCoef(:, :, iatom), mCoef(:, :, iatom), &
          & nRadFun, iloTable, nlo_atom, sJuPhonTemp, hJuPhonTemp )
        ! Symmetrize kinetic energy of MT matrix element again to be consistent with FLEUR
        lmp = 0
        do oqn_l = 0, atoms%lmax(itype)
          do mqn_m = -oqn_l, oqn_l
            ! go to acof of next lm
            lmp = lmp + 1
            do nBand = 1, nobdNe(ikpt, 1)
              do pBand = 1, ne(ikpt)
                hJuPhonTemp(pBand, nBand) = hJuPhonTemp(pBand, nBand) + 0.5 * &
                  & (conjg(mCoef(pband, lmp + 1, iatom)) * mCoef(nBand, lmp, iatom) - conjg(mCoef(pband, lmp, iatom)) * mCoef(nBand, lmp + 1, iatom))
!            write(2128, '(5i8, 2f15.8)') oqn_l, mqn_m, lmp, nBand, pBand, 0.5 * (conjg(mCoef(pband, lmp, iatom)) * mCoef(nBand, lmp + 1, iatom) + conjg(mCoef(pband, lmp + 1, iatom)) * mCoef(nBand, lmp, iatom ))
!            write(2129, '(5i8, 2f15.8)') oqn_l, mqn_m, lmp, nBand, pBand, conjg(mCoef(pband, lmp, iatom)) * mCoef(nBand, lmp + 1, iatom)! + conjg(mCoefBp(pband, lmlo + 1)) * mCoefKb(nBand, lmlo ))

              end do ! pband
            end do ! nband
            lmp = lmp + 1
          end do ! mqn_m
        end do ! oqn_l
        sJuPhon = sJuPhon + sJuPhonTemp
        hJuPhon = hJuPhon + hJuPhonTemp
      end do
    end do

!    do pBand = 1, ne(ikpt)
!      do nBand = 1, pband
!        if ( abs(real(sJuPhon(pBand, nBand))) < 1e-9 ) then
!          write (950, '(i8, i8, 2(f24.8))') pBand, nBand, 0.0
!        else
!          write (950, '(i8, i8, 2(f24.8))') pBand, nBand, real(sJuPhon(pBand, nBand))
!        end if
!        if ( abs(aimag(sJuPhon(pBand, nBand))) < 1e-9 ) then
!          write (950, '(i8, i8, 2(f24.8))') pBand, nBand, 0.0
!        else
!          write (950, '(i8, i8, 2(f24.8))') pBand, nBand, aimag(sJuPhon(pBand, nBand))
!        end if
!        if ( abs(real(hJuPhon(pBand, nBand))) < 1e-9 ) then
!          write (960, '(i8, i8, 2(f24.8))') pBand, nBand, 0.0
!        else
!          write (960, '(i8, i8, 2(f24.8))') pBand, nBand, real(hJuPhon(pBand, nBand))
!        end if
!        if ( abs(aimag(hJuPhon(pBand, nBand))) < 1e-9 ) then
!          write (960, '(i8, i8, 2(f24.8))') pBand, nBand, 0.0
!        else
!          write (960, '(i8, i8, 2(f24.8))') pBand, nBand, aimag(hJuPhon(pBand, nBand))
!        end if
!!       ! write (*, *) 'bari'
!!       ! write (*, *) s(pBand, nBand)
!        write (951, '(i8, i8, 2(f24.8))') pBand, nBand, hJuPhon(pBand, nBand)
!        write (952, '(i8, i8, 2(f24.8))') pBand, nBand, sJuPhon(pBand, nBand)
!      end do
!    end do

    ! Load Hamiltonian matrix and overlap for certain k-point !todo which k-point?
    allocate(hFleur(dimens%neigd, dimens%neigd), sFleur(dimens%neigd, dimens%neigd))
    call fopen( 1000, name='hmatTest', status='old', action='read', form='unformatted' )
    read(1000) hFleur
    call fclose(1000)
    ! Load overlap matrix and overlap for certain k-point !todo which k-point?
    call fopen( 1000, name='smatTest', status='old', action='read', form='unformatted' )
    read(1000) sFleur
    call fclose(1000)
    ! Transform it into Kohn--Sham wavefunction space
    ! Compare it !todo really for every band combination?
    passedS = .true.
    passedH = .true.
    do pBand = 1, ne(ikpt)
      do nBand = 1, pband
        if ( abs( hFleur(pBand, nBand) - conjg(hJuPhon(pBand, nBand)) ) > 9e-6 ) then
          write (logUnit, '(a,i4,a,i4, 2(es15.8))') 'Problems with Hamiltonian at band indices ', pBand, 'and', nBand, hFleur(pBand, nBand) - hJuPhon(pBand, nBand)
          passedH = passedH .and. .false.
        end if
        if ( abs( sFleur(pBand, nBand) - conjg(sJuPhon(pBand, nBand)) ) > 9e-6 ) then
          write (logUnit, '(a,i4,a,i4,2(es15.8))') 'Problems with overlap matrix at band indices ', pBand, 'and', nBand, hFleur(pBand, nBand) - hJuPhon(pBand, nBand)
          passedS = passedS .and. .false.
        end if
      end do
    end do
    if ( .not.passedH ) then
      call juDFT_warn('Test of Hamilton matrix elements within Sternheimer equation failed', calledby='TestHSME', hint='Check log file .')
    else
      write (logUnit, '(a)') 'Test comparison of muffin-tin Hamilton matrix in Sternheimer equation and Fleur with 1e-7 accuracy!'
      write (logUnit, '(a)') '---------------------------------------------------------------------------------------------------'
      write (logUnit, '(a)') '                                                                                                  |__&
                                                                                                                            &passed'
      write (logUnit, *)
    end if
    if ( .not.passedS ) then
      call juDFT_warn('Test of Overlap matrix elements within Sternheimer equation failed', calledby='TestHSME', hint='Check log file.')
    else
      write (logUnit, '(a)') 'Test comparison of muffin-tin overlap matrix in Sternheimer equation and Fleur with 1e-7 accuracy!'
      write (logUnit, '(a)') '--------------------------------------------------------------------------------------------------'
      write (logUnit, '(a)') '                                                                                                 |__p&
                                                                                                                             &assed'
      write (logUnit, *)
    end if

  end subroutine TestHSME

  !--------------------------------------------------------------------------------------------------------------------------------------
  !> @author
  !> Christian-Roman Gerhorst, Gregor Michaliczek
  !>
  !> @brief
  !> Test of warping the IR potential before putting it into the Sternheimer equation.
  !>
  !> @details
  !> Test warping routine in Sternheimer module by reading in warped and unwarped IR potential and warping unwarped IR potential and
  !> compare it to the read-in warped potential.
  !>
  !> @note
  !> The dimensions stars%k1d, stars%k2d and stars%k3d have to be chosen big enough so that the potentials can be compared with the Fleur
  !> potentials. In Fleur an averaging is being made when remapping the Fourier mesh to the stars again. If the mesh is mapped to a
  !> G-vector representation the averaging is omitted. Therefore aliasing errors weigh heavier! Alising errors are reduced when enlarging
  !> the FFT-mesh within the fl7para file (i.e. k1d, k2d, k3d). One can determine the required dimensions with the Nyquist-Shannon
  !> sampling theorem to optimize performance later that is to only enlarge as much as required. Aliasing errors cannot be balanced out
  !> they have to be avoided! Calculations with forces have given indications that a higher Gmax-cutoff might also solve the problem in
  !> a way so that the aliasing error does not weigh as heavy as before.
  !>
  !> @attention
  !> This routine can require manipulation of k1d, k2d, k3d within fl7para file! This should be taken care of when integrating in Fleur.
  !--------------------------------------------------------------------------------------------------------------------------------------
  subroutine TestWarping(stars, ngdp, logUnit, gdp)

    use m_types
    use m_juDFT_NOstopNO, only : juDFT_warn
    use mod_juPhonUtils, only : fopen, fclose
    use m_jpPotDensHelper, only : WarpIRPot, convertStar2G

    implicit none

    ! Type parameters
    type(t_stars),             intent(in) :: stars

    ! Scalar parameters
    integer,                   intent(in) :: ngdp
    integer,                   intent(in) :: logUnit

    ! Array parameters
    integer,                   intent(in) :: gdp(:, :)

    ! Local array variables
    complex,      allocatable             :: vpw_eff(:, :)
    complex,      allocatable             :: vpw_eff_w(:, :)
    complex,      allocatable             :: vpwG(:)
    complex,      allocatable             :: vpwG3D(:, :)
    complex,      allocatable             :: vpw_wFleur(:, :)
    complex,      allocatable             :: vpw_wFleurG( :)
    complex,      allocatable             :: vpw_wGjP(:, :)
    real,         allocatable             :: ufft(:)
    
integer :: ii

    integer :: i
    write ( logUnit, '(a)' ) 'Test of Sternheimer warping routine with 1e-6 accuracy'
    write ( logUnit, '(a)' ) '------------------------------------------------------'

    ! Load warped and unwarped potential from Fleur
    ! todo WARNING: This is not the effective but only the Coulomb potential
    allocate( vpw_eff( stars%ng3, 1 ), vpw_eff_w( stars%ng3, 1 ) )
    vpw_eff = 0
    vpw_eff_w = 0
    call Fopen( 1000, name='vpw_wt', status='old', form='unformatted' )
    read( 1000 ) vpw_eff
    call Fclose( 1000 )

    call Fopen( 1000, name='vpw_wt_w', status='old', form='unformatted' )
    read( 1000 ) vpw_eff_w
    call Fclose( 1000 )



    ! Transform Fleur potentials into G-vector representation
    allocate( vpwG( ngdp ), vpwG3D( ngdp, 3 ), vpw_wFleurG( ngdp ) )
    vpwG = 0
    vpwG3D = 0
    vpw_wFleurG = 0
    call convertStar2G( vpw_eff(:, 1), vpwG, stars, ngdp, gdp )
    vpwG3D (:, 1) = vpwG
    vpwG3D (:, 2) = vpwG
    vpwG3D (:, 3) = vpwG

    call convertStar2G( vpw_eff_w(:, 1), vpw_wFleurG, stars, ngdp, gdp )
  !  open (1000, file='veffG', status='replace', form='formatted')
  !  do ii = 1, ngdp
  !    write (1000, *) stars%ig(gdp(1, ii), gdp(2, ii), gdp(3, ii)), vpwG3D(ii, 1)
  !  end do
  !  close(1000)

  !  open (1000, file='veffWG', status='replace', form='formatted')
  !  do ii = 1, ngdp
  !    write (1000, *) stars%ig(gdp(1, ii), gdp(2, ii), gdp(3, ii)), vpw_wFleurG(ii)
  !  end do
  !  close(1000)
    ! Warp unwarped potential of FLEUR
    allocate( vpw_wGjP( ngdp, 1 ) )
    vpw_wGjP = 0
    !todo this can be wrong now, without direction from outside it worked, the louter llop over directions only went from 1 to 1
    call warpIRPot( stars, ngdp, 1, gdp, vpwG3D, vpw_wGjP(:, 1) )

!    write (11111, *) vpw_wGjP(:, 1)
!    write (22222, *) vpw_wFleurG(:)
!    write (33333, *) vpw_wGjP(:, 1) - vpw_wFleurG
!    do ii = 1, ngdp
!      write (44444, '(i5, 3(i3), 5(f20.11))') ii, gdp(:, ii), vpw_wGjP(ii, 1), vpw_wFleurG(ii), abs(vpw_wGjP(ii, 1)- vpw_wFleurG(ii))
!    end do
    ! Compare the result to warped potential read in from Fleur
    if ( any( abs( vpw_wGjP(:, 1) - vpw_wFleurG(:) ) >= 1e-6 ) ) then
      write ( logUnit, '(a)' ) '                                   |__failed'
      call JuDFT_warn( 'Sternheimer warping routine for interstitial potential inconsistent', calledby='TestWarping', &
        & hint='Increase the size of the FFT box k1d, k2d, k3d within fl7para file!' )
    else
      write (logUnit, '(a)')   '                                                     |__passed'
    end if
    write (logUnit, *)
  end subroutine TestWarping

  subroutine Testz1ikGz0q0Diff( atoms, dimens, kpts, cell, nv, nobd, gbas, mapGbas, kveclo, z0)

    ! move this routine to readInz1
    use m_jpSetupDynMatHelper, only : readInz1
    use m_types
    use m_jpConstants, only : iu

    implicit none

    ! Type parameter
    type(t_atoms),                  intent(in)  :: atoms
    type(t_dimension),              intent(in)  :: dimens
    type(t_kpts),                   intent(in)  :: kpts
    type(t_cell),                   intent(in)  :: cell

    ! Array parameter
    integer,                        intent(in)  :: nv(:, :)
    integer,                        intent(in)  :: nobd(:, :)
    integer,                        intent(in)  :: gbas(:, :)
    integer,                        intent(in)  :: mapGbas(:, :, :)
    integer,                        intent(in)  :: kveclo(:, :)
    MCOMPLEX,                       intent(in)  :: z0(:,:,:,:)

    ! Scalar variables
    integer                                     :: ikpt
    integer                                     :: iband
    integer                                     :: iDtype
    integer                                     :: iDeqat
    integer                                     :: iDatom
    integer                                     :: idir
    integer                                     :: iBas
    real                                        :: invsAtNr
    complex                                     :: z1Bench


    ! Array variables
    complex,           allocatable              :: z1nG(:, :, :, :)
    real                                        :: kExt(3)
    real                                        :: GExt(3)

    !todo for k /= 0 this test still does not work
    ! but somehow if we multiply it with the other z and sum up the density comes out

    allocate( z1nG(dimens%nvd, 3, atoms%nat, maxval(nobd(:, :))) )
    z1nG = cmplx(0., 0.)

    invsAtNr = 1. / real(atoms%nat)

    do ikpt = 1, 5!kpts%nkpt
      z1nG = cmplx(0., 0.)
      call readInz1( atoms, ikpt, 1, ikpt, nobd, nv, z1nG )
      kExt(1:3) = matmul( cell%bmat(1:3, 1:3), kpts%bk(1:3, ikpt) )

      do iband = 1, nobd(ikpt, 1)
        iDatom = 0
        do iDtype = 1, atoms%ntype
          do iDeqat = 1, atoms%neq(iDtype)
            iDatom = iDatom + 1
            do idir = 1, 3
              do iBas = 1, nv(1, ikpt) + atoms%nlotot
                Gext(1:3) = matmul( cell%bmat(1:3, 1:3), gbas(1:3, iBas) )
                z1Bench = - iu * cmplx(invsAtNr, 0) * ( Gext(idir) + kExt(idir) ) * z0(iBas, iband, ikpt, 1)
                write(2002, '(7(i8),2(f15.8))') iDatom, iband, idir, ikpt, gbas(1, iBas), gbas(2, iBas), gbas(3, iBas), z1Bench
                write(2003, '(7(i8),2(f15.8))') iDatom, iband, idir, ikpt, gbas(1, iBas), gbas(2, iBas), gbas(3, iBas), z1nG(iBas, idir, iDatom, iband)
              end do ! iBas
              do iBas = nv(1, ikpt) + 1, nv(1, ikpt) + atoms%nlotot
                Gext(1:3) = matmul( cell%bmat(1:3, 1:3), gbas(1:3, mapGbas(kveclo(iBas - nv(1, ikpt), ikpt), ikpt, 1)) )
                z1Bench =  - iu * cmplx(invsAtNr, 0) * ( Gext(idir) + kExt(idir)) * z0(iBas, iband, ikpt, 1)
                write(2002, '(7(i8),4(f15.8))') iDatom, iband, idir, ikpt, gbas(1, mapGbas(kveclo(iBas - nv(1, ikpt), ikpt), ikpt, 1)),   &
                  & gbas(2, mapGbas(kveclo(iBas - nv(1, ikpt), ikpt), ikpt, 1)),                                                    &
                  & gbas(3, mapGbas(kveclo(iBas - nv(1, ikpt), ikpt), ikpt, 1)), z1Bench
                write(2003, '(7(i8),4(f15.8))') iDatom, iband, idir, ikpt, gbas(1, mapGbas(kveclo(iBas - nv(1, ikpt), ikpt), ikpt, 1)),   &
                  & gbas(2, mapGbas(kveclo(iBas - nv(1, ikpt), ikpt), ikpt, 1)),                                                    &
                  & gbas(3, mapGbas(kveclo(iBas - nv(1, ikpt), ikpt), ikpt, 1)), z1nG(iBas, idir, iDatom, iband)
              end do ! iBas
            end do ! idir
          end do ! iDeqat
        end do ! iDtype
      end do ! iband
    end do ! ikpt
  end subroutine Testz1ikGz0q0Diff

  subroutine testInstq0Conv(atoms, stars, cell, kpts, qpts, dimens, input, results, lathar, sym, enpara, uds, V0Fleur, ngdp, nobd, kveclo,&
      & memd_atom, rbas1, rbas2, uuilon, duilon, ilo2p, ulouilopn, z0, nv, gbas, mapGbas, gdp2iLim, kpq2kPrVec, mapKpq2K, gdp2Ind, &
      & clnu_atom, mlh_atom, nmem_atom, gdp, rho0IR, rho0MT, ne, eig, El, nRadFun, iloTable, noPtsCon, testrho0ContSw, logUnit, vEff0MTsh, vEff0pwUw)

    use m_jpConstants, only : iu
    use m_jpDens1stVar, only : calcRho1IRValDS, calcPsDensMT, calcVarCTcorr, calcPDCinAlph
    use m_jpPotDensHelper, only : WarpIRPot,  convertStar2G
    use mod_juPhonUtils, only :fopen, fclose, calcGrR2FinLH
    use m_jpPotDensHelper, only : calcIRdVxcKern, calcMTdVxcKern
    use m_jpGrVeff0, only : GenGrVeff0
    use m_jpSternheimer, only : solveSternheimerEq, GenRecEdifME
    use m_jpIOnMixing, only : storeZ1nG
    use m_jpSternhHF, only : tlmplm4V
    use m_jpSternhPulaySurface, only : tlmplm4H0, IRcoeffVeffUV
    use m_jpTestInput, only : checkDensnPot
    use m_jpTestPotential, only : checkjuPhPots
    use m_jpVeff1
    use m_types
    use m_jpSternhPulaySurface, only : calcSfVeffFast, calcSintKinEps
    use m_jpPlotObservables, only : plot1stDensVar
    use m_jpIOnMixing, only : storeDensity

    implicit none

    ! Type parameters
    type(t_atoms),                  intent(in) :: atoms
    type(t_cell),                   intent(in) :: cell
    type(t_kpts),                   intent(in) :: kpts
    type(t_kpts),                   intent(in) :: qpts
    type(t_dimension),              intent(in) :: dimens
    type(t_input),                  intent(in) :: input
    type(t_results),                intent(in) :: results
    type(t_sphhar),                 intent(in) :: lathar
    type(t_stars),                  intent(in) :: stars
    type(t_sym),                    intent(in) :: sym
    type(t_enpara),                 intent(in) :: enpara
    type(t_usdus),                  intent(in) :: uds
    type(t_potential),              intent(in) :: V0Fleur


    ! Scalar parameters
    integer,                        intent(in) :: ngdp
    integer,                        intent(in) :: memd_atom
    integer,                        intent(in) :: noPtsCon
    logical,                        intent(in) :: testrho0ContSw
    integer,                        intent(in) :: logUnit

    ! Array parameters
    integer,                        intent(in) :: mlh_atom(:, 0:, :)
    integer,                        intent(in) :: nmem_atom(0:, :)
    complex,                        intent(in) :: clnu_atom(:, 0:, :)
    integer,                        intent(in) :: nobd(:, :)
    integer,                        intent(in) :: kveclo(:, :)
    MCOMPLEX,                       intent(in) :: z0(:, :, :, :)
    integer,                        intent(in) :: nv(:, :)
    integer,                        intent(in) :: kpq2kPrVec(:, :, :)
    integer,                        intent(in) :: mapGbas(:, :, :)
    integer,                        intent(in) :: mapKpq2K(:, :)
    integer,                        intent(in) :: gdp2Ind(:, :, :)
    integer,                        intent(in) :: gbas(:, :)
    complex,                        intent(in) :: rho0IR(:, :)
    integer,                        intent(in) :: gdp(:, :)
    real,                           intent(in) :: rho0MT(:, 0:, :, :)
    real,                           intent(in) :: rbas1(:, :, 0:, :, :)
    real,                           intent(in) :: rbas2(:, :, 0:, :, :)
    real,                           intent(in) :: uuilon(:, :)
    real,                           intent(in) :: duilon(:, :)
    real,                           intent(in) :: ulouilopn(:, :, :)
    integer,                        intent(in) :: ilo2p(:, :)
    integer,                        intent(in) :: ne(:)
    real,                           intent(in) :: eig(:,:,:)
    real,                           intent(in) :: El(:, 0:, :, :)
    integer,                        intent(in) :: nRadFun(0:, :)
    integer,                        intent(in) :: iloTable(:, 0:, :)
    complex,                        intent(in) :: vEff0MTsh(:, :, :, :)
    integer,                        intent(in) :: gdp2iLim(2, 3)
    complex,                        intent(in) :: vEff0pwUw(:, :)

    ! Type variables
    type(t_tlmplm)                             :: tdHS0
    type(t_tlmplm)                             :: tdVx
    type(t_tlmplm)                             :: tdVy
    type(t_tlmplm)                             :: tdVz
    type(t_tlmplm)                             :: tdVx2
    type(t_tlmplm)                             :: tdVy2
    type(t_tlmplm)                             :: tdVz2

    ! Scalar variables
    integer                                    :: iDtype
    integer                                    :: iDatom
    integer                                    :: iDeqat
    integer                                    :: ikpt
    integer                                    :: iG
    integer                                    :: iqpt
    integer                                    :: idir
    integer                                    :: iBas
    integer                                    :: iband
    integer                                    :: iatom
    integer                                    :: itype
    integer                                    :: ieqat
    integer                                    :: ptsym
    integer                                    :: ilh
    integer                                    :: oqn_l
    integer                                    :: lm_pre
    integer                                    :: imem
    integer                                    :: lm
    integer                                    :: lm2
    integer                                    :: imesh
    integer                                    :: mqn_m
    real                                       :: invsAtNr
    logical                                    :: harSw
    logical                                    :: extSw
    logical                                    :: xcSw
    logical                                    :: testGoldstein
    logical                                    :: vExtFull
    integer                                    :: ispin
    logical                                    :: densityNotPotential
    integer                                    :: lmpMax
    logical                                    :: grRhoTermSw
    logical                                    :: vHarNum
    integer                                    :: maxlmp

    ! Array variables
    complex,           allocatable             :: z1Special(:, :, :, :)
    complex,           allocatable             :: vEff1IR(:, :)
    complex,           allocatable             :: vEff1MT(:, :, :, :)
    complex,           allocatable             :: rho1IR(:, :, :)
    complex,           allocatable             :: rho1IRReset(:, :, :)
    complex,           allocatable             :: rho1IRAfter(:, :, :)
    complex,           allocatable             :: rho1MT(:, :, :, :, :)
    complex,           allocatable             :: grRho0MT(:, :, :, :)
    real,              allocatable             :: r2Rho0MT(:, :, :, :)
    complex,           allocatable             :: rho1MTCoreDispAt(:, :, :, :)
    complex,           allocatable             :: rho1MTDummy(:, :, :, :)
    complex,           allocatable             :: r2GrRho0MT(:, :, :, :)
    complex,           allocatable             :: qpwcG(:, :)
    real,              allocatable             :: acoff(:)
    real,              allocatable             :: alpha(:)
    complex,           allocatable             :: rho0IRPW(:, :)
    complex,           allocatable             :: rho1IRctC(:, :, :)
    complex,           allocatable             :: rho1MTctC(:, :, :, :, :)
    complex,           allocatable             :: grVxcIRKern(:)
    real,              allocatable             :: dKernMTGPts(:, :, :)
    complex,           allocatable             :: ylm(:, : )
    real,              allocatable             :: gausWts(:) ! gaussian weights belonging to gausPts
    complex,           allocatable             :: grVeff0MT(:, :, :, :)
    complex,           allocatable             :: grRho0MTDummy(:, :, :, :) ! Dummy quantity at the moment
    complex,           allocatable             :: grVeff0IR(:, :)
    complex,           allocatable             :: sumVMTs(:, :, :, :)
    complex,           allocatable             :: sumVMTs2(:, :, :, :)
    complex,           allocatable             :: w_vEff1IR(:, :)
    integer,           allocatable             :: nlo_atom(:)
    real,              allocatable             :: recEdiffME(:, :)
    real,              allocatable             :: cutContr(:, :)
    complex,           allocatable             :: z1nG(:, :, :)
!    complex,           allocatable             :: vEff0IRG( :)
    complex,           allocatable             :: veffUvIR(:, :)
    complex,           allocatable             :: surfIntVFast(:, :, :)
    complex,           allocatable             :: rho0MTsh(:, :, :, :)
    complex,           allocatable             :: grRho0IR(:, :)

    complex,   allocatable                      :: hFullNoAbcofBK(:, :, :, :)
    complex,   allocatable                      :: overlapNoAbcofBK(:, :, :, :)

    real                                       :: kExt(3)
    real                                       :: Gext(3)

    integer :: ii, jj

    complex,           allocatable              :: zDummy(:, :, :, :)
    complex, allocatable :: vEff0IRG_w(:)

    complex              :: vpwStar(stars%n3d, 1)
    complex              :: vpwStar1(stars%n3d, 1)
    complex              :: vpwStarDiff(stars%n3d, 1)
    complex              :: vpwStarDiffG(ngdp)
    character(len=22)                            :: filename

    write(*, *)
    write(*, *)
    write(*, *)
    write(*, *)

    allocate( z1Special(dimens%nbasfcn, dimens%neigd, 3, kpts%nkpt) )
    allocate( rho1IR( ngdp, 3, atoms%nat ) )
    allocate( rho1IRAfter( ngdp, 3, atoms%nat ) )
    allocate( rho1IRReset( ngdp, 3, atoms%nat ) )
    allocate( rho1MT( atoms%jmtd, (atoms%lmaxd + 1)**2, atoms%nat, 3, atoms%nat ) )
    allocate( sumVMTs(atoms%jmtd, (atoms%lmaxd + 1 )**2, 3, atoms%nat) )
    allocate( sumVMTs2(atoms%jmtd, (atoms%lmaxd + 1 )**2, 3, atoms%nat) )

    write(*, *) 'remove zDummy'
    allocate( zDummy(dimens%nbasfcn, dimens%neigd, kpts%nkpt, 1) )
    zDummy = cmplx(0., 0.)
!    zDummy = cmplx(real(z0), 0.)
    zDummy = z0

!    do ikpt = 1, kpts%nkpt
!      write(*, *) z0(:, :, ikpt, 1)
!      read(*, *)
!    end do ! ikpt
!    NOstopNO

    iqpt = 1
    ispin = 1
    rho1IR( :, :, : ) = cmplx(0., 0.)
    rho1IRAfter( :, :, : ) = cmplx(0., 0.)
    rho1IRReset( :, :, : ) = cmplx(0., 0.)
    rho1MT = cmplx(0., 0.)

    ! todo If we loop over displaced atoms, take care of rho1MT and where the gradient is subtracted
    iDatom = 1
    iDtype = 1
    invsAtNr = 1.! / atoms%nat
          z1Special(:, :, :, :) = cmplx(0.0, 0.0)
        ! Set up Sternheimer q = 0 solution artificially and calculate the interstitial potential
        do ikpt = 1, kpts%nkpt

          kExt(1:3) = matmul( cell%bmat(1:3, 1:3), kpts%bk(1:3, ikpt) )
          do idir = 1, 3
            do iBas = 1, nv(1, ikpt)
              Gext(1:3) = matmul( cell%bmat(1:3, 1:3), gbas(1:3, mapGbas(iBas, ikpt, 1)) )
              do iband = 1, nobd(ikpt, 1) !todo LOs mit berücksichtigen, auch mit Gs?, + or - q?
                z1Special(iBas, iband, idir, ikpt) = z1Special(iBas, iband, idir, ikpt) - iu * cmplx(invsAtNr, 0) * ( Gext(idir) + kExt(idir) )&
                  & * zDummy(iBas, iband, ikpt, 1)

                !write(2007, '(4(i8),2(f15.8))') ikpt, idir, iBas, iband, z1Special(iBas, iband, idir)
              end do ! iband
            end do ! iBas
            do iBas = nv(1, ikpt) + 1, nv(1, ikpt) + atoms%nlotot
              Gext(1:3) = matmul( cell%bmat(1:3, 1:3), gbas(1:3, mapGbas(kveclo(iBas - nv(1, ikpt), ikpt), ikpt, 1)) )
              do iband = 1, 1!nobd(ikpt, 1) !todo LOs mit berücksichtigen, auch mit Gs?, + or - q?
                z1Special(iBas, iband, idir, ikpt) = z1Special(iBas, iband, idir, ikpt) - iu * cmplx(invsAtNr, 0) * ( Gext(idir) + kExt(idir) )&
                  & * zDummy(iBas, iband, ikpt, 1)
                write(2007, '(4(i8),2(f15.8))') ikpt, idir, iBas, iband, aimag(z1Special(iBas, iband, idir, ikpt))
              end do ! iband
            end do ! iBas
            ! calculate interstitial density with the special z1
            call calcRho1IRValDS( cell, results, nobd(ikpt, 1), nv(1, :), ikpt, iqpt, mapKpq2K(ikpt,iqpt), idir, gbas(:, :),&
              & zDummy(:, :, ikpt, 1), z1Special(:, :, :, ikpt), rho1IRReset(:, :, iDatom), gdp2Ind, mapGbas, gdp2iLim, kpq2kPrVec)
            call calcRho1IRValDS( cell, results, nobd(ikpt, 1), nv(1, :), ikpt, iqpt, mapKpq2K(ikpt,iqpt), idir, gbas(:, :),&
              & zDummy(:, :, ikpt, 1), z1Special(:, :, :, ikpt), rho1IR(:, :, iDatom), gdp2Ind, mapGbas, gdp2iLim, kpq2kPrVec)
          end do ! idir

      !  if (ikpt ==56) then
      !!    call plot1stDensVar(sym, cell, input, lathar, atoms, rho1IRReset(:, :, iDatom), ngdp, gdp, mlh_atom, nmem_atom, clnu_atom, rho0MT, qpts%bk(:, iqpt), rho1MT(:, :, :, :, iDatom), 1., 1., 1., iDatom, 1, 1) !todo remove counter!
      !    do idir = 1, 3
      !      do iBas = 1, nv(1, ikpt)
      !        do iband = 1, nobd(ikpt, 1)
      !          write(2007, '(4(i8),2(f15.8))') ikpt, idir, iBas, iband, z1Special(iBas, iband, idir, ikpt)
      !        end do ! iband
      !      end do ! iBas
      !    end do ! idir
       !   do idir = 1, 3
       !     do iG = 1, ngdp
       !       write(2098, '(4(i8),2(f15.8))') iDatom, ikpt, idir, iG, rho1IR(iG, idir, iDatom)
       !     end do ! iG
       !   end do ! idir
       ! endif
        rho1IRReset = cmplx(0., 0.)
      end do ! ikpt


        ! Add core contribution
        call calcPsDensMT( atoms, cell, sym, stars, dimens, input, ngdp, acoff, alpha, qpwcG, gdp, logUnit )

        allocate( rho1IRctC(ngdp, atoms%nat, 3), rho1MTctC(atoms%jmtd, (atoms%lmaxd + 1)**2, atoms%nat, 3, atoms%nat) )
        rho1IRctC = cmplx(0.0, 0.0)
        rho1MTctC = cmplx(0.0, 0.0)

        ! For neon we have no core contribution
        call calcVarCTcorr( atoms, cell, ngdp, gdp(:, :ngdp), [0., 0., 0.], qpwcG, rho1IRctC, rho1MTctC )
        do idir = 1, 3
          do iG = 1,ngdp
            ! IR core contribution correction to the density variation
            ! We have a minus in the calculation of rho1IRctC
            rho1IR(iG, idir, iDatom) = rho1IR(iG, idir, iDatom) + rho1IRctC(iG, iDatom, idir)
            write(2122, '(2(i8),1(f15.8))') idir, iG, aimag(rho1IR(iG, idir, iDatom))
          end do
        end do

        call plot1stDensVar(sym, cell, input, lathar, atoms, rho1IR(:, :, iDatom), ngdp, gdp, mlh_atom, nmem_atom, clnu_atom, rho0MT, qpts%bk(:, iqpt), rho1MT(:, :, :, :, iDatom), 1., 1., 1., iDatom, 1, 1) !todo remove counter!
        do idir = 1, 3
          write(filename, '(a10,i1,a4,i1,a4,i2)') 'JPcdn1_Dat', iDatom, 'Ddir', idir, 'qInd', iqpt
          call storeDensity(atoms, ngdp, filename, rho1IR(:, idir, iDatom), rho1MT(:, :, :, idir, iDatom))
        end do ! idir

        ! Check if interstitial density is correct
        allocate( rho0IRpw(ngdp, 1) )
        rho0IRpw(:, :) = cmplx(.0, .0)
        call convertStar2G( rho0IR(:, 1), rho0IRpw(:, 1), stars, ngdp, gdp )

        do iG = 1, ngdp
          Gext(1:3) = matmul( cell%bmat(1:3, 1:3), gdp(1:3, iG) )
          do idir = 1, 3
            !write(2000, '(i8,i8,2(f15.8))') iG, idir, -iu * Gext(idir) * rho0IRpw(iG)
            !write(2001, '(i8,i8,2(f15.8))') iG, idir, rho1IR(iG, idir, iDatom)
            if ( abs( -iu * Gext(idir) * rho0IRpw(iG, 1) - rho1IR(iG, idir, iDatom) ) > 1e-9 ) then
              write(*, '(a,i8,i8)') 'q=0 Sternheimer IR density broken at G and direction', iG, idir
            end if
          end do ! idir
        end do ! iG


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
                  end do ! imesh
                end do ! mqn_m
              end do ! oqn_l
            end do ! ieqat
          end do ! itype
        end do ! idir

        ! We have checked whether all aother terms are zero in another place, todo do this also in productio later
        ! Numerical stable muffin-tin gradient routine
        rho1MT = cmplx(0., 0.)
        do idir = 1, 3
          do oqn_l = 0, atoms%lmax(iDtype) !+ 1
            do mqn_m = -oqn_l, oqn_l
              lm = oqn_l * (oqn_l + 1) + mqn_m + 1! todo does everythink start at 1?
              do imesh = 1, atoms%jri(iDtype)
                ! If we are at the displaced atom we have to add the gradient of the density
                rho1MT(imesh, lm, iDatom, idir, iDatom) =  rho1MT(imesh, lm, iDatom, idir, iDatom) &
                  &- grRho0MT(imesh, lm, iDatom, idir)
!                write(2004, '(2(f15.8))') rho1MT(imesh, lm, iDatom, idir, iDatom)
              end do ! imesh
            end do ! mqn_m
          end do ! oqn_l
        end do ! idir

        ! Coretail corrections are zero because there are no coretails reaching out of the muffin-tin which could be corrected
        allocate( rho1MTCoreDispAt(atoms%jmtd, 1 : 4, atoms%nat, 3) )
        rho1MTCoreDispAt = cmplx(0.0, 0.0)
        call calcPDCinAlph( atoms, alpha, acoff, rho1MTCoreDispAt )
        iatom = 0
        do itype = 1, atoms%ntype
          do ieqat = 1, atoms%neq(itype)
            iatom = iatom + 1
            do idir = 1, 3
              if ( iDatom == iatom ) then
                do lm = 2, 4
                  do imesh = 1, atoms%jri(itype)
                    ! The correction of the core correction at the displaced atom. We need not to correct the core contribution of
                    ! the displaced atom as it is the reason for the core corrections in the non-displaced atoms
!                    rho1MT( imesh, lm, iatom, idir, iDatom) = rho1MT( imesh, lm, iatom, idir, iDatom ) -                           &
!                      & rho1MTCoreDispAt( imesh, lm, iDatom, idir )
                    !write(2005, '(2(f15.8))') rho1MTCoreDispAt(imesh, lm, iDatom, idir)
                  end do
                end do
              end if
              do oqn_l = 0, atoms%lmax(itype)
                do mqn_m = -oqn_l, oqn_l
                  lm = oqn_l * (oqn_l + 1) + mqn_m + 1
                  do imesh = 1, atoms%jri(itype)
                    ! MT core contribution correction to the density variation
                    ! we have a minus in the calculation of rho1MTctC
!                    rho1MT(imesh, lm, iatom, idir, iDatom) = rho1MT(imesh, lm, iatom, idir, iDatom) +                              &
!                                                                                        & rho1MTctC( imesh, lm, iatom, idir, iDatom)
                    ! it is okay that the coretail correction is zero because there is no electron loss from the core for Neon!
                    !write(2006, '(2(f15.8))') rho1MTctC(imesh, lm, iatom, idir, iDatom)
                  end do ! imesh
                end do ! mqn_m
              end do ! oqn_l
            end do ! idir
          end do ! ieqat
        end do ! itype

        if (.true.) then
          write(logUnit, '(a)') 'Check and compare continuity of electron density and density variation.'
          ! Validate the continuity of the density and the density variation
          if (.not.testrho0ContSw) then
            call CheckDensnPot( atoms, sym, stars, cell, lathar, noPtsCon, logUnit, densityNotPotential, ispin, rho0IR, rho0MT)
          else
            write(logUnit, '(a)') 'Continuity of density already checked. Please, scroll up!'
          end if
          !todo this dummy variable is necessary because the order of indices is not the same in potential and density
          allocate(rho1MTDummy(atoms%jmtd, (atoms%lmaxd + 1)**2, 3, atoms%nat))
          rho1MTDummy = cmplx(0., 0.)
          iatom = 0
          do itype = 1, atoms%ntype
            do ieqat = 1, atoms%neq(itype)
              iatom = iatom + 1
              do idir = 1, 3
                do oqn_l = 0, atoms%lmax(itype)
                  do mqn_m = -oqn_l, oqn_l
                    lm = oqn_l * (oqn_l + 1) + 1 + mqn_m
                    do imesh = 1, atoms%jri(itype)
                      rho1MTDummy(imesh, lm, idir, iatom) = rho1MT(imesh, lm, iatom, idir, iDatom)
                    end do ! imesh
                  end do ! mqn_m
                end do ! oqn_l
              end do ! idir
            end do ! ieqat
          end do ! itype

          ! This at the currecnt stage only expands until lmax and not lmax + 1
          call checkjuPhPots(noPtsCon, atoms, ngdp, gdp, cell, lathar, sym, rho1IR(:, :, iDatom), rho1MTDummy, logUnit)
        end if

        allocate( grRho0IR(ngdp, 3) )
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

        call calcIRdVxcKern(stars, gdp, ngdp, rho0IR(:, 1), grVxcIRKern) ! add spin for this and next line
        call calcMTdVxcKern(atoms, dimens, lathar, rho0MT(:, :, :, 1), nmem_atom, clnu_atom, mlh_atom, gausWts, ylm,        &
          & dKernMTGPts)

        harSw = .true.
        extSw = .true.
        xcSw = .true.
        testGoldstein = .false.
        grRhoTermSw = .false. ! this has to be equal to vExtFull going into Veff1
        call GenGrVeff0( atoms, cell, dimens, stars, ngdp, harSw, extSw, xcSw, gdp, rho0IRpw, rho0MTsh, grRho0IR, grRho0MT, gausWts, ylm,&
          & dKernMTGPts, grVxcIRKern, testGoldstein, grRhoTermSw, grVeff0IR, grVeff0MT ) ! add spin

        if (.false.) then
          do idir = 1, 3
            do iG = 1, ngdp
              write(2011, '(2(i8),2(f15.8))') idir, iG, -grVeff0IR(iG, idir)
            end do ! iG
          end do
        end if

        if (.false.) then
          ! vpw_eff and vpw_coul have to be present for this to work!
          ! The fluctuations of the interstitial potential only stem from calculating the hartree pseudo potential with an unprecise
          ! interal method. The xc and the external potential deliver the exact analytical gradient However this has no influence yet
          ! on the matrix element used in Sternheimer for Neon
          vpwStar = cmplx(0., 0.)
          vpwStar1 = cmplx(0., 0.)
          vpwStarDiff= cmplx(0., 0.)
          call fopen( 1000, name='vpw_eff', status='old', action='read', form='unformatted' )
            read( 1000 ) vpwStar
          call fclose( 1000 )
          call fopen( 1000, name='vpw_coul', status='old', action='read', form='unformatted' )
            read( 1000 ) vpwStar1
          call fclose( 1000 )
          vpwStarDiff = vpwStar1 - vpwStar
          vpwStarDiffG = cmplx(0., 0.)
          call convertStar2G(vpwStarDiff(:, 1), vpwStarDiffG, stars, ngdp, gdp)
          do idir = 1, 3
            do iG = 1, ngdp
              Gext(1:3) = matmul(cell%bmat(1:3, 1:3), gdp(1:3, iG))
              write(2094, '(2(i8),2(f15.8))') idir, iG, -iu * Gext(idir) * vpwStarDiffG(iG)
            end do
          end do ! idir
        end if


        harSw = .true.
        extSw = .true.
        xcSw = .true.
        vExtFull = .false.
        vHarNum = .false.
        ! The density only consists of - the gradient. Within the routine genveff1 the gradient is added so to be compatible with
        ! the Sternheimer production routine we do not insert it here
        rho1MT(:, :, :, :, :) = cmplx(0., 0.)
        call GenVeff1( stars, cell, atoms, dimens, harSw, extSw, xcSw, VextFull, ngdp, [0., 0., 0.], rho0IRpw, rho0MTsh, &
          & rho1IR(:, :, iDatom), rho1MT(:, :, :, :, iDatom), grRho0MT, gdp, vEff1IR, vEff1MT, grVxcIRKern, ylm, dKernMTGPts, gausWts, iDatom, iDtype, &
          & iqpt, ngdp, gdp, vHarNum )

        if (.false.) then
          do idir = 1, 3
            do iG = 1, ngdp
              write(2012, '(2(i8),2(f15.8))') idir, iG, vEff1IR(iG, idir)
            end do ! iG
          end do
        end if

        if (any(abs(vEff1IR - (-grVeff0IR)) > 1e-7)) NOstopNO'Veff1 and grVeff0 give not the same'

        ! These quantities are not required for Sternheimer routines
        deallocate ( grVeff0IR)!, grRho0MTDummy )

        sumVMTs = cmplx(0.0, 0.0)
        sumVMTs2 = cmplx(0.0, 0.0)
        iatom = 0
        do itype = 1, atoms%ntype
          do ieqat = 1, atoms%neq(itype)
            iatom = iatom + 1
            do idir = 1, 3
               do oqn_l = 0, atoms%lmax(itype)
                 do mqn_m = -oqn_l, oqn_l
                   lm = oqn_l * (oqn_l + 1) + 1 + mqn_m
                   lm2 = oqn_l * (oqn_l + 1) + 1 - mqn_m
                   do imesh = 1, atoms%jri(itype)
                   !todo really grVeff0MT_init
                     sumVMTs(imesh, lm, idir, iatom) = vEff1MT(imesh, lm, iatom, idir) + grVeff0MT(imesh, lm, idir, iatom)
                     sumVMTs2(imesh, lm2, idir, iatom) = (-1)**(mqn_m) * conjg(sumVMTs(imesh, lm, idir, iatom))
                     write(2009, '(5(i8),2(f15.8))') iatom, idir, oqn_l, mqn_m, imesh, sumVMTs(imesh, lm, idir, iatom)
                     write(2010, '(5(i8),2(f15.8))') iatom, idir, oqn_l, mqn_m, imesh, sumVMTs(imesh, lm, idir, iatom)
                   end do
                end do
              end do
            end do
          end do
        end do
        if (any( abs(sumVMTs) >1e-9) .or. any( abs(sumVMTs2) >1e-9)) NOstopNO'not vanishing'



        if (.false.) then
          ! Setting up the effective potential variation like this does not make a difference in the matrix element
          allocate(vEff0IRG_w(ngdp))
          vEff0IRG_w = cmplx(0.0, 0.0)
          vpwStar(:, :) = cmplx(0., 0.)
          call fopen( 1000, name='vpw_eff', status='old', action='read', form='unformatted' )
            read( 1000 ) vpwStar
          call fclose( 1000 )
          !call convertStar2G(V0Fleur%vpw_uw(:, 1), vEff0IRG_w, stars, ngdp, gdp)
          call convertStar2G(vpwStar(:, 1), vEff0IRG_w, stars, ngdp, gdp)
          vEff1IR(:, :) = cmplx(0., 0.)
          do idir = 1, 3
            do iG = 1, ngdp
              Gext(1:3) = matmul( cell%bmat(1:3, 1:3), gdp(1:3, iG) )
              ! Remember to first apply the gradient then warp not vice versa!
              vEff1IR(iG, idir) = -iu * Gext(idir) * vEff0IRG_w(iG)
            end do ! iG
          end do ! idir
        end if

        allocate(w_vEff1IR(ngdp, 3))
        do idir = 1, 3
          call warpIRPot(stars, ngdp, idir, gdp, vEff1IR, w_vEff1IR(:, idir))
        end do

        if (.false.) then
          do idir = 1, 3
            do iG = 1, ngdp
              write(2102, '(2i8,2f15.8)') idir, iG, w_vEff1IR(iG, idir)
            end do ! iG
          end do ! idir
        end if

        call tlmplm4V( atoms, lathar, dimens, enpara, uds, input, tdVx, input%jspins, input%jspins, 1, sumVMTs, rbas1, rbas2, uuilon, &
          &duilon, ulouilopn, ilo2p, nlo_atom )
        call tlmplm4V( atoms, lathar, dimens, enpara, uds, input, tdVy, input%jspins, input%jspins, 2, sumVMTs, rbas1, rbas2, uuilon, &
          &duilon, ulouilopn, ilo2p, nlo_atom )
        call tlmplm4V( atoms, lathar, dimens, enpara, uds, input, tdVz, input%jspins, input%jspins, 3, sumVMTs, rbas1, rbas2, uuilon, &
          &duilon, ulouilopn, ilo2p, nlo_atom )

        call tlmplm4V( atoms, lathar, dimens, enpara, uds, input, tdVx2, input%jspins, input%jspins, 1, sumVMTs2, rbas1, rbas2, uuilon, &
          &duilon, ulouilopn, ilo2p, nlo_atom )
        call tlmplm4V( atoms, lathar, dimens, enpara, uds, input, tdVy2, input%jspins, input%jspins, 2, sumVMTs2, rbas1, rbas2, uuilon, &
          &duilon, ulouilopn, ilo2p, nlo_atom )
        call tlmplm4V( atoms, lathar, dimens, enpara, uds, input, tdVz2, input%jspins, input%jspins, 3, sumVMTs2, rbas1, rbas2, uuilon, &
          &duilon, ulouilopn, ilo2p, nlo_atom )

        ! Integrals for calculating the Sternheimer equation matrix element with the unperturbed Hamiltonian for the displaced atom α
        ! The form which is read out from pottot is used in tlmplm. However, the spherical component of the potential is not used at
        ! all therefore, we can manipulate them as much we want. For radfun and radflo we use vr0 which is not manipulated before
        ! begin used. rbas1 and rbas2 have to be multiplied by r to form the Jacobi determinant for the MT integrals
        call tlmplm4H0( atoms, dimens, enpara, uds, input, tdHS0, 1, logUnit, rbas1, rbas2, uuilon, duilon, ulouilopn, ilo2p, &
                                                                                                           & vEff0MTsh(:, :, :, 1) ) ! add spin

!        allocate(vEff0IRG(ngdp))
!        allocate(nlo_atom(atoms%nat))
!        vEff0IRG = cmplx(0.0, 0.0)
!        call convertStar2G(V0Fleur%vpw(:, 1), vEff0IRG, stars, ngdp, gdp)

!    lmpMax = maxval( (/ (sum( (/ ((2*oqn_l+1)* nRadFun(oqn_l,itype), oqn_l = 0,atoms%lmax(itype)) /) ),itype=1,atoms%ntype) /) )
!        allocate( hFullNoAbcofBK(lmpMax, lmpMax, 3, atoms%nat) )
!        allocate( overlapNoAbcofBK(lmpMax, lmpMax, 3, atoms%ntype) )
!        hFullNoAbcofBK(:, :, :, :) = cmplx(0., 0.)
!        overlapNoAbcofBK(:, :, :, :) = cmplx(0., 0.)
!        call PrepareMTSurfInt( atoms, dimens, kpts, lathar, V0Fleur, iDtype, iDatom, nRadFun, El, rbas1(:, :, :, :, 1), rbas2(:, :, :, :, 1), nmem_atom, mlh_atom,&
!          & clnu_atom, hFullNoAbcofBK, overlapNoAbcofBK, zDummy)
       allocate(surfIntVFast(dimens%neigd, maxval(nobd), 3))

    allocate( veffUvIR(3, (1 * atoms%lmaxd + 1)**2) )
    veffUvIR(:, :) = cmplx(0., 0.)
    call IRcoeffVeffUv( atoms, stars, cell, iDtype, iDatom, ngdp, 1, gdp, veffUvIR, vEff0pwUw )

          allocate (z1nG(dimens%nbasfcn, dimens%neigd, 3)) !segfault in abcof if second argument not like this, this can be optimized
        maxlmp = maxval( (/ (sum( (/ ((2 * oqn_l + 1) * nRadFun(oqn_l,iDtype), oqn_l = 0, atoms%lmax(iDtype)) /) ), iDtype = 1,&
                                                                                                                  &atoms%ntype) /) )
        do ikpt = 1 , kpts%nkpt ! additional k-points still have to be added

          write(*, *) 'ikpt = ', ikpt
          ! The energy differences are calculated for every atom which is redundant. But this is no bottleneck so keeping in memor is worse.
          call GenRecEdifME( ne(ikpt), nobd(ikpt, 1), eig(:, ikpt, 1), eig(:, ikpt, 1), recEdiffME, cutContr )

        !  do ii = 1, ne(1)
        !    do jj = 1, nobd(1, 1)
        !      write(2030, '(f15.8)') recEdiff(ii, jj)
        !      write(2031, '(f15.8)') cutContr(ii, jj)
        !    end do ! jj
        !  end do ! ii

          ! todo does an allocation make sense here?
          z1nG = cmplx(0.0, 0.0)
          surfIntVFast = cmplx(0., 0.)

          call calcSfVeffFast( atoms, kpts, qpts, cell, dimens, ikpt, ikpt, iqpt, kpq2kPrVec, 1, gbas, veffUvIR, nv, ne, nobd, mapGbas, z0, iDtype, iDatom, surfIntVFast)
!          !todo el für l > 3 ist einfach el von l = 3 man kann in radfun und radflo ein array reinstecken, dass um ein l erweitertet ist
!          ! Due to performance reasons we do not make a loop over idir (array of types, type of arrays)
!ORI    GINAL VERSION
          call solveSternheimerEq( atoms, sym, cell, kpts, qpts, dimens, uds, tdHS0, tdVx, stars, gdp, ne, nv, w_vEff1IR, &
            & eig(:, ikpt, ispin), eig(:, mapKpq2K(ikpt, iqpt), ispin), El, nRadFun, iloTable, mapGbas, gbas, zDummy(:, :, ikpt, ispin), &
            & zDummy(:, :, ikpt, ispin), kveclo(:, :), iDtype, iDatom, ikpt, mapKpq2K(ikpt, iqpt), iqpt, 1, &
            & ngdp, nobd, z1nG(:, :, 1), nlo_atom, recEdiffME, kpq2kPrVec, tdVx2, cutContr, surfIntVFast, ngdp, gdp, maxlmp )
          call solveSternheimerEq( atoms, sym, cell, kpts, qpts, dimens, uds, tdHS0, tdVy, stars, gdp, ne, nv, w_vEff1IR, &
            & eig(:, ikpt, ispin), eig(:, mapKpq2K(ikpt, iqpt), ispin), El, nRadFun, iloTable, mapGbas, gbas, zDummy(:, :, ikpt, ispin), &
            & zDummy(:, :, ikpt, ispin), kveclo(:, :), iDtype, iDatom, ikpt, mapKpq2K(ikpt, iqpt), iqpt, 2, &
            & ngdp, nobd, z1nG(:, :, 2), nlo_atom, recEdiffME, kpq2kPrVec, tdVy2, cutContr, surfIntVFast, ngdp, gdp, maxlmp )
          !todo rename last argument in all calls of solveSternheimerEq!
          call solveSternheimerEq( atoms, sym, cell, kpts, qpts, dimens, uds, tdHS0, tdVz, stars, gdp, ne, nv, w_vEff1IR, &
            & eig(:, ikpt, ispin), eig(:, mapKpq2K(ikpt, iqpt), ispin), El, nRadFun, iloTable, mapGbas, gbas, zDummy(:, :, ikpt, ispin), &
            & zDummy(:, :, ikpt, ispin), kveclo(:, :), iDtype, iDatom, ikpt, mapKpq2K(ikpt, iqpt), iqpt, 3, &
            & ngdp, nobd, z1nG(:, :, 3), nlo_atom, recEdiffME, kpq2kPrVec, tdVz2, cutContr, surfIntVFast, ngdp, gdp, maxlmp )

          deallocate(recEdiffME, cutContr)

          call storeZ1nG(atoms, ikpt, 1, mapKpq2K, iDatom, nobd, nv, z1nG)
!          do idir = 1, 3
!            do iBas = 1, nv(1, ikpt)
!              do iband = 1, nobd(ikpt, 1)
!                write(2008, '(4(i8),2(f15.8))') ikpt, idir, iBas, iband, z1nG(iBas, iband, idir)
!                write(2121, '(4(i8),2(es15.8,3x))') ikpt, idir, iBas, iband, z1nG(iBas, iband, idir) - z1Special(iBas, iband, idir, ikpt)
!              end do ! iband
!            end do ! iBas
!            do iBas = nv(1, ikpt) + 1, nv(1, ikpt) + atoms%nlotot
!              do iband = 1, nobd(ikpt, 1)
!                write(2008, '(4(i8),2(f15.8))') ikpt, idir, iBas, iband, z1nG(iBas, iband, idir)
!                write(2121, '(4(i8),2(es15.8,3x))') ikpt, idir, iBas, iband, z1nG(iBas, iband, idir) - z1Special(iBas, iband, idir, ikpt)
!              end do ! iband
!            end do ! iBas
!          end do ! idir

          do idir = 1, 3
            call calcRho1IRValDS( cell, results, nobd(ikpt, 1), nv(1, :), ikpt, iqpt, mapKpq2K(ikpt,iqpt), idir, gbas(:, :),&
              & zDummy(:, :, ikpt, 1), z1nG, rho1IRAfter(:, :, iDatom), gdp2Ind, mapGbas, gdp2iLim, kpq2kPrVec)
          end do ! idir


!        if (ikpt ==56) then
!          call plot1stDensVar(sym, cell, input, lathar, atoms, rho1IR(:, :, iDatom), ngdp, gdp, mlh_atom, nmem_atom, clnu_atom, rho0MT, qpts%bk(:, iqpt), rho1MT(:, :, :, :, iDatom), 1., 1., 1., iDatom, 2, 1) !todo remove counter!
!          do idir = 1, 3
!            do iG = 1, ngdp
!              write(2099, '(4(i8),2(f15.8))') iDatom, ikpt, idir, iG, rho1IR(iG, idir, iDatom)
!            end do ! iG
!          end do ! idir
!        endif
        end do ! ikpt
        deallocate (z1nG)
        do idir = 1, 3
          do iG = 1,ngdp
            ! IR core contribution correction to the density variation
            ! We have a minus in the calculation of rho1IRctC
            rho1IRAfter(iG, idir, iDatom) = rho1IRAfter(iG, idir, iDatom) + rho1IRctC(iG, iDatom, idir)
            write(2123, '(2(i8),2(f15.8))') idir, iG, rho1IRAfter(iG, idir, iDatom)
            write(2124, '(2(i8),1(es15.8))') idir, iG, aimag(rho1IRAfter(iG, idir, iDatom)) - aimag(rho1IR(iG, idir, iDatom))
            if (aimag(rho1IR(iG, idir, iDatom)) >1e-10) then
              write(2125, '(2(i8),1(es15.8))') idir, iG, (aimag(rho1IRAfter(iG, idir, iDatom)) - aimag(rho1IR(iG, idir, iDatom))) / aimag(rho1IR(iG, idir, iDatom))
            end if

          end do
        end do
        call plot1stDensVar(sym, cell, input, lathar, atoms, rho1IRAfter(:, :, iDatom), ngdp, gdp, mlh_atom, nmem_atom, clnu_atom, rho0MT, qpts%bk(:, iqpt), rho1MT(:, :, :, :, iDatom), 1., 1., 1., iDatom, 2, 1) !todo remove counter!

    write(*, *)
    write(*, *)
    write(*, *)
    write(*, *)
  end subroutine testInstq0Conv

  ! Compare the surfinthepsmt and surfinthepsIR and discuss in phd thesis why this is important
  subroutine SurfIntHepsMT( atoms, lathar, kpts, dimens, sym, V0Fleur, cell, input, El, rbas1, rbas2, nRadFun, mlh_atom, nmem_atom, clnu_atom, ilst, nv, gbas, ne, z, ud, kveclo, eig, iloTable, nobd, surfIntMT, onlyPotentialSw )

#include "cppmacro.h"

    use m_abcof
    use m_types, only : t_atoms, t_sphhar, t_potential, t_kpts, t_sym, t_dimension, t_cell, t_input, t_usdus, t_noco
     
    use m_jpConstants, only : c_im, iu, fpi
    use m_gaunt, only : gaunt1

    implicit none

    ! Type parameters
    type(t_atoms),                  intent(in)  :: atoms
    type(t_sphhar),                 intent(in)  :: lathar
    type(t_potential),              intent(in)  :: V0Fleur
    type(t_sym),                    intent(in)  :: sym
    type(t_dimension),              intent(in)  :: dimens
    type(t_cell),                   intent(in)  :: cell
    type(t_kpts),                   intent(in)  :: kpts
    type(t_input),                  intent(in)  :: input
    type(t_usdus),                  intent(in)  :: ud

    logical,                        intent(in)  :: onlyPotentialSw

    ! Array parameters
    real,                           intent(in)  :: eig(:, :, :)
    real,                           intent(in)  :: El(:, 0:, :, :)
    real,                           intent(in)  :: rbas1(:, :, 0:, :)
    real,                           intent(in)  :: rbas2(:, :, 0:, :)
    integer,                        intent(in)  :: nRadFun(0:, :)
    integer,                        intent(in)  :: mlh_atom(:,0:,:)
    integer,                        intent(in)  :: nmem_atom(0:, :)
    complex,                        intent(in)  :: clnu_atom(:,0:,:)
    integer,                        intent(in)  :: gbas(:, :)
    integer,                        intent(in)  :: ilst(:, :, :)
    integer,                        intent(in)  :: kveclo(:,:)
    integer,                        intent(in)  :: nv(:, :)
    integer,                        intent(in)  :: ne(:)
    integer,                        intent(in)  :: nobd(:, :)
    MCOMPLEX,                       intent(in)  :: z(:,:,:,:)
    integer,                        intent(in)  :: iloTable(:, 0:, :)
    complex,           allocatable, intent(out) :: surfIntMT(:, :, :, :)


    ! Type variable
    type(od_inp)                                :: odi
    type(od_sym)                                :: ods
    type(t_noco)                                :: noco

    ! Scalar variable
    integer                                     :: nmat
    integer                                     :: lmp
    integer                                     :: oqn_l
    integer                                     :: mqn_m
    integer                                     :: imesh
    integer                                     :: itype
    integer                                     :: iradf
    integer                                     :: rMt
    integer                                     :: iatom
    integer                                     :: ptsym
    integer                                     :: ilh
    integer                                     :: oqn_l3p
    integer                                     :: oqn_l2p
    integer                                     :: mqn_m4p
    integer                                     :: mqn_m3p
    integer                                     :: mqn_m2p
    real                                        :: gauntVeffNv
    real                                        :: gauntComplete
    integer                                     :: idir
    integer                                     :: oqn_l1p
    integer                                     :: mqn_m1p
    integer                                     :: ieqat
    integer                                     :: imem
    integer                                     :: lm
    integer                                     :: lm1p
    integer                                     :: lmp1p
    integer                                     :: iradf1p
    integer                                     :: lmpMax
    integer                                     :: ikpt
    integer                                     :: ilo
    integer                                     :: ibandK
    integer                                     :: ibandB
    complex                                     :: hsphSymLmp1pLmp
    real                                        :: rbasOverlap
    complex                                     :: veffnSphLmp1pLmp
    real                                        :: epsSym

    ! Array variable
    integer,           allocatable              :: ngoprI(:)
    real,              allocatable              :: hSphUlmRmt(:, :, :)
!    real,              allocatable              :: hepsSphUlmRmtK(:, :, :)
!    real,              allocatable              :: hepsSphUlmRmtB(:, :, :)
    complex,           allocatable              :: vEfflmLm1P(:, :, :)
    complex,           allocatable              :: vEfflmpLmp1P(:, :, :, :)
    complex,           allocatable              :: gauntVec(:, :, :, :)
    complex,           allocatable              :: acof(:, :, :)
    complex,           allocatable              :: bcof(:, :, :)
    complex,           allocatable              :: ccof(:, :, :, :)
    complex,           allocatable              :: mCofK(:, :)
    complex,           allocatable              :: mCofB(:, :)
!    complex,           allocatable              :: hepsMatR2lmp1pK(:, :)
!    complex,           allocatable              :: hepsMatR2lmp1pB(:, :)
    complex,           allocatable              :: zDummy(:, :, :, :)
    complex,           allocatable              :: veffDummy(:, :, :, :)
    complex,           allocatable              :: hFullNoAbcofBK(:, :, :, :)
    complex,           allocatable              :: overlapNoAbcofBK(:, :, :, :)
    complex,           allocatable              :: hFullNoAbcofK(:, :)
    complex,           allocatable              :: hFull(:, :)
    complex,           allocatable              :: overlapNoAbcofK(:, :)
    complex,           allocatable              :: overlap(:, :, :, :)

    complex                                     :: sumTest

    ! Check that abcof
    allocate(acof(dimens%neigd, 0:dimens%lmd, atoms%nat), bcof(dimens%neigd, 0:dimens%lmd, atoms%nat), &
      &ccof(-atoms%llod:atoms%llod, dimens%neigd, atoms%nlod, atoms%nat))
    allocate( noco%alph(atoms%ntype), noco%beta(atoms%ntype) ) !Up to now those variables are only of dummy character
    allocate(ngoprI(atoms%nat))
    ngoprI(:) = 1

    allocate( zDummy(dimens%nbasfcn, dimens%neigd, kpts%nkpt, 1) )
    zDummy = cmplx(0., 0.)
    zDummy = z!cmplx(real(z), 0.)

    ! Test potential
    allocate( veffDummy(atoms%jmtd, (atoms%lmaxd + 1)**2, 3, atoms%nat) )
    vEffDummy(:, :, :, :)  = cmplx(0., 0.)
    vEffDummy(atoms%jri(1), 1, 1:3, 1)  = cmplx(sqrt(fpi), 0.)


    ! Spherical part only contributes until lmax and not lmax + 2 as we have only the ket here expanded until lmax and in the overlap we
    ! cut there at lmax.
    lmpMax = maxval( (/ (sum( (/ ((2*oqn_l+1)* nRadFun(oqn_l,itype), oqn_l = 0,atoms%lmax(itype)) /) ),itype=1,atoms%ntype) /) )
    allocate( hSphUlmRmt(2, lmpMax, atoms%ntype) )
    hSphUlmRmt(:, :, :) = 0.

    allocate( vEfflmLm1P((atoms%lmaxd + 1)**2, (atoms%lmaxd + 1)**2, 3) )
    allocate( vEfflmpLmp1P(lmpMax, lmpMax, 3, atoms%nat) )
    allocate( gauntVec( (atoms%lmaxd + 1)**2, (atoms%lmaxd + 1)**2, 3, atoms%ntype ) )
    allocate( hFullNoAbcofBK(lmpMax, lmpMax, 3, atoms%nat) )
    allocate( overlapNoAbcofBK(lmpMax, lmpMax, 3, atoms%ntype) )

    vEfflmpLmp1P = cmplx(0., 0.)
    gauntVec = cmplx(0., 0.)
    hFullNoAbcofBK = cmplx(0., 0.)
    overlapNoAbcofBK = cmplx(0., 0.)

    do itype = 1, atoms%ntype
      rMt = atoms%jri(itype)
      lmp = 0
      ! Calculate H_sph u_lmp(R_Mt) * R_mT (rbas arrays are already multiplied with the mesh once)
      do oqn_l = 0, atoms%lmax(itype)
        do mqn_m = -oqn_l, oqn_l
          lmp = lmp + 1
          hSphUlmRmt(1, lmp, itype) = El(1, oqn_l, itype, 1) * rbas1(rMt, 1, oqn_l, itype)
          hSphUlmRmt(2, lmp, itype) = El(1, oqn_l, itype, 1) * rbas2(rMt, 1, oqn_l, itype)

          lmp = lmp + 1
          hSphUlmRmt(1, lmp, itype) = rbas1(rMt, 1, oqn_l, itype) + El(1, oqn_l, itype, 1) * rbas1(rMt, 2, oqn_l, itype)
          hSphUlmRmt(2, lmp, itype) = rbas2(rMt, 1, oqn_l, itype) + El(1, oqn_l, itype, 1) * rbas2(rMt, 2, oqn_l, itype)

          do iradf = 3, nRadFun(oqn_l, itype)
            lmp = lmp + 1
            hSphUlmRmt(1, lmp, itype) = El(iradf - 1, oqn_l, itype, 1) * rbas1(rMt, iradf, oqn_l, itype)
            hSphUlmRmt(2, lmp, itype) = El(iradf - 1, oqn_l, itype, 1) * rbas2(rMt, iradf, oqn_l, itype)
          end do ! p
        end do ! mqn_m
      end do ! oqn_l


      iatom = 0
      do ieqat = 1, atoms%neq(itype)
        iatom = iatom + 1
        ptsym = atoms%ntypsy(iatom)
        vEfflmLm1P = cmplx(0., 0.)
        ! Do all the triangular rules math for the effective potential and calculate the product of the the effective potential with
        ! the normal vector and two Gaunt coefficients
        ! We use only the non-spherical part of the potential
        ! NOTE:  For only Veff0 set ilh = 0 to lathar%nlh, and uncomment if that oqn_l2 not smaller than 0, and sum hsph to veff
        if (onlyPotentialSw) then
          do ilh = 0, lathar%nlh(ptsym)
            oqn_l3p = lathar%llh(ilh, ptsym)
            ! As we only treat the non-spherical part oqn_l2p never is smaller than 0
            do oqn_l2p = abs(oqn_l3p - 1), min(oqn_l3p + 1, atoms%lmax(itype)), 2
                        if (oqn_l2p < 0) cycle ! if we want only Veff complete
              do mqn_m4p = -1, 1
                do imem = 1, nmem_atom(ilh, iatom)
                  mqn_m3p = mlh_atom(imem, ilh, iatom)
                  mqn_m2p = mqn_m3p + mqn_m4p ! tested
                  if (abs(mqn_m2p) > oqn_l2p) cycle
                  gauntVeffNv = gaunt1( oqn_l2p, oqn_l3p, 1, mqn_m2p, mqn_m3p, mqn_m4p, atoms%lmaxd )
                  do idir = 1, 3
                    do oqn_l = 0, atoms%lmax(itype)
                      do oqn_l1p = abs(oqn_l2p - oqn_l), min(oqn_l2p + oqn_l, atoms%lmax(itype))
                        do mqn_m = -oqn_l, oqn_l
                          mqn_m1p = mqn_m3p + mqn_m4p + mqn_m
                          if ( abs(mqn_m1p) > oqn_l1p ) cycle
                          lm = oqn_l * (oqn_l + 1) + 1 + mqn_m
                          lm1p = oqn_l1p * (oqn_l1p + 1) + 1 + mqn_m1p
                          gauntComplete = gaunt1( oqn_l1p, oqn_l2p, oqn_l, mqn_m1p, mqn_m2p, mqn_m, atoms%lmaxd )
                      ! Constant potential for testing
                      !    vEfflmLm1P(lm1p, lm, idir) = vEfflmLm1P(lm1p, lm, idir) + c_im(idir, mqn_m4p + 2) * vEffDummy(rMt, ilh + 1, idir, iatom)&
                      !      & * gauntVeffNv * gauntComplete
                          vEfflmLm1P(lm1p, lm, idir) = vEfflmLm1P(lm1p, lm, idir) + c_im(idir, mqn_m4p + 2) * V0Fleur%vr(rMt, ilh, itype, 1)                      &
                            & * clnu_atom(imem, ilh, iatom) * gauntVeffNv * gauntComplete
                        end do ! mqn_m
                      end do ! oqn_l1p
                    end do ! oqn_l
                  end do ! idir
                end do ! imem
              end do ! mqn_m4p
            end do ! oqn_l2p
          end do ! ilh
        else
          do ilh = 1, lathar%nlh(ptsym)
            oqn_l3p = lathar%llh(ilh, ptsym)
            ! As we only treat the non-spherical part oqn_l2p never is smaller than 0
            do oqn_l2p = abs(oqn_l3p - 1), min(oqn_l3p + 1, atoms%lmax(itype)), 2
                   !     if (oqn_l2p < 0) cycle ! if we want only Veff complete
              do mqn_m4p = -1, 1
                do imem = 1, nmem_atom(ilh, iatom)
                  mqn_m3p = mlh_atom(imem, ilh, iatom)
                  mqn_m2p = mqn_m3p + mqn_m4p ! tested
                  if (abs(mqn_m2p) > oqn_l2p) cycle
                  gauntVeffNv = gaunt1( oqn_l2p, oqn_l3p, 1, mqn_m2p, mqn_m3p, mqn_m4p, atoms%lmaxd )
                  do idir = 1, 3
                    do oqn_l = 0, atoms%lmax(itype)
                      do oqn_l1p = abs(oqn_l2p - oqn_l), min(oqn_l2p + oqn_l, atoms%lmax(itype))
                        do mqn_m = -oqn_l, oqn_l
                          mqn_m1p = mqn_m3p + mqn_m4p + mqn_m
                          if ( abs(mqn_m1p) > oqn_l1p ) cycle
                          lm = oqn_l * (oqn_l + 1) + 1 + mqn_m
                          lm1p = oqn_l1p * (oqn_l1p + 1) + 1 + mqn_m1p
                          gauntComplete = gaunt1( oqn_l1p, oqn_l2p, oqn_l, mqn_m1p, mqn_m2p, mqn_m, atoms%lmaxd )
                      ! Constant potential for testing
                      !    vEfflmLm1P(lm1p, lm, idir) = vEfflmLm1P(lm1p, lm, idir) + c_im(idir, mqn_m4p + 2) * vEffDummy(rMt, ilh + 1, idir, iatom)&
                      !      & * gauntVeffNv * gauntComplete
                          vEfflmLm1P(lm1p, lm, idir) = vEfflmLm1P(lm1p, lm, idir) + c_im(idir, mqn_m4p + 2) * V0Fleur%vr(rMt, ilh, itype, 1)                      &
                            & * clnu_atom(imem, ilh, iatom) * gauntVeffNv * gauntComplete
                        end do ! mqn_m
                      end do ! oqn_l1p
                    end do ! oqn_l
                  end do ! idir
                end do ! imem
              end do ! mqn_m4p
            end do ! oqn_l2p
          end do ! ilh
        end if

        ! Sum up the Gaunt coefficients of the spherical Hamiltonian
        do idir = 1, 3
          do oqn_l = 0, atoms%lmax(itype)
            do mqn_m = -oqn_l, oqn_l
              lm = oqn_l * (oqn_l + 1) + 1 + mqn_m
              do oqn_l1p = abs(oqn_l - 1), min(oqn_l + 1, atoms%lmax(itype)), 2
                do mqn_m2p = -1, 1
                  mqn_m1p = mqn_m + mqn_m2p
                  if ( abs(mqn_m1p) > oqn_l1p ) cycle
                  lm1p = oqn_l1p * (oqn_l1p + 1) + 1 + mqn_m1p
                  gauntComplete = gaunt1( oqn_l1p, oqn_l, 1, mqn_m1p, mqn_m, mqn_m2p, atoms%lmaxd )
                  gauntVec(lm1p, lm, idir, itype) = gauntVec(lm1p, lm, idir, itype) + c_im(idir, mqn_m2p + 2) * gauntComplete
                end do ! mqn_m1Pr
              end do ! oqn_l1Pr
            end do ! mqn_m
          end do ! oqn_l
        end do ! idir

        ! Calculate the matrixelement with the potential and the overlap potential matrix element so far that only the matching
        ! coefficients have to be multiplied by a matrix multiplication
        do idir = 1, 3
          lmp = 0
          do oqn_l = 0, atoms%lmax(itype)
            do mqn_m = -oqn_l, oqn_l
              lm = oqn_l * ( oqn_l + 1 ) + 1 + mqn_m
              do iradf = 1, nRadFun(oqn_l, itype)
                lmp = lmp + 1
                lmp1p = 0
                do oqn_l1p = 0, atoms%lmax(itype)
                  do mqn_m1p = -oqn_l1p, oqn_l1p
                    lm1p = oqn_l1p * ( oqn_l1p + 1 ) + 1 + mqn_m1p
                    do iradf1p = 1, nRadFun(oqn_l1p, itype)
                      lmp1p = lmp1p + 1
!                      hsphSymLmp1pLmp = 0.5 * gauntVec(lm1p, lm, idir, itype) * &
                      hsphSymLmp1pLmp = gauntVec(lm1p, lm, idir, itype) * &
                        & ( rbas1(rMT, iradf1p, oqn_l1p, itype) * hSphUlmRmt(1, lmp, itype) &
                        & + rbas2(rMT, iradf1p, oqn_l1p, itype) * hSphUlmRmt(2, lmp, itype) )!&
!                        & + hSphUlmRmt(1, lmp1p, itype)         * rbas1(rMT, iradf, oqn_l, itype) &
!                        & + hSphUlmRmt(2, lmp1p, itype)         * rbas2(rMT, iradf, oqn_l, itype) )
                      rbasOverlap = rbas1(rMt, iradf1p, oqn_l1p, itype) * rbas1(rMt, iradf, oqn_l, itype) &
                                & + rbas2(rMt, iradf1p, oqn_l1p, itype) * rbas2(rMt, iradf, oqn_l, itype)
                      veffnSphLmp1pLmp = rbasOverlap * vEfflmLm1P(lm1p, lm, idir)
                      if (onlyPotentialSw) then
                        hFullNoAbcofBK(lmp1p, lmp, idir, iatom) = veffnSphLmp1pLmp
                      else
                        hFullNoAbcofBK(lmp1p, lmp, idir, iatom) = hsphSymLmp1pLmp + veffnSphLmp1pLmp
                      end if
!vEfflmpLmp1P(lmp, lmp1p, idir, iatom) = rbasOverlap * vEfflmLm1P(lm1p, lm, idir)
                      overlapNoAbcofBK(lmp1p, lmp, idir, itype) = rbasOverlap * gauntVec(lm1p, lm, idir, itype)
                    !write(4002, '(5i8,2f15.8)') idir, oqn_l1p, mqn_m1p, iradf1p, lmp, hFullNoAbcofBK(lmp1p, lmp, idir, iatom)
                    end do ! iradf1p
                  end do ! mqn_m1p
                end do ! oqn_l1p
              end do ! iradf
            end do ! mqn_m
          end do ! oqn_l
        end do ! idir

      end do ! ieqat
    end do ! itype

    allocate(mCofK(lmpMax, maxval(nobd(:, :))) )
    allocate(mCofB(dimens%neigd, lmpMax ))
!    allocate(hepsSphUlmRmtK(2, lmpMax, maxval(nobd(:, :))))
!    allocate(hepsSphUlmRmtB(2, lmpMax, dimens%neigd))
!    allocate(hepsMatR2lmp1pK(lmpMax, maxval(nobd(:, :))))
!    allocate(hepsMatR2lmp1pB(lmpMax, lmpMax)
    allocate( hFullNoAbcofK(dimens%neigd, lmpMax) )
    allocate( hFull(dimens%neigd, maxval(nobd(:, :))) )
    allocate( overlapNoAbcofK( dimens%neigd, lmpMax) )
    allocate( overlap(dimens%neigd, maxval(nobd(:, :)), 3, atoms%nat) )
    allocate(surfIntMT(dimens%neigd, maxval(nobd(:, :)), 3, atoms%nat))
    do ikpt = 1, kpts%nkpt

      surfIntMT = cmplx(0., 0.)
      hFullNoAbcofK = cmplx(0., 0.)
      hFull = cmplx(0., 0.)
      overlapNoAbcofK = cmplx(0., 0.)
      overlap = cmplx(0., 0.)
      acof(:, :, :) = cmplx(0., 0.)
      bcof(:, :, :) = cmplx(0., 0.)
      ccof(:, :, :, :) = cmplx(0., 0.)
      nmat = nv(1, ikpt) + atoms%nlotot
      call abcof ( atoms%lmaxd, atoms%ntype, dimens%neigd, dimens%neigd, atoms%nat, sym%nop, dimens%nvd, dimens%jspd, dimens%lmd, &
        & dimens%nbasfcn, atoms%llod, atoms%nlod, atoms%nlotot, sym%invtab, atoms%ntype, sym%mrot, ngoprI, atoms%taual, atoms%neq, &
        & atoms%lmax, atoms%rmt, cell%omtil, cell%bmat, cell%bbmat, kpts%bk(:, ikpt), &
        & gbas(1, ilst(:nv(input%jspins, ikpt), ikpt, input%jspins)), gbas(2, ilst(:nv(input%jspins, ikpt), ikpt, input%jspins)), &
        & gbas(3, ilst(:nv(input%jspins, ikpt), ikpt, input%jspins)), nv(:, ikpt),  nmat, ne(ikpt), zDummy(:, :, ikpt, 1), ud%us(:, :, 1), &
        & ud%dus(:, :, 1), ud%uds, ud%duds(:, :, 1), ud%ddn(:, :, 1), atoms%invsat, sym%invsatnr, ud%ulos(:, :, 1), &
        & ud%uulon(:, :, 1), ud%dulon(:, :, 1),  ud%dulos(:, :, 1), atoms%llo, atoms%nlo, atoms%l_dulo, &
        & atoms%lapw_l, noco%l_noco, noco%l_ss, 1, noco%alph, noco%beta, noco%qss, kveclo(:, ikpt), odi, ods, &
        & acof, bcof, ccof)

      ! Print wavefunction and complex wavefunction of muffin-tin to compare it with the interstitial wavefunction
      if (.false.) then
        do oqn_l = 0, atoms%lmax(1)
          do mqn_m = -oqn_l, oqn_l
            lm = oqn_l * (oqn_l + 1) + mqn_m
            do ibandK = 1, nobd(ikpt, 1)
              write(2043, '(3(i8),2(f15.8))') ikpt, ibandK, lm + 1, iu**oqn_l * (acof(ibandK, lm, 1) * rbas1(atoms%jmtd, 1, oqn_l, 1) + bcof(ibandK, lm, 1) * rbas1(atoms%jmtd, 2, oqn_l, 1)) / atoms%rmt(1)
            end do ! ibandK
          end do ! mqn_m
        end do ! oqn_l

        do oqn_l = 0, atoms%lmax(1)
          do mqn_m = -oqn_l, oqn_l
            lm = oqn_l * (oqn_l + 1) + mqn_m
            do ibandK = 1, ne(ikpt)!nobd(ikpt, 1)
              write(2044, '(3(i8),2(f15.8))') ikpt, ibandK, lm + 1, iu**oqn_l * conjg((acof(ibandK, lm, 1) * rbas1(atoms%jmtd, 1, oqn_l, 1) + bcof(ibandK, lm, 1) * rbas1(atoms%jmtd, 2, oqn_l, 1)) / atoms%rmt(1))
            end do ! ibandK
          end do ! mqn_m
        end do ! oqn_l
      end if

      ! Reorganize matching coefficients to include the LO contributions more naturally
      iatom = 0
      mCofB(:, :) = cmplx(0., 0.)
      mCofK(:, :) = cmplx(0., 0.)
      do itype = 1, atoms%ntype
        do ieqat = 1, atoms%neq(itype)
          iatom = iatom + 1
          lmp   =   0
          lm    =  -1
          do oqn_l = 0, atoms%lmax(itype)
            do mqn_m = -oqn_l, oqn_l
              lm = lm + 1
              !p = 1
              lmp = lmp + 1
              do ibandB = 1, ne(ikpt)
                mCofB(ibandB, lmp) = acof(ibandB, lm, iatom) * iu**oqn_l
              end do ! ibandB
              !p = 2
              lmp = lmp + 1
              do ibandB = 1, ne(ikpt)
                mCofB(ibandB, lmp) = bcof(ibandB, lm, iatom)* iu**oqn_l
              end do ! ibandB
              !LOs
              do iradf = 3, nRadFun(oqn_l, itype)
                ilo = iloTable(iradf, oqn_l, itype)
                lmp = lmp + 1
                do ibandB = 1, ne(ikpt)
                  mCofB(ibandB, lmp) = ccof(mqn_m, ibandB, ilo, iatom)* iu**oqn_l
                end do ! ibandB
              end do ! iradf
            end do ! mqn_m
          end do ! oqn_l

          do ibandK = 1, nobd(ikpt, 1)
            lmp   =   0
            lm    =  -1
            do oqn_l = 0, atoms%lmax(itype)
              do mqn_m = -oqn_l, oqn_l
                lm = lm + 1
                !p = 1
                lmp = lmp + 1
                mCofK(lmp, ibandK) = acof(ibandK, lm, iatom) * iu**oqn_l

                !p = 2
                lmp = lmp + 1
                mCofK(lmp, ibandK) = bcof(ibandK, lm, iatom)* iu**oqn_l

                !LOs
                do iradf = 3, nRadFun(oqn_l, itype)
                  ilo = iloTable(iradf, oqn_l, itype)
                  lmp = lmp + 1
                  mCofK(lmp, ibandK) = ccof(mqn_m, ibandK, ilo, iatom)* iu**oqn_l
                end do ! iradf
              end do ! mqn_m
            end do ! oqn_l
          end do ! ibandK

      !    hepsSphUlmRmtK(:, :, :) = cmplx(0., 0.)
      !    do ibandK = 1, nobd(ikpt, 1)
      !      lmp = 0
      !      ! Add the eigen energies to the spherical hamiltonian part
      !      do oqn_l = 0, atoms%lmax(itype)
      !        do mqn_m = -oqn_l, oqn_l
      !          lmp = lmp + 1
      !          hepsSphUlmRmtK(1, lmp, ibandK) = hSphUlmRmt(1, lmp, itype) !- eig(ibandK, ikpt, 1) * rbas1(rMt, 1, oqn_l, itype)
      !          hepsSphUlmRmtK(2, lmp, ibandK) = hSphUlmRmt(2, lmp, itype) !- eig(ibandK, ikpt, 1) * rbas2(rMt, 1, oqn_l, itype)

      !          lmp = lmp + 1
      !          hepsSphUlmRmtK(1, lmp, ibandK) = hSphUlmRmt(1, lmp, itype) !- eig(ibandK, ikpt, 1) * rbas1(rMt, 2, oqn_l, itype)
      !          hepsSphUlmRmtK(2, lmp, ibandK) = hSphUlmRmt(2, lmp, itype) !- eig(ibandK, ikpt, 1) * rbas2(rMt, 2, oqn_l, itype)

      !          do iradf = 3, nRadFun(oqn_l, itype)
      !            lmp = lmp + 1
      !            hepsSphUlmRmtK(1, lmp, ibandK) = hSphUlmRmt(1, lmp, itype)! - eig(ibandK, ikpt, 1) * rbas1(rMt, iradf, oqn_l, itype)
      !            hepsSphUlmRmtK(2, lmp, ibandK) = hSphUlmRmt(2, lmp, itype)! - eig(ibandK, ikpt, 1) * rbas2(rMt, iradf, oqn_l, itype)
      !          end do ! iradf
      !        end do ! mqn_m
      !      end do ! oqn_l
      !    end do ! ibandK

      !    hepsSphUlmRmtB(:, :, :) = cmplx(0., 0.)
      !    do ibandB = 1, ne(ikpt)
      !      lmp = 0
      !      ! Add the eigen energies to the spherical hamiltonian part
      !      do oqn_l = 0, atoms%lmax(itype)
      !        do mqn_m = -oqn_l, oqn_l
      !          lmp = lmp + 1
      !          hepsSphUlmRmtB(1, lmp, ibandB) = hSphUlmRmt(1, lmp, itype)! - eig(ibandB, ikpt, 1) * rbas1(rMt, 1, oqn_l, itype)
      !          hepsSphUlmRmtB(2, lmp, ibandB) = hSphUlmRmt(2, lmp, itype)! - eig(ibandB, ikpt, 1) * rbas2(rMt, 1, oqn_l, itype)

      !          lmp = lmp + 1
      !          hepsSphUlmRmtB(1, lmp, ibandB) = hSphUlmRmt(1, lmp, itype)! - eig(ibandB, ikpt, 1) * rbas1(rMt, 2, oqn_l, itype)
      !          hepsSphUlmRmtB(2, lmp, ibandB) = hSphUlmRmt(2, lmp, itype)! - eig(ibandB, ikpt, 1) * rbas2(rMt, 2, oqn_l, itype)

      !          do iradf = 3, nRadFun(oqn_l, itype)
      !            lmp = lmp + 1
      !            hepsSphUlmRmtB(1, lmp, ibandB) = hSphUlmRmt(1, lmp, itype)! - eig(ibandB, ikpt, 1) * rbas1(rMt, iradf, oqn_l, itype)
      !            hepsSphUlmRmtB(2, lmp, ibandB) = hSphUlmRmt(2, lmp, itype)! - eig(ibandB, ikpt, 1) * rbas2(rMt, iradf, oqn_l, itype)
      !          end do ! iradf
      !        end do ! mqn_m
      !      end do ! oqn_l
      !    end do ! ibandB

          ! Print out only the action of the spherical Hamiltonian.
!          if (.false.) then
!            do ibandK = 1, nobd(ikpt, 1)
!              lmp = 0
!              do oqn_l = 0, atoms%lmax(itype)
!                do mqn_m = -oqn_l, oqn_l
!
!                  lmp = lmp + 1
!                  sumTest = hepsSphUlmRmtK(1, lmp, ibandK) * mCofK(lmp, ibandK)
!
!                  lmp = lmp + 1
!                  sumTest = sumTest + hepsSphUlmRmtK(1, lmp, ibandK) * mCofK(lmp, ibandK)
!                  !write(2047, '(4(i8),2(f15.8))') ikpt, ibandK, oqn_l, mqn_m, sumTest
!                  if ( ( abs(real(sumTest)) < 1e-8 ) .and. ( abs(aimag(sumTest)) >= 1e-8 )) then
!                    write(2047, '(4(i8),2(f15.8))') ikpt, ibandK, oqn_l, mqn_m,  0., aimag(sumTest)
!                  else if ( ( abs( real(sumTest)) >= 1e-8 ) .and. ( abs(aimag(sumTest)) < 1e-8 )) then
!                    write(2047, '(4(i8),2(f15.8))') ikpt, ibandK, oqn_l, mqn_m, real(sumTest), 0.
!                  else if ( ( abs(real(sumTest)) < 1e-8 ) .and. ( abs(aimag(sumTest)) < 1e-8 )) then
!                    write(2047, '(4(i8),2(f15.8))') ikpt, ibandK, oqn_l, mqn_m, 0., 0.
!                  else if ( ( abs(real(sumTest)) >= 1e-8 ) .and. ( abs(aimag(sumTest)) >= 1e-8 )) then
!                    write(2047, '(4(i8),2(f15.8))') ikpt, ibandK, oqn_l, mqn_m, real(sumTest), aimag(sumTest)
!                  end if
!
!                end do ! mqn_m
!              end do ! oqn_l
!            end do ! ibandK
!          end if

          do idir = 1, 3
!            hepsMatR2lmp1pK(:, :) = cmplx(0., 0.)
!            do ibandK = 1, nobd(ikpt, 1)
!              lmp1p = 0
!              do oqn_l1p = 0, atoms%lmax(itype)
!                do mqn_m1p = -oqn_l1p, oqn_l1p
!                  lm1p = oqn_l1p * (oqn_l1p + 1) + 1 + mqn_m1p
!                  do iradf1p = 1, nRadFun(oqn_l1p, itype)
!                    lmp1p = lmp1p + 1
!                    lmp = 0
!                    do oqn_l = 0, atoms%lmax(itype)
!                      do mqn_m = -oqn_l, oqn_l
!                        lm = oqn_l * (oqn_l + 1) + 1 + mqn_m
!                        do iradf = 1, nRadFun(oqn_l, itype)
!                          lmp = lmp + 1
!!                          ! Put everything together except for the bra matching coefficients.
!!                          hepsMatR2lmp1pK(lmp1p, ibandK) = hepsMatR2lmp1pK(lmp1p, ibandK) +                                          &
!!                            & ( ( rbas1(rMt, iradf1p, oqn_l1p, itype) * hepsSphUlmRmtK(1, lmp, ibandK) &
!!                            &   + rbas2(rMt, iradf1p, oqn_l1p, itype) * hepsSphUlmRmtK(2, lmp, ibandK) ) &
!!                            & * gauntVec(lm1p, lm, idir, itype) * 0.5 + VefflmpLmp1P(lmp1p, lmp, idir, iatom) ) * mCofK(lmp, ibandK)
!                          ! this one is only the potential, which has to be run from ilh so change this if you uncomment this
!                           hepsMatR2lmp1pK(lmp1p, ibandK) = hepsMatR2lmp1pK(lmp1p, ibandK) + VefflmpLmp1P(lmp, lmp1p, idir, iatom) * mCofK(lmp, ibandK) ! working version without kinetic energy
!                        end do ! iradf
!                      end do ! mqn_m
!                    end do ! oqn_l
!                  end do ! iradf1p
!                end do ! mqn_m1p
!              end do ! oqn_l1p
!            end do ! ibandK

             ! delete surfIntMT for every ikpt
            hFullNoAbcofK(1:ne(ikpt), 1:lmpMax) = &
                                   & matmul( conjg(mCofB(1:ne(ikpt), 1:lmpMax)), hFullNoAbcofBK(1:lmpMax, 1:lmpMax, idir, iatom) )
            hFull(1:ne(ikpt), 1:nobd(ikpt, 1)) = &
                                   & matmul( hFullNoAbcofK(1:ne(ikpt), 1:lmpMax), mCofK(1:lmpMax, 1:nobd(ikpt, 1)) )

            overlapNoAbcofK(1:ne(ikpt), 1:lmpMax) = &
                                   & matmul( conjg(mCofB(1:ne(ikpt), 1:lmpMax)), overlapNoAbcofBK(1:lmpMax, 1:lmpmax, idir, itype) )
            ! has to be iatom because of mCofs(iatom)
            overlap(1:ne(ikpt), 1:nobd(ikpt, 1), idir, iatom) = &
                                   & matmul( overlapNoAbcofK(1:ne(ikpt), 1:lmpMax), mCofK(1:lmpMax, 1:nobd(ikpt, 1)) )

            do ibandK = 1, nobd(ikpt, 1)
              do ibandB = 1, ne(ikpt)
                !epsSym = 0.5 * (eig(ibandK, ikpt, 1) + eig(ibandB, ikpt, 1))
                epsSym = eig(ibandK, ikpt, 1)
                surfIntMT(ibandB, ibandK, idir, iatom) = hFull(ibandB, ibandK) !- epsSym * overlap(ibandB, ibandK, idir, iatom)
              end do ! ibandB
            end do ! ibandK

            ! todo check whether lmpmax for various types is required
            ! Complete the surface integral by the bra matching coefficients.
            ! for performance reasons mCofK instead of mCofB
!           surfIntMT(:ne(ikpt), :nobd(ikpt, 1), idir, iatom) = surfIntMT(:ne(ikpt), :nobd(ikpt, 1), idir, iatom) +                              &
!              & matmul(conjg(mCofB(:ne(ikpt), :lmpMax)), hepsMatR2lmp1pK(:lmpMax, :nobd(ikpt, 1)))

            !mat(:ne(ikpt), :lmpMax) = matmul(conjg(mCofB(:ne(ikpt), :lmpMax)), foo(:lmpMax, :lmpMax) )
            !matFin(:ne(ikpt), :nobd(ikpt, 1), idir) = matmul(mat(:ne(ikpt), :lmpMax), mCofK(:lmpMax, nobd(ikpt, 1)))
!            lmp1p = 0
!            hepsMatR2lmp1pB(:, :) = cmplx(0., 0.)
!            do oqn_l1p = 0, atoms%lmax(itype)
!              do mqn_m1p = -oqn_l1p, oqn_l1p
!                lm1p = oqn_l1p * (oqn_l1p + 1) + 1 + mqn_m1p
!                do iradf1p = 1, nRadFun(oqn_l1p, itype)
!                  lmp1p = lmp1p + 1
!                  lmp = 0
!                  do oqn_l = 0, atoms%lmax(itype)
!                    do mqn_m = -oqn_l, oqn_l
!                      lm = oqn_l * (oqn_l + 1) + 1 + mqn_m
!                      do iradf = 1, nRadFun(oqn_l, itype)
!                        lmp = lmp + 1
!                        ! Put everything together except for the bra matching coefficients.
!                        hepsMatR2lmp1pB() = hepsMatR2lmp1pB + 0.5 * conjg(mCofB(ibandB, lmp1p)) *            &
!                          &   ( hepsSphUlmRmtB(1, lmp1p, ibandB) * rbas1(rMt, iradf, oqn_l, itype)   &
!                          &   + hepsSphUlmRmtB(2, lmp1p, ibandB) * rbas2(rMt, iradf, oqn_l, itype) ) &
!                          & * gauntVec(lm1p, lm, idir, itype)   * mCofK(ibandB, lmp)
!                      end do ! iradf
!                    end do ! mqn_m
!                  end do ! oqn_l
!                end do ! iradf1p
!              end do ! mqn_m1p
!            end do ! oqn_l1p
          end do ! idir
        end do ! ieqat
      end do ! itype

           !NOTE: We want only compare the potential in the surface integrals. Because the kinetic energy shows differences because of the discontinuity, which is in sum zero again which can be seen in the tests of the surface integrals of the dynamical matrix.
      if (onlyPotentialSw) then
        do idir = 1, 3
          do ibandK = 1, nobd(ikpt, 1)
            do ibandB = 1, ne(ikpt)
              if ( ( abs(real(surfIntMT(ibandB, ibandK, idir, 1))) < 1e-8 ) .and. ( abs(aimag(surfIntMT(ibandB, ibandK, idir, 1))) >= 1e-8 )) then
                write(2033, '(4(i8),2(f15.8))') ikpt, idir, ibandK, ibandB, 0., -aimag(surfIntMT(ibandB, ibandK, idir, 1))
              else if ( ( abs( real(surfIntMT(ibandB, ibandK, idir, 1))) >= 1e-8 ) .and. ( abs(aimag(surfIntMT(ibandB, ibandK, idir, 1))) < 1e-8 )) then
                write(2033, '(4(i8),2(f15.8))') ikpt, idir, ibandK, ibandB, -real(surfIntMT(ibandB, ibandK, idir, 1)), 0.
              else if ( ( abs(real(surfIntMT(ibandB, ibandK, idir, 1))) < 1e-8 ) .and. ( abs(aimag(surfIntMT(ibandB, ibandK, idir, 1))) < 1e-8 )) then
                write(2033, '(4(i8),2(f15.8))') ikpt, idir, ibandK, ibandB, 0., 0.
              else if ( ( abs(real(surfIntMT(ibandB, ibandK, idir, 1))) >= 1e-8 ) .and. ( abs(aimag(surfIntMT(ibandB, ibandK, idir, 1))) >= 1e-8 )) then
                write(2033, '(4(i8),2(f15.8))') ikpt, idir, ibandK, ibandB, -real(surfIntMT(ibandB, ibandK, idir, 1)), -aimag(surfIntMT(ibandB, ibandK, idir, 1))
              end if
              write(2111, '(4(i8),1(f15.8))') ikpt, idir, ibandK, ibandB, -aimag(overlap(ibandB, ibandK, idir, iatom))
            end do ! ibandB
          end do ! ibandK
        end do ! idir
      else
        do idir = 1, 3
          do ibandK = 1, nobd(ikpt, 1)
            do ibandB = 1, ne(ikpt)
              if ( ( abs(real(surfIntMT(ibandB, ibandK, idir, 1))) < 1e-8 ) .and. ( abs(aimag(surfIntMT(ibandB, ibandK, idir, 1))) >= 1e-8 )) then
                write(2333, '(4(i8),2(f15.8))') ikpt, idir, ibandK, ibandB, 0., -aimag(surfIntMT(ibandB, ibandK, idir, 1))
              else if ( ( abs( real(surfIntMT(ibandB, ibandK, idir, 1))) >= 1e-8 ) .and. ( abs(aimag(surfIntMT(ibandB, ibandK, idir, 1))) < 1e-8 )) then
                write(2333, '(4(i8),2(f15.8))') ikpt, idir, ibandK, ibandB, -real(surfIntMT(ibandB, ibandK, idir, 1)), 0.
              else if ( ( abs(real(surfIntMT(ibandB, ibandK, idir, 1))) < 1e-8 ) .and. ( abs(aimag(surfIntMT(ibandB, ibandK, idir, 1))) < 1e-8 )) then
                write(2333, '(4(i8),2(f15.8))') ikpt, idir, ibandK, ibandB, 0., 0.
              else if ( ( abs(real(surfIntMT(ibandB, ibandK, idir, 1))) >= 1e-8 ) .and. ( abs(aimag(surfIntMT(ibandB, ibandK, idir, 1))) >= 1e-8 )) then
                write(2333, '(4(i8),2(f15.8))') ikpt, idir, ibandK, ibandB, -real(surfIntMT(ibandB, ibandK, idir, 1)), -aimag(surfIntMT(ibandB, ibandK, idir, 1))
              end if
              write(2311, '(4(i8),1(f15.8))') ikpt, idir, ibandK, ibandB, -aimag(overlap(ibandB, ibandK, idir, iatom))
            end do ! ibandB
          end do ! ibandK
        end do ! idir
      end if
    end do ! ikpt

    ! Note: If we only switch on the potential, this integral compares well to the interstitial version. The same is true for a
    !       constant potential. Comparing directly the wavefunction from the interstitial and the muffin-tin version, the same comes
    !       out. But if we compare the full action of the spherical Hamiltonian, there is a difference in the values which is due to
    !       the fact, that we compare the second derivative of the LAPW basis function which is not continious any more, necessarily
    !       Therefore, we belive this discrepance to come from this fact.

  end subroutine SurfIntHepsMT


  ! Implementation of SurfIntHepsMT without triangular rules
  subroutine SurfIntHepsMTSafe( atoms, lathar, kpts, dimens, sym, V0Fleur, cell, input, El, rbas1, rbas2, nRadFun, mlh_atom, nmem_atom, clnu_atom, ilst, nv, gbas, ne, z, ud, kveclo, eig, iloTable, nobd )

#include "cppmacro.h"

    use m_abcof
    use m_types, only : t_atoms, t_sphhar, t_potential, t_kpts, t_sym, t_dimension, t_cell, t_input, t_usdus, t_noco
     
    use m_jpConstants, only : c_im
    use m_gaunt, only : gaunt1

    implicit none

    ! Type parameters
    type(t_atoms),                  intent(in)  :: atoms
    type(t_sphhar),                 intent(in)  :: lathar
    type(t_potential),              intent(in)  :: V0Fleur
    type(t_sym),                    intent(in)  :: sym
    type(t_dimension),              intent(in)  :: dimens
    type(t_cell),                   intent(in)  :: cell
    type(t_kpts),                   intent(in)  :: kpts
    type(t_input),                  intent(in)  :: input
    type(t_usdus),                  intent(in)  :: ud

    ! Array parameters
    real,                           intent(in)  :: eig(:, :, :)
    real,                           intent(in)  :: El(:, 0:, :, :)
    real,                           intent(in)  :: rbas1(:, :, 0:, :)
    real,                           intent(in)  :: rbas2(:, :, 0:, :)
    integer,                        intent(in)  :: nRadFun(0:, :)
    integer,                        intent(in)  :: mlh_atom(:,0:,:)
    integer,                        intent(in)  :: nmem_atom(0:, :)
    complex,                        intent(in)  :: clnu_atom(:,0:,:)
    integer,                        intent(in)  :: gbas(:, :)
    integer,                        intent(in)  :: ilst(:, :, :)
    integer,                        intent(in)  :: kveclo(:,:)
    integer,                        intent(in)  :: nv(:, :)
    integer,                        intent(in)  :: ne(:)
    integer,                        intent(in)  :: nobd(:, :)
    MCOMPLEX,                       intent(in)  :: z(:,:,:,:)
    integer,                        intent(in)  :: iloTable(:, :, :)


    ! Type variable
    type(od_inp)                                :: odi
    type(od_sym)                                :: ods
    type(t_noco)                                :: noco

    ! Scalar variable
    integer                                     :: nmat
    integer                                     :: lmp
    integer                                     :: oqn_l
    integer                                     :: mqn_m
    integer                                     :: imesh
    integer                                     :: itype
    integer                                     :: iradf
    integer                                     :: rMt
    integer                                     :: iatom
    integer                                     :: ptsym
    integer                                     :: ilh
    integer                                     :: oqn_l3p
    integer                                     :: oqn_l2p
    integer                                     :: mqn_m4p
    integer                                     :: mqn_m3p
    integer                                     :: mqn_m2p
    real                                        :: gauntVeffNv
    real                                        :: gauntComplete
    integer                                     :: idir
    integer                                     :: oqn_l1p
    integer                                     :: mqn_m1p
    integer                                     :: ieqat
    integer                                     :: imem
    integer                                     :: lm
    integer                                     :: lm1p
    integer                                     :: lmp1p
    integer                                     :: iradf1p
    integer                                     :: lmpMax
    integer                                     :: ikpt
    integer                                     :: ilo
    integer                                     :: ibandK
    integer                                     :: ibandB

    ! Array variable
    integer,           allocatable              :: ngoprI(:)
    complex,           allocatable              :: hSphUlmRmt(:, :, :)
    complex,           allocatable              :: vEfflmLm1P(:, :, :)
    complex,           allocatable              :: vEfflmpLmp1P(:, :, :, :)
    complex,           allocatable              :: gauntVec(:, :, :, :)
    complex,           allocatable              :: acof(:, :, :)
    complex,           allocatable              :: bcof(:, :, :)
    complex,           allocatable              :: ccof(:, :, :, :)
    complex,           allocatable              :: mCof(:, :)
    complex,           allocatable              :: hepsSphUlmRmt(:, :, :)
    complex,           allocatable              :: hepsMatR2lmp1p(:, :)
    complex,           allocatable              :: surfIntMT(:, :, :)

    ! Check that abcof
    allocate(acof(dimens%neigd, 0:dimens%lmd, atoms%nat), bcof(dimens%neigd, 0:dimens%lmd, atoms%nat), &
      &ccof(-atoms%llod:atoms%llod, dimens%neigd, atoms%nlod, atoms%nat))
    allocate( noco%alph(atoms%ntype), noco%beta(atoms%ntype) ) !Up to now those variables are only of dummy character
    allocate(ngoprI(atoms%nat))
    ngoprI(:) = 1


    ! Spherical part only contributes until lmax and not lmax + 2 as we have only the ket here expanded until lmax and in the overlap we
    ! cut there at lmax.
    lmpMax = maxval( (/ (sum( (/ ((2*oqn_l+1)* nRadFun(oqn_l,itype), oqn_l = 0,atoms%lmax(itype)) /) ),itype=1,atoms%ntype) /) )
    allocate( hSphUlmRmt(2, lmpMax, atoms%ntype) )
    hSphUlmRmt(:, :, :) = cmplx(0., 0.)

    allocate( vEfflmLm1P((atoms%lmaxd + 1)**2, (atoms%lmaxd + 1)**2, 3) )
    allocate( vEfflmpLmp1P(lmpMax, lmpMax, 3, atoms%nat) )
    allocate( gauntVec( (atoms%lmaxd + 1)**2, (atoms%lmaxd + 1)**2, 3, atoms%ntype ) )

    vEfflmpLmp1P = cmplx(0., 0.)
    gauntVec = cmplx(0., 0.)

    do itype = 1, atoms%ntype
      rMt = atoms%jri(itype)
      lmp = 0
      ! Calculate H_sph u_lmp(R_Mt) * R_mT (rbas arrays are already multiplied with the mesh once)
      do oqn_l = 0, atoms%lmax(itype)
        do mqn_m = -oqn_l, oqn_l
          lmp = lmp + 1
          hSphUlmRmt(1, lmp, itype) = El(1, oqn_l, itype, 1) * rbas1(rMt, 1, oqn_l, itype)
          hSphUlmRmt(2, lmp, itype) = El(1, oqn_l, itype, 1) * rbas2(rMt, 1, oqn_l, itype)

          lmp = lmp + 1
          hSphUlmRmt(1, lmp, itype) = rbas1(rMt, 1, oqn_l, itype) + El(1, oqn_l, itype, 1) * rbas1(rMt, 2, oqn_l, itype)
          hSphUlmRmt(2, lmp, itype) = rbas2(rMt, 1, oqn_l, itype) + El(1, oqn_l, itype, 1) * rbas2(rMt, 2, oqn_l, itype)

          do iradf = 3, nRadFun(oqn_l, itype)
            lmp = lmp + 1
            hSphUlmRmt(1, lmp, itype) = El(iradf - 1, oqn_l, itype, 1) * rbas1(rMt, iradf, oqn_l, itype)
            hSphUlmRmt(2, lmp, itype) = El(iradf - 1, oqn_l, itype, 1) * rbas2(rMt, iradf, oqn_l, itype)
          end do ! p
        end do ! mqn_m
      end do ! oqn_l

      iatom = 0
      do ieqat = 1, atoms%neq(itype)
        iatom = iatom + 1
        ptsym = atoms%ntypsy(iatom)
        vEfflmLm1P = cmplx(0., 0.)
        ! Do all the triangular rules math for the effective potential and calculate the product of the the effective potential with
        ! the normal vector and two Gaunt coefficients
        ! We use only the non-spherical part of the potential
        do ilh = 1, lathar%nlh(ptsym)
          oqn_l3p = lathar%llh(ilh, ptsym)
          do oqn_l2p = 0, atoms%lmax(itype)
            do mqn_m4p = -1, 1
              do imem = 1, nmem_atom(ilh, iatom)
                mqn_m3p = mlh_atom(imem, ilh, iatom)
                do mqn_m2p = -oqn_l2p, oqn_l2p
                  gauntVeffNv = gaunt1( oqn_l2p, oqn_l3p, 1, mqn_m2p, mqn_m3p, mqn_m4p, atoms%lmaxd )
                  do idir = 1, 3
                    do oqn_l = 0, atoms%lmax(itype)
                      do oqn_l1p = 0, atoms%lmax(itype)
                        do mqn_m = -oqn_l, oqn_l
                          do mqn_m1p = -oqn_l1p, oqn_l1p
                            lm = oqn_l * (oqn_l + 1) + 1 + mqn_m
                            lm1p = oqn_l1p * (oqn_l1p + 1) + 1 + mqn_m1p
                            gauntComplete = gaunt1( oqn_l1p, oqn_l2p, oqn_l, mqn_m1p, mqn_m2p, mqn_m, atoms%lmaxd )
                            vEfflmLm1P(lm, lm1p, idir) = vEfflmLm1P(lm, lm1p, idir) + c_im(idir, mqn_m4p + 2) * V0Fleur%vr(rMt, ilh, itype, 1)                      &
                              & * clnu_atom(imem, ilh, iatom) * gauntVeffNv * gauntComplete
                          end do ! mqn_m1p
                        end do ! mqn_m
                      end do ! oqn_l1p
                    end do ! oqn_l
                  end do ! idir
                end do ! mqn_m2p
              end do ! imem
            end do ! mqn_m4p
          end do ! oqn_l2p
        end do ! ilh

        ! Calculate the matrix element with basis function where the spherical part of the potential is within
        do idir = 1, 3
          lmp1p = 0
          do oqn_l1p = 0, atoms%lmax(itype)
            do mqn_m1p = -oqn_l1p, oqn_l1p
              lm1p = oqn_l1p * ( oqn_l1p + 1 ) + 1 + mqn_m1p
              do iradf1p = 1, nRadFun(oqn_l1p, itype)
                lmp1p = lmp1p + 1
                lmp = 0
                do oqn_l = 0, atoms%lmax(itype)
                  do mqn_m = -oqn_l, oqn_l
                    lm = oqn_l * ( oqn_l + 1 ) + 1 + mqn_m
                    do iradf = 1, nRadFun(oqn_l, itype)
                      lmp = lmp + 1
                      VefflmpLmp1P(lmp, lm1p, idir, iatom) = ( rbas1(rMt, iradf1p, oqn_l1p, itype) * rbas1(rMt, iradf, oqn_l, itype)      &
                        & + rbas2(rMt, iradf1p, oqn_l1p, itype) * rbas2(rMt, iradf, oqn_l, itype) ) * vEfflmLm1P(lm, lm1p, idir)
                    end do ! iradf
                  end do ! mqn_m
                end do ! oqn_l
              end do ! iradf1p
            end do ! mqn_m1p
          end do ! oqn_l1p
        end do ! idir

      end do ! ieqat

      ! Sum up the Gaunt coefficients of the spherical Hamiltonian
      do idir = 1, 3
        do oqn_l = 0, atoms%lmax(itype)
          do mqn_m = -oqn_l, oqn_l
            lm = oqn_l * (oqn_l + 1) + 1 + mqn_m
            do oqn_l1p = max(oqn_l - 1, 0), min(oqn_l + 1, atoms%lmax(itype)), 2
              do mqn_m2p = -1, 1
                mqn_m1p = mqn_m + mqn_m2p
                if ( abs(mqn_m1p) > oqn_l1p ) cycle
                lm1p = oqn_l1p * (oqn_l1p + 1) + 1 + mqn_m1p
                gauntComplete = gaunt1( oqn_l1p, oqn_l, 1, mqn_m1p, mqn_m, mqn_m2p, atoms%lmaxd )
                gauntVec(lm1p, lm, idir, itype) = gauntVec(lm1p, lm, idir, itype) + c_im(idir, mqn_m2p + 2) * gauntComplete
              end do ! mqn_m1Pr
            end do ! oqn_l1Pr
          end do ! mqn_m
        end do ! oqn_l
      end do ! idir
    end do ! itype

    allocate(mCof(dimens%neigd, lmpMax))
    allocate(hepsSphUlmRmt(2, lmpMax, maxval(nobd(:, :))))
    allocate(hepsMatR2lmp1p(lmpMax, maxval(nobd(:, :))))
    allocate(surfIntMT(dimens%neigd, maxval(nobd(:, :)), 3))
    do ikpt = 1, kpts%nkpt
      ! We want the k-point resolved surfIntMT
      surfIntMT = cmplx(0., 0.)

      acof(:, :, :) = cmplx(0., 0.)
      bcof(:, :, :) = cmplx(0., 0.)
      ccof(:, :, :, :) = cmplx(0., 0.)
      nmat = nv(1, ikpt) + atoms%nlotot
      call abcof ( atoms%lmaxd, atoms%ntype, dimens%neigd, dimens%neigd, atoms%nat, sym%nop, dimens%nvd, dimens%jspd, dimens%lmd, &
        & dimens%nbasfcn, atoms%llod, atoms%nlod, atoms%nlotot, sym%invtab, atoms%ntype, sym%mrot, ngoprI, atoms%taual, atoms%neq, &
        & atoms%lmax, atoms%rmt, cell%omtil, cell%bmat, cell%bbmat, kpts%bk(:, ikpt), &
        & gbas(1, ilst(:nv(input%jspins, ikpt), ikpt, input%jspins)), gbas(2, ilst(:nv(input%jspins, ikpt), ikpt, input%jspins)), &
        & gbas(3, ilst(:nv(input%jspins, ikpt), ikpt, input%jspins)), nv(:, ikpt),  nmat, ne(ikpt), z(:, :, ikpt, 1), ud%us(:, :, 1), &
        & ud%dus(:, :, 1), ud%uds, ud%duds(:, :, 1), ud%ddn(:, :, 1), atoms%invsat, sym%invsatnr, ud%ulos(:, :, 1), &
        & ud%uulon(:, :, 1), ud%dulon(:, :, 1),  ud%dulos(:, :, 1), atoms%llo, atoms%nlo, atoms%l_dulo, &
        & atoms%lapw_l, noco%l_noco, noco%l_ss, 1, noco%alph, noco%beta, noco%qss, kveclo(:, ikpt), odi, ods, &
        & acof, bcof, ccof)

      iatom = 0
      mCof(:, :) = cmplx(0., 0.)
      do itype = 1, atoms%ntype
        do ieqat = 1, atoms%neq(itype)
          iatom = iatom + 1
          lmp   =   0
          lm    =  -1
          do oqn_l = 0, atoms%lmax(itype)
            do mqn_m = -oqn_l, oqn_l
              lm = lm + 1
              !p = 1
              lmp = lmp + 1
              do ibandK = 1, ne(ikpt)
                mCof(ibandK, lmp) = acof(ibandK, lm, iatom)
              end do ! ibandK
              !p = 2
              lmp = lmp + 1
              do ibandK = 1, ne(ikpt)
                mCof(ibandK, lmp) = bcof(ibandK, lm, iatom)
              end do ! ibandK
              !LOs
              do iradf = 3, nRadFun(oqn_l, itype)
                ilo = iloTable(iradf, oqn_l, itype)
                lmp = lmp + 1
                do ibandK = 1, ne(ikpt)
                  mCof(ibandK, lmp) = ccof(mqn_m, ibandK, ilo, iatom)
                end do ! ibandK
              end do ! iradf
            end do ! mqn_m
          end do ! oqn_l

          hepsSphUlmRmt(:, :, :) = cmplx(0., 0.)
          do ibandK = 1, nobd(ikpt, 1)
            lmp = 0
            ! Add the eigen energies to the spherical hamiltonian part
            do oqn_l = 0, atoms%lmax(itype)
              do mqn_m = -oqn_l, oqn_l
                lmp = lmp + 1
                hepsSphUlmRmt(1, lmp, ibandK) = hSphUlmRmt(1, lmp, itype) - eig(ibandK, ikpt, 1) * rbas1(rMt, 1, oqn_l, itype)
                hepsSphUlmRmt(2, lmp, ibandK) = hSphUlmRmt(2, lmp, itype) - eig(ibandK, ikpt, 1) * rbas2(rMt, 1, oqn_l, itype)

                lmp = lmp + 1
                hepsSphUlmRmt(1, lmp, ibandK) = hSphUlmRmt(1, lmp, itype) - eig(ibandK, ikpt, 1) * rbas1(rMt, 2, oqn_l, itype)
                hepsSphUlmRmt(2, lmp, ibandK) = hSphUlmRmt(2, lmp, itype) - eig(ibandK, ikpt, 1) * rbas2(rMt, 2, oqn_l, itype)

                do iradf = 3, nRadFun(oqn_l, itype)
                  lmp = lmp + 1
                  hepsSphUlmRmt(1, lmp, ibandK) = hSphUlmRmt(1, lmp, itype) - eig(ibandK, ikpt, 1) * rbas1(rMt, iradf, oqn_l, itype)
                  hepsSphUlmRmt(2, lmp, ibandK) = hSphUlmRmt(2, lmp, itype) - eig(ibandK, ikpt, 1) * rbas2(rMt, iradf, oqn_l, itype)
                end do ! iradf
              end do ! mqn_m
            end do ! oqn_l
          end do ! ibandK


          do idir = 1, 3
            hepsMatR2lmp1p(:, :) = cmplx(0., 0.)
            do ibandK = 1, nobd(ikpt, 1)
              lmp1p = 0
              do oqn_l1p = 0, atoms%lmax(itype)
                do mqn_m1p = -oqn_l1p, oqn_l1p
                  lm1p = oqn_l1p * (oqn_l1p + 1) + 1 + mqn_m1p
                  do iradf1p = 1, nRadFun(oqn_l1p, itype)
                    lmp1p = lmp1p + 1
                    lmp = 0
                    do oqn_l = 0, atoms%lmax(itype)
                      do mqn_m = -oqn_l, oqn_l
                        do iradf = 1, nRadFun(oqn_l, itype)
                          lmp = lmp + 1
                          hepsMatR2lmp1p(lmp1p, ibandK) = hepsMatR2lmp1p(lmp1p, ibandK) + ((rbas1(rMt, iradf1p, oqn_l1p, itype)    &
                            &* hepsSphUlmRmt(1, lmp, ibandK) + rbas2(rMt, iradf1p, oqn_l1p, itype) * hepsSphUlmRmt(2, lmp, ibandK))&
                            & + VefflmpLmp1P(lmp, lm1p, idir, iatom)) * gauntVec(lm1p, lm, idir, itype) * mCof(ibandK, lmp)
                        end do ! iradf
                      end do ! mqn_m
                    end do ! oqn_l
                  end do ! iradf1p
                end do ! mqn_m1p
              end do ! oqn_l1p
            end do ! ibandK
            ! todo check whether lmpmax for various types is required
            surfIntMT(:ne(ikpt), :nobd(ikpt, 1), idir) = surfIntMT(:ne(ikpt), :nobd(ikpt, 1), idir) +                              &
              & matmul(conjg(mCof(:ne(ikpt), :lmpMax)), hepsMatR2lmp1p(:lmpMax, :nobd(ikpt, 1)))
          end do ! idir
        end do ! ieqat
      end do ! itype
      do idir = 1, 3
        do ibandK = 1, nobd(ikpt, 1)
          do ibandB = 1, ne(ikpt)
            write(2034, '(4(i8),2(f15.8))') ikpt, idir, ibandK, ibandB, surfIntMT(ibandB, ibandK, idir)
          end do ! ibandB
        end do ! ibandK
      end do ! idir
    end do ! ikpt

  end subroutine SurfIntHepsMTSafe

  subroutine CompareSurfInt( atoms, lathar, kpts, qpts, dimens, sym, V0Fleur, cell, input, ud, results, stars, ngdp, gdp, El, rbas1,&
      & rbas2, nRadFun, mlh_atom, nmem_atom, clnu_atom, ilst, nv, gbas, ne, z, kveclo, eig, iloTable, nobd, kpq2kPrVec, vEff0IRpwUw )

    use m_types
    use m_jpSternhPulaySurface, only : CalcSintKinEps
    use mod_juPhonUtils, only : fopen, fclose
    use m_jpPotDensHelper, only : convertStar2G
    use m_jpSternhPulaySurface, only : calcSfVeffFast, IRcoeffVeffUV

    implicit none

    ! Type parameters
    type(t_atoms),     intent(in) :: atoms
    type(t_sphhar),    intent(in) :: lathar
    type(t_kpts),      intent(in) :: kpts
    type(t_kpts),      intent(in) :: qpts
    type(t_dimension), intent(in) :: dimens
    type(t_sym),       intent(in) :: sym
    type(t_potential), intent(in) :: V0Fleur
    type(t_cell),      intent(in) :: cell
    type(t_input),     intent(in) :: input
    type(t_usdus),     intent(in) :: ud
    type(t_results),   intent(in) :: results
    type(t_stars),     intent(in) :: stars

    ! Scalar parameters
    integer,           intent(in) :: ngdp

    ! Array parameters
    integer,           intent(in) :: gdp(:, :)
    real,              intent(in) :: El(:, 0:, :, :)
    real,              intent(in) :: rbas1(:, :, 0:, :, :)
    real,              intent(in) :: rbas2(:, :, 0:, :, :)
    integer,           intent(in) :: nRadFun(0:, :)
    integer,           intent(in) :: mlh_atom(:,0:,:)
    integer,           intent(in) :: nmem_atom(0:, :)
    complex,           intent(in) :: clnu_atom(:,0:,:)
    integer,           intent(in) :: ilst(:, :, :)
    integer,           intent(in) :: nv(:, :)
    integer,           intent(in) :: gbas(:, :)
    integer,           intent(in) :: ne(:)
    MCOMPLEX,          intent(in) :: z(:,:,:,:)
    integer,           intent(in) :: kveclo(:,:)
    real,              intent(in) :: eig(:, :, :)
    integer,           intent(in) :: iloTable(:, 0:, :)
    integer,           intent(in) :: nobd(:, :)
    integer,           intent(in) :: kpq2kPrVec(:, :, :)
    complex,           intent(in) :: vEff0IRpwUw(:, :)

    ! Scalar variables
    integer                       :: iDatom
    integer                       :: iqpt
    integer                       :: iDtype
    integer                       :: iDeqat
    integer                       :: ikpt
    integer                       :: idir
    integer                       :: ibandK
    integer                       :: ibandB
    integer                       :: coScale
    integer                       :: ikpq
    complex                       :: surfIntRes

    ! Array variables
    !todo read in this array
    complex,           allocatable             :: rho0IRnCtG(:)
    complex,           allocatable             :: surfIntTeps(:, :)
    complex,           allocatable             :: surfIntVFast(:, :, :)
    complex,           allocatable             :: surfIntVeff(:, :)
    complex,           allocatable             :: surfInt(:, :)
    complex,           allocatable             :: veffUvIR(:, :)
    complex,           allocatable             :: zDummy(:, :, :, :)
    complex,           allocatable             :: surfIntMT(:, :, :, :)


    ! This is for the last matrix element the qpw without coretail correction.
    !todo remember that the rho0IRnCTG is warped
    allocate(rho0IRnCTG(ngdp))
    rho0IRnCTG = cmplx(0.0, 0.0)
    call convertStar2G(V0Fleur%vpw(:, 1), rho0IRnCTG, stars, ngdp, gdp)

    ! NOTE: This routine is not the routine used in jpSternhPUlaySurface
    ! todo make the routines consistent
    call SurfIntHepsMT( atoms, lathar, kpts, dimens, sym, V0Fleur, cell, input, El, rbas1(:, :, :, :, 1), rbas2(:, :, :, :, 1),    &
      & nRadFun, mlh_atom, nmem_atom, clnu_atom, ilst, nv, gbas, ne, z, ud, kveclo, eig, iloTable, nobd, surfIntMT, .true. )

    if (.false.) then
      ! This is a safe implementation of surfIntHepsMT which has been served as auxiliary routine for testing surfinthepsmt but is
      ! not up to date any more concerning bugs which have been found later. This routine does not implement the Gaunt triangular
      ! rules. This routine is deprecated and will be outsourced later
      call SurfIntHepsMTSafe( atoms, lathar, kpts, dimens, sym, V0Fleur, cell, input, El, rbas1(:, :, :, :, 1),                    &
        & rbas2(:, :, :, :, 1), nRadFun, mlh_atom, nmem_atom, clnu_atom, ilst, nv, gbas, ne, z, ud, kveclo, eig, iloTable, nobd )
    end if

    ! This is the implementation of the surface integral within the Sternheimer equation with IR wavefunctions. Where all entities
    ! of the integral are summarized and then Rayleigh expansion is performed. Therefore 3 G summs are needed performing as
    ! Gmax * Kmax * Kmax, the faster version is calcsfveffFast
    if (.false.) then
      call calcSfVeffSums( atoms, kpts, cell, stars, nv, ne, nobd, gdp, ngdp, gbas, ilst, z )
    end if

    ! This is only for calcsintkineps, the other routine do this internaly, remember to put this away later
    allocate( zDummy(dimens%nbasfcn, dimens%neigd, kpts%nkpt, 1) )
    zDummy = cmplx(0., 0.)
    zDummy = z !cmplx(real(z), 0.)

    if (.false.) then
      call testIRkinEnerg( atoms, stars, kpts, cell, ngdp, gbas, gdp, ilst, nv, nobd, zDummy)
    end if

    iDtype = 1
    iDatom = 1
    coScale = 1
    !TODO ATTENTION: zDummy is not valid in general!!!
    write(*, *) 'remove zDummy'
    iqpt = 1
    allocate( surfIntTeps(dimens%neigd, maxval(nobd(:, 1))) )
    allocate( surfInt(dimens%neigd, maxval(nobd(:, 1))) )
    !allocate( veffUvIR(3, (coScale * (atoms%lmaxd + 1) + 1)**2) )
    allocate( veffUvIR(3, (coScale * atoms%lmaxd + 1)**2) )
    veffUvIR(:, :) = cmplx(0., 0.)
    call IRcoeffVeffUv( atoms, stars, cell, iDtype, iDatom, ngdp, coScale, gdp, veffUvIR, vEff0IRpwUw )
    allocate( surfIntVFast(dimens%neigd, maxval(nobd), 3) )
    do ikpt = 1, kpts%nkpt
      write(*, *) ikpt
      ikpq = ikpt ! todo change it later
      surfIntVFast(:, :, :) = cmplx(0., 0.)
      call calcSfVeffFast( atoms, kpts, qpts, cell, dimens, ikpt, ikpq, iqpt, kpq2kPrVec, coScale, gbas, veffUvIR, nv, ne, nobd, ilst, z, iDtype, iDatom, surfIntVFast)
      do idir = 1, 3
        surfIntTeps = cmplx(0.0, 0.0)
        surfInt = cmplx(0.0, 0.0)
        call calcSintKinEps(atoms, cell, kpts, qpts, iDtype, iDatom, nobd(ikpt, 1), ne(ikpt), ikpt, ikpq, iqpt, idir, &
          & nv(1, :), gbas, ilst, zDummy(:, :, ikpt, 1), &
          & zDummy(:, :, ikpt, 1), surfIntTeps, surfInt, kpq2kPrVec)
       do ibandK = 1, nobd(ikpt, 1)
         do ibandB = 1, ne(ikpt)
           !surfIntRes = surfIntTeps(ibandB, ibandK) +  surfIntVFast(ibandB, ibandK, idir) 
           !surfIntRes = surfIntTeps(ibandB, ibandK) +  surfIntVFast(ibandB, ibandK, idir) - eig(ibandK, ikpt, 1) * surfInt(ibandB, ibandK)
           !NOTE: We want only compare the potential in the surface integrals. Because the kinetic energy shows differences because of the discontinuity, which is in sum zero again which can be seen in the tests of the surface integrals of the dynamical matrix.
           surfIntRes = surfIntVFast(ibandB, ibandK, idir)
           if ( ( abs(real(surfIntRes)) < 1e-8 ) .and. ( abs(aimag(surfIntRes)) >= 1e-8 )) then
             write(2036, '(4(i8),2(f15.8))') ikpt, idir, ibandK, ibandB, 0., aimag(surfIntRes)
           else if ( ( abs( real(surfIntRes)) >= 1e-8 ) .and. ( abs(aimag(surfIntRes)) < 1e-8 )) then
             write(2036, '(4(i8),2(f15.8))') ikpt, idir, ibandK, ibandB, real(surfIntRes), 0.
           else if ( ( abs(real(surfIntRes)) < 1e-8 ) .and. ( abs(aimag(surfIntRes)) < 1e-8 )) then
             write(2036, '(4(i8),2(f15.8))') ikpt, idir, ibandK, ibandB, 0., 0.
           else if ( ( abs(real(surfIntRes)) >= 1e-8 ) .and. ( abs(aimag(surfIntRes)) >= 1e-8 )) then
             write(2036, '(4(i8),2(f15.8))') ikpt, idir, ibandK, ibandB, real(surfIntRes), aimag(surfIntRes)
           end if
           write(2112, '(4(i8),1(f15.8))') ikpt, idir, ibandK, ibandB, aimag(surfInt(ibandB, ibandK))
!           if ( ( abs(real(surfIntVFast(ibandB, ibandK, idir))) < 1e-8 ) .and. ( abs(aimag(surfIntVFast(ibandB, ibandK, idir))) >= 1e-8 )) then
!             write(2036, '(4(i8),2(f15.8))') ikpt, idir, ibandK, ibandB, 0., aimag(surfIntVFast(ibandB, ibandK, idir))
!           else if ( ( abs( real(surfIntVFast(ibandB, ibandK, idir))) >= 1e-8 ) .and. ( abs(aimag(surfIntVFast(ibandB, ibandK, idir))) < 1e-8 )) then
!             write(2036, '(4(i8),2(f15.8))') ikpt, idir, ibandK, ibandB, real(surfIntVFast(ibandB, ibandK, idir)), 0.
!           else if ( ( abs(real(surfIntVFast(ibandB, ibandK, idir))) < 1e-8 ) .and. ( abs(aimag(surfIntVFast(ibandB, ibandK, idir))) < 1e-8 )) then
!             write(2036, '(4(i8),2(f15.8))') ikpt, idir, ibandK, ibandB, 0., 0.
!           else if ( ( abs(real(surfIntVFast(ibandB, ibandK, idir))) >= 1e-8 ) .and. ( abs(aimag(surfIntVFast(ibandB, ibandK, idir))) >= 1e-8 )) then
!             write(2036, '(4(i8),2(f15.8))') ikpt, idir, ibandK, ibandB, real(surfIntVFast(ibandB, ibandK, idir)), aimag(surfIntVFast(ibandB, ibandK, idir))
!           end if
          write(2114, '(4i8,2f15.8)') ikpt, idir, ibandK, ibandB, surfIntMT(ibandB, ibandK, idir, 1) + surfIntRes
         end do ! ibandB
       end do ! ibandK
     end do ! idir
    end do ! ikpt

    !NOTE: Compare to 2033 and 2111

!    allocate( surfIntTeps(dimens%neigd, maxval(nobd(:, 1))) )
!    allocate( surfInt(dimens%neigd, maxval(nobd(:, 1))) )
!    allocate( surfIntVeff(dimens%neigd, maxval(nobd(:, 1))) )
!    iDatom = 0
!    iqpt = 1
!    do iDtype = 1, atoms%ntype
!      do iDeqat = 1, atoms%neq(iDtype)
!        iDatom = iDatom + 1
!        do ikpt = 1, kpts%nkpt
!          do idir = 1, 3
!            surfIntTeps = cmplx(0.0, 0.0)
!            surfInt = cmplx(0.0, 0.0)
!            surfIntVeff = cmplx(0.0, 0.0)
!!            call calcSintKinEps(atoms, cell, kpts, qpts, results, iDtype, iDatom, nobd(ikpt, 1), ne(ikpt), ikpt, ikpt, iqpt, idir, &
!!              & nv(1, :), eig(:, ikpt, 1), eig(:, ikpt, 1), gbas, ilst, z(:, :, ikpt, 1), &
!!              & z(:, :, ikpt, 1), surfIntTeps, surfInt, kpq2kPrVec)
!            call calcSfVeff( atoms, cell, stars, dimens, qpts, ikpt, ikpt, iqpt, ngdp, ne(ikpt), nobd(ikpt, 1), iDtype, iDatom,    &
!              & idir, nv, gbas(:, ilst(:nv(1, ikpt), ikpt, 1)), gbas(:, ilst(:nv(1, ikpt), ikpt, 1)), gdp, rho0IRncTG,             &
!              & z(:, :, ikpt, 1), z(:, :, ikpt, 1), surfIntVeff, kpq2kPrVec )
!            do ibandK = 1, nobd(ikpt, 1)
!              do ibandB = 1, ne(ikpt)
!                !write(2035, '(4(i8),2(f15.8))') ikpt, idir, ibandK, ibandB, -(surfIntTeps(ibandB, ibandK) + (eig(ibandB, ikpt, 1) - eig(ibandK, ikpt, 1)) * surfInt(ibandB, ibandK)) !+ surfIntVeff(ibandB, ibandK))
!           if ( ( abs(real(surfIntVeff(ibandB, ibandK))) < 1e-8 ) .and. ( abs(aimag(surfIntVeff(ibandB, ibandK))) >= 1e-8 )) then
!             write(2035, '(4(i8),2(f15.8))') ikpt, idir, ibandK, ibandB, 0., aimag(surfIntVeff(ibandB, ibandK))
!           else if ( ( abs( real(surfIntVeff(ibandB, ibandK))) >= 1e-8 ) .and. ( abs(aimag(surfIntVeff(ibandB, ibandK))) < 1e-8 )) then
!             write(2035, '(4(i8),2(f15.8))') ikpt, idir, ibandK, ibandB, real(surfIntVeff(ibandB, ibandK)), 0.
!           else if ( ( abs(real(surfIntVeff(ibandB, ibandK))) < 1e-8 ) .and. ( abs(aimag(surfIntVeff(ibandB, ibandK))) < 1e-8 )) then
!             write(2035, '(4(i8),2(f15.8))') ikpt, idir, ibandK, ibandB, 0., 0.
!           else if ( ( abs(real(surfIntVeff(ibandB, ibandK))) >= 1e-8 ) .and. ( abs(aimag(surfIntVeff(ibandB, ibandK))) >= 1e-8 )) then
!             write(2035, '(4(i8),2(f15.8))') ikpt, idir, ibandK, ibandB, real(surfIntVeff(ibandB, ibandK)), aimag(surfIntVeff(ibandB, ibandK))
!           end if
!              end do ! ibandB
!            end do ! ibandK
!          end do ! idir
!        end do ! ikpt
!      end do ! iDeqat
!    end do ! iDtype

  end subroutine CompareSurfInt

  ! Seperate routine which has been writen to only compare the kinetic energy + the spherical potential in the intersitital, i.e.
  ! the spherical Hamiltonian which can be compared to the MT action.
  subroutine testIRkinEnerg( atoms, stars, kpts, cell, ngdp, gbas, gdp, ilst, nv, nobd, z)

    use mod_juPhonUtils, only : fopen, fclose
    use m_jpPotDensHelper, only : convertStar2G
    use m_jpConstants, only : iu, fpi, tpi
    use m_types
    use m_ylm_old
    use m_sphbes

    implicit none

    ! Type parameter
    type(t_atoms),             intent(in) :: atoms
    type(t_stars),             intent(in) :: stars
    type(t_kpts),              intent(in) :: kpts
    type(t_cell),              intent(in) :: cell

    ! Scalar parameter
    integer,                   intent(in) :: ngdp

    ! Array parameter
    integer,                   intent(in) :: gdp(:, :)
    integer,                   intent(in) :: gbas(:, :)
    integer,                   intent(in) :: ilst(:, :, :)
    integer,                   intent(in) :: nv(:, :)
    integer,                   intent(in) :: nobd(:, :)
    complex,                   intent(in) :: z(:, :, :, :)

    ! Scalar variables
    integer                               :: iDtype
    integer                               :: iDatom
    integer                               :: iG
    integer                               :: oqn_l
    integer                               :: mqn_m
    integer                               :: lm_pre
    integer                               :: lm
    integer                               :: ikpt
    integer                               :: ibandK
    complex                               :: kinEffect
    complex                               :: V0

    ! Array variables
    real                                   :: gpk(3)
    real                                   :: gpkCart(3)
    complex,       allocatable             :: vpw_effSt(:, :)
    complex,       allocatable             :: vpw_effPw(:)
    complex,       allocatable             :: ylm(:)
    complex,       allocatable             :: kinEnergPsi(:, :)
    real,          allocatable             :: sbes(:)

    ! Load unwarped potential from Fleur
    allocate( vpw_effSt( stars%ng3, 1 ) )
    vpw_effSt = 0
    call Fopen( 1000, name='vpw_eff', status='old', form='unformatted' )
    read( 1000 ) vpw_effSt
    call Fclose( 1000 )

    allocate( vpw_effPw( ngdp ) )
    vpw_effPw = 0
    call convertStar2G( vpw_effSt(:, 1), vpw_effPw, stars, ngdp, gdp )

    iDtype = 1
    iDatom = 1
    kinEffect = cmplx(0., 0.)

    allocate( ylm( (atoms%lmax(1) + 1)**2 ) )
    allocate( sbes(0: atoms%lmax(1)) )
    allocate( kinEnergPsi( ( atoms%lmax(1) + 1)**2, maxval(nobd(:, :)) ) )

    ! todo ask Gustav if this is the right way to go
    V0 = 0
    do iG = 1, ngdp
!      if (all(gdp(1:3, iG) == 0)) cycle
      gpkCart(1:3) = matmul( cell%bmat(1:3, 1:3), gdp(1:3, iG) )
      sbes(:) = 0.
      call sphbes(1, norm2(gpkCart(1:3)) * atoms%rmt(iDtype), sbes)

      V0 = V0 + vpw_effPw(iG) * sbes(0)
    end do ! iG

    do ikpt = 1, kpts%nkpt
      kinEnergPsi(:, :) = cmplx(0., 0.)
      do iG = 1, nv(1, ikpt)
        gpk(1:3) = gbas(1:3, ilst(iG, ikpt, 1)) + kpts%bk(1:3, ikpt)
        gpkCart(1:3) = matmul( cell%bmat(1:3, 1:3), gpk(:) )
        kinEffect = norm2(gpkCart(1:3))**2 / 2.

        call ylm4(atoms%lmax(iDtype), gpkCart, ylm)
        call sphbes(atoms%lmax(iDtype), norm2(gpkCart) * atoms%rmt(iDtype), sbes)

        do ibandK = 1, nobd(ikpt, 1)
          do oqn_l = 0, atoms%lmax(iDtype)
            lm_pre = oqn_l * (oqn_l + 1) + 1
            do mqn_m = -oqn_l, oqn_l
              lm = lm_pre + mqn_m
              kinEnergPsi(lm, ibandK) = kinEnergPsi(lm, ibandK) + z(iG, ibandK, ikpt, 1) * (kinEffect + V0) * &
                & exp( iu * tpi * dot_product(gpk(:), atoms%taual(1:3, iDatom)) ) &
                & * sbes(oqn_l) * conjg(ylm(lm)) * fpi * iu**oqn_l / sqrt(cell%omtil) * atoms%rmt(iDtype)
            end do ! mqn_m
          end do ! oqn_l
        end do ! ibandK
      end do ! iG
      do ibandK = 1, nobd(ikpt, 1)
        do oqn_l = 0, atoms%lmax(iDtype)
          lm_pre = oqn_l * (oqn_l + 1) + 1
          do mqn_m = -oqn_l, oqn_l
            lm = lm_pre + mqn_m
           if ( ( abs(real(kinEnergPsi(lm, ibandK))) < 1e-8 ) .and. ( abs(aimag(kinEnergPsi(lm, ibandK))) >= 1e-8 )) then
             write(2046, '(4(i8),2(f15.8))') ikpt, ibandK, oqn_l, mqn_m, 0., aimag(kinEnergPsi(lm, ibandK))
           else if ( ( abs( real(kinEnergPsi(lm, ibandK))) >= 1e-8 ) .and. ( abs(aimag(kinEnergPsi(lm, ibandK))) < 1e-8 )) then
             write(2046, '(4(i8),2(f15.8))')  ikpt, ibandK, oqn_l, mqn_m, real(kinEnergPsi(lm, ibandK)), 0.
           else if ( ( abs(real(kinEnergPsi(lm, ibandK))) < 1e-8 ) .and. ( abs(aimag(kinEnergPsi(lm, ibandK))) < 1e-8 )) then
             write(2046, '(4(i8),2(f15.8))')  ikpt, ibandK, oqn_l, mqn_m, 0., 0.
           else if ( ( abs(real(kinEnergPsi(lm, ibandK))) >= 1e-8 ) .and. ( abs(aimag(kinEnergPsi(lm, ibandK))) >= 1e-8 )) then
             write(2046, '(4(i8),2(f15.8))') ikpt, ibandK, oqn_l, mqn_m, real(kinEnergPsi(lm, ibandK)), aimag(kinEnergPsi(lm, ibandK))
           end if
           ! write(2046, '(4(i8),2(f15.8))') ikpt, ibandK, oqn_l, mqn_m, kinEnergPsi(lm, ibandK)
          end do ! mqn_m
        end do ! oqn_l
      end do ! ibandK
    end do ! ikpt

  end subroutine testIRkinEnerg

  ! IR surface integrals without FFT, where all G-summations are performed
  subroutine calcSfVeffSums( atoms, kpts, cell, stars, nv, ne, nobd, gdp, ngdp, gbas, ilst, z )

    use m_types
    use m_jpConstants, only : iu, fpi, tpi, c_im
    use m_sphbes
    use m_ylm_old
    use mod_juPhonUtils, only : fopen, fclose
    use m_jpPotDensHelper, only : convertStar2G

    implicit none

    ! Type parameters
    type(t_atoms),             intent(in) :: atoms
    type(t_kpts),              intent(in) :: kpts
    type(t_cell),              intent(in) :: cell
    type(t_stars),             intent(in) :: stars

    ! Scalar parameter
    integer,                   intent(in) :: ngdp

    ! Array parameters
    integer,                   intent(in) :: nv(:, :)
    integer,                   intent(in) :: nobd(:, :)
    integer,                   intent(in) :: ne(:)
    integer,                   intent(in) :: gdp(:, :)
    integer,                   intent(in) :: gbas(:, :)
    integer,                   intent(in) :: ilst(:, :, :)
    MCOMPLEX,                  intent(in) :: z(:, :, :, :)

    ! Scalar variables
    integer                               :: ikpt
    integer                               :: iatom
    integer                               :: itype
    complex                               :: pref
    integer                               :: ibasB
    integer                               :: ibasK
    integer                               :: iGdp
    integer                               :: ieqat
    integer                               :: ibandK
    integer                               :: ibandB
    integer                               :: idir
    integer                               :: mqn_m
    integer                               :: ikpq
    complex                               :: phaseFac
    complex                               :: surfIntNonVecBK
    complex                               :: surfIntNonVecB
    complex                               :: surfIntNonVec
    complex                               :: tSummedCitY1t


    ! Array variables
    integer                               :: gSum(3)
    real                                  :: gSumCart(3)
    complex                               :: ylm(4)
    real                                  :: sbes(0:1)
    complex,       allocatable            :: surfInt(:, :, :)
    complex,       allocatable            :: vpw_eff(:, :)
    complex,       allocatable            :: vpwG(:)


    ! Load unwarped potential from Fleur
    allocate( vpw_eff( stars%ng3, 1 ) )
    vpw_eff = 0
    call Fopen( 1000, name='vpw_eff', status='old', form='unformatted' )
    read( 1000 ) vpw_eff
    call Fclose( 1000 )

    ! Transform Fleur potentials into G-vector representation
    allocate( vpwG( ngdp ) )
    vpwG = 0
    call convertStar2G( vpw_eff(:, 1), vpwG, stars, ngdp, gdp )

!    vpwG = 0
!    do iGdp = 1, ngdp
!      if (all(gdp(:, iGdp) == 0)) then
!        vpwG(iGdp) = 1
!      end if
!    end do ! iGdp

    write(*, *) 'calcsfveffsums'
    allocate( surfInt( maxval(ne(:)), maxval(nobd(:, :)), 3))
    do ikpt = 1, kpts%nkpt
      ikpq = ikpt
      iatom = 0
      surfInt(:, :, :) = cmplx(0., 0.)
      do itype = 1, atoms%ntype
        pref = -fpi * iu * atoms%rmt(itype)**2 / cell%omtil
        do ieqat = 1, atoms%neq(itype)
          iatom = iatom + 1
          do ibasB = 1, nv(1, ikpq)
            write(*, *) ibasB, nv(1, ikpq)
            do ibasK = 1, nv(1, ikpt)
              do iGdp = 1, ngdp
                !if (all(gdp(:, iGdp) == 0)) cycle
                !gSum(1:3) = gdp(1:3, iGdp) + Gbas(1:3, ilst(ibasK, ikpt, 1)) - Gbas(1:3, ilst(ibasB, ikpq, 1)) - qpts%bk(1:3, iqpt) - kpq2kPrVec(1:3, ikpt, iqpt)
                gSum(1:3) = gdp(1:3, iGdp) + gbas(1:3, ilst(ibasK, ikpt, 1)) - gbas(1:3, ilst(ibasB, ikpq, 1))
                gSumCart(1:3) = matmul( cell%bmat(1:3, 1:3), gSum(1:3) )

                ylm(:) = cmplx(0., 0.)
                call ylm4( 1, gSumCart, ylm )

                sbes(:) = 0.
                call sphbes(1, norm2(gSumCart) * atoms%rmt(itype), sbes)

                phaseFac = exp( iu * tpi * dot_product(gSum(:), atoms%taual(1:3, iatom)) )


                surfIntNonVecBK = pref * phaseFac * sbes(1) * vpwG(iGdp)
                do ibandK = 1, nobd(ikpt, 1)
                  surfIntNonVecB = surfIntNonVecBK * real(z(ibasK, ibandK, ikpt, 1))
                  do ibandB = 1, ne(ikpq)
                    surfIntNonVec = surfIntNonVecB * real(conjg( z(ibasB, ibandB, ikpq, 1) ))
                    do idir = 1, 3
                      tSummedCitY1t = cmplx(0., 0.)
                      do mqn_m = -1, 1
                        tSummedCitY1t = tSummedCitY1t + conjg(c_im(idir, mqn_m + 2)) * conjg(ylm(mqn_m + 3))
                      end do ! mqn_m
                      surfInt(ibandB, ibandK, idir) = surfInt(ibandB, ibandK, idir) + surfIntNonVec * tSummedCitY1t
                    end do ! idir
                  end do ! ibandK
                end do ! ibandB
              end do ! iGdp
            end do ! ibasB
          end do ! ibasK
        end do ! ieqat
      end do ! itype
      do idir = 1, 3
        do ibandK = 1, nobd(ikpt, 1)
          do ibandB = 1, ne(ikpq)
           if ( ( abs(real(surfInt(ibandB, ibandK, idir))) < 1e-8 ) .and. ( abs(aimag(surfInt(ibandB, ibandK, idir))) >= 1e-8 )) then
             write(2040, '(4(i8),2(f15.8))') ikpt, idir, ibandK, ibandB, 0., aimag(surfInt(ibandB, ibandK, idir))
           else if ( ( abs( real(surfInt(ibandB, ibandK, idir))) >= 1e-8 ) .and. ( abs(aimag(surfInt(ibandB, ibandK, idir))) < 1e-8 )) then
             write(2040, '(4(i8),2(f15.8))') ikpt, idir, ibandK, ibandB, real(surfInt(ibandB, ibandK, idir)), 0.
           else if ( ( abs(real(surfInt(ibandB, ibandK, idir))) < 1e-8 ) .and. ( abs(aimag(surfInt(ibandB, ibandK, idir))) < 1e-8 )) then
             write(2040, '(4(i8),2(f15.8))') ikpt, idir, ibandK, ibandB, 0., 0.
           else if ( ( abs(real(surfInt(ibandB, ibandK, idir))) >= 1e-8 ) .and. ( abs(aimag(surfInt(ibandB, ibandK, idir))) >= 1e-8 )) then
             write(2040, '(4(i8),2(f15.8))') ikpt, idir, ibandK, ibandB, real(surfInt(ibandB, ibandK, idir)), aimag(surfInt(ibandB, ibandK, idir))
           end if
          end do ! ibandK
        end do ! ibandB
      end do ! idir
    end do ! ikpt

  end subroutine calcSfVeffSums

  subroutine CheckSplittedMTSurfIntRoutine( atoms, dimens, kpts, lathar, V0Fleur, sym, cell, input, ud, kveclo, nobd, ne, iloTable, eig, nv, gbas, ilst, nRadFun, El, rbas1, rbas2, nmem_atom, mlh_atom, clnu_atom, z )

    use m_types

    implicit none

    type(t_atoms),                  intent(in)  :: atoms
    type(t_dimension),              intent(in)  :: dimens
    type(t_kpts),                   intent(in)  :: kpts
    type(t_sphhar),                 intent(in)  :: lathar
    type(t_potential),              intent(in)  :: V0Fleur
    type(t_sym),                    intent(in)  :: sym
    type(t_cell),                   intent(in)  :: cell
    type(t_input),                  intent(in)  :: input
    type(t_usdus),                  intent(in)  :: ud


    ! Array parameters
    integer,                        intent(in)  :: kveclo(:,:)
    integer,                        intent(in)  :: nobd(:, :)
    integer,                        intent(in)  :: ne(:)
    integer,                        intent(in)  :: iloTable(:, 0:, :)
    real,                           intent(in)  :: eig(:, :, :)
    integer,                        intent(in)  :: nv(:, :)
    integer,                        intent(in)  :: gbas(:, :)
    integer,                        intent(in)  :: ilst(:, :, :)
    integer,                        intent(in)  :: nRadFun(0:, :)
    real,                           intent(in)  :: El(:, 0:, :, :)
    real,                           intent(in)  :: rbas1(:, :, 0:, :)
    real,                           intent(in)  :: rbas2(:, :, 0:, :)
    integer,                        intent(in)  :: nmem_atom(0:, :)
    integer,                        intent(in)  :: mlh_atom(:,0:,:)
    complex,                        intent(in)  :: clnu_atom(:,0:,:)
    MCOMPLEX,                       intent(in)  :: z(:,:,:,:)


    ! Scalar variables
    integer                                     :: ikpt
    integer                                     :: idir
    integer                                     :: iDtype
    integer                                     :: iDatom
    integer                                     :: lmpMax
    integer                                     :: oqn_l
    integer                                     :: iDeqat
    integer                                     :: iDdir
    integer                                     :: ibandK
    integer                                     :: ibandB

    complex,   allocatable                      :: hFullNoAbcofBK(:, :, :, :)
    complex,   allocatable                      :: overlapNoAbcofBK(:, :, :, :)
    complex,   allocatable                      :: overlap(:, :)
    complex,   allocatable                      :: surfIntMT(:, :)
    complex, allocatable                    :: surfIntDummy(:, :, :, :)

    lmpMax = maxval( (/ (sum( (/ ((2*oqn_l+1)* nRadFun(oqn_l,iDtype), oqn_l = 0,atoms%lmax(iDtype)) /) ),iDtype=1,atoms%ntype) /) )

    allocate( hFullNoAbcofBK(lmpMax, lmpMax, 3, atoms%nat) )
    allocate( overlapNoAbcofBK(lmpMax, lmpMax, 3, atoms%ntype) )

    allocate( overlap(dimens%neigd, maxval(nobd(:, :)) ) )
    allocate( surfIntMT(dimens%neigd, maxval(nobd(:, :))) )


    iDatom = 0
    do iDtype = 1, atoms%ntype
      do iDeqat = 1, atoms%neq(iDtype)
        iDatom = iDatom + 1
        hFullNoAbcofBK(:, :, :, :) = cmplx(0., 0.)
        overlapNoAbcofBK(:, :, :, :) = cmplx(0., 0.)
        call PrepareMTSurfInt( atoms, dimens, kpts, lathar, V0Fleur, iDtype, iDatom, nRadFun, El, rbas1, rbas2, nmem_atom, mlh_atom,&
          & clnu_atom, hFullNoAbcofBK, overlapNoAbcofBK, z)

        do ikpt = 1, kpts%nkpt
          do idir = 1, 3
            overlap(:, :) = cmplx(0., 0.)
            surfIntMT(:, :) = cmplx(0., 0.)
            call CalcSurfIntMT(atoms, dimens, sym, cell, kpts, input, ud, ikpt, ikpt, idir, iDatom, iDtype, nv, gbas, ilst, z(:, :, ikpt, 1), z(:, :, ikpt, 1), nobd, &
              & ne, kveclo, nRadFun, iloTable, hFullNoAbcofBK, overlapNoAbcofBK, overlap, surfIntMT)

            do ibandK = 1, nobd(ikpt, 1)
              do ibandB = 1, ne(ikpt)
                surfIntMT(ibandB, ibandK) = surfIntMT(ibandB, ibandK)! - 0.5 * (eig(ibandK, ikpt, 1) + eig(ibandB, ikpt, 1)) * overlap(ibandB, ibandK) 
                if ( ( abs(real(surfIntMT(ibandB, ibandK))) < 1e-8 ) .and. ( abs(aimag(surfIntMT(ibandB, ibandK))) >= 1e-8 )) then
                  write(2037, '(4(i8),2(f15.8))') ikpt, idir, ibandK, ibandB, 0., -aimag(surfIntMT(ibandB, ibandK))
                else if ( ( abs( real(surfIntMT(ibandB, ibandK))) >= 1e-8 ) .and. ( abs(aimag(surfIntMT(ibandB, ibandK))) < 1e-8 )) then
                  write(2037, '(4(i8),2(f15.8))') ikpt, idir, ibandK, ibandB, -real(surfIntMT(ibandB, ibandK)), 0.
                else if ( ( abs(real(surfIntMT(ibandB, ibandK))) < 1e-8 ) .and. ( abs(aimag(surfIntMT(ibandB, ibandK))) < 1e-8 )) then
                  write(2037, '(4(i8),2(f15.8))') ikpt, idir, ibandK, ibandB, 0., 0.
                else if ( ( abs(real(surfIntMT(ibandB, ibandK))) >= 1e-8 ) .and. ( abs(aimag(surfIntMT(ibandB, ibandK))) >= 1e-8 )) then
                  write(2037, '(4(i8),2(f15.8))') ikpt, idir, ibandK, ibandB, -real(surfIntMT(ibandB, ibandK)), -aimag(surfIntMT(ibandB, ibandK))
                end if
                  write(2113, '(4(i8),2(f15.8))') ikpt, idir, ibandK, ibandB, overlap(ibandB, ibandK)
              end do ! ibandB
            end do ! ibandK
          end do ! idir
        end do ! ikpt
      end do ! iDeqat
    end do ! iDtype

!NOTE: COMPARE 2037 to the output of the surfinthepsmt (should be 2333 and 2311)
    call SurfIntHepsMT( atoms, lathar, kpts, dimens, sym, V0Fleur, cell, input, El, rbas1(:, :, :, :), rbas2(:, :, :, :),    &
      & nRadFun, mlh_atom, nmem_atom, clnu_atom, ilst, nv, gbas, ne, z, ud, kveclo, eig, iloTable, nobd, surfIntDummy, .false. )

  end subroutine CheckSplittedMTSurfIntRoutine

  ! Tests two ways of calculating <Psi|Veff1|Psi>_IR by using the normal potential Veff0
  subroutine TestIRV1ME( atoms, kpts, cell, lathar, input, stars, dimens, V0Fleur, results, memd_atom, logUnit, rho0MT, clnu_atom, nmem_atom, mlh_atom, &
      &gbas, ilst, nv, ne, nobd, z, rho0IR, kpq2kPrVec)

    use m_types
    use m_jpPotDensHelper, only : CalcIRdVxcKern, CalcMTdVxcKern, WarpIRPot, convertStar2G, genPotDensGvecs
    use m_jpSternhHF, only : CalcMEPotIR
    use m_juDFT_NOstopNO, only : juDFT_warn

    implicit none

    ! Type parameter
    type(t_atoms),                   intent(in) :: atoms
    type(t_kpts),                    intent(in) :: kpts
    type(t_cell),                    intent(in) :: cell
    type(t_sphhar),                  intent(in) :: lathar
    type(t_input),                   intent(in) :: input
    type(t_stars),                   intent(in) :: stars
    type(t_dimension),               intent(in) :: dimens
    type(t_results),                 intent(in) :: results
    type(t_potential),               intent(in) :: V0Fleur

    ! Scalar parameter
    integer,                         intent(in) :: memd_atom
    integer,                         intent(in) :: logUnit

    ! Array parameter
    real,                            intent(in) :: rho0MT(:, 0:, :, :)
    complex,                         intent(in) :: clnu_atom(:, 0:, :)
    integer,                         intent(in) :: nmem_atom(0:, :)
    integer,                         intent(in) :: mlh_atom(:, 0:, :)
    integer,                         intent(in) :: gbas(:, :)
    integer,                         intent(in) :: ilst(:, :, :)
    integer,                         intent(in) :: nv(:, :)
    integer,                         intent(in) :: ne(:)
    integer,                         intent(in) :: nobd(:, :)
    MCOMPLEX,                        intent(in) :: z(:, :, :, :) !Attention this is every time complex for inversion symmetry there must be another solution
    complex,                         intent(in) :: rho0IR(:,:) !Interstitial density
    integer,                         intent(in) :: kpq2kPrVec(:, :, :)

    ! Scalar variable
    integer                                     :: nmat
    integer                                     :: ngdp
    integer                                     :: ikpt
    integer                                     :: iqpt
    integer                                     :: idir
    integer                                     :: iband
    integer                                     :: nrBndBra
    integer                                     :: nrBndKet
    complex                                     :: vEff0IRMEInt
    logical                                     :: testGoldstein
    integer                                     :: ngdp2km
    complex                                     :: vEff0IRMEFFTSum
    integer                                     :: iG
    logical                                     :: passed = .true.

    ! Array variable
    integer,           allocatable              :: gdp(:, :)
    integer                                     :: gdp2iLim(2, 3)
    integer,           allocatable              :: gdp2Ind(:, :,:)
    MCOMPLEX,          allocatable              :: vEff1IR(:, :)
    complex,           allocatable              :: grVeff0IR(:, :)
    complex,           allocatable              :: grRho0MT(:, :, :, :)
    complex,           allocatable              :: grVeff0MT(:, :, :, :)
    complex,           allocatable              :: w_vEff1IR(:, :)
    complex,           allocatable              :: vEff0IRMEEFFT(:, :)
    real,              allocatable              :: gWghts(:) ! gaussian weights belonging to gausPts
    complex,           allocatable              :: ylm1(:, : )
    complex,           allocatable              :: ylm2(:, :)
    real,              allocatable              :: dKernMTGPts(:, :, :)
    complex,           allocatable              :: grVxcIRKern(:)
    complex,           allocatable              :: rho0IRPw(:)
    complex,           allocatable              :: vEff0IRG(:)
    complex,           allocatable              :: w_vEff0IRG(:)
    complex,           allocatable              :: z1Dummy(:, :, :)
    complex,           allocatable              :: rho0IRn(:, :)
    ! todo reduce number of variables


    !todo The G-vector set has to be built again because if gdp2Ind contains G-vectors until Gmax the Fleur density is not
    ! reproduced any more, this does only work if it goes until 2kmax, Check this further!!!
    call genPotDensGvecs(stars, cell, input, ngdp, ngdp2km, gdp, gdp2Ind, gdp2iLim, .true.)

    write(logUnit, '(a)') 'Testing different ways to calculate <Psi|Veff|Psi>_IR'
    write(logUnit, '(a)') '-----------------------------------------------------'
    ! This test only works for q = 0
    iqpt = 1

    ! We have to take the unpeturbed effective potential instead of the first-order effective potential, because due to symmetry,
    ! Veff1 is an odd function while the product of two wavefunctions with two odd and even bands are always even!
    ! Using the first-order effective potential delivers zero.


    allocate( vEff0IRMEEFFT(dimens%neigd, maxval(nobd(:, :))) )

    ! this is already warped
    allocate(vEff0IRG(ngdp))
    vEff0IRG = cmplx(0.0, 0.0)
    call convertStar2G(V0Fleur%vpw_uw(:, 1), vEff0IRG, stars, ngdp, gdp)
    allocate(w_vEff0IRG(ngdp))
    w_vEff0IRG = cmplx(0.0, 0.0)
    call convertStar2G(V0Fleur%vpw(:, 1), w_vEff0IRG, stars, ngdp, gdp)

    allocate( rho0IRn(ngdp, 1) )
    allocate(z1Dummy(dimens%nbasfcn, maxval(nobd(:, :)), 3))
    do ikpt = 1, kpts%nkpt
      z1Dummy = cmplx(0.0, 0.0)
      z1Dummy(:, :nobd(ikpt, 1), 1) = z(:, :nobd(ikpt, 1), ikpt, 1)
      nrBndBra = ne(ikpt)
      nrBndKet = nobd(ikpt, 1)
      vEff0IRMEEFFT(:, :) = cmplx(0., 0.)
      idir = 1
      nmat = nv(1, ikpt) + atoms%nlotot
      call calcMEPotIR( stars, dimens, gbas(:, ilst(:nv(input%jspins, ikpt), ikpt, input%jspins)),                    &
        & gbas(:, ilst(:nv(input%jspins, ikpt), ikpt, input%jspins)), nv, vEff0IRG, z(:, :, ikpt, 1), z(:, :, ikpt, 1), gdp,       &
        & nmat, nrBndBra, nrBndKet, ikpt, iqpt, ikpt, ngdp, vEff0IRMEEFFT, kpq2kPrVec ) !todo spin-relevant

      do iband = 1, nrBndKet
        rho0IRn = cmplx(0., 0.)
        call CalcRho1IRValDSn( cell, input, results, nobd(ikpt, 1), nv(1, :), ikpt, iqpt, ikpt, 1, gbas, z(:, :, ikpt, 1), z1Dummy,&
          & rho0IRn, gdp2Ind, ilst, gdp2iLim, kpq2kPrVec, iband )
        vEff0IRMEInt = cmplx(0., 0.)
        ! We have 1 / omtil in the density stemming from the z. We still need a omtil coming from the volume integral in the IR.
        ! Within the FFT the 1 / omtil from the zs and from the integral can cancel away because we never really calculate the
        ! density
        vEff0IRMEInt = vEff0IRMEInt + cell%omtil * dot_product(conjg(w_vEff0IRG(:)), rho0IRn(:, 1))
        if (.false.) then
          write(2050, '(2(i8),2(f15.8))') ikpt, iband, real(vEff0IRMEEFFT(iband, iband))
          write(2051, '(2(i8),2(f15.8))') ikpt, iband, real(vEff0IRMEInt )
        end if
        if (abs(vEff0IRMEEFFT(iband, iband) - vEff0IRMEInt) > 9e-7) passed = .false.
      end do ! iband
    end do ! ikpt

    if (passed) then
      write(logUnit, '(a)') '                                                    |_passed'

    else
      write(logUnit, '(a)') '                                                    |_failed'
      call juDFT_warn( 'Testing different ways to calculate <Psi|Veff|Psi>_IR! failed', &
        &calledby='TestIRV1ME', hint='Debug.' )
    end if

  end subroutine TestIRV1ME

  subroutine calcRho1IRValDSn( cell, input, results, nobd, nv, ikpt, iqpt, ikpq, idir, GbasVec, z, z1nG, rho1IR, gdp2Ind, mapGbas,   &
      & gdp2iLimi, kpq2kPrVec, iband )

#include "cppmacro.h"

    use m_types

    implicit none

    ! Type parameter
    type(t_cell),    intent(in)    :: cell
    type(t_input),   intent(in)    :: input
    type(t_results), intent(in)    :: results

    ! Scalar parameter
    integer,         intent(in)    :: nobd
    integer,         intent(in)    :: ikpt
    integer,         intent(in)    :: iqpt
    integer,         intent(in)    :: ikpq
    integer,         intent(in)    :: idir
    integer,         intent(in)    :: iband

    ! Array parameter
    integer,         intent(in)    :: nv(:)
    integer,         intent(in)    :: GbasVec(:, :)
    MCOMPLEX,        intent(in)    :: z(:, :)
    complex,         intent(in)    :: z1nG(:, :, :)
    integer,         intent(in)    :: gdp2iLimi(2, 3)
    integer,         intent(in)    :: gdp2Ind(gdp2iLimi(1, 1):gdp2iLimi(2, 1), gdp2iLimi(1, 2):gdp2iLimi(2, 2),                     &
                                             &gdp2iLimi(1, 3):gdp2iLimi(2, 3))
    integer,         intent(in)    :: mapGbas(:, :, :) !todo rename ilst to this everywhere
    integer,         intent(in)    :: kpq2kPrVec(:, :, :)
    complex,         intent(inout) :: rho1IR(:, :)

    ! Scalar variables
    integer                        :: iGppp
    integer                        :: iGp
    integer                        :: Gind
    real                           :: prf

    ! Array variables
    integer                        :: Gvec(3)

    ! This routine is a fork of the original rho1IRValDS only creating the unperturbed density for every band. Since this is for
    ! the first Sternheimer matrix element on the right side, we need no spin-degeneracy factor 2, no 2 from product rule, no
    ! point dependent weights-occupation number products w_iks
    ! TODO omtil still under question, because the z's are already normalized
    prf = 1. / cell%omtil
    do iGppp = 1, nv(ikpq)
      do iGp = 1, nv(ikpt)
        Gvec(:) = GbasVec(:, mapGbas(iGppp, ikpq, 1)) - GbasVec(:, mapGbas(iGp, ikpt, 1)) +  kpq2kPrVec(:, ikpt, iqpt) ! G = G" + G' - G_f
        Gind = gdp2Ind(Gvec(1), Gvec(2), Gvec(3))
        if ( Gind /= 0 ) then
          rho1IR(Gind, idir) = rho1IR(Gind, idir) + prf * conjg( z(iGp, iband) ) * z1nG(iGppp, iband, idir)
        end if
      end do
    end do

  end subroutine calcRho1IRValDSn

  ! Prints analytical z1 = -i (k + G) z0 into fort.2063 for k = 0 in band band representation
  subroutine TestSternhBandBandAnalyt( atoms, dimens, cell, kpts, nv, ne, gbas, nobd, mapGbas, eig, z0, kveclo )

    use m_jpConstants, only : iu
    use m_types
    use mod_juPhonUtils, only : inv

    implicit none

    ! Type parameters
    type(t_atoms),                  intent(in) :: atoms
    type(t_dimension),              intent(in) :: dimens
    type(t_cell),                   intent(in) :: cell
    type(t_kpts),                   intent(in) :: kpts

    ! Array parameters
    integer,                        intent(in) :: nv(:, :)
    integer,                        intent(in) :: ne(:)
    integer,                        intent(in) :: gbas(:, :)
    integer,                        intent(in) :: nobd(:, :)
    integer,                        intent(in) :: mapGbas(:, :, :)
    real,                           intent(in) :: eig(:,:,:)
    MCOMPLEX,                       intent(in) :: z0(:, :, :, :)
    integer,                        intent(in) :: kveclo(:, :)

    ! Scalar variables
    integer                                    :: isp
    integer                                    :: ikpt
    integer                                    :: idir
    integer                                    :: iBas
    integer                                    :: ibandK
    integer                                    :: ibandB
    real                                       :: invsAtNr

    ! Array variables
    complex,           allocatable             :: z1Special(:, :, :)
    complex,           allocatable             :: z1SpecialTest(:, :, :)
    complex,           allocatable             :: zDummy(:, :, :, :)
    complex,           allocatable             :: z0k0Inv(:, :)
    complex,           allocatable             :: idMTest(:, :)
    complex,           allocatable             :: z1SpecialBand(:, :, :)
    complex,           allocatable             :: z1BandResult(:, :, :)
    real                                       :: kExt(3)
    real                                       :: Gext(3)

    allocate( z1Special(dimens%nbasfcn, dimens%neigd, 3) )
    allocate( z1SpecialTest(dimens%nbasfcn, dimens%neigd, 3) )

    write(*, *) 'remove zDummy from TestSternhBandBandAnalyt'
    allocate( zDummy(dimens%nbasfcn, dimens%neigd, kpts%nkpt, 1) )
    zDummy = cmplx(0., 0.)
    !zDummy = cmplx(real(z0), 0.)
    zDummy = z0

    ikpt = 1
    write(*, *) 'foo', ikpt
    invsAtNr = 1 / atoms%nat

    allocate( z1BandResult(ne(ikpt), nobd(ikpt, 1), 3) )
    z1BandResult = cmplx(0., 0.)

    ! Calculate the analytical z1 = -i(k + G) z0 for q = 0 and k = 0
    z1Special(:, :, :) = cmplx(0.0, 0.0)
    kExt(1:3) = matmul( cell%bmat(1:3, 1:3), kpts%bk(1:3, ikpt) )
    do idir = 1, 3
      do iBas = 1, nv(1, ikpt)
        Gext(1:3) = matmul( cell%bmat(1:3, 1:3), gbas(1:3, mapGbas(iBas, ikpt, 1)) )
        do ibandK = 1, nobd(ikpt, 1)
          z1Special(iBas, ibandK, idir) = z1Special(iBas, ibandK, idir) - iu * cmplx(invsAtNr, 0) * ( Gext(idir) + kExt(idir) )&
            & * zDummy(iBas, ibandK, ikpt, 1)
        end do ! ibandK
      end do ! iBas
      do iBas = nv(1, ikpt) + 1, nv(1, ikpt) + atoms%nlotot
        Gext(1:3) = matmul( cell%bmat(1:3, 1:3), gbas(1:3, mapGbas(kveclo(iBas - nv(1, ikpt), ikpt), ikpt, 1)) )
        do ibandK = 1, nobd(ikpt, 1)
          z1Special(iBas, ibandK, idir) = z1Special(iBas, ibandK, idir) - iu * cmplx(invsAtNr, 0) * ( Gext(idir) + kExt(idir) )&
            & * zDummy(iBas, ibandK, ikpt, 1)
        end do ! ibandK
      end do ! iBas

    end do ! idir
    if (.false.) then
      do idir = 1, 3
        do iBas = 1, nv(1, ikpt)
          do ibandB = 1, nobd(ikpt, 1)
            write(2061, '(4(i8),2(f15.8))') ikpt, idir, iBas, ibandB, z1Special(iBas, ibandB, idir)
          end do ! iband
        end do ! iBas
      end do ! idir
    end if

    rewind(2115)
    do idir = 1, 3
      do ibandK = 1, nobd(ikpt, 1)
        do ibandB = 1, ne(ikpt)
          read(2115) z1BandResult(ibandB, ibandK, idir)
        end do ! ibandB
      end do ! ibandK
    end do ! idir


!    ! Invert the z0 for k = 0
     ! NOTE: This only works for systems with inversion symmetry
     allocate( z0k0Inv(nv(1, ikpt), ne(ikpt)) )
     allocate( idMTest(nv(1, ikpt), ne(ikpt)) )
     z0K0Inv(:, :) = 0.
     z0k0Inv(:nv(1, ikpt), :ne(ikpt)) = inv(z0(:nv(1, ikpt), :ne(ikpt), ikpt, 1))


     if (.false. ) then
       ! Inverse is really inverse, this was tested!
       idMTest(:nv(1, ikpt), :ne(ikpt)) = matmul(z0k0Inv(:nv(1, ikpt), :ne(ikpt)), z0(:nv(1, ikpt), :ne(ikpt), ikpt, 1) )
       do ibandB = 1, ne(ikpt)
         do iBas = 1, nv(1, ikpt)
           write(2062, '(2(i8), 2(f15.8))') iBas, ibandB, idMTest(iBas, ibandB)
           if (ibandB == iBas) then
             if ( abs(idMTest(iBas, ibandB) - 1.) > 1e-8 ) NOstopNO'matrix inversion failure'
           else
             if ( abs(idMTest(iBas, ibandB))  > 1e-8 ) NOstopNO'matrix inversion failure'
           end if
         end do ! iBas
       end do ! ibandB
     end if

    ! Transform z1
    allocate( z1SpecialBand(ne(ikpt), nobd(ikpt, 1), 3) )
    z1SpecialBand = cmplx(0., 0.)
    do idir = 1, 3
      ! We have to keep this in this dimension order to work!
      z1SpecialBand(:ne(ikpt), :nobd(ikpt, 1), idir) = matmul( z0k0Inv(:ne(ikpt), :nv(1, ikpt)), z1Special(:nv(1, ikpt), :nobd(ikpt, 1), idir) )
    end do ! idir

    if (.false.) then
    ! Transformation from basis-function to wave-function space
      z1SpecialTest(:, :, :) = cmplx(0., 0.)
      do idir = 1, 3
        z1SpecialTest(:nv(1, ikpt), :nobd(ikpt, 1), idir) = matmul(z0(:nv(1, ikpt), :ne(ikpt), ikpt, 1), z1SpecialBand(:ne(ikpt), :nobd(ikpt, 1), idir))
      end do ! idir
      do idir = 1, 3
        do iBas = 1, nv(1, ikpt)
          do ibandB = 1, nobd(ikpt, 1)
            ! Have to be equal to fort.2061
            write(2090, '(4(i8),2(f15.8))') ikpt, idir, iBas, ibandB, z1SpecialTest(iBas, ibandB, idir)
            if (abs(z1SpecialTest(iBas, ibandB, idir) - z1Special(iBas, ibandB, idir)) > 1e-9) NOstopNO'matrix and trafo failure'
          end do ! iband
        end do ! iBas
      end do ! idir
    end if

    ! Print z1 in band band representation
    do idir = 1, 3
      do ibandK = 1, nobd(ikpt, 1)
        do ibandB = 1, ne(ikpt)
          ! Does only hold for k = 0
          if (abs(eig(ibandB, ikpt, 1) - eig(ibandK, ikpt, 1)) < 1e-10) then
            write(2063, '(3(i8),2(f15.8))') idir, ibandK, ibandB, cmplx(0.,0.)
            write(2117, '(3(i8),2(es15.8))') idir, ibandK, ibandB, cmplx(0.,0.) - z1BandResult(ibandB, ibandK, 1)
          else
            write(2063, '(3(i8),2(f15.8))') idir, ibandK, ibandB, (eig(ibandB, ikpt, 1) - eig(ibandK, ikpt, 1)) * z1SpecialBand(ibandB, ibandK, idir)
            write(2117, '(3(i8),2(es15.8))') idir, ibandK, ibandB, (eig(ibandB, ikpt, 1) - eig(ibandK, ikpt, 1)) * z1SpecialBand(ibandB, ibandK, idir) - z1BandResult(ibandB, ibandK, idir)
          end if
        end do ! ibandB
      end do ! iBandK
    end do ! idir

  end subroutine TestSternhBandBandAnalyt

  ! Tests whether the first order eigen energy vanishes for q = 0 and a monotatomic system
  subroutine TestKSEnerg1stVar( atoms, cell, lathar, enpara, stars, dimens, sym, input, usdus, V0Fleur, kpts, qpts, results,       &
      & memd_atom, ngdp, logUnit, gdp, rho0IR, rho0MT, mlh_atom, nmem_atom, clnu_atom, rbas1, rbas2, uuilon, duilon, ulouilopn, ilo2p,&
      & ne, El, nv, nobd, GbasVec, ilst, kveclo, iloTable, z, kpq2kPrVec, eig, nRadFun, vEff0MTsh, vEff0IRpwUw )

    use m_juDFT_NOstopNO, only : juDFT_warn
    use m_types

    implicit none

    type(t_atoms),                  intent(in)  :: atoms
    type(t_cell),                   intent(in)  :: cell
    type(t_sphhar),                 intent(in)  :: lathar
    type(t_enpara),                 intent(in)  :: enpara
    type(t_stars),                  intent(in)  :: stars
    type(t_dimension),              intent(in)  :: dimens
    type(t_sym),                    intent(in)  :: sym
    type(t_input),                  intent(in)  :: input
    type(t_usdus),                  intent(in)  :: usdus
    type(t_potential),              intent(in)  :: V0Fleur
    type(t_kpts),                   intent(in)  :: kpts
    type(t_kpts),                   intent(in)  :: qpts
    type(t_results),                intent(in) :: results

    ! Scalar parameters
    integer,                        intent(in)  :: memd_atom
    integer,                        intent(in)  :: ngdp
    integer,                        intent(in)  :: logUnit

    ! Array parameters
    integer,                        intent(in)  :: gdp(:, :)
    complex,                        intent(in)  :: rho0IR(:, :) !n3d !todo  with spin
    real,                           intent(in)  :: rho0MT(:, 0:, :, :) ! todo with spin
    integer,                        intent(in)  :: mlh_atom(:,0:,:)
    integer,                        intent(in)  :: nmem_atom(0:, :)
    complex,                        intent(in)  :: clnu_atom(:,0:,:)
    real,                           intent(in)  :: rbas1(:, :, 0:, :, :)
    real,                           intent(in)  :: rbas2(:, :, 0:, :, :)
    real,                           intent(in)  :: uuilon(:, :)
    real,                           intent(in)  :: duilon(:, :)
    real,                           intent(in)  :: ulouilopn(:, :, :)
    integer,                        intent(in)  :: ilo2p(:, :)
    integer,                        intent(in)  :: ne(:)
    real,                           intent(in)  :: El(:, 0:, :, :)
    integer,                        intent(in)  :: nv(:, :)
    integer,                        intent(in)  :: nobd(:, :)
    integer,                        intent(in)  :: GbasVec(:, :)
    integer,                        intent(in)  :: ilst(:, :, :)
    integer,                        intent(in)  :: kveclo(:,:)
    integer,                        intent(in)  :: iloTable(:, 0:, :)
    MCOMPLEX,                       intent(in)  :: z(:, :, :, :) ! Attention: see zBar
    integer,                        intent(in)  :: kpq2kPrVec(:, :, :)
    real,                           intent(in)  :: eig(:, :, :)
    integer,                        intent(in)  :: nRadFun(0:, :)
    complex,                        intent(in)  :: vEff0IRpwUw(:, :)
    complex,                        intent(in)  :: vEff0MTsh(:, :, :, :)

    ! Scalar variables
    integer                                     :: ikpt
    integer                                     :: idir
    integer                                     :: iband
    integer                                     :: iDtype
    integer                                     :: iDatom
    logical                                     :: passed = .true.

    ! Array variables
    complex,           allocatable              :: eig1(:, :, :)

    write(logUnit, '(a)') 'Test vanishing of epsilon1 for q = 0 for monoatomic system.'
    write(logUnit, '(a)') '-----------------------------------------------------------'
    allocate ( eig1( maxval(nobd(:, :)), kpts%nkpt, 3) )
    eig1 = cmplx(0., 0.)

    iDtype = 1
    iDatom = 1

    call CalcKSEnerg1stVar( atoms, cell, lathar, enpara, stars, dimens, sym, input, usdus, V0Fleur, kpts, qpts, results,       &
      & memd_atom, ngdp, logUnit, gdp, rho0IR, rho0MT, mlh_atom, nmem_atom, clnu_atom, rbas1, rbas2, uuilon, duilon, ulouilopn, ilo2p,&
      & ne, El, nv, nobd, GbasVec, ilst, kveclo, iloTable, z, kpq2kPrVec, eig, nRadFun, iDtype, iDatom, eig1, vEff0MTsh, vEff0IRpwUw )

    if (.false.) then
      do ikpt = 1, kpts%nkpt
        do iband = 1, nobd(ikpt, 1)
          do idir = 1, 3
            write(2100, '(3i8, 2f15.8)') ikpt, idir, iband, eig1(iband, ikpt, idir)
          end do ! iband
        end do ! ikpt
      end do ! idir
    end if

    if ( any(abs(eig1(:, :, :)) > 9e-6) ) passed = .false.

    if (passed) then
      write(logUnit, '(a)') '                                                          |_passed'
    else
      write(logUnit, '(a)') '                                                          |_failed'
      call juDFT_warn( 'Test of vanishing eps1 for q = 0 for monoatomic system failed', &
        &calledby='TestKSEnerg1stVar', hint='Debug.' )
    end if

  end subroutine TestKSEnerg1stVar

  subroutine CalcKSEnerg1stVar( atoms, cell, lathar, enpara, stars, dimens, sym, input, usdus, V0Fleur, kpts, qpts, results,       &
      & memd_atom, ngdp, logUnit, gdp, rho0IR, rho0MT, mlh_atom, nmem_atom, clnu_atom, rbas1, rbas2, uuilon, duilon, ulouilopn, ilo2p,&
      & ne, El, nv, nobd, GbasVec, ilst, kveclo, iloTable, z, kpq2kPrVec, eig, nRadFun, iDtype, iDatom, eig1, vEff0MTsh, vEff0IRpwUw )

    use m_jpGrVeff0, only : GenGrVeff0
     
    use m_types, only : t_atoms, t_cell, t_sphhar, t_enpara, t_stars, t_dimension, t_sym, t_input, t_usdus, t_potential, t_kpts, &
                                                                                                       & t_results, t_noco, t_tlmplm
    use m_abcof
    use m_jpSternhHF, only : calcMEPotIR
    use m_jpSternhPulaySurface, only : tlmplm4H0, calcHS0MT, calcSfVeffFast, calcSintKinEps, IRcoeffVeffUV
    use m_jpConstants, only : iu
    use m_jpPotDensHelper, only : calcIRdVxcKern, calcMTdVxcKern, WarpIRPot, ConvertStar2G
    use mod_juPhonUtils, only : CalcGrR2FinLh

    implicit none

    ! Type parameters
    type(t_atoms),                  intent(in)  :: atoms
    type(t_cell),                   intent(in)  :: cell
    type(t_sphhar),                 intent(in)  :: lathar
    type(t_enpara),                 intent(in)  :: enpara
    type(t_stars),                  intent(in)  :: stars
    type(t_dimension),              intent(in)  :: dimens
    type(t_sym),                    intent(in)  :: sym
    type(t_input),                  intent(in)  :: input
    type(t_usdus),                  intent(in)  :: usdus
    type(t_potential),              intent(in)  :: V0Fleur
    type(t_kpts),                   intent(in)  :: kpts
    type(t_kpts),                   intent(in)  :: qpts
    type(t_results),                intent(in) :: results

    ! Scalar parameters
    integer,                        intent(in)  :: memd_atom
    integer,                        intent(in)  :: ngdp
    integer,                        intent(in)  :: iDtype
    integer,                        intent(in)  :: iDatom
    integer,                        intent(in)  :: logUnit

    ! Array parameters
    integer,                        intent(in)  :: gdp(:, :)
    complex,                        intent(in)  :: rho0IR(:, :) !n3d !todo  with spin
    real,                           intent(in)  :: rho0MT(:, 0:, :, :) ! todo with spin
    integer,                        intent(in)  :: mlh_atom(:,0:,:)
    integer,                        intent(in)  :: nmem_atom(0:, :)
    complex,                        intent(in)  :: clnu_atom(:,0:,:)
    real,                           intent(in)  :: rbas1(:, :, 0:, :, :)
    real,                           intent(in)  :: rbas2(:, :, 0:, :, :)
    real,                           intent(in)  :: uuilon(:, :)
    real,                           intent(in)  :: duilon(:, :)
    real,                           intent(in)  :: ulouilopn(:, :, :)
    integer,                        intent(in)  :: ilo2p(:, :)
    integer,                        intent(in)  :: ne(:)
    real,                           intent(in)  :: El(:, 0:, :, :)
    integer,                        intent(in)  :: nv(:, :)
    integer,                        intent(in)  :: nobd(:, :)
    integer,                        intent(in)  :: GbasVec(:, :)
    integer,                        intent(in)  :: ilst(:, :, :)
    integer,                        intent(in)  :: kveclo(:,:)
    integer,                        intent(in)  :: iloTable(:, 0:, :)
    MCOMPLEX,                       intent(in)  :: z(:, :, :, :) ! Attention: see zBar
    integer,                        intent(in)  :: kpq2kPrVec(:, :, :)
    real,                           intent(in)  :: eig(:, :, :)
    integer,                        intent(in)  :: nRadFun(0:, :)
    complex,                        intent(in)  :: vEff0MTsh(:, :, :, :)
    complex,                        intent(in)  :: vEff0IRpwUw(:, :)
    complex,                        intent(out) :: eig1(:, :, :)

    ! Type variables !todo beware maybe not take them from fleur_init might be dangerous
    type(t_noco)                               :: noco
    type(od_inp)                               :: odi
    type(od_sym)                               :: ods
    type(t_tlmplm)                             :: td4HS0

    ! Scalar variables
    logical                                     :: harSw =.true.
    logical                                     :: extSw =.true.
    logical                                     :: xcSw =.true.
    logical                                     :: testGoldstein
    integer                                     :: itype
    integer                                     :: ilh
    integer                                     :: imesh
    integer                                     :: idir
    integer                                     :: iatom
    integer                                     :: lm
    integer                                     :: lmpMax
    integer                                     :: nmat
    integer                                     :: ieqat
    integer                                     :: oqn_l
    integer                                     :: mqn_m
    integer                                     :: lm_pre
    integer                                     :: ikpt
    integer                                     :: nrBnd
    integer                                     :: maxNrBnd
    integer                                     :: ieig
    integer                                     :: iBas
    integer                                     :: lmp
    integer                                     :: iband
    integer                                     :: ilo
    integer                                     :: coScale
    integer                                     :: p
    integer                                     :: ibandB
    integer                                     :: ibandK
    integer                                     :: ptsym
    integer                                     :: imem

    integer                                     :: iG
    logical                                     :: vExtFull

    ! Array variables
    complex,           allocatable              :: grVeff0IR(:, :)
    complex,           allocatable              :: grVeff0MT(:, :, :, :)
    complex,           allocatable              :: grRho0IR(:, :) ! Dummy quantity at the moment
    complex,           allocatable              :: grRho0MT(:, :, :, :) ! Dummy quantity at the moment
    real,              allocatable              :: gausWts(:) ! gaussian weights belonging to gausPts
    complex,           allocatable              :: ylm(:, : )
    real,              allocatable              :: dKernMTGPts(:, :, :)
    complex,           allocatable              :: grVxcIRKern(:)
    real,              allocatable              :: r2Rho0MT(:, :, :, :)
    complex,           allocatable              :: r2GrRho0MT(:, :, :, :)
    complex,           allocatable              :: rho0MTsh(:, :, :, :)
    complex,           allocatable              :: acof(:, :, :)
    complex,           allocatable              :: bcof(:, :, :)
    complex,           allocatable              :: ccof(:, :, :, :)
    complex,           allocatable              :: acofBar(:, :, :)
    complex,           allocatable              :: bcofBar(:, :, :)
    complex,           allocatable              :: ccofBar(:, :, :, :)
    integer,           allocatable              :: ngoprI(:)
    complex,           allocatable              :: zBar(:, :)
    complex,           allocatable              :: mCoef(:, :, :)
    complex,           allocatable              :: mCoefBar(:, :, :)
    integer,           allocatable              :: nlo_atom(:)
    MCOMPLEX,          allocatable              :: vEff1IR(:, :)
    MCOMPLEX,          allocatable              :: h0MTBv(:, :)
    MCOMPLEX,          allocatable              :: s0MTBv(:, :)
    MCOMPLEX,          allocatable              :: h0MTKv(:, :)
    MCOMPLEX,          allocatable              :: s0MTKv(:, :)
    complex,           allocatable              :: w_grVeff0IR(:, :)
    real                                        :: Gext(3)
    complex,           allocatable              :: veffUvIR(:, :)
    complex,           allocatable              :: surfIntTeps(:, :)
    complex,           allocatable              :: surfInt(:, :)
    complex,           allocatable              :: rho0IRpw(:, :)
    complex,           allocatable             :: surfIntVFast(:, :, :)


    ! Generate gradient of unperturbed density
    ! The factor r^2 has beeen divided out so that rho0MT is only rho0MT and not r^2 rho0MT as it is done in Fleur. The factor
    ! sqrt(4pi) for the zeroth component was already dividied out when constructing rho0MT in cdnmt routine in Fleur. Here to
    ! improve stability of the gradient routine we derive r2Rho0MT and divide out the r^2 again later. Doing so avoids the
    ! subtraction of small numbers close to the core.
    allocate( r2Rho0MT( atoms%jmtd, 0:lathar%nlhd, atoms%ntype, 1) )
    allocate( grRho0MT( atoms%jmtd, ( atoms%lmaxd + 1 )**2, atoms%nat, 3) )
    r2Rho0MT(:, :, :, :) = 0.
    grRho0MT(:, :, :, :) = cmplx(0., 0.)

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

    do itype = 1, atoms%ntype
      do ilh = 0, lathar%nlhd
        do imesh = 1, atoms%jri(itype)
          r2Rho0MT(imesh, ilh, itype, 1) = rho0MT(imesh, ilh, itype, 1) * atoms%rmsh(imesh, itype) * atoms%rmsh(imesh, itype)
        end do
      end do
    end do

    call calcGrR2FinLH( atoms, lathar, clnu_atom, nmem_atom, mlh_atom, r2Rho0MT(:, :, :, 1), r2GrRho0MT )

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
              end do ! imesh
            end do ! mqn_m
          end do ! oqn_l
        end do ! ieqat
      end do ! itype
    end do ! idir

    ! Generate gradient of unperturbed effective potential, switches have been set in declaration.
    !todo the switches have to placed still
    call calcIRdVxcKern(stars, gdp, ngdp, rho0IR(:, 1), grVxcIRKern) ! add spin for this and next line
    call calcMTdVxcKern(atoms, dimens, lathar, rho0MT(:, :, :, 1), nmem_atom, clnu_atom, mlh_atom, gausWts, ylm, dKernMTGPts)
    !TODO REVIEW WHICH TERMS ARE USED FOR FIRST ITERATION OF STERNHEIMMER
    vExtFull = .false.
    call GenGrVeff0( atoms, cell, dimens, stars, ngdp, harSw, extSw, xcSw, gdp, rho0IRpw, rho0MTsh, grRho0IR, grRho0MT, gausWts, &
      & ylm, dKernMTGPts, grVxcIRKern, testGoldstein, vExtFull, grVeff0IR, grVeff0MT ) ! add spin

!    allocate(w_grVeff0IR(ngdp, 3))
!    w_grVeff0IR = cmplx(0., 0.)
!    do idir = 1, 3
!      call warpIRPot(stars, ngdp, idir, gdp, grVeff0IR, w_grVeff0IR(:, idir))
!    end do ! idir

!    do idir = 1, 3
!      do iG = 1, ngdp
!        write(2101, '(2i8,2f15.8)') idir, iG, -w_grVeff0IR(iG, idir)
!      end do ! iG
!    end do ! idir

    ! Call of abcof for standard kets of MT potentials matrix element
    maxNrBnd = maxval(nobd(:, :))
    allocate(noco%alph(atoms%ntype), noco%beta(atoms%ntype)) ! Dummy variables to run abcof
    allocate(acof(maxNrBnd, 0:dimens%lmd, atoms%nat), bcof(maxNrBnd, 0:dimens%lmd, atoms%nat), &
      &ccof(-atoms%llod:atoms%llod, maxNrBnd, atoms%nlod, atoms%nat))

    allocate(acofBar(maxNrBnd, 0:dimens%lmd, atoms%nat), bcofBar(maxNrBnd, 0:dimens%lmd, atoms%nat), &
      &ccofBar(-atoms%llod:atoms%llod, maxNrBnd, atoms%nlod, atoms%nat))
    allocate(ngoprI(atoms%nat))
    ngoprI(:) = 1
    allocate( vEff1IR(dimens%neigd, maxNrBnd) )

    coScale = 1.
    allocate( veffUvIR(3, (coScale * atoms%lmaxd + 1)**2) )
    allocate( surfIntVFast(dimens%neigd, maxval(nobd), 3) )

    allocate( surfIntTeps(dimens%neigd, maxval(nobd(:, 1))) )
    allocate( surfInt(dimens%neigd, maxval(nobd(:, 1))) )
    allocate(zBar(dimens%nbasfcn, maxNrBnd))

    allocate( h0MTBv(dimens%neigd, maxNrBnd), s0MTBv(dimens%neigd, maxNrBnd) )
    allocate( h0MTKv(dimens%neigd, maxNrBnd), s0MTKv(dimens%neigd, maxNrBnd) )

    lmpMax = maxval( (/ (sum( (/ ((2*oqn_l+1)* nRadFun(oqn_l,itype), oqn_l = 0,atoms%lmax(itype)) /) ),itype=1,atoms%ntype) /) )
    ! We dimension the array until dimens%neigd, however only fill them up the occupied bands for reasons of performance
    allocate( mCoef(dimens%neigd, lmpMax, atoms%nat), mCoefBar(dimens%neigd, lmpMax, atoms%nat) )
    mCoef = cmplx(0., 0.)
    mCoefBar = cmplx(0., 0.)

    ! Integrals for calculating the Sternheimer equation matrix element with the unperturbed Hamiltonian for the displaced atom α
    ! The form which is read out from pottot is used in tlmplm. However, the spherical component of the potential is not used at
    ! all therefore, we can manipulate them as much we want. For radfun and radflo we use vr0 which is not manipulated before
    ! begin used. rbas1 and rbas2 have to be multiplied by r to form the Jacobi determinant for the MT integrals
    call tlmplm4H0( atoms, dimens, enpara, usdus, input, td4HS0, 1, logUnit, rbas1, rbas2, uuilon, duilon, ulouilopn, ilo2p, &
                                                                                                           & Veff0MTsh(:, :, :, 1) )

    allocate(nlo_atom(atoms%nat))
    nlo_atom = 0

    do ikpt = 1, kpts%nkpt

      ! Generate abcofs
      nrBnd = nobd(ikpt, 1) ! todo spin-relevant
      nmat = nv(1, ikpt) + atoms%nlotot
      acof(:, :, :) = cmplx(0.0, 0.0)
      bcof(:, :, :) = cmplx(0.0, 0.0)
      ccof(:, :, :, :) = cmplx(0.0, 0.0)
      call abcof ( atoms%lmaxd, atoms%ntype, dimens%neigd, maxNrBnd, atoms%nat, sym%nop, dimens%nvd, dimens%jspd, dimens%lmd, &
        & dimens%nbasfcn, atoms%llod, atoms%nlod, atoms%nlotot, sym%invtab, atoms%ntype, sym%mrot, ngoprI, atoms%taual, atoms%neq, &
        & atoms%lmax, atoms%rmt, cell%omtil, cell%bmat, cell%bbmat, kpts%bk(:, ikpt), &
        & GbasVec(1, ilst(:nv(input%jspins, ikpt), ikpt, input%jspins)), GbasVec(2, ilst(:nv(input%jspins, ikpt), ikpt, input%jspins)), &
        & GbasVec(3, ilst(:nv(input%jspins, ikpt), ikpt, input%jspins)), nv(:, ikpt),  nmat, nrBnd, z(:, :, ikpt, 1), usdus%us(:, :, 1), &
        & usdus%dus(:, :, 1), usdus%uds, usdus%duds(:, :, 1), usdus%ddn(:, :, 1), atoms%invsat, sym%invsatnr, usdus%ulos(:, :, 1), &
        & usdus%uulon(:, :, 1), usdus%dulon(:, :, 1),  usdus%dulos(:, :, 1), atoms%llo, atoms%nlo, atoms%l_dulo, &
        & atoms%lapw_l, noco%l_noco, noco%l_ss, 1, noco%alph, noco%beta, noco%qss, kveclo(:, ikpt), odi, ods, &
        & acof, bcof, ccof)

        ! Calculate surface integral
        veffUvIR(:, :) = cmplx(0., 0.)
        call IRcoeffVeffUv( atoms, stars, cell, iDtype, iDatom, ngdp, coScale, gdp, veffUvIR, vEff0IRpwUw )
        surfIntVFast(:, :, :) = cmplx(0., 0.)
        call calcSfVeffFast( atoms, kpts, qpts, cell, dimens, ikpt, ikpt, ikpt, kpq2kPrVec, coScale, GbasVec, veffUvIR, nv, ne, nobd, ilst,                     &
          & z(:, :, :, :), iDtype, iDatom, surfIntVFast)
      do idir = 1, 3
        ! Calculating zPrime for the bar acof bcof and ccof of the MT Hamiltonian matrix element
        zBar = 0
        do ieig = 1, nrBnd
          do iBas = 1, nv(1, ikpt)
           ! Gext(:) = matmul( cell%bmat, GbasVec(:, ilst(iBas, ikpt, 1)) )
           ! todo To increase performance, we can source out the matrix multiplication of Gext as we calculate a lot redundantly.
            Gext(:) = matmul( cell%bmat, GbasVec(:, ilst(iBas, ikpt, 1)) + kpts%bk(:, ikpt))
            zBar(iBas, ieig) = iu * Gext(idir) * z(iBas, ieig, ikpt, 1)
          end do
          do iBas = nv(1, ikpt) + 1, nv(1, ikpt) + atoms%nlotot
           ! Gext(:) = matmul( cell%bmat, GbasVec(:, ilst(kveclo(iBas - nv(1, ikpt), ikpt), ikpt, 1)) )
            Gext(:) = matmul( cell%bmat, GbasVec(:, ilst(kveclo(iBas - nv(1, ikpt), ikpt), ikpt, 1)) + kpts%bk(:, ikpt))
            zBar(iBas, ieig) = iu * Gext(idir) * z(iBas, ieig, ikpt, 1)
          end do
        end do

        ! Calculating acofBar, bcofBar, ccofBar using zBar
        nmat = nv(1, ikpt) + atoms%nlotot
        call abcof ( atoms%lmaxd, atoms%ntype, maxNrbnd, maxNrBnd, atoms%nat, sym%nop, dimens%nvd, dimens%jspd, dimens%lmd, &
          & dimens%nbasfcn, atoms%llod, atoms%nlod, atoms%nlotot, sym%invtab, atoms%ntype, sym%mrot, ngoprI, atoms%taual, atoms%neq, &
          & atoms%lmax, atoms%rmt, cell%omtil, cell%bmat, cell%bbmat, kpts%bk(:, ikpt), &
          & GbasVec(1, ilst(:nv(input%jspins, ikpt), ikpt, input%jspins)), GbasVec(2, ilst(:nv(input%jspins, ikpt), ikpt, input%jspins)), &
          & GbasVec(3, ilst(:nv(input%jspins, ikpt), ikpt, input%jspins)), nv(:, ikpt),  nmat, nrBnd, zBar(:, :), usdus%us(:, :, 1), &
          & usdus%dus(:, :, 1), usdus%uds, usdus%duds(:, :, 1), usdus%ddn(:, :, 1), atoms%invsat, sym%invsatnr, usdus%ulos(:, :, 1), &
          & usdus%uulon(:, :, 1), usdus%dulon(:, :, 1),  usdus%dulos(:, :, 1), atoms%llo, atoms%nlo, atoms%l_dulo, &
          & atoms%lapw_l, noco%l_noco, noco%l_ss, 1, noco%alph, noco%beta, noco%qss, kveclo(:, ikpt), odi, ods, &
          & acofBar, bcofBar, ccofBar)

        ! Rearrange potential acofs, bcofs, ccofs ensuring an consistent handling of LO terms within the loop the structure compared to LAPWs
        iatom  = 0
        mCoef =  0
        mCoefBar = 0
        do itype = 1, atoms%ntype
          do ieqat = 1, atoms%neq(itype)
            iatom = iatom + 1
            lmp   =   0
            lm    =  -1
            do oqn_l = 0, atoms%lmax(itype)
              do mqn_m = -oqn_l, oqn_l
                lm = lm + 1
                !p = 1
                lmp = lmp + 1
                do iband = 1, nrBnd
                  mCoef(iband, lmp, iatom) = acof(iband, lm, iatom)
                  mCoefBar(iband, lmp, iatom) = acofBar(iband, lm, iatom)
                end do ! iband
                !p = 2
                lmp = lmp + 1
                do iband = 1, nrBnd
                  mCoef(iband, lmp, iatom) = bcof(iband, lm, iatom)
                  mCoefBar(iband, lmp, iatom) = bcofBar(iband, lm, iatom)
                end do ! iband
                !LOs
                do p = 3, nRadFun(oqn_l, itype)
                  ilo = iloTable(p, oqn_l, itype)
                  lmp = lmp + 1
                  do iband = 1, nrBnd
                    mCoef(iband, lmp, iatom) = ccof(mqn_m, iband, ilo, iatom)
                    mCoefBar(iband, lmp, iatom) = ccofBar(mqn_m, iband, ilo, iatom)
                  end do ! iband
                end do
              end do
            end do
          end do
        end do

!        lmp = 0
!        do oqn_l = 0, atoms%lmax(1)
!          do mqn_m = -oqn_l, oqn_l
!            do p = 1, nRadFun(oqn_l, 1)
!              lmp = lmp + 1
!              do iband = 1, nrBnd
!                write(2106, '(4i8,2f15.8)') idir, ikpt, lmp, iband, mCoef(iband, lmp, 1)
!                write(2107, '(4i8,2f15.8)') idir, ikpt, lmp, iband, mCoefBar(iband, lmp, 1)
!              end do ! iband
!            end do ! p
!          end do ! mqn_m
!        end do ! oqn_l

        ! Calculate the IR matrix element of the first-variation effective potential.
        vEff1IR = cmplx(0., 0.)
        ! nmat has to be the dimension of the bras due to the concept of calcMEPotIR
        nmat = nv(1, ikpt) + atoms%nlotot
        call calcMEPotIR( stars, dimens, GbasVec(:, ilst(:nv(input%jspins, ikpt), ikpt, input%jspins)), &
          & GbasVec(:, ilst(:nv(input%jspins, ikpt), ikpt, input%jspins)), nv, grVeff0IR(:, idir), z(:, :, ikpt, 1),        &
          & z(:, :, ikpt, 1), gdp, nmat, ne(ikpt), nrBnd, ikpt, 1, ikpt, ngdp, vEff1IR, kpq2kPrVec ) !todo spin-relevant

        if (.false.) then
          do ibandK = 1, nrBnd 
            do ibandB = 1, ne(ikpt) ! pBand >> nBand
              write(2103, '(3(i8),2(f15.8))') idir, ibandK, ibandB, -aimag(vEff1IR(ibandB, ibandK))
            end do ! pBand
          end do ! nBand
        end if

        ! Calculate Pulay terms
        h0MTBv = cmplx(0.0, 0.0)
        s0MTBv = cmplx(0.0, 0.0)
        h0MTKv = cmplx(0.0, 0.0)
        s0MTKv = cmplx(0.0, 0.0)
        call calcHS0MT( atoms, usdus, td4HS0, ikpt, ikpt, iDtype, iDatom, ne, nobd, El, mCoefBar(:, :, iDatom), mCoef(:, :, iDatom), nRadFun, &
          & iloTable, nlo_atom, s0MTBv, h0MTBv )
        call calcHS0MT( atoms, usdus, td4HS0, ikpt, ikpt, iDtype, iDatom, ne, nobd, El, mCoef(:, :, iDatom), mCoefBar(:, :, iDatom), nRadFun, &
          & iloTable, nlo_atom, s0MTKv, h0MTKv )

        do ibandK = 1, nrBnd
          do ibandB = 1, nrBnd ! ibandB >> ibandK
!            write(2104, '(3(i8),2(f15.8))') idir, iBandK, iBandB, aimag(h0MTBv(ibandB, ibandK) - eig(ibandK, ikpt, 1) * s0MTBv(ibandB, ibandK))
!            write(2105, '(3(i8),2(f15.8))') idir, iBandK, iBandB, aimag(h0MTKv(ibandB, ibandK) - eig(ibandK, ikpt, 1) * s0MTKv(ibandB, ibandK))
!            write(2104, '(3(i8),2(f15.8))') idir, iBandK, iBandB, aimag(s0MTBv(ibandB, ibandK))
!            write(2105, '(3(i8),2(f15.8))') idir, iBandK, iBandB, aimag(s0MTKv(ibandB, ibandK))
!            write(2104, '(3(i8),2(f15.8))') idir, iBandK, iBandB, aimag(h0MTBv(ibandB, ibandK))
!            write(2105, '(3(i8),2(f15.8))') idir, iBandK, iBandB, aimag(h0MTKv(ibandB, ibandK))
          end do ! ibandB
        end do ! ibandK
        ! Calculate Surface integral

        surfIntTeps = cmplx(0.0, 0.0)
        surfInt = cmplx(0.0, 0.0)
        call calcSintKinEps(atoms, cell, kpts, qpts, iDtype, iDatom, nobd(ikpt, 1), ne(ikpt), ikpt, ikpt, 1, idir, &
          & nv(1, :), GbasVec, ilst, z(:, :, ikpt, 1),&
          & z(:, :, ikpt, 1), surfIntTeps, surfInt, kpq2kPrVec)
!todo   actually we should the IR surface integral here
!        allocate( overlap(dimens%neigd, maxval(nobd(:, :)) ) )
!        allocate( surfIntMT(dimens%neigd, maxval(nobd(:, :))) )
!        overlap = cmplx(0., 0.)
!        surfIntMT = cmplx(0., 0.)
!        call CalcSurfIntMT(atoms, dimens, sym, cell, kpts, input, usdus, ikpt, idir, dispAtInd, dispType, nv, Gbasvec, ilst, zKet, nobd, &
!          & ne, kveclo, nRadFun, iloTable, hFullNoAbcofBK, overlapNoAbcofBK, overlap, surfIntMT)
        do ibandK = 1, nrBnd
          do ibandB = 1, nrBnd ! ibandB >> ibandK
         !   write(2108, '(3(i8),2(f15.8))') idir, iBandK, iBandB, (surfIntTeps(ibandB, ibandK)                        &
         !   & + eig(ibandK, ikpt, 1) * surfInt(ibandB, ibandK) - surfIntVFast(ibandB, ibandK, idir))
            !write(2108, '(4(i8),2(f15.8))') ikpt, idir, iBandK, iBandB, (surfIntTeps(ibandB, ibandK) + surfIntVFast(ibandB, ibandK, idir))
!            write(2108, '(3(i8),2(f15.8))') idir, iBandK, iBandB, aimag(surfIntTeps(ibandB, ibandK) + surfIntVFast(ibandB, ibandK, idir))
!            write(2110, '(3(i8),2(f15.8))') idir, iBandK, iBandB, aimag(surfIntVFast(ibandK, ibandK, idir))
!            write(2109, '(3(i8),2(f15.8))') idir, iBandK, iBandB, aimag(h0MTKv(ibandB, ibandK) + h0MTBv(ibandB, ibandK) + surfIntTeps(ibandB, ibandK) + surfIntVFast(ibandB, ibandK, idir))
            !write(2108, '(4(i8),2(f15.8))') ikpt, idir, iBandK, iBandB, surfIntVFast(ibandB, ibandK, idir)
          end do ! ibandB
        end do ! ibandK
        !TODO Actually one should do that with Veff1 not with grad Veff1 to check Sternheimer better!

!        ! Sum up contributions
        do iband = 1, nrBnd
          !eig1(iband, ikpt, idir) = -vEff1IR(iband, iband) - h0MTBv(iband, iband) + eig(iband, ikpt, 1) * s0MTBv(iband, iband)      &
          !  & - h0MTKv(iband, iband) + eig(iband, ikpt, 1) * s0MTKv(iband, iband) - surfIntTeps(iband, iband)                        &
          !  & + eig(iband, ikpt, 1) * surfInt(iband, iband) - surfIntVFast(iband, iband, idir)
          eig1(iband, ikpt, idir) = -vEff1IR(iband, iband) + h0MTBv(iband, iband) - eig(iband, ikpt, 1) * s0MTBv(iband, iband)      &
            & + h0MTKv(iband, iband) - eig(iband, ikpt, 1) * s0MTKv(iband, iband) + surfIntTeps(iband, iband)                        &
            & - eig(iband, ikpt, 1) * surfInt(iband, iband) + surfIntVFast(iband, iband, idir)
!          eig1(iband, ikpt, idir) = -vEff1IR(iband, iband) + h0MTBv(iband, iband) - eig(iband, ikpt, 1) * s0MTBv(iband, iband)      &
!            & + h0MTKv(iband, iband) - eig(iband, ikpt, 1) * s0MTKv(iband, iband) + surfIntTeps(iband, iband)                        &
!            & - eig(iband, ikpt, 1) * surfInt(iband, iband) + surfIntVFast(iband, iband, idir)
        end do ! iband
      end do ! idir
    end do ! ikpt

  end subroutine CalcKSEnerg1stVar

  subroutine testIRbackFolding(atoms, cell, lathar, input, stars, dimens, Veff0, kpts, sym, ngdp, ne, nv, nobd, GbasVec, ilst, gdp, kpq2kPrVec, mapKpq2K, z, rho0IR, rho0MT, nmem_atom, clnu_atom, mlh_atom, memd_atom)

    use m_types
    use m_jpSternhHF, only : calcMEPotIR
    use m_jpConstants, only : iu
    use m_jpPotDensHelper, only : WarpIRPot, genPertPotDensGvecs, convertStar2G
    use m_jpGrVeff0, only : GenGrVeff0
    use mod_juPhonUtils, only : calcGrR2FinLH
    use m_jpVeff1, only : GenVeff1

    implicit none

    ! Type variable
    type(t_atoms),                  intent(in) :: atoms
    type(t_cell),                   intent(in) :: cell
    type(t_sphhar),                 intent(in) :: lathar
    type(t_input),                  intent(in) :: input
    type(t_stars),                  intent(in) :: stars
    type(t_dimension),              intent(in) :: dimens
    type(t_potential),              intent(in) :: Veff0
    type(t_kpts),                   intent(in) :: kpts
    type(t_sym),                    intent(in) :: sym

    ! Scalar variable
    integer,                        intent(in) :: ngdp
    integer,                        intent(in) :: memd_atom

    ! Array variable
    integer,                        intent(in) :: ne(:)
    integer,                        intent(in) :: nv(:, :)
    integer,                        intent(in) :: nobd(:, :)
    integer,                        intent(in) :: GbasVec(:, :)
    integer,                        intent(in) :: ilst(:, :, :)
    integer,                        intent(in) :: gdp(:, :)
    integer,                        intent(in) :: kpq2kPrVec(:, :, :)
    integer,                        intent(in) :: mapKpq2K(:, :) ! todo possible error here
    MCOMPLEX,                       intent(in) :: z(:,:,:,:)
    complex,                        intent(in) :: rho0IR(:, :)
    real,                           intent(in) :: rho0MT(:, 0:, :, :)
    integer,                        intent(in) :: nmem_atom(0:, :)
    complex,                        intent(in) :: clnu_atom(:, 0:, :)
    integer,                        intent(in) :: mlh_atom(:, 0:, :)

    ! Type parameters
    type(t_dimension)                          :: dimensSpecial

    ! Scalar parameters
    integer                                    :: nmat
    integer                                    :: ikpt
    integer                                    :: ikpq
    integer                                    :: iqpt
    integer                                    :: idir
    integer                                    :: iG
    integer                                    :: irecl
    integer                                    :: nkpt
    integer                                    :: nvd
    integer                                    :: nbasfcn
    integer                                    :: neigd
    integer                                    :: ikptSpecial
    integer                                    :: inband
    integer                                    :: ipband
    logical                                    :: harSw
    logical                                    :: extSw
    logical                                    :: xcSw
    logical                                    :: vExt1FullSw
    logical                                    :: vHarNum
    integer                                    :: iDatom
    integer                                    :: iDtype
    integer                                    :: ngpqdp
    integer                                    :: ngpqdp2km
!    logical                                    :: testGoldstein
!    logical                                    :: grRhoTermSw
!    integer                                    :: itype
!    integer                                    :: ilh
!    integer                                    :: imesh
!    integer                                    :: iatom
!    integer                                    :: ieqat
!    integer                                    :: oqn_l
!    integer                                    :: lm_pre
!    integer                                    :: mqn_m
!    integer                                    :: lm

    ! Array parameters
    MCOMPLEX,          allocatable             :: vEff1IR(:, :)
    complex,           allocatable             :: w_vEff1IR(:, :)
    complex,           allocatable             :: vEff0IRG(:)
    complex,           allocatable             :: vEff1Mat(:, :)
    real,              allocatable             :: el0Special(:, :, :)
    real,              allocatable             :: ello0Special(:, :, :)
    integer,           allocatable             :: neSpecial(:)
    integer,           allocatable             :: nvSpecial(:, :)
    integer,           allocatable             :: nvExtended(:, :)
    integer,           allocatable             :: GbasVecSpecial(:, :)
    real,              allocatable             :: eigSpecial(:,:,:)
    integer,           allocatable             :: ilstSpecial(:, :, :)
    integer,           allocatable             :: kvecloSpecial(:,:)
    MCOMPLEX,          allocatable             :: zSpecial(:,:,:,:)
    integer,           allocatable             :: kpq2kPrVecSpecial(:, :, :)
!    complex,           allocatable             :: grVeff0IR(:, :)
!    real,              allocatable             :: r2Rho0MT(:, :, :, :)
!    complex,           allocatable             :: r2GrRho0MT(:, :, :, :)
!    complex,           allocatable             :: grVeff0MT(:, :, :, :)
!    complex,           allocatable             :: grRho0MT(:, :, :, :)
!    real,              allocatable             :: gausWts(:) ! gaussian weights belonging to gausPts
!    complex,           allocatable             :: ylm2(:, :)
!    complex,           allocatable             :: ylm1(:, : )
!    real,              allocatable             :: dKernMTGPts(:, :, :)
!    complex,           allocatable             :: grVxcIRKern(:)
    real                                       :: Gext(3)
    real                                       :: qpoint(3)
    complex,           allocatable            :: grRho0MT(:, :, :, :)
    complex,           allocatable            :: vxc1IRKern(:)
    complex,           allocatable            :: ylm(:, :)
    real,              allocatable            :: gWghts(:) ! gaussian weights belonging to gausPts
    real,              allocatable            :: dKernMTGPts(:, :, :)
    complex,           allocatable            :: vExt1MTtemp(:, :, :, :)
    complex,           allocatable            :: rho1IR(:, :, :)
    complex,           allocatable            :: rho1MT(:, :, :, :, :)
    integer,  allocatable                :: gpqdp(:, :)
    integer,  allocatable                :: gpqdp2Ind(:, :, :)
      complex,           allocatable              :: rho0IRDummy(:, :)
      complex,           allocatable              :: rho0MTDummy(:, :, :, :)
    integer                              :: gpqdp2iLim(2, 3)


!    ! The factor r^2 has beeen divided out so that rho0MT is only rho0MT and not r^2 rho0MT as it is done in Fleur. The factor
!    ! sqrt(4pi) for the zeroth component was already dividied out when constructing rho0MT in cdnmt routine in Fleur. Here to
!    ! improve stability of the gradient routine we derive r2Rho0MT and divide out the r^2 again later. Doing so avoids the
!    ! subtraction of small numbers close to the core.
!    allocate( r2Rho0MT( atoms%jmtd, 0:lathar%nlhd, atoms%ntype, 1) )
!    allocate( grRho0MT( atoms%jmtd, ( atoms%lmaxd + 2 )**2, atoms%nat, 3) )
!    r2Rho0MT(:, :, :, :) = 0.
!    grRho0MT(:, :, :, :) = cmplx(0., 0.)
!
!    do itype = 1, atoms%ntype
!      do ilh = 0, lathar%nlhd
!        do imesh = 1, atoms%jri(itype)
!          r2Rho0MT(imesh, ilh, itype, 1) = rho0MT(imesh, ilh, itype, 1) * atoms%rmsh(imesh, itype) * atoms%rmsh(imesh, itype)
!        end do
!      end do
!    end do
!
!    call calcGrR2FinLH( atoms, lathar, clnu_atom, nmem_atom, mlh_atom, r2Rho0MT(:, :, :, 1), r2GrRho0MT )
!
!    do idir = 1, 3
!      iatom = 0
!      do itype = 1, atoms%ntype
!        do ieqat = 1, atoms%neq(itype)
!          iatom = iatom + 1
!          do oqn_l = 0, atoms%lmax(itype)
!            lm_pre = oqn_l * (oqn_l + 1) + 1
!            do mqn_m = -oqn_l, oqn_l
!              lm = lm_pre + mqn_m
!              do imesh = 1, atoms%jri(itype)
!                grRho0MT(imesh, lm, iatom, idir) = r2GrRho0MT(imesh, lm, iatom, idir) / atoms%rmsh(imesh, itype)**2
!              end do ! imesh
!            end do ! mqn_m
!          end do ! oqn_l
!        end do ! ieqat
!      end do ! itype
!    end do ! idir
!
!    harSw = .true.
!    extSw = .true.
!    xcSw = .false.
!    testGoldstein = .false.
!    grRhoTermSw = .false. ! this has to be equal to vExtFull going into Veff1
!    call GenGrVeff0( atoms, cell, lathar, stars, dimens, memd_atom, ngdp, harSw, extSw, xcSw, gdp, &
!      & rho0IR( :, 1 ), rho0MT(:, :, :, 1), nmem_atom, mlh_atom, clnu_atom, grVeff0IR, grVeff0MT, grRho0MT, gausWts, ylm1,&
!      & ylm2, dKernMTGPts, grVxcIRKern, testGoldstein, grRhoTermSw ) ! add spin

    ! Calculate the IR matrix element of the first-variation effective potential.
    allocate( vEff1IR(ngdp, 3) )
    allocate(vEff0IRG(ngdp))
    allocate( rho1IR( ngdp, 3, atoms%nat ) )
    allocate( rho1MT( atoms%jmtd, (atoms%lmaxd + 1)**2, atoms%nat, 3, atoms%nat ) )
    allocate( grRho0MT( atoms%jmtd, ( atoms%lmaxd + 1 )**2, atoms%nat, 3) )
    allocate( vxc1IRKern(ngdp) )
    allocate( gWghts(dimens%nspd), ylm(dimens%nspd, ( atoms%lmaxd + 1 )**2), dKernMTGPts(dimens%nspd, atoms%jmtd, atoms%nat))
     allocate(  rho0IRDummy(ngdp, 1) )
     allocate( rho0MTDummy( atoms%jmtd, (atoms%lmaxd + 1)**2, atoms%nat, 1) )
     rho0IRDummy(:, :) = cmplx(0., 0.)
     rho0MTDummy(:, :, :, :) = cmplx(0., 0.)
    vxc1IRKern(:) = cmplx(0., 0.)
    grRho0MT(:, :, :, :) = cmplx(0., 0.)
    vEff1IR(:, :) = 0.
    ylm(:, :) = cmplx(0., 0.)
    dKernMTGPts(:, :, :) = 0.
    gWghts(:) = 0.
    harSw = .false.
    extSw = .true.
    iDatom = 1
    iDtype = 1
    xcSw = .false.
    qpoint(:) = [0., 0., 0.5]
    iqpt = 3
    ! There is nothing canceling away as in Sternheimer therefore we need the full contribution of the external potential
    vExt1FullSw = .true.
    vHarNum=.false.
    ! todo we don't need rho0IR in the argument list
    ! todo attention we do not have a factor r^2 here
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !TODO ATTENTION REACTIVATE GENVEFF1 BEFORE COMMIT
    write(*, *) 'BUG for this test~~~~~~ shifted G set really correct?'
    call genPertPotDensGvecs( stars, cell, input, ngpqdp, ngpqdp2km, qpoint, gpqdp, gpqdp2Ind, gpqdp2iLim )
    call GenVeff1( stars, cell, atoms, dimens, harSw, extSw, xcSw, vExt1FullSw, ngdp, qpoint, rho0IRDummy, rho0MTDummy, rho1IR(:, :, iDatom), &
      & rho1MT(:, :, :, :, iDatom), grRho0MT, gdp, vEff1IR, vExt1MTtemp, vxc1IRKern, ylm, dKernMTGPts, gWghts, iDatom, iDtype, iqpt,&
      & ngpqdp, gpqdp, vHarNum ) ! add spin
    !vEff0IRG(:) = 0
    !call convertStar2G( Veff0%vpw(:, 1), vEff0IRG, stars, ngdp, gdp )
    !do idir = 1, 3
    !  do iG = 1, ngdp
    !    Gext(1:3) = matmul( cell%bmat(1:3, 1:3), gdp(1:3, iG) )
    !    ! Remember to first apply the gradient then warp not vice versa!
    !    vEff1IR(iG, idir) = -iu * Gext(idir) * vEff0IRG(iG)
    !  end do ! iG
    !end do ! idir
!    allocate(w_vEff1IR(ngdp, 3))
!    do idir = 1, 3
!      call warpIRPot(stars, ngdp, idir, gdp, vEff1IR, w_vEff1IR(:, idir))
!    end do
    write(2152, '(2f15.8)') w_vEff1IR

    do idir = 1, 3
      ! Test it with k + q = k' and backfolding
      ikpt = 4
      ikpq = mapKpq2K(ikpt, iqpt) ! todo should be 2
      write(*, *) 'ikpq', ikpq
      allocate( vEff1Mat(ne(ikpq), nobd(ikpt, 1)) )
      ! nmat has to be the dimension of the bras due to the concept of calcMEPotIR
      nmat = nv(1, ikpq) + atoms%nlotot
      call calcMEPotIR( stars, dimens, GbasVec(:, ilst(:nv(input%jspins, ikpq), ikpq, input%jspins)),                         &
        & GbasVec(:, ilst(:nv(input%jspins, ikpt), ikpt, input%jspins)), nv, vEff1IR(:,idir), z(:, :, mapKpq2K(ikpt, iqpt), 1),    &
        & z(:, :, ikpt, 1), gdp, nmat, ne(ikpq), nobd(ikpt, 1), ikpt, iqpt, ikpq, ngdp, vEff1Mat, kpq2kPrVec ) !todo spin-relevant

      do inband = 1, nobd(ikpt, 1)
        do ipband = 1, ne(ikpq)
          write(2149, '(3i8,2f15.8)') idir, ipband, inband, vEff1Mat(ipband, inband)
        end do ! ipband
      end do ! inband
      deallocate(vEff1Mat)

      !todo set them correctly
      irecl = 307616
      nkpt = 1
      nvd = 137
      nbasfcn = nvd + atoms%nat * atoms%nlod * ( 2 * atoms%llod + 1)
      neigd = 137
      ikptSpecial = 1
      ikpq = kpts%nkpt + 1

      !irecl = atoms%ntype * (atoms%lmaxd + 1 + atoms%nlod) + 2
      !irecl = 8 * (irecl + 4 + dimens%neigd * (2 * dimens%nbasfcn + 1))
      !irecl = irecl + 4 * atoms%nat * atoms%nlod * (2 * atoms%llod + 1)
      write(*, *) 307616, irecl
      call readSpecialEig(atoms, irecl, nkpt, nvd, nbasfcn, neigd, el0Special, ello0Special, neSpecial, nvSpecial, GbasVecSpecial,&
        & eigSpecial, ilstSpecial, kvecloSpecial, zSpecial )

      ! Test with k + q and no backfolding
      allocate( vEff1Mat(neSpecial(ikptSpecial), nobd(ikpt, 1)) )
      vEff1Mat = cmplx(0., 0.)
      ! nmat has to be the dimension of the bras due to the concept of calcMEPotIR
      nmat = nvSpecial(1, ikptSpecial) + atoms%nlotot
      allocate( kpq2kPrVecSpecial(3, kpts%nkpt, kpts%nkpt) )
      kpq2kPrVecSpecial = cmplx(0., 0.)
      allocate( nvExtended( input%jspins, kpts%nkpt + 1 ) )
      nvExtended(1, 1:kpts%nkpt) = nv(1, 1:kpts%nkpt)
      nvExtended(1, kpts%nkpt + 1:kpts%nkpt + nkpt) = nvSpecial(1, 1:nkpt)

      dimensSpecial%jspd = dimens%jspd
      dimensSpecial%nspd = dimens%nspd
      dimensSpecial%nvd  = nvd
      dimensSpecial%nv2d = dimens%nv2d
      dimensSpecial%neigd= neigd
      dimensSpecial%neigd2 = dimens%neigd2
      dimensSpecial%ncvd = dimens%ncvd
      dimensSpecial%nn2d = dimens%nn2d
      dimensSpecial%nn3d = dimens%nn3d
      dimensSpecial%nstd = dimens%nstd
      dimensSpecial%msh = dimens%msh
      dimensSpecial%lmd = dimens%lmd
      dimensSpecial%lmplmd = dimens%lmplmd
!      dimensSpecial%nbasfcn = nbasfcn
      dimensSpecial%nbasfcn = dimens%nbasfcn

      call calcMEPotIR( stars, dimensSpecial, GbasVecSpecial(:, ilstSpecial(:nvSpecial(input%jspins, ikptSpecial), ikptSpecial, input%jspins)),  &
        & GbasVec(:, ilst(:nv(input%jspins, ikpt), ikpt, input%jspins)), nvExtended, vEff1IR(:,idir), zSpecial(:, :, ikptSpecial, 1),       &
        & z(:, :, ikpt, 1), gdp, nmat, neSpecial(ikptSpecial), nobd(ikpt, 1), ikpt, iqpt, ikpq, ngdp, vEff1Mat, kpq2kPrVecSpecial ) !todo spin-relevant

      do inband = 1, nobd(ikpt, 1)
        do ipband = 1, neSpecial(ikptSpecial)
          write(2150, '(3i8,2f15.8)') idir, ipband, inband, vEff1Mat(ipband, inband)
        end do ! ipband
      end do ! inband
      deallocate(vEff1Mat, kpq2kPrVecSpecial, nvExtended)
    end do ! idir



  end subroutine testIRbackFolding

  subroutine readSpecialEig(atoms, irecl, nkpt, nvd, nbasfcn, neigd, el0Special, ello0Special, neSpecial, nvSpecial, GbasVecSpecial,&
      & eigSpecial, ilstSpecial, kvecloSpecial, zSpecial )

    use mod_juPhonUtils, only : fopen, fclose
    use m_juDFT_NOstopNO, only : juDFT_error, juDFT_warn
    use m_types, only : t_atoms

    implicit none

    ! Type parameters
    type(t_atoms),         intent(in) :: atoms

    ! Scalar parameters
    integer,               intent(in)  :: irecl
    integer,               intent(in)  :: nkpt
    integer,               intent(in)  :: nvd
    integer,               intent(in)  :: nbasfcn
    integer,               intent(in)  :: neigd

    ! Array parameters
    real,     allocatable, intent(out) :: el0Special(:, :, :)
    real,     allocatable, intent(out) :: ello0Special(:, :, :)
    integer,  allocatable, intent(out) :: neSpecial(:)
    integer,  allocatable, intent(out) :: nvSpecial(:, :)
    integer,  allocatable, intent(out) :: GbasVecSpecial(:, :)
    real,     allocatable, intent(out) :: eigSpecial(:,:,:)
    integer,  allocatable, intent(out) :: ilstSpecial(:, :, :)
    integer,  allocatable, intent(out) :: kvecloSpecial(:,:)
    MCOMPLEX, allocatable, intent(out) :: zSpecial(:,:,:,:)

    ! Scalar variables
    integer                            :: ierror
    real                               :: rdum
    integer                            :: intDump                 ! integer temporary variable
    integer                            :: ikpt
    integer                            :: indexG                  ! accumulated index
    integer                            :: isp
    integer                            :: iG
    logical                            :: addG                    ! Helping logical for compressed storage of basis vectors G
    integer                            :: ifind                    ! Helping logical for compressed storage of basis vectors G
    integer                            :: irec                   ! Helping logical for compressed storage of basis vectors G
    integer                            :: ii

    ! Array variables
    real                               :: bk(3)
    integer, allocatable               :: GbasVec_temp(:, :)      ! Helping array for compressed storage of basis vectors G
    integer, allocatable               :: GbasVec_eig(:, :, :, :) ! temporary storage for reading basis vectors G from eig
    integer, allocatable               :: k1(:)                   ! temporary x coordinates of basis vectors for a k-point
    integer, allocatable               :: k2(:)                   ! temporary y coordinates of basis vectors for a k-point
    integer, allocatable               :: k3(:)                   ! temporary z coordinates of basis vectors for a k-point

    write(*, *) 'irecl', irecl
    call fopen(1000, name='eigSpecial', status='old', action='read', form='unformatted', access='direct', recl=irecl)

    ! Allocate quantities to be read out from eig file
    allocate( nvSpecial(1, nkpt), neSpecial(nkpt), k1(nvd), k2(nvd), k3(nvd) )

    allocate( zSpecial(nbasfcn, neigd, nkpt, 1), stat=ierror )
    if( ierror /=  0 ) then
      call juDFT_error( 'Error in allocation of eigenvectors z', calledby='Read_eig', hint='Check for having enough storage' )
    end if

    allocate( eigSpecial(neigd, nkpt, 1), stat=ierror )
    if( ierror /=  0 ) then
      call juDFT_error( 'Error in allocation of eigen energies eig', calledby='Read_eig', hint='Check for having enough storage' )
    end if

    allocate( el0Special(0:atoms%lmaxd, atoms%ntype, 1), ello0Special(atoms%nlod, atoms%ntype, 1), &
      & stat=ierror )
    if( ierror /=  0 ) then
      call juDFT_error( 'Error in allocation of energy parameters el0 ello0', calledby='Read_eig', &
        & hint='Check for having enough storage' )
    end if

    allocate( kvecloSpecial(atoms%nlotot, nkpt), stat=ierror )
    if( ierror /=  0 ) then
      call juDFT_error( 'Error in allocation of G-basis vectors for LOs kveclo', calledby='Read_eig', &
        & hint='Check for having enough storage' )
    end if

    allocate( ilstSpecial(nvd, nkpt, 1), GbasVec_eig(3, nvd, nkpt, 1), GbasVec_temp(3, 0), GbasVecSpecial(3, 0) )

    ! Actually reading out eig file line after line.
    indexG = 0
!    GbasVec_eig = 0
    do isp = 1, 1
      do ikpt = 1, nkpt
        !irec = nkpt * (isp - 1) + ikpt
        irec = 1
        read(1000, rec=irec) el0Special(:,:,isp), (rdum,ii=1,2), ello0Special(:,:,isp), bk, rdum, neSpecial(ikpt), nvSpecial(isp, ikpt), intDump, &
          & (eigSpecial(ii, ikpt, isp), ii=1, neigd), k1(:), k2(:), k3(:), kvecloSpecial(:, ikpt), zSpecial(:, :, ikpt, isp)
        GbasVec_eig(1, :nvSpecial(1, ikpt), ikpt, isp) = k1(:nvSpecial(1, ikpt))
        GbasVec_eig(2, :nvSpecial(1, ikpt), ikpt, isp) = k2(:nvSpecial(1, ikpt))
        GbasVec_eig(3, :nvSpecial(1, ikpt), ikpt, isp) = k3(:nvSpecial(1, ikpt))

      end do ! ikpt
    end do ! isp

    deallocate( k1, k2, k3 )

    if ( neigd /= nvd + atoms%nlotot ) then
      call juDFT_warn( 'Number of eigenvalues is not equals number of basis functions. You might miss relevant occupied bands!', &
        &calledby='Read_eig', hint='Set neigd = nvd + nlotot within fl7para file' )
    end if

    ! Store basis vectors G in a compressed way and not redudantly as they occur multiple times for different k-points. The array ilst
    ! stores the index where the basis vector is stored if it has occured for a different k-point or spin yet.
    indexG = 0
    ! Initialize ilstSpecial generally to -1 to indicate unvalid value
    ilstSpecial = -1
    do isp = 1, 1
      do ikpt = 1, nkpt
        do iG = 1, nvSpecial(isp, ikpt)
          if ( uBound( GbasVecSpecial, 2 ) == 0 ) then
            addG = .true.
          else
            do ifind =  1, uBound( GbasVecSpecial, 2 )

              if ( .not. all( abs( GbasVecSpecial(:, ifind) - GbasVec_eig(:, iG, ikpt, isp) ) <= 1E-8 ) ) then
                addG = .true.
              else
                ilstSpecial(iG, ikpt, isp) = ifind
                addG = .false.
                exit
              endif
            enddo ! ifind
          endif

          if( addG ) then
            if ( uBound( GbasVecSpecial, 2 ) == 0 ) then
              deallocate( GbasVecSpecial, GbasVec_temp )
              allocate( GbasVecSpecial(3, 1), GbasVec_temp(3, 1) )
            else
              GbasVec_temp = GbasVecSpecial
              deallocate( GbasVecSpecial)
              allocate( GbasVecSpecial(3, indexG + 1) )
              GbasVecSpecial(:, :indexG) = GbasVec_temp
            endif

            indexG = indexG + 1
            ilstSpecial(iG, ikpt, isp) = indexG
            GbasVecSpecial(1, indexG) = GbasVec_eig(1, iG, ikpt, isp)
            GbasVecSpecial(2, indexG) = GbasVec_eig(2, iG, ikpt, isp)
            GbasVecSpecial(3, indexG) = GbasVec_eig(3, iG, ikpt, isp)
            if( uBound( GbasVecSpecial, 2 ) /= 0 ) then
              ! resize GbasVec_temp to current size of GbasVec
              deallocate( GbasVec_temp )
              allocate( GbasVec_temp(3, indexG) )
            end if
          end if
        end do ! iG
      end do ! ikpt
    end do ! isp

    call fclose(1000)

  end subroutine readSpecialEig

  ! NOTE: The Gbas vectors have to be already sorted by the array ilst.
  !---------------------------------------------------------------------------------------------------------------------------------
  !> @author
  !> Christian-Roman Gerhorst, Gregor Michaliczek
  !>
  !> @brief
  !> Calls routines to set-up Sternheimer equation and solves it for first-order wavefunction expansion coefficients.
  !>
  !> @details
  !>
  !> @note
  !> We have to consider all (occupied and unoccupied) bands p in the bra and only the occupied bands n in the kets.
  !> Using the OEP approach of Markus Betzinger, it should be possible to also only consider the occupied bands in the bras leading to
  !> significant runtime enhancements.
  !---------------------------------------------------------------------------------------------------------------------------------
  subroutine calcSfVeff( atoms, cell, stars, dimens, qpts, ikpt, ikpq, iqpt, ngdp, nrBraBands, nrKetBands, iDtype, iDatom, iDdir, nv, &
      &GbasBra, GbasKet, gdp, vEff0IR, zBra, zKet, surfInt, kpq2kPrVec )

#include "recycledRoutines/cpp_arch.h"

    use m_types
    use m_ylm_old
    use m_sphbes
    use m_JPConstants, only : iu, fpi, tpi, c_im

    implicit none

    include "fftw3.f" ! include FFTW constants

    ! Type parameter
    type(t_atoms),                  intent(in) :: atoms
    type(t_cell),                   intent(in) :: cell
    type(t_stars),                  intent(in) :: stars
    type(t_dimension),              intent(in) :: dimens
    type(t_kpts),                   intent(in) :: qpts

    ! Scalar parameter
    integer,                        intent(in) :: ikpt
    integer,                        intent(in) :: ikpq
    integer,                        intent(in) :: iqpt
    integer,                        intent(in) :: ngdp
    integer,                        intent(in) :: nrBraBands
    integer,                        intent(in) :: nrKetBands
    integer,                        intent(in) :: iDtype
    integer,                        intent(in) :: iDatom
    integer,                        intent(in) :: iDdir

    ! Array parameter
    integer,                        intent(in)  :: nv(:, :)
    integer,                        intent(in)  :: GbasKet(:, :)
    integer,                        intent(in)  :: GbasBra(:, :)
    integer,                        intent(in)  :: gdp(:, :)
    complex,                        intent(in)  :: vEff0IR(:) ! thetaV(G) in stars expansion
    MCOMPLEX,                       intent(in)  :: zBra(dimens%nbasfcn,dimens%neigd) ! basis functions in PW expansion
    MCOMPLEX,                       intent(in)  :: zKet(dimens%nbasfcn,dimens%neigd) ! basis functions in PW expansion
    integer,                        intent(in)  :: kpq2kPrVec(:, :, :)
    complex,                        intent(out) :: surfInt(:, :)

    ! Local Scalars
  !  integer                                     :: iv
    integer                                     :: il
    integer                                     :: im
    integer                                     :: in
    integer                                     :: sk1d
    integer                                     :: sk2d
    integer                                     :: sk3d
    integer                                     :: ifft1ds
    integer                                     :: ifft2ds
    integer                                     :: ifft3ds
    integer                                     :: ifftds
  !  integer                                     :: iG1
    integer                                     :: iG
    integer                                     :: iBas
    integer                                     :: idir
    integer                                     :: xIndex
    integer                                     :: xis
    integer                                     :: xil
    integer                                     :: yIndex
    integer                                     :: yis
    integer                                     :: yil
    integer                                     :: zIndex
    integer                                     :: zis
    integer                                     :: zil
    integer                                     :: smallIndex
    integer                                     :: largeIndex
    integer                                     :: iBband
    integer                                     :: iKband
    integer                                     :: imesh
    integer                                     :: gridIndex
    integer                                     :: lm
    real                                        :: raf
    real                                        :: testG
    integer*8                                      backwardPlan
    integer*8                                      backwardPlanK
    integer*8                                      backwardPlanB
    complex                                     :: prf
    complex                                     :: factor
    complex                                     :: zeta
    complex                                     :: sumt
    integer, parameter :: boxSizesMaxIndex = 200

    ! Local Arrays
    integer                                     :: maxGk(3)
    integer                                     :: maxGb(3)
    integer                                     :: maxGp(3)
    integer                                     :: nfft(3)
    integer                                     :: gMirr(3)
    real                                        :: Gpqvec(3)
    real                                        :: Greal(3)
    complex                                     :: ylm(4)
    real                                        :: fj(0 : 1)
    complex,           allocatable              :: tempGrid(:)
    integer,           allocatable              :: igfft(:)
    complex,           allocatable              :: VFFTBox(:)
    integer,           allocatable              :: iv1dBra(:,:)
    integer,           allocatable              :: iv1dKet(:,:)
    complex,           allocatable              :: zFFTkBox(:)
    complex,           allocatable              :: zFFTbBox(:)
    complex,           allocatable              :: zVzFFTBox(:, :, :)
    complex,           allocatable              :: preRaileigh(:)
    integer(kind=8),   allocatable              :: forwardPlanR(:, :) ! todo test IR routine whether this form of writing the kind is correct
    integer                                     :: fftBoxSizes(1:boxSizesMaxIndex) !array with optimal FFT box sizes larger or equal to array index.

    intrinsic isign,real,cmplx,aimag,conjg

    external dfftw_plan_dft_3d
    external dfftw_execute_dft
    external dfftw_destroy_plan


    !***********************
    !* Initialization step *
    !***********************

    ! Optimization of the FFT mesh sizes (Measured by Gregor Michaliczek) so that the FFT from fftw runs faster
    fftBoxSizes = [  1,  2,  3,  4,  5,  6,  7,  8,  9, 10, 11, 12, 13, 14, 15, 16, 20, 20, 20, 20, 21, 25, 25, 25, 25, 26, 27, 28, &
                   & 30, 30, 32, 32, 36, 36, 36, 36, 40, 40, 40, 40, 42, 42, 44, 44, 45, 48, 48, 48, 50, 50, 54, 54, 54, 54, 56, 56, &
                   & 60, 60, 60, 60, 64, 64, 64, 64, 66, 66, 72, 72, 72, 72, 72, 72, 75, 75, 75, 80, 80, 80, 80, 80, 84, 84, 84, 84, &
                   & 90, 90, 90, 90, 90, 90, 96, 96, 96, 96, 96, 96, 98, 98,100,100,104,104,104,104,105,108,108,108,112,112, 112,112,&
                   &120,120,120,120,120,120,120,120, 128,128,128,128,128,128,128,128,130,130, 132,132,135,135,135,144,144,144,144,144, &
                   &144,144,144,144,150,150,150,150,150,150, 160,160,160,160,160,160,160,160,160,160, 162,162,168,168,168,168,168,168, &
                   & 180,180, 180,180,180,180,180,180,180,180,180,180, 192,192,192,192,192,192,192,192,192,192, 192,192,200,200,200,200, &
                   & 200,200,200,200 ]

    surfInt = cmplx(0.0, 0.0)

    ! Measure maximal absolute coordinate of G-basis vector and G-potential vector
    ! todo it might be that we do not hit the same basis vector only max coordinates from different basis vectors.
    maxGk = 0
    do iBas = 1, nv(1,ikpt)
      il = GbasKet(1, iBas)
      im = GbasKet(2, iBas)
      in = GbasKet(3, iBas)
      testG = abs(il)
      if (testG > maxGk(1)) then
        maxGk(1) = ceiling(testG)
      end if
      testG = abs(im)
      if (testG > maxGk(2)) then
        maxGk(2) = ceiling(testG)
      end if
      testG = abs(in)
      if (testG > maxGk(3)) then
        maxGk(3) = ceiling(testG)
      end if
    end do

    maxGb = 0
    do iBas = 1, nv(1,ikpq)
    ! todo do we really need the kpq2kPrVec here?
      il = GbasBra(1, iBas) + kpq2kPrVec(1, ikpt, iqpt)
      im = GbasBra(2, iBas) + kpq2kPrVec(2, ikpt, iqpt)
      in = GbasBra(3, iBas) + kpq2kPrVec(3, ikpt, iqpt)
      testG = abs(il)
      if (testG > maxGb(1)) then
        maxGb(1) = ceiling(testG)
      end if
      testG = abs(im)
      if (testG > maxGb(2)) then
        maxGb(2) = ceiling(testG)
      end if
      testG = abs(in)
      if (testG > maxGb(3)) then
        maxGb(3) = ceiling(testG)
      end if
    end do

    maxGp = 0
    do iG = 1, ngdp
      il = gdp(1, iG)
      im = gdp(2, iG)
      in = gdp(3, iG)
      testG = abs(il)
      if (testG > maxGp(1)) then
        maxGp(1) = ceiling(testG)
      end if
      testG = abs(im)
      if (testG > maxGp(2)) then
        maxGp(2) = ceiling(testG)
      end if
      testG = abs(in)
      if (testG > maxGp(3)) then
        maxGp(3) = ceiling(testG)
      end if
    end do

    ! We have the product of two wavefunctions and a potential so we need a Fourier mesh of 2 * kmax + Gmax...
    ! todo the kpq2kprvec is used for maxGK but all in all it should not be used for this sum here
    sk1d = maxGb(1) + maxGk(1) + maxGp(1)! bis 2 kmax(=maxG) + Gmax
    sk2d = maxGb(2) + maxGk(2) + maxGp(2)
    sk3d = maxGb(3) + maxGk(3) + maxGp(3)

    ! ...and we have to account for both the positive and negative axis and zero.
    ifft1ds=(2*sk1d+1)
    ifft2ds=(2*sk2d+1)
    ifft3ds=(2*sk3d+1)

    ! Optimize the size of the FFT mesh to run the FFT quicker
    !todo comment this out, this might be a reason that the calculation fails
!    if (ifft1ds <= boxSizesMaxIndex) then ! kann ich drin lassen
!       ifft1ds = fftBoxSizes(ifft1ds)
!    end if
!    if (ifft2ds <= boxSizesMaxIndex) then
!       ifft2ds = fftBoxSizes(ifft2ds)
!    end if
!    if (ifft3ds <= boxSizesMaxIndex) then
!       ifft3ds = fftBoxSizes(ifft3ds)
!    end if

    ! We reduce 3D to a one-dimensional array of size ifftds
    ifftds=ifft1ds*ifft2ds*ifft3ds

    !******************************************************************************************
    !* Step 1: Put the unwarped interstitial potential onto the main FFT mesh and perform a FFT
    !******************************************************************************************

    ! Generate mapping array storing the index in the FFT mesh for every G-vector
    !todo Reasons of convention we should check whether a -G instead of a G is more apropriate here
    ! todo it is strange that a gvector like (3 , 3, -3) is only partly mirrored
    ! todo if we construct this mapping array igfft, we use k1d k2d, k3d which we change in fl7para is this a mistake?
    allocate( igfft(ngdp) )
    nfft = [ 3 * stars%k1d, 3 * stars%k2d, 3 * stars%k3d ]
    gMirr = 0
    igfft = 0
    do iG = 1, ngdp
      do idir = 1, 3
        gMirr(idir) = gdp(idir, iG)
        if ( gMirr(idir) < 0 ) then
          gMirr(idir) = gMirr(idir) + nfft(idir)
        end if
      end do
      igfft(iG) = gMirr(1) + nfft(1) * gMirr(2) + nfft(1) * nfft(2) * gMirr(3)
      !write(1013, '(4(i8))')  igfft(iG), gdp(:, iG)
      gMirr = 0
    end do

    ! Initialize FFT from fftw library
    allocate (VFFTBox(0:ifftds - 1))
    VFFTBox = cmplx(0.0,0.0)
    CALL dfftw_plan_dft_3d(backwardPlan, ifft1ds, ifft2ds, ifft3ds, VFFTBox, VFFTBox, FFTW_BACKWARD, FFTW_ESTIMATE)
    VFFTBox = cmplx(0.0,0.0)

    ! The mapping indices from igfft are only valid for a FFT mesh of the size the variable tempGrid features (using k1d, k2d, k3d) so,
    ! in a first step, we put it onto this grid.
    allocate (tempGrid(0:27 * stars%k1d * stars%k2d * stars%k3d - 1))
    tempGrid = (0.0,0.0)
    do iG = 1, ngdp
       tempGrid(igfft(iG)) = vEff0IR(iG)
    end do

    ! Migrate potential from temporary grid tempGrid to a grid of the final size VFFTBox
    ! The coordinates of G are sorted that the positives come first starting from the beginning and the negative start in reverse order
    ! from ifftNds
    !todo mapping correctly?
    do xIndex = -min(sk1d, stars%k1d), min(sk1d, stars%k1d) ! aufpassen bzgl. indices !TODO Loop indices anpassen!!!!
       xis = xIndex
       xil = xIndex
       if (xIndex < 0) then
          xis = xis + ifft1ds
          xil = xil + 3 * stars%k1d
       end if
       do yIndex = -min(sk2d, stars%k2d), min(sk2d, stars%k2d)
          yis = yIndex
          yil = yIndex
          if (yIndex < 0) then
             yis = yis + ifft2ds
             yil = yil + 3 * stars%k2d
          end if
          do zIndex = -min(sk3d, stars%k3d), min(sk3d, stars%k3d)
             zis = zIndex
             zil = zIndex
             if (zIndex < 0) then
                zis = zis + ifft3ds
                zil = zil + 3 * stars%k3d
             end if
             smallIndex = (xis + ifft1ds * yis + ifft1ds * ifft2ds * zis)
             largeIndex = (xil + 3 * stars%k1d * yil + 9 * stars%k1d * stars%k2d * zil)
             VFFTBox(smallIndex) = tempGrid(largeIndex)
          end do
       end do
    end do

    deallocate (tempGrid)

    ! perform a complex 3D FFT on "VFFTBox" from reciprocal space to real space
    call dfftw_execute_dft(backwardPlan, VFFTBox, VFFTBox)
    call dfftw_destroy_plan(backwardPlan)


    !*********************************************************************
    !* Step 2: Generate basis functions on real-space grid by performing *
    !*         a FFT on z. Then calculate matrix elements H_ij.          *
    !*********************************************************************
    !todo beware of the minus sign in front of k_i, nv from ket
    ! frage der konsistenz der definition von ebenen wellen von fleur und fourier bibliothek und wie ist vorwärts und rückwärts trafo definiert
    ! Create mapping array to put ket-wavefunction coefficient, which is dependent on G, onto FFT mesh
    allocate(iv1dKet(dimens%nvd, 1))
    do iBas = 1, nv(1, ikpt)
       il = GbasKet(1, iBas)
       im = GbasKet(2, iBas)
       in = GbasKet(3, iBas)
       if (il < 0) then
         il = il + ifft1ds ! is ifft1d correct form for ket?
       end if
       if (im < 0) then
          im = im + ifft2ds
       end if
       if (in < 0) then
          in = in + ifft3ds
       end if
       iv1dKet(iBas,1) = il+ifft1ds*im+ifft1ds*ifft2ds*in ! spin
       !write(1018, '(4(i8))') iv1dKet(iBas, 1), GbasKet(1, iBas), GbasKet(2, iBas), GbasKet(3, iBas)
    end do

    ! Create mapping array to put bra-wavefunction coefficient, which is dependent on G, onto FFT mesh
    allocate(iv1dBra(dimens%nvd, 1))
    !todo is kpq2kPrVec correct
    do iBas = 1, nv(1, ikpq)
       il = GbasBra(1, iBas) + kpq2kPrVec(1, ikpt, iqpt)
       im = GbasBra(2, iBas) + kpq2kPrVec(2, ikpt, iqpt)
       in = GbasBra(3, iBas) + kpq2kPrVec(3, ikpt, iqpt)
       if (il < 0) then
          il = il + ifft1ds
       end if
       if (im < 0) then
          im = im + ifft2ds
       end if
       if (in < 0) then
          in = in + ifft3ds
       end if
       iv1dBra(iBas,1) = il+ifft1ds*im+ifft1ds*ifft2ds*in ! spin
!       write(1023, '(4(i8))') iv1dBra(iBas, 1), GbasBra(1, iBas), GbasBra(2, iBas), GbasBra(3, iBas)
    end do
!    NOstopNO

    ! Initialize FFT with the fftw libary from reciprocal space to real space for ket- and bra-wavefunction coefficients
    allocate (zFFTkBox(0: ifftds - 1))
    zFFTkBox = (0.0,0.0)
    call dfftw_plan_dft_3d(backwardPlanK, ifft1ds, ifft2ds, ifft3ds, zFFTkBox, zFFTkBox, FFTW_BACKWARD, FFTW_MEASURE)

    allocate (zFFTbBox(0:ifftds-1))
    zFFTbBox = (0.0,0.0)
    call dfftw_plan_dft_3d(backwardPlanB, ifft1ds, ifft2ds, ifft3ds, zFFTbBox, zFFTbBox, FFTW_BACKWARD, FFTW_MEASURE)

    ! Initialize FFT with the fftw libary from real space to reciprocal space for every band combination possible within the Sternheimer
    ! equation
    allocate ( forwardPlanR(nrBraBands, nrKetBands) )
    allocate ( zVzFFTBox(0 : ifftds - 1, nrBraBands, nrKetBands) )
    zVzFFTBox = (0.0,0.0)
    do iKband = 1, nrKetBands
      do iBband = 1, nrBraBands
        call dfftw_plan_dft_3d(forwardPlanR(iBband, iKband), ifft1ds, ifft2ds, ifft3ds, zVzFFTBox(:, iBband, iKband), zVzFFTBox(:, iBband, iKband), FFTW_FORWARD, FFTW_MEASURE)
      end do
    end do

    zVzFFTBox = (0.0,0.0)
    do iKband = 1, nrKetBands
      zFFTkBox = (0.0,0.0)
      ! Put ket-wavefunction expansion coefficients onto main FFT mesh and perform FFT with the fftw library
      do iBas = 1, nv(1, ikpt)
        zFFTkBox(iv1dKet(iBas,1)) = zKet(iBas,iKband)
      end do
      call dfftw_execute_dft(backwardPlanK, zFFTkBox, zFFTkBox)

      do iBband = 1, nrBraBands

        zFFTbBox = (0.0,0.0)
        ! Put bra-wavefunction expansion coefficients onto main FFT mesh and perform FFT with the fftw library
        do iBas = 1, nv(1, ikpq)
          ! is conjugated below
          !zFFTbBox(iv1dBra(iBas,1)) = conjg(zBra(iBas,iBband))
          zFFTbBox(iv1dBra(iBas,1)) = zBra(iBas,iBband)
        end do
        CALL dfftw_execute_dft(backwardPlanB, zFFTbBox, zFFTbBox)

        ! In real space the convolution of the reciprocal space is just the product of all quantities
        do imesh = 0, ifftds-1
          ! Bra was conjugated before FFT
          zVzFFTBox(imesh, iBband, iKband ) = conjg(zFFTbBox(imesh)) * VFFTBox(imesh) * zFFTkBox(imesh)
        end do

        ! Perform a FFT from real space to reciprocal space of the convoluted quantity and delete the configuration data of the fftw
        ! library for this kind of transformation
        call dfftw_execute_dft(forwardPlanR(iBband, iKband), zVzFFTBox(:, iBband, iKband), zVzFFTBox(:, iBband, iKband))
        call dfftw_destroy_plan(forwardPlanR(iBband, iKband))
      end do
    end do
    ! Correct quantity with factor obtained during the FFT
    zVzFFTBox(:, :, : ) = zVzFFTBox(:, :, :) / real(ifftds) !todo also check this in other calculation!

    ! Destroy remaining FFT configuration data to clear up memory
    call dfftw_destroy_plan(backwardPlanK)
    call dfftw_destroy_plan(backwardPlanB)

    ! As the following calculation is independent of the bands we optimize runtime for sake of storage by running over the FFT grid
    ! indices which also match to the grid where the convoluted z V z in reciprocal space is given for every G. We can determine the
    ! G-vector for every grid index and assume that the convoluted quantity is 1 for the moment.
    ! the just 

    prf = fpi * iu / cell%omtil
    raf = atoms%rmt(iDtype)**2

    allocate(preRaileigh(0:ifftds -1))
    preRaileigh = cmplx(0.0, 0.0)
    do xIndex = -sk1d, sk1d
      xis = xIndex
      if (xIndex < 0) then
        xis = xis + ifft1ds
      end if
      do yIndex = -sk2d, sk2d
        yis = yIndex
        if (yIndex < 0) then
          yis = yis + ifft2ds
        end if
        do zIndex = -sk3d, sk3d
          zis = zIndex
          if (zIndex < 0) then
            zis = zis + ifft3ds
          end if
          gridIndex = (xis+ifft1ds*yis+ifft1ds*ifft2ds*zis)
          Gpqvec(:) = [xIndex - qpts%bk(1, iqpt), yIndex - qpts%bk(2, iqpt), zIndex - qpts%bk(3, iqpt)]
          Greal(:) = matmul(cell%bmat(:, :), Gpqvec(:))
          call ylm4(1, Greal, ylm)
          fj = 0.0
          call sphbes(1, norm2(Greal) * atoms%rmt(iDtype), fj)
          sumt = cmplx(0.0,0.0)
          do lm = 2, 4
            ! We do not have the conjugated natural coordinate expansion coefficients because only the Y_lm were conjugated
            ! by multiplying them with a factor (-1)^m canceling away with the same factor stemming from ylm not being
            ! conjugated, although it is claimed in the Rayleigh expansion
            sumt = sumt + c_im(iDdir, lm - 1) * ylm(lm)
          end do ! t
          preRaileigh(gridIndex) = preRaileigh(gridIndex) +  prf * raf * fj(1) * sumt * exp( iu * tpi * dot_product(real(Gpqvec(:)), atoms%taual(:, iDatom)))
         end do
      end do
    end do

    do iKband = 1, nrKetBands
      do iBband = 1, nrBraBands
        surfInt(iBband, iKband) = surfInt(iBband, iKband) + dot_product(conjg(preRaileigh(:)), zVZFFTBox(:, iBband, iKband)) 
      end do
    end do

  end subroutine calcSfVeff

  ! NOTE: The Gbas vectors have to be already sorted by the array ilst.
  !---------------------------------------------------------------------------------------------------------------------------------
  !> @author
  !> Christian-Roman Gerhorst, Gregor Michaliczek
  !>
  !> @brief
  !> Calls routines to set-up Sternheimer equation and solves it for first-order wavefunction expansion coefficients.
  !>
  !> @details
  !>
  !> @note
  !> We have to consider all (occupied and unoccupied) bands p in the bra and only the occupied bands n in the kets.
  !> Using the OEP approach of Markus Betzinger, it should be possible to also only consider the occupied bands in the bras leading to
  !> significant runtime enhancements.
  !---------------------------------------------------------------------------------------------------------------------------------
  subroutine PrepareMTSurfInt( atoms, dimens, kpts, lathar, V0Fleur, iDtype, iDatom, nRadFun, El, rbas1, rbas2, nmem_atom, mlh_atom,&
      & clnu_atom, hFullNoAbcofBK, overlapNoAbcofBK, z)

    use m_types
    use m_jpConstants, only : c_im
    use m_gaunt, only : gaunt1

    implicit none

    ! Type parameters
    type(t_atoms),                  intent(in)  :: atoms
    type(t_dimension),              intent(in)  :: dimens
    type(t_kpts),                   intent(in)  :: kpts
    type(t_sphhar),                 intent(in)  :: lathar
    type(t_potential),              intent(in)  :: V0Fleur

    ! Scalar parameters
    integer,                        intent(in)  :: iDtype
    integer,                        intent(in)  :: iDatom

    ! Array parameters
    integer,                        intent(in)  :: nRadFun(0:, :)
    real,                           intent(in)  :: El(:, 0:, :, :)
    real,                           intent(in)  :: rbas1(:, :, 0:, :)
    real,                           intent(in)  :: rbas2(:, :, 0:, :)
    integer,                        intent(in)  :: nmem_atom(0:, :)
    integer,                        intent(in)  :: mlh_atom(:,0:,:)
    complex,                        intent(in)  :: clnu_atom(:,0:,:)
    MCOMPLEX,                       intent(in)  :: z(:,:,:,:)
    complex,                        intent(out) :: hFullNoAbcofBK(:, :, :, :)
    complex,                        intent(out) :: overlapNoAbcofBK(:, :, :, :)

    ! Scalar variables
    integer                                     :: lmpMax
    integer                                     :: oqn_l
    integer                                     :: rMt
    integer                                     :: lmp
    integer                                     :: mqn_m
    integer                                     :: iradf
    integer                                     :: ptsym
    integer                                     :: ilh
    integer                                     :: oqn_l3p
    integer                                     :: oqn_l2p
    integer                                     :: mqn_m4p
    integer                                     :: mqn_m2p
    integer                                     :: oqn_l1p
    integer                                     :: mqn_m1p
    integer                                     :: mqn_m3p
    integer                                     :: lm
    integer                                     :: lmp1p
    integer                                     :: lm1p
    integer                                     :: imem
    real                                        :: gauntVeffNv
    integer                                     :: idir
    real                                        :: gauntComplete
    complex                                     :: hsphSymLmp1pLmp
    complex                                     :: veffnSphLmp1pLmp
    real                                        :: rbasOverlap
    integer                                     :: iradf1p

    ! Array variables
    complex,           allocatable              :: zDummy(:, :, :, :)
    complex,           allocatable              :: vEfflmLm1P(:, :, :)
    complex,           allocatable              :: vEfflmpLmp1P(:, :, :, :)
    complex,           allocatable              :: gauntVec(:, :, :, :)
    real,              allocatable              :: hSphUlmRmt(:, :, :)

    !todo write(*, *) 'remove zDummy'
    allocate( zDummy(dimens%nbasfcn, dimens%neigd, kpts%nkpt, 1) )
    zDummy = cmplx(0., 0.)
    zDummy = z

    ! Spherical part only contributes until lmax and not lmax + 2 as we have only the ket here expanded until lmax and in the overlap we
    ! cut there at lmax.
    lmpMax = maxval( (/ (sum( (/ ((2*oqn_l+1)* nRadFun(oqn_l,iDtype), oqn_l = 0,atoms%lmax(iDtype)) /) ),iDtype=1,atoms%ntype) /) )
    allocate( hSphUlmRmt(2, lmpMax, atoms%ntype) )
    hSphUlmRmt(:, :, :) = 0.

    allocate( vEfflmLm1P((atoms%lmaxd + 1)**2, (atoms%lmaxd + 1)**2, 3) )
    allocate( vEfflmpLmp1P(lmpMax, lmpMax, 3, atoms%nat) )
    allocate( gauntVec( (atoms%lmaxd + 1)**2, (atoms%lmaxd + 1)**2, 3, atoms%ntype ) )
    gauntVec(:, :, :, :) = cmplx(0., 0.)

    rMt = atoms%jri(iDtype)
    lmp = 0
    ! Calculate H_sph u_lmp(R_Mt) * R_mT (rbas arrays are already multiplied with the mesh once)
    do oqn_l = 0, atoms%lmax(iDtype)
      do mqn_m = -oqn_l, oqn_l
        lmp = lmp + 1
        hSphUlmRmt(1, lmp, iDtype) = El(1, oqn_l, iDtype, 1) * rbas1(rMt, 1, oqn_l, iDtype)
        hSphUlmRmt(2, lmp, iDtype) = El(1, oqn_l, iDtype, 1) * rbas2(rMt, 1, oqn_l, iDtype)

        lmp = lmp + 1
        hSphUlmRmt(1, lmp, iDtype) = rbas1(rMt, 1, oqn_l, iDtype) + El(1, oqn_l, iDtype, 1) * rbas1(rMt, 2, oqn_l, iDtype)
        hSphUlmRmt(2, lmp, iDtype) = rbas2(rMt, 1, oqn_l, iDtype) + El(1, oqn_l, iDtype, 1) * rbas2(rMt, 2, oqn_l, iDtype)

        do iradf = 3, nRadFun(oqn_l, iDtype)
          lmp = lmp + 1
          hSphUlmRmt(1, lmp, iDtype) = El(iradf - 1, oqn_l, iDtype, 1) * rbas1(rMt, iradf, oqn_l, iDtype)
          hSphUlmRmt(2, lmp, iDtype) = El(iradf - 1, oqn_l, iDtype, 1) * rbas2(rMt, iradf, oqn_l, iDtype)
        end do ! p
      end do ! mqn_m
    end do ! oqn_l

    ptsym = atoms%ntypsy(iDatom)
    vEfflmLm1P = cmplx(0., 0.)
    ! Do all the triangular rules math for the effective potential and calculate the product of the the effective potential with
    ! the normal vector and two Gaunt coefficients
    ! We use only the non-spherical part of the potential
    do ilh = 1, lathar%nlh(ptsym)
      oqn_l3p = lathar%llh(ilh, ptsym) ! this probably is not required
      if (oqn_l3p == 0) cycle
      ! As we only treat the non-spherical part oqn_l2p never is smaller than 0
      do oqn_l2p = abs(oqn_l3p - 1), min(oqn_l3p + 1, atoms%lmax(iDtype)), 2
        do mqn_m4p = -1, 1
          do imem = 1, nmem_atom(ilh, iDatom)
            mqn_m3p = mlh_atom(imem, ilh, iDatom)
            mqn_m2p = mqn_m3p + mqn_m4p ! tested
            if (abs(mqn_m2p) > oqn_l2p) cycle
            gauntVeffNv = gaunt1( oqn_l2p, oqn_l3p, 1, mqn_m2p, mqn_m3p, mqn_m4p, atoms%lmaxd )
            do idir = 1, 3
              do oqn_l = 0, atoms%lmax(iDtype)
                do oqn_l1p = abs(oqn_l2p - oqn_l), min(oqn_l2p + oqn_l, atoms%lmax(iDtype))
                  do mqn_m = -oqn_l, oqn_l
                    mqn_m1p = mqn_m3p + mqn_m4p + mqn_m
                    if ( abs(mqn_m1p) > oqn_l1p ) cycle
                    lm = oqn_l * (oqn_l + 1) + 1 + mqn_m
                    lm1p = oqn_l1p * (oqn_l1p + 1) + 1 + mqn_m1p
                    gauntComplete = gaunt1( oqn_l1p, oqn_l2p, oqn_l, mqn_m1p, mqn_m2p, mqn_m, atoms%lmaxd )
                ! Constant potential for testing
                !    vEfflmLm1P(lm1p, lm, idir) = vEfflmLm1P(lm1p, lm, idir) + c_im(idir, mqn_m4p + 2) * vEffDummy(rMt, ilh + 1, idir, iDatom)&
                !      & * gauntVeffNv * gauntComplete
                    vEfflmLm1P(lm1p, lm, idir) = vEfflmLm1P(lm1p, lm, idir) + c_im(idir, mqn_m4p + 2) * V0Fleur%vr(rMt, ilh, iDtype, 1)                      &
                      & * clnu_atom(imem, ilh, iDatom) * gauntVeffNv * gauntComplete
                  end do ! mqn_m
                end do ! oqn_l1p
              end do ! oqn_l
            end do ! idir
          end do ! imem
        end do ! mqn_m4p
      end do ! oqn_l2p
    end do ! ilh

    ! Sum up the Gaunt coefficients of the spherical Hamiltonian
    do idir = 1, 3
      do oqn_l = 0, atoms%lmax(iDtype)
        do mqn_m = -oqn_l, oqn_l
          lm = oqn_l * (oqn_l + 1) + 1 + mqn_m
          do oqn_l1p = abs(oqn_l - 1), min(oqn_l + 1, atoms%lmax(iDtype)), 2
            do mqn_m2p = -1, 1
              mqn_m1p = mqn_m + mqn_m2p
              if ( abs(mqn_m1p) > oqn_l1p ) cycle
              lm1p = oqn_l1p * (oqn_l1p + 1) + 1 + mqn_m1p
              gauntComplete = gaunt1( oqn_l1p, oqn_l, 1, mqn_m1p, mqn_m, mqn_m2p, atoms%lmaxd )
              gauntVec(lm1p, lm, idir, iDtype) = gauntVec(lm1p, lm, idir, iDtype) + c_im(idir, mqn_m2p + 2) * gauntComplete
            end do ! mqn_m1Pr
          end do ! oqn_l1Pr
        end do ! mqn_m
      end do ! oqn_l
    end do ! idir

    ! Calculate the matrixelement with the potential and the overlap potential matrix element so far that only the matching
    ! coefficients have to be multiplied by a matrix multiplication
    do idir = 1, 3
      lmp = 0
      do oqn_l = 0, atoms%lmax(iDtype)
        do mqn_m = -oqn_l, oqn_l
          lm = oqn_l * ( oqn_l + 1 ) + 1 + mqn_m
          do iradf = 1, nRadFun(oqn_l, iDtype)
            lmp = lmp + 1
            lmp1p = 0
            do oqn_l1p = 0, atoms%lmax(iDtype)
              do mqn_m1p = -oqn_l1p, oqn_l1p
                lm1p = oqn_l1p * ( oqn_l1p + 1 ) + 1 + mqn_m1p
                do iradf1p = 1, nRadFun(oqn_l1p, iDtype)
                  lmp1p = lmp1p + 1
!                  hsphSymLmp1pLmp = 0.5 * gauntVec(lm1p, lm, idir, iDtype) * &
                  hsphSymLmp1pLmp = gauntVec(lm1p, lm, idir, iDtype) * &
                    & ( rbas1(rMT, iradf1p, oqn_l1p, iDtype) * hSphUlmRmt(1, lmp, iDtype) &
                    & + rbas2(rMT, iradf1p, oqn_l1p, iDtype) * hSphUlmRmt(2, lmp, iDtype) )!&
!                    & + hSphUlmRmt(1, lmp1p, iDtype)         * rbas1(rMT, iradf, oqn_l, iDtype) &
!                    & + hSphUlmRmt(2, lmp1p, iDtype)         * rbas2(rMT, iradf, oqn_l, iDtype) )
                  rbasOverlap = rbas1(rMt, iradf1p, oqn_l1p, iDtype) * rbas1(rMt, iradf, oqn_l, iDtype) &
                            & + rbas2(rMt, iradf1p, oqn_l1p, iDtype) * rbas2(rMt, iradf, oqn_l, iDtype)
                  veffnSphLmp1pLmp = rbasOverlap * vEfflmLm1P(lm1p, lm, idir)
                  hFullNoAbcofBK(lmp1p, lmp, idir, iDatom) = hsphSymLmp1pLmp + veffnSphLmp1pLmp
                  overlapNoAbcofBK(lmp1p, lmp, idir, iDtype) = rbasOverlap * gauntVec(lm1p, lm, idir, iDtype)
                end do ! iradf1p
              end do ! mqn_m1p
            end do ! oqn_l1p
          end do ! iradf
        end do ! mqn_m
      end do ! oqn_l
    end do ! idir

  end subroutine PrepareMTSurfInt

  ! NOTE: The Gbas vectors have to be already sorted by the array ilst.
  !---------------------------------------------------------------------------------------------------------------------------------
  !> @author
  !> Christian-Roman Gerhorst, Gregor Michaliczek
  !>
  !> @brief
  !> Calls routines to set-up Sternheimer equation and solves it for first-order wavefunction expansion coefficients.
  !>
  !> @details
  !>
  !> @note
  !> We have to consider all (occupied and unoccupied) bands p in the bra and only the occupied bands n in the kets.
  !> Using the OEP approach of Markus Betzinger, it should be possible to also only consider the occupied bands in the bras leading to
  !> significant runtime enhancements.
  !---------------------------------------------------------------------------------------------------------------------------------
  subroutine CalcSurfIntMT(atoms, dimens, sym, cell, kpts, input, ud, ikpq, ikpt, iDdir, iDatom, iDtype, nv, gbas, ilst, zkpq, zkpt, nobd, &
      & ne, kveclo, nRadFun, iloTable, hFullNoAbcofBK, overlapNoAbcofBK, overlap, surfIntMT)

    use m_types, only : t_atoms, t_kpts, t_sym, t_dimension, t_cell, t_input, t_usdus, t_noco
     
    use m_jpConstants, only : iu
    use m_abcof

    implicit none

    ! Type parameters
    type(t_atoms),                  intent(in)  :: atoms
    type(t_dimension),              intent(in)  :: dimens
    type(t_sym),                    intent(in)  :: sym
    type(t_cell),                   intent(in)  :: cell
    type(t_kpts),                   intent(in)  :: kpts
    type(t_input),                  intent(in)  :: input
    type(t_usdus),                  intent(in)  :: ud

    integer,                        intent(in)  :: ikpt
    integer,                        intent(in)  :: ikpq
    integer,                        intent(in)  :: iDdir
    integer,                        intent(in)  :: iDatom
    integer,                        intent(in)  :: iDtype

    integer,                        intent(in)  :: nv(:, :)
    integer,                        intent(in)  :: gbas(:, :)
    integer,                        intent(in)  :: ilst(:, :, :)
    MCOMPLEX,                       intent(in)  :: zkpt(:,:)
    MCOMPLEX,                       intent(in)  :: zkpq(:,:)
    integer,                        intent(in)  :: nobd(:, :)
    integer,                        intent(in)  :: kveclo(:,:)
    integer,                        intent(in)  :: ne(:)
    integer,                        intent(in)  :: iloTable(:, 0:, :)
    integer,                        intent(in)  :: nRadFun(0:, :)
    complex,                        intent(in)  :: hFullNoAbcofBK(:, :, :, :)
    complex,                        intent(in)  :: overlapNoAbcofBK(:, :, :, :)
    complex,                        intent(out) :: overlap(:, :)
    complex,                        intent(out) :: surfIntMT(:, :)

    ! Type variables
    type(od_inp)                                :: odi
    type(od_sym)                                :: ods
    type(t_noco)                                :: noco

    ! Scalar variables
    integer                                     :: nmat
    integer                                     :: lmp
    integer                                     :: ilo
    integer                                     :: oqn_l
    integer                                     :: iradf
    integer                                     :: ibandK
    integer                                     :: ibandB
    integer                                     :: mqn_m
    integer                                     :: lm
    integer                                     :: lmpMax

    ! Array variables
    integer,           allocatable              :: ngoprI(:)
    complex,           allocatable              :: mCofK(:, :)
    complex,           allocatable              :: mCofB(:, :)
    complex,           allocatable              :: hFullNoAbcofK(:, :)
    complex,           allocatable              :: overlapNoAbcofK(:, :)
    complex,           allocatable              :: hFull(:, :)
    complex,           allocatable              :: acof(:, :, :)
    complex,           allocatable              :: bcof(:, :, :)
    complex,           allocatable              :: ccof(:, :, :, :)
    complex,           allocatable              :: acofkpq(:, :, :)
    complex,           allocatable              :: bcofkpq(:, :, :)
    complex,           allocatable              :: ccofkpq(:, :, :, :)

    ! Spherical part only contributes until lmax and not lmax + 2 as we have only the ket here expanded until lmax and in the overlap we
    ! cut there at lmax.
    lmpMax = maxval( (/ (sum( (/ ((2*oqn_l+1)* nRadFun(oqn_l,iDtype), oqn_l = 0,atoms%lmax(iDtype)) /) ),iDtype=1,atoms%ntype) /) )

    allocate( mCofK(lmpMax, maxval(nobd(:, :))) )
    allocate( mCofB(dimens%neigd, lmpMax ) )
    allocate( hFullNoAbcofK(dimens%neigd, lmpMax) )
    allocate( hFull(dimens%neigd, maxval(nobd(:, :))) )
    allocate( overlapNoAbcofK( dimens%neigd, lmpMax) )

    ! Check that abcof
    allocate(acof(dimens%neigd, 0:dimens%lmd, atoms%nat), bcof(dimens%neigd, 0:dimens%lmd, atoms%nat), &
      &ccof(-atoms%llod:atoms%llod, dimens%neigd, atoms%nlod, atoms%nat))
    allocate(acofkpq(dimens%neigd, 0:dimens%lmd, atoms%nat), bcofkpq(dimens%neigd, 0:dimens%lmd, atoms%nat), &
      &ccofkpq(-atoms%llod:atoms%llod, dimens%neigd, atoms%nlod, atoms%nat))
    allocate( noco%alph(atoms%ntype), noco%beta(atoms%ntype) ) !Up to now those variables are only of dummy character
    allocate(ngoprI(atoms%nat))
    ngoprI(:) = 1


    surfIntMT = cmplx(0., 0.)
    hFullNoAbcofK = cmplx(0., 0.)
    hFull = cmplx(0., 0.)
    overlapNoAbcofK = cmplx(0., 0.)
    overlap = cmplx(0., 0.)
    acof(:, :, :) = cmplx(0., 0.)
    bcof(:, :, :) = cmplx(0., 0.)
    ccof(:, :, :, :) = cmplx(0., 0.)
    nmat = nv(1, ikpt) + atoms%nlotot
    call abcof ( atoms%lmaxd, atoms%ntype, dimens%neigd, dimens%neigd, atoms%nat, sym%nop, dimens%nvd, dimens%jspd, dimens%lmd, &
      & dimens%nbasfcn, atoms%llod, atoms%nlod, atoms%nlotot, sym%invtab, atoms%ntype, sym%mrot, ngoprI, atoms%taual, atoms%neq, &
      & atoms%lmax, atoms%rmt, cell%omtil, cell%bmat, cell%bbmat, kpts%bk(:, ikpt), &
      & gbas(1, ilst(:nv(input%jspins, ikpt), ikpt, input%jspins)), gbas(2, ilst(:nv(input%jspins, ikpt), ikpt, input%jspins)), &
      & gbas(3, ilst(:nv(input%jspins, ikpt), ikpt, input%jspins)), nv(:, ikpt),  nmat, ne(ikpt), zkpt(:, :), ud%us(:, :, 1), &
      & ud%dus(:, :, 1), ud%uds, ud%duds(:, :, 1), ud%ddn(:, :, 1), atoms%invsat, sym%invsatnr, ud%ulos(:, :, 1), &
      & ud%uulon(:, :, 1), ud%dulon(:, :, 1),  ud%dulos(:, :, 1), atoms%llo, atoms%nlo, atoms%l_dulo, &
      & atoms%lapw_l, noco%l_noco, noco%l_ss, 1, noco%alph, noco%beta, noco%qss, kveclo(:, ikpt), odi, ods, &
      & acof, bcof, ccof)

    acofkpq(:, :, :) = cmplx(0., 0.)
    bcofkpq(:, :, :) = cmplx(0., 0.)
    ccofkpq(:, :, :, :) = cmplx(0., 0.)
    nmat = nv(1, ikpq) + atoms%nlotot
    call abcof ( atoms%lmaxd, atoms%ntype, dimens%neigd, dimens%neigd, atoms%nat, sym%nop, dimens%nvd, dimens%jspd, dimens%lmd, &
      & dimens%nbasfcn, atoms%llod, atoms%nlod, atoms%nlotot, sym%invtab, atoms%ntype, sym%mrot, ngoprI, atoms%taual, atoms%neq, &
      & atoms%lmax, atoms%rmt, cell%omtil, cell%bmat, cell%bbmat, kpts%bk(:, ikpq), &
      & gbas(1, ilst(:nv(input%jspins, ikpq), ikpq, input%jspins)), gbas(2, ilst(:nv(input%jspins, ikpq), ikpq, input%jspins)), &
      & gbas(3, ilst(:nv(input%jspins, ikpq), ikpq, input%jspins)), nv(:, ikpq),  nmat, ne(ikpq), zkpq(:, :), ud%us(:, :, 1), &
      & ud%dus(:, :, 1), ud%uds, ud%duds(:, :, 1), ud%ddn(:, :, 1), atoms%invsat, sym%invsatnr, ud%ulos(:, :, 1), &
      & ud%uulon(:, :, 1), ud%dulon(:, :, 1),  ud%dulos(:, :, 1), atoms%llo, atoms%nlo, atoms%l_dulo, &
      & atoms%lapw_l, noco%l_noco, noco%l_ss, 1, noco%alph, noco%beta, noco%qss, kveclo(:, ikpq), odi, ods, &
      & acofkpq, bcofkpq, ccofkpq)

    ! Reorganize matching coefficients to include the LO contributions more naturally
    mCofB(:, :) = cmplx(0., 0.)
    mCofK(:, :) = cmplx(0., 0.)

    lmp   =   0
    lm    =  -1
    do oqn_l = 0, atoms%lmax(iDtype)
      do mqn_m = -oqn_l, oqn_l
        lm = lm + 1
        !p = 1
        lmp = lmp + 1
        do ibandB = 1, ne(ikpq)
          mCofB(ibandB, lmp) = acofkpq(ibandB, lm, iDatom) * iu**oqn_l
        end do ! ibandB
        !p = 2
        lmp = lmp + 1
        do ibandB = 1, ne(ikpq)
          mCofB(ibandB, lmp) = bcofkpq(ibandB, lm, iDatom)* iu**oqn_l
        end do ! ibandB
        !LOs
        do iradf = 3, nRadFun(oqn_l, iDtype)
          ilo = iloTable(iradf, oqn_l, iDtype)
          lmp = lmp + 1
          do ibandB = 1, ne(ikpq)
            mCofB(ibandB, lmp) = ccofkpq(mqn_m, ibandB, ilo, iDatom)* iu**oqn_l
          end do ! ibandB
        end do ! iradf
      end do ! mqn_m
    end do ! oqn_l

    do ibandK = 1, nobd(ikpt, 1)
      lmp   =   0
      lm    =  -1
      do oqn_l = 0, atoms%lmax(iDtype)
        do mqn_m = -oqn_l, oqn_l
          lm = lm + 1
          !p = 1
          lmp = lmp + 1
          mCofK(lmp, ibandK) = acof(ibandK, lm, iDatom) * iu**oqn_l

          !p = 2
          lmp = lmp + 1
          mCofK(lmp, ibandK) = bcof(ibandK, lm, iDatom)* iu**oqn_l

          !LOs
          do iradf = 3, nRadFun(oqn_l, iDtype)
            ilo = iloTable(iradf, oqn_l, iDtype)
            lmp = lmp + 1
            mCofK(lmp, ibandK) = ccof(mqn_m, ibandK, ilo, iDatom)* iu**oqn_l
          end do ! iradf
        end do ! mqn_m
      end do ! oqn_l
    end do ! ibandK

    ! Just from here on iDdir is relevant!
    ! delete surfIntMT for every ikpt
    hFullNoAbcofK(1:ne(ikpq), 1:lmpMax) = &
                           & matmul( conjg(mCofB(1:ne(ikpq), 1:lmpMax)), hFullNoAbcofBK(1:lmpMax, 1:lmpMax, iDdir, iDatom) )
    surfIntMT(1:ne(ikpq), 1:nobd(ikpt, 1)) = &
                           & matmul( hFullNoAbcofK(1:ne(ikpq), 1:lmpMax), mCofK(1:lmpMax, 1:nobd(ikpt, 1)) )

    overlapNoAbcofK(1:ne(ikpq), 1:lmpMax) = &
                           & matmul( conjg(mCofB(1:ne(ikpq), 1:lmpMax)), overlapNoAbcofBK(1:lmpMax, 1:lmpmax, iDdir, iDtype) )
    ! has to be iDatom because of mCofs(iDatom)
    overlap(1:ne(ikpq), 1:nobd(ikpt, 1)) = &
                           & matmul( overlapNoAbcofK(1:ne(ikpq), 1:lmpMax), mCofK(1:lmpMax, 1:nobd(ikpt, 1)) )


!    if (.false.) then
!      do ibandK = 1, nobd(ikpt, 1)
!        do ibandB = 1, ne(ikpt)
!          if ( ( abs(real(surfIntMT(ibandB, ibandK, iDdir, iDatom))) < 1e-8 ) .and. ( abs(aimag(surfIntMT(ibandB, ibandK, iDdir, iDatom))) >= 1e-8 )) then
!            write(2033, '(4(i8),2(f15.8))') ikpt, iDdir, ibandK, ibandB, 0., -aimag(surfIntMT(ibandB, ibandK, iDdir, iDatom))
!          else if ( ( abs( real(surfIntMT(ibandB, ibandK, iDdir, iDatom))) >= 1e-8 ) .and. ( abs(aimag(surfIntMT(ibandB, ibandK, iDdir, iDatom))) < 1e-8 )) then
!            write(2033, '(4(i8),2(f15.8))') ikpt, iDdir, ibandK, ibandB, -real(surfIntMT(ibandB, ibandK, iDdir, iDatom)), 0.
!          else if ( ( abs(real(surfIntMT(ibandB, ibandK, iDdir, iDatom))) < 1e-8 ) .and. ( abs(aimag(surfIntMT(ibandB, ibandK, iDdir, iDatom))) < 1e-8 )) then
!            write(2033, '(4(i8),2(f15.8))') ikpt, iDdir, ibandK, ibandB, 0., 0.
!          else if ( ( abs(real(surfIntMT(ibandB, ibandK, iDdir, iDatom))) >= 1e-8 ) .and. ( abs(aimag(surfIntMT(ibandB, ibandK, iDdir, iDatom))) >= 1e-8 )) then
!            write(2033, '(4(i8),2(f15.8))') ikpt, iDdir, ibandK, ibandB, -real(surfIntMT(ibandB, ibandK, iDdir, iDatom)), -aimag(surfIntMT(ibandB, ibandK, iDdir, iDatom))
!          end if
!        end do ! ibandB
!      end do ! ibandK
!    end if

    ! Note: If we only switch on the potential, this integral compares well to the interstitial version. The same is true for a
    !       constant potential. Comparing directly the wavefunction from the interstitial and the muffin-tin version, the same comes
    !       out. But if we compare the full action of the spherical Hamiltonian, there is a difference in the values which is due to
    !       the fact, that we compare the second derivative of the LAPW basis function which is not continious any more, necessarily
    !       Therefore, we belive this discrepance to come from this fact.
  end subroutine CalcSurfIntMT

  !todo code snippet of #123, do not delete
  ! Maintenance
 ! if (.false.) then
 !   allocate( rho1IRDStest(ngpqdp, 3, atoms%nat) )
 !   allocate( rho1MTTest(atoms%jmtd, (atoms%lmaxd + 1)**2, atoms%nat, 3, atoms%nat))
 !   rho1IRDStest(:, :, :) = cmplx(0., 0.)
 !   rho1MTTest(:, :, :, :, :) = cmplx(0., 0.)
 !   ! Test whether the loading and the storage of the density works well.
 !   write(filename, '(a10,i1,a4,i1,a4,i1)') 'JPcdn1_Dat', iDatom, 'Ddir', idir, 'qInd', iqpt
 !   call loadDensity(atoms, ngdp, filename, rho1IRDSTest(:, idir, iDatom), rho1MTTest(:, :, :, idir, iDatom))
 !   if (any(abs(rho1IRDSTest(:, idir ,iDatom) - rho1IRDS(:,idir, iDatom)) > 1e-9)) then
 !     write(*, *) 'IR inconsistent for dir', idir
 !   end if
 !   if (any(abs(rho1MTTest(:, :, :, idir ,iDatom) - rho1MT(:, :, :, idir, iDatom)) > 1e-9)) then
 !     write(*, *) 'MT inconsistent for dir', idir
 !   end if
 ! end if

end module m_jpTestSternheimer
