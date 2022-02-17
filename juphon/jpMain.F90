!>----------------------------------------------------------------------------------------------------------------------------------
!> IAS 1 / PGI 1 Forschungszentrum Jülich, juPhon Plugin for the FLEUR program
!>----------------------------------------------------------------------------------------------------------------------------------
!>
!> @author
!> Christian-Roman Gerhorst
!>
!> @brief
!> Main program operating juPhon.
!>
!> @details
!> This program calls subroutines fulfilling the major tasks/milestones in aquiring the Dynamical Matrix within the DFPT method.
!>
!> @param atoms      : Atoms type, see types.f90.
!> @param cell       : Unit cell type, see types.f90.
!> @param sym        : Symmetries type, see types.f90.
!> @param stars      : Stars type, see types.f90.
!> @param lathar     : Lattice harmonics type, see types.f90.
!> @param input      : Input type, see types.f90.
!> @param dimens     : Dimension type, see types.f90.
!> @param enpara     : Energy parameter type, see types.f90.
!> @param xcpot      : Exchange potential type, see types.f90.
!> @param kpts       : K-points type, see types.f90.
!> @param qpts       : Q-points type, see types.f90.
!> @param Veff0      : Type containing unperturbed output potentials of Fleur, see types.f90.
!> @param usdus      : Type containing quantities consisting of the radial solutions, see types.f90.
!> @param results    : Results type, see types.f90.
!> @param logUnit    : Unit number for juPhon.log.
!> @param noPotsCon  : Number of points for which continuity test should be performed.
!> @param ngdp       : Number of G-vectors for potentials and densities.
!> @param ngdp2km    : Number of G-vectors for potentials and densities which are smaller than 2 kmax.
!> @param memd_atom  : Maximal number of members in all lattice harmonics for all atoms.
!> @param paPoX      : X-coordinate of path's direction for potential's gradient test.
!> @param paPoY      : Y-coordinate of path's direction for potential's gradient test.
!> @param paPoZ      : Z-coordinate of path's direction for potential's gradient test.
!> @param harSw      : Switch to activate Hartree potential.
!> @param extSw      : Switch to activate external potential.
!> @param xcSw       : Switch to activate exchange-correlation potential.
!> @param numaddQ    : Number of additional q-vectors
!> @param startTime  : Starting time of juPhon.
!> @param El         : Contains LAPW and LO energy parameters.
!> @param rbas1      : Large components of radial solution, its energy derivative and u_LO
!> @param rbas2      : Small components of radial solution, its energy derivative and u_LO
!> @param kveclo     : Basis G-vectors of local orbitals.
!> @param eig        : Contains Kohn\--Sham eigenvalues.
!> @param z          : Kohn\--Sham eigenvectors.
!> @param nRadFun    : Number of radial functions per orbital quantum number l and atom type.
!> @param GbasVec    : G-basis vectors
!> @param ilst       : For various k-points G-basis vectors occur more than once, thus they are only stored once in juPhon. This pointer
!>                     array contains the right index for GbasVec array to "unfold" G-basis vectors again.
!> @param iloTable   : Number of local orbital if orbital quantum number l, atom type and index p > 2 from nRadFun is given.
!> @param nv         : Number of LAPW G-basis vectors for given k-point.
!> @param ne         : Number of eigenvalues per k-point.
!> @param nobd       : Number of occupied bands per k-point and spin
!> @param rho0IR     : Star coefficients of the unperturbed and converged interstitial density parsed from Fleur.
!> @param rho0MT     : Radial coefficients of the unperturbed and converged muffin-tin densities parsed from Fleur.
!> @param gridf      : Initialized grid quantity to use intgrf routine.
!> @param gdp        : G-vectors of potentials and densities.
!> @param mlh_atom   : Magnetic quantum number m of lattice harmonic members for every atom.
!> @param nmem_atom  : Number of lattice harmonic members for every atom.
!> @param clnu_atom  : Phase mediating between stars and plane waves.
!> @param vrCoul     : Unperturbed radial coefficients of the converged Coulomb potential parsed from FLEUR.
!> @param vpwCoul_uw : Star coefficients of the unwarped (no convolution with step function) and converged interstitial Coulomb
!> @param mapKpq2K   : For a given k-point index and q-point index, this array gives the k-point set index of the result k + q
!> @param uuilon     : overlap integral between the radial functions of the integral (multiplied by ulo_der) of a local orbital and the
!>                     flapw radial function with the same l
!> @param duilon     : overlap integral between the radial functions of the integral of a local orbital and the energy derivative of the
!>                     flapw radial function with the same l
!> @param ulouilopn  : overlap integral between the radial functions of the integral of a local orbital and another local orbital with
!>                     the same l.
!> @param ilo2p      : mapping array giving the p value for given number of LO and itype
!> @param gdp2Ind    : Stores the index of a potential and density G-Vector. The dimensions are the G-vector components.
!> @param   gdp2iLim : Stores the min and maxvals of gdp2Ind
!> @param   vExt2IR  : Stores the interstitial coefficients of the 2nd order external potential.
!> @param E2ndOrdII  : Stores the interstitial coefficients of the 2nd order external potential.
!>
!> @attention        : Do not use q = 0 for optical phonons. Some terms for that are still missing.
!>---------------------------------------------------------------------------------------------------------------------------------------
program juPhon

#include "cppmacro.h"

  use m_jpInit,            only : InitializeJuPhon
  use m_types
  use m_jpSternheimer,     only : initSternheimerSCC
  use m_jpSternheimer,     only : solveSternheimerSCC
  use m_juDFT_time,        only : TimeStart, TimeNOstopNO, WriteTimes
  use m_jp2ndOrdQuant,     only : GenVext2, CalcIIEnerg2
  use m_jpSetupDynMat,     only : SetupDynamicMatrix
  use m_jpProcessDynMat,   only : DiagonalizeDynMat, CalculateFrequencies
  use m_jpPlotObservables, only : PlotGradV0, plotConvergedQuantities
  use mod_juPhonUtils,     only : Fclose
  use m_jpLog,             only : FinishLogFile
  use m_jpConstants,       only : compPhon

#ifdef DEBUG_MODE

  use m_jpTestInput,       only : DetailedInitializationTests
  use m_jpTestSternheimer, only : TestSternheimerSetup

  use jpTest1stVarDens,    only : test3DplotPotDens
  use m_jpTestDynMatHF, only : TestDynMatHF
  use m_jpTestDynMatPulay, only : TestDynMatPulay
  use m_jpTestDynMatSurf, only : TestDynMatSurf
  use m_jpTestGoldsteinCond, only : TestGoldsteinCond

  use m_jpTestConverged, only : TestConvergedQuantities

!  use m_jpTestDynMatDeprecated,  only : TestDynMatDeprecated

#endif

  implicit none

  ! Type variables
  type(t_atoms)                   :: atoms
  type(t_cell)                    :: cell
  type(t_sym)                     :: sym
  type(t_stars)                   :: stars
  type(t_sphhar)                  :: lathar
  type(t_input)                   :: input
  type(t_dimension)               :: dimens
  type(t_enpara)                  :: enpara
  type(t_xcpot)                   :: xcpot
  type(t_kpts)                    :: kpts
  type(t_kpts)                    :: qpts
  type(t_potential)               :: Veff0
  type(t_usdus)                   :: usdus
  type(t_results)                 :: results
  type(t_tlmplm)                  :: tdHS0
  type(t_oneD) :: oneD
  type(t_vacuum) :: vacuum

  ! Scalar variables
  integer                         :: logUnit = 100
  integer                         :: iqpt
  integer                         :: noPtsCon
  integer                         :: ngdp
  integer                         :: ngpqdp
  integer                         :: ngdp2km
  integer                         :: memd_atom
  real                            :: paPoX
  real                            :: paPoY
  real                            :: paPoZ
  logical                         :: harSw
  logical                         :: extSw
  logical                         :: xcSw
  logical                         :: testCompareGrVeff0FleurSw
  logical                         :: testVeff1Sw
  logical                         :: testUnfoldStarsSw
  logical                         :: testRadDerivativeSw
  logical                         :: testGauntCoeffSw
  logical                         :: testGradLhExpandFuncSw
  logical                         :: testContGrVeff0Sw
  logical                         :: testWarpingSw
  logical                         :: testSternhHSMEtestSw
  logical                         :: testSternhSchroedPertTheoSw
  logical                         :: testz1Phi0ContSw
  logical                         :: testRho1BasCorrSw
  logical                         :: testPlotRho03Dsw
  logical                         :: testRadSolSw
  logical                         :: testKptsWeightSw
  logical                         :: testCountValElecSw
  logical                         :: testVeff0ContSw
  logical                         :: testrho0ContSw
  logical                         :: testBackRotMTCoordSysSw
  logical                         :: testPsi0ContSw
  logical                         :: testOverlapSw
  logical                         :: testGradRho0PathSw
  logical                         :: testEii2LatPeriodQSw
  logical                         :: testVarphiHepsVarphiSw
  logical                         :: testRho1IRsw
  logical                         :: testRho1MTsw
  integer                         :: numAddQs
  logical                         :: calcEigenVec
  logical                         :: oneSternhCycle
  logical                         :: recipLengthUnit
  logical                         :: onlyTests
  logical                         :: testsActivated
  logical                         :: test1st2ndPulDynMatEps1
  logical                         :: test1st2ndPulDynMatCancel
  logical                         :: test3rdPulDynMatCancel
  logical                         :: testIntVeff1Rho1Val
  logical                         :: testGrHepsGrtMatElem
  logical                         :: testGrPsiPsiMatElem
  logical                         :: testCompareSurfInt
  logical                         :: testSplitMTSurfIntSterh
  logical                         :: testVeff1IRMESternh
  logical                         :: testEps1q0
  logical                         :: testVeff1IRMatqBackFold
  logical                         :: testVeff1IRqLatPeriod
  logical                         :: testGrMatElemPsiHepsPsiGaussTheo
  logical                         :: testPsiHepsTildePsi
  logical                         :: testGoldsteinRemaining
  logical                         :: testR2orNotWfMtGradNgrNTensGrOvls
  logical                         :: testComp3ArgSFIntsSw
  logical                         :: testComp2ArgSFIntsSw
  logical                         :: testGoldsteinSurfSw
  logical                         :: testComp2ArgGrSFIntsSw
  logical                         :: testIRIntegralSw
  logical                         :: testIR3rdMatElemSw
  logical                         :: testActionHgrPhiSw
  logical                         :: testXCintegrals
  logical                         :: testEii2PsDens
  logical                         :: plot3DgrVhar0
  logical                         :: plot3DgrVext0
  logical                         :: plot3DgrVxc0
  logical                         :: plot3DgrVeff0
  logical                         :: plotVhar1fccXZSw
  logical                         :: plotVext1fccXZSw
  logical                         :: plotVxc1fccXZSw
  logical                         :: plotVeff1fccXZSw
  logical                         :: plotRho1fccXZSw

  ! Array variables
  real,               allocatable :: El(:, :, :, :)
  real,               allocatable :: rbas1(:,:,:,:,:)
  real,               allocatable :: rbas2(:,:,:,:,:)
  integer,            allocatable :: kveclo(:,:)
  real,               allocatable :: eig(:,:,:)
  MCOMPLEX,           allocatable :: z(:,:,:,:)
  integer,            allocatable :: nRadFun(:, :)
  integer,            allocatable :: GbasVec(:, :)
  integer,            allocatable :: ilst(:, :, :)
  integer,            allocatable :: iloTable(:, :, :)
  integer,            allocatable :: nv(:, :)
  integer,            allocatable :: ne(:)
  integer,            allocatable :: nobd(:, :)
  complex,            allocatable :: rho0IR(:,:)
  real,               allocatable :: rho0MT(:,:,:,:)
  real,               allocatable :: gridf(:, :)
  integer,            allocatable :: gdp(:,:)
  integer,            allocatable :: mlh_atom(:,:,:)
  integer,            allocatable :: nmem_atom(:, :)
  complex,            allocatable :: clnu_atom(:,:,:)
  real,               allocatable :: vrCoul(:, :, :, :)
  complex,            allocatable :: vpwCoul_uw(:, :)
  integer,            allocatable :: mapKpq2K(:, :)
  real,               allocatable :: uuilon(:, :)
  real,               allocatable :: duilon(:, :)
  real,               allocatable :: ulouilopn(:, :, :)
  integer,            allocatable :: ilo2p(:, :)
  integer,            allocatable :: gdp2Ind(:, :, :)
  integer,            allocatable :: kpq2kPrVec(:, :, :)
  complex,            allocatable :: vExt2IR(:, :, :, :)
  complex,            allocatable :: vExt2MT(:, :, :, :)
  complex,            allocatable :: E2ndOrdII(:, :)
  complex,            allocatable :: dynMat(:, :)
  complex,            allocatable :: qpwcG(:, :)
  complex,            allocatable :: rho1MTCoreDispAt(:, :, :, :)
  complex,            allocatable :: grVxcIRKern(:)
  real,               allocatable :: dKernMTGPts(:, :, :)
  complex,            allocatable :: grVeff0MT_init(:, :, :, :)
  complex,            allocatable :: grVeff0MT_main(:, :, :, :)
  complex,            allocatable :: grVext0IR_DM(:, :)
  complex,            allocatable :: grVext0MT_DM(:, :, :, :)
  complex,            allocatable :: grVCoul0IR_DM_SF(:, :)
  complex,            allocatable :: grVCoul0MT_DM_SF(:, :, :, :)
  complex,            allocatable :: grVeff0IR_DM(:, :)
  complex,            allocatable :: grVeff0MT_DM(:, :, :, :)
  complex,            allocatable :: grVeff0MT_DMhxc(:, :, :, :)
  complex,            allocatable :: grRho0IR(:, :)
  complex,            allocatable :: grRho0MT(:, :, :, :)
  complex,            allocatable :: ylm(:, :)
  real,               allocatable :: gausWts(:) ! gaussian weights belonging to gausPts
  complex,            allocatable :: rho1IR(:, :, :)
  complex,            allocatable :: rho1MT(:, :, :, :, :)
  complex,            allocatable :: rho1MTz0(:, :, :, :)
  complex,            allocatable :: vExt1MT(:, :, :, :, :)
  complex,            allocatable :: vExt1MTnoVol(:, :, :, :, :)
  complex,            allocatable :: vEff1MTnoVol(:, :, :, :, :)
  complex,            allocatable :: vExt1MTnoVolnoq(:, :, :, :, :)
  complex,            allocatable :: vH1MTnoVol(:, :, :, :, :)
  complex,            allocatable :: vEff1IR(:, :, :)
  complex,            allocatable :: vEff1MT(:, :, :, :, :)
  real,               allocatable :: eigenVals(:)
  complex,            allocatable :: eigenFreqs(:)
  complex,            allocatable :: eigenVecs(:, :)
  integer,            allocatable :: gpqdp(:, :)
  complex,            allocatable :: vXC0IR(:, :)
  complex,            allocatable :: eXCIR(:)
  real,               allocatable :: vXC0MT(:, :, :, :)
  real,               allocatable :: eXCMT(:, :, :)
  complex,            allocatable :: vExt1IR_final(:, :, :)
  complex,            allocatable :: vExt1noqIR_final(:, :, :)
  complex,            allocatable :: vHar1IR_final(:, :, :)
  complex,            allocatable :: vHar1MT_final(:, :, :, :, :)
  complex,            allocatable :: rho1MTDelta(:, :, :, :, :)
  complex,            allocatable :: vExt1MTDelta(:, :, :, :, :)
  complex,            allocatable :: vExt1MTq0(:, :, :, :, :)
  complex,            allocatable :: vHar1MTDelta(:, :, :, :, :)
  complex,            allocatable :: vHar1MTq0(:, :, :, :, :)
  complex,            allocatable :: vXc1MTDelta(:, :, :, :, :)
  complex,            allocatable :: vXc1MTq0(:, :, :, :, :)
  complex,            allocatable :: rho0IRpw(:, :)
  complex,            allocatable :: vEff0IRpw(:, :)
  complex,            allocatable :: vEff0IRpwUw(:, :)
  complex,            allocatable :: rho0MTsh(:, :, :, :)
  complex,            allocatable :: vEff0MTsh(:, :, :, :)
  complex,            allocatable :: vCoul1IRtempNoVol(:, :)
  complex,            allocatable :: vCoul1MTtempNoVol(:, :, :, :)
  integer                         :: gdp2iLim(2, 3)
  integer                         :: startTime(8)



  call date_and_time( VALUES=startTime )

  ! Read in all quantities required for juPhon calculation from converged Fleur calculation or recalculate them, respectively.
  call TimeStart('Init JuPhon')
  call InitializeJuPhon( atoms, cell, sym, stars, lathar, input, dimens, enpara, xcpot, kpts, qpts, Veff0, usdus, results, logUnit,&
    & iqpt, paPoX, paPoY, paPoZ, harSw, extSw, xcSw, noPtsCon, ngdp, ngdp2km, memd_atom, numAddQs, startTime, kveclo, nv, nRadFun, &
    & iloTable, GbasVec, ilst, ne, mapKpq2K, eig, rbas1, rbas2, El, vrCoul, z, vpwCoul_uw, rho0IR, rho0MT, gridf, gdp, mlh_atom,   &
    & nmem_atom, clnu_atom, uuilon, duilon, ulouilopn, nobd, ilo2p, gdp2Ind, gdp2iLim, kpq2kPrVec, calcEigenVec, oneSternhCycle,   &
    & recipLengthUnit, onlyTests, testCompareGrVeff0FleurSw, testVeff1Sw, testUnfoldStarsSw, testRadDerivativeSw, testGauntCoeffSw,&
    & testGradLhExpandFuncSw, testContGrVeff0Sw, testWarpingSw, testSternhHSMEtestSw, testSternhSchroedPertTheoSw,                 &
    & testz1Phi0ContSw, testRho1BasCorrSw, testPlotRho03Dsw, testRadSolSw, testKptsWeightSw, testCountValElecSw, testVeff0ContSw,  &
    & testrho0ContSw, testBackRotMTCoordSysSw, testPsi0ContSw, testOverlapSw, testGradRho0PathSw, testEii2LatPeriodQSw,                     &
    & testVarphiHepsVarphiSw, testRho1IRsw, testRho1MTsw, testsActivated, oneD, vacuum,          &
    & test1st2ndPulDynMatEps1, test1st2ndPulDynMatCancel, test3rdPulDynMatCancel, testIntVeff1Rho1Val,                             &
    & testGrPsiPsiMatElem, testCompareSurfInt, testSplitMTSurfIntSterh, testVeff1IRMESternh,                   &
    & testEps1q0, testVeff1IRMatqBackFold, testVeff1IRqLatPeriod, testGrMatElemPsiHepsPsiGaussTheo, testPsiHepsTildePsi,           &
    & testGoldsteinRemaining, testR2orNotWfMtGradNgrNTensGrOvls, testComp3ArgSFIntsSw, testComp2ArgSFIntsSw, testGoldsteinSurfSw, testComp2ArgGrSFIntsSw, testIRIntegralSw, testIR3rdMatElemSw, testActionHgrPhiSw, testXCintegrals, testEii2PsDens, vXC0IR, eXCIR, vXC0MT, eXCMT, rho0IRpw, rho0MTsh, vEff0IRpw, vEff0IRpwUw, vEff0MTsh )
  call TimeNOstopNO('Init JuPhon')
  write(*, '(a)') 'JuPhon initialized!'

#ifdef DEBUG_MODE
  if (testsActivated) then
    ! Additional more elaborated tests are performed to further test quantities created by InitializeJuPhon routine.
    ! todo attention with bmat, this should not be transposed back for init tests!!!!! Different lattice constants might not work
    call TimeStart('Test: Initialization')
    call DetailedInitializationTests( atoms, stars, cell, dimens, lathar, sym, input, kpts, usdus, Veff0, testRadSolSw,            &
      & testKptsWeightSw, testCountValElecSw, testVeff0ContSw, testrho0ContSw, testBackRotMTCoordSysSw, testPsi0ContSw,            &
      & testOverlapSw, testXCintegrals, logUnit, noPtsCon, rbas1, rbas2, iloTable, nRadFun, gridf, kveclo, nv, GbasVec, ilst, ne, z, rho0IR, rho0MT,&
      & mlh_atom, nmem_atom, clnu_atom, vXC0IR, eXCIR, vXC0MT, eXCMT, ngdp, gdp )
    call TimeNOstopNO('Test: Initialization')

    ! Tests of the single terms which setup the Sternheimer equation.
    call TimeStart('Test: Sternheimer')
    call TestSternheimerSetup( atoms, cell, sym, stars, lathar, enpara, kpts, qpts, input, dimens, usdus, Veff0, results, iqpt,    &
      & ngdp, ngdp2km, paPoX, paPoY, paPoZ, harSw, extSw, xcSw, testCompareGrVeff0FleurSw, testVeff1Sw, testUnfoldStarsSw,         &
      & testRadDerivativeSw, testGauntCoeffSw, testGradLhExpandFuncSw, testContGrVeff0Sw, testWarpingSw, testSternhHSMEtestSw,     &
      & testSternhSchroedPertTheoSw, testz1Phi0ContSw, testRho1BasCorrSw, testPlotRho03Dsw, testGradRho0PathSw, testRho1IRsw,      &
      & testRho1MTsw, logUnit, gdp, mlh_atom, nmem_atom, clnu_atom, rho0IR, rho0MT, ne, nobd, eig, rbas1, rbas2, uuilon, duilon,   &
      & ulouilopn, GbasVec, ilst, nv, z, kveclo, iloTable, nRadFun, El, mapKpq2K, ilo2p, memd_atom, noPtsCon, gdp2Ind, gdp2iLim,   &
      & kpq2kPrVec, testrho0ContSw, testCompareSurfInt, testSplitMTSurfIntSterh,                &
      & testVeff1IRMESternh, testEps1q0, testVeff1IRMatqBackFold, testVeff1IRqLatPeriod, vEff0MTsh, vEff0IRpwUw )
    call TimeNOstopNO('Test: Sternheimer')

    call TimeStart('Test: Hellmann-Feynman dynamical matrix')
    call TestDynMatHF( atoms, lathar, input, stars, cell, dimens, sym, qpts, harSw, extSw, xcSw, testEii2PsDens, testEii2LatPeriodQSw, ngdp, paPoX,         &
      & paPoY, paPoZ, memd_atom, gdp, mlh_atom, nmem_atom, clnu_atom, vExt2IR, vExt2MT, rho0IR, rho0MT, logUnit )
    call TimeNOstopNO('Test: Hellmann-Feynman dynamical matrix')

    call TimeStart('Test: Pulay dynamical matrix')
    call TestDynMatPulay( atoms, dimens, cell, sym, stars, lathar, input, kpts, qpts, enpara, usdus, Veff0, results, ngdp, &
      & memd_atom, logUnit, testGrMatElemPsiHepsPsiGaussTheo, testPsiHepsTildePsi, testR2orNotWfMtGradNgrNTensGrOvls,              &
      & testVarphiHepsVarphiSw, testGrPsiPsiMatElem, testIntVeff1Rho1Val, test3rdPulDynMatCancel, testIRIntegralSw, testIR3rdMatElemSw, testActionHgrPhiSw, clnu_atom, nmem_atom, mlh_atom,  &
      & rho0IR, rho0MT, gdp, gBasVec, ilst, nv, ne, &
      & nobd, nRadFun, kpq2kPrVec, z, rbas1, rbas2, El, eig, uuilon, duilon, ulouilopn, ilo2p, iloTable, kveclo, mapKpq2K, vEff0MTsh )
    call TimeNOstopNO('Test: Pulay dynamical matrix')

    call TimeStart('Test: Surface dynamical matrix')
    call TestDynMatSurf( atoms, dimens, stars, sym, cell, kpts, qpts, input, lathar, usdus, Veff0, results, testComp3ArgSFIntsSw, testComp2ArgSFIntsSw, testComp2ArgGrSFIntsSw, nRadFun, El, eig,      &
                       & gdp, rbas1, rbas2, nmem_atom, mlh_atom, clnu_atom, z, nv, GbasVec, ilst, nobd, kveclo, iloTable, ngdp,    &
                       & kpq2kPrVec, mapKpq2K, logUnit, memd_atom, rho0IR, rho0MT )
    call TimeNOstopNO('Test: Surface dynamical matrix')

    call TimeStart('Test: Goldstein dynamical matrix')
    call TestGoldsteinCond( atoms, cell, lathar, stars, dimens, kpts, qpts, results, Veff0, sym, usdus, enpara, input, ngdp, &
      & memd_atom, logUnit, testGoldsteinRemaining, test1st2ndPulDynMatEps1, test1st2ndPulDynMatCancel, testGoldsteinSurfSw, rho0IR, rho0MT, clnu_atom, &
      & nmem_atom, mlh_atom, gdp, GbasVec, ilst, kveclo, nv, kpq2kPrVec, eig, nobd, &
      & mapKpq2K, nRadFun, rbas1, rbas2, El, iloTable, z, uuilon, duilon, ulouilopn, ilo2p )
    call TimeNOstopNO('Test: Goldstein dynamical matrix')

!    if (.false.) then
!      call TestDynMatDeprecated( atoms, enpara, lathar, sym, cell, kpts, dimens, usdus, input, results, qpts, stars, Veff0,        &
!        &  logUnit, ngdp, memd_atom, GbasVec, gdp2iLim, kpq2kPrVec, gdp2Ind, gdp, mapKpq2K,                          &
!        & rho0IR, rho0MT, nv, ilst, z, kveclo, nRadFun, rbas1, rbas2, iloTable, El, eig, nobd, clnu_atom, nmem_atom, mlh_atom,     &
!        & uuilon, duilon, ulouilopn, ilo2p, oneD, vacuum, ne )
!    end if

    write(*, '(a)') 'Tests finished'

  else
      write ( logUnit, * )
      write ( logUnit, * )
      write(*, *)
      write(*, '(a)') '----------------------------'
      write(*, '(a)') 'DISABLED input data test(s)!'
      write(*, '(a)') '----------------------------'
      write(logUnit, '(a)') 'DISABLED input data tests!'
      write(logUnit, '(a)') '**************************'
      write ( logUnit, * )
      write(*, '(a)') '-----------------------------------'
      write(*, '(a)') 'DISABLED density variation test(s)!'
      write(*, '(a)') '-----------------------------------'
      write(logUnit, '(a)') 'DISABLED density variation tests!'
      write(logUnit, '(a)') '*********************************'
      write(*, '(a)') '----------------------------'
      write(*, '(a)') 'DISABLED potential test(s)!'
      write(*, '(a)') '----------------------------'
      write ( logUnit, * )
      write(logUnit, '(a)') 'DISABLED potential tests!'
      write(logUnit, '(a)') '*************************'
      write ( logUnit, * )
      write ( logUnit, * )
      write(*, '(a)') '----------------------------'
      write(*, '(a)') 'DISABLED Sternheimer test(s)!'
      write(*, '(a)') '----------------------------'
      write(logUnit, '(a)') 'DISABLED Sternheimer tests!'
      write(logUnit, '(a)') '**************************'
      write(*, '(a)') '----------------------------------------------------------------'
      write(*, '(a)') 'DISABLED dynamical matrix Hellmann-Feynman contribution test(s)!'
      write(*, '(a)') '----------------------------------------------------------------'
      write(logUnit, '(a)')
      write(logUnit, '(a)') 'DISABLED dynamical matrix Hellmann-Feynman contribution test(s)!'
      write(logUnit, '(a)') '****************************************************************'
      write(*, '(a)') '-----------------------------------------------------'
      write(*, '(a)') 'DISABLED dynamical matrix Pulay contribution test(s)!'
      write(*, '(a)') '-----------------------------------------------------'
      write(logUnit, '(a)')
      write(logUnit, '(a)') 'DISABLED dynamical matrix Pulay contribution test(s)!'
      write(logUnit, '(a)') '*****************************************************'
      write(*, '(a)') '---------------------------------------------------'
      write(*, '(a)') 'DISABLED dynamical matrix surface integral test(s)!'
      write(*, '(a)') '---------------------------------------------------'
      write(logUnit, '(a)')
      write(logUnit, '(a)') 'DISABLED dynamical matrix surface integral test(s)!'
      write(logUnit, '(a)') '***************************************************'
      write(*, '(a)') '------------------------------------------------------'
      write(*, '(a)') 'DISABLED dynamical matrix Goldstein condition test(s)!'
      write(*, '(a)') '------------------------------------------------------'
      write(logUnit, *)
      write(logUnit, '(a)') 'DISABLED dynamical matrix Goldstein condition test(s)!'
      write(logUnit, '(a)') '******************************************************'
  end if

  if ( onlyTests )then
    call WriteTimes()
    call finishLogFile( startTime, logUnit )
    call fclose(logUnit, status='keep')
    NOstopNO
  end if
#endif

  write (*, '(a,1x,i3)') 'Performing calculation for q-point with index:', iqpt

  call TimeStart('Init Sternheimer SCC')
  ! We determine all quantities not dependent on the q-points
  call initSternheimerSCC( atoms, sym, cell, stars, dimens, lathar, enpara, usdus, input, tdHS0, qpts, logUnit, ngdp, iqpt, &
    & gdp, rho0IR, rho0MT, mlh_atom, nmem_atom, clnu_atom, rbas1, rbas2, uuilon, duilon, ulouilopn, ilo2p, rho0IRpw, rho0MTsh, &
    & vEff0MTsh, qpwcG, rho1MTCoreDispAt, grVxcIRKern, dKernMTGPts, grVeff0MT_init, grVeff0MT_main, grRho0IR, grRho0MT, gausWts, &
    & ylm, grVext0IR_DM, grVext0MT_DM, grVCoul0IR_DM_SF, grVCoul0MT_DM_SF, grVeff0IR_DM, grVeff0MT_DM, grVeff0MT_DMhxc )
  call TimeNOstopNO('Init Sternheimer SCC')
  write(*, '(a)') 'Sternheimer self-consistency cycle initialized!'

  plot3DgrVhar0 = .false.
  plot3DgrVext0 = .false.
  plot3DgrVxc0  = .false.
  plot3DgrVeff0 = .false.
  call PlotGradV0( atoms, cell, lathar, stars, dimens, input, sym, memd_atom, ngdp, logUnit, plot3DgrVhar0, &
    & plot3DgrVext0, plot3DgrVxc0, plot3DgrVeff0, gdp, rho0IR, rho0MT, nmem_atom, mlh_atom, clnu_atom, grRho0IR, grRho0MT, gausWts,&
    & ylm, dKernMTGPts, grVxcIRKern )

  ! Solve Sternheimer self-consistency cycle and obtain first order densities and first order potentials for all atoms
  write(*, '(a,1x,i3,a)') 'Sternheimer self-consistency cycle for q-point', iqpt, '...'
  call TimeStart('Sternheimer SCC')
  if (compPhon) then
    call solveSternheimerSCC( atoms, sym, stars, lathar, dimens, cell, enpara, usdus, input, kpts, qpts, results, usdus,      &
      & logUnit, ngdp, rbas1, rbas2, kveclo, uuilon, duilon, ulouilopn, &
      & gdp, mapKpq2K, ne, eig, GbasVec, ilst, z, nv, El, nradFun, iloTable, nobd, ilo2p, gdp2Ind,     &
      & gdp2iLim, kpq2kPrVec, qpwcG, iqpt, tdHS0, ylm, grRho0IR, grRho0MT, grVeff0MT_DM, grVeff0MT_main, dKernMTGPts,       &
      & grVxcIRKern, rho1MTCoreDispAt, gausWts, rho1IR, rho1MT, vExt1MT, vEff1IR, vEff1MT, oneSternhCycle, ngpqdp, gpqdp, &
      & vExt1IR_final, vHar1IR_final, vHar1MT_final, rho1MTDelta, vExt1MTDelta, vExt1MTq0, vHar1MTDelta, vHar1MTq0, vXc1MTDelta, &
      & vXc1MTq0, rho0IRpw, rho0MTsh, vEff0IRpwUw, noPtsCon, vExt1MTnoVol, vEff1MTnoVol, vH1MTnoVol, vExt1MTnoVolnoq, grVeff0MT_DMhxc, vExt1noqIR_final, rho1MTz0,  vCoul1IRtempNoVol, vCoul1MTtempNoVol, vEff0MTsh)
  else
    call solveSternheimerSCC( atoms, sym, stars, lathar, dimens, cell, enpara, usdus, input, kpts, qpts, results, usdus,      &
      & logUnit, ngdp, rbas1, rbas2, kveclo, uuilon, duilon, ulouilopn, &
      & gdp, mapKpq2K, ne, eig, GbasVec, ilst, z, nv, El, nradFun, iloTable, nobd, ilo2p, gdp2Ind,     &
      & gdp2iLim, kpq2kPrVec, qpwcG, iqpt, tdHS0, ylm, grRho0IR, grRho0MT, grVeff0MT_init, grVeff0MT_main, dKernMTGPts,       &! main --> DM for Alex
      & grVxcIRKern, rho1MTCoreDispAt, gausWts, rho1IR, rho1MT, vExt1MT, vEff1IR, vEff1MT, oneSternhCycle, ngpqdp, gpqdp,&
      & vExt1IR_final, vHar1IR_final, vHar1MT_final, rho1MTDelta, vExt1MTDelta, vExt1MTq0, vHar1MTDelta, vHar1MTq0, vXc1MTDelta, &
      & vXc1MTq0, rho0IRpw, rho0MTsh, vEff0IRpwUw, noPtsCon, vExt1MTnoVol, vEff1MTnoVol, vH1MTnoVol, vExt1MTnoVolnoq, grVeff0MT_DMhxc, vExt1noqIR_final, rho1MTz0, vCoul1IRtempNoVol, vCoul1MTtempNoVol )
  end if
  call TimeNOstopNO('Sternheimer SCC')
  write(*, '(a)') '1st-order quantities determined!'

  plotVhar1fccXZSw = .false.
  plotVext1fccXZSw = .false.
  plotVxc1fccXZSw = .false.
  plotVeff1fccXZSw = .false.
  plotRho1fccXZSw = .false.
  call plotConvergedQuantities( atoms, stars, cell, lathar, dimens, sym, input, qpts, ngdp, iqpt, memd_atom, logUnit, &
    & plotVhar1fccXZSw, plotVext1fccXZSw, plotVxc1fccXZSw, plotVeff1fccXZSw, gdp, rho0IR, rho0MT, nmem_atom, mlh_atom, clnu_atom, rho1IR, &
    & rho1MT, gausWts, ylm, dKernMTGPts, grVxcIRKern, rho1MTDelta, plotRho1fccXZSw )
#ifdef DEBUG_MODE
  !todo write interface to turn on and off this test environment and to NOstopNOafter this without calculating the dynamical matrix
  call TestConvergedQuantities( atoms, sym, lathar, dimens, stars, input, cell, kpts, qpts, enpara, usdus, Veff0, results, ngdp, ngdp2km, logUnit, memd_atom, numAddQs, gdp, mlh_atom, nmem_atom, clnu_atom, rho0IR, rho0MT, rbas1, rbas2, uuilon, duilon, ulouilopn, ilo2p, kveclo, mapKpq2K, ne, eig, gbasVec, ilst, z, nv, El, nRadFun, iloTable, nobd, kpq2kPrVec, gdp2Ind, rho0IRpw, rho0MTsh, vEff0MTsh, vEff0IRpwUw, noPtsCon )
#endif

  ! Calculate the second-order external potential and the second-order response of the ion-ion energy. Test-mode set false!
  call TimeStart('2nd order quantities')
  ! todo we leave it on true for production until we found out why we need to subtract the trace
!  call GenVext2(atoms, cell, dimens, ngdp, gdp, vExt2IR, vExt2MT, .false.)
  call CalcIIEnerg2(atoms, cell, dimens, qpts, stars, input, iqpt, ngdp, gdp, E2ndOrdII)
  call TimeNOstopNO('2nd order quantities')
  write(*, '(a)') '2nd-order quantities determined!'


  ! Set up Dynamical Matrix comprising of a basic Hellmann-Feynman contribution and a Pulay contribution
  call TimeStart('Setup Dynamical Matrix')
  call SetupDynamicMatrix( atoms, input, sym, dimens, cell, lathar, stars, kpts, qpts, usdus, results, Veff0, iqpt, ngdp, ngpqdp, gdp, mlh_atom, nmem_atom,&
    & clnu_atom, rho0IR, rho1IR, rho1MT, vExt1MT, vEff1IR, vEff1MT, Veff0%vpw, Veff0%vr(:, :, :, 1),&
    & rho0MT, E2ndOrdII, El, eig, rbas1, rbas2, iloTable, nv, nobd, ilst, GbasVec, z, kveclo, nRadFun, mapKpq2K, kpq2kPrVec,       &
    & gpqdp, memd_atom, logUnit, vXC0IR, eXCIR, vXC0MT, eXCMT, vExt1IR_final, vHar1IR_final, vHar1MT_final, grRho0IR, grRho0MT, &
    & grVext0IR_DM, grVext0MT_DM, grVeff0IR_DM, grVeff0MT_DM, dynMat, rho1MTDelta, vExt1MTDelta, vExt1MTq0, vHar1MTDelta, vHar1MTq0, &
    & vXc1MTDelta, vXc1MTq0, vExt1MTnoVol, grVeff0MT_DMhxc, vEff1MTnoVol, vH1MTnoVol, vExt1MTnoVolnoq, vExt1noqIR_final, rho1MTz0, &
    & grVCoul0IR_DM_SF, grVCoul0MT_DM_SF, vCoul1IRtempNoVol, vCoul1MTtempNoVol)
  call TimeNOstopNO('Setup Dynamical Matrix')
  write(*, '(a)') 'Setup of Dynamical Matrix done!'



  ! Wrapper for LAPACK Diagonalization routine
  call TimeStart('Diagonalize Dynamical Matrix')
  call DiagonalizeDynMat(atoms, qpts, calcEigenVec, dynMat, eigenVals, eigenVecs, iqpt)
  call TimeNOstopNO('Diagonalize Dynamical Matrix')
  write(*, '(a)') 'Dynamical Matrix diagonalized!'

  ! Transform the eigenvalues of the Dynamical Matrix to Phonon frequencies in units of THz or 1 / cm
  call TimeStart('Postprocessing')
  call CalculateFrequencies(atoms, iqpt, eigenVals, eigenFreqs)
  call TimeNOstopNO('Postprocessing')
  write(*, '(a)') 'Post-processing finished!'

  ! Finish juphon writing out the computation time, completing the logfile and closing it
  call WriteTimes()
  call finishLogFile( startTime, logUnit )
  call fclose(logUnit, status='keep')
  write(*, '(a)') 'juPhon terminated!'

end program juPhon
