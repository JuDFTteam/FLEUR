!-----------------------------------------------------------------------------------------------------------------------------------
! Forschungszentrum Jülich, juPhon Plugin for the FLEUR program
!-----------------------------------------------------------------------------------------------------------------------------------
!
! module m_jpInput: Initializes juPhon
!
!> @author
!> Christian-Roman Gerhorst and Markus Betzinger
!>
!> @brief
!> Initializes juPhon generating the input data for the linear response calculation.
!>
!> @details
!> This module contains routines which provide all needed input quantities for the calculation with juPhon. They are either parsed
!> from a preceeding FLEUR calculation or recalculated here for sake of better performance. Some routines also calculate input
!> quantities for juPhon which have not been calculated by FLEUR as they were not needed. The main routine operating the input is
!> m_input::initializejuphon. Some routines were recycled from former programs of Markus Betzinger.
!>
!> @note
!> Additional information and formulas pointing out the routines of this module can be found within this
!> <a href='jpInit.pdf'>document</a>.
!-----------------------------------------------------------------------------------------------------------------------------------
module m_jpInit

  implicit none

  contains

  !>--------------------------------------------------------------------------------------------------------------------------------
  !> @author
  !> Christian-Roman Gerhorst and Markus Betzinger
  !>
  !> @brief
  !> Main input method calling other routines for generating input data for juPhon.
  !>
  !> @details
  !> This routine operates the generation of the input data and calls other helping routines. Furthermore input data is adjusted to
  !> avoid confusion, i.e. bmat is transposed to not calculate with the transposed matrix but with the original reciprocal Bravais
  !> matrix.
  !>
  !> @param[out] atoms       : Atoms type, see types.f90.
  !> @param[out] cell        : Unit cell type, see types.f90.
  !> @param[out] sym         : Symmetries type, see types.f90.
  !> @param[out] stars       : Stars type, see types.f90.
  !> @param[out] lathar      : Lattice harmonics type, see types.f90.
  !> @param[out] input       : Input type, see types.f90.
  !> @param[out] dimens      : Dimension type, see types.f90.
  !> @param[out] enpara      : Energy parameter type, see types.f90.
  !> @param[out] xcpot       : Exchange potential type, see types.f90.
  !> @param[out] kpts        : K-points type, see types.f90.
  !> @param[out] qpts        : Q-points type, see types.f90.
  !> @param[out] Veff0       : Type containing unperturbed output potentials of Fleur, see types.f90.
  !> @param[out] usdus       : Type containing quantities consisting of the radial solutions, see types.f90.
  !> @param[out] results     : Results type, see types.f90.
  !> @param[in]  logUnit     : Unit number for juPhon.log.
  !> @param[out] paPoX       : X-coordinate of path's direction for potential's gradient test.
  !> @param[out] paPoY       : Y-coordinate of path's direction for potential's gradient test.
  !> @param[out] paPoZ       : Z-coordinate of path's direction for potential's gradient test.
  !> @param[out] harSw       : Switch to activate Hartree potential.
  !> @param[out] extSw       : Switch to activate external potential.
  !> @param[out] xcSw        : Switch to activate exchange-correlation potential.
  !> @param[out] noPtsCon   : Number of points for which continuity test should be performed.
  !> @param[out] ngdp        : Number of G-vectors for potentials and densities.
  !> @param[out] ngdp2km     : Number of G-vectors for potentials and densities which are smaller than 2 kmax.
  !> @param[out] memd_atom   : Maximal number of members in all lattice harmonics for all atoms.
  !> @param[out] numaddQs    : Number of additional q-vectors
  !> @param[in]  startTime   : Starting time of juPhon.
  !> @param[out] kveclo      : Basis G-vectors of local orbitals.
  !> @param[out] nv          : Number of LAPW G-basis vectors for given k-point.
  !> @param[out] nRadFun     : Number of radial functions per orbital quantum number l and atom type.
  !> @param[out] iloTable    : Number of local orbital if orbital quantum number l, atom type and index p > 2 from nRadFun is given.
  !> @param[out] GbasVec     : G-basis vectors
  !> @param[out] ilst        : For various k-points G-basis vectors occur more than once, thus they are only stored once in juPhon.
  !>                           This pointer array contains the right index for GbasVec array to "unfold" G-basis vectors again.
  !> @param[out] ne          : Number of eigenvalues per k-point.
  !> @param[out] mapKpq2K    : For a given k-point index and q-point index, this array gives the k-point set index of the result
  !>                           k + q (already mapped back to Brillouin zone).
  !> @param[out] eig         : Contains Kohn\--Sham eigenvalues.
  !> @param[out] rbas1       : Large components of radial solution, its energy derivative and u_LO
  !> @param[out] rbas2       : Small components of radial solution, its energy derivative and u_LO
  !> @param[out] El          : Contains LAPW and LO energy parameters.
  !> @param[out] vrCoul      : Unperturbed radial coefficients of the converged Coulomb potential parsed from FLEUR.
  !> @param[out] z           : Kohn-Sham eigenvectors.
  !> @param[out] vpwCoul_uw  : Star coefficients of the unwarped (no convolution with step function) and converged interstitial
  !>                           Coulomb potential parsed from Fleur.
  !> @param[out] rho0IR      : Star coefficients of the unperturbed and converged interstitial density parsed from Fleur.
  !> @param[out] rho0MT      : Radial coefficients of the unperturbed and converged muffin-tin densities parsed from Fleur.
  !> @param[out] gridf       : Initialized grid quantity to use intgrf routine.
  !> @param[out] gdp         : G-vectors of potentials and densities.
  !> @param[out] mlh_atom    : Magnetic quantum number m of lattice harmonic members for every atom.
  !> @param[out] nmem_atom   : Number of lattice harmonic members for every atom.
  !> @param[out] clnu_atom   : Phase mediating between stars and plane waves.
  !> @param[out] uuilon      : overlap integral between the radial functions of the integral (multiplied by ulo_der) of a local
  !>                           orbital and the flapw radial function with the same l
  !> @param[out] duilon      : overlap integral between the radial functions of the integral of a local orbital and the energy
  !>                           derivative of the flapw radial function with the same l
  !> @param[out] ulouilopn   : overlap integral between the radial functions of the integral of a local orbital and another local
  !>                           orbital with the same l.
  !> @param[out] nobd        : Number of occupied bands per k-point and spin
  !> @param[out] ilo2p       : mapping array giving the p value for given number of LO and itype
  !> @param[out] gdp2Ind     : Stores the index of a potential and density G-Vector. The dimensions are the G-vector components.
  !> @param[out] gdp2iLim    : Stores the min and maxvals of gdp2Ind
  !>
  !>--------------------------------------------------------------------------------------------------------------------------------
  subroutine InitializeJuPhon( atoms, cell, sym, stars, lathar, input, dimens, enpara, xcpot, kpts, qpts, Veff0, usdus, results,   &
      & logUnit, iqpt, paPoX, paPoY, paPoZ, harSw, extSw, xcSw, noPtsCon, ngdp, ngdp2km, memd_atom, numAddQs, startTime, kveclo,   &
      & nv, nRadFun, iloTable, GbasVec, ilst, ne, mapKpq2K, eig, rbas1, rbas2, El, vrCoul, z, vpwCoul_uw, rho0IR, rho0MT, gridf,   &
      & gdp, mlh_atom, nmem_atom, clnu_atom, uuilon, duilon, ulouilopn, nobd, ilo2p, gdp2Ind, gdp2iLim, kpq2kPrVec, calcEigenVec,  &
      & oneSternhCycle, recipLengthUnit, onlyTests, testCompareGrVeff0FleurSw, testVeff1Sw, testUnfoldStarsSw, testRadDerivativeSw,&
      & testGauntCoeffSw, testGradLhExpandFuncSw, testContGrVeff0Sw, testWarpingSw, testSternhHSMEtestSw,                          &
      & testSternhSchroedPertTheoSw, testz1Phi0ContSw, testRho1BasCorrSw, testPlotRho03Dsw, testRadSolSw, testKptsWeightSw,        &
      & testCountValElecSw, testVeff0ContSw, testrho0ContSw, testBackRotMTCoordSysSw, testPsi0ContSw, testOverlapSw,               &
      & testGradRho0PathSw, testEii2LatPeriodQSw, testVarphiHepsVarphiSw, testRho1IRsw, testRho1MTsw, testsActivated,   vacuum,&
      & test1st2ndPulDynMatEps1, test1st2ndPulDynMatCancel, test3rdPulDynMatCancel, testIntVeff1Rho1Val, testGrPsiPsiMatElem,      &
      & testCompareSurfInt, testSplitMTSurfIntSterh, testVeff1IRMESternh, testEps1q0, testVeff1IRMatqBackFold,                     &
      & testVeff1IRqLatPeriod, testGrMatElemPsiHepsPsiGaussTheo, testPsiHepsTildePsi, testGoldsteinRemaining,                      &
      & testR2orNotWfMtGradNgrNTensGrOvls, testComp3ArgSFIntsSw, testComp2ArgSFIntsSw, testGoldsteinSurfSw, testComp2ArgGrSFIntsSw,&
      & testIRIntegralSw, testIR3rdMatElemSw, testActionHgrPhiSw, testXCintegrals, testEii2PsDens, vXC0IR, eXCIR, vXC0MT, eXCMT,   &
      & rho0IRpw, rho0MTsh, vEff0IRpw, vEff0IRpwUw, vEff0MTsh )

    use m_types
    use m_firstglance
    use m_jpParamInpF, only : ReadInpFile
    use m_jpLog, only : LogFleurInit, InpParam2Log, LogRotLh, LogBandsBFcn, InitLogFile, FinishLogFile
    use m_fleur_init
    use m_rwinp
    use mod_juPhonUtils, only : Intgrf_init
    use m_jpPotDensHelper, only : genPotDensGvecs, genRotLh
    use m_juDFT_NOstopNO, only : juDFT_error, juDFT_warn
    use mod_juPhonUtils,     only : fopen, fclose

#include "cppmacro.h"

    implicit none

    ! Type parameters
    type(t_atoms),                    intent(out) :: atoms
    type(t_cell),                     intent(out) :: cell
    type(t_sym),                      intent(out) :: sym
    type(t_stars),                    intent(out) :: stars
    type(t_sphhar),                   intent(out) :: lathar
    type(t_input),                    intent(out) :: input
    type(t_dimension),                intent(out) :: dimens
    type(t_enpara),                   intent(out) :: enpara
    type(t_xcpot),                    intent(out) :: xcpot
    type(t_kpts),                     intent(out) :: kpts
    type(t_kpts),                     intent(out) :: qpts
    type(t_potential),                intent(out) :: Veff0
    type(t_usdus),                    intent(out) :: usdus
    type(t_results),                  intent(out) :: results

    ! Scalar parameters
    integer,                          intent(in)  :: logUnit !
    integer,                          intent(in)  :: startTime(8)
    integer,                          intent(out) :: iqpt
    real,                             intent(out) :: paPoX
    real,                             intent(out) :: paPoY
    real,                             intent(out) :: paPoZ
    logical,                          intent(out) :: harSw
    logical,                          intent(out) :: extSw
    logical,                          intent(out) :: xcSw
    logical,                          intent(out) :: calcEigenVec
    logical,                          intent(out) :: oneSternhCycle
    logical,                          intent(out) :: recipLengthUnit
    logical,                          intent(out) :: onlyTests
    logical,                          intent(out) :: testsActivated
    integer,                          intent(out) :: noPtsCon
    integer,                          intent(out) :: ngdp
    integer,                          intent(out) :: ngdp2km
    integer,                          intent(out) :: memd_atom
    integer,                          intent(out) :: numAddQs
    logical,                          intent(out) :: testCompareGrVeff0FleurSw
    logical,                          intent(out) :: testVeff1Sw
    logical,                          intent(out) :: testUnfoldStarsSw
    logical,                          intent(out) :: testRadDerivativeSw
    logical,                          intent(out) :: testGauntCoeffSw
    logical,                          intent(out) :: testGradLhExpandFuncSw
    logical,                          intent(out) :: testContGrVeff0Sw
    logical,                          intent(out) :: testWarpingSw
    logical,                          intent(out) :: testSternhHSMEtestSw
    logical,                          intent(out) :: testSternhSchroedPertTheoSw
    logical,                          intent(out) :: testz1Phi0ContSw
    logical,                          intent(out) :: testRho1BasCorrSw
    logical,                          intent(out) :: testPlotRho03Dsw
    logical,                          intent(out) :: testRadSolSw
    logical,                          intent(out) :: testKptsWeightSw
    logical,                          intent(out) :: testCountValElecSw
    logical,                          intent(out) :: testVeff0ContSw
    logical,                          intent(out) :: testrho0ContSw
    logical,                          intent(out) :: testBackRotMTCoordSysSw
    logical,                          intent(out) :: testRho1IRsw
    logical,                          intent(out) :: testRho1MTsw
    logical,                          intent(out) :: testPsi0ContSw
    logical,                          intent(out) :: testOverlapSw
    logical,                          intent(out) :: testGradRho0PathSw
    logical,                          intent(out) :: testEii2LatPeriodQSw
    logical,                          intent(out) :: testVarphiHepsVarphiSw
    logical,                          intent(out) :: test1st2ndPulDynMatEps1
    logical,                          intent(out) :: test1st2ndPulDynMatCancel
    logical,                          intent(out) :: test3rdPulDynMatCancel
    logical,                          intent(out) :: testIntVeff1Rho1Val
    logical,                          intent(out) :: testGrPsiPsiMatElem
    logical,                          intent(out) :: testCompareSurfInt
    logical,                          intent(out) :: testSplitMTSurfIntSterh
    logical,                          intent(out) :: testVeff1IRMESternh
    logical,                          intent(out) :: testEps1q0
    logical,                          intent(out) :: testVeff1IRMatqBackFold
    logical,                          intent(out) :: testVeff1IRqLatPeriod
    logical,                          intent(out) :: testGrMatElemPsiHepsPsiGaussTheo
    logical,                          intent(out) :: testPsiHepsTildePsi
    logical,                          intent(out) :: testGoldsteinRemaining
    logical,                          intent(out) :: testR2orNotWfMtGradNgrNTensGrOvls
    logical,                          intent(out) :: testComp3ArgSFIntsSw
    logical,                          intent(out) :: testComp2ArgSFIntsSw
    logical,                          intent(out) :: testGoldsteinSurfSw
    logical,                          intent(out) :: testComp2ArgGrSFIntsSw
    logical,                          intent(out) :: testIRIntegralSw
    logical,                          intent(out) :: testIR3rdMatElemSw
    logical,                          intent(out) :: testActionHgrPhiSw
    logical,                          intent(out) :: testXCintegrals
    logical,                          intent(out) :: testEii2PsDens

    ! Array parameters
    integer,            allocatable,  intent(out) :: kveclo(:,:)
    integer,            allocatable,  intent(out) :: nv(:, :)
    integer,            allocatable,  intent(out) :: nRadFun(:,:)
    integer,            allocatable,  intent(out) :: iloTable(:, :, :)
    integer,            allocatable,  intent(out) :: GbasVec(:, :)
    integer,            allocatable,  intent(out) :: ilst(:, :, :)
    integer,            allocatable,  intent(out) :: ne(:)
    integer,            allocatable,  intent(out) :: mapKpq2K(:, :)
    real,               allocatable,  intent(out) :: eig(:,:,:)
    real,               allocatable,  intent(out) :: rbas1(:,:,:,:,:)
    real,               allocatable,  intent(out) :: rbas2(:,:,:,:,:)
    real,               allocatable,  intent(out) :: El(:, :, :, :)
    real,               allocatable,  intent(out) :: vrCoul(:, :, :, :)
    MCOMPLEX,           allocatable,  intent(out) :: z(:,:,:,:)
    complex,            allocatable,  intent(out) :: vpwCoul_uw(:, :)
    complex,            allocatable,  intent(out) :: rho0IR(:,:)
    real,               allocatable,  intent(out) :: rho0MT(:,:,:,:)
    real,               allocatable,  intent(out) :: gridf(:, :)
    integer,            allocatable,  intent(out) :: gdp(:,:)
    integer,            allocatable,  intent(out) :: mlh_atom(:,:,:)
    integer,            allocatable,  intent(out) :: nmem_atom(:, :)
    complex,            allocatable,  intent(out) :: clnu_atom(:,:,:)
    real,               allocatable,  intent(out) :: uuilon(:, :)
    real,               allocatable,  intent(out) :: duilon(:, :)
    real,               allocatable,  intent(out) :: ulouilopn(:, :, :)
    integer,            allocatable,  intent(out) :: nobd(:, :)
    integer,            allocatable,  intent(out) :: ilo2p(:, :)
    integer,            allocatable,  intent(out) :: gdp2Ind(:, :, :)
    integer,            allocatable,  intent(out) :: kpq2kPrVec(:, :, :)
    complex,            allocatable,  intent(out) :: vXC0IR(:, :)
    complex,            allocatable,  intent(out) :: eXCIR(:)
    real,               allocatable,  intent(out) :: vXC0MT(:, :, :, :)
    real,               allocatable,  intent(out) :: eXCMT(:, :, :)
    integer,                          intent(out) :: gdp2iLim(2, 3)
    complex,            allocatable,  intent(out) :: rho0IRpw(:, :)
    complex,            allocatable,  intent(out) :: vEff0IRpw(:, :)
    complex,            allocatable,  intent(out) :: vEff0IRPwUw(:, :)
    complex,            allocatable,  intent(out) :: rho0MTsh(:, :, :, :)
    complex,            allocatable,  intent(out) :: vEff0MTsh(:, :, :, :)


    ! Scalar Variables
    !
    ! kMode_sw        : if true, a Gamma including equidistant k-point set is created
    ! writeKpqArraySw : if true, write kpqMappingArray to file
    ! ivers           : Specifies version of FLEUR at top of FLEUR out file
    ! l_opti          : Irrelevant for juPhon
    ! nw              : Irrelevant for juPhon
    ! imixJP          : Index of mixing method. (0 means simple mixing)
    ! numAddQs        : Number of additional phonon vectors q
    ! namex           : Irrelevant for juPhon
    ! relcor          : Irrelevant for juPhon
    ! mixAlpha        : Coefficient alpha determing the grade of mixing
    ! latVecScal      : Irrelevant for juPhon
    ! log_dump        : Irrelevant for juPhon
    ! cdn1Exist       : Logical, whether density file is existing
    logical                                       :: kMode_sw
    logical                                       :: writeKpqArraySw
    character(len=9),   parameter                 :: ivers='fleurInit'
    logical                                       :: l_opti
    integer                                       :: nw = 1
    integer                                       :: imixJP
    character(len=4)                              :: namex
    character(len=12)                             :: relcor
    real                                          :: mixAlpha
    real                                          :: latVecScal
    logical                                       :: log_dump
    logical                                       :: cdn1Exist

    ! Type Variables
    ! All these type variables are irrelevant for juPhon
    type(t_mpi)                                   :: mpi
    type(t_noco)                                  :: noco
    type(t_vacuum)     ,intent(out)                           :: vacuum
    type(t_sliceplot)                             :: sliceplot
    type(t_banddos)                               :: banddos
    type(t_obsolete)                              :: obsolete
    type(t_jij)                                   :: jij
    type(t_hybrid)                                :: hybrid
     

    ! Array Variables
    ! uds        : energy derivative of radial solution u at muffin-tin radius
    ! us         : radial solution at muffin-tin radius
    ! duds       : energy derivative of radial derivative of radial solution u at muffin-tin radius
    ! dus        : radial derivative of radial solution at muffin-tin radius
    ! ddn        : overlap integral of udots
    ! ulos       : the value of the radial function of a local orbital at the muffin tin radius
    ! dulos      : the value of the radial derivative of the radial function of a local orbital at the muffin tin radius
    ! uulon      : overlap integral between the radial functions of a local obital and the flapw radial function with the same l
    ! dulon      : overlap integral between the radial functions of a local obital and the energy derivative of the flapw radial
    !              function with the same l
    ! addQs      : sequential storage of additional phonon q-vectors
    ! kSetDim    : dimensions of additional k-point set
    ! qSetDim    : dimensions of additioanl q-point set
    ! kSetShift  : shift of origianl k-point set
    ! nkptiAqs   : number of irreducible additional k-points due to new phonon wave-vector q
    ! nkptfAqs   : total number of additional k-points due to new phonon wave-vector q
    ! noel       : irrelevant for juPhon calculation
    ! a1         : lattice vector 1
    ! a2         : lattice vector 2
    ! a3         : lattice vector 3
    ! intArr_dump: irrelevant for juPhon calculation
    real,               allocatable               :: uds(:, :)
    real,               allocatable               :: us(:, :)
    real,               allocatable               :: duds(:, :)
    real,               allocatable               :: dus(:, :)
    real,               allocatable               :: ddn(:, :)
    real,               allocatable               :: ulos( :, :)
    real,               allocatable               :: dulos( :, :)
    real,               allocatable               :: uulon( :, :)
    real,               allocatable               :: dulon( :, :)
    real,               allocatable               :: addQs(:)
    integer,            allocatable               :: kSetDim(:)
    integer,            allocatable               :: qSetDim(:)
    real,               allocatable               :: kSetShift(:)
    integer,            allocatable               :: nkptiAqs(:)
    integer,            allocatable               :: nkptfAqs(:)
    character(len=3),   allocatable               :: noel(:)
    real                                          :: a1(3)
    real                                          :: a2(3)
    real                                          :: a3(3)
    integer                                       :: intArr_dump(3)
    character(len=15)                             :: filenameTemp

    !integer :: iband, iG
    ! Parsing juPhon input file juPhon.inp and writing relevant parameters to logfile juPhon.log
    call TimeStart('Init:InpFile')
    call ReadInpFile( addQs, kSetDim, qSetDim, kMode_sw, kSetShift, iqpt, imixJP, mixAlpha, noPtsCon, paPoX, paPoY, paPoZ, harSW,  &
      & extSw, xcSw, writeKpqArraySw, calcEigenVec, oneSternhCycle, recipLengthUnit, onlyTests, testCompareGrVeff0FleurSw,         &
      & testVeff1Sw, testUnfoldStarsSw, testRadDerivativeSw, testGauntCoeffSw, testGradLhExpandFuncSw, testContGrVeff0Sw,          &
      & testWarpingSw, testSternhHSMEtestSw, testSternhSchroedPertTheoSw, testz1Phi0ContSw, testRho1BasCorrSw, testPlotRho03Dsw,   &
      & testRadSolSw, testKptsWeightSw, testCountValElecSw, testVeff0ContSw, testrho0ContSw, testBackRotMTCoordSysSw,              &
      & testPsi0ContSw, testOverlapSw, testGradRho0PathSw, testEii2LatPeriodQSw, testVarphiHepsVarphiSw, testRho1IRsw,             &
      & testRho1MTsw, testsActivated, test1st2ndPulDynMatEps1, test1st2ndPulDynMatCancel, test3rdPulDynMatCancel,                  &
      & testIntVeff1Rho1Val, testGrPsiPsiMatElem, testCompareSurfInt, testSplitMTSurfIntSterh, testVeff1IRMESternh, testEps1q0,    &
      & testVeff1IRMatqBackFold, testVeff1IRqLatPeriod, testGrMatElemPsiHepsPsiGaussTheo, testPsiHepsTildePsi,                     &
      & testGoldsteinRemaining, testR2orNotWfMtGradNgrNTensGrOvls, testComp3ArgSFIntsSw, testComp2ArgSFIntsSw, testGoldsteinSurfSw,&
      & testComp2ArgGrSFIntsSw, testIRIntegralSw, testIR3rdMatElemSw, testActionHgrPhiSw, testXCintegrals, testEii2PsDens )

    ! Open logfile and init with start time and hostname
    if (iqpt < 10) then
      write(filenameTemp, '("juPhon_q00",i1,".log")') iqpt
    else if (iqpt > 9 .and. iqpt < 100) then
      write(filenameTemp, '("juPhon_q0",i2".log")') iqpt
    else if (iqpt > 99 .and. iqpt < 1000) then
      write(filenameTemp, '("juPhon_q",i3".log")') iqpt
    end if
    call fopen( logUnit, name=filenameTemp, status='replace', action='write', form='formatted', access='sequential' )
    call initLogFile( logUnit, startTime )
    ! Obtain Input Data for juPhon Calculation
    write( logUnit, * )
    write( logUnit, * )
    write( logUnit, '(a)' ) 'juPhon.inp read! Start juPhon initialization...'
    write( logUnit, '(a)' ) '***********************************************'
    call InpParam2Log( kMode_sw, kSetDim, qSetDim, kSetShift, addQs, noPtsCon, [paPoX, paPoY, PaPoZ], harSw, extSw, xcSw, &
      & writeKpqArraySw, logUnit )
    call TimeNOstopNO('Init:InpFile')

    ! If k-points generation mode activated, obtain required input information and generate Γ-point centered k-point mesh
    if ( kMode_sw ) then
      call TimeStart('Init:kpts-genPrep')
      write( logUnit, * )
      write( logUnit, '(a)' ) 'Entering kpts-generation mode...'
      write( logUnit, '(a)' ) '================================'
      write( logUnit, * )
      write( logUnit, '(a)' ) 'Collecting minimal symmetry and cell information for k-point generation...'
      write( logUnit, '(a)' ) '--------------------------------------------------------------------------'

      ! This check catches the case that kpts generator is started before starting density is created by Fleur.
      inquire ( file='cdn1', exist=cdn1Exist )
      if ( .not.cdn1Exist ) then
        call juDFT_error( 'Files needed for generation of kpts set missing', calledby='operateInput',                              &
          & hint='Please run FLEUR first.' )
      end if

      ! Call of first_glance is necessary, otherwise noel cannot be allocated which is needed for rw_inp
      call First_glance( atoms%ntype,sym%nop,atoms%natd,obsolete%nwdd,atoms%nlod,vacuum%layerd, input%itmax,log_dump,log_dump,     &
        & log_dump, kpts%nkpt,kpts%nmop,jij%nqpt,intArr_dump )

      ! Allocate and initiate further quantities needed for the call of rw_inp
      atoms%ntypd = atoms%ntype
      atoms%nlod = max( atoms%nlod, 1 )
      allocate( atoms%lmax(atoms%ntype),atoms%ntypsy(atoms%natd),atoms%neq(atoms%ntype),atoms%nlhtyp(atoms%ntype), &
        & atoms%rmt(atoms%ntype),atoms%zatom(atoms%ntype),atoms%jri(atoms%ntype),atoms%dx(atoms%ntype), atoms%nlo(atoms%ntype), &
        & atoms%llo(atoms%nlod,atoms%ntype),atoms%nflip(atoms%ntype),atoms%bmu(atoms%ntype), noel(atoms%ntype), &
        & vacuum%izlay(vacuum%layerd,2),atoms%ncst(atoms%ntype),atoms%lnonsph(atoms%ntype), atoms%taual(3,atoms%natd), &
        & atoms%pos(3,atoms%natd), atoms%nz(atoms%ntype),atoms%relax(3,atoms%ntype), atoms%l_geo(atoms%ntype), &
        & noco%soc_opt(atoms%ntype+2),noco%alph(atoms%ntype),noco%beta(atoms%ntype), atoms%lda_u(atoms%ntype), &
        & noco%l_relax(atoms%ntype),jij%l_magn(atoms%ntype),jij%M(atoms%ntype), jij%magtype(atoms%ntype),jij%nmagtype(atoms%ntype),&
        & noco%b_con(2,atoms%ntype), lathar%clnu(1,1,1),lathar%nlh(1),lathar%llh(1,1),lathar%nmem(1,1),lathar%mlh(1,1,1), &
        & hybrid%select1(4,atoms%ntype),hybrid%lcutm1(atoms%ntype),hybrid%select2(4,atoms%ntype),hybrid%lcutm2(atoms%ntype),&
        & hybrid%lcutwf(atoms%ntype) )

      ! read in inp file to aquire Bravais matrix
      a1(:) = 0
      a2(:) = 0
      a3(:) = 0
      call Rw_inp( 'r', atoms, obsolete, vacuum, input, stars, sliceplot, banddos, cell, sym, xcpot, noco, jij,   hybrid, kpts,&
        & noel, namex, relcor, a1, a2, a3, latVecScal )
      cell%amat(:, 1) =  a1(:) * latVecScal
      cell%amat(:, 2) =  a2(:) * latVecScal
      cell%amat(:, 3) =  a3(:) * latVecScal

      ! Collect remaining symmetry information for the generation of k-points
      call PrepareKpointGeneration( sym, cell, logUnit )
      call DetermineNlotot( atoms, logUnit )
      call TimeNOstopNO('Init:kpts-genPrep')
      call TimeNOstopNO('Init JuPhon') ! is this really necessary?
    else
      write( logUnit, * )
      write( logUnit, * )
      write( logUnit, '(a)' ) 'Entering regular juPhon calculation mode...'
      write( logUnit, '(a)' ) '==========================================='
      write( logUnit, * )

      ! Recycle Fleur initialization routine. It might be modernized in the meanwhile.
      call TimeStart('Init:fleur-init')
      call fleur_init(ivers, mpi, input, dimens, atoms, lathar, cell, stars, sym, &
                     &noco, vacuum, sliceplot, banddos, obsolete, enpara, xcpot, results, &
                     &jij, kpts, hybrid,   l_opti)
      call TimeNOstopNO('Init:fleur-init')


      if ( input%jspins > 1 ) then
        call juDFT_error('Spins are not implemented', calledby='operateInput', hint='Implement spins before running again with this&
          & setting.' )
      end if
      input%spinf = 1
      input%imix = imixJP
      input%alpha = mixAlpha
      atoms%n_u = 0 ! no LDA + U

      if ( any ( abs( atoms%lmax - atoms%lnonsph ) /= 0 ) ) then
        call juDFT_warn('Non-spherical l is not equals lmax', calledby='operateInput', hint='Set non-spherical l to lmax')
      end if

      kpts%nkptf = kSetDim(1) * kSetDim(2) * kSetDim(3)
      kpts%nkpt3 = kSetDim
      if ( atoms%natd /= atoms%nat ) then
        call juDFT_error('natd not equals nat', calledby='operateInput', hint='This can cause dimension problems in read_eig' )
      end if

      call logFleurInit(atoms, input, dimens, cell, lathar, stars, sym, xcpot, logUnit)

      ! Generate inverse symmetry operations
      call TimeStart('Init:kpts-genPrep')
      call GenerateRecSymOp( sym )

      write ( logUnit, '(a)' ) 'Generation of inverse symmetry operations completed'
      write ( logUnit, '(a)' ) '===================================================='
      write ( logUnit, *)
      write ( logUnit, *)
    endif

    ! If kpts-generation mode kpts are generated otherwise they are read in and checked.
    call GenNCheckKnQPointInfo(qSetDim, logUnit, kSetDim, addQs, kMode_sw, writeKpqArraySw, sym, nkptiAqs, numAddQs, startTime, &
     & nkptfAqs, kpts, qpts, mapKpq2K, kpq2kPrVec)
    call TimeNOstopNO('Init:kpts-genPrep')

    ! Read relevant information from Fleur eig file and write some quantities to logfile
    call TimeStart('Init:readEigPotDens')
    call Read_eig( atoms, dimens, kpts, input, enpara, cell, eig, kveclo, z, nv, ilst, GbasVec, ne, El, logUnit )
!    if (sym%invs) then
!      z=0.5*(z+conjg(z))
!    end if

    ! Read converged potentials from Fleur and write some quantities to logfile
    call ReadPotFleur( atoms, stars, dimens, lathar, input, sym, vacuum,   Veff0, vrCoul, vpwCoul_uw, logUnit )

    ! Read in xc-potential and xc-energy density both in the IR and the MT
    call ReadXCquantities( atoms, stars, dimens, lathar, vXC0IR, eXCIR, vXC0MT, eXCMT )

    ! Read converged density from Fleur and write some quantities to logfile
    call ReadDensFleur( atoms, sym, lathar, dimens, stars, input, cell, vacuum,   logUnit, rho0IR, rho0MT )
    call TimeNOstopNO('Init:readEigPotDens')

    ! Generate radial solutions u, udot and u_LO and rearrange them in rbas1 and rbas2. Furthermore, generate overlap integrals of
    ! them.
    call TimeStart('Init:jPSpecial')
    call GenRadSolnInt( atoms, input, dimens, enpara, Veff0, mpi, usdus, logUnit, rbas1, rbas2, nRadFun, iloTable, uuilon, duilon, &
      & ulouilopn, ilo2p )


    !Transposing bmat due to practical reaasons
    cell%bmat =  transpose(cell%bmat)
    write ( logUnit, * )
    write ( logUnit, * )
    write ( logUnit, '(a)' ) 'Transposed reciprocal transposed Bravais matrix to avoid confusion within juPhon calculations'
    write ( logUnit, '(a)' ) '============================================================================================='

    ! Initialize fast numerical integration intgrf
    call Intgrf_init( atoms%ntype, atoms%jmtd, atoms%jri, atoms%dx, atoms%rmsh, gridf )
    write ( logUnit, * )
    write ( logUnit, '(a)' ) 'Initialized intgrf integration routine'
    write ( logUnit, '(a)' ) '======================================'

    ! generate G-vectors which are needed for calculating potentials and densities from now on
    write ( logUnit, * )
    write ( logUnit, '(a)' ) 'Creating G-vector set for potentials and densities'
    write ( logUnit, '(a)' ) '=================================================='
    call genPotDensGvecs(stars, cell, input, ngdp, ngdp2km, gdp, gdp2Ind, gdp2iLim, .false.)
    write ( logUnit, '(a,1x,i6,1x,a)' ) 'Found', ngdp, 'G-vectors.'

    ! For atoms of equal atom types, create coefficients for back-rotated local coordinate systems to use unperturbed densities and
    ! potentials converged in Fleur, which are expanded in lattice harmonics, for calculations in juPhon lacking any symmetry (see
    ! details in documenation of this routine).
    call GenRotLh( atoms, sym, cell, lathar, clnu_atom, memd_atom, nmem_atom, mlh_atom )
    call LogRotLh( atoms, lathar, memd_atom, logUnit, clnu_atom, nmem_atom, mlh_atom )

    call UnfoldRho0nVeff0Symmetry( atoms, stars, lathar, ngdp, gdp, clnu_atom, nmem_atom, mlh_atom, rho0IR, Veff0%vpw,             &
                                           & Veff0%vpw_uw, rho0MT, Veff0%vr, rho0IRpw, vEff0IRpw, vEff0IRPwUw, rho0MTsh, vEff0MTsh, cell)

    ! Determine number of occupied bands for each k-point
    call DetermineOccBands( atoms, cell, dimens, input, kpts, results, logUnit, nobd )
    call LogBandsBFcn( kpts, logUnit, nv, ne, nobd )

    write (logUnit, * )
    write (logUnit, * )
    write (logUnit, * )
    write (logUnit, '(a)') 'Acquiring input data completed!'
    write (logUnit, '(a)') '*******************************'
    call TimeNOstopNO('Init:jPSpecial')


  end subroutine InitializeJuPhon

  !>--------------------------------------------------------------------------------------------------------------------------------
  !> @author
  !> Christian-Roman Gerhorst
  !>
  !> @brief
  !> Unfolds the effective potential and the charge density read in from FLEUR from star to plane-wave and form lattice harmonic to
  !> spherical harmonic expansion.
  !>
  !> @details
  !> Factors that are multiplied to the potential and the density are already neutralized in the routine that reads them in
  !>
  !> @param[in]    atoms      : Atoms type, see types.f90.
  !> @param[in]    logUnit    : Unit number for juPhon.log.
  !>
  !>--------------------------------------------------------------------------------------------------------------------------------
  subroutine UnfoldRho0nVeff0Symmetry( atoms, stars, lathar, ngdp, gdp, clnu_atom, nmem_atom, mlh_atom, rho0IRst, vEff0IRst,       &
      & vEff0IRstUw, rho0MTlh, vEff0MTlh, rho0IRpw, vEff0IRpw, vEff0IRPwUw, rho0MTsh, vEff0MTsh, cell )

    use m_types, only : t_atoms, t_stars, t_sphhar, t_cell
    use m_jpPotDensHelper, only : ConvertStar2G, ConvertLH2SphHarm
    use mod_juPhonUtils, only : calcGrFinLH
    use m_jpConstants, only : compPhon

    implicit none

    ! Type arguments
    type(t_atoms),               intent(in)  :: atoms
    type(t_stars),               intent(in)  :: stars
    type(t_sphhar),              intent(in)  :: lathar
    type(t_cell), optional, intent(in)       :: cell

    ! Scalar arguments
    integer,                     intent(in)  :: ngdp

    ! Array Arguments
    integer,                     intent(in)  :: gdp(:, :)
    integer,                     intent(in)  :: mlh_atom(:,0:,:)
    integer,                     intent(in)  :: nmem_atom(0:, :)
    complex,                     intent(in)  :: clnu_atom(:,0:,:)
    complex,                     intent(in)  :: rho0IRst(:, :)
    complex,                     intent(in)  :: vEff0IRst(:, :)
    complex,                     intent(in)  :: vEff0IRstUw(:, :)
    real,                        intent(in)  :: rho0MTlh(:, 0:, :, :)
    real,                        intent(in)  :: vEff0MTlh(:, 0:, :, :)
    complex,  allocatable,       intent(out) :: rho0IRpw(:, :)
    complex,  allocatable,       intent(out) :: vEff0IRpw(:, :)
    complex,  allocatable,       intent(out) :: vEff0IRPwUw(:, :)
    complex,  allocatable,       intent(out) :: rho0MTsh(:, :, :, :)
    complex,  allocatable,       intent(out) :: vEff0MTsh(:, :, :, :)

    complex, allocatable                     :: GrvEff0MT(:, :, :, :)
    ! Array Arguments
    allocate( rho0IRpw(ngdp, 1) )
    allocate( vEff0IRpw(ngdp, 1) )
    allocate( vEff0IRpwUw(ngdp, 1) )
    allocate( rho0MTsh( atoms%jmtd, (atoms%lmaxd + 1)**2, atoms%nat, 1) )
    allocate( vEff0MTsh( atoms%jmtd, (atoms%lmaxd + 1)**2, atoms%nat, 1) )

    rho0IRpw(:, :) = cmplx(0., 0.)
    vEff0IRpw(:, :) = cmplx(0., 0.)
    rho0MTsh(:, :, :, :) = cmplx(0., 0.)
    vEff0MTsh(:, :, :, :) = cmplx(0., 0.)

    if (compPhon) then
      ! Conversion of the interstitial charge density
      open(108,file='000_Gvec',form='FORMATTED',action='WRITE',status='replace')
      open(109,file='000_rho_pw',form='FORMATTED',action='WRITE',status='replace')
      call convertStar2G( rho0IRst(:, 1), rho0IRpw(:, 1), stars, ngdp, gdp, cell )
      close(108)
      close(109)
      ! Conversion of the interstitial warped effective potential
      open(109,file='000_V_pw_w',form='FORMATTED',action='WRITE',status='replace')
      call convertStar2G( vEff0IRst(:, 1), vEff0IRpw(:, 1), stars, ngdp, gdp )
      close(109)
      ! Conversion of the interstitial unwarped effective potential
      open(109,file='000_V_pw',form='FORMATTED',action='WRITE',status='replace')
      call convertStar2G( vEff0IRstUw(:, 1), vEff0IRpwUw(:, 1), stars, ngdp, gdp )
      close(109)
      ! Conversion of the muffin-tin charge density
      open(109,file='000_rho_mt',form='FORMATTED',action='WRITE',status='replace')
      call convertLH2SphHarm( atoms, lathar, clnu_atom, nmem_atom, mlh_atom, rho0MTlh, rho0MTsh )
      close(109)
      ! Conversion of the muffin-tin effective potential
      open(109,file='000_V_mt',form='FORMATTED',action='WRITE',status='replace')
      call convertLH2SphHarm( atoms, lathar, clnu_atom, nmem_atom, mlh_atom, vEff0MTlh, vEff0MTsh )
      close(109)

      open(109,file='000_grad_V_mt_num',form='FORMATTED',action='WRITE',status='replace')
      call calcGrFinLH( atoms, lathar, clnu_atom, nmem_atom, mlh_atom, vEff0MTlh(:, :, :, 1), GrvEff0MT )
      close(109)
    else
      ! Conversion of the interstitial charge density
      call convertStar2G( rho0IRst(:, 1), rho0IRpw(:, 1), stars, ngdp, gdp )

      ! Conversion of the interstitial warped effective potential
      call convertStar2G( vEff0IRst(:, 1), vEff0IRpw(:, 1), stars, ngdp, gdp )

      ! Conversion of the interstitial unwarped effective potential
      call convertStar2G( vEff0IRstUw(:, 1), vEff0IRpwUw(:, 1), stars, ngdp, gdp )

      ! Conversion of the muffin-tin charge density
      call convertLH2SphHarm( atoms, lathar, clnu_atom, nmem_atom, mlh_atom, rho0MTlh, rho0MTsh )

      ! Conversion of the muffin-tin effective potential
      call convertLH2SphHarm( atoms, lathar, clnu_atom, nmem_atom, mlh_atom, vEff0MTlh, vEff0MTsh )

    end if
  end subroutine UnfoldRho0nVeff0Symmetry

  !>--------------------------------------------------------------------------------------------------------------------------------
  !> @author
  !> Christian-Roman Gerhorst
  !>
  !> @brief
  !> Determines nlotot dimension to be printed out in juPhon.log after creating k-points
  !>
  !> @details
  !> nlotot is determined so early because it has to be added to nvd and to get the number of bands for the calculation and to be
  !> changed in the fl7para file.
  !>
  !> @param[in]    atoms      : Atoms type, see types.f90.
  !> @param[in]    logUnit    : Unit number for juPhon.log.
  !>
  !>--------------------------------------------------------------------------------------------------------------------------------
  subroutine DetermineNlotot( atoms, logUnit )

    use m_types

    implicit none

    type(t_atoms), intent(inout) :: atoms

    integer,       intent(in)    :: logUnit

    integer                      :: itype ! Loops over atom types
    integer                      :: ilo   ! Loops over local orbitals of one atom type

    atoms%nlotot = 0
    do itype = 1, atoms%ntype
      do ilo = 1,atoms%nlo(itype)
        atoms%nlotot = atoms%nlotot + atoms%neq(itype) * ( 2 * atoms%llo(ilo, itype) + 1 )
      end do
    end do

    write (logUnit, '(a,1x,i6)') 'nlotot =', atoms%nlotot

  end subroutine DetermineNlotot

  !>--------------------------------------------------------------------------------------------------------------------------------
  !> @author
  !> Christian-Roman Gerhorst
  !>
  !> @brief
  !> Determing number of occupied bands.
  !>
  !> @details
  !> The idea of the routine has been adapted from Fleur. It calls the routine fermie to acquire the occupation number array w_iks.
  !> After that this array is scanned to determine the number of occupied bands.
  !>
  !> @param[in]    atoms      : Atoms type, see types.f90.
  !> @param[in]    cell       : Unit cell type, see types.f90.
  !> @param[in]    dimens     : Dimension type, see types.f90.
  !> @param[in]    input      : Input type, see types.f90.
  !> @param[in]    kpts       : K-points type, see types.f90.
  !> @param[inout] results    : Results type, see types.f90.
  !> @param[in]    logUnit    : Unit number for juPhon.log.
  !> @param[out]   nobd       : Number of occupied bands per k-point and spin
  !>
  !> @todo : check that occupied bands are determined correctly, see todo within the routine
  !>--------------------------------------------------------------------------------------------------------------------------------
  subroutine DetermineOccBands( atoms, cell, dimens, input, kpts, results, logUnit, nobd )

    use m_types
    use m_fermie
    use m_JPConstants, only : tpi, compPhon
    use mod_juPhonUtils, only : fopen, fclose

    implicit none

    ! Type parameters
    type(t_atoms),                  intent(in)    :: atoms
    type(t_cell),                   intent(in)    :: cell
    type(t_dimension),              intent(in)    :: dimens
    type(t_input),                  intent(in)    :: input
    type(t_kpts),                   intent(in)    :: kpts
    type(t_results),                intent(inout) :: results

    ! Scalar parmaeters
    integer,                        intent(in)    :: logUnit

    ! Array parameters
    integer,           allocatable, intent(out)   :: nobd(:, :)

    ! Local type variables
    type(t_jij)                                   :: jij
    type(t_noco)                                  :: noco

    ! Local scalar variables
    integer                                       :: irecl
    integer                                       :: int_len
    integer                                       :: iTemp
    integer                                       :: iunit = 203
    real                                          :: ef
    real                                          :: seigscv
    real                                          :: ts
    integer                                       :: ikpt
    integer                                       :: noccbd
    integer                                       :: nband

    ! Local array variables
    integer                                       :: nkptOthFormat(1)
    real                                          :: zelecOthFormat(1)

    ! This is required due to fermie which assumes nkpt and zelec to be dependent of energy windows
    nkptOthFormat(1) = kpts%nkpt
    zelecOthFormat(1) = input%zelec

    ! Allocation of variables to execute fermie
    allocate( results%w_iks(dimens%neigd, kpts%nkpt, dimens%jspd) )
    allocate( jij%eig_l(dimens%neigd + 6, jij%nkpt_l) )

    ! Calculating record length for reading the eig file
    irecl = atoms%ntype * (atoms%lmaxd + 1 + atoms%nlod) + 2
#ifdef CPP_INVERSION
    irecl = 8 * (irecl + 4 + dimens%neigd * (dimens%nbasfcn + 1))
#else
    irecl = 8 * (irecl + 4 + dimens%neigd * (2 * dimens%nbasfcn + 1))
#endif
    irecl = irecl + 4 * atoms%nat * atoms%nlod * (2 * atoms%llod + 1)

    iTemp = size(transfer(tpi,(/1/)))
    if (iTemp == 1) then
      int_len = 8                 ! iTemp = 1 if real*8
    else if (iTemp == 2) then
      int_len = 4                 ! iTemp = 2 if real*4
    else
      NOstopNO'read_eig: strange length of integers'
    endif
    irecl = irecl + int_len * (3 + 3 * dimens%nvd)
    !Open eigfile for fermie
    call fopen(iunit, name='eig', status='old', action='read', form='unformatted', access='direct', recl=irecl)

    ! Determine w_iks
    call fermie( dimens%neigd, kpts%nkptd, dimens%jspd, 1, atoms%ntype, atoms%lmaxd, atoms%nlod, input%jspins, 1, atoms%ntype, &
      & nkptOthFormat, iunit, irecl, .false., .false., input%gauss, input%delgau, zelecOthFormat, input%tkb, input%tria,           &
      & noco%l_noco, noco%l_ss, ef, seigscv, ts, results%w_iks, noco%qss, jij%l_J, jij%l_disp, cell%bmat, jij%nkpt_l, jij%eig_l, &
      & 0, 1)
    call fclose(iunit)


    ! Count number of occupied bands. Having found the first occupied band the next k-point is scanned.

    if (compPhon) then
      open(111,file='000_occs',form='FORMATTED',action='WRITE',status='replace')
    end if

    allocate ( nobd(kpts%nkpt, input%jspins) )
    nobd = 0
    do ikpt = 1, kpts%nkpt
      noccbd = 0
      do nband = 1, dimens%neigd ! todo not until ne(ikpt?)

        if (compPhon) then
          write(111,*) ikpt, nband
          write(111,*) results%w_iks(nband,ikpt,1)
        end if

        if ( results%w_iks( nband, ikpt, 1 )  >= 1.e-8 ) then
          noccbd = noccbd + 1
        else
          exit
        end if
      end do
      nobd(ikpt, 1) = noccbd
    end do

    if (compPhon) then
      close(111)
    end if

  end subroutine DetermineOccBands

  !>--------------------------------------------------------------------------------------------------------------------------------
  !> @author
  !> Christian-Roman Gerhorst and Markus Betzinger
  !>
  !> @brief
  !> Acquire input which is needed for the minimal run of juPhon, the k-point generation mode.
  !>
  !> @details
  !> For generating a Gamma-point centered k-point set not all information of the fleur init routine is needed. Hence, this routine
  !> constitutes a minimal version of the fleur init routine.
  !>
  !> @param[out] sym        : Symmetries type, see types.f90.
  !> @param[out] cell       : Unit cell type, see types.f90.
  !> @param[in]  logUnit    : Unit number for juPhon.log.
  !>
  !>--------------------------------------------------------------------------------------------------------------------------------
  subroutine PrepareKpointGeneration( sym, cell, logFUnit )

    use m_types
    use mod_juPhonUtils
    use m_JPConstants, only : tpi
    use m_juDFT_NOstopNO, only : juDFT_error
    use m_inv3
    use m_rwsymfile

    implicit none

    ! Scalar parameter
    integer,                        intent(in)  :: logFUnit

    ! Type parameters
    type(t_sym),                    intent(out) :: sym
    type(t_cell),                   intent(out) :: cell

    ! Scalar variables
    integer                                     :: FLEURinpUnit=200   ! global unit for FLEUR inp file
    real                                        :: latVecScal         ! Scaling of the lattice vectors
    character(len=1)                            :: rw='rR'            ! switch whehter sym.out file is read or is written in
    character(len=7)                            :: symfn='sym.out'    ! file name of sym.out file
    integer                                     :: nopd=48            ! maximal number of possible symmetries
    integer                                     :: symfhUnit=201      ! global unit for sym.out file
    integer                                     :: isym, jsym, ii, jj ! loop variables

    ! Array variables
    real,             allocatable               :: rtau(:,:)          ! Temporary array to store sym%tau
    integer,          allocatable               :: rot(:,:,:)         ! Temporary array to store sym%mrot

    ! Create reciprocal Bravais matrix
    call Inv3( cell%amat, cell%bmat, cell%omtil )
    cell%bmat(:, :) = tpi * cell%bmat(:, :)
    write ( logFUnit, '(a)' ) 'Bravais matrix in real lattice read in from inp:'
    do ii = 1, 3
      write ( logFUnit, '(5x,3(f10.6))' ) ( cell%amat(ii,jj), jj=1, 3 )
    enddo
    write ( logFUnit, '(a)' )
    write ( logFUnit, '(a)' ) 'Transposed Bravais matrix in reciprocal lattice:'
    do ii = 1, 3
      write ( logFUnit, '(5x,3(f10.6))' ) ( cell%bmat(ii,jj), jj=1, 3 )
    enddo

    ! Find symmetry transformation
    allocate( sym%mrot(3, 3, 48), sym%tau(3, 48) )
    call Rw_symfile( rw, symfhUnit, symfn, nopd, cell%bmat, sym%mrot, sym%tau, sym%nop, sym%nop2, sym%symor )
    write( logFUnit, '(a,x,i2,x,a)' )
    if ( sym%nop == 1 ) then
      write( logFUnit, '(a,x,i1,x,a)' ) 'Extracted', sym%nop, 'operation (identity) from sym.out file.'
    else
      write( logFUnit, '(a,x,i2,x,a)' ) 'Extracted', sym%nop, 'operations from sym.out file.'
    end if

    ! Check symmetries
    allocate( rot(3, 3, sym%nop), rtau(3, sym%nop) )
    do isym=1,sym%nop
      rot(:, :, isym) = sym%mrot(:, :, isym)
      rtau(  :, isym) = sym%tau(:, isym)
    end do
    if(any( rot(:, :, 1) - reshape((/1,0,0,0,1,0,0,0,1/), (/3,3/)) /= 0) ) then
      call juDFT_error( 'First symmetry operation is not the identity.', calledby='prepareKpointGeneration', &
        & hint='Review FLEUR calculation' )
    endif

    !Generate and check inverse symmetry transformations
    allocate( sym%invtab(sym%nop) )
    sym%invtab = 0
    do isym = 1,sym%nop
      do jsym = 1,sym%nop
        if( all( matmul( rot(:,:,isym), rot(:,:,jsym) ) == reshape((/1,0,0,0,1,0,0,0,1/), (/3,3/))) &
          &.and. all( modulo( matmul( rot(:,:,isym), rtau(:,jsym)) + rtau(:,isym), 1.0) < 1d-10 ) ) then
          if( sym%invtab(isym) /= 0 ) then
            call juDFT_error( 'inverse operation already defined.', calledby='prepareKpointGeneration', &
              & hint='Review FLEUR calculation' )
          endif
          sym%invtab(isym) = jsym
        end if
      end do
      if( sym%invtab(isym) == 0 ) then
        call juDFT_error( 'inverse operation not found.', calledby='prepareKpointGeneration', hint='Review FLEUR calculation' )
      endif
    end do
    write( logFUnit, * )
    call GenerateRecSymOp( sym )

  end subroutine PrepareKpointGeneration

  !>--------------------------------------------------------------------------------------------------------------------------------
  !> @author
  !> Christian-Roman Gerhorst and Markus Betzinger
  !>
  !> @brief
  !> Parses the eig file which has been generated by a preceeding FLEUR calculation.
  !>
  !> @details
  !> Several unperturbed quantities which have been already acquired by a FLEUR calculation are stored into the eig file. This routine
  !> recycles the stored quantities by parsing the eig file.
  !>
  !> @note
  !> For this routine to work there must be an eig file.
  !> @note GbasVec, ilst, kveclo, z are tested in test of wavefunction continuity, among others. Eig should be right as other quantities
  !>
  !> @param[in] atoms        : Atoms type, see types.f90.
  !> @param[in] dimens       : Dimension type, see types.f90.
  !> @param[in] kpts         : K-points type, see types.f90.
  !> @param[in] input        : Input type, see types.f90.
  !> @param[out] enpara      : Energy parameter type, see types.f90.
  !> @param[out] eig         : Contains Kohn\--Sham eigenvalues.
  !> @param[out] kveclo      : Basis G-vectors of local orbitals.
  !> @param[out] nv          : Number of LAPW G-basis vectors for given k-point.
  !> @param[out] ilst        : For various k-points G-basis vectors occur more than once, thus they are only stored once in juPhon. This
  !>                           pointer array contains the right index for GbasVec array to "unfold" G-basis vectors again.
  !> @param[out] GbasVec     : G-basis vectors
  !> @param[out] ne          : Number of eigenvalues per k-point.
  !> @param[out] El          : Contains LAPW and LOs energy parameters.
  !> @param[in]  logUnit     : Unit number for juPhon.log.
  !>
  !> @todo kveclo is treated seperately maybe this can be changed in future
  !>  read out with certain record_length are tested so eig lying between these quantities in eig file is as well.
  !>
  !>--------------------------------------------------------------------------------------------------------------------------------
  subroutine Read_eig( atoms, dimens, kpts, input, enpara, cell, eig, kveclo, z, nv, ilst, GbasVec, ne, El, logUnit )

#include "cppmacro.h"

    use m_JPConstants
    use mod_juPhonUtils
    use m_types
    use m_juDFT_NOstopNO, only : juDFT_error, juDFT_warn
    use m_jpLog, only : LogReadFromEigFile

    implicit none

    ! Type parameters
    type(t_atoms),         intent(in)  :: atoms
    type(t_dimension),     intent(in)  :: dimens
    type(t_kpts),          intent(in)  :: kpts
    type(t_input),         intent(in)  :: input
    type(t_enpara),        intent(out) :: enpara
    type(t_cell),          intent(in)  :: cell

    ! Scalar parameters
    integer,               intent(in)  :: logUnit

    ! Array parameters
    real,     allocatable, intent(out) :: eig(:,:,:)
    integer,  allocatable, intent(out) :: kveclo(:,:)
    MCOMPLEX, allocatable, intent(out) :: z(:,:,:,:)
    integer,  allocatable, intent(out) :: nv(:, :)
    integer,  allocatable, intent(out) :: ilst(:, :, :)
    integer,  allocatable, intent(out) :: GbasVec(:, :)
    integer,  allocatable, intent(out) :: ne(:)
    real,     allocatable, intent(out) :: El(:, :, :, :)


    ! Scalar variables
    integer                             :: ii, ikpt, isp, n1, n2   ! loop variables
    integer                             :: intDump                 ! integer temporary variable
    integer                             :: irecl, irec             ! variables for record length
    integer                             :: int_len                 ! integer length to determine correct record length
    integer                             :: ierror                  ! error nubmer
    integer                             :: iunit = 203             ! unit for eig file
    real                                :: rdum                    ! real temporary variable
    logical                             :: eigExists               ! logical temporary variable

    ! Array variables
    real                                :: bkpt(3), kvec(3)                 ! temporary k-point read in from eig
    integer ,allocatable                :: k1(:)                   ! temporary x coordinates of basis vectors for a k-point
    integer ,allocatable                :: k2(:)                   ! temporary y coordinates of basis vectors for a k-point
    integer ,allocatable                :: k3(:)                   ! temporary z coordinates of basis vectors for a k-point
    integer                             :: iG                      ! loop index
    integer                             :: indexG                  ! accumulated index
    integer, allocatable                :: GbasVec_eig(:, :, :, :) ! temporary storage for reading basis vectors G from eig
    integer, allocatable                :: GbasVec_temp(:, :)      ! Helping array for compressed storage of basis vectors G
    character(len=1024)                 :: errorMessage            ! String for error messages
    integer                             :: ifind                   ! loop Variable
    logical                             :: addG                    ! Helping logical for compressed storage of basis vectors G
    integer                             :: kindIndicat             ! Helping integer for determing right record length

    ! Find out whether eig file exists
    inquire(file='eig', exist=eigExists)
    if( .not.eigExists ) then
      call juDFT_error('File eig does not exist', calledby='Read_eig', hint='Perform new Fleur calculation.')
    end if

    ! calculate record length of the eig-file to read it out properly
    irecl = atoms%ntype * (atoms%lmaxd + 1 + atoms%nlod) + 2
#ifdef CPP_INVERSION
    irecl = 8 * (irecl + 4 + dimens%neigd * (dimens%nbasfcn + 1))
#else
    irecl = 8 * (irecl + 4 + dimens%neigd * (2 * dimens%nbasfcn + 1))
#endif
    irecl = irecl + 4 * atoms%nat * atoms%nlod * (2 * atoms%llod + 1)

    kindIndicat = size(transfer(tpi,(/1/)))
    if ( kindIndicat == 1 ) then
      int_len = 8                 ! idum = 1 if real*8
    else if ( kindIndicat == 2 ) then
      int_len = 4                 ! idum = 2 if real*4
    else
      NOstopNO'read_eig: strange length of integers'
    endif
    irecl = irecl + int_len * (3 + 3 * dimens%nvd)


    call fopen(iunit, name='eig', status='old', action='read', form='unformatted', access='direct', recl=irecl)

    ! Allocate quantities to be read out from eig file
    allocate( nv(input%jspins, kpts%nkpt), ne(kpts%nkpt), k1(dimens%nvd), k2(dimens%nvd), k3(dimens%nvd) )

    allocate( z(dimens%nbasfcn, dimens%neigd, kpts%nkpt, input%jspins), stat=ierror )
    if( ierror /=  0 ) then
      call juDFT_error( 'Error in allocation of eigenvectors z', calledby='Read_eig', hint='Check for having enough storage' )
    end if

    allocate( eig(dimens%neigd, kpts%nkpt, input%jspins), stat=ierror )
    if( ierror /=  0 ) then
      call juDFT_error( 'Error in allocation of eigen energies eig', calledby='Read_eig', hint='Check for having enough storage' )
    end if

    allocate( enpara%el0(0:atoms%lmaxd, atoms%ntype, input%jspins), enpara%ello0(atoms%nlod, atoms%ntype, input%jspins), &
      & stat=ierror )
    if( ierror /=  0 ) then
      call juDFT_error( 'Error in allocation of energy parameters el0 ello0', calledby='Read_eig', &
        & hint='Check for having enough storage' )
    end if

    allocate( kveclo(atoms%nlotot, kpts%nkpt), stat=ierror )
    if( ierror /=  0 ) then
      call juDFT_error( 'Error in allocation of G-basis vectors for LOs kveclo', calledby='Read_eig', &
        & hint='Check for having enough storage' )
    end if

    allocate( ilst(dimens%nvd, kpts%nkpt, input%jspins), GbasVec_eig(3, dimens%nvd, kpts%nkpt, input%jspins), &
      & GbasVec_temp(3, 0), GbasVec(3, 0) )

    ! Actually reading out eig file line after line.
    indexG = 0
    GbasVec_eig = 0

    if (compPhon) then
      open(109,file='000_evals',form='FORMATTED',action='WRITE',status='replace')
      open(110,file='000_kgs',form='FORMATTED',action='WRITE',status='replace')
    end if

    do isp = 1, input%jspins
      do ikpt = 1, kpts%nkpt
        irec = kpts%nkpt * (isp - 1) + ikpt
        read(iunit, rec=irec) enpara%el0(:,:,isp), (rdum,ii=1,2), enpara%ello0(:,:,isp), bkpt, rdum, ne(ikpt), nv(isp, ikpt), &
          & intDump, (eig(ii, ikpt, isp), ii=1, dimens%neigd), k1(:), k2(:), k3(:), kveclo(:, ikpt), z(:, :, ikpt, isp)

        if (compPhon) then
          do ii=1, dimens%neigd
            write(109,*) ikpt, ii
            write(109,*) eig(ii,ikpt,1)
          end do
        end if

        GbasVec_eig(1, :nv(1, ikpt), ikpt, isp) = k1(:nv(1, ikpt))
        GbasVec_eig(2, :nv(1, ikpt), ikpt, isp) = k2(:nv(1, ikpt))
        GbasVec_eig(3, :nv(1, ikpt), ikpt, isp) = k3(:nv(1, ikpt))

        if (compPhon) then
          do n1=1, nv(1, ikpt)
            kvec(:) = matmul(cell%bmat(1:3, 1:3), GbasVec_eig(:, n1, ikpt, isp) - bkpt(:))
            write(110,*) ikpt, n1
            write(110,*) kvec(1), kvec(2), kvec(3)
          end do
        end if

        ! Consistency check of k-points
        if(maxval(abs(bkpt - kpts%bk(:,ikpt))) >=  1E-8) then !idum, i=1, 3*dimens%nvd G-vektoren
          write (errorMessage, '("k-point",1x,"(",3(f15.8),")",1x,"from eig file differs from internal k-point",1x,"(",3(f15.8),")"&
            &,1x,"with index",1x,i3)') bkpt, kpts%bk(:, ikpt), ikpt
          call juDFT_error(errorMessage, calledby='Read_eig', hint='Start complete calculation procedure from beginning')
        end if

      end do ! ikpt
    end do ! isp

    if (compPhon) then
      close(109)
      close(110)
    end if

    deallocate( k1, k2, k3 )

    if ( dimens%neigd /= dimens%nvd + atoms%nlotot ) then
      call juDFT_warn( 'Number of eigenvalues is not equals number of basis functions. You might miss relevant occupied bands!', &
        &calledby='Read_eig', hint='Set neigd = nvd + nlotot within fl7para file' )
    end if

    ! Store basis vectors G in a compressed way and not redundantly, as they occur multiple times for different k-points. The array
    ! ilst stores the index where the basis vector is stored if it has occured for a different k-point or spin yet.
    indexG = 0
    ! Initialize ilst generally to -1 to indicate unvalid value
    ilst = -1
    do isp = 1, input%jspins
      do ikpt = 1, kpts%nkpt
        do iG = 1, nv(isp, ikpt)
          if ( uBound( GbasVec, 2 ) == 0 ) then
            addG = .true.
          else
            do ifind =  1, uBound( GbasVec, 2 )

              if ( .not. all( abs( GbasVec(:, ifind) - GbasVec_eig(:, iG, ikpt, isp) ) <= 1E-8 ) ) then
                addG = .true.
              else
                ilst(iG, ikpt, isp) = ifind
                addG = .false.
                exit
              endif
            enddo ! ifind
          endif

          if( addG ) then
            if ( uBound( GbasVec, 2 ) == 0 ) then
              deallocate( GbasVec, GbasVec_temp )
              allocate( GbasVec(3, 1), GbasVec_temp(3, 1) )
            else
              GbasVec_temp = GbasVec
              deallocate( GbasVec )
              allocate( GbasVec(3, indexG + 1) )
              GbasVec(:, :indexG) = GbasVec_temp
            endif

            indexG = indexG + 1
            ilst(iG, ikpt, isp) = indexG
            GbasVec(1, indexG) = GbasVec_eig(1, iG, ikpt, isp)
            GbasVec(2, indexG) = GbasVec_eig(2, iG, ikpt, isp)
            GbasVec(3, indexG) = GbasVec_eig(3, iG, ikpt, isp)
            if( uBound( GbasVec, 2 ) /= 0 ) then
              ! resize GbasVec_temp to current size of GbasVec
              deallocate( GbasVec_temp )
              allocate( GbasVec_temp(3, indexG) )
            end if
          end if
        end do ! iG
      end do ! ikpt
    end do ! isp

    !if (compPhon) then
      !open(110,file='000_z0',form='FORMATTED',action='WRITE',status='replace')
      open(111,file='000_kqpts',form='FORMATTED',action='WRITE',status='replace')
      do ikpt = 1, kpts%nkpt
        write(111,*) ikpt
        write(111,*) kpts%bk(1, ikpt), kpts%bk(2, ikpt), kpts%bk(3, ikpt)
      !  do n2=1, size(z,2)
      !    do n1=1, nv(1, ikpt)
      !      write(110,*) ikpt, n1, n2
      !      write(110,*) real(z(n1, n2, ikpt, 1)), aimag(z(n1, n2, ikpt, 1))
      !    end do
      !  end do
      end do
      !close(110)
      close(111)
      !NOstopNO
    !end if

    call fclose(iunit)

    call uniteEnergyParameters( atoms, enpara, dimens, El )

    call LogReadFromEigFile(atoms, kpts, dimens, El, nv, ne, logUnit)

  end subroutine Read_eig

  !>-------------------------------------------------------------------------------------------------------------------------------------
  !> @author
  !> Christian-Roman Gerhorst and Markus Betzinger
  !>
  !> @brief
  !> Reads in potential files of unperturbed potentials generated by FLEUR
  !>
  !> @details
  !>
  !> @note
  !> For this routine to work there must be the three potential files vpw, pottot_unwarped and pottot.
  !>
  !> @param[in] atoms       : Atoms type, see types.f90.
  !> @param[in] stars       : Stars type, see types.f90.
  !> @param[in] dimens      : Dimension type, see types.f90.
  !> @param[in] lathar      : Lattice harmonics type, see types.f90.
  !> @param[in] input       : Input type, see types.f90.
  !> @param[in] sym         : Symmetries type, see types.f90.
  !> @param[in] vacuum      : Vacuum calculations type, see types.f90
  !> @param[in] oneD        : Type 1D-calculations
  !> @param[out] Veff0      : Type containing unperturbed output potentials of Fleur, see types.f90.
  !> @param[out] vrCoul     : Unperturbed radial coefficients of the converged Coulomb potential parsed from FLEUR.
  !> @param[out] vpwCoul_uw : Star coefficients of the unwarped (no convolution with step function) and converged interstitial
  !>                          Coulomb potential parsed from Fleur.
  !> @param[in]  logUnit    : Unit number for juPhon.log.
  !>
  !> @todo : Potentials for Debugging reasons should be read in within the debugging methods, after that this routine have to be
  !>         cleaned
  !>--------------------------------------------------------------------------------------------------------------------------------
  subroutine ReadPotFleur( atoms, stars, dimens, lathar, input, sym, vacuum,   V0Fleur, vrCoul, vpwCoul_uw, logUnit )

    use m_types
    use m_loddop
    use mod_juPhonUtils
    use m_JPConstants, only: fpi, compPhon

    implicit none

    ! Type parameter
    type(t_atoms),                  intent(in)  :: atoms
    type(t_stars),                  intent(in)  :: stars
    type(t_dimension),              intent(in)  :: dimens
    type(t_sphhar),                 intent(in)  :: lathar
    type(t_input),                  intent(in)  :: input
    type(t_sym),                    intent(in)  :: sym
    type(t_vacuum),                 intent(in)  :: vacuum
     
    type(t_potential),              intent(out) :: V0Fleur

    ! Array parameter
    real,              allocatable, intent(out) :: vrCoul(:, :, :, :)
    complex,           allocatable, intent(out) :: vpwCoul_uw(:, :)

    ! Scalar parameter
    integer,                        intent(in)  :: logUnit

    ! Scalar variables
    integer(4)                                  :: iunit = 204        ! unit number for potential files
    integer                                     :: isp                ! loop variable
    integer                                     :: itype              ! loop variable
    integer                                     :: imesh              ! loop variable
    integer                                     :: iter               ! irrelevant for juPhon

    ! Array variables
    character(8)                                :: dop, iop, name(10) ! irrlevant for juPhon

    ! todo review this code after it is clear what this routine should finally read in (only standard effective potential or more)
    allocate(V0Fleur%vpw_uw(stars%ng3, input%jspins), V0Fleur%vpw(stars%ng3,input%jspins), &
      & V0Fleur%vzxy(vacuum%nmzxyd, stars%ng2-1,2, input%jspins), V0Fleur%vz(vacuum%nmzd, 2, input%jspins), &
      & V0Fleur%vr(atoms%jmtd, 0:lathar%nlhd, atoms%ntype, input%jspins),V0Fleur%vr0(atoms%jmtd, atoms%ntype, input%jspins))
    allocate(vpwCoul_uw(stars%ng3, input%jspins), vrCoul(atoms%jmtd, 0:lathar%nlhd, atoms%ntype, input%jspins))


    call fopen(iunit, name='pottot_unwarped', status='old', action='read', form='unformatted')

    rewind (iunit)

    call loddop(input%jspins, stars%ng3, stars%ng2, vacuum%nmzxyd, vacuum%nmzd, atoms%jmtd, lathar%nlhd, atoms%ntype,           &
      & input%jspins, stars%ng3, stars%ng2, vacuum%nvac, atoms%ntype, sym%invs, sym%invs2, input%film, lathar%nlh, atoms%jri,   &
      & lathar%ntypsd, atoms%ntypsy, iunit, atoms%nat, atoms%neq, iop, dop, iter, V0Fleur%vr, V0Fleur%vpw_uw, V0Fleur%vz,          &
      & V0Fleur%vzxy, name)
    call fclose(iunit)

    call fopen(iunit, name='pottot', status='old', action='read', form='unformatted')

    rewind (iunit)

    call loddop(input%jspins, stars%ng3, stars%ng2, vacuum%nmzxyd, vacuum%nmzd, atoms%jmtd, lathar%nlhd, atoms%ntype,           &
      & input%jspins, stars%ng3, stars%ng2, vacuum%nvac, atoms%ntype, sym%invs, sym%invs2, input%film, lathar%nlh, atoms%jri,   &
      & lathar%ntypsd, atoms%ntypsy, iunit, atoms%nat, atoms%neq, iop, dop, iter,V0Fleur%vr, V0Fleur%vpw, V0Fleur%vz, V0Fleur%vzxy,&
      & name)

    ! store spherical part in V0Fleur%vr0. This is needed for correctly calculating the radial basis functions.
    V0Fleur%vr0 = V0Fleur%vr(:,0,:,:)

    ! remove r/sqrt(4*pi) from V0Fleur%vr

    if (compPhon) then
      open(109,file='000_r',form='FORMATTED',action='WRITE',status='replace')
    end if

    do isp = 1,input%jspins
      do itype = 1,atoms%ntype
        do imesh = 1,atoms%jri(itype)
          V0Fleur%vr(imesh,0,itype,isp) = V0Fleur%vr(imesh,0,itype,isp)*sqrt(fpi)/atoms%rmsh(imesh,itype)

          if (compPhon) then
            write(109,*) imesh, itype
            write(109,*) atoms%rmsh(imesh,itype)
          end if

        end do
      end do
    end do

    if (compPhon) then
      close(109)
    end if

    call fclose(iunit)
    ! end todo

    write ( logUnit, * )
    write ( logUnit, * )
    write ( logUnit, '(a)' ) 'Read in converged and unperturbed Coulomb and effective potentials for interstitial and muffin-tin re&
      &gions from files pottot_uw and potcoul_uw'
    write ( logUnit, '(a)' ) '=====================================================================================================&
      &========================================='
    write ( logUnit, * )

  end subroutine ReadPotFleur

  !>--------------------------------------------------------------------------------------------------------------------------------
  !> @author
  !> Christian-Roman Gerhorst and Markus Betzinger
  !>
  !> @brief
  !> Read in the unperturbed density file generated by FLEUR.
  !>
  !> @note
  !> For this routine to work, there must be a density file called cdn1.
  !> @details
  !>
  !> @param[in] atoms      : Atoms type, see types.f90.
  !> @param[in] sym        : Symmetries type, see types.f90.
  !> @param[in] lathar     : Lattice harmonics type, see types.f90.
  !> @param[in] dimens     : Dimension type, see types.f90.
  !> @param[in] stars      : Stars type, see types.f90.
  !> @param[in] input      : Input type, see types.f90.
  !> @param[in] cell       : Unit cell type, see types.f90.
  !> @param[in] vacuum     : Vacuum calculations type, see types.f90.
  !> @param[in] oneD       : Type 1D-calculationst type, see types.f90.
  !> @param[in]  logUnit   : Unit number for juPhon.log.
  !> @param[out] rho0IR    : Star coefficients of the unperturbed and converged interstitial density parsed from Fleur.
  !> @param[out] rho0MT    : Radial coefficients of the unperturbed and converged muffin-tin densities parsed from Fleur.
  !>
  !> @todo : Put rhos into general density type which is similiar to the potential type
  !>--------------------------------------------------------------------------------------------------------------------------------
  subroutine ReadDensFleur( atoms, sym, lathar, dimens, stars, input, cell, vacuum,   logUnit, rho0IR, rho0MT )

    use m_types
    use m_loddop
    use mod_juPhonUtils
    use m_cdntot
    use m_jpLog, only : LogReadcdn1

    implicit none

    ! Type parameters
    type(t_atoms),                  intent(in)  :: atoms
    type(t_sym),                    intent(in)  :: sym
    type(t_sphhar),                 intent(in)  :: lathar
    type(t_dimension),              intent(in)  :: dimens
    type(t_stars),                  intent(in)  :: stars
    type(t_input),                  intent(in)  :: input
    type(t_cell),                   intent(in)  :: cell
    type(t_vacuum),                 intent(in)  :: vacuum
     

    ! Scalar parameters
    integer,                        intent(in)  :: logUnit
    ! Array parameters
    !todo write them into general potential and density type!
    complex,           allocatable, intent(out) :: rho0IR(:,:)
    real,              allocatable, intent(out) :: rho0MT(:,:,:,:)
    !todo end

    ! Scalar variables
    integer(4)                                  :: iunit = 205        ! unit number for density
    integer                                     :: iter               ! irrelevant for juPhon
    integer                                     :: ispin              ! loop variable
    integer                                     :: ilh                ! loop variable
    integer                                     :: imesh              ! loop variable
    integer                                     :: itype              ! loop variable
    real                                        :: qtotIR             ! total charge in interstitial
    real                                        :: qtot               ! total charge in unit cell

    ! Array variables
    character(8)                                :: dop, iop, name(10) ! irrelevant for juPhon calculation
    real,              allocatable              :: nz(:,:,:)          ! irrelevant for juPhon calculation
    complex,           allocatable              :: nzxy(:,:,:,:)      ! irrelevant for juPhon calcualtion
    real                                        :: qMTs(atoms%ntype)  ! charge of muffin-tins

    allocate( rho0IR(stars%ng3,input%jspins), nzxy(vacuum%nmzxyd, stars%ng2-1,2, input%jspins), &
      & nz(vacuum%nmzd, 2, input%jspins), rho0MT(atoms%jmtd, 0:lathar%nlhd, atoms%ntype, input%jspins) )


    call Fopen(iunit, name='cdn1', status='old', action='read', form='unformatted')

    rewind (iunit)

    call loddop( input%jspins, stars%ng3, stars%ng2, vacuum%nmzxyd, vacuum%nmzd, atoms%jmtd, lathar%nlhd, atoms%ntype,          &
      & input%jspins, stars%ng3, stars%ng2, vacuum%nvac, atoms%ntype, sym%invs, sym%invs2, input%film, lathar%nlh, atoms%jri,   &
      & lathar%ntypsd, atoms%ntypsy, iunit, atoms%natd, atoms%neq, iop, dop, iter, rho0MT, rho0IR, nz, nzxy, name )

    call Fclose(iunit)

    ! Calculates charge in IR and complete unit cell
    call cdntot( stars, atoms, sym, vacuum, input, cell,   rho0IR, rho0MT, qtot, qtotIR, qMTs )

    call LogReadcdn1( atoms, qtotIR, qMTs, qtot, logUnit )

    ! Divide out factor r^2 so that rho0MT is only rho0MT and not r^2 rho0MT as it is done in Fleur. The factor sqrt(4pi) for the
    ! zeroth component was already dividied out when constructing rho0MT in cdnmt routine in Fleur.
    do ispin = 1, input%jspins
      do itype = 1, atoms%ntype
        do ilh = 0, lathar%nlhd
          do imesh = 1, atoms%jri(itype)
            rho0MT(imesh, ilh, itype, ispin) = rho0MT(imesh, ilh, itype, ispin) / atoms%rmsh(imesh, itype)                         &
              & / atoms%rmsh(imesh, itype)
        end do
      end do
    end do
  end do

  end subroutine ReadDensFleur

  !>--------------------------------------------------------------------------------------------------------------------------------
  !> @author
  !> Christian-Roman Gerhorst and Markus Betzinger
  !>
  !> @brief
  !> Generates radial solutions u, udot, u_LOs and overlap intergrals of them.
  !>
  !> @details
  !>
  !> @param[in] atoms      : Atoms type, see types.f90.
  !> @param[in] input      : Input type, see types.f90.
  !> @param[in] dimens     : Dimension type, see types.f90.
  !> @param[in] enpara     : Energy parameter type, see types.f90.
  !> @param[in] Veff0      : Type containing unperturbed output potentials of Fleur, see types.f90.
  !> @param[in] mpi        : MPI type, see types.f90
  !> @param[out] usdus     : Type containing quantities consisting of the radial solutions, see types.f90.
  !> @param[in]  logUnit   : Unit number for juPhon.log.
  !> @param[out] rbas1     : Large components of radial solution, its energy derivative and u_LO
  !> @param[out] rbas2     : Small components of radial solution, its energy derivative and u_LO
  !> @param[out] nRadFun   : Number of radial functions per orbital quantum number l and atom type.
  !> @param[out] iloTable  : Number of local orbital if orbital quantum number l, atom type and index p > 2 from nRadFun is given.
  !> @param[out] uuilon    : overlap integral between the radial functions of the integral (multiplied by ulo_der) of a local
  !>                         oribtal and the flapw radial function with the same l
  !> @param[out] duilon    : overlap integral between the radial functions of the integral of a local orbital and the energy
  !>                         derivative of the flapw radial function with the same l
  !> @param[out] ulouilopn : overlap integral between the radial functions of the integral of a local orbital and another local
  !>                         orbital with the same l.
  !> @param[out] ilo2p     : mapping array giving the p value for given number of LO and itype
  !>
  !>--------------------------------------------------------------------------------------------------------------------------------
  subroutine GenRadSolnInt( atoms, input, dimens, enpara, V0sFleur, mpi, usdus, logUnit, rbas1, rbas2, nRadFun, iloTable, uuilon,  &
      & duilon, ulouilopn, ilo2p )

    use m_radfun
    use m_radflo
    use m_types
    use m_jpLog, only : LogRadSolutions
    use m_jpConstants, only : compPhon

    implicit none

    ! Type parameters
    type(t_atoms),                  intent(in)  :: atoms
    type(t_input),                  intent(in)  :: input
    type(t_dimension),              intent(in)  :: dimens
    type(t_enpara),                 intent(in)  :: enpara
    type(t_potential),              intent(in)  :: V0sFleur
    type(t_mpi),                    intent(in)  :: mpi
    type(t_usdus),                  intent(out) :: usdus

    ! Scalar paramter
    integer,                        intent(in)  :: logUnit

    ! Array parameters
    real,              allocatable, intent(out) :: rbas1(:,:,:,:,:)
    real,              allocatable, intent(out) :: rbas2(:,:,:,:,:)
    integer,           allocatable, intent(out) :: nRadFun(:,:)
    integer,           allocatable, intent(out) :: iloTable(:, :, :)
    real,              allocatable, intent(out) :: uuilon(:, :)
    real,              allocatable, intent(out) :: duilon(:, :)
    real,              allocatable, intent(out) :: ulouilopn(:, :, :)
    integer,           allocatable, intent(out) :: ilo2p(:, :)

    ! Scalar variables
    real                                        :: wronk             ! irrelevant for juPhon
    integer                                     :: nodeu             ! irrelevant for juPhon
    integer                                     :: noded             ! irrelevant for juPhon
    integer                                     :: oqn_l             ! loop variable
    integer                                     :: itype             ! loop variable
    integer                                     :: isp               ! loop variable
    integer                                     :: igrid             ! loop variable
    integer                                     :: p                 ! loop variable
    integer                                     :: ilo               ! loop variable

    ! Array variables
    ! us        : radial solution at muffin-tin radius
    ! uds       : energy derivative of radial solution u at muffin-tin radius
    ! dus       : radial derivative of radial solution at muffin-tin radius
    ! duds      : energy derivative of radial derivative of radial solution u at muffin-tin radius
    ! ddn       : overlap integral of udots
    ! ulos      : the value of the radial function of a local orbital at the muffin tin radius
    ! dulos     : the value of the radial derivative of the radial function of a local orbital at the muffin tin radius
    ! uulon     : overlap integral between the radial functions of a local obital and the flapw radial function with the same l
    ! dulon     : overlap integral between the radial functions of a local obital and the energy derivative of the flapw radial
    !              function with the same l
    ! uloulopn  : overlap integral between the radial functions of two different local orbitals only needed if they have the same l
    ! f         : stores radial solution u coming out of radfun
    ! g         : stores energy derivative of radial solution udot coming out of radfun
    ! flo       : stores radial solution of local orbital coming out of radflo
    !
    real,          allocatable                  :: us(:, :)
    real,          allocatable                  :: uds(:, :)
    real,          allocatable                  :: dus(:, :)
    real,          allocatable                  :: duds(:, :)
    real,          allocatable                  :: ddn(:, :)
    real,          allocatable                  :: ulos(:, :)
    real,          allocatable                  :: dulos(:, :)
    real,          allocatable                  :: uulon(:, :)
    real,          allocatable                  :: dulon(:, :)
    real,          allocatable                  :: uloulopn(:, :, :)
    real,          allocatable                  :: f(:, :, :)
    real,          allocatable                  :: g(:, :, :)
    real,          allocatable                  :: flo(:, :, :)


    allocate( uuilon(atoms%nlod, atoms%ntype) )
    allocate( duilon(atoms%nlod, atoms%ntype) )
    allocate( ulouilopn(atoms%nlod, atoms%nlod, atoms%ntype) )
    allocate( f(atoms%jmtd, 2, 0:atoms%lmaxd), g(atoms%jmtd, 2, 0:atoms%lmaxd) ) ! second index runs over vector-component
    allocate( flo(atoms%jmtd, 2, atoms%nlod) )                                   ! second index runs over vector-component
    allocate( iloTable(2 + atoms%nlod, 0: atoms%lmaxd, atoms%ntype) )
    allocate( uds(0:atoms%lmaxd, atoms%ntype), us(0:atoms%lmaxd, atoms%ntype) )
    allocate( duds(0:atoms%lmaxd, atoms%ntype), dus(0:atoms%lmaxd, atoms%ntype) )
    allocate( ddn(0:atoms%lmaxd, atoms%ntype) )
    allocate( ulos(atoms%nlod, atoms%ntype), dulos(atoms%nlod, atoms%ntype) )
    allocate( uulon(atoms%nlod, atoms%ntype), dulon(atoms%nlod, atoms%ntype) )
    allocate( uloulopn(atoms%nlod, atoms%nlod, atoms%ntype) )
    allocate( rbas1( atoms%jmtd, 2 + atoms%nlod, 0:atoms%lmaxd, atoms%ntype, dimens%jspd) ) !maybe samller dimension see Markus
    allocate( rbas2(atoms%jmtd, 2 + atoms%nlod, 0:atoms%lmaxd, atoms%ntype, dimens%jspd) )
    allocate( nRadFun(0:atoms%lmaxd, atoms%ntype) )

    ! Initialize to avoid segfaults or unprecise multiplication errors
    uuilon = 0.0
    duilon = 0.0
    ulouilopn = 0.0
    uulon = 0.0
    uloulopn = 0.0
    ulos = 0.0
    dulos = 0.0
    f = 0.0
    g = 0.0
    flo = 0.0
    uds = 0.0
    us = 0.0
    duds = 0.0
    dus = 0.0
    ddn = 0.0
    dulon = 0.0
    uloulopn = 0.0

    ! This is initialized to the neutral element of multiplication to avoid division by zero
    rbas1 = 1
    rbas2 = 1

    ! The number of radial functions p for a given atom type and orbital quantum number l is at least 2 ( u and udot ).
    nRadFun=2

    ! Mapping array gives zero if for a given atom type and orbital quantum number l there is no local orbital radial solution u_LO
    iloTable = 0

    ! Generate Arrays which count the number of radial solutions for a given atom type and orbital quantum number l and give the
    ! number of the local orbital if the index which counts the radial solutions for a given atom type and orbital quantum number l
    ! is given as well as the orbital quantum number l and the atom type themselves.
    allocate(ilo2p(maxval(atoms%nlo), atoms%ntype))
    ilo2p = 0
    do isp = 1, input%jspins
      do itype = 1, atoms%ntype
        do ilo = 1, atoms%nlo(itype)
          oqn_l = atoms%llo(ilo, itype)
          nRadFun(oqn_l, itype) = nRadFun(oqn_l, itype) + 1
          ilo2p(ilo, itype) = nRadFun(oqn_l, itype)
          iloTable(nRadFun(oqn_l, itype), oqn_l, itype) = ilo
        end do ! ilo
      end do ! itype
    end do ! isp


    ! Generate radial solutions with from Fleur recycled routines radfun and radflo. Rearrange them into rbas1 and rbas2 to not give
    ! LOs a special treatment but to only loop over all radial functions and treat LOs implicitely
    do isp = 1, input%jspins
      do itype = 1, atoms%ntype
        do oqn_l = 0, atoms%lmax(itype)

          ! The El from oqn_l = 3 upwards are the same, so we just use El(oqn_l = 4) here.
          call radfun( oqn_l, enpara%el0(oqn_l, itype, isp), V0sFleur%vr0(1, itype, 1), atoms%jri(itype), atoms%rmsh(1, itype), &
            & atoms%dx(itype), atoms%jmtd, f(1, 1, oqn_l), g(1, 1, oqn_l), us(oqn_l, itype), dus(oqn_l, itype), uds(oqn_l, itype), &
            & duds(oqn_l, itype), ddn(oqn_l, itype), nodeu, noded, wronk )

          call radflo( atoms%ntype, atoms%nlod, dimens%jspd, atoms%jmtd, atoms%lmaxd, itype, isp, enpara%ello0(1, 1, isp), &
            & V0sFleur%vr0(1, itype, 1), atoms%jri(itype), atoms%rmsh(1, itype), atoms%dx(itype), f, g, atoms%llo, atoms%nlo, &
            & atoms%l_dulo(1, itype), mpi%irank, atoms%ulo_der, ulos, dulos, uulon, dulon, uloulopn, uuilon, duilon, ulouilopn, flo)

          do igrid = 1, atoms%jri(itype) ! was jmtd in former times, which is not necessary

            rbas1(igrid, 1, oqn_l, itype, isp) = f(igrid, 1, oqn_l)
            rbas2(igrid, 1, oqn_l, itype, isp) = f(igrid, 2, oqn_l)
            rbas1(igrid, 2, oqn_l, itype, isp) = g(igrid, 1, oqn_l)
            rbas2(igrid, 2, oqn_l, itype, isp) = g(igrid, 2, oqn_l)

          end do
          do p = 3, nRadFun(oqn_l, itype)
              ! p > 2 indexes LOs
            do igrid = 1, atoms%jri(itype)
             rbas1(igrid,p, oqn_l, itype, isp) = flo(igrid, 1,iloTable(p, oqn_l, itype))
             rbas2(igrid,p, oqn_l, itype, isp) = flo(igrid, 2,iloTable(p, oqn_l, itype))
            end do
          end do
        end do
      end do
    end do

    if (compPhon) then
      open(109,file='000_u',form='FORMATTED',position='append',action='WRITE',status='REPLACE')
      do itype = 1, atoms%ntype
        do oqn_l = 0, atoms%lmax(itype)
          do igrid = 1, atoms%jri(itype)
            write(109,*) igrid, 1, oqn_l, itype, 1
            write(109,*) rbas1(igrid, 1, oqn_l, itype, 1)
            write(109,*) igrid, 1, oqn_l, itype, 2
            write(109,*) rbas2(igrid, 1, oqn_l, itype, 1)
            write(109,*) igrid, 2, oqn_l, itype, 1
            write(109,*) rbas1(igrid, 2, oqn_l, itype, 1)
            write(109,*) igrid, 2, oqn_l, itype, 2
            write(109,*) rbas2(igrid, 2, oqn_l, itype, 1)
          end do
        end do
      end do
      close(109)
    end if

    ! Write radial solutions into modern usdus type
    allocate(usdus%us(0:atoms%lmaxd, atoms%ntype, input%jspins), usdus%dus(0:atoms%lmaxd, atoms%ntype, input%jspins), &
      & usdus%uds(0:atoms%lmaxd, atoms%ntype, input%jspins), usdus%duds(0:atoms%lmaxd, atoms%ntype, input%jspins), &
      & usdus%ddn(0:atoms%lmaxd, atoms%ntype, input%jspins), usdus%ulos(atoms%nlod, atoms%ntype, input%jspins), &
      & usdus%dulos(atoms%nlod, atoms%ntype, input%jspins), usdus%uulon(atoms%nlod, atoms%ntype, input%jspins), &
      & usdus%dulon(atoms%nlod, atoms%ntype, input%jspins), usdus%uloulopn(atoms%nlod, atoms%nlod, atoms%ntype, input%jspins))

    usdus%us   (0:atoms%lmaxd, 1:atoms%ntype, 1) = us  (0:atoms%lmaxd, 1:atoms%ntype)
    usdus%dus  (0:atoms%lmaxd, 1:atoms%ntype, 1) = dus (0:atoms%lmaxd, 1:atoms%ntype)
    usdus%uds  (0:atoms%lmaxd, 1:atoms%ntype, 1) = uds (0:atoms%lmaxd, 1:atoms%ntype)
    usdus%duds (0:atoms%lmaxd, 1:atoms%ntype, 1) = duds(0:atoms%lmaxd, 1:atoms%ntype)
    usdus%ddn  (0:atoms%lmaxd, 1:atoms%ntype, 1) = ddn (0:atoms%lmaxd, 1:atoms%ntype)
    usdus%ulos (:, :, 1) = ulos(:, :)
    usdus%dulos(:, :, 1) = dulos(:, :)
    usdus%uulon(:, :, 1) = uulon(:, :)
    usdus%dulon(:, :, 1) = dulon(:, :)
    usdus%uloulopn(:, :, :, 1) = uloulopn(:, :, :)

    call LogRadSolutions ( atoms, usdus, nRadFun, iloTable, uuilon, duilon, ulouilopn, logUnit )

  end subroutine GenRadSolnInt

  !>--------------------------------------------------------------------------------------------------------------------------------
  !> @author
  !> Christian-Roman Gerhorst and Markus Betzinger
  !>
  !> @brief
  !> Unites the energy parameter of the FLAPW basis and the LOs into one energy parameter array.
  !>
  !> @details
  !>
  !>
  !> @param[in] atoms      : Atoms type, see types.f90.
  !> @param[in] enpara     : Energy parameter type, see types.f90.
  !> @param[in] dimens     : Dimension type, see types.f90.
  !> @param[out] El        : Contains LAPW and LOs energy parameters.
  !>
  !>--------------------------------------------------------------------------------------------------------------------------------
  subroutine uniteEnergyParameters( atoms, enpara, dimens, El )

    use m_types

    implicit none

    type(t_atoms),               intent(in)  :: atoms
    type(t_enpara),              intent(in)  :: enpara
    type(t_dimension),           intent(in)  :: dimens
    real,           allocatable, intent(out) :: El(:, :, :, :)

    integer                                  :: oqn_l !loop variable
    integer                                  :: itype !loop variable
    integer                                  :: isp   !loop variable
    integer                                  :: ilo   !loop variable

    ! Fill in regular energy parameter to index p = 1 similiar to radial solutions
    allocate(El(1 + atoms%nlod, 0:atoms%lmaxd, atoms%ntype, dimens%jspd))
    do isp = 1, dimens%jspd
      do itype = 1, atoms%ntype
        do oqn_l = 0, atoms%lmaxd
          El(1, oqn_l, itype, isp) = enpara%el0(oqn_l, itype, isp)
        end do
      end do
    end do

    ! Fill in LO energy parameter to index p (indices LOs for p > 1)
    do isp = 1, dimens%jspd
      do itype = 1, atoms%ntype
        do ilo = 1, atoms%nlo(itype)
          El(1 + ilo, atoms%llo(ilo, itype), itype, isp) = enpara%ello0(ilo, itype, isp)
        end do
      end do
    end do

  end subroutine uniteEnergyParameters

  !>--------------------------------------------------------------------------------------------------------------------------------
  !> @author
  !> Christian-Roman Gerhorst and Markus Betzinger
  !>
  !> @brief
  !> Generates reciprocal symmetry operations.
  !>
  !> @details
  !>
  !> @param[out] sym : Symmetries type, see types.f90.
  !>
  !>--------------------------------------------------------------------------------------------------------------------------------
  subroutine GenerateRecSymOp(sym)

    use m_types

    implicit none

    ! Type parameters
    type(t_sym), intent(inout) :: sym

    ! Scalar variables
    integer                    :: isym ! loop variable

    allocate(sym%mrrot(3, 3, sym%nop))
    do isym = 1, sym%nop
      sym%mrrot(:, :, isym) = transpose(sym%mrot(:, :, sym%invtab(isym)))
    end do

  end subroutine GenerateRecSymOp

  !>--------------------------------------------------------------------------------------------------------------------------------
  !> @author
  !> Christian-Roman Gerhorst and Markus Betzinger
  !>
  !> @brief
  !> k- or q-point set generator.
  !>
  !> @details
  !>
  !> @param[in] nkpt3   : number of k-points per direction
  !> @param[in] kptadd  : shift of kpoints
  !> @param[in] nop     : Number of symmetry operations
  !> @param[in] mrrot   : symmetry operations in reciprocal space
  !> @param[out] nkpt   : number of k-points in IBZ (irreducible Brillouin zone)
  !> @param[out] nkptf  : number of k-kpoints in full brillouin zone
  !> @param[out] bkp    : stores index of parent's k-point
  !> @param[out] bksym  : symmetry operation to map parent to child kpoint
  !> @param[out] bk     : holds all k-points, at first all irreducible then the remaining of the full Brillouin zone.
  !> @param[out] wtkpt  : weight of k-point, i.e. to how many points can a irreducible k-point be mapped to.
  !>
  !> @todo removal of multiple Gamma points if additional qpoints not yet working
  !>
  !>--------------------------------------------------------------------------------------------------------------------------------
  subroutine GenKnQPointSet( nkpt3, kptadd, nop, mrrot, nkpt, nkptf, bkp, bksym, bk, wtkpt )

    use mod_juPhonUtils

    implicit none

    ! Scalar parameters
    integer,                        intent(in)   :: nop
    integer,                        intent(out)  :: nkpt
    integer,                        intent(out)  :: nkptf

    ! Array parameters
    integer,                        intent(in)   :: nkpt3(:)
    real,                           intent(in)   :: kptadd(3)
    integer,                        intent(in)   :: mrrot(:, :, :)
    real,    allocatable,           intent(out)  :: bk(:, :)
    integer, allocatable,           intent(out)  :: bkp(:)
    integer, allocatable,           intent(out)  :: bksym(:)
    real,    allocatable, optional, intent(out)  :: wtkpt(:)

    ! Scalar variables
    integer                                      :: i1, i2       ! temporary indices to aquire number of (irreducibel) k-points
    integer                                      :: k1, k2, k3   ! temporary k points
    integer                                      :: isym         ! loop variable for symmetry operation
    integer                                      :: ikpt         ! loop variable for k-points
    integer                                      :: maxkdim      ! maximal k-set dimension

    ! Array variables
    logical, allocatable                         :: def(:, :, :) ! stores whether k-point is defined or not
    integer                                      :: k(3)         ! stores current k-point


    ! Calculate total number of k-points
    nkptf = nkpt3(1) * nkpt3(2) * nkpt3(3)

    ! Allocate and initialize parameters and variables
    allocate( bk(3, nkptf), bkp(nkptf), bksym(nkptf) )
    allocate( def(0:nkpt3(1)-1, 0:nkpt3(2)-1, 0:nkpt3(3)-1) )
    def = .false.
    i1  = 0
    i2  = nkptf + 1

    ! Loop over all possible k-point. If a new k-point is found (.false. in def array) then find all k-points related to it by
    ! symmetry and mark them found in the def array.
    do k1 = 0, nkpt3(1) - 1
      do k2 = 0, nkpt3(2) - 1
        do k3 = 0, nkpt3(3) - 1
          if(.not.def(k1, k2, k3)) then
            ! irreducible k-points stretch from 1 to i1
            i1              = i1 + 1
            ! project into 1st brillouin zone again
            bk(:, i1)      = modulo1r( 1d0 * (/k1, k2, k3/) / nkpt3 + kptadd )
            bkp(i1)         = i1
            bksym(i1)       = 1
            def(k1, k2, k3) = .true.
            do isym = 1, nop
              ! skip operations that do not leave kptadd invariant
              if( any( modulo1r( matmul( mrrot(:, :, isym), kptadd ) - kptadd ) >= 1d-12 ) ) cycle
              k = modulo( matmul( mrrot(:, :, isym), (/k1, k2, k3/) ), nkpt3 )
              if( .not.def(k(1), k(2), k(3)) ) then
                i2                    = i2 - 1
                bk(:,i2)              = modulo1r (1d0 * k / nkpt3 + kptadd)
                bkp(i2)               = i1
                bksym(i2)             = isym
                def(k(1), k(2), k(3)) = .true.
              end if
            end do ! isym
          end if
        end do ! k3
      end do ! k2
    end do ! k1

    !Todo this has to be refined and reviewed also the condition if additional k-points are added
    if ( sum(kptadd) > 1e-12 ) then
      bk(:,i2:)  = bk(:, (/(i1, i1=nkptf, i2, -1)/))
      bkp(i2:)   = bkp((/(i1, i1=nkptf, i2,-1)/))
      bksym(i2:) = bksym((/(i1, i1=nkptf, i2, -1)/))
    end if

    nkpt = i1

    ! get weight of k-points
    if (present(wtkpt)) then
      allocate(wtkpt(nkptf))
      wtkpt = 0
      do ikpt = 1, nkptf
        wtkpt(bkp(ikpt)) = wtkpt(bkp(ikpt)) + 1
      end do
    endif

  end subroutine GenKnQPointSet

  !>--------------------------------------------------------------------------------------------------------------------------------
  !> @author
  !> Christian-Roman Gerhorst and Markus Betzinger
  !>
  !> @brief
  !> Creates and checks the kpts file.
  !>
  !> @details
  !> This routine also supports adding additional phonon wave-vectors q and to shift the whole k-point set
  !>
  !> @param[in] qSetDim         : Dimension of the q-point set
  !> @param[in] logUnit         : Globel unit for the log file juPhon.log
  !> @param[in] kSetDim         : Dimension of the k-point set
  !> @param[in] addQs           : Additional phonon wave-vectors added to the k-point set
  !> @param[in] kModeSw         : Switch for k-point generation mode
  !> @param[in] writeKpqArraySw : kpqArray is written into file kpqMapArray for debugging.
  !> @param[out] sym            : Symmetries type, see types.f90.
  !> @param[out] addQnkptis     : number of irreducible k-points which are added due to additional q-point inducing an additional
  !>                              shifted k-point set
  !> @param[out] nrAddQs        : Number of additional phonon wave-vectors
  !> @param[in]  startTime      : Starting time of juPhon.
  !> @param[out] addQnkptfs     : number of all k-points which are added due to additional q-point inducing an additional shifted
  !>                              k-point set
  !> @param[inout] kpts         : K-points type, see types.f90.
  !> @param[inout] qpts         : Q-points type, see types.f90.
  !> @param[out] mapKpq2K       : For a given k-point index and q-point index, this array gives the k-point set index of the result
  !>                              k + q (already mapped back to Brillouin zone).
  !>
  !> @todo : implement that regular k-point set is shifted from the beginning
  !> @todo : ensure that multiple Gamma Points are thrown out from kpts set, until now there were not such problems!
  !>--------------------------------------------------------------------------------------------------------------------------------
  subroutine GenNCheckKnQPointInfo( qSetDim, logUnit, kSetDim, addQs, kModeSw, writeKpqArraySw, sym, addQnkptis, nrAddQs, &
      & startTime, addQnkptfs, kpts, qpts, mapKpq2K, kpq2kPrVec )

    use mod_juPhonUtils, only : kgv, fopen, fclose
    use m_juDFT_NOstopNO, only : juDFT_error
    use m_types
    use m_jpLog, only : finishLogFile
    implicit none

    ! type parameter
    type(t_sym),               intent(in)    :: sym
    type(t_kpts),              intent(inout) :: kpts
    type(t_kpts),              intent(inout) :: qpts

    ! scalar parameter
    integer,                   intent(in)    :: logUnit
    logical,                   intent(in)    :: kModeSw
    logical,                   intent(in)    :: writeKpqArraySw
    integer,                   intent(out)   :: nrAddQs

    ! array parameter
    integer,                   intent(in)    :: startTime(8)
    integer,      allocatable, intent(in)    :: qSetDim(:)
    integer,      allocatable, intent(in)    :: kSetDim(:)
    real,         allocatable, intent(in)    :: addQs(:)
    integer,      allocatable, intent(out)   :: addQnkptis(:)
    integer,      allocatable, intent(out)   :: addQnkptfs(:)
    integer,      allocatable, intent(out)   :: mapKpq2K(:, :)
    integer,      allocatable, intent(out)   :: kpq2kPrVec(:, :, :)

    ! scalar variable
    integer                                  :: cumulativeSumI            ! Helping index to create (e.g. kpts%)bk
    integer                                  :: cumulativeSumF            ! Helping index to create (e.g. kpts%)bk
    integer                                  :: ikpt                      ! loop variable
    integer                                  :: indexOfZero               ! index in bk where zero vector occurs multiple times
    integer                                  :: nkpt_proof                ! used to check nkpt in kpts file
    real                                     :: lcm_proof                 ! used to check lcm in kpts file
    real                                     :: bk_proof(3)               ! used to check k-points in kpts file
    real                                     :: wtkpt_proof               ! used to check weight in kpts file
    real                                     :: qDuplicate(3)             ! used to check for duplicates of phonon vectors q
    real                                     :: lcm                       ! factor with which k-points in kpts file are multiplied
    integer                                  :: iqadd                     ! loop variable
    integer                                  :: iqpt                      ! loop variable
    real                                     :: sumKptW                   ! stores sum of k-point weights to determine correct wtkpt

    ! array variable
    ! Help to handle additional q-points
    integer,      allocatable                :: bkp_temp(:)               ! helps to setup (e.g. kpts%)kpts%bkp
    integer,      allocatable                :: bksym_temp(:)             ! helps to setup (e.g. kpts%)kpts%bksym
    real,         allocatable                :: bk_temp(:, :)             ! helps to setup (e.g. kpts%)kpts%bk
    integer,      allocatable                :: wtkpt_temp(:)             ! helps to setupt (e.g. kpts%)wtkpt
    integer,      allocatable                :: bkp_tempContainer(:, :)   ! helps to setup (e.g. kpts%)kpts%bkp
    integer,      allocatable                :: bksym_tempContainer(:, :) ! helps to setup (e.g. kpts%)kpts%bksym
    real,         allocatable                :: bk_tempContainer(:, :, :) ! helps to setup (e.g. kpts%)kpts%bk
    character(len=1024)                      :: errorMessage              ! stores error message for error output

    ! set pointer array addQnkptis to minimal size, i.e. no additional q-points are assumed and initialize
    allocate( addQnkptis(1) )
    addQnkptis = 0
    nrAddQs = 0

    ! generate Γ-centered k-point set, k-point type arrays are cautiously stored in temp arrays, to handle additional k-points
    ! todo if normal k-set should be shifted, variable from input should be build in here
    if ( kModeSw) then
      call GenKnQPointSet( kSetDim, [0., 0., 0.], sym%nop, sym%mrrot, kpts%nkpt, kpts%nkptf, bkp_temp, bksym_temp, bk_temp, &
        & wtkpt=kpts%wtkpt )
    else
      ! nop is set to 1 to avoid rotating wavefunctions
      call GenKnQPointSet( kSetDim, [0., 0., 0.], 1, sym%mrrot, kpts%nkpt, kpts%nkptf, bkp_temp, bksym_temp, bk_temp, &
        & wtkpt=kpts%wtkpt )
    end if

    ! Γ-point is already considered and would be a duplicate
    if ( any(addQs(:) /= 0) ) then
      ! q-points are stored in sequencing tripel bunches (x, y, z)
      nrAddQs = size(addQs) / 3
      ! adjust final k- and q-point arrays and scalars according to additional q-points and initialize
      deallocate( addQnkptis )
      allocate( addQnkptis(nrAddQs) )
      allocate( addQnkptfs(nrAddQs) )
      addQnkptis = 0
      addQnkptfs = 0
      allocate( kpts%bkp(kpts%nkptf + nrAddQs * kpts%nkptf) )
      allocate( kpts%bksym(kpts%nkptf + nrAddQs * kpts%nkptf) )
      if( .not.kModeSw ) then
        deallocate(kpts%bk)
      endif
      allocate( kpts%bk(3, kpts%nkptf + nrAddQs * kpts%nkptf) )
      kpts%bkp = 0
      kpts%bksym = 0
      kpts%bk = 0

      ! Container to store all k-point sets, regular k-point set is stored first at index 0, for every additional q-point a
      ! respective k-point set has to be generated which are stored in the subsequent indices
      allocate( bkp_tempContainer(0:nrAddQs, kpts%nkptf) )
      allocate( bksym_tempContainer(0:nrAddQs, kpts%nkptf) )
      allocate( bk_tempContainer(0:nrAddQs, 3, kpts%nkptf) )
      bkp_tempContainer(0, :) = bkp_temp(:)
      bksym_tempContainer(0, :) = bksym_temp(:)
      bk_tempContainer(0, :, :) =  bk_temp(:, :)

      ! Iterate over all additional q-points and generate their k-point sets and store them to just mentioned container arrays
      do iqadd = 1, nrAddQs
        if ( kModeSw ) then
          call GenKnQPointSet(kSetDim, addQs((iqadd - 1) * 3 + 1 : iqadd * 3), sym%nop, sym%mrrot, addQnkptis(iqadd),              &
            & addQnkptfs(iqadd), bkp_temp, bksym_temp, bk_temp)
        else
          ! nop is set to 0 to avoid rotation of wavefunctions
          call GenKnQPointSet(kSetDim, addQs((iqadd - 1) * 3 + 1 : iqadd * 3), 1, sym%mrrot, addQnkptis(iqadd), addQnkptfs(iqadd), &
            & bkp_temp, bksym_temp, bk_temp)
        end if
        bkp_tempContainer(iqadd, :) = bkp_temp(:)
        bksym_tempContainer(iqadd, :) = bksym_temp(:)
        bk_tempContainer(iqadd, :, :) = bk_temp(:, :)
        deallocate(bkp_temp, bksym_temp, bk_temp)
      end do

      ! Collaps temp arrays into final k- and q-point arrays and determine all relevant scalars
      ! Begin with irreducible regular k-point set
      kpts%bkp(1:kpts%nkpt) = bkp_tempContainer(0, 1:kpts%nkpt)
      kpts%bksym(1:kpts%nkpt) = bksym_tempContainer(0, 1:kpts%nkpt)
      kpts%bk(:, 1:kpts%nkpt) = bk_tempContainer(0, :, 1:kpts%nkpt)

      ! Add irreducible additional shifted k-point sets
      cumulativeSumI = 0
      do iqadd = 1, nrAddQs
        kpts%bksym(kpts%nkpt + cumulativeSumI + 1:kpts%nkpt + cumulativeSumI + addQnkptis(iqadd)) = &
          & bksym_tempContainer(iqadd, 1:addQnkptis(iqadd))
        kpts%bkp(kpts%nkpt + cumulativeSumI + 1:kpts%nkpt + cumulativeSumI + addQnkptis(iqadd)) = &
          & bkp_tempContainer(iqadd, 1:addQnkptis(iqadd))
        kpts%bk(:, kpts%nkpt + cumulativeSumI + 1:kpts%nkpt + cumulativeSumI + addQnkptis(iqadd)) = &
          & bk_tempContainer(iqadd, :, 1:addQnkptis(iqadd))
        cumulativeSumI = cumulativeSumI + addQnkptis(iqadd)
      end do

      ! Apend child k-points after all (regular plus additional q-point induced) irreducible k-point sets; start with regular      &
      ! k-point set In order to use nkptf we have to insert sum(addQnkptis) points which are between end of nkpt and beginning of  &
      ! nkptf points size(array(i+1:i)) = 0 in Fortran!
      kpts%bkp(kpts%nkpt + sum(addQnkptis) + 1: sum(addQnkptis) + kpts%nkptf) = bkp_tempContainer(0, kpts%nkpt+1:kpts%nkptf)
      kpts%bksym(kpts%nkpt + sum(addQnkptis) + 1: sum(addQnkptis) + kpts%nkptf) = bksym_tempContainer(0, kpts%nkpt+1:kpts%nkptf )
      kpts%bk(:, kpts%nkpt + sum(addQnkptis) + 1: sum(addQnkptis) + kpts%nkptf) = bk_tempContainer(0, :, kpts%nkpt+1:kpts%nkptf)

      ! Continue with additional shifted k-point sets
      ! The difference addQnkptfs(iqadd) - addQnkptis(iqadd)) as there are additional q-points which have child k-points and not
      cumulativeSumF = 0
      do iqadd = 1, nrAddQs
        kpts%bksym(kpts%nkptf + sum(addQnkptis) + cumulativeSumF + 1:kpts%nkptf + sum(addQnkptis) + cumulativeSumF &
          & + (addQnkptfs(iqadd) - addQnkptis(iqadd))) = bksym_tempContainer(iqadd, addQnkptis(iqadd) + 1:addQnkptfs(iqadd))
        kpts%bkp(kpts%nkptf + sum(addQnkptis) + cumulativeSumF + 1:kpts%nkptf + sum(addQnkptis) + cumulativeSumF &
          & + (addQnkptfs(iqadd) - addQnkptis(iqadd))) = bkp_tempContainer(iqadd, addQnkptis(iqadd) + 1:addQnkptfs(iqadd))
        kpts%bk(:, kpts%nkptf + sum(addQnkptis) + cumulativeSumF + 1:kpts%nkptf + sum(addQnkptis) + cumulativeSumF &
          & + (addQnkptfs(iqadd)- addQnkptis(iqadd))) = bk_tempContainer(iqadd, :, addQnkptis(iqadd) + 1:addQnkptfs(iqadd))

        cumulativeSumF = cumulativeSumF + addQnkptfs(iqadd)
      end do

    else
      ! copy temp arrays to final k-point array
      allocate( kpts%bkp(kpts%nkptf) )
      allocate( kpts%bksym(kpts%nkptf) )
      if (kModeSw) then
        allocate( kpts%bk(3, kpts%nkptf) )
      endif
      kpts%bkp = bkp_temp
      kpts%bksym = bksym_temp
      kpts%bk = bk_temp
    endif

    ! Adding of additional k-points not yet implemented.
    if (.false.) then
      !kpts unit is 202
      ! Γ-point might be multiple times in k-point set so duplicates have to be removed
      ! todo Problem of this algorithm arrays which are scanned are shortened during the scanning
      if (any(addQs(:) /= 0)) then
        allocate( bkp_temp(size(kpts%bkp)), bksym_temp(size(kpts%bksym)), bk_temp(3, size(kpts%bk, 2)), &
                                                                                                    & wtkpt_temp(size(kpts%wtkpt)) )
        indexOfZero = 0
        do ikpt = 2, kpts%nkptf + sum(addQnkptfs)
          if ( all( abs( kpts%bk(:, ikpt) - 0.0 ) < 1e-10 ) ) then
            deallocate( bk_temp, bksym_temp, bkp_temp, wtkpt_temp )
            allocate( bk_temp(3, size( kpts%bk, 2 )), bksym_temp(size( kpts%bksym )), bkp_temp(size( kpts%bkp )), &
              & wtkpt_temp(size( kpts%wtkpt )) )
            bkp_temp = 0
            bksym_temp = 0
            bk_temp = 0
            wtkpt_temp = 0
            indexOfZero = ikpt
            bkp_temp = kpts%bkp
            bksym_temp = kpts%bksym
            bk_temp = kpts%bk
            wtkpt_temp = kpts%wtkpt
            deallocate(kpts%bkp, kpts%bk, kpts%bksym, kpts%wtkpt)
            allocate(kpts%bkp(size(bkp_temp) - 1), kpts%bksym(size(bksym_temp) - 1), kpts%bk(3, size(bk_temp, 2) - 1), &
              & kpts%wtkpt(size(wtkpt_temp) - 1))
            kpts%bkp(:indexOfZero - 1) = bkp_temp(:indexOfZero - 1)
            kpts%bkp(indexOfZero:) = bkp_temp(indexOfZero + 1 :)
            kpts%bksym(:indexOfZero - 1) = bksym_temp(:indexOfZero - 1)
            kpts%bksym(indexOfZero:) = bksym_temp(indexOfZero + 1 :)
            kpts%bk(:, :indexOfZero - 1) = bk_temp(:, :indexOfZero - 1)
            kpts%bk(:, indexOfZero:) = bk_temp(:, indexOfZero + 1 :)
            kpts%wtkpt( :indexOfZero - 1) = wtkpt_temp(:indexOfZero - 1)
            kpts%wtkpt(indexOfZero:) = wtkpt_temp(indexOfZero + 1 :)

            ! fix pointers
            indexOfZero = indexOfZero - kpts%nkpt
            if ( indexOfZero < 0 ) then
              kpts%nkpt = kpts%nkpt - 1
              kpts%nkptf = kpts%nkptf - 1
            else
              do iqadd = 1, nrAddQs
                indexOfZero =  indexOfZero - addQnkptis(iqadd)
                if (indexOfZero < 0) then
                  addQnkptis(iqadd) = addQnkptis(iqadd) - 1
                  addQnkptfs(iqadd) = addQnkptfs(iqadd) - 1
                  exit
                end if
              end do
            end if
            indexOfZero = indexOfZero - ( kpts%nkptf - kpts%nkpt )
            if ( indexOfZero < 0 ) then
              kpts%nkptf = kpts%nkptf - 1
            else
              do iqadd = 1, nrAddQs
                indexOfZero =  indexOfZero - ( addQnkptfs(iqadd) - addQnkptis(iqadd) )
                if (indexOfZero < 0) then
                  addQnkptfs(iqadd) = addQnkptfs(iqadd) - 1
                  exit
                end if
              end do
            end if
          end if
        end do
      end if
    end if

    ! Generation of q-point set
    if ( modulo( kSetDim(1), qSetDim(1) ) == 0 .and. modulo( kSetDim(2), qSetDim(2) ) == 0 .and. &
      & modulo( kSetDim(3), qSetDim(3) ) == 0 ) then

      ! TODO kptadd into general type?
      ! If in the end the mode of operation is clear q-point set will only be generated if not in k-point mode!
      if ( kModeSw ) then
        call GenKnQPointSet( qSetDim, (/ 0., 0., 0./), sym%nop, sym%mrrot, qpts%nkpt, qpts%nkptf, qpts%bkp, qpts%bksym, qpts%bk, &
          & wtkpt=qpts%wtkpt )
      else
        ! nop is set to 0 to avoid rotation of wavefunctions
        call GenKnQPointSet( qSetDim, (/ 0., 0., 0./), 1, sym%mrrot, qpts%nkpt, qpts%nkptf, qpts%bkp, qpts%bksym, qpts%bk, &
          & wtkpt=qpts%wtkpt )
      end if

      !add additional q-points to q-set
      if ( any( addQs(:) /= 0 ) ) then
        ! Find duplicates of q, awkward indexing comes from the fact that additional q-points are stored in subsequent triple
        ! bunches
        do iqadd = 1, nrAddQs
          ! Is q-vector in first Brillouin zone?
          qDuplicate(1) = addQs((iqadd - 1) * 3 + 1)
          qDuplicate(2) = addQs((iqadd - 1) * 3 + 2)
          qDuplicate(3) = addQs((iqadd - 1) * 3 + 3)
          if ( any( qDuplicate > 1 ) ) then
            write( errorMessage, '(i2,a,3(f14.6),a)' ) iqadd, '. q-vector ', qDuplicate, ' is not in  first Brillouin zone.'
            call juDFT_error( trim(errorMessage), calledby='GenNCheckKnQPointInfo', &
              & hint='Please change this q-vector in input file and generate new kpts file.' )
          endif
          ! Is additional q-point already part of q-point set?
          do ikpt = 1, qpts%nkptf
            if ( all( abs( qpts%bk(:, ikpt) - qDuplicate(:) ) < 1E-6 ) ) then
              write(errorMessage, '(a,3(f14.6), a)') 'Please remove ', qDuplicate, ' from input file.'
              call juDFT_error( 'Duplicates in q-points detected.', calledby='GenNCheckKnQPointInfo', &
                & hint=errorMessage )
            endif
          end do
        end do

        allocate( bk_temp(3, qpts%nkptf) )
        allocate( bkp_temp(qpts%nkptf) )
        allocate( bksym_temp(qpts%nkptf) )
        allocate( wtkpt_temp(qpts%nkptf) )
        bk_temp = qpts%bk
        bkp_temp = qpts%bkp
        bksym_temp = qpts%bksym
        wtkpt_temp = qpts%wtkpt
        deallocate( qpts%bk, qpts%bkp, qpts%bksym, qpts%wtkpt )
        allocate( qpts%bk(3, qpts%nkptf + nrAddQs) )
        allocate( qpts%bkp(qpts%nkptf + nrAddQs) )
        allocate( qpts%bksym(qpts%nkptf + nrAddQs) )
        allocate( qpts%wtkpt(qpts%nkptf + nrAddQs) )
        qpts%bk(:, :qpts%nkpt) = bk_temp(:, :qpts%nkpt)
        qpts%bkp(:qpts%nkpt) = bkp_temp(:qpts%nkpt)
        qpts%bksym(:qpts%nkpt) = bksym_temp(:qpts%nkpt)
        qpts%wtkpt(:qpts%nkpt) = wtkpt_temp(:qpts%nkpt)
        do iqadd = 1, nrAddQs
          qpts%bk(1, qpts%nkpt + iqadd) = addQs((iqadd - 1) * 3 + 1)
          qpts%bk(2, qpts%nkpt + iqadd) = addQs((iqadd - 1) * 3 + 2)
          qpts%bk(3, qpts%nkpt + iqadd) = addQs((iqadd - 1) * 3 + 3)
          qpts%bkp(qpts%nkpt + iqadd) = 0
          qpts%bksym(qpts%nkpt + iqadd) = 0
          qpts%wtkpt(qpts%nkpt + iqadd) = 0
        enddo
        qpts%bk(:, qpts%nkpt + nrAddQs + 1:) =  bk_temp(:, qpts%nkpt + 1:)
        qpts%bkp(qpts%nkpt + nrAddQs + 1:)   =  bkp_temp(qpts%nkpt + 1:)
        qpts%bksym(qpts%nkpt + nrAddQs + 1:) =  bksym_temp(qpts%nkpt + 1:)
        qpts%wtkpt(qpts%nkpt + nrAddQs + 1:) =  wtkpt_temp(qpts%nkpt + 1:)
      endif

    else
      call juDFT_error( 'q-point set is not part of k-point set.', calledby='GenNCheckKnQPointInfo', &
        & hint='Please ensure dimensions k being multiples of dimensions of q.')
    endif

    ! Some logfile output about results of k-point generator
    if (.not.kModeSw) then
      write ( logUnit, '(a)' ) "Checking read-in kpts-file and generating remaining k-point and q-point information"
      write ( logUnit, '(a)' ) "==================================================================================="
    end if
    write ( logUnit, * )
    write ( logUnit, '(a)' ) 'k-point generator statistics'
    write ( logUnit, '(a)' ) '----------------------------'
    write ( logUnit, '(a,1x,i4)' ) 'Number of regular    irreducible k-points:', kpts%nkpt
    if (nrAddQs > 0) then
      write ( logUnit, '(a,1x,i4)' ) 'Number of additional irreducible k-points:', sum(addQnkptis)
    end if
    write ( logUnit, '(a,1x,i4)' ) 'Number of regular    total       k-points:', kpts%nkptf
    if (nrAddQs > 0) then
      write ( logUnit, '(a,1x,i4)' ) 'Number of additional total       k-points:', sum(addQnkptfs)
    end if
    write ( logUnit, '(a,1x,i4)' ) 'Number of            total       q-points:', qpts%nkptf + nrAddQs

    write ( logUnit, * )
    write ( logUnit, '(a)' ) 'q-points to be evaluated:'
    write ( logUnit, '(a)' ) '-------------------------'
    do iqpt = 1, size( qpts%bk, 2 )
      write ( logUnit, '("("1x,3(f10.8,1x)")")' ) qpts%bk(:, iqpt)
    end do
    write ( logUnit, * )

    if (.false.) then
      ! Add k-points to kpts file not yet implemented
      deallocate(bk_temp)
      allocate(bk_temp(3, kpts%nkpt + 1))
      bk_temp = 0
      bk_temp(:, 1:4) = kpts%bk(:, 1:4)
      bk_temp(:, 5) = [0.,0.,1.25]
      bk_temp(:, 6:kpts%nkpt+1) = kpts%bk(:, 5:kpts%nkpt)
      deallocate(kpts%bk)
      allocate(kpts%bk(3, kpts%nkpt + 1))
      kpts%bk = bk_temp
      kpts%nkpt = kpts%nkpt + 1
      kpts%nkptf = kpts%nkpt
      !deallocate(wtkpt_temp)
      allocate(wtkpt_temp(kpts%nkpt))
      wtkpt_temp(1:4) = kpts%wtkpt(1:4)
      wtkpt_temp(5) = 0
      wtkpt_temp(6:kpts%nkpt) = kpts%wtkpt(5:kpts%nkpt - 1)
      kpts%wtkpt = wtkpt_temp
    end if

    ! If k-mode activated a kpts file will be created, if not the kpts file will be read in and checked
    if ( kModeSw ) then

      ! Write information to kpts file
      lcm = real( kgv(kSetDim, 3) )
      call fopen( 202, name='kpts', status='replace', action='write', form='formatted', access='sequential' )
      write( 202, '(i5,f20.10)' ) kpts%nkpt + sum(addQnkptis), lcm
      do ikpt=1, kpts%nkpt
        write (202, '(4(f10.5))') kpts%bk(:, ikpt) * lcm, kpts%wtkpt(ikpt)
      end do
      if ( any(addQs(:) /= 0) ) then
        do ikpt = kpts%nkpt + 1, kpts%nkpt + sum(addQnkptis)
          write( 202, '(4(f10.5))' ) kpts%bk(:, ikpt) * lcm, 0.
        end do
      endif
      call fclose( 202 )

      ! Write information to juphon.log file
      write( logUnit, '(a)' ) 'k-point set with high-symmetry point Γ generated and written to kpts file.'
      write( logUnit, '(a)' ) '--------------------------------------------------------------------------'
      if ( sym%nop > 1 ) then
      write( logUnit, * )
      write( logUnit, '(a)' ) 'As a next step, please converge out Fleur calculation.'
      else if ( sym%nop == 1 ) then
      write( logUnit, * )
      write( logUnit, '(a)' ) 'As a next step, please run the converged Fleur calculation once again to gain all wavefunctions requ&
        &ired.'
      end if
      write( logUnit, * )
      write( logUnit, '(a)' ) 'Terminating juPhon...'
      write( logUnit, '(a)' ) '====================='

      ! Terminates k-point generation mode
      call finishLogFile( startTime, logUnit )
      call fclose( logUnit, status='keep' )
      NOstopNO

    else
      ! Read in existing kpts file again and compare to kpts which have been created internally to ensure consistency
      call fopen( 202, name='kpts', status='old', action='read', form='formatted', access='sequential' )
      read( 202, '(i5,f20.10)' ) nkpt_proof, lcm_proof
      !write(*, *) nkpt_proof, kpts%nkpt+ sum(addQnkptis)
      if (nkpt_proof /= kpts%nkpt + sum(addQnkptis)) then
        call juDFT_error( 'Number of kpoints incorrect in kpts file.', calledby='GenNCheckKnQPointInfo', &
          & hint='Please generate new kpts file.' )
      endif
      lcm = real( kgv(kSetDim, 3) )
      if ( abs(lcm_proof -  lcm) > 1E-6 ) then
        call juDFT_error( 'Denominator is not correct in first line.', calledby='GenNCheckKnQPointInfo', &
          & hint='Please generate new kpts file.' )
      endif
      do ikpt = 1, kpts%nkpt
        read ( 202, '(4(f10.5))' ) bk_proof, wtkpt_proof
        bk_proof = bk_proof / lcm
        if ( any( abs( kpts%bk(:, ikpt) - bk_proof(:) ) > 1E-7) ) then
          errorMessage = ''
          write (errorMessage,'("k-point",1x,i4,1x,"wrong.")') ikpt
          call juDFT_error( errorMessage, calledby='GenNCheckKnQPointInfo', hint='Please generate new kpts file.' )
        endif

        if ( kpts%wtkpt(ikpt) /= wtkpt_proof ) then
          errorMessage=''
          write (errorMessage,'("Weight of k-point",1x,i4,1x,"wrong.")') ikpt
          call juDFT_error( errorMessage, calledby='GenNCheckKnQPointInfo', hint='Please generate new kpts file.' )
        endif
      end do
      if ( any(addQs(:) /= 0) ) then
        do ikpt = kpts%nkpt + 1, kpts%nkpt +  sum(addQnkptis)
          read ( 202, '(4(f10.5))' ) bk_proof(:), wtkpt_proof
          bk_proof = bk_proof / lcm
          if ( any( abs( kpts%bk(:, ikpt) - bk_proof(:) ) > 1E-7) ) then
            errorMessage=''
            write (errorMessage,'("k-point",1x,i4,1x,"wrong.")') ikpt
            call juDFT_error( errorMessage, calledby='GenNCheckKnQPointInfo', hint='Please generate new kpts file.' )
          endif
          if ( wtkpt_proof /= 0 ) then
            errorMessage=''
            write (errorMessage,'("Weight of k-point",1x,i4,1x," in shifted k-point set wrong.")') ikpt
            call juDFT_error( errorMessage, calledby='GenNCheckKnQPointInfo', hint='Please generate new kpts file.' )
          endif
        end do
      endif
      call fclose( 202 )

      call createkqMapArrays( kpts, qpts, nrAddQs, lcm, kSetDim, addQnkptis, mapKpq2K, kpq2kPrVec )

      write ( logUnit, '(a)' ) 'Created k + q => k mapping array'
      write ( logUnit, '(a)' ) '================================'

#ifdef DEBUG_MODE

      if ( writeKpqArraySw ) then
        ! Write out of k + q mapping array
        call fopen( 1000, name='kpqMapArray', status='replace', action='write', form='formatted', access='sequential' )
        write ( 1000, '("k-index",16x,"k-coordinates",19x,"q-index",16x,"q-coordinates",21x,"kpq-index",14x,"kpq-coordinates",/)' )
        do iqpt = 1, qpts%nkpt + nrAddQs
          do ikpt = 1, kpts%nkpt
            write ( 1000, '(i3,2x,3(f15.8),7x,i3,3(f15.8),10x,i3,3(f15.8))') ikpt, kpts%bk(:, ikpt), iqpt, qpts%bk(:, iqpt), &
              &mapKpq2K(ikpt,iqpt), kpts%bk(:, mapKpq2K(ikpt,iqpt))
          end do
        end do
        call fclose( 1000 )
      end if

      ! Check k + q -> k mapping array
      do iqpt = 1, qpts%nkpt + nrAddQs
        do ikpt = 1, kpts%nkpt
          if ( norm2( mod ( kpts%bk(:, ikpt) + qpts%bk(:, iqpt) , 1.0 ) - kpts%bk(:, mapKpq2K(ikpt, iqpt)) ) > 1e-7 ) then
            errorMessage=''
            write (errorMessage,'("Inconsistency in determing index of k + q with k-index",1x,i3,1x,"and q-index",1x,i3)') ikpt,   &
              & iqpt
            call juDFT_error( errorMessage, calledby='GenNCheckKnQPointInfo', hint='Please check mapKpq2K array!' )
          end if
        end do
      end do

#endif
      write ( logUnit, * ) ''
      write ( logUnit, '(a)' ) 'Check of k-point and q-point information passed!'
      write ( logUnit, '(a)' ) '================================================'
      write ( logUnit, * ) ''
    endif

    ! Correct kpts to be not an integer but the integer of the kpts file divided by the number of irreducible k-points
    sumKptW = 0
    sumKptW = sum( kpts%wtkpt(:kpts%nkpt) )
    if ( sumKptW == 0 ) then
      call juDFT_error( 'Weights of k-points damaged.', calledby='GenNCheckKnQPointInfo', hint='Please check generation of kpts%wkp&
        &ts!' )
    end if
    kpts%wtkpt = kpts%wtkpt /sumKptW


  end subroutine GenNCheckKnQPointInfo

  !---------------------------------------------------------------------------------------------------------------------------------
  !> @author
  !> Christian-Roman Gerhorst and Markus Betzinger
  !>
  !> @brief
  !> Determine indices of k-vectors which are equals to the respective k + q combination
  !>
  !> @details
  !> For every k and q vector, there a k + q vector can be found which is again a k-vector of the k-point set (regular or shifted in
  !> addition). This routine determines the mapping between every k and q to their respective k + q vector within the Brillouin zone
  !> and determines connects the indices of the k and q-vectors to the index of the respective k + q vector which is again element
  !> of the k-point set.
  !>
  !> @param[in] kpts       : K-points type, see types.f90.
  !> @param[in] qpts       : Q-points type, see types.f90.
  !> @param[in] nrAddQs    : Number of additional phonon wave-vectors
  !> @param[in] lcm        : Factor with which k-points in k-point file are multiplied
  !> @param[in] kSetDim    : Dimension of the k-point set
  !> @param[in] addQnkptis : Number of irreducible k-points which are added due to additional q-point inducing an additional shifted
  !>                         k-point set
  !> @param[out] mapKpq2K  : For a given k-point index and q-point index, this array gives the k-point set index of the result k + q
  !>                         (already mapped back to Brillouin zone).
  !>
  !---------------------------------------------------------------------------------------------------------------------------------
  subroutine createkqMapArrays( kpts, qpts, nrAddQs, lcm, kSetDim, addQnkptis, mapKpq2K, kpq2kPrVec )

    use m_types
    use mod_juPhonUtils
    use m_juDFT_NOstopNO, only : juDFT_error

    implicit none

    ! Type parameters
    type(t_kpts),              intent(in)  :: kpts
    type(t_kpts),              intent(in)  :: qpts

    ! Array parameter
    integer,                   intent(in)  :: addQnkptis(:)
    integer,                   intent(in)  :: kSetDim(:)
    integer,      allocatable, intent(out) :: mapKpq2K(:, :)
    integer,      allocatable, intent(out) :: kpq2kPrVec(:, :, :)

    ! Scalar parameter
    integer,                   intent(in)  :: nrAddQs
    real,                      intent(in)  :: lcm

    ! Array variable
    integer,      allocatable              :: mapK2Ind(:, :, :) ! takes components of kpts%bk * lcm and gives its index
    integer,      allocatable              :: mapK2mK(:)
    integer                                :: kpqNomin(3)       ! helps to find k + q mapped back to Brillouin zone
    real                                   :: kpqTemp(3)        ! helps to find k + q mapped back to Brillouin zone
    character(len=1024)                    :: errorMessage      ! stores error message for error output

    ! Scalar variable
    integer                                :: ikpt              ! loop variable
    integer                                :: maxKcomp(3)       ! stores maximal k-point component * lcm
    integer                                :: idir              ! loop variable
    integer                                :: iqpt              ! loop variable
    integer                                :: nkptShift         ! stores shift of shifted k-point set
    integer                                :: ikptSh            ! loop variable
    logical                                :: matchFound        ! is true if index for k + q is found

    ! Determine maximal value of k-vector per direction in internal representation for allocation of mapKpq2K array and allocate it
    maxKcomp = 0
    do idir = 1, 3
      maxKcomp(idir) = maxval( kpts%bk(idir, :kpts%nkpt) * lcm )
    end do
    allocate( mapK2Ind(0:maxKcomp(1), 0:maxKcomp(2), 0:maxKcomp(3)) )
    mapK2Ind = 0
    allocate( mapKpq2K(kpts%nkpt, qpts%nkptf + nrAddQs) )

    ! Fill up array which stores the index of a given k-point
    do ikpt = 1, kpts%nkpt
      mapK2Ind(nint(kpts%bk(1, ikpt) * lcm), nint(kpts%bk(2, ikpt) * lcm), nint(kpts%bk(3, ikpt) * lcm)) = ikpt
    end do

    ! Determine k-point on which k + q can be folded back and determine the respective reciprocal lattice vector.
    ! The absolute value of every coordinate of the reciprocal lattice vector can be 1 maximally.
    allocate( kpq2kPrVec(3, kpts%nkpt, qpts%nkpt) )
    kpq2kPrVec = 0
    do iqpt = 1, qpts%nkpt
      do ikpt = 1, kpts%nkpt
        kpqNomin = 0
        do idir = 1, 3
          kpqNomin(idir) = nint( mod( kpts%bk(idir, ikpt) + qpts%bk(idir, iqpt), 1. ) * lcm )
          !kpqNomin(idir) = nint(lcm*kpts%bk(idir, ikpt) + lcm*qpts%bk(idir, iqpt))
          ! Is in 1st Brillouin zone
          if (abs(real(kpqNomin(idir)) / real(lcm) - (kpts%bk(idir, ikpt) + qpts%bk(idir, iqpt))) < 1e-5) then
          !if (kpqNomin(idir).lt.int(lcm)) then
            kpq2kPrVec(idir, ikpt, iqpt) = 0
          ! Has to be backfolded
          else
            kpq2kPrVec(idir, ikpt, iqpt) = -1
          end if
        end do
        if(.false.) then
          write(1005, '(i5, i5, 3(i5))') iqpt, ikpt, kpq2kPrVec(:, ikpt, iqpt)
        end if
        mapKpq2K(ikpt, iqpt) = mapK2Ind( kpqNomin(1), kpqNomin(2), kpqNomin(3) )
        !mapKpq2K(ikpt, iqpt) = mapK2Ind( mod(kpqNomin(1),int(lcm)), mod(kpqNomin(2),int(lcm)), mod(kpqNomin(3),int(lcm)) )
      end do
    end do

    ! For this found k-vector equals to k + q, determine the index and fill up array which connects the kpts indices of the k and q
    ! qpoint with the index of the k-vector equals to k + q.
    nkptShift = 0
    matchFound = .false.
    do iqpt = qpts%nkpt + 1, qpts%nkpt + nrAddQs
      do ikpt = 1, kpts%nkpt
        kpqTemp(:) = modulo1r( kpts%bk(:, ikpt) + qpts%bk(:, iqpt) )
        do ikptSh = 1, addQnkptis(iqpt - qpts%nkpt)
          if ( norm2( kpts%bk(:, kpts%nkpt + ikptSh + nkptShift) - kpqTemp(:) ) < 1e-7 ) then
            mapKpq2K(ikpt, iqpt) = kpts%nkpt + nkptShift + ikptSh
            matchFound = .true.
            exit
          end if
        end do
        if ( .not.matchFound ) then
          write (errorMessage, '(a,1x,i3)') 'no match for k+q-point', kpqTemp
          call juDFT_error( errorMessage, calledby='createkqMapArrays', hint='Please check k-points, q-points and mapKpq2K array!' )
        else
          matchFound = .false.
        end if
      end do
      nkptShift = nkptShift + addQnkptis(iqpt - qpts%nkpt)
    end do

    if (.false.) then
      ! Finds out which k' results when sending k to -k and then backfolding into 1st Brillouin zone.
      do ikpt = 1, kpts%nkpt
        do idir = 1, 3
          if ( kpts%bk(idir, ikpt)==0 ) then
            kpqTemp(idir) = kpts%bk(idir, ikpt)
          else
            kpqTemp(idir) = -kpts%bk(idir, ikpt) + 1
          end if
        end do ! idir
        mapK2mK(ikpt) = mapK2Ind(int(kpqTemp(1) * lcm), int(kpqTemp(2) * lcm), int(kpqTemp(3) * lcm))
      end do ! ikpt
    end if

  end subroutine createkqMapArrays

  subroutine ReadXCquantities( atoms, stars, dimens, lathar, vXC0IR, eXCIR, vXC0MT, eXCMT )

    use mod_juPhonUtils, only : fopen, fclose
    use m_types

    implicit none

    type(t_atoms),                  intent(in)  :: atoms
    type(t_stars),                  intent(in)  :: stars
    type(t_dimension),              intent(in)  :: dimens
    type(t_sphhar),                 intent(in)  :: lathar

    complex,           allocatable, intent(out) :: vXC0IR(:, :)
    complex,           allocatable, intent(out) :: eXCIR(:)
    real,              allocatable, intent(out) :: vXC0MT(:, :, :, :)
    real,              allocatable, intent(out) :: eXCMT(:, :, :)

    allocate( vXC0IR(stars%n3d, dimens%jspd), vXC0MT(atoms%jmtd, 0:lathar%nlhd, atoms%ntype, dimens%jspd) )
    allocate( eXCIR(stars%n3d), eXCMT( atoms%jmtd, 0:lathar%nlhd, atoms%ntype ) )
    vXC0IR(:, :) = cmplx(0., 0.)
    eXCIR(:) = cmplx(0., 0.)
    vXC0MT(:, :, :, :) = 0.
    eXCMT(:, :, : ) = 0.

    ! XC potential unwarped in the IR read in from FLEUR
    call fopen(1000, name='vxc0_uw_fleur', status='old', action='read', form='unformatted')
    read(1000) vXC0IR(:, :)
    call fclose(1000)

    ! XC energy density unwarped in the IR read in from FLEUR
    call fopen(1000, name='excpw_uw_fleur', status='old', action='read', form='unformatted')
    read(1000) eXCIR(:)
    call fclose(1000)

    ! XC potential in the MT read in from FLEUR
    call fopen(1000, name='v0MTFLEUR_xc', status='old', action='read', form='unformatted')
    read(1000) vXC0MT(:, :, :, :)
    call fclose(1000)

    ! XC energy density in the MT read in from FLEUR
    call fopen(1000, name='excr_fleur', status='old', action='read', form='unformatted')
    read(1000) eXCMT(:, :, :)
    call fclose(1000)

  end subroutine ReadXCquantities

end module m_jpInit
