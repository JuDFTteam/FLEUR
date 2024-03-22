module m_jpTestDynMatPulay

  implicit none

  contains

  ! tests only for Pulay contribution to the dynamical matrix, that can be performed before the Sternheimer SCC
  subroutine TestDynMatPulay( atoms, dimens, cell, sym, stars, lathar, input, kpts, qpts, enpara, usdus, Veff0, results, ngdp, &
      & memd_atom, logUnit, testGrMatElemPsiHepsPsiGaussTheoSw, testPsiHepsTildePsiSw, testR2orNotWfMtGradNgrNTensGrOvlsSw, testVarphiHepsVarphiSw, testGrPsiPsiMatElem, testIntVeff1Rho1Val, test3rdPulDynMatCancel, testIRIntegralSw, testIR3rdMatElemSw, testActionHgrPhiSw, clnu_atom, nmem_atom,&
      & mlh_atom, rho0IR, rho0MT, gdp, gBas, mapGbas, nv, ne, nobd, nRadFun, kpq2kPrVec, z0, rbas1, rbas2, El, eig, uuilon, duilon,&
      & ulouilopn, ilo2p, iloTable, kveclo, mapKpq2K, vEff0MTsh )


#include "cppmacro.h"
    use m_types

    implicit none

    ! Type parameters
    type(t_atoms),                  intent(in) :: atoms
    type(t_dimension),              intent(in) :: dimens
    type(t_cell),                   intent(in) :: cell
    type(t_sym),                    intent(in) :: sym
    type(t_stars),                  intent(in) :: stars
    type(t_sphhar),                 intent(in) :: lathar
    type(t_input),                  intent(in) :: input
    type(t_kpts),                   intent(in) :: kpts
    type(t_kpts),                   intent(in) :: qpts
    type(t_enpara),                 intent(in) :: enpara
    type(t_usdus),                  intent(in) :: usdus
    type(t_potential),              intent(in) :: Veff0
    type(t_results),                intent(in) :: results

    ! Scalar parameters
    integer,                        intent(in) :: ngdp
    integer,                        intent(in) :: memd_atom
    integer,                        intent(in) :: logUnit
    logical,                        intent(in) :: testGrMatElemPsiHepsPsiGaussTheoSw
    logical,                        intent(in) :: testPsiHepsTildePsiSw
    logical,                        intent(in) :: testR2orNotWfMtGradNgrNTensGrOvlsSw
    logical,                        intent(in) :: testVarphiHepsVarphiSw
    logical,                        intent(in) :: testGrPsiPsiMatElem
    logical,                        intent(in) :: testIntVeff1Rho1Val
    logical,                        intent(in) :: test3rdPulDynMatCancel
    logical,                        intent(in) :: testIRIntegralSw
    logical,                        intent(in) :: testIR3rdMatElemSw
    logical,                        intent(in) :: testActionHgrPhiSw

    ! Type parameters
    complex,                        intent(in) :: clnu_atom(:, 0:, :)
    integer,                        intent(in) :: nmem_atom(0:, :)
    integer,                        intent(in) :: mlh_atom(:, 0:, :)
    complex,                        intent(in) :: rho0IR(:,:)
    real,                           intent(in) :: rho0MT(:, 0:, :, :)
    integer,                        intent(in) :: gdp(:, :)
    integer,                        intent(in) :: gBas(:, :)
    integer,                        intent(in) :: mapGbas(:, :, :)
    integer,                        intent(in) :: nv(:, :)
    integer,                        intent(in) :: ne(:)
    integer,                        intent(in) :: nobd(:, :)
    integer,                        intent(in) :: nRadFun(0:, :)
    integer,                        intent(in) :: kpq2kPrVec(:, :, :)
    MCOMPLEX,                       intent(in) :: z0(:,:,:,:)
    real,                           intent(in) :: rbas1(:, :, 0:, :, :)
    real,                           intent(in) :: rbas2(:, :, 0:, :, :)
    real,                           intent(in) :: El(:, 0:, :, :)
    real,                           intent(in) :: eig(:, :, :)
    real,                           intent(in) :: uuilon(:, :)
    real,                           intent(in) :: duilon(:, :)
    real,                           intent(in) :: ulouilopn(:, :, :)
    integer,                        intent(in) :: ilo2p(:, :)
    integer,                        intent(in) :: iloTable(:, 0:, :)
    integer,                        intent(in) :: kveclo(:, :)
    integer,                        intent(in) :: mapKpq2K(:, :)
    complex,                        intent(in) :: vEff0MTsh(:, :, :, :)

    if ( testGrMatElemPsiHepsPsiGaussTheoSw .or. testPsiHepsTildePsiSw .or. testR2orNotWfMtGradNgrNTensGrOvlsSw .or. &
      & testIRIntegralSw .or. testIntVeff1Rho1Val .or. testIR3rdMatElemSw .or. testVarphiHepsVarphiSw .or. testGrPsiPsiMatElem &
      & .or. test3rdPulDynMatCancel .or. testActionHgrPhiSw ) then
      write(*, *)
      write(*, '(a)') 'Initiating dynamical matrix Pulay contribution test(s)...'
      write(*, '(a)') '---------------------------------------------------------'
      write(logUnit, *)
      write(logUnit, '(a)') 'Dynamical matrix Pulay contribution test(s)'
      write(logUnit, '(a)') '*******************************************'
      write(logUnit, *)
    else
      write(*, '(a)') '-----------------------------------------------------'
      write(*, '(a)') 'DISABLED dynamical matrix Pulay contribution test(s)!'
      write(*, '(a)') '-----------------------------------------------------'
      write(logUnit, '(a)')
      write(logUnit, '(a)') 'DISABLED dynamical matrix Pulay contribution test(s)!'
      write(logUnit, '(a)') '*****************************************************'
      return
    end if
    if (testGrMatElemPsiHepsPsiGaussTheoSw) then
      write(*, '(2x,a)') 'Performing TestGrMatElemPsiHepsPsiGaussTheo...'
      call TestGrMatElemPsiHepsPsiGaussTheo( atoms, lathar, stars, dimens, cell, input, kpts, enpara, usdus, sym, Veff0, ngdp, memd_atom, &
        & logUnit, clnu_atom, nmem_atom, mlh_atom, gdp, rho0IR, rho0MT, gBas, mapGbas, nv, nobd, z0, rbas1, rbas2, &
        & uuilon, duilon, ulouilopn, ilo2p, nRadFun, El, iloTable, kveclo, eig )
    else
      write(*, '(2x,a)') 'DISABLED TestGrMatElemPsiHepsPsiGaussTheo!'
      write(logUnit, '(a)') 'Testing nabla <Psi_kn|H - eps|Psi_kn> - oint dS Psi_kn^* (H - eps) Psi_kn = 0'
      write(logUnit, '(a)') '-----------------------------------------------------------------------------'
      write (logUnit, '(a)')'                                                                             |__ DISABLED!'
    end if

    if (testPsiHepsTildePsiSw) then
      write(*, '(2x,a)') 'Performing testPsiHepsTildePsi...'
      call testPsiHepsTildePsi( atoms, cell, lathar, dimens, kpts, results, Veff0, sym, usdus, enpara, input, clnu_atom, nmem_atom,&
        & mlh_atom, gBas, mapGbas, kveclo, nv, eig, nobd, z0, mapKpq2K, nRadFun, rbas1, rbas2, El, iloTable, uuilon, duilon, &
        & ulouilopn, ilo2p, logUnit, vEff0MTsh )
    else
      write(*, '(2x,a)') 'DISABLED testPsiHepsTildePsi!'
      write(logUnit, '(a)') 'Testing <BarPsi_kn|H - eps|Psi_kn> and <Psi_kn|H - eps|BarPsi_kn>'
      write(logUnit, '(a)') '-----------------------------------------------------------------'
      write (logUnit, '(a)')'                                                                 |__ DISABLED!'
    end if

    if (testR2orNotWfMtGradNgrNTensGrOvlsSw) then
      write(*, '(2x,a)') 'Performing TestR2orNotWfMtGradNgrNTensGrOvls...'
      call TestR2orNotWfMtGradNgrNTensGrOvls(atoms, usdus, kpts, cell, sym, dimens, results, lathar, logUnit, clnu_atom, nRadFun, z0, nv, kveclo, gBas, nobd, mapGbas, rbas1, rbas2, iloTable, mlh_atom, nmem_atom)
    else
      write(*, '(2x,a)') 'DISABLED TestR2orNotWfMtGradNgrNTensGrOvls!'
    write(logUnit, '(a)') 'Testing <grPsi|Psi>, <grGrtPsi|Psi> and <grPsi|grTPsi>'
    write(logUnit, '(a)') '------------------------------------------------------'
    write (logUnit, '(a)')'                                                      |__ DISABLED!'
    end if

    if (testIRIntegralSw) then
      !IR integral test and test of warping by using the reciprocal representation of the step function to compare with the routine for warping
      write(*, '(2x,a)') 'Performing TestIRIntegral...'
      call testIRIntegral(atoms, cell, stars, Veff0, dimens, kpts, input, results, ngdp, gdp, nv, mapGbas, gBas, nobd, z0, kpq2kPrVec, logUnit)
    else
      write(*, '(2x,a)') 'DISABLED TestIRIntegral!'
    write(logUnit, '(a)')   'Compare int dV grRho_IR grVeff_IR with <grPsi|grVeff|Psi>_IR'
    write(logUnit, '(a)')   '------------------------------------------------------------'
      write( logUnit, '(a)' ) '                                                           |__ DISABLED!'
    end if

    if (testIntVeff1Rho1Val) then
      ! Tests two representations of int Veff1 rho1 this is the mt test of the upper routine
      write(*, '(2x,a)') 'Performing TestDynMatPulInt...'
      call TestDynMatPulInt( atoms, sym, dimens, lathar, stars, input, usdus, Veff0, kpts, cell, results, ngdp, logUnit, rho0IR, rho0MT, Veff0%vr, gdp,     &
        & clnu_atom, nmem_atom, mlh_atom, nRadFun, rbas1, rbas2, nv, mapGbas, gBas, nobd, z0, kpq2kPrVec, El, kveclo, iloTable )
      !todo does the same as the routine above so we have to unite them
    !  call CalcGrVeffGrtRhoInt( atoms, cell, lathar, dimens, stars, Veff0, input, ngdp, memd_atom, clnu_atom, nmem_atom, mlh_atom, &
    !    & rho0IR, gdp, rho0MT )
  else
      write(*, '(2x,a)') 'DISABLED TestDynMatPulInt!'
      write(logUnit, '(a)')   'Compare of MT integral grRho grVeff with MT restricted <grPsi|grVeff|Psi>'
      write(logUnit, '(a)')     '-------------------------------------------------------------------------'
      write( logUnit, '(a)' ) '                                                                        |__ failed!'
    end if

    if (testIR3rdMatElemSw) then
      write(*, '(2x,a)') 'Performing TestIR3rdMatElem...'
      call testIR3rdMatElem( atoms, stars, lathar, dimens, usdus, sym, kpts, qpts, cell, input, logUnit, mlh_atom, clnu_atom, Veff0%vr(:, :, :, 1), rbas1, rbas2,&
      & ilo2p, nmem_atom, nRadFun, gBas, mapGbas, nv, kveclo, nobd, z0, El, iloTable, eig, kpq2kPrVec, Veff0%vpw, gdp, ngdp, Veff0%vpw_uw )
  else
      write(*, '(2x,a)') 'DISABLED TestIR3rdMatElem!'
      write(logUnit, '(a)')   'Compare of MT integral grRho grVeff with MT restricted <grPsi|grVeff|Psi>'
      write(logUnit, '(a)')     '-------------------------------------------------------------------------'
      write( logUnit, '(a)' ) '                                                                        |__ DISABLED!'
    end if

    if ( testVarphiHepsVarphiSw ) then
      write(*, '(2x,a)') 'Performing TestHeps1Phi...'
    call testHeps1Phi( atoms, lathar, dimens, enpara, usdus, input, sym, kpts, cell, logUnit, mlh_atom, clnu_atom,      &
      & Veff0%vr(:, :, :, 1), rbas1, rbas2, uuilon, duilon, ulouilopn, ilo2p, nmem_atom, nRadFun, gBas, mapGbas, nv, kveclo, nobd, z0, El,     &
      & iloTable, vEff0MTsh )
  else
      write(*, '(2x,a)') 'DISABLED TestHeps1Phi!'
    write( logUnit, '(a)' ) 'PsiHepsPsi test'
    write( logUnit, '(a)' ) '---------------'
      write( logUnit, '(a)' ) '              |__ DISABLED!'
    end if

    if (testGrPsiPsiMatElem) then
      ! lasts long!
      ! Tests whether 2 <grad Psi| Psi> = grad Rho for certain lm! There is a formula in Aarons thesis for that
      write(*, '(2x,a)') 'Performing TestGrVarphi...'
      call testGrVarphi( atoms, sym, dimens, lathar, stars, Veff0, kpts, cell, input, usdus, results, ngdp, logUnit, rho0IR, rho0Mt, gdp,      &
        & clnu_atom, nmem_atom, mlh_atom, nRadFun, rbas1, rbas2, nv, mapGbas, gBas, nobd, z0, kpq2kPrVec, El, kveclo, iloTable )
    else
      write(*, '(2x,a)') 'DISABLED TestGrVarphi!'
      write( logUnit, '(a)' ) 'Testing the lm components of grRho0MT with Gaunt coefficient'
      write( logUnit, '(a)' ) '------------------------------------------------------------'
      write( logUnit, '(a)' ) '                                                           |__ DISABLED!'
    end if

    if (test3rdPulDynMatCancel) then
      ! Tests whether for q = 0 all contributions except for <grad Psi | H - eps | gradT Psi> vanish
      write(*, '(2x,a)') 'Performing Test3rdBraKetHeps...'
      call Test3rdBraKetHeps( atoms, lathar, kpts, qpts, sym, dimens, cell, usdus, stars, Veff0, results, logUnit, nRadFun, nv, gBas, mapGbas, kveclo,     &
        & mapKpq2K, nobd, z0, iloTable, eig, kpq2kPrVec, El, Veff0%vr(:, :, :, 1), nmem_atom, mlh_atom, clnu_atom, rbas1, rbas2 )
    else
      write(*, '(2x,a)') 'DISABLED Test3rdBraKetHeps!'
      write( logUnit, '(a)' ) 'Check cancelling within 3rd braket in Pulay contribution to the dynamical matrix.'
      write( logUnit, '(a)' ) '---------------------------------------------------------------------------------'
      write( logUnit, '(a)' ) '                                                                                |__ DISABLED!'
    end if

    if (testActionHgrPhiSw) then
      ! tests the application of H on grVarphi in the same style as H on Varphi and it tests the projection onto the gradient of the basis function
      ! by using another algorithm which does the same but not modularized
      write(*, '(2x,a)') 'Performing TestActionHgrPhi...'
      call TestActionHgrPhi( atoms, dimens, lathar, Veff0, logUnit, rbas1(:, :, :, :, 1), rbas2(:, :, :, :, 1), nRadFun, El, mlh_atom, nmem_atom, clnu_atom )
    else
      write(*, '(2x,a)') 'DISABLED TestActionHgrPhi!'
      write(logUnit, '(a)') 'Test consistency of Hamiltonian action onto gradVarphi.'
      write(logUnit, '(a)') '-------------------------------------------------------'
      write (logUnit, '(a)')'                                                      |__ DISABLED!'
      write(logUnit, *)
    end if

  end subroutine TestDynMatPulay

  !>-------------------------------------------------------------------------------------------------------------------------------------
  !> @author
  !> Christian-Roman Gerhorst
  !>
  !> @brief
  !> Tests <grPsi|H - eps|Psi>_MT and <Psi|H - eps|grPsi>_MT using Gauss theorem
  !>
  !> @details
  !> According to Gauss theorem nabla <Psi_kn|H - eps|Psi_kn>_MT - oint_MT dS Psi_kn^* (H - eps_kn) Psi_kn = 0 holds and is checked
  !> here numerically. The deeper reason is that <grVarphi|H - eps|Varphi>_MT and <Varphi|H - eps|grVarphi>_MT should be checked;
  !> they are used in the matrix contributions for the Pulay dynamical matrix. Although this is only a test with the conventional
  !> matching coefficients it is better than nothing.
  !>
  !> @note
  !> If setting internal variable output true, the upper sum can be visualized for every k-point and band in the file grBKPsiHePsiGausRelKbnd.ssv. Debug output is hard-coded.
  !>
  !> @param[in] atoms        : Atoms type, see types.f90.
  !> @param[in] dimens       : Dimension type, see types.f90.
  !> @param[in] kpts         : K-points type, see types.f90.
  !> @param[in] input        : Input type, see types.f90.
  !> @param[out] enpara      : Energy parameter type, see types.f90.
  !> @param[out] eig         : Contains Kohn\--Sham eigenvalues.
  !> @param[out] kveclo      : Basis G-vectors of local orbitals.
  !> @param[out] nv          : Number of LAPW G-basis vectors for given k-point.
  !> @param[out] mapGbas     : For various k-points G-basis vectors occur more than once, thus they are only stored once in juPhon. This
  !>                           pointer array contains the right index for GbasVec array to "unfold" G-basis vectors again.
  !> @param[out] gBas        : G-basis vectors
  !> @param[out] ne          : Number of eigenvalues per k-point.
  !> @param[out] El          : Contains LAPW and LOs energy parameters.
  !> @param[in]  logUnit     : Unit number for juPhon.log.
  !>
  !> todo : sum passed varibales together for all tests
  !> todo : nicer output in files and logfile
  !> todo : rename for nicer variables
  !> todo : Check test for more than one atom
  !> todo : Nicer formatting
  !> todo : LOs not checked
  !>
  !>-------------------------------------------------------------------------------------------------------------------------------------
  !TODO Can also be performed with arbitrary abcofs, especially the abcofs used later in the test the condition should still hold
  !     and one would also check other channels accessed by the abcof1 for example.
  subroutine TestGrMatElemPsiHepsPsiGaussTheo( atoms, lathar, stars, dimens, cell, input, kpts, enpara, uds, sym, Veff0, &
      & ngdp, memd_atom, logUnit, clnu_atom, nmem_atom, mlh_atom, gdp, rho0IR, rho0MT, gBas, mapGbas, nv, nobd, z0, rbas1, &
      & rbas2, uuilon, duilon, ulouilopn, ilo2p, nRadFun, El, iloTable, kveclo, eig )

    use m_types, only : t_atoms, t_sphhar, t_stars, t_dimension, t_cell, t_input, t_kpts, t_enpara, t_usdus, t_sym, t_tlmplm, t_potential
    use m_jpPotDensHelper, only : calcIRdVxcKern, calcMTdVxcKern, ConvertStar2G
    use m_jpGrVeff0, only : GenGrVeff0
    use m_jpSternhHF, only : tlmplm4V, CalcVsumMT
    use m_jpConstants, only : iu, Tmatrix
    use m_abcof3
     
    use mod_juPhonUtils, only : Derivative, CalcChannelsGrFlpNat, Fopen, Fclose, calcGrR2FinLh
    use m_jpSetupDynMatHelper, only : CalcHnGrV0Varphi
    use m_jpSetupDynMatSF, only : CalcHGrVarphi, PrepareMTSurfIntDM, CalcSintMT, CalcSintKinEnergOvl
    use m_jpSetupDynMat, only : CalcVecBasfMatElems, CalcScalBasfMatElems
    use m_juDFT_NOstopNO, only : juDFT_warn

    implicit none

    ! Type parameter
    type(t_atoms),                  intent(in) :: atoms
    type(t_sphhar),                 intent(in) :: lathar
    type(t_stars),                  intent(in) :: stars
    type(t_dimension),              intent(in) :: dimens
    type(t_cell),                   intent(in) :: cell
    type(t_input),                  intent(in) :: input
    type(t_kpts),                   intent(in) :: kpts
    type(t_enpara),                 intent(in) :: enpara
    type(t_usdus),                  intent(in) :: uds
    type(t_sym),                    intent(in) :: sym
    type(t_potential),              intent(in) :: Veff0

    ! Scalar parameters
    integer,                        intent(in) :: ngdp
    integer,                        intent(in) :: memd_atom
    integer,                        intent(in) :: logUnit

    ! Array parameters
    complex,                        intent(in) :: clnu_atom(:, 0:, :) !member, 0,lathar, nat
    integer,                        intent(in) :: nmem_atom(0:, :) !0,lathar, nat
    integer,                        intent(in) :: mlh_atom(:, 0:, :) ! memd_atom, 0:lathar, nat
    integer,                        intent(in) :: gdp(:, :)
    complex,                        intent(in) :: rho0IR(:,:)
    real,                           intent(in) :: rho0MT(:, 0:, :, :)
    integer,                        intent(in) :: gBas(:, :)
    integer,                        intent(in) :: mapGbas(:, :, :)
    integer,                        intent(in) :: nv(:, :)
    integer,                        intent(in) :: nobd(:, :)
    MCOMPLEX,                       intent(in) :: z0(:, :, :, :) ! Attention: see zBar
    real,                           intent(in) :: rbas1(:, :, 0:, :, :)
    real,                           intent(in) :: rbas2(:, :, 0:, :, :)
    real,                           intent(in) :: uuilon(:, :)
    real,                           intent(in) :: duilon(:, :)
    real,                           intent(in) :: ulouilopn(:, :, :)
    integer,                        intent(in) :: ilo2p(:, :)
    integer,                        intent(in) :: nRadFun(0:, :)
    real,                           intent(in) :: El(:, 0:, :, :)
    integer,                        intent(in) :: iloTable(:, 0:, :)
    integer,                        intent(in) :: kveclo(:,:)
    real,                           intent(in) :: eig(:, :, :)

    ! Type parameters
    type(t_tlmplm)                             :: tdVx
    type(t_tlmplm)                             :: tdVy
    type(t_tlmplm)                             :: tdVz
    type(t_tlmplm)                             :: tdVx2
    type(t_tlmplm)                             :: tdVy2
    type(t_tlmplm)                             :: tdVz2
    type(od_inp)                               :: odi
    type(od_sym)                               :: ods

    ! Scalar parameters
    integer                                    :: itype
    integer                                    :: iG
    integer                                    :: ilh
    integer                                    :: idir
    integer                                    :: iatom
    integer                                    :: ieqat
    integer                                    :: oqn_l
    integer                                    :: mqn_m
    integer                                    :: lm
    integer                                    :: lm2
    integer                                    :: lmp
    integer                                    :: imesh
    integer                                    :: lm_pre
    integer                                    :: nmat
    integer                                    :: nrBnd
    logical                                    :: harSw =.true.
    logical                                    :: extSw =.true.
    logical                                    :: xcSw =.true.
    logical                                    :: vExtFull =.true.
    logical                                    :: testGoldstein =.false.
    integer                                    :: maxNrBnd
    integer                                    :: ikpt
    integer                                    :: iband
    integer                                    :: pMaxLocal
    integer                                    :: iradf
    integer                                    :: lmpMax
    integer                                    :: ptsym
    integer                                    :: imem
    integer                                    :: lmaxBra
    integer                                    :: mqn_m2Pr
    integer                                    :: jj
    integer                                    :: ii
    integer                                    :: iBas
    integer                                    :: nRadFunMax
    integer                                    :: chanMaxBra
    integer                                    :: chanMaxKet
    integer                                    :: mqn_m2PrBra
    integer                                    :: mqn_m2PrKet
    integer                                    :: coScale
    complex                                    :: gaussRelSum
    logical                                    :: output = .true.
    logical                                    :: passed = .true.

    ! Array parameters
    real,              allocatable             :: r2Rho0MT(:, :, :, :)
    complex,           allocatable             :: grRho0MT(:, :, :, :)
    complex,           allocatable             :: r2GrRho0MT(:, :, :, :)
    complex,           allocatable             :: grVxcIRKern(:)
    real,              allocatable             :: gaussWts(:) ! gaussian weights belonging to gausPts
    complex,           allocatable             :: ylm(:, : )
    real,              allocatable             :: dKernMTGPts(:, :, :)
    complex,           allocatable             :: grVeff0IR(:, :)
    complex,           allocatable             :: grVeff0MT(:, :, :, :)
    integer,           allocatable             :: nlo_atom(:)
    complex,           allocatable             :: ab0cofScl(:, :)
    complex,           allocatable             :: ab0cofSclNoIu(:, :, :)
    complex,           allocatable             :: a(:, :, :)
    complex,           allocatable             :: b(:, :, :)
    complex,           allocatable             :: bascof_lo(:, :, :, :, :)
    complex ,          allocatable             :: vSumMT(:, :, :)
    complex,           allocatable             :: sumVMTs(:, :, :, :)
    complex,           allocatable             :: sumVMTs2(:, :, :, :)
    integer,           allocatable             :: ngoprI(:)
    integer,           allocatable             :: lmpT(:)
    real,              allocatable             :: r2(:)
    real,              allocatable             :: varphi1(:, :, :)
    real,              allocatable             :: varphi2(:, :, :)
    real,              allocatable             :: delrVarphi1(:, :, :)
    real,              allocatable             :: delrVarphi2(:, :, :)
    integer,           allocatable             :: grVarphiChLout(:, :)
    integer,           allocatable             :: grVarphiChMout(:, :)
    real,              allocatable             :: grVarphiCh1(:, :, :, :)
    real,              allocatable             :: grVarphiCh2(:, :, :, :)
    complex,           allocatable             :: vEff0MtSpH(:, :)
    complex,           allocatable             :: hVarphi(:, :, :, :)
    complex,           allocatable             :: vEff0NsphGrVarphi(:, :, :, :, :)
    complex,           allocatable             :: hGrVarphi(:, :, :, :, :)
    complex,           allocatable             :: r2grVeff0SphVarphi(:, :, :, :, :)
    complex,           allocatable             :: r2grVeff0SphVarphiDummy(:, :, :, :, :)
    real,              allocatable             :: grVarphiVarphi(:, :, :, :)
    complex,           allocatable             :: grVarphiHVarphi(:, :, :, :)
    real,              allocatable             :: varphiGrVarphi(:, :, :, :)
    real,              allocatable             :: varphiVarphiDummy(:, :, :)
    complex,           allocatable             :: varphiHGrvarphi(:, :, :, :)
    complex,           allocatable             :: varphiGrVeffvarphi(:, :, :, :)
    real,              allocatable             :: varphiVarphi(:, :, :)
    complex,           allocatable             :: varphiHvarphi(:, :, :)
    complex,           allocatable             :: psiHepsGrPsi(:, :, :)
    complex,           allocatable             :: grPsiHepsPsi(:, :, :)
    complex,           allocatable             :: varphiHGrPsi(:)
    complex,           allocatable             :: varphiGrPsi(:)
    complex,           allocatable             :: grVarphiHPsi(:)
    complex,           allocatable             :: grVarphiPsi(:)
    complex,           allocatable             :: varphiGrVeffPsi(:)
    complex,           allocatable             :: hFullNoAbcofBK(:, :, :, :)
    complex,           allocatable             :: overlapNoAbcofBK(:, :, :, :)
    integer,           allocatable             :: muKet(:, :)
    integer,           allocatable             :: muBra(:, :)
    integer,           allocatable             :: lambdaKet(:, :)
    integer,           allocatable             :: lambdaBra(:, :)
    real,              allocatable             :: varphiBra1(:, :, :, :)
    real,              allocatable             :: varphiKet1(:, :, :, :)
    real,              allocatable             :: varphiBra2(:, :, :, :)
    real,              allocatable             :: varphiKet2(:, :, :, :)
    complex,           allocatable             :: overlap(:, :, :)
    complex,           allocatable             :: surfIntMT(:, :, :)
    complex,           allocatable             :: grRho0IR(:, :) ! Dummy quantity at the moment
    complex,           allocatable              :: rho0IRpw(:, :)
    complex,           allocatable              :: rho0MTsh(:, :, :, :)
    complex                                    :: grPsiHepsPsiNat(-1:1)
    complex                                    :: psiHepsGrPsiNat(-1:1)
    real                                       :: kExt(3)
    complex                                    :: psiGrVeffPsi(3)
    complex                                    :: grPsiPsiSum(3)
    complex                                    :: psiGrPsiSum(3)
    complex                                    :: surfIntSum(3)
    real                                       :: Gext(3)

    write(logUnit, *)
    write(logUnit, *)
    write(logUnit, '(a)') 'Testing nabla <Psi_kn|H - eps|Psi_kn> - oint dS Psi_kn^* (H - eps) Psi_kn = 0'
    write(logUnit, '(a)') '-----------------------------------------------------------------------------'

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
    ! Quantities for initialization
    allocate( lmpT(atoms%ntype) )
    lmpT(:) = 0
    do itype = 1, atoms%ntype
      lmpT(itype) = sum( [ ( (2 * oqn_l + 1)* nRadFun(oqn_l, itype), oqn_l = 0, atoms%lmax(itype) ) ] )
    end do ! itype
    lmpMax = maxval( lmpT(:) )

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
    ! Generate gradient of unperturbed density
    ! The factor r^2 has beeen divided out so that rho0MT is only rho0MT and not r^2 rho0MT as it is done in Fleur. The factor
    ! sqrt(4pi) for the zeroth component was already dividied out when constructing rho0MT in cdnmt routine in Fleur. Here to
    ! improve stability of the gradient routine we derive r2Rho0MT and divide out the r^2 again later. Doing so avoids the
    ! subtraction of small numbers close to the core.
    allocate( r2Rho0MT( atoms%jmtd, 0:lathar%nlhd, atoms%ntype, 1) )
    allocate( grRho0MT( atoms%jmtd, ( atoms%lmaxd + 1 )**2, atoms%nat, 3) )
    r2Rho0MT(:, :, :, :) = 0.
    grRho0MT(:, :, :, :) = cmplx(0., 0.)

    ! TODO Spin is missing
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

    ! Generate gradient of unperturbed effective potential, switches have been set in declaration. Start with xc-kernel calculation
    ! TODO Spin is missing
    call calcIRdVxcKern(stars, gdp, ngdp, rho0IR(:, 1), grVxcIRKern) ! add spin for this and next line
    call calcMTdVxcKern(atoms, dimens, lathar, rho0MT(:, :, :, 1), nmem_atom, clnu_atom, mlh_atom, gaussWts, ylm, dKernMTGPts)
    harSw = .true.
    extSw = .true.
    xcSw = .true.
    vExtFull = .true.
    testGoldstein = .false.
    call GenGrVeff0( atoms, cell, dimens, stars, ngdp, harSw, extSw, xcSw, gdp, rho0IRpw, rho0MTsh, grRho0IR, grRho0MT, gaussWts, ylm, &
      & dKernMTGPts, grVxcIRKern, testGoldstein, vExtFull, grVeff0IR, grVeff0MT ) ! add spin

    ! Allocations for basis function arrays
    maxNrBnd = maxval(nobd(:, :))
    nRadFunMax = maxval( nRadFun(:, :) )
    allocate( sumVMTs(atoms%jmtd, (atoms%lmaxd + 1 )**2, 3, atoms%nat) )
    allocate( sumVMTs2(atoms%jmtd, (atoms%lmaxd + 1 )**2, 3, atoms%nat) )
    allocate( varphiGrVeffVarphi(lmpMax, lmpMax, atoms%nat, 3), varphiHGrvarphi(lmpMax, lmpMax, atoms%ntype, -1:1) )
    allocate( varphiVarphiDummy(lmpMax, lmpMax, atoms%ntype))
    allocate( varphiGrVarphi(lmpMax, lmpMax, -1:1, atoms%ntype))
    allocate( grVarphiVarphi(lmpMax, lmpMax, -1:1, atoms%ntype), grVarphiHVarphi(lmpMax, lmpMax, -1:1, atoms%nat ) )
    allocate( r2(atoms%jmtd) )
    allocate( varphi1(atoms%jmtd, nRadFunMax, 0:atoms%lmaxd), varphi2(atoms%jmtd, nRadFunMax, 0:atoms%lmaxd) )
    allocate( vEff0MtSpH( atoms%jmtd, 0:dimens%lmd), vEff0NsphGrVarphi(2, atoms%jmtd, 0:(atoms%lmaxd + 3)**2 - 1, lmpMax, -1:1), &
            & r2grVeff0SphVarphi(2, atoms%jmtd, 0:(atoms%lmaxd + 3)**2 - 1, lmpMax, 3) )
    allocate( r2grVeff0SphVarphiDummy(2, atoms%jmtd, 0:(atoms%lmaxd + 3)**2 - 1, lmpMax, 3) )
    allocate( hGrVarphi(2, atoms%jmtd, 0:(atoms%lmaxd + 3)**2 - 1, lmpMax, -1:1))
    allocate( hVarphi(2, atoms%jmtd, 0:(atoms%lmaxd + 3)**2 - 1, lmpMax)  )
    allocate( delrVarphi1( atoms%jmtd, nRadFunMax, 0:atoms%lmaxd), delrVarphi2( atoms%jmtd, nRadFunMax, 0:atoms%lmaxd) )
    allocate( grVarphiCh1(atoms%jmtd, 2, lmpMax, -1:1), grVarphiCh2(atoms%jmtd, 2, lmpMax, -1:1), &
            & grVarphiChLout(2, 0:atoms%lmaxd), grVarphiChMout(-atoms%lmaxd:atoms%lmaxd, -1:1) )
    allocate( hFullNoAbcofBK(lmpMax, lmpMax, 3, atoms%nat), overlapNoAbcofBK(lmpMax, lmpMax, 3, atoms%ntype) )
    allocate( muKet(-atoms%lmaxd:atoms%lmaxd, -1:1), muBra(-atoms%lmaxd:atoms%lmaxd, -1:1) )
    allocate( lambdaKet(2, 0:atoms%lmaxd), lambdaBra(2, 0:atoms%lmaxd) )
    allocate( varphiBra1(atoms%jmtd, 2, lmpMax, -1:1), varphiBra2(atoms%jmtd, 2, lmpMax, -1:1), &
            & varphiKet1(atoms%jmtd, 2, lmpMax, -1:1), varphiKet2(atoms%jmtd, 2, lmpMax, -1:1) )

    hFullNoAbcofBK(:, :, :, :) = cmplx(0., 0.)
    overlapNoAbcofBK(:, :, :, :) = cmplx(0., 0.)
    varphiGrVeffVarphi(:, :, :, :) = cmplx(0., 0.)
    varphiHGrvarphi(:, :, :, :) = cmplx(0., 0.)
    varphiVarphiDummy(:, :, :) = cmplx(0., 0.)
    varphiGrVarphi(:, :, :, :) = 0.
    grVarphiVarphi(:, :, :, :) = cmplx(0., 0.)
    grVarphiHVarphi(:, :, :, :) = cmplx(0., 0.)
    r2(:) = 0.

    ! Fill the sum variable analogously to the Sternheimer algorithm, only minus the MT Weinert gradient of the unperturbed
    ! effective potential is used
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
                sumVMTs(imesh, lm, idir, iatom) = -grVeff0MT(imesh, lm, idir, iatom)
                sumVMTs2(imesh, lm2, idir, iatom) = (-1)**(mqn_m) * conjg(sumVMTs(imesh, lm, idir, iatom))
              end do
            end do
          end do
        end do
      end do
    end do

    ! Construct analogously to FLEUR the upper and lower diagonal of the non-spherical potential integrals
    ! todo Spin is missing
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

    deallocate(sumVMTs, sumVMTs2)

    iatom = 0
    do itype = 1, atoms%ntype

      r2(:) = 0.
      ! Jacobi determinant of radial part
      do imesh = 1, atoms%jri(itype)
       r2(imesh) = atoms%rmsh(imesh, itype) * atoms%rmsh(imesh, itype)
      end do ! imesh

      do ieqat = 1, atoms%neq(itype)
        iatom = iatom + 1

        varphi1(:, :, :) = 0.
        varphi2(:, :, :) = 0.
        delrVarphi1(:, :, :) = 0.
        delrVarphi2(:, :, :) = 0.
        do oqn_l = 0, atoms%lmax(itype)
          do iradf = 1, nRadFun(oqn_l, itype)
            do imesh = 1, atoms%jri(itype)
              ! In Fleur the radial solutions u_lp are multiplied by a factor r by default to avoid an additional multiplication of the
              ! Jacobi determinant r^2 in radial integrals given spherical coordinates.
              varphi1(imesh, iradf, oqn_l) = rbas1(imesh, iradf, oqn_l, itype, 1) / atoms%rmsh(imesh, itype)
              varphi2(imesh, iradf, oqn_l) = rbas2(imesh, iradf, oqn_l, itype, 1) / atoms%rmsh(imesh, itype)
            end do ! imesh
            ! Precalculate partial derivatives of varphis in r-direction since it is needed twice
            call Derivative( varphi1(1:atoms%jri(itype), iradf, oqn_l), itype, atoms, delrVarphi1(1:atoms%jri(itype), iradf, oqn_l) )
            call Derivative( varphi2(1:atoms%jri(itype), iradf, oqn_l), itype, atoms, delrVarphi2(1:atoms%jri(itype), iradf, oqn_l) )
          end do ! iradf
        end do ! oqn_l

        ! Calculate the muffin-tin gradients of the radial basis functions multiplied with spherical harmonics
        grVarphiChLout(:, :) = 0
        grVarphiChMout(:, :) = 0
        grVarphiCh1(:, :, :, :) = 0.
        grVarphiCh2(:, :, :, :) = 0.
        call CalcChannelsGrFlpNat( atoms, itype, nRadFun, varphi1, varphi2, delrVarphi1, delrVarphi2, grVarphiChLout, grVarphiChMout, &
                                                                                                        & grVarphiCh1, grVarphiCh2 )

        ! Expand the coefficients of the lattice-harmonic given potential into spherical-harmonic coefficients for the given atom.
        ! TODO PUT THIS INTO calcHnGrV0VArphi
        vEff0MtSpH(:, :) = cmplx(0.0, 0.0)
        ptsym = atoms%ntypsy(iatom)
        do ilh = 0, lathar%nlh(ptsym)
          oqn_l = lathar%llh(ilh, ptsym)
          lm_pre = oqn_l * (oqn_l + 1)
          do imem = 1, nmem_atom(ilh, iatom)
            mqn_m = mlh_atom(imem, ilh, iatom)
            lm = lm_pre + mqn_m
            !todo one could only evaluate the vEff0MtSpH which have a contribution, i.e. the oqn_l and mqn_m which are in llh and mlh_atom
            ! TODO Spin is missing
            ! maybe construct a pointer and run only over them to make it memory efficient.
            do imesh = 1, atoms%jri(itype)
              vEff0MtSpH(imesh, lm) = vEff0MtSpH(imesh, lm) + Veff0%vr(imesh, ilh, itype, 1) * clnu_atom(imem, ilh, iatom)
            end do ! imesh
          end do ! imem
        end do ! ilh


        ! Calculate the action of the full Hamiltonian on a basis function varphi, calculate the product of a non-spherical
        ! effective potential, a ket basis-function and a Gaunt coefficient expecting the bra up to lmax + 1 and calculate the
        ! product of the radial Jacobi determiant r^2 times the spherical part of the effective potential times a MT basis function
        ! times a Gaunt coefficient expecting a bra up to lmax + 1
        ! TODO Spin is missing
        vEff0NsphGrVarphi(:, :, :, :, :) = cmplx(0., 0.)
        r2grVeff0SphVarphi(:, :, :, :, :) = cmplx(0., 0.)
        hVarphi(:, :, :, :) = cmplx(0., 0.)
        call CalcHnGrV0Varphi( atoms, lathar, itype, iatom, lmpMax, El, varphi1, varphi2, nRadFun, vEff0MtSpH, Veff0%vr(:, :, :, 1), clnu_atom, &
          & nmem_atom, mlh_atom, grVarphiCh1, grVarphiCh2, grVarphiChLout, grVarphiChMout, hVarphi, vEff0NsphGrVarphi, r2grVeff0SphVarphi, r2grVeff0SphVarphiDummy)
        deallocate(vEff0MtSpH, r2grVeff0SphVarphiDummy)

        ! Call routine to calculate H grVarphi
        lmaxBra = atoms%lmax(itype)
        hGrVarphi(:, :, :, :, :) = cmplx(0., 0.)
        do mqn_m2Pr = -1, 1
          call CalcHGrVarphi( atoms, itype, mqn_m2Pr, lmpMax, lmaxBra, grVarphiChMout, nRadFun, grVarPhiCh1, grVarPhiCh2,        &
                                                      & grVarphiChLout, vEff0NsphGrVarphi, El, lmpT, hGrVarphi )
        end do ! mqn_m2Pr
        deallocate(vEff0NsphGrVarphi)

        ! Calculate <grVarphi|varphi> and <grVarphi|H|varphi> using the action of the full Hamiltonian on varphi
        call CalcVecBasfMatElems( atoms, itype, 2, nRadFun, r2, grVarphiChLout, grVarphiChMout, varphi1, varphi2, &
           & grVarPhiCh1, grVarphiCh2, hVarphi, grVarphiVarphi(:, :, :, itype), grVarphiHVarphi(:, :, :, iatom) )

        ! Transpose <grVarphi|Varphi> to <varphi|grVarphi> for testing reasons.
        ! NOTE: Actually this is not necessary in the production code
        do mqn_m2Pr = -1, 1
          do jj = 1, lmpMax
            do ii = 1, lmpMax
              varphiGrVarphi(ii, jj, mqn_m2Pr, itype) = grVarphiVarphi(jj, ii, mqn_m2Pr, itype)
            end do ! ii
          end do ! jj

          ! Calculate <varphi| H | grVarphi> using the action of H onto grVarphi
          call CalcScalBasfMatElems( atoms, itype, iatom, nRadFun, r2, varphi1, varphi2, hGrVarphi(:, :, :, :, mqn_m2Pr), &
            & varphiVarphiDummy, varphiHGrvarphi(:, :, :, mqn_m2Pr) )
        end do ! mqn_m2Pr

        ! Calculate <varphi|Veff0Sph|varphi>, remember that the Jacobi determinant is already present, therefore, we do not need to
        ! multiply it anymore which is why r2 is set onto 1
        r2(:) = 1.
        do idir = 1, 3
          call CalcScalBasfMatElems( atoms, itype, iatom, nRadFun, r2, varphi1, varphi2, r2grVeff0SphVarphi(:, :, :, :, idir), &
            & varphiVarphiDummy, varphigrVeffvarphi(:, :, :, idir) )
        end do ! idir

      end do ! ieqat

      ! Calculate the MT surface integral using the action of H onto varphi and the general routine for that
      chanMaxBra = 1
      chanMaxKet = 1
      lmaxBra = atoms%lmax(itype)
      muBra(:, :) = 0
      muKet(:, :) = 0
      mqn_m2PrBra = -1
      mqn_m2PrKet = -1
      varphiBra1 = 0.
      varphiKet1 = 0.
      varphiBra2 = 0.
      varphiKet2 = 0.
      do mqn_m = -atoms%lmax(itype), atoms%lmax(itype)
        muKet(mqn_m, -1) = mqn_m
        muBra(mqn_m, -1) = mqn_m
      end do ! mqn_m
      do oqn_l = 0, atoms%lmax(itype)
        lambdaKet(1, oqn_l) = oqn_l
        lambdaBra(1, oqn_l) = oqn_l
      end do ! oqn_l
      lmp = 0
      do oqn_l = 0, atoms%lmax(itype)
        do mqn_m = -oqn_l, oqn_l
          do iradf = 1, nRadFun(oqn_l, itype)
            lmp = lmp + 1
            varphiBra1(atoms%jri(itype), 1, lmp, -1) = varphi1(atoms%jri(itype), iradf, oqn_l )
            varphiKet1(atoms%jri(itype), 1, lmp, -1) = varphi1(atoms%jri(itype), iradf, oqn_l )
            varphiBra2(atoms%jri(itype), 1, lmp, -1) = varphi2(atoms%jri(itype), iradf, oqn_l )
            varphiKet2(atoms%jri(itype), 1, lmp, -1) = varphi2(atoms%jri(itype), iradf, oqn_l )
          end do ! iradf
        end do ! mqn_m
      end do ! oqn_l

      call PrepareMTSurfIntDM( atoms, itype, chanMaxBra, chanMaxKet, lmaxBra, lmpMax, mqn_m2PrBra, mqn_m2PrKet, muKet, muBra,     &
        & lambdaKet, lambdaBra, nRadFun, varphiKet1, varphiKet2, varphiBra1, varphiBra2, hVarphi, atoms%rmt(itype)**2, hFullNoAbcofBK, overlapNoAbcofBK )

    end do ! itype

    ! Deallocate helping arrays for the basis function arrays
    deallocate(varphiVarphiDummy)
    deallocate( r2grVeff0SphVarphi, hGrVarphi, hVarphi, delrVarphi1, delrVarphi2, grVarphiCh1, grVarphiCh2, grVarphiChLout, grVarphiChMout )
    deallocate( varphi1, varphi2, muKet, muBra, lambdaKet, lambdaBra, varphiBra1, varphiBra2, varphiKet1, varphiKet2 )

    ! Allocate arrays to calcualte the final required quantities
    allocate( ngoprI(atoms%nat) )
    allocate( vSumMT(maxNrBnd, maxNrBnd, 3) )
    allocate( ab0cofScl(lmpMax, maxNrBnd))
    allocate( ab0cofSclNoIu(maxNrBnd, lmpMax, atoms%nat))
    allocate( a( dimens%nvd, 0:dimens%lmd, atoms%nat), b(dimens%nvd, 0:dimens%lmd, atoms%nat), &
            & bascof_lo(3, -atoms%llod:atoms%llod, 4 * atoms%llod + 2, atoms%nlod, atoms%nat) )
    allocate( overlap(maxNrBnd, 3, atoms%ntype), surfIntMT(maxNrBnd, 3, atoms%ntype) )
    allocate( varphiGrVeffPsi(lmpMax))
    allocate( psiHepsGrPsi(3, maxNrBnd, atoms%nat), grPsiHepsPsi(3, maxNrBnd, atoms%nat))
    allocate( varphiHGrPsi(lmpMax), varphiGrPsi(lmpMax), grVarphiHPsi(lmpMax), grVarphiPsi(lmpMax) )
    psiHepsGrPsi(:, :, :) = cmplx(0., 0.)
    grPsiHepsPsi(:, :, :) = cmplx(0., 0.)
    varphiHGrPsi(:) = cmplx(0., 0.)
    varphiGrPsi(:) = 0.
    grVarphiHPsi(:) = cmplx(0., 0.)
    grVarphiPsi(:) = 0.
    varphiGrVeffPsi(:) = cmplx(0., 0.)
    ngoprI(:) = 1
    grPsiPsiSum(:) = cmplx(0., 0.)
    psiGrPsiSum(:) = cmplx(0., 0.)
    surfIntSum(:) = cmplx(0., 0.)

    if (output) then
      call Fopen(1000, name='grBKPsiHePsiGausRelKbnd.ssv', status='replace', action='write', form='formatted')
      write(1000, '(a)') 'ikpt     kx  ky kz   direction    band  Sum of MT volume integrals with gradients and surface integral without gradients'
    end if

    do ikpt = 1, kpts%nkpt

      ! Calculate once the small abcofs which can be used twice for the different abcofs
      nmat = nv(1, ikpt) + atoms%nlotot
      a(:, :, :) = cmplx(0.0, 0.0)
      b(:, :, :) = cmplx(0.0, 0.0)
      bascof_lo(:, :, :, :, :) = cmplx(0.0, 0.0)
      ! TODO : Spin missing
      call abcof3( atoms%lmaxd, atoms%ntype, atoms%nat, sym%nop, dimens%nvd, dimens%jspd, 1, dimens%lmd, dimens%nbasfcn, &
        & atoms%llod, atoms%nlod, atoms%nlotot, sym%invtab, atoms%ntype, sym%mrot, ngoprI, atoms%taual, atoms%neq, atoms%lmax, &
        & atoms%rmt, cell%omtil, cell%bmat, cell%bbmat, kpts%bk(:, ikpt), gBas(1, mapGbas(:nv(1, ikpt), ikpt, 1)), &
        & gbas(2, mapGbas(:nv(1, ikpt), ikpt, 1)), gbas(3, mapGbas(:nv(1, ikpt), ikpt, 1)), nv(:, ikpt), nmat, &
        & uds%us, uds%dus, uds%uds, uds%duds, uds%ddn, atoms%invsat, sym%invsatnr, uds%ulos, uds%uulon, uds%dulon, &
        & uds%dulos, atoms%llo, atoms%nlo, atoms%l_dulo, atoms%lapw_l, kveclo(:,ikpt), odi, ods, a, b, bascof_lo )

      ! Calcuate the abcofs which would be the result of the abcof routine from FLEUR times i^l.
      ab0cofSclNoIu(:, :, :) = cmplx(0., 0.)
      psiHepsGrPsi(:, :, :) = cmplx(0., 0.)
      grPsiHepsPsi(:, :, :) = cmplx(0., 0.)
      iatom = 0
      do itype = 1, atoms%ntype
        do ieqat = 1, atoms%neq(itype)
          iatom = iatom + 1
          ab0cofScl(:, :) = cmplx(0., 0.)
          do iband = 1, nobd(ikpt, 1)
            lmp = 0
            lm  = 0
            do oqn_l = 0, atoms%lmax(itype)
              do mqn_m = -oqn_l, oqn_l
                pMaxLocal = nRadFun(oqn_l, itype)
                ! p = 1
                ab0cofScl(lmp + 1, iband) = iu**oqn_l * dot_product( conjg(z0(:nv(1, ikpt), iband, ikpt, 1)), a(:nv(1, ikpt), lm, iatom) )
                ! p = 2
                ab0cofScl(lmp + 2, iband) = iu**oqn_l * dot_product( conjg(z0(:nv(1, ikpt), iband, ikpt, 1)), b(:nv(1, ikpt), lm, iatom) )
                ! Add LO contributions
                do iradf = 3, pMaxLocal
                  ! p = 1
                  ab0cofScl(lmp + 1, iband) = ab0cofScl(lmp + 1, iband) + iu**oqn_l * &
                    & dot_product( conjg(z0(nv(1, ikpt) + 1:nv(1, ikpt) + atoms%nlotot, iband, ikpt, 1)), &
                                                      & bascof_lo(1, mqn_m, :atoms%nlotot, iloTable(iradf, oqn_l, itype), iatom) )
                  ! p = 2
                  ab0cofScl(lmp + 2, iband) = ab0cofScl(lmp + 2, iband) + iu**oqn_l * &
                    & dot_product( conjg(z0(nv(1, ikpt) + 1:nv(1, ikpt) + atoms%nlotot, iband, ikpt, 1)), &
                                                      & bascof_lo(2, mqn_m, :atoms%nlotot, iloTable(iradf, oqn_l, itype), iatom) )
                  ! 2 < p < LOs for that l and that atom type
                  ab0cofScl(lmp + iradf, iband) = iu**oqn_l * dot_product( conjg(z0(nv(1, ikpt) + 1:nv(1, ikpt) + atoms%nlotot, iband, ikpt, 1)), &
                                                      & bascof_lo(3, mqn_m, :atoms%nlotot, iloTable(iradf, oqn_l, itype), iatom) )
                end do ! iradf

                ! This is a precalculation of the 1st and 3rd line of A.50. Actually the resulting lmp should have been primed or at
                ! least the p as lm1Pr = lm. But for sake of performance we place it here.
                ! sake of performance
                lm = lm + 1
                lmp = lmp + pMaxLocal
              end do ! mqn_m
            end do !oqn_l

            ! Multiply the matching coefficients times the wave function coefficients to the basis function arrays.
            grPsiHepsPsiNat(:) = cmplx(0., 0.)
            psiHepsGrPsiNat(:) = cmplx(0., 0.)
            psiGrVeffPsi(:) = cmplx(0., 0.)
            do mqn_m2Pr = -1, 1
              ! <Psi|GrVeffSph|Psi> which lacks as contribution to <varphi|H|grVarphi>
              varphiGrVeffPsi(:) = cmplx(0., 0.)
              varphiGrVeffPsi(1:lmpT(itype)) = matmul(varphigrVeffvarphi(1:lmpT(itype), 1:lmpT(itype), iatom, mqn_m2Pr + 2), ab0cofScl(1:lmpT(itype), iband))
              psiGrVeffPsi(mqn_m2Pr + 2) = dot_product(ab0cofScl(1:lmpT(itype), iband), varphiGrVeffPsi(1:lmpT(itype)))

              ! <grVarphi|H|Psi> and <grVarphi|varphi>
              grVarphiPsi(:) = cmplx(0., 0.)
              grVarphiHPsi(:) = cmplx(0., 0.)
              grVarphiPsi(1:lmpT(itype)) = matmul(grVarphiVarphi(1:lmpT(itype), 1:lmpT(itype), mqn_m2Pr, itype), ab0cofScl(1:lmpT(itype), iband))
              grVarphiHPsi(1:lmpT(itype)) = matmul(grVarphiHVarphi(1:lmpT(itype), 1:lmpT(itype), mqn_m2Pr, iatom), ab0cofScl(1:lmpT(itype), iband))

              ! <varphi|GrPsi>
              ! <varphi|H|GrPsi>
              varphiGrPsi(:) = cmplx(0., 0.)
              varphiHGrPsi(:) = cmplx(0., 0.)
              varphiGrPsi(1:lmpT(itype)) = matmul(varphiGrvarphi(1:lmpT(itype), 1:lmpT(itype), mqn_m2Pr, itype), ab0cofScl(1:lmpT(itype), iband))
              varphiHGrPsi(1:lmpT(itype)) = matmul(varphiHGrvarphi(1:lmpT(itype), 1:lmpT(itype), itype, mqn_m2Pr), ab0cofScl(1:lmpT(itype), iband))

              ! <grPsi|H - eps|Psi> and <psi|H - eps|grPsi> both following and the origin arrays are given in natural coordinates
              grPsiHepsPsiNat(mqn_m2Pr) = dot_product(ab0cofScl(1:lmpT(itype), iband), grVarphiHPsi(1:lmpT(itype))) &
                      &   - eig(iband, ikpt, 1) * dot_product(ab0cofScl(1:lmpT(itype), iband), grVarphiPsi(1:lmpT(itype)))
              psiHepsGrPsiNat(mqn_m2Pr) = dot_product(ab0cofScl(1:lmpT(itype), iband), varphiHGrPsi(1:lmpT(itype))) &
                      &   - eig(iband, ikpt, 1) * dot_product(ab0cofScl(1:lmpT(itype), iband), varphiGrPsi(1:lmpT(itype)))
            end do ! mqn_m2Pr

            ! Transforming <grPsi|H - eps|Psi> and <psi|H - eps|grPsi> to cartesian coordinates
            grPsiHepsPsi(1:3, iband, iatom) = matmul(conjg(Tmatrix(1:3, 1:3)), grPsiHepsPsiNat(-1:1))
            psiHepsGrPsi(1:3, iband, iatom) = matmul(Tmatrix(1:3, 1:3), psiHepsGrPsiNat(-1:1)) - psiGrVeffPsi(1:3)

          end do ! iband

          ! Calculating the large matching coefficients as they would have been the output in abcof from FLEUR (no i^l multiplied)
          lmp = 0
          lm  = 0
          do oqn_l = 0, atoms%lmax(itype)
            do mqn_m = -oqn_l, oqn_l
              do iband = 1, nobd(ikpt, 1)
                pMaxLocal = nRadFun(oqn_l, itype)
                ! p = 1
                ab0cofSclNoIu(iband, lmp + 1, iatom) = dot_product( conjg(z0(:nv(1, ikpt), iband, ikpt, 1)), a(:nv(1, ikpt), lm, iatom) )
                ! p = 2
                ab0cofSclNoIu(iband, lmp + 2, iatom) = dot_product( conjg(z0(:nv(1, ikpt), iband, ikpt, 1)), b(:nv(1, ikpt), lm, iatom) )
                ! Add LO contributions
                do iradf = 3, pMaxLocal
                  ! p = 1
                  ab0cofSclNoIu(iband, lmp + 1, iatom) = ab0cofSclNoIu(iband, lmp + 1, iatom) +  &
                    & dot_product( conjg(z0(nv(1, ikpt) + 1:nv(1, ikpt) + atoms%nlotot, iband, ikpt, 1)), &
                                                      & bascof_lo(1, mqn_m, :atoms%nlotot, iloTable(iradf, oqn_l, itype), iatom) )
                  ! p = 2
                  ab0cofSclNoIu(iband, lmp + 2, iatom) = ab0cofSclNoIu(iband, lmp + 2, iatom) +  &
                    & dot_product( conjg(z0(nv(1, ikpt) + 1:nv(1, ikpt) + atoms%nlotot, iband, ikpt, 1)), &
                                                      & bascof_lo(2, mqn_m, :atoms%nlotot, iloTable(iradf, oqn_l, itype), iatom) )
                  ! 2 < p < LOs for that l and that atom type
                  ab0cofSclNoIu(iband, lmp + iradf, iatom) =  dot_product( conjg(z0(nv(1, ikpt) + 1:nv(1, ikpt) + atoms%nlotot, iband, ikpt, 1)), &
                                                      & bascof_lo(3, mqn_m, :atoms%nlotot, iloTable(iradf, oqn_l, itype), iatom) )
                end do ! iradf
              end do ! iband

              ! This is a precalculation of the 1st and 3rd line of A.50. Actually the resulting lmp should have been primed or at
              ! least the p as lm1Pr = lm. But for sake of performance we place it here.
              ! sake of performance
              lm = lm + 1
              lmp = lmp + pMaxLocal
            end do ! mqn_m
          end do !oqn_l

          ! todo probably wrong for more than one atom
          ! Calculate the MT surface integral with the integrand e_r Psi^* (H - eps) Psi
          ! TODO Spin missing
          overlap(:, :, :) = cmplx(0., 0.)
          surfIntMT(:, :, :) = cmplx(0., 0.)
          do idir = 1, 3
            call CalcSintMT( ikpt, ikpt, idir, iatom, itype, lmpMax, nobd, hFullNoAbcofBK, overlapNoAbcofBK, ab0cofScl, ab0cofScl,   &
                                                                              & overlap(:, idir, itype), surfIntMT(:, idir, itype) )
          end do ! idir
        end do ! ieqat
      end do ! itype

      ! Us the Sternheimer routine to calculate <Psi_pk'|V_eff^(1)|Psi_nk> to calculate <Psi_nk|grad V_eff^(0)|Psi_nk>
      !TODO Here, we have to think, when dealing with more than one atom, because this routine differentiates between displaced atoms and all atoms, this routine intrinsically performs a loop over all atoms, we have to think whether this is correct
      !TODO Spin missing
      vSumMT = cmplx(0.0, 0.0)
      call calcVsumMT( atoms, tdVx, tdVx2, ikpt, ikpt, nobd(:, 1), nobd, ab0cofSclNoIu, ab0cofSclNoIu, nRadFun, iloTable, nlo_atom, vSumMT(:, :, 1) )
      call calcVsumMT( atoms, tdVy, tdVy2, ikpt, ikpt, nobd(:, 1), nobd, ab0cofSclNoIu, ab0cofSclNoIu, nRadFun, iloTable, nlo_atom, vSumMT(:, :, 2) )
      call calcVsumMT( atoms, tdVz, tdVz2, ikpt, ikpt, nobd(:, 1), nobd, ab0cofSclNoIu, ab0cofSclNoIu, nRadFun, iloTable, nlo_atom, vSumMT(:, :, 3) )

      ! Sum up contributions
      iatom = 0
      do itype = 1, atoms%ntype
        do ieqat = 1, atoms%neq(itype)
          iatom = iatom + 1
          do idir = 1, 3
            do iband = 1, nobd(ikpt, 1)
              ! Todo The handling of vSumMT is incorrect for more than one atom
              ! TODO Spin missing
              gaussRelSum = vSumMT(iband, iband, idir) + surfIntMT(iband, idir, iatom) - eig(iband, ikpt, 1) * overlap(iband, idir, iatom) - grPsiHepsPsi(idir, iband, iatom) - psiHepsGrPsi(idir, iband, iatom) + surfIntMT(iband, idir, iatom)
              if (abs(real(gaussRelSum)) > 1e-7) then
                write(logUnit, '(a, 3i5)') 'Error: Relation contains real contributions at (k-point index   band index   displacement direction)', ikpt, iband, idir
                passed = passed .and. .false.
              end if

              if (output) then
                write(1000, '(i8,3f15.8,2i8,2f15.8)') ikpt, kpts%bk(1, ikpt), kpts%bk(2, ikpt), kpts%bk(3, ikpt), idir, iband, gaussRelSum
              end if
              if (.false.) then
                write(4000, '(3i8,2f15.8)') idir, ikpt, iband, surfIntMT(iband, idir, iatom)
                write(4001, '(3i8,2f15.8)') idir, ikpt, iband, overlap(iband, idir, iatom)
                write(4002, '(3i8,2f15.8)') ikpt, idir, iband, vSumMT(iband, iband, idir)
                write(4003, '(3i8,2f15.8)') ikpt, idir, iband, grPsiHepsPsi(idir, iband, iatom)
                write(4004, '(3i8,2f15.8)') ikpt, idir, iband, psiHepsGrPsi(idir, iband, iatom)
              end if

              ! Sum up the matrix elements which have imaginary contributions, at least in the sum over the k-points and the bands,
              ! they should give zero. The symmetric k-points have no contribution here, but it seems that we have to add pairs of
              ! unsymmetric k-points with its partner, i.e. 0 0 0.25 with 0 0 0.75 for example. This does not work completely,
              ! probably because we have degenerated s-bands for Neon.
              ! But in sum, everything is okay, furthermore, such terms only occur under sums of k-points and bands
              ! If we only have potential instead of using the E_l, then grPsiHepsPsi and psiHepsGrPsi cancel
              ! Also in the surface integral only when using any term of spehrical Hamiltonian does not cancel, potentials solely do
              grPsiPsiSum(idir) = grPsiPsiSum(idir) + grPsiHepsPsi(idir, iband, iatom)
              psiGrPsiSum(idir) = psiGrPsiSum(idir) + psiHepsGrPsi(idir, iband, iatom)
              !TODO : Spin missing
              surfIntSum(idir) = surfIntSum(idir) + surfIntMT(iband, idir, iatom) - eig(iband, ikpt, 1) * overlap(iband, idir, iatom)
            end do ! iband
          end do ! idir
        end do ! ieqat
      end do ! itype

    end do ! ikpt

    if (output) then
      call Fclose(1000)
    end if

    ! Write out the sum over k-points and band for <grPsi|H - eps|Psi>, <Psi|H - eps|grPsi> and the MT surface integral with the
    ! integrand e_r Psi^* (H - eps) Psi
    if ( any(abs(grPsiPsiSum(:))> 1e-7) ) then
      passed = passed .and. .false.
      write(logUnit, '(a)') 'Error: k-point and band sum does not cancel for grPsiPsiSum'
      write(logUnit, '(3(2f15.8,2x))') grPsiPsiSum(1), grPsiPsiSum(2), grPsiPsiSum(3)
    end if

    if ( any(abs(psiGrPsiSum(:)) > 1e-7) ) then
      passed = passed .and. .false.
      write(logUnit, '(a)') 'Error: k-point and band sum does not cancel for psiGrPsiSum'
      write(logUnit, '(3(2f15.8,2x))') psiGrPsiSum(1), psiGrPsiSum(2), psiGrPsiSum(3)
    end if

    if ( any(abs(surfIntSum(:)) > 1e-7) ) then
      passed = passed .and. .false.
      write(logUnit, '(a)') 'Error: k-point and band sum does not cancel for surfIntSum'
      write(logUnit, '(3(2f15.8,2x))') surfIntSum(1), surfIntSum(2), surfIntSum(3)
    end if

    ! TODO : sum passed varibales together for all tests
    ! TODO : nicer output in files and logfile
    ! TODO : rename for nicer variables
    ! TODO : Check test for more than one atom
    ! TODO : Nicer formatting
    ! todo : LOs not checked
    if (passed) then
      write (logUnit, '(a)')'                                                                             |__ passed!'
    else
      write (logUnit, '(a)')'                                                                             |__ failed!'
      call juDFT_warn('Test of nabla <Psi_kn|H - eps|Psi_kn> - oint dS Psi_kn^* (H - eps) Psi_kn = 0 failed', &
        & calledby='TestGrMatElemPsiHepsPsiGaussTheo', hint='Check log file for hints where to debug.')
    end if

  end subroutine TestGrMatElemPsiHepsPsiGaussTheo

  !>-------------------------------------------------------------------------------------------------------------------------------------
  !> @author
  !> Christian-Roman Gerhorst
  !>
  !> @brief
  !> Tests <Psi|H - eps|BarPsi>_MT and <BarPsi|H - eps|Psi>_MT methods
  !>
  !> @details
  !> We need the multiplication of <varphi | H - eps | varphi > with scalar and vectorial abcofs to bra or ket in the Sternheimer
  !> equation but also in the Pulay contributions to the dynamical matrix. We have two methods, the old one is used in the old
  !> Fleur code and a new one used in the Pulay contributions to the dynamical matrix. One method is to use the tlmplm routines and
  !> the second variation (Sternheimer), the other is to setup <varphi|H - eps|varphi> and then calculate the small matching
  !> coefficients contrust different types of large matching coefficients and then mulitply to <varphi|H - eps|varphi>.
  !> Due to the fact, that the action of the spherical Hamiltonian is non-hermitian in the current implementation, there is a
  !> difference between <Psi|H - eps|BarPsi> and <BarPsi|H - eps|Psi>.
  !> This routine tets the two methods for calculating such matrix elements by comparing them.
  !>
  !> @param[in] atoms        : Atoms type, see types.f90.
  !> @param[in] dimens       : Dimension type, see types.f90.
  !> @param[in] kpts         : K-points type, see types.f90.
  !> @param[in] input        : Input type, see types.f90.
  !> @param[out] enpara      : Energy parameter type, see types.f90.
  !> @param[out] eig         : Contains Kohn\--Sham eigenvalues.
  !> @param[out] kveclo      : Basis G-vectors of local orbitals.
  !> @param[out] nv          : Number of LAPW G-basis vectors for given k-point.
  !> @param[out] mapGbas     : For various k-points G-basis vectors occur more than once, thus they are only stored once in juPhon. This
  !>                           pointer array contains the right index for GbasVec array to "unfold" G-basis vectors again.
  !> @param[out] gBas        : G-basis vectors
  !> @param[out] ne          : Number of eigenvalues per k-point.
  !> @param[out] El          : Contains LAPW and LOs energy parameters.
  !> @param[in]  logUnit     : Unit number for juPhon.log.
  !>
  !>-------------------------------------------------------------------------------------------------------------------------------------
  subroutine testPsiHepsTildePsi( atoms, cell, lathar, dimens, kpts, results, Veff0, sym, usdus, enpara, input,    &
      & clnu_atom, nmem_atom, mlh_atom, gBas, gBasUnwrap, kveclo, nv, eig,   &
      & nobd, z, mapKpq2K, nRadFun, rbas1, rbas2, El, iloTable, uuilon, duilon, ulouilopn, ilo2p, logUnit, vEff0MTsh )

#include "cppmacro.h"

    use m_types, only : t_atoms, t_cell, t_sphhar, t_stars, t_dimension, t_kpts, t_results, t_potential, t_sym, t_usdus, t_enpara, t_input, t_tlmplm, t_noco
    use m_jpConstants, only : iu
    use m_jpSetupDynMatHelper, only : CalcHnGrV0Varphi
    use m_jpSetupDynMat, only : CalcScalBasfMatElems
    use m_abcof3, only : abcof3
     
    use mod_juPhonUtils, only : Derivative
    use m_jpSternhPulaySurface, only : tlmplm4H0, calcHS0MT
    use m_abcof
    use m_juDFT_NOstopNO, only : juDFT_warn

    implicit none

    ! Type parameters
    type(t_atoms),                  intent(in) :: atoms
    type(t_cell),                   intent(in) :: cell
    type(t_sphhar),                 intent(in) :: lathar
    type(t_dimension),              intent(in) :: dimens
    type(t_kpts),                   intent(in) :: kpts
    type(t_results),                intent(in) :: results
    type(t_potential),              intent(in) :: Veff0
    type(t_sym),                    intent(in) :: sym
    type(t_usdus),                  intent(in) :: usdus
    type(t_enpara),                 intent(in) :: enpara
    type(t_input),                  intent(in) :: input

    ! Scalar parameters
    integer,                        intent(in) :: logUnit

    ! Array parameters
    complex,                        intent(in) :: clnu_atom(:, 0:, :)
    integer,                        intent(in) :: nmem_atom(0:, :)
    integer,                        intent(in) :: mlh_atom(:, 0:, :)
    integer,                        intent(in) :: gBas(:, :)
    integer,                        intent(in) :: gBasUnwrap(:, :, :)
    integer,                        intent(in) :: kveclo(:,:)
    integer,                        intent(in) :: nv(:, :)
    real,                           intent(in) :: eig(:, :, :)
    integer,                        intent(in) :: nobd(:, :)
    integer,                        intent(in) :: mapKpq2K(:, :)
    integer,                        intent(in) :: nRadFun(0:, :)
    real,                           intent(in) :: rbas1(:, :, 0:, :, :)
    real,                           intent(in) :: rbas2(:, :, 0:, :, :)
    real,                           intent(in) :: El(:, 0:, :, :)
    integer,                        intent(in) :: iloTable(:, 0:, :)
    MCOMPLEX,                       intent(in) :: z(:, :, :, :)
    real,                           intent(in) :: uuilon(:, :)
    real,                           intent(in) :: duilon(:, :)
    real,                           intent(in) :: ulouilopn(:, :, :)
    integer,                        intent(in) :: ilo2p(:, :)
    complex,                        intent(in) :: vEff0MTsh(:, :, :, :)

    ! Type variable
    type(od_inp)                               :: odi
    type(od_sym)                               :: ods
    type(t_tlmplm)                             :: tdHS0
    type(t_noco)                               :: noco

    ! Scalar variables
    integer                                    :: itype
    integer                                    :: ilh
    integer                                    :: imesh
    integer                                    :: idir
    integer                                    :: iatom
    integer                                    :: ieqat
    integer                                    :: oqn_l
    integer                                    :: lm_pre
    integer                                    :: mqn_m
    integer                                    :: iband
    integer                                    :: nmat
    integer                                    :: ikpt
    integer                                    :: iBas
    integer                                    :: idirC
    integer                                    :: idirR
    integer                                    :: ikpq
    integer                                    :: iDatomA
    integer                                    :: iDtypeA
    integer                                    :: iDeqatA
    integer                                    :: lm
    integer                                    :: lmp
    integer                                    :: pMaxLocal
    integer                                    :: iDatomB
    integer                                    :: iDtypeB
    integer                                    :: iDeqatB
    integer                                    :: lmpMax
    integer                                    :: nRadFunMax
    integer                                    :: iradf
    integer                                    :: ptsym
    integer                                    :: imem
    integer                                    :: maxlmp ! todo should be the same as lmpT(:)
    integer                                    :: p
    integer                                    :: ilo
    logical                                    :: passed = .true.


!    ! rho1IR     : interstitial first variation of density
!    ! rho1MT     : muffin-tin first variation of density
!    ! w_vEff1IR  : warped first variation of effective potential in interstitial
!    ! vEff1MT    : first variation of effective potential in muffin-tin
!    ! Array variables
    complex,           allocatable             :: z1nG(:, :, :, :)
    integer,           allocatable             :: lmpT(:)
    real,              allocatable             :: gbasExtKpq(:, :)
    complex,           allocatable             :: a(:, :, :)
    complex,           allocatable             :: b(:, :, :)
    complex,           allocatable             :: bascof_lo(:, :, :, :, :)
    integer,           allocatable             :: ngoprI(:)
    complex,           allocatable             :: ab0cofScl(:)
    complex,           allocatable             :: ab1cofVec(:, :, :)
    complex,           allocatable             :: abcofSumMat(:, :, :, :)
    complex,           allocatable             :: abcofMat(:)
    real,              allocatable             :: varphiVarphi(:, :, :)
    complex,           allocatable             :: varphiHvarphi(:, :, :)
    real,              allocatable             :: r2(:)
    real,              allocatable             :: varphi1(:, :, :)
    real,              allocatable             :: varphi2(:, :, :)
    complex,           allocatable             :: vEff0MtSpH(:, :)
    complex,           allocatable             :: hVarphi(:, :, :, :)
    integer,           allocatable             :: grVarphiChLout(:, :)
    integer,           allocatable             :: grVarphiChMout(:, :)
    real,              allocatable             :: grVarphiCh1(:, :, :, :)
    real,              allocatable             :: grVarphiCh2(:, :, :, :)
    complex,           allocatable             :: vEff0NsphGrVarphi(:, :, :, :, :)
    complex,           allocatable             :: r2grVeff0SphVarphi(:, :, :, :, :)
    complex,           allocatable             :: r2grVeff0SphVarphiDummy(:, :, :, :, :)
    complex,           allocatable             :: aKpq(:, :, :)
    complex,           allocatable             :: bKpq(:, :, :)
    complex,           allocatable             :: bascof_loKpq(:, :, :, :, :)
    complex,           allocatable             :: z1Gext(:)
    complex,           allocatable             :: varphiPsi(:)
    complex,           allocatable             :: varphiHPsi(:)
    complex,           allocatable             :: varphiPsiCjg(:)
    complex,           allocatable             :: varphiHPsiCjg(:)
    complex,           allocatable             :: acofKet(:, :, :)
    complex,           allocatable             :: bcofKet(:, :, :)
    complex,           allocatable             :: ccofKet(:, :, :, :)
    complex,           allocatable             :: mCoefK(:, :, :)
    complex,           allocatable             :: acofBar(:, :, :)
    complex,           allocatable             :: bcofBar(:, :, :)
    complex,           allocatable             :: ccofBar(:, :, :, :)
    complex,           allocatable             :: mCoefKv(:, :)
    MCOMPLEX,          allocatable             :: h0MTBv(:, :)
    MCOMPLEX,          allocatable             :: s0MTBv(:, :)
    MCOMPLEX,          allocatable             :: h0MTKv(:, :)
    MCOMPLEX,          allocatable             :: s0MTKv(:, :)
    complex,           allocatable             :: zBar(:, :)
    integer,           allocatable             :: nlo_atom(:)
    complex                                    :: psiHepsPsi(3, 3)
    complex                                    :: psiHepsPsiCjg(3, 3)
    real                                       :: gExt(3)
    real                                       :: kext(3)
    real                                       :: kpqExt(3)
    complex                                    :: dynMatPuMt1(3, 3)
    complex                                    :: dynMatPuMt2(3, 3)

    write(logUnit, '(a)') 'Testing <BarPsi_kn|H - eps|Psi_kn> and <Psi_kn|H - eps|BarPsi_kn>'
    write(logUnit, '(a)') '-----------------------------------------------------------------'


    ! Calculate contributions of the 1st and the 2nd braket for the MT with a z1 and i (k + G), z1 is set to the q=0 solution again

    ! ==============================================================================================================================
    ! METHOD 1: Calculate the overlap of the basis functions and the action of the Hamiltonian onto the basis function and project
    !           it onto other basis function. Then, construct the wished kinds of abcof that can be multiplied to latter arrays.
    ! ==============================================================================================================================
    ! Quantities for initialization
    allocate( lmpT(atoms%ntype) )
    lmpT(:) = 0
    do itype = 1, atoms%ntype
      lmpT(itype) = sum( [ ( (2 * oqn_l + 1)* nRadFun(oqn_l, itype), oqn_l = 0, atoms%lmax(itype) ) ] )
    end do ! itype
    lmpMax     = maxval( lmpT(:) )
    nRadFunMax = maxval( nRadFun(:, :) )

    allocate( varphiVarphi(lmpMax, lmpMax, atoms%ntype), varphiHvarphi(lmpMax, lmpMax, atoms%nat))
    allocate( r2(atoms%jmtd) )
    allocate( varphi1(atoms%jmtd, nRadFunMax, 0:atoms%lmaxd), varphi2(atoms%jmtd, nRadFunMax, 0:atoms%lmaxd) )
    allocate( vEff0MtSpH( atoms%jmtd, 0:dimens%lmd), vEff0NsphGrVarphi(2, atoms%jmtd, 0:(atoms%lmaxd + 3)**2 - 1, lmpMax, -1:1), &
            & r2grVeff0SphVarphi(2, atoms%jmtd, 0:(atoms%lmaxd + 3)**2 - 1, lmpMax, 3) )
    allocate( r2grVeff0SphVarphiDummy(2, atoms%jmtd, 0:(atoms%lmaxd + 3)**2 - 1, lmpMax, 3) )
    allocate( hVarphi(2, atoms%jmtd, 0:(atoms%lmaxd + 3)**2 - 1, lmpMax)  )
    allocate( grVarphiCh1(atoms%jmtd, 2, lmpMax, -1:1), grVarphiCh2(atoms%jmtd, 2, lmpMax, -1:1), &
            & grVarphiChLout(2, 0:atoms%lmaxd), grVarphiChMout(-atoms%lmaxd:atoms%lmaxd, -1:1) )

    varphiVarphi(:, :, :) = cmplx(0., 0.)
    varphiHvarphi(:, :, :) = cmplx(0., 0.)
    vEff0NsphGrVarphi(:, :, :, :, :) = cmplx(0., 0.)
    r2grVeff0SphVarphi(:, :, :, :, :) = cmplx(0., 0.)
    r2grVeff0SphVarphiDummy(:, :, :, :, :) = cmplx(0., 0.)
    hVarphi(:, :, :, :) = cmplx(0., 0.)
    grVarphiCh1(:, :, :, :) = 0.
    grVarphiCh2(:, :, :, :) = 0.
    grVarphiChLout(:, :) = 0
    grVarphiChMout(:, :) = 0

    iatom = 0
    do itype = 1, atoms%ntype

      ! Precalculate radial Jacobi determinant for later integrals
      r2(:) = 0.
      do imesh = 1, atoms%jri(itype)
        r2(imesh) = atoms%rmsh(imesh, itype) * atoms%rmsh(imesh, itype)
      end do ! imesh

      varphi1(:, :, :) = 0.
      varphi2(:, :, :) = 0.
      do oqn_l = 0, atoms%lmax(itype)
        do iradf = 1, nRadFun(oqn_l, itype)
          do imesh = 1, atoms%jri(itype)
            ! In Fleur the radial solutions u_lp are multiplied by a factor r by default to avoid an additional multiplication of the
            ! Jacobi determinant r^2 in radial integrals given spherical coordinates.
            varphi1(imesh, iradf, oqn_l) = rbas1(imesh, iradf, oqn_l, itype, 1) / atoms%rmsh(imesh, itype)
            varphi2(imesh, iradf, oqn_l) = rbas2(imesh, iradf, oqn_l, itype, 1) / atoms%rmsh(imesh, itype)
          end do ! imesh
        end do ! iradf
      end do ! oqn_l
      ! Calculate H |varphi>

      do ieqat = 1, atoms%neq(itype)
        iatom = iatom + 1
        !todo block begin                             place this block into calchngrv0varphi
        ! Expand the coefficients of the lattice-harmonic given potential into spherical-harmonic coefficients for the given atom.
        vEff0MtSpH(:, :) = cmplx(0.0, 0.0)
        ptsym = atoms%ntypsy(iatom)
        do ilh = 0, lathar%nlh(ptsym)
          oqn_l = lathar%llh(ilh, ptsym)
          lm_pre = oqn_l * (oqn_l + 1)
          do imem = 1, nmem_atom(ilh, iatom)
            mqn_m = mlh_atom(imem, ilh, iatom)
            lm = lm_pre + mqn_m
            !todo one could only evaluate the vEff0MtSpH which have a contribution, i.e. the oqn_l and mqn_m which are in llh and mlh_atom
            ! maybe construct a pointer and run only over them to make it memory efficient.
            do imesh = 1, atoms%jri(itype)
              vEff0MtSpH(imesh, lm) = vEff0MtSpH(imesh, lm) + Veff0%vr(imesh, ilh, itype, 1) * clnu_atom(imem, ilh, iatom)
            end do ! imesh
          end do ! imem
        end do ! ilh

        !todo block end.............................................

        ! Calculate the action of the Hamiltonian onto the basis function, of the non-spherical potential onto the gradient of
        ! the basis function and of the gradient of the spherical effective potential onto the basis function.
        hVarphi = cmplx(0.0, 0.0)
        call CalcHnGrV0Varphi( atoms, lathar, itype, iatom, lmpMax, El, varphi1, varphi2, nRadFun, vEff0MtSpH, Veff0%vr(:, :, :, 1), clnu_atom, &
          & nmem_atom, mlh_atom, grVarphiCh1, grVarphiCh2, grVarphiChLout, grVarphiChMout, hVarphi, vEff0NsphGrVarphi, r2grVeff0SphVarphi, r2grVeff0SphVarphiDummy )

        ! Calculate scalar basis function matrix elements
        call CalcScalBasfMatElems( atoms, itype, iatom, nRadFun, r2, varphi1, varphi2, hVarphi, varphiVarphi, varphiHvarphi )
      end do ! ieqat
    end do ! itype

    deallocate( r2, varphi1, varphi2, vEff0MtSpH, vEff0NsphGrVarphi, r2grVeff0SphVarphi, hVarphi, grVarphiCh1, grVarphiCh2, r2grVeff0SphVarphiDummy )

    ! Allocate the k-dependent and band-dependent arrays
    allocate(ngoprI(atoms%nat))
    ngoprI(:) = 1
    allocate( gbasExtKpq( dimens%nbasfcn, 3) )
    allocate( a( dimens%nvd, 0:dimens%lmd, atoms%nat), b(dimens%nvd, 0:dimens%lmd, atoms%nat), &
      & bascof_lo(3, -atoms%llod:atoms%llod, 4 * atoms%llod + 2, atoms%nlod, atoms%nat) )
    allocate( aKpq( dimens%nvd, 0:dimens%lmd, atoms%nat), bKpq(dimens%nvd, 0:dimens%lmd, atoms%nat), &
      & bascof_loKpq(3, -atoms%llod:atoms%llod, 4 * atoms%llod + 2, atoms%nlod, atoms%nat) )
    allocate( ab0cofScl(lmpMax), ab1cofVec( lmpMax, atoms%nat, 3), abcofMat(nRadFunMax), &
      & abcofSumMat(lmpMax, atoms%nat, 3, 3))
    allocate( z1Gext(dimens%nvd) )
    allocate( varphiPsi(lmpMax), varphiHPsi(lmpMax) )
    allocate( varphiPsiCjg(lmpMax), varphiHPsiCjg(lmpMax) )
    allocate( z1nG(dimens%nbasfcn, 3, atoms%nat, maxval(nobd(:, :))) )

    dynMatPuMt1(:, :) = cmplx(0., 0.)

    do ikpt = 1, kpts%nkpt
      !todo This is only now to test q = 0, but in principle we can also test another qs.
      ikpq = ikpt

      ! Setup the z1 for q = 0
      z1nG(:, :, :, :) = cmplx(0., 0.)
      kExt(1:3) = matmul(cell%bmat, kpts%bk(1:3, ikpt))
      gExt(:) = 0.
      do iBas = 1, nv(1, ikpt) !+ atoms%nlotot
        Gext(1:3) = matmul( cell%bmat(1:3, 1:3), gBas(1:3, gBasUnwrap(iBas, ikpt, 1)))
        do iband = 1, nobd(ikpt, 1)
          do idir = 1, 3
            z1nG(iBas, idir, 1, iband) = -iu * ( kExt(idir) + Gext(idir) ) * z(iBas, iband, ikpt, 1)
          end do ! iBas
        end do ! idir
      end do ! iband

      ! We leave the k + q infrastructure, so that it is as close as possible to the actual implementation in the production code.
      ! Furthermore, we are able to extend the functionality later to other qs different from q = 0.
      ! We have terms with factors of k, k + q or G, where a transformation into cartesian coordinates is required.
      kpqExt(1:3) = matmul(cell%bmat, kpts%bk(1:3, ikpq))

      gbasExtKpq(:, :) = 0.
      ! The nv has to be taken at ikpq. We assign the G-basis vectors once so later we can just loop without jumps over iBas.
      do iBas = 1, nv(1, ikpq)
        gbasExtKpq(iBas, 1:3) = matmul( cell%bmat, gbas(1:3, gBasUnwrap(iBas, ikpq, 1)) )
      end do ! iBas
      do iBas = nv(1, ikpq) + 1, nv(1, ikpq) + atoms%nlotot
        gbasExtKpq(iBas, 1:3) = matmul( cell%bmat, gbas(1:3, gBasUnwrap(kveclo(iBas, ikpq), ikpq, 1)) )
      end do ! iBas

      ! Calculate the basis matching coefficients at k + q = k' and the matching coefficients at k.
      nmat = nv(1, ikpq) + atoms%nlotot
      aKpq(:, :, :) = cmplx(0.0, 0.0)
      bKpq(:, :, :) = cmplx(0.0, 0.0)
      bascof_loKpq(:, :, :, :, :) = cmplx(0.0, 0.0)
      call abcof3( atoms%lmaxd, atoms%ntype, atoms%nat, sym%nop, dimens%nvd, dimens%jspd, 1, dimens%lmd, dimens%nbasfcn, &
        & atoms%llod, atoms%nlod, atoms%nlotot, sym%invtab, atoms%ntype, sym%mrot, ngoprI, atoms%taual, atoms%neq, atoms%lmax, &
        & atoms%rmt, cell%omtil, cell%bmat, cell%bbmat, kpts%bk(:, ikpq), gbas(1, gBasUnwrap(:nv(1, ikpq), ikpq, 1)), &
        & gbas(2, gBasUnwrap(:nv(1, ikpq), ikpq, 1)), gbas(3, gBasUnwrap(:nv(1, ikpq), ikpq, 1)), nv(:, ikpq), nmat, &
        & usdus%us, usdus%dus, usdus%uds, usdus%duds, usdus%ddn, atoms%invsat, sym%invsatnr, usdus%ulos, usdus%uulon, usdus%dulon, &
        & usdus%dulos, atoms%llo, atoms%nlo, atoms%l_dulo, atoms%lapw_l, kveclo(:,ikpq), odi, ods, aKpq, bKpq, bascof_loKpq )

      nmat = nv(1, ikpt) + atoms%nlotot
      a(:, :, :) = cmplx(0.0, 0.0)
      b(:, :, :) = cmplx(0.0, 0.0)
      bascof_lo(:, :, :, :, :) = cmplx(0.0, 0.0)
      call abcof3( atoms%lmaxd, atoms%ntype, atoms%nat, sym%nop, dimens%nvd, dimens%jspd, 1, dimens%lmd, dimens%nbasfcn, &
        & atoms%llod, atoms%nlod, atoms%nlotot, sym%invtab, atoms%ntype, sym%mrot, ngoprI, atoms%taual, atoms%neq, atoms%lmax, &
        & atoms%rmt, cell%omtil, cell%bmat, cell%bbmat, kpts%bk(:, ikpt), gbas(1, gBasUnwrap(:nv(1, ikpt), ikpt, 1)), &
        & gbas(2, gBasUnwrap(:nv(1, ikpt), ikpt, 1)), gbas(3, gBasUnwrap(:nv(1, ikpt), ikpt, 1)), nv(:, ikpt), nmat, &
        & usdus%us, usdus%dus, usdus%uds, usdus%duds, usdus%ddn, atoms%invsat, sym%invsatnr, usdus%ulos, usdus%uulon, usdus%dulon, &
        & usdus%dulos, atoms%llo, atoms%nlo, atoms%l_dulo, atoms%lapw_l, kveclo(:,ikpt), odi, ods, a, b, bascof_lo )

      do iband = 1, nobd(ikpt, 1)
        iDatomA = 0
        do iDtypeA = 1, atoms%ntype
          do iDeqatA = 1, atoms%neq(iDtypeA)
            iDatomA = iDatomA + 1

            ! Conventional matching coefficients as we know them from FLEUR.
            ab0cofScl(:) = cmplx(0., 0.)
            lmp = 0
            lm  = 0
            do oqn_l = 0, atoms%lmax(iDtypeA)
              do mqn_m = - oqn_l, oqn_l
                pMaxLocal = nRadFun(oqn_l, iDtypeA)
                ! p = 1
                ab0cofScl(lmp + 1) = iu**oqn_l * dot_product( conjg(z(:nv(1, ikpt), iband, ikpt, 1)), a(:nv(1, ikpt), lm, iDatomA) )
                ! p = 2
                ab0cofScl(lmp + 2) = iu**oqn_l * dot_product( conjg(z(:nv(1, ikpt), iband, ikpt, 1)), b(:nv(1, ikpt), lm, iDatomA) )
                ! Add LO contributions
                do iradf = 3, pMaxLocal
                  ! p = 1
                  ab0cofScl(lmp + 1) = ab0cofScl(lmp + 1) + iu**oqn_l * &
                    & dot_product( conjg(z(nv(1, ikpt) + 1:nv(1, ikpt) + atoms%nlotot, iband, ikpt, 1)), &
                                                         & bascof_lo(1, mqn_m, :atoms%nlotot, iloTable(iradf, oqn_l, iDtypeA), iDatomA) )
                  ! p = 2
                  ab0cofScl(lmp + 2) = ab0cofScl(lmp + 2) + iu**oqn_l * &
                    & dot_product( conjg(z(nv(1, ikpt) + 1:nv(1, ikpt) + atoms%nlotot, iband, ikpt, 1)), &
                                                         & bascof_lo(2, mqn_m, :atoms%nlotot, iloTable(iradf, oqn_l, iDtypeA), iDatomA) )
                  ! 2 < p < LOs for that l and that atom type
                  ab0cofScl(lmp + iradf) = iu**oqn_l * dot_product( conjg(z(nv(1, ikpt) + 1:nv(1, ikpt) + atoms%nlotot, iband, ikpt, 1)), &
                                                         & bascof_lo(3, mqn_m, :atoms%nlotot, iloTable(iradf, oqn_l, iDtypeA), iDatomA) )
                end do ! iradf


                ! This is a precalculation of the 1st and 3rd line of A.50. Actually the resulting lmp should have been primed or at
                ! least the p as lm1Pr = lm. But for sake of performance we place it here.
                ! sake of performance
                lm = lm + 1
                lmp = lmp + pMaxLocal
              end do ! mqn_m
            end do !oqn_l


            ! Calculate the vector-like large matching coefficients using the first-order wave-function expansion coefficients gained
            ! from solving the Sternheimer equation.
            ab1cofVec(:, :, :) = cmplx(0.0, 0.0)
            do idirR = 1, 3
              iDatomB = 0
              do iDtypeB = 1, atoms%ntype
                do iDeqatB = 1, atoms%neq(iDtypeB)
                  iDatomB = iDatomB + 1
                  lm  = 0
                  lmp = 0
                  do oqn_l = 0, atoms%lmax(iDtypeB)
                    do mqn_m = -oqn_l, oqn_l
                      pMaxLocal = nRadFun(oqn_l, iDtypeB)
                      ! p = 1
                      ab1cofVec(lmp + 1, iDatomB, idirR) = iu**oqn_l * &
                                     & dot_product( conjg( z1nG(:nv(1, ikpq), idirR, iDatomA, iband) ), aKpq(:nv(1, ikpq), lm, iDatomB) )
                      ! p = 2
                      ab1cofVec(lmp + 2, iDatomB, idirR) = iu**oqn_l * &
                                     & dot_product( conjg( z1nG(:nv(1, ikpq), idirR, iDatomA, iband) ), bKpq(:nv(1, ikpq), lm, iDatomB) )
                      do iradf = 3, pMaxLocal
                        ! p = 1
                        ab1cofVec(lmp + 1, iDatomB, idirR) = ab1cofVec(lmp + 1, iDatomB, idirR) &
                          & + iu**oqn_l * dot_product( conjg( z1nG(nv(1, ikpq) + 1:nv(1, ikpq) + atoms%nlotot, idirR, iDatomA, iband) ), &
                                                      & bascof_loKpq(1, mqn_m, :atoms%nlotot, iloTable(iradf, oqn_l, iDtypeB), iDatomB) )
                        ! p = 2
                        ab1cofVec(lmp + 2, iDatomB, idirR) = ab1cofVec(lmp + 2, iDatomB, idirR) &
                          & + iu**oqn_l * dot_product(conjg(z1nG(nv(1, ikpq) + 1:nv(1, ikpq) + atoms%nlotot, idirR, iDatomA, iband)), &
                                                      & bascof_loKpq(2, mqn_m, :atoms%nlotot, iloTable(iradf, oqn_l, iDtypeB), iDatomB) )
                        ! 2 < p < LOs for that l and that atom type
                        ab1cofVec(lmp + iradf, iDatomB, idirR) = &
                          & iu**oqn_l * dot_product( conjg( z1nG(nv(1, ikpq) + 1:nv(1, ikpq) + atoms%nlotot, idirR, iDatomA, iband) ),&
                                                      & bascof_loKpq(3, mqn_m, :atoms%nlotot, iloTable(iradf, oqn_l, iDtypeB), iDatomB) )
                      end do ! iradf
                      lm = lm + 1
                      lmp = lmp + nRadFun(oqn_l, iDtypeB)
                    end do ! mqn_m
                  end do ! oqn_l
                end do ! iDeqatB
              end do ! iDtypeB
            end do ! idirR

            ! Matrix abcof
            abcofSumMat(:, :, :, :) = cmplx(0., 0.)
            do idirC = 1, 3
              do idirR = 1, 3
                lm  = 0
                lmp = 0
                ! Precalculation of A.49. PhD thesis A. Klueppelberg
                z1Gext(:) = cmplx(0., 0.)
                do iBas = 1, nv(1, ikpq) + atoms%nlotot
                  z1Gext(iBas) = z1nG(iBas, idirR, iDatomA, iband) * gbasExtKpq(iBas, idirC)
                end do ! iBas

                iDatomB = 0
                do iDtypeB = 1, atoms%ntype
                  do iDeqatB = 1, atoms%neq(iDtypeB)
                    iDatomB = iDatomB + 1
                    lm  = 0
                    lmp = 0
                    do oqn_l = 0, atoms%lmax(iDtypeB)
                      do mqn_m = -oqn_l, oqn_l
                        abcofMat(:) = cmplx(0.0, 0.0)
                        pMaxLocal = nRadFun(oqn_l, iDtypeB)

                        ! p = 1
                        abcofMat(1) = iu**oqn_l * dot_product( conjg(z1Gext(:nv(1, ikpq)) ), aKpq(:nv(1, ikpq), lm, iDatomB) )
                        ! p = 2
                        abcofMat(2) = iu**oqn_l * dot_product( conjg(z1Gext(:nv(1, ikpq)) ), bKpq(:nv(1, ikpq), lm, iDatomB) )
                        do iradf = 3, pMaxLocal
                          ! p = 1
                          abcofMat(1) = abcofMat(1) &
                            & + iu**oqn_l * dot_product( conjg( z1Gext(nv(1, ikpq) + 1:nv(1, ikpq) + atoms%nlotot)), &
                                                      & bascof_loKpq(1, mqn_m, :atoms%nlotot, iloTable(iradf, oqn_l, iDtypeB), iDatomB) )
                          ! p = 2
                          abcofMat(2) = abcofMat(2) &
                            & + iu**oqn_l * dot_product( conjg( z1Gext(nv(1, ikpq) + 1:nv(1, ikpq) + atoms%nlotot)), &
                                                     & bascof_loKpq(2, mqn_m, :atoms%nlotot, iloTable(iradf, oqn_l, iDtypeB), iDatomB) )
                          ! 2 < p < LOs for that l and that atom type
                          abcofMat(iradf) = iu**oqn_l * dot_product( conjg(z1Gext(nv(1, ikpq) + 1:nv(1, ikpq) + atoms%nlotot) ) &
                                                    &, bascof_loKpq(3, mqn_m, :atoms%nlotot, iloTable(iradf, oqn_l, iDtypeB), iDatomB) )
                        end do ! p

                        do iradf = 1, nRadFun(oqn_l, iDtypeB)
                          ! Precalculation of the 1st line in A.50 PhD thesis A. Klueppelberg. Here a plus sign stand before the 2
                          ! because abcofSumMat is complex conjugated later!
                          ! todo review it according to correct atom assignment B should be A here
                          abcofSumMat(lmp + iradf, iDatomB, idirR, idirC) = &
                            & iu * ( ab1cofVec(lmp + iradf, iDatomB, idirR) * kpqExt(idirC) + abcofMat(iradf) )
                        end do ! iradf
                        lm = lm + 1
                        lmp = lmp + pMaxLocal
                      end do ! mqn_m
                    end do ! oqn_l
                  end do ! iDeqatB
                end do ! iDtypeB
              end do ! idirR
            end do ! idirC

            ! Calculation of the final matrix elements.
            iDatomB = 0
            do iDtypeB = 1, atoms%ntype
              do iDeqatB = 1, atoms%neq(iDtypeB)
                iDatomB = iDatomB + 1
                ! Precalculation of the 1st and 3rd line in A.50 PhD thesis A. Klueppelberg.
                varphiPsi(:) = cmplx(0., 0.)
                varphiPsi(1:lmpT(iDtypeA)) = matmul( varphiVarphi(1:lmpT(iDtypeA), 1:lmpT(iDtypeA), iDtypeA), ab0cofScl(1:lmpT(iDtypeA)) )
                varphiHPsi(:) = cmplx(0.0, 0.0)
                varphiHPsi(1:lmpT(iDtypeA)) = matmul( varphiHvarphi(1:lmpT(iDtypeA), 1:lmpT(iDtypeA), iDatomA), ab0cofScl(1:lmpT(iDtypeA)) )
                psiHepsPsi(:, :) = cmplx(0., 0.)
                psiHepsPsiCjg(:, :) = cmplx(0., 0.)
                do idirC = 1, 3
                  do idirR = 1, 3
                    varphiPsiCjg(:) = cmplx(0., 0.)
                    varphiPsiCjg(1:lmpT(iDtypeA)) = matmul( varphiVarphi(1:lmpT(iDtypeA), 1:lmpT(iDtypeA), iDtypeA), abcofSumMat(:lmpT(iDtypeB), iDatomB, idirR, idirC) )
                    varphiHPsiCjg(:) = cmplx(0.0, 0.0)
                    varphiHPsiCjg(1:lmpT(iDtypeA)) = matmul( varphiHvarphi(1:lmpT(iDtypeA), 1:lmpT(iDtypeA), iDatomA), abcofSumMat(:lmpT(iDtypeB), iDatomB, idirR, idirC) )
                    ! Calculation of the 2nd and 4th line of A.50 PhD thesis A. Klueppelberg. The distinction between atom alpha and
                    ! beta in abcofSumMat has already be performed upwards.
                    psiHepsPsiCjg(idirR, idirC) = dot_product( ab0cofScl(1:lmpT(iDtypeA)), varphiHPsiCjg(:lmpT(iDtypeB)) ) &
                                             & - eig(iband, ikpt, 1) * dot_product( ab0cofScl(1:lmpT(iDtypeA)), varphiPsiCjg(:lmpT(iDtypeB)) )
                    psiHepsPsi(idirR, idirC) = dot_product(abcofSumMat(:lmpT(iDtypeB), iDatomB, idirR, idirC), varphiHPsi(:lmpT(iDtypeB)) ) &
                                             & - eig(iband, ikpt, 1) * dot_product(abcofSumMat(:lmpT(iDtypeB), iDatomB, idirR, idirC), varphiPsi(:lmpT(iDtypeB)) )
                    !todo this does not work for more than one atom
                    dynMatPuMt1(idirR, idirC) = dynMatPuMt1(idirR, idirC)  &
                      & + 2 * results%w_iks(iband, ikpt, 1) * (psiHepsPsi(idirR,idirC) + psiHepsPsiCjg(idirR,idirC))
                  end do ! idirR
                end do ! idirC
              end do ! iDeqatB
            end do ! iDtypeB
          end do ! iDeqatA
        end do ! iDtypeA
      end do ! iband
    end do ! ikpt

    if (.false.) then
      write(*, '(a)') 'BarPsi (H - eps) Psi + Psi (H - eps) BarPsi 2nd method'
      write(*, '(3(2(es16.8,1x),3x))') dynMatPuMt1(1, :)
      write(*, '(3(2(es16.8,1x),3x))') dynMatPuMt1(2, :)
      write(*, '(3(2(es16.8,1x),3x))') dynMatPuMt1(3, :)
    end if

    deallocate(varphiVarphi, varphiHvarphi, a, b, bascof_lo, aKpq, bKpq, bascof_loKpq, ab0cofScl, ab1cofVec, abcofMat, abcofSumMat, varphiPsi, varphiHPsi, varphiPsiCjg, varphiHPsiCjg, z1Gext )
    ! ==============================================================================================================================
    ! METHOD 2: Calculate the integrals of the non-spherical potential and the basis functions in the tlmplm routine. Then use the
    !           abcof routine with different z to generate the matching coefficients which are used in the second-variation routine.
    ! ==============================================================================================================================

    ! Calculate the integrals of the non-spherical potential and the basis functions
    call tlmplm4H0( atoms, dimens, enpara, usdus, input, tdHS0, 1, logUnit, rbas1, rbas2, uuilon, duilon, ulouilopn, ilo2p, &
                                                                                                           & vEff0MTsh(:, :, :, 1) )

    maxlmp = maxval( (/ (sum( (/ ((2*oqn_l+1)* nRadFun(oqn_l,itype), oqn_l = 0,atoms%lmax(itype)) /) ),itype=1,atoms%ntype) /) )
    allocate(noco%alph(atoms%ntype), noco%beta(atoms%ntype)) ! Dummy variables to run abcof
    allocate(acofKet(maxval(nobd), 0:dimens%lmd, atoms%nat), bcofKet(maxval(nobd), 0:dimens%lmd, atoms%nat), &
      &ccofKet(-atoms%llod:atoms%llod, maxval(nobd), atoms%nlod, atoms%nat))
    allocate(mCoefK(maxval(nobd), maxlmp, atoms%nat))
    !todo this could be an error when passing the variable with maxval nobd for non metals
    allocate(acofBar(maxval(nobd), 0:dimens%lmd, atoms%nat), bcofBar(maxval(nobd), 0:dimens%lmd, atoms%nat), &
     &ccofBar(-atoms%llod:atoms%llod, maxval(nobd), atoms%nlod, atoms%nat))
    allocate( mCoefKv(maxval(nobd), maxlmp) )
    allocate( h0MTBv(maxval(nobd), maxval(nobd)), s0MTBv(maxval(nobd), maxval(nobd)) )
    allocate( h0MTKv(maxval(nobd), maxval(nobd)), s0MTKv(maxval(nobd), maxval(nobd)) )
    allocate(zBar(dimens%nbasfcn, maxval(nobd)))
    allocate(nlo_atom(atoms%nat))
    nlo_atom = 0

    dynMatPuMt2(:, :) = cmplx(0., 0.)
    do ikpt = 1, kpts%nkpt

      ! Conventional matching coefficients.
      nmat = nv(1, ikpt) + atoms%nlotot
      acofKet(:, :, :) = cmplx(0.0, 0.0)
      bcofKet(:, :, :) = cmplx(0.0, 0.0)
      ccofKet(:, :, :, :) = cmplx(0.0, 0.0)
      call abcof ( atoms%lmaxd, atoms%ntype, dimens%neigd, maxval(nobd), atoms%nat, sym%nop, dimens%nvd, dimens%jspd, dimens%lmd, &
        & dimens%nbasfcn, atoms%llod, atoms%nlod, atoms%nlotot, sym%invtab, atoms%ntype, sym%mrot, ngoprI, atoms%taual, atoms%neq, &
        & atoms%lmax, atoms%rmt, cell%omtil, cell%bmat, cell%bbmat, kpts%bk(:, ikpt), &
        & gbas(1, gbasUnwrap(:nv(input%jspins, ikpt), ikpt, input%jspins)), gbas(2, gbasUnwrap(:nv(input%jspins, ikpt), ikpt, input%jspins)), &
        & gbas(3, gbasUnwrap(:nv(input%jspins, ikpt), ikpt, input%jspins)), nv(:, ikpt),  nmat, nobd(ikpt, 1), z(:, :, ikpt, 1), usdus%us(:, :, 1), &
        & usdus%dus(:, :, 1), usdus%uds, usdus%duds(:, :, 1), usdus%ddn(:, :, 1), atoms%invsat, sym%invsatnr, usdus%ulos(:, :, 1), &
        & usdus%uulon(:, :, 1), usdus%dulon(:, :, 1),  usdus%dulos(:, :, 1), atoms%llo, atoms%nlo, atoms%l_dulo, &
        & atoms%lapw_l, noco%l_noco, noco%l_ss, 1, noco%alph, noco%beta, noco%qss, kveclo(:, ikpt), odi, ods, &
        & acofKet, bcofKet, ccofKet)

      ! Rearrange potential acofs, bcofs, ccofs ensuring an consistent handling of LO terms within the loop the structure compared to LAPWs
      iatom  = 0
      mCoefK =  0
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
              mCoefK(:nobd(ikpt, 1), lmp, iatom) = acofKet(:nobd(ikpt, 1), lm, iatom)
              !p = 2
              lmp = lmp + 1
              mCoefK(:nobd(ikpt, 1), lmp, iatom) = bcofKet(:nobd(ikpt, 1), lm, iatom)
              !LOs
              do p = 3, nRadFun(oqn_l, itype)
                ilo = iloTable(p, oqn_l, itype)
                lmp = lmp + 1
                mCoefK(:nobd(ikpt, 1), lmp, iatom) = ccofKet(mqn_m, :nobd(ikpt, 1), ilo, iatom)
              end do
            end do
          end do
        end do
      end do

      do idirC = 1, 3
        do idirR = 1, 3
          ! Calculating zPrime for the bar acof bcof and ccof of the MT Hamiltonian matrix element
          zBar = 0
          ! is that correct that it oes to ne(ikpt and not ne(ikpq)
          do iband = 1, nobd(ikpt, 1)
            do iBas = 1, nv(1, ikpt)
             ! todo To increase performance, we can source out the matrix multiplication of Gext as we calculate a lot redundantly.
              Gext(1:3) = matmul( cell%bmat(1:3, 1:3), gbas(1:3, gBasUnwrap(iBas, ikpt, 1)) + kpts%bk(1:3, ikpt))
              zBar(iBas, iband) = Gext(idirR) * Gext(idirC) * z(iBas, iband, ikpt, 1)
            end do
          end do

          ! Calculating acofBar, bcofBar, ccofBar using zBar
          ! When rewriting abcof dimensions of acof, bcof, and ccof and be made even smaller
          acofBar(:, :, :) = cmplx(0., 0.)
          bcofBar(:, :, :) = cmplx(0., 0.)
          ccofBar(:, :, :, :) = cmplx(0., 0.)
          nmat = nv(1, ikpt) + atoms%nlotot
          call abcof ( atoms%lmaxd, atoms%ntype, maxval(nobd), maxval(nobd), atoms%nat, sym%nop, dimens%nvd, dimens%jspd, dimens%lmd, &
            & dimens%nbasfcn, atoms%llod, atoms%nlod, atoms%nlotot, sym%invtab, atoms%ntype, sym%mrot, ngoprI, atoms%taual, atoms%neq, &
            & atoms%lmax, atoms%rmt, cell%omtil, cell%bmat, cell%bbmat, kpts%bk(:, ikpt), &
            & gbas(1, gBasUnwrap(:nv(input%jspins, ikpt), ikpt, input%jspins)), gbas(2, gBasUnwrap(:nv(input%jspins, ikpt), ikpt, input%jspins)), &
            & gbas(3, gBasUnwrap(:nv(input%jspins, ikpt), ikpt, input%jspins)), nv(:, ikpt),  nmat, nobd(ikpt, 1), zBar(:, :), usdus%us(:, :, 1), &
            & usdus%dus(:, :, 1), usdus%uds, usdus%duds(:, :, 1), usdus%ddn(:, :, 1), atoms%invsat, sym%invsatnr, usdus%ulos(:, :, 1), &
            & usdus%uulon(:, :, 1), usdus%dulon(:, :, 1),  usdus%dulos(:, :, 1), atoms%llo, atoms%nlo, atoms%l_dulo, &
            & atoms%lapw_l, noco%l_noco, noco%l_ss, 1, noco%alph, noco%beta, noco%qss, kveclo(:, ikpt), odi, ods, &
            & acofBar, bcofBar, ccofBar)

          ! Rearrange Hamiltonian acofs, bcofs, ccofs ensuring a consistent handling of LOs within the loop the structure compared to LAPWs
          mCoefKv = 0
          lmp = 0
          lm = -1
          do oqn_l = 0, atoms%lmax(1)
            do mqn_m = -oqn_l, oqn_l
              lm = lm + 1
              !p = 1
              lmp = lmp + 1
              mCoefKv(:nobd(ikpt, 1), lmp) = acofBar(:nobd(ikpt, 1), lm, 1)
              !p = 2
              lmp = lmp + 1
              mCoefKv(:nobd(ikpt, 1), lmp) = bcofBar(:nobd(ikpt, 1), lm, 1)
              !LOs
              do p = 3, nRadFun(oqn_l, 1)
                ilo = iloTable(p, oqn_l, 1)
                lmp = lmp + 1
                mCoefKv(:nobd(ikpt, 1), lmp) = ccofBar(mqn_m, :nobd(ikpt, 1), ilo, 1)
              end do
            end do
          end do

          ! Calculate Hamiltonian MT matrix element and overlap matrix with bra or ket varied
          h0MTBv = cmplx(0.0, 0.0)
          s0MTBv = cmplx(0.0, 0.0)
          h0MTKv = cmplx(0.0, 0.0)
          s0MTKv = cmplx(0.0, 0.0)
          call calcHS0MT( atoms, usdus, tdHS0, ikpt, ikpt, 1, 1, nobd(:, 1), nobd, El, mCoefKv, mCoefK(:, :, 1), nRadFun, &
            & iloTable, nlo_atom, s0MTBv, h0MTBv )
          call calcHS0MT( atoms, usdus, tdHS0, ikpt, ikpt, 1, 1, nobd(:, 1), nobd, El, mCoefK(:, :, 1), mCoefKv, nRadFun, &
            & iloTable, nlo_atom, s0MTKv, h0MTKv )

          do iband = 1, nobd(ikpt, 1)
            dynMatPuMt2(idirR, idirC) = dynMatPuMt2(idirR, idirC) + 2 * results%w_iks(iband, ikpt, 1) * ( &
              & h0MTBv(iband, iband) - eig(iband, ikpt, 1) * s0MTBv(iband, iband) &
              & + h0MTKv(iband, iband) - eig(iband, ikpt, 1) * s0MTKv(iband, iband) &
            & )
          end do ! iband
        end do ! idirR
      end do ! idirC
    end do ! ikpt

    if (.false.) then
      write(*, '(a)') 'BarPsi (H - eps) Psi + Psi (H - eps) BarPsi 2nd method'
      write(*, '(3(2(es16.8,1x),3x))') dynMatPuMt2(1, :)
      write(*, '(3(2(es16.8,1x),3x))') dynMatPuMt2(2, :)
      write(*, '(3(2(es16.8,1x),3x))') dynMatPuMt2(3, :)
    end if

    if ( any(abs(dynMatPuMT2(:, :) - dynMatPuMT1(:, :)) > 1e-8) ) passed = .false.

    if (passed) then
      write (logUnit, '(a)')'                                                                 |__ passed!'
    else
      write (logUnit, '(a)')'                                                                 |__ failed!'
      call juDFT_warn('Testing <BarPsi_kn|H - eps|Psi_kn> and <Psi_kn|H - eps|BarPsi_kn> failed', &
        & calledby='testPsiHepsTildePsi', hint='Debug.')
    end if
  end subroutine testPsiHepsTildePsi

  !>-------------------------------------------------------------------------------------------------------------------------------------
  !> @author
  !> Christian-Roman Gerhorst
  !>
  !> @brief
  !> Tests <grPsi|Psi> in the muffin-tin by default. Analyzes the neccesity of a r2-multiplication onto the the wave function in the
  !> MT before applying the gradient. The muffin-tin matrix elements <grGrtPsi|Psi> and <grPsi|GrtPsi> can also be tested.
  !>
  !> @details
  !> This compares in the regular mode <grPsi|Psi>_MT to the numerical MT gradient of the unperturbed density .
  !> This test can also be used to analyzse the impact of multiplying a r2 to the wavefunction in the MT before applying the MT
  !> gradient both before integration, i.e. mesh-resolved and after the integration over the MT mesh.
  !> The sum of <grGrtPsi|Psi> and <grPsi|GrtPsi> in the MT can be compared to the numerical tensor gradient of the unperturbed MT
  !> density but we can avoid the tensor gradient in the production code at the moment.
  !>
  !> @param[in] atoms        : Atoms type, see types.f90.
  !> @param[in] dimens       : Dimension type, see types.f90.
  !> @param[in] kpts         : K-points type, see types.f90.
  !> @param[out] kveclo      : Basis G-vectors of local orbitals.
  !> @param[out] nv          : Number of LAPW G-basis vectors for given k-point.
  !> @param[out] mapGbas     : For various k-points G-basis vectors occur more than once, thus they are only stored once in juPhon. This
  !>                           pointer array contains the right index for GbasVec array to "unfold" G-basis vectors again.
  !> @param[out] gBas        : G-basis vectors
  !> @param[in]  logUnit     : Unit number for juPhon.log.
  !>
  !>-------------------------------------------------------------------------------------------------------------------------------------
  subroutine TestR2orNotWfMtGradNgrNTensGrOvls(atoms, usdus, kpts, cell, sym, dimens, results, lathar, logUnit, clnu_atom, nRadFun,&
      & z0, nv, kveclo, gBas, nobd, mapGbas, rbas1, rbas2, iloTable, mlh_atom, nmem_atom)

    use m_jPConstants, only : iu, fpi
     
    use m_abcof3
    use m_types, only : t_atoms, t_usdus, t_kpts, t_cell, t_sym, t_dimension, t_results, t_sphhar
    use m_intgr, only : intgr3LinIntp
    use m_jpSetupDynMat, only : CalcScalBasfMatElems, CalcVecBasfMatElems
    use m_jpConstants, only : Tmatrix, Tmatrix_transposed
    use mod_juPhonUtils, only : Derivative, CalcChannelsGrFlpNat, CalcChannelsGrGrtFlpNat, CalcChannelsGrR2FlpNat, CalcChannelsGrGrtR2FlpNat, calcGrR2FinLH, CalcGrR2FinSH
    use m_gaunt, only : gaunt1
    use m_juDFT_NOstopNO, only : juDFT_warn

    implicit none

    ! Type parameters
    type(t_atoms),                  intent(in) :: atoms
    type(t_usdus),                  intent(in) :: usdus
    type(t_kpts),                   intent(in) :: kpts
    type(t_cell),                   intent(in) :: cell
    type(t_sym),                    intent(in) :: sym
    type(t_dimension),              intent(in) :: dimens
    type(t_results),                intent(in) :: results
    type(t_sphhar),                 intent(in) :: lathar

    ! Scalar parameter
    integer,                        intent(in) :: logUnit

    ! Array parameters
    complex,                        intent(in) :: clnu_atom(:,0:,:)
    integer,                        intent(in) :: mlh_atom(:, 0:, :)
    integer,                        intent(in) :: nmem_atom(0:, :)
    integer,                        intent(in) :: nRadFun(0:, :)
    MCOMPLEX,                       intent(in) :: z0(:, :, :, :)
    integer,                        intent(in) :: nv(:, :)
    integer,                        intent(in) :: kveclo(:,:)
    integer,                        intent(in) :: gbas(:, :)
    integer,                        intent(in) :: nobd(:, :)
    integer,                        intent(in) :: mapGbas(:, :, :)
    real,                           intent(in) :: rbas1(:, :, 0:, :, :)
    real,                           intent(in) :: rbas2(:, :, 0:, :, :)
    integer,                        intent(in) :: iloTable(:, 0:, :)

    ! Scalar variables
    real                                       :: intgrReal
    real                                       :: intgrImag
    integer                                    :: imesh
    integer                                    :: iqpt
    integer                                    :: oqn_l
    integer                                    :: mqn_m
    integer                                    :: iradf
    integer                                    :: pMaxLocal
    integer                                    :: lm
    integer                                    :: lmp
    integer                                    :: nmat
    integer                                    :: iatom
    integer                                    :: itype
    integer                                    :: ieqat
    integer                                    :: ikpt
    integer                                    :: iband
    integer                                    :: lmpMax
    integer                                    :: nRadFunMax
    integer                                    :: idir
    integer                                    :: mqn_m2Pr
    integer                                    :: lm3Pr
    integer                                    :: oqn_l3Pr
    integer                                    :: mqn_m3Pr
    integer                                    :: mqn_m4Pr
    integer                                    :: iradf3Pr
    integer                                    :: lmp3Pr
    integer                                    :: lmp1Pr
    integer                                    :: oqn_l1Pr
    integer                                    :: mqn_m1Pr
    integer                                    :: iradf1Pr
    integer                                    :: ichan
    integer                                    :: oqn_l4Pr
    real                                       :: gauntCoeff
    integer                                    :: mqn_m2PrR
    integer                                    :: mqn_m2PrC
    integer                                    :: mqn_m5Pr
    integer                                    :: oqn_l5Pr
    integer                                    :: ichanPr
    integer                                    :: idirR
    integer                                    :: idirC
    logical                                    :: simpGradient = .true. ! default: true. only simple gradient (if one wants to check additionaly before integration then activate beforeIntegraton
    logical                                    :: tensorGradient = .false. ! default: false, only tensorGradient(one of the switches beforeIntegration and after Integration has to be activated
    logical                                    :: beforeIntegration = .false. ! default: false. check a mesh resolved quantity
    logical                                    :: afterIntegration = .true. ! default: true. check the quantitiy after the integration has been done.
    logical                                    :: r2Multiplied = .false. ! default: false. multiply r2 to the wavefunctions before they are multiplied with the tensor gradient or do it not
    logical                                    :: output =.false. ! default:false. the output before integration is activated by default because one does only want to have a look here for debugging, this steers the output for the branch after Integration
    logical                                    :: passed = .true. !should be set to true
    integer                                    :: lMesh
    integer                                    :: uMesh

    ! Type variable
    type(od_inp)                               :: odi
    type(od_sym)                               :: ods

    ! Array variables
    real,              allocatable             :: fReal(:)
    real,              allocatable             :: fImag(:)
    real,              allocatable             :: r2Rho0MTlhVal(:, :, :, :)
    complex,           allocatable             :: PsiPsi(:)
    complex,           allocatable             :: varphiPsi(:)
    real,              allocatable             :: varphiVarphi(:, :, :)
    complex,           allocatable             :: varphiHVarphiDummy(:, :, :)
    real,              allocatable             :: varphi1(:, :, :)
    complex,           allocatable             :: hVarphiDummy(:, :, :, :)
    real,              allocatable             :: varphi2(:, :, :)
    integer,           allocatable             :: lmpT(:)
    real,              allocatable             :: r2(:)
    integer,           allocatable             :: ngoprI(:)
    complex,           allocatable             :: a(:, :, :)
    complex,           allocatable             :: b(:, :, :)
    complex,           allocatable             :: bascof_lo(:, :, :, :, :)
    complex,           allocatable             :: ab0cofScl(:)
    real,              allocatable             :: delrVarphi1(:, :, :)
    real,              allocatable             :: delrVarphi2(:, :, :)
    real,              allocatable             :: r2delrVarphi1(:, :, :)
    real,              allocatable             :: r2delrVarphi2(:, :, :)
    real,              allocatable             :: grVarphiCh1(:, :, :, :)
    real,              allocatable             :: grVarphiCh2(:, :, :, :)
    integer,           allocatable             :: grVarphiChLout(:, :)
    integer,           allocatable             :: grVarphiChMout(:, :)
    complex,           allocatable             :: grPsiPsi(:, :)
    complex,           allocatable             :: grVarphiPsi(:)
    complex,           allocatable             :: r2GrRho0MtVal(:, :, :, :)
    complex,           allocatable             :: grVarphiHvarphiDummy(:, :, :)
    real,              allocatable             :: grVarphivarphi(:, :, :)
    complex,           allocatable             :: grPsiPsiMesh(:, :, :)
    complex,           allocatable             :: grPsiPsiNatMesh(:, :, :)
    real,              allocatable             :: grVarphiVarphiMesh(:, :, :, :, :)
    real,              allocatable             :: grVArphiVarphiTest(:, :, :, :)
    complex,           allocatable             :: grPsiPsiTestNat(:, :)
    complex,           allocatable             :: grPsiPsiTest(:, :)
    real,              allocatable             :: intgrdR(:)
    real,              allocatable             :: delrGrVarphiCh1(:, :, :, :)
    real,              allocatable             :: delrGrVarphiCh2(:, :, :, :)
    real,              allocatable             :: r2delrGrVarphiCh1(:, :, :, :)
    real,              allocatable             :: r2delrGrVarphiCh2(:, :, :, :)
    integer,           allocatable             :: grGrtVarphiChLout(:, :)
    integer,           allocatable             :: grGrtVarphiChMout(:, :, :)
    real,              allocatable             :: grGrtVarphiCh1(:, :, :, :, :)
    real,              allocatable             :: grGrtVarphiCh2(:, :, :, :, :)
    real,              allocatable             :: r2grGrtVarphiCh1(:, :, :, :, :)
    real,              allocatable             :: r2grGrtVarphiCh2(:, :, :, :, :)
    real,              allocatable             :: grGrtVarphiVarphiMesh(:, :, :, :, :, :)
    real,              allocatable             :: grVarphiGrtVarphiMesh(:, :, :, :, :, :)
    real,              allocatable             :: grGrtVarphiVarphiTest(:, :, :, :, :)
    real,              allocatable             :: grVarphiGrtVarphiTest(:, :, :, :, :)
    complex,           allocatable             :: grGrtPsiPsiMeshNatNat(:, :, :, :)
    complex,           allocatable             :: grPsiGrtPsiMeshNatNat(:, :, :, :)
    complex,           allocatable             :: grGrtPsiPsiTestNatNat(:, :, :)
    complex,           allocatable             :: grPsiGrtPsiTestNatNat(:, :, :)
    complex,           allocatable             :: varphiPsiTemp(:)
    complex,           allocatable             :: grGrtPsiPsiMesh(:, :, :, :)
    complex,           allocatable             :: grPsiGrtPsiMesh(:, :, :, :)
    complex,           allocatable             :: grGrtPsiPsiTest(:, :, :)
    complex,           allocatable             :: grPsiGrtPsiTest(:, :, :)
    complex,           allocatable             :: r2GrGrtRho0MTVal(:, :, :, :)
    real,              allocatable             :: r2grVarphiCh1(:, :, :, :)
    real,              allocatable             :: r2grVarphiCh2(:, :, :, :)
    real,              allocatable             :: r2varphi1(:, :, :)
    real,              allocatable             :: r2varphi2(:, :, :)
    complex                                    :: grPsiPsiNat(-1:1)
    complex                                    :: grGrtPsiPsiMeshNat(1:3, -1:1)
    complex                                    :: grPsiGrtPsiMeshNat(1:3, -1:1)
    complex                                    :: grGrtPsiPsiTestNat(1:3, -1:1)
    complex                                    :: grPsiGrtPsiTestNat(1:3, -1:1)

    ! The test can also anlalyze the impact of an r2-multiplication or test the matrix elements incorporating a tensor gradient.
    ! Then, the default parameters have to be changed.
    if (.not.simpGradient .or. tensorGradient .or. beforeIntegration .or. .not.afterIntegration .or. r2Multiplied .or. output) then
      call juDFT_warn('Internal control switches have not a default value', &
        & calledby='TestR2orNotWfMtGradNgrNTensGrOvls', hint='Check the switches in the code.')
    end if

    write(logUnit, '(a)') 'Testing <grPsi|Psi>, <grGrtPsi|Psi> and <grPsi|grTPsi>'
    write(logUnit, '(a)') '------------------------------------------------------'

    if (tensorGradient) then
      ! Lower and upper boundary of mesh, when examining quantities before integration for the tensor gradient
      lmesh = 1
      umesh = 100
      if (lmesh > umesh .or. umesh - lmesh > 100) then
        NOstopNO'Mesh window to big, that causes overflow errors'
      end if
    end if


    ! Quantities for initialization
    allocate( lmpT(atoms%ntype) )
    lmpT(:) = 0
    do itype = 1, atoms%ntype
      lmpT(itype) = sum( [ ( (2 * oqn_l + 1)* nRadFun(oqn_l, itype), oqn_l = 0, atoms%lmax(itype) ) ] )
    end do ! itype
    lmpMax     = maxval( lmpT(:) )
    nRadFunMax = maxval( nRadFun(:, :) )

    allocate( r2(atoms%jmtd) )
    allocate( a( dimens%nvd, 0:dimens%lmd, atoms%nat), b(dimens%nvd, 0:dimens%lmd, atoms%nat), &
            & bascof_lo(3, -atoms%llod:atoms%llod, 4 * atoms%llod + 2, atoms%nlod, atoms%nat) )
    allocate( ab0cofScl(lmpMax) )
    allocate(ngoprI(atoms%nat))
    allocate( varphi1(atoms%jmtd, nRadFunMax, 0:atoms%lmaxd), varphi2(atoms%jmtd, nRadFunMax, 0:atoms%lmaxd) )
    allocate( varphiVarphi(lmpMax, lmpMax, atoms%ntype), varphiHvarphiDummy(lmpMax, lmpMax, atoms%nat))
    allocate( varphiPsi(lmpMax), psiPsi(atoms%nat) )
    allocate( hVarphiDummy(2, atoms%jmtd, 0:(atoms%lmaxd + 3)**2 - 1, lmpMax)  )
    allocate( grVarphiCh1(atoms%jmtd, 2, lmpMax, -1:1), grVarphiCh2(atoms%jmtd, 2, lmpMax, -1:1), &
            & grVarphiChLout(2, 0:atoms%lmaxd), grVarphiChMout(-atoms%lmaxd:atoms%lmaxd, -1:1) )
    allocate( grVarphiHvarphiDummy(lmpMax, lmpMax, -1:1) )
    allocate( grVarphiVarphi(lmpMax, lmpMax, -1:1) )
    allocate( grPsiPsi(3, atoms%nat))
    allocate( grVarphiPsi(lmpMax))
    allocate( grPsiPsiNatMesh(atoms%jmtd, 4, -1:1))
    allocate( grVarphiVarphiMesh(atoms%jmtd, 4, lmpMax, lmpMax, -1:1))
    if (tensorGradient) then
      allocate( grGrtVarphiVarphiMesh(lmesh:umesh, 9, lmpMax, lmpMax, -1:1, -1:1))
      allocate( grVarphiGrtVarphiMesh(lmesh:umesh, 9, lmpMax, lmpMax, -1:1, -1:1))
      allocate( grGrtVarphiChLout(4, 0:atoms%lmaxd), grGrtVarphiChMout(-atoms%lmaxd:atoms%lmaxd, -1:1, -1:1), &
              & grGrtVarphiCh1(atoms%jmtd, 4, lmpMax, -1:1, -1:1), grGrtVarphiCh2(atoms%jmtd, 4, lmpMax, -1:1, -1:1) )
      allocate( grGrtVarphiVarphiTest(9, lmpMax, lmpMax, -1:1, -1:1))
      allocate( grVarphiGrtVarphiTest(9, lmpMax, lmpMax, -1:1, -1:1))
      allocate( grGrtPsiPsiMeshNatNat(atoms%jmtd, 9, -1:1, -1:1), grPsiGrtPsiMeshNatNat(atoms%jmtd, 9, -1:1, -1:1))
      allocate( grGrtPsiPsiTestNatNat( 9, -1:1, -1:1), grPsiGrtPsiTestNatNat( 9, -1:1, -1:1) )
      allocate( grGrtPsiPsiMesh(atoms%jmtd, 9, 3, 3), grPsiGrtPsiMesh(atoms%jmtd, 9, 3, 3))
      allocate( grGrtPsiPsiTest( 9, 3, 3), grPsiGrtPsiTest(9, 3, 3))
      allocate( r2GrGrtRho0MTVal(atoms%jmtd, (atoms%lmaxd + 2)**2, atoms%ntype, 3) )
    end if

    allocate( grVarphiVarphiTest(4, lmpMax, lmpMax, -1:1))
    allocate( grPsiPsiTestNat(4, -1:1))
    allocate( grPsiPsiTest(4, 1:3))
    allocate( varphiPsiTemp(lmpMax))
    if (r2Multiplied) then
      allocate( r2grVarphiCh1(atoms%jmtd, 2, lmpMax, -1:1), r2grVarphiCh2(atoms%jmtd, 2, lmpMax, -1:1) )
      allocate( r2delrVarphi1( atoms%jmtd, nRadFunMax, 0:atoms%lmaxd), r2delrVarphi2( atoms%jmtd, nRadFunMax, 0:atoms%lmaxd) )
      allocate( r2varphi1(atoms%jmtd, nRadFunMax, 0:atoms%lmaxd), r2varphi2(atoms%jmtd, nRadFunMax, 0:atoms%lmaxd) )
      if (tensorGradient) then
        allocate( r2delrGrVarphiCh1(atoms%jmtd, 2, lmpMax, -1:1), r2delrGrVarphiCh2(atoms%jmtd, 2, lmpMax, -1:1) )
        allocate( r2grGrtVarphiCh1(atoms%jmtd, 4, lmpMax, -1:1, -1:1), r2grGrtVarphiCh2(atoms%jmtd, 4, lmpMax, -1:1, -1:1) )
      end if
    else
      allocate( delrVarphi1( atoms%jmtd, nRadFunMax, 0:atoms%lmaxd), delrVarphi2( atoms%jmtd, nRadFunMax, 0:atoms%lmaxd) )
      if (tensorGradient) then
        allocate( delrGrVarphiCh1(atoms%jmtd, 2, lmpMax, -1:1), delrGrVarphiCh2(atoms%jmtd, 2, lmpMax, -1:1) )
      end if
    end if

    varphiVarphi(:, :, :) = cmplx(0., 0.)
    varphiHvarphiDummy(:, :, :) = cmplx(0., 0.)
    psiPsi(:) = cmplx(0., 0.)
    hVarphiDummy(:, :, :, :) = cmplx(0., 0.)
    grVarphiHvarphiDummy(:, :, :) = cmplx(0., 0.)
    grPsiPsi(:, :) = cmplx(0., 0.)

    iatom = 0
    do itype = 1, atoms%ntype
      do ieqat = 1, atoms%neq(itype)
        iatom = iatom + 1

        grVarphiVarphi(:, :, :) = cmplx(0., 0.)
        grPsiPsiNat(:) = cmplx(0., 0.)
        grPsiPsiNatMesh(:, :, :) = cmplx(0., 0.)
        grVarphiVarphiTest(:, :, :, :) = cmplx(0., 0.)
        grPsiPsiTestNat(:, :) = cmplx(0., 0.)
        grPsiPsiTest(:, :) = cmplx(0., 0.)
        if (tensorGradient) then
          grGrtPsiPsiMeshNatNat(:, :, :, :) = cmplx(0., 0.)
          grPsiGrtPsiMeshNatNat(:, :, :, :) = cmplx(0., 0.)
          grGrtPsiPsiTestNatNat( :, :, :) = cmplx(0., 0.)
          grPsiGrtPsiTestNatNat( :, :, :) = cmplx(0., 0.)
        end if

        varphi1(:, :, :) = 0.
        varphi2(:, :, :) = 0.
        if (r2Multiplied) then
          r2varphi1(:, :, :) = 0.
          r2varphi2(:, :, :) = 0.
          r2delrVarphi1(:, :, :) = 0.
          r2delrVarphi2(:, :, :) = 0.
        else
          delrVarphi1(:, :, :) = 0.
          delrVarphi2(:, :, :) = 0.
        end if
        do oqn_l = 0, atoms%lmax(itype)
          do iradf = 1, nRadFun(oqn_l, itype)
            do imesh = 1, atoms%jri(itype)
              ! In Fleur the radial solutions u_lp are multiplied by a factor r by default to avoid an additional multiplication of the
              ! Jacobi determinant r^2 in radial integrals given spherical coordinates. We divide this out here so that we do not
              ! have overhead (additional terms) when deriving.
              varphi1(imesh, iradf, oqn_l) = rbas1(imesh, iradf, oqn_l, itype, 1) / atoms%rmsh(imesh, itype)
              varphi2(imesh, iradf, oqn_l) = rbas2(imesh, iradf, oqn_l, itype, 1) / atoms%rmsh(imesh, itype)
              if (r2Multiplied) then
                r2varphi1(imesh, iradf, oqn_l) = rbas1(imesh, iradf, oqn_l, itype, 1) * atoms%rmsh(imesh, itype)
                r2varphi2(imesh, iradf, oqn_l) = rbas2(imesh, iradf, oqn_l, itype, 1) * atoms%rmsh(imesh, itype)
              end if
            end do ! imesh
            ! Calculate the radial derivative of varphi1 and varphi2
            if (r2Multiplied) then
              call Derivative( r2varphi1(1:atoms%jri(itype), iradf, oqn_l), itype, atoms, r2delrVarphi1(1:atoms%jri(itype), iradf, oqn_l) )
              call Derivative( r2varphi2(1:atoms%jri(itype), iradf, oqn_l), itype, atoms, r2delrVarphi2(1:atoms%jri(itype), iradf, oqn_l) )
            else
              call Derivative( varphi1(1:atoms%jri(itype), iradf, oqn_l), itype, atoms, delrVarphi1(1:atoms%jri(itype), iradf, oqn_l) )
              call Derivative( varphi2(1:atoms%jri(itype), iradf, oqn_l), itype, atoms, delrVarphi2(1:atoms%jri(itype), iradf, oqn_l) )
            end if
          end do ! iradf
        end do ! oqn_l

       ! Precalculate radial Jacobi determinant for later integrals
       r2(:) = 0.
       do imesh = 1, atoms%jri(itype)
         r2(imesh) = atoms%rmsh(imesh, itype) * atoms%rmsh(imesh, itype)
       end do ! imesh

      ! Calculate the application of the gradient and the gradient's dyadic product onto the MT basis functions (matching coefficients
      ! have no spatial dependence) and determing its scattering channels.
      grVarphiChLout(:, :) = 0
      grVarphiChMout(:, :) = 0
      grVarphiCh1(:, :, :, :) = 0.
      grVarphiCh2(:, :, :, :) = 0.

    if (r2Multiplied) then
      r2grVarphiCh1(:, :, :, :) = cmplx(0., 0.)
      r2grVarphiCh2(:, :, :, :) = cmplx(0., 0.)
      call CalcChannelsGrR2FlpNat( atoms, itype, nRadFun, r2varphi1, r2varphi2, r2delrVarphi1, r2delrVarphi2, grVarphiChLout, grVarphiChMout, &
                                                                                                        & r2grVarphiCh1, r2grVarphiCh2 )
      do mqn_m2PrR = -1, 1
        lmp = 0
        do oqn_l = 0, atoms%lmax(itype)
          do mqn_m = -oqn_l, oqn_l
            do iradf = 1, nRadFun(oqn_l, itype)
              lmp = lmp + 1
              do ichan = 1, 2
                do imesh = 1, atoms%jmtd
                  grVarphiCh1(imesh, ichan, lmp, mqn_m2PrR) = r2grVarphiCh1(imesh, ichan, lmp, mqn_m2PrR) / atoms%rmsh(imesh, itype)**2
                  grVarphiCh2(imesh, ichan, lmp, mqn_m2PrR) = r2grVarphiCh2(imesh, ichan, lmp, mqn_m2PrR) / atoms%rmsh(imesh, itype)**2
                end do ! imesh
              end do ! ichan
            end do !iradf
          end do ! mqn_m
        end do ! oqn_l
      end do ! mqn_m2PrR
    else
      call CalcChannelsGrFlpNat( atoms, itype, nRadFun, varphi1, varphi2, delrVarphi1, delrVarphi2, grVarphiChLout, grVarphiChMout, &
                                                                                                        & grVarphiCh1, grVarphiCh2 )
    end if

    if (tensorGradient) then

      ! Precalculate partial derivatives of the varphis' gradients in r-direction since it is required two times.
      if (r2Multiplied) then
        r2delrGrVarphiCh1(:, :, :, :) = 0.
        r2delrGrVarphiCh2(:, :, :, :) = 0.
      else
        delrGrVarphiCh1(:, :, :, :) = 0.
        delrGrVarphiCh2(:, :, :, :) = 0.
      end if
      do mqn_m2PrR = -1, 1
        lmp = 0
        do oqn_l = 0, atoms%lmax(itype)
          do mqn_m = -oqn_l, oqn_l
            do iradf = 1, nRadFun(oqn_l, itype)
              lmp = lmp + 1
              do ichan = 1, 2
                if (r2Multiplied) then
                  call Derivative( r2grVarphiCh1(:, ichan, lmp, mqn_m2PrR), itype, atoms, r2delrGrVarphiCh1(:, ichan, lmp, mqn_m2PrR) )
                  call Derivative( r2grVarphiCh2(:, ichan, lmp, mqn_m2PrR), itype, atoms, r2delrGrVarphiCh2(:, ichan, lmp, mqn_m2PrR) )
                else
                  call Derivative( grVarphiCh1(:, ichan, lmp, mqn_m2PrR), itype, atoms, delrGrVarphiCh1(:, ichan, lmp, mqn_m2PrR) )
                  call Derivative( grVarphiCh2(:, ichan, lmp, mqn_m2PrR), itype, atoms, delrGrVarphiCh2(:, ichan, lmp, mqn_m2PrR) )
                end if
              end do ! ichan
            end do !iradf
          end do ! mqn_m
        end do ! oqn_l
      end do ! mqn_m2PrR

      grGrtVarphiChLout(:, :) = 0
      grGrtVarphiChMout(:, :, :) = 0
      grGrtVarphiCh1(:, :, :, :, :) = 0.
      grGrtVarphiCh2(:, :, :, :, :) = 0.
      if (r2Multiplied) then
        r2grGrtVarphiCh1(:, :, :, :, :) = 0.
        r2grGrtVarphiCh2(:, :, :, :, :) = 0.
        call CalcChannelsGrGrtR2FlpNat( atoms, itype, lmpMax, nRadFun, grVarphiChLout, grVarphiChMout, r2grVarphiCh1, r2grVarphiCh2, &
                            & r2delrGrVarphiCh1, r2delrGrVarphiCh2, grGrtVarphiChLout, grGrtVarphiChMout, r2grGrtVarphiCh1, r2grGrtVarphiCh2 )
        do mqn_m2PrC = -1, 1
          do mqn_m2PrR = -1, 1
            lmp = 0
            do oqn_l = 0, atoms%lmax(itype)
              do mqn_m = -oqn_l, oqn_l
                do iradf = 1, nRadFun(oqn_l, itype)
                  lmp = lmp + 1
                  do ichan = 1, 4
                    do imesh = 1, atoms%jri(itype)
                      grGrtVarphiCh1(imesh, ichan, lmp, mqn_m2PrR, mqn_m2PrC) = r2grGrtVarphiCh1(imesh, ichan, lmp, mqn_m2PrR, mqn_m2PrC) / atoms%rmsh(imesh, itype)**2
                      grGrtVarphiCh2(imesh, ichan, lmp, mqn_m2PrR, mqn_m2PrC) = r2grGrtVarphiCh2(imesh, ichan, lmp, mqn_m2PrR, mqn_m2PrC) / atoms%rmsh(imesh, itype)**2
                    end do ! imesh
                  end do ! ichan
                end do !iradf
              end do ! mqn_m
            end do ! oqn_l
          end do ! mqn_m2PrR
        end do ! mqn_m2PrC
      else

        call CalcChannelsGrGrtFlpNat( atoms, itype, lmpMax, nRadFun, grVarphiChLout, grVarphiChMout, grVarphiCh1, grVarphiCh2, &
                          & delrGrVarphiCh1, delrGrVarphiCh2, grGrtVarphiChLout, grGrtVarphiChMout, grGrtVarphiCh1, grGrtVarphiCh2 )
      end if

    end if

      ! These functions are for the case in which the quantities are not mesh or lm resolved.
      ! Project varphi on varphi, the rest are dummy arguments
      call CalcScalBasfMatElems( atoms, itype, iatom, nRadFun, r2, varphi1, varphi2, hVarphiDummy, varphiVarphi,    &
                                                                                                     & varphiHvarphiDummy(:, :, :) )

      ! Project grVarphi onto varphi, the rest are dummy arguments
      call CalcVecBasfMatElems( atoms, itype, 2, nRadFun, r2, grVarphiChLout, grVarphiChMout, varphi1, varphi2, &
          & grVarPhiCh1, grVarphiCh2, hVarphiDummy, grVarphiVarphi, grVarphiHvarphiDummy )

      allocate( intgrdR(atoms%jri(itype)) )

      ! As in calcscalbasfmatelems and calcvecbasfmatelems the integral over the MT mesh is already performed, we cannot use it to
      ! analyze the impact of r2 for <grVarphi|varphi> before integration. Therefore the integration takes here place, again.
      ! Especially in contrast to CalcScalBasfMatElems and CalcVecBasfMatElems, we here have a lm resolved matrix element of
      ! grVarphiVarphi that we can compare with the gradient of the density and it is mesh resolved.
      !todo maybe we do things here redundantly
      if (simpGradient) then
        grVarphiVarphiMesh = cmplx(0., 0.)
        do mqn_m2Pr = -1, 1
          lmp3Pr = 0
          do oqn_l3Pr = 0, atoms%lmax(itype)
            do mqn_m3Pr = -oqn_l3Pr, oqn_l3Pr
              mqn_m4Pr = grVarphiChMout(mqn_m3Pr, mqn_m2Pr)
              do iradf3Pr = 1, nRadFun(oqn_l3Pr, itype)
                lmp3Pr = lmp3Pr + 1
                lmp1Pr = 0
                do oqn_l1Pr = 0, atoms%lmax(itype)
                  do mqn_m1Pr = -oqn_l1Pr, oqn_l1Pr
                    do iradf1Pr = 1, nRadFun(oqn_l1Pr, itype)
                      lmp1Pr = lmp1Pr + 1
                      do ichan = 1, 2
                        oqn_l4Pr = grVarphiChLout(ichan, oqn_l3Pr)
                        if ( ( abs(mqn_m4Pr) > oqn_l4Pr ) .or. ( oqn_l4Pr < 0 ) .or. ( oqn_l4Pr > atoms%lmax(itype) ) ) cycle
                        oqn_l = 1
                        do mqn_m = -oqn_l, oqn_l
                          lm = oqn_l * (oqn_l  + 1) + 1 + mqn_m
                          intgrdR(:) = 0.
                          gauntCoeff = gaunt1(oqn_l1Pr, oqn_l, oqn_l4Pr, mqn_m1Pr, mqn_m, mqn_m4Pr, atoms%lmax(itype))
                          do imesh = 1, atoms%jri(itype)
                            intgrdR(imesh) = r2(imesh) * ( grVarPhiCh1(imesh, ichan, lmp3Pr, mqn_m2Pr) * varphi1(imesh, iradf1Pr, oqn_l1Pr)       &
                              & + grVarPhiCh2(imesh, ichan, lmp3Pr, mqn_m2Pr) * varphi2(imesh, iradf1Pr, oqn_l1Pr) ) * gauntCoeff
                          if (beforeIntegration) then
                            grVarphiVarphiMesh(imesh, lm, lmp3Pr, lmp1Pr, mqn_m2Pr) = grVarphiVarphiMesh(imesh, lm, lmp3Pr, lmp1Pr, mqn_m2Pr) + ( grVarPhiCh1(imesh, ichan, lmp3Pr, mqn_m2Pr) * varphi1(imesh, iradf1Pr, oqn_l1Pr)       &
                              & + grVarPhiCh2(imesh, ichan, lmp3Pr, mqn_m2Pr) * varphi2(imesh, iradf1Pr, oqn_l1Pr) ) * gauntCoeff
                          end if
                          end do ! imesh
                          call intgr3LinIntp(intgrdR(:), atoms%rmsh(1, itype), atoms%dx(itype), atoms%jri(itype), intgrReal, 1)
                          ! < grVarphi | varphi >
                          grVarphiVarphiTest(lm, lmp3Pr, lmp1Pr, mqn_m2Pr) = grVarphiVarphiTest(lm, lmp3Pr, lmp1Pr, mqn_m2Pr) +intgrReal
                        end do ! mqn_m
                      end do ! ichan
                    end do ! iradf1Pr
                  end do ! mqn_m1Pr
                end do ! oqn_l1Pr
              end do ! iradf3Pr
            end do ! mqn_m3Pr
          end do ! oqn_l3Pr
        end do ! mqn_m2Pr
      end if

      ! The same argumentation as the last comment for the singelGradient. Manual lm-resoveld mesh resolved setup of <grGrtPsi|Psi>
      ! integrand.
      if (tensorGradient .and. beforeIntegration) then
        grGrtVarphiVarphiMesh(:, :, :, :, :, :) = cmplx(0., 0.)
        do mqn_m2PrC = -1, 1
          do mqn_m2PrR = -1, 1
            lmp1Pr = 0
            do oqn_l1Pr = 0, atoms%lmax(itype)
              do mqn_m1Pr = -oqn_l1Pr, oqn_l1Pr
                do iradf1Pr = 1, nRadFun(oqn_l1Pr, itype)
                  lmp1Pr = lmp1Pr + 1
                  lmp3Pr = 0
                  do oqn_l3Pr = 0, atoms%lmax(itype)
                    do mqn_m3Pr = -oqn_l3Pr, oqn_l3Pr
                      mqn_m4Pr = grGrtVarphiChMout(mqn_m3Pr, mqn_m2PrR, mqn_m2PrC)
                      do iradf3Pr = 1, nRadFun(oqn_l3Pr, itype)
                        lmp3Pr = lmp3Pr + 1
                        do ichan = 1, 4
                          oqn_l4Pr = grGrtVarphiChLout(ichan, oqn_l3Pr)
                          !todo is lmax cutoff correct?
                          if ( ( abs(mqn_m4Pr) > oqn_l4Pr ) .or. ( oqn_l4Pr < 0 ) .or. ( oqn_l4Pr > atoms%lmax(itype) ) ) cycle
                          do oqn_l = 0, 2, 2
                            do mqn_m = -oqn_l, oqn_l
                              lm = oqn_l * (oqn_l  + 1) + 1 + mqn_m
                              gauntCoeff = gaunt1(oqn_l1Pr, oqn_l, oqn_l4Pr, mqn_m1Pr, mqn_m, mqn_m4Pr, atoms%lmax(itype))
                              do imesh = lMesh, uMesh
                                grGrtVarphiVarphiMesh(imesh, lm, lmp3Pr, lmp1Pr, mqn_m2PrR, mqn_m2PrC) = grGrtVarphiVarphiMesh(imesh, lm, lmp3Pr, lmp1Pr, mqn_m2PrR, mqn_m2PrC) + ( grGrtVarPhiCh1(5*imesh, ichan, lmp3Pr, mqn_m2PrR, mqn_m2PrC) * varphi1(5*imesh, iradf1Pr, oqn_l1Pr)       &
                                  & + grGrtVarPhiCh2(5*imesh, ichan, lmp3Pr, mqn_m2PrR, mqn_m2PrC) * varphi2(5*imesh, iradf1Pr, oqn_l1Pr) ) * gauntCoeff
                              end do ! imesh
                            end do ! mqn_m
                          end do ! oqn_l
                        end do ! ichan
                      end do ! iradf1Pr
                    end do ! mqn_m1Pr
                  end do ! oqn_l1Pr
                end do ! iradf3Pr
              end do ! mqn_m3Pr
            end do ! oqn_l3Pr
          end do ! mqn_m2PrR
        end do ! mqn_m2PrC
      end if ! beforeIntegration


      ! The same argumentation as the last comment for the singelGradient. Manual lm-resolved integration of grGrtVarphiVarphi to
      ! compare it later with the lm-resolved density gradient integration.
      if (tensorGradient .and. afterIntegration) then
        grGrtVarphiVarphiTest( :, :, :, :, :) = cmplx(0., 0.)
        do mqn_m2PrC = -1, 1
          do mqn_m2PrR = -1, 1
            lmp1Pr = 0
            do oqn_l1Pr = 0, atoms%lmax(itype)
              do mqn_m1Pr = -oqn_l1Pr, oqn_l1Pr
                do iradf1Pr = 1, nRadFun(oqn_l1Pr, itype)
                  lmp1Pr = lmp1Pr + 1
                  lmp3Pr = 0
                  do oqn_l3Pr = 0, atoms%lmax(itype)
                    do mqn_m3Pr = -oqn_l3Pr, oqn_l3Pr
                      mqn_m4Pr = grGrtVarphiChMout(mqn_m3Pr, mqn_m2PrR, mqn_m2PrC)
                      do iradf3Pr = 1, nRadFun(oqn_l3Pr, itype)
                        lmp3Pr = lmp3Pr + 1
                        do ichan = 1, 4
                          oqn_l4Pr = grGrtVarphiChLout(ichan, oqn_l3Pr)
                          !todo is lmax cutoff correct?
                          if ( ( abs(mqn_m4Pr) > oqn_l4Pr ) .or. ( oqn_l4Pr < 0 ) .or. ( oqn_l4Pr > atoms%lmax(itype) ) ) cycle
                          do oqn_l = 0, 2, 2
                            do mqn_m = -oqn_l, oqn_l
                              lm = oqn_l * (oqn_l  + 1) + 1 + mqn_m
                              intgrdR(:) = 0.
                              gauntCoeff = gaunt1(oqn_l1Pr, oqn_l, oqn_l4Pr, mqn_m1Pr, mqn_m, mqn_m4Pr, atoms%lmax(itype))
                              do imesh = 1, atoms%jri(itype)
                                intgrdR(imesh) = r2(imesh) * ( grGrtVarPhiCh1(imesh, ichan, lmp3Pr, mqn_m2PrR, mqn_m2PrC) * varphi1(imesh, iradf1Pr, oqn_l1Pr)       &
                                  & + grGrtVarPhiCh2(imesh, ichan, lmp3Pr, mqn_m2PrR, mqn_m2PrC) * varphi2(imesh, iradf1Pr, oqn_l1Pr) ) * gauntCoeff
                              end do ! imesh
                              call intgr3LinIntp(intgrdR(:), atoms%rmsh(1, itype), atoms%dx(itype), atoms%jri(itype), intgrReal, 1)
                              ! < grVarphi | varphi >
                              grGrtVarphiVarphiTest(lm, lmp3Pr, lmp1Pr, mqn_m2PrR, mqn_m2PrC) = grGrtVarphiVarphiTest(lm, lmp3Pr, lmp1Pr, mqn_m2PrR, mqn_m2PrC) + intgrReal
                            end do ! mqn_m
                          end do ! oqn_l
                        end do ! ichan
                      end do ! iradf1Pr
                    end do ! mqn_m1Pr
                  end do ! oqn_l1Pr
                end do ! iradf3Pr
              end do ! mqn_m3Pr
            end do ! oqn_l3Pr
          end do ! mqn_m2PrR
        end do ! mqn_m2PrC

      end if ! afterIntegration


      ! The same argumentation as the last comment for the singelGradient. Manual lm-resolved and mesh resolved setup of
      ! <grPsi|GrtPsi> integrand
      if (tensorGradient .and. beforeIntegration) then
        grVarphiGrtVarphiMesh(:, :, :, :, :, :) = cmplx(0., 0.)
        do mqn_m2PrC = -1, 1
          do mqn_m2PrR = -1, 1
          lmp1Pr = 0
          do oqn_l1Pr = 0, atoms%lmax(itype)
            do mqn_m1Pr = -oqn_l1Pr, oqn_l1Pr
              mqn_m5Pr = grVarphiChMout(mqn_m1Pr, mqn_m2PrC)
              do iradf1Pr = 1, nRadFun(oqn_l1Pr, itype)
                lmp1Pr = lmp1Pr + 1
                lmp3Pr = 0
                do oqn_l3Pr = 0, atoms%lmax(itype)
                  do mqn_m3Pr = -oqn_l3Pr, oqn_l3Pr
                    mqn_m4Pr = grVarphiChMout(mqn_m3Pr, mqn_m2PrR)
                    do iradf3Pr = 1, nRadFun(oqn_l3Pr, itype)
                      lmp3Pr = lmp3Pr + 1
                        do ichanPr = 1, 2
                          oqn_l4Pr = grVarphiChLout(ichanPr, oqn_l3Pr)
                          if ( ( abs(mqn_m4Pr) > oqn_l4Pr ) .or. ( oqn_l4Pr < 0 ) .or. ( oqn_l4Pr > atoms%lmax(itype) ) ) cycle
                          do ichan = 1, 2
                            oqn_l5Pr = grVarphiChLout(ichan, oqn_l1Pr)
                            if ( ( abs(mqn_m5Pr) > oqn_l5Pr ) .or. ( oqn_l5Pr < 0 ) .or. ( oqn_l5Pr > atoms%lmax(itype) ) ) cycle
                            !todo is lmax cutoff correct?
                            do oqn_l = 0, 2, 2
                              do mqn_m = -oqn_l, oqn_l
                                lm = oqn_l * (oqn_l  + 1) + 1 + mqn_m
                                intgrdR(:) = 0.
                                gauntCoeff = gaunt1(oqn_l5Pr, oqn_l, oqn_l4Pr, mqn_m5Pr, mqn_m, mqn_m4Pr, atoms%lmax(itype))
                                do imesh = lMesh, uMesh
                                  grVarphiGrtVarphiMesh(imesh, lm, lmp3Pr, lmp1Pr, mqn_m2PrR, mqn_m2PrC) = grVarphiGrtVarphiMesh(imesh, lm, lmp3Pr, lmp1Pr, mqn_m2PrR, mqn_m2PrC) + ( grVarPhiCh1(imesh, ichanPr, lmp3Pr, mqn_m2PrR) *  grVarPhiCh1(imesh, ichan, lmp1Pr, mqn_m2PrC)      &
                                    & + grVarPhiCh2(imesh, ichanPr, lmp3Pr, mqn_m2PrR) * grVarPhiCh2(imesh, ichan, lmp1Pr, mqn_m2PrC) ) * gauntCoeff
                                end do ! imesh
                              end do ! mqn_m
                            end do ! oqn_l
                          end do ! ichan
                        end do ! ichanPr
                      end do ! iradf1Pr
                    end do ! mqn_m1Pr
                  end do ! oqn_l1Pr
                end do ! iradf3Pr
              end do ! mqn_m3Pr
            end do ! oqn_l3Pr
          end do ! mqn_m2PrR
        end do ! mqn_m2PrC
      end if

      ! The same argumentation as the last comment for the singelGradient. Manual lm-resolved integration of <grPsi|GrtPsi>
      if (tensorGradient .and. afterIntegration) then
        grVarphiGrtVarphiTest( :, :, :, :, :) = cmplx(0., 0.)
        do mqn_m2PrC = -1, 1
          do mqn_m2PrR = -1, 1
          lmp1Pr = 0
          do oqn_l1Pr = 0, atoms%lmax(itype)
            do mqn_m1Pr = -oqn_l1Pr, oqn_l1Pr
              mqn_m5Pr = grVarphiChMout(mqn_m1Pr, mqn_m2PrC)
              do iradf1Pr = 1, nRadFun(oqn_l1Pr, itype)
                lmp1Pr = lmp1Pr + 1
                lmp3Pr = 0
                do oqn_l3Pr = 0, atoms%lmax(itype)
                  do mqn_m3Pr = -oqn_l3Pr, oqn_l3Pr
                    mqn_m4Pr = grVarphiChMout(mqn_m3Pr, mqn_m2PrR)
                    do iradf3Pr = 1, nRadFun(oqn_l3Pr, itype)
                      lmp3Pr = lmp3Pr + 1
                        do ichanPr = 1, 2
                          oqn_l4Pr = grVarphiChLout(ichanPr, oqn_l3Pr)
                          if ( ( abs(mqn_m4Pr) > oqn_l4Pr ) .or. ( oqn_l4Pr < 0 ) .or. ( oqn_l4Pr > atoms%lmax(itype) ) ) cycle
                          do ichan = 1, 2
                            oqn_l5Pr = grVarphiChLout(ichan, oqn_l1Pr)
                            if ( ( abs(mqn_m5Pr) > oqn_l5Pr ) .or. ( oqn_l5Pr < 0 ) .or. ( oqn_l5Pr > atoms%lmax(itype) ) ) cycle
                            !todo is lmax cutoff correct?
                            do oqn_l = 0, 2, 2
                              do mqn_m = -oqn_l, oqn_l
                                lm = oqn_l * (oqn_l  + 1) + 1 + mqn_m
                                intgrdR(:) = 0.
                                gauntCoeff = gaunt1(oqn_l5Pr, oqn_l, oqn_l4Pr, mqn_m5Pr, mqn_m, mqn_m4Pr, atoms%lmax(itype))
                                do imesh = 1, atoms%jri(itype)
                                  intgrdR(imesh) = r2(imesh) * ( grVarPhiCh1(imesh, ichanPr, lmp3Pr, mqn_m2PrR) * grVarPhiCh1(imesh, ichan, lmp1Pr, mqn_m2PrC)       &
                                    & + grVarPhiCh2(imesh, ichanPr, lmp3Pr, mqn_m2PrR) * grVarPhiCh2(imesh, ichan, lmp1Pr, mqn_m2PrC) ) * gauntCoeff
                                end do ! imesh
                                !call intgr3NoIntp(intgrdR(:), atoms%rmsh(1, itype), atoms%dx(itype), atoms%jri(itype), integralR)
                                call intgr3LinIntp(intgrdR(:), atoms%rmsh(1, itype), atoms%dx(itype), atoms%jri(itype), intgrReal, 1)
                                !call intgr3(intgrdR(:), atoms%rmsh(1, itype), atoms%dx(itype), atoms%jri(itype), integralR)
                                ! < grVarphi | varphi >
                                grVarphiGrtVarphiTest(lm, lmp3Pr, lmp1Pr, mqn_m2PrR, mqn_m2PrC) = grVarphiGrtVarphiTest(lm, lmp3Pr, lmp1Pr, mqn_m2PrR, mqn_m2PrC) + intgrReal
                              end do ! mqn_m
                            end do ! oqn_l
                          end do ! ichan
                        end do ! ichanPr
                      end do ! iradf1Pr
                    end do ! mqn_m1Pr
                  end do ! oqn_l1Pr
                end do ! iradf3Pr
              end do ! mqn_m3Pr
            end do ! oqn_l3Pr
          end do ! mqn_m2PrR
        end do ! mqn_m2PrC
      end if ! afterIntegration


      PsiPsi(:) = cmplx(0., 0.)
      ngoprI(:) = 1
      iqpt = 1
      do ikpt = 1, kpts%nkpt


        ! Generate large matching coefficients
        nmat = nv(1, ikpt) + atoms%nlotot
        a(:, :, :) = cmplx(0.0, 0.0)
        b(:, :, :) = cmplx(0.0, 0.0)
        bascof_lo(:, :, :, :, :) = cmplx(0.0, 0.0)
        call abcof3( atoms%lmaxd, atoms%ntype, atoms%nat, sym%nop, dimens%nvd, dimens%jspd, 1, dimens%lmd, dimens%nbasfcn, &
          & atoms%llod, atoms%nlod, atoms%nlotot, sym%invtab, atoms%ntype, sym%mrot, ngoprI, atoms%taual, atoms%neq, atoms%lmax, &
          & atoms%rmt, cell%omtil, cell%bmat, cell%bbmat, kpts%bk(:, ikpt), gbas(1, mapGbas(:nv(1, ikpt), ikpt, 1)), &
          & gbas(2, mapGbas(:nv(1, ikpt), ikpt, 1)), gbas(3, mapGbas(:nv(1, ikpt), ikpt, 1)), nv(:, ikpt), nmat, &
          & usdus%us, usdus%dus, usdus%uds, usdus%duds, usdus%ddn, atoms%invsat, sym%invsatnr, usdus%ulos, usdus%uulon, usdus%dulon, &
          & usdus%dulos, atoms%llo, atoms%nlo, atoms%l_dulo, atoms%lapw_l, kveclo(:,ikpt), odi, ods, a, b, bascof_lo )

        do iband = 1, nobd(ikpt, 1)
          ab0cofScl(:) = cmplx(0., 0.)
          lmp = 0
          lm  = 0
          do oqn_l = 0, atoms%lmax(itype)
            do mqn_m = - oqn_l, oqn_l
              pMaxLocal = nRadFun(oqn_l, itype)
              ! p = 1
              ab0cofScl(lmp + 1) = iu**oqn_l * dot_product( conjg(z0(:nv(1, ikpt), iband, ikpt, 1)), a(:nv(1, ikpt), lm, iatom) )
              ! p = 2
              ab0cofScl(lmp + 2) = iu**oqn_l * dot_product( conjg(z0(:nv(1, ikpt), iband, ikpt, 1)), b(:nv(1, ikpt), lm, iatom) )
              ! Add LO contributions
              do iradf = 3, pMaxLocal
                ! p = 1
                ab0cofScl(lmp + 1) = ab0cofScl(lmp + 1) + iu**oqn_l * &
                  & dot_product( conjg(z0(nv(1, ikpt) + 1:nv(1, ikpt) + atoms%nlotot, iband, ikpt, 1)), &
                                                       & bascof_lo(1, mqn_m, :atoms%nlotot, iloTable(iradf, oqn_l, itype), iatom) )
                ! p = 2
                ab0cofScl(lmp + 2) = ab0cofScl(lmp + 2) + iu**oqn_l * &
                  & dot_product( conjg(z0(nv(1, ikpt) + 1:nv(1, ikpt) + atoms%nlotot, iband, ikpt, 1)), &
                                                       & bascof_lo(2, mqn_m, :atoms%nlotot, iloTable(iradf, oqn_l, itype), iatom) )
                ! 2 < p < LOs for that l and that atom type
                ab0cofScl(lmp + iradf) = iu**oqn_l * dot_product( conjg(z0(nv(1, ikpt) + 1:nv(1, ikpt) + atoms%nlotot, iband, ikpt, 1)), &
                                                       & bascof_lo(3, mqn_m, :atoms%nlotot, iloTable(iradf, oqn_l, itype), iatom) )
              end do ! iradf
              lm = lm + 1
              lmp = lmp + pMaxLocal
            end do ! mqn_m
          end do !oqn_l

          ! Calculate Psi* Psi overlap with lm summed and mesh integrated
          varphiPsi(:) = cmplx(0., 0.)
          varphiPsi(:) = matmul( varphiVarphi(1:lmpT(itype), 1:lmpT(itype), itype), ab0cofScl(1:lmpT(itype)) )
          PsiPsi(iatom) = PsiPsi(iatom) + results%w_iks(iband, ikpt, 1) * dot_product(ab0cofScl(1:lmpT(itype)), varphiPsi(1:lmpT(itype)))

          ! Calculate grPsi* Psi overlap with lm summed and mesh integrated
          do mqn_m2Pr = -1, 1
            grVarphiPsi(:) = cmplx(0., 0.)
            grVarphiPsi(1:lmpT(itype)) = matmul( grVarphiVarphi(1:lmpT(itype), 1:lmpT(itype), itype), ab0cofScl(1:lmpT(itype)) )
            grPsiPsiNat(mqn_m2Pr) = grPsiPsiNat(mqn_m2Pr) + results%w_iks(iband, ikpt, 1) &
                                                             & * dot_product( ab0cofScl(1:lmpT(itype)), grVarphiPsi(1:lmpT(itype)) )

            ! Calculate grPsi* Psi overlap integrand (mesh-resolved) and lm resolved for largest l = 1 component
            ! Calculate grPsi* Psi overlap integrand lm resolved for largest l = 1 component
            if (simpGradient) then
              oqn_l = 1
              do mqn_m = -oqn_l, oqn_l
                lm = oqn_l * (oqn_l + 1) + 1 + mqn_m
                if (beforeIntegration) then
                  do imesh = 1, atoms%jri(itype)
                    grVarphiPsi(:) = cmplx(0., 0.)
                    grVarphiPsi(1:lmpT(itype)) = matmul( grVarphiVarphiMesh(imesh, lm, 1:lmpT(itype), 1:lmpT(itype), mqn_m2Pr), ab0cofScl(1:lmpT(itype)) )
                    grPsiPsiNatMesh(imesh, lm, mqn_m2Pr) = grPsiPsiNatMesh(imesh, lm, mqn_m2Pr) + results%w_iks(iband, ikpt, 1) * dot_product( ab0cofScl(1:lmpT(itype)), grVarphiPsi(1:lmpT(itype)) )
                  end do ! imesh
                end if
                grVarphiPsi(:) = cmplx(0., 0.)
                grVarphiPsi(1:lmpT(itype)) = matmul( grVarphiVarphiTest(lm, 1:lmpT(itype), 1:lmpT(itype), mqn_m2Pr), ab0cofScl(1:lmpT(itype)) )
                grPsiPsiTestNat(lm, mqn_m2Pr) = grPsiPsiTestNat(lm, mqn_m2Pr) + results%w_iks(iband, ikpt, 1) * dot_product( ab0cofScl(1:lmpT(itype)), grVarphiPsi(1:lmpT(itype)) )
              end do ! mqn_m
            end if
          end do ! mqn_m2Pr

          ! Calculate lm-resolved and mesh resolved integrand of <grGrtPsi|Psi>
          if (tensorGradient .and. beforeIntegration) then
            do mqn_m2PrC = -1, 1
              do mqn_m2PrR = -1, 1
                do oqn_l = 0, 2, 2
                  do mqn_m = -oqn_l, oqn_l
                    lm = oqn_l * (oqn_l + 1) + 1 + mqn_m
                    do imesh = lMesh, uMesh
                      varphiPsiTemp(:) = cmplx(0., 0.)
                      varphiPsiTemp(1:lmpT(itype)) = matmul(grGrtVarphiVarphiMesh(imesh, lm, 1:lmpT(itype), 1:lmpT(itype), mqn_m2PrR, mqn_m2PrC), ab0cofScl(1:lmpT(itype)))
                      grGrtPsiPsiMeshNatNat(imesh, lm, mqn_m2PrR, mqn_m2PrC) = grGrtPsiPsiMeshNatNat(imesh, lm, mqn_m2PrR, mqn_m2PrC) + results%w_iks(iband, ikpt, 1) * dot_product(ab0cofScl(1:lmpT(itype)), varphiPsiTemp(1:lmpT(itype)))
                      varphiPsiTemp(:) = cmplx(0., 0.)
                      varphiPsiTemp(1:lmpT(itype)) = matmul(grVarphiGrtVarphiMesh(imesh, lm, 1:lmpT(itype), 1:lmpT(itype), mqn_m2PrR, mqn_m2PrC), ab0cofScl(1:lmpT(itype)))
                      grPsiGrtPsiMeshNatNat(imesh, lm, mqn_m2PrR, mqn_m2PrC) = grPsiGrtPsiMeshNatNat(imesh, lm, mqn_m2PrR, mqn_m2PrC) + results%w_iks(iband, ikpt, 1) * dot_product(ab0cofScl(1:lmpT(itype)), varphiPsiTemp(1:lmpT(itype)))
                    end do ! imesh
                  end do ! mqn_m
                end do ! oqn_l
              end do ! mqn_m2PrR
            end do ! mqn_m2PrC
          end if ! beforeIntegration

          ! Calculate lm-resolved integral of <grGrtPsi|Psi>
          if (tensorGradient .and. afterIntegration) then
            do mqn_m2PrC = -1, 1
              do mqn_m2PrR = -1, 1
                do oqn_l = 0, 2, 2
                  do mqn_m = -oqn_l, oqn_l
                    lm = oqn_l * (oqn_l + 1) + 1 + mqn_m
                    varphiPsiTemp(:) = cmplx(0., 0.)
                    varphiPsiTemp(1:lmpT(itype)) = matmul(grGrtVarphiVarphiTest(lm, 1:lmpT(itype), 1:lmpT(itype), mqn_m2PrR, mqn_m2PrC), ab0cofScl(1:lmpT(itype)))
                    grGrtPsiPsiTestNatNat(lm, mqn_m2PrR, mqn_m2PrC) = grGrtPsiPsiTestNatNat(lm, mqn_m2PrR, mqn_m2PrC) + results%w_iks(iband, ikpt, 1) * dot_product(ab0cofScl(1:lmpT(itype)), varphiPsiTemp(1:lmpT(itype)))
                    varphiPsiTemp(:) = cmplx(0., 0.)
                    varphiPsiTemp(1:lmpT(itype)) = matmul(grVarphiGrtVarphiTest(lm, 1:lmpT(itype), 1:lmpT(itype), mqn_m2PrR, mqn_m2PrC), ab0cofScl(1:lmpT(itype)))
                    grPsiGrtPsiTestNatNat(lm, mqn_m2PrR, mqn_m2PrC) = grPsiGrtPsiTestNatNat(lm, mqn_m2PrR, mqn_m2PrC) + results%w_iks(iband, ikpt, 1) * dot_product( ab0cofScl(1:lmpT(itype)), varphiPsiTemp(1:lmpT(itype)) )
                  end do ! mqn_m
                end do ! oqn_l
              end do ! mqn_m2PrR
            end do ! mqn_m2PrC
          end if ! afterIntegration

        end do ! iband

      end do ! ikpt

      ! Transform into cartesian coordinates
      if (simpGradient) then
        allocate(grPsiPsiMesh(atoms%jmtd, 4, 3) )
        grPsiPsiMesh(:, :, :) = cmplx(0., 0.)
        oqn_l = 1
        do mqn_m = -oqn_l, oqn_l
          lm = oqn_l * (oqn_l + 1) + 1 + mqn_m
          if (beforeIntegration) then
            do imesh = 1, atoms%jri(itype)
              grPsiPsiMesh(imesh, lm, 1:3) = matmul(conjg(Tmatrix(1:3, 1:3)), grPsiPsiNatMesh(imesh, lm, -1:1))
            end do ! imesh
          end if
          grPsiPsiTest(lm, 1:3) = matmul(conjg(Tmatrix(1:3, 1:3)), grPsiPsiTestNat(lm, -1:1))
        end do ! mqn_m
      end if

      grPsiPsi(1:3, iatom) = matmul(conjg(Tmatrix(1:3, 1:3)), grPsiPsiNat(-1:1))


      if (tensorGradient) then
        grGrtPsiPsiMesh(:, :, :, :) = cmplx(0., 0.)
        grPsiGrtPsiMesh(:, :, :, :) = cmplx(0., 0.)
        grPsiGrtPsiTest(:, :, :) = cmplx(0., 0.)
        grGrtPsiPsiTest(:, :, :) = cmplx(0., 0.)
        grGrtPsiPsiMeshNat(:, :) = cmplx(0., 0.)
        grPsiGrtPsiMeshNat(:, :) = cmplx(0., 0.)
        grPsiGrtPsiTestNat(:, :) = cmplx(0., 0.)
        grGrtPsiPsiTestNat(:, :) = cmplx(0., 0.)
        if (beforeIntegration) then
          do oqn_l = 0, 2, 2
            do mqn_m = -oqn_l, oqn_l
              lm = oqn_l * (oqn_l + 1) + 1 + mqn_m
              do imesh = lMesh, uMesh
                grGrtPsiPsiMeshNat(1:3, -1:1) = matmul( conjg(Tmatrix(1:3, 1:3)), grGrtPsiPsiMeshNatNat(imesh, lm, -1:1, -1:1))
                grGrtPsiPsiMesh(imesh, lm, 1:3, 1:3) = matmul( grGrtPsiPsiMeshNat(1:3, -1:1), conjg(Tmatrix_transposed(1:3, 1:3)) )
                grPsiGrtPsiMeshNat(1:3, -1:1) = matmul( conjg(Tmatrix(1:3, 1:3)), grPsiGrtPsiMeshNatNat(imesh, lm, -1:1, -1:1))
                grPsiGrtPsiMesh(imesh, lm, 1:3, 1:3) = matmul( grPsiGrtPsiMeshNat(1:3, -1:1), Tmatrix_transposed(1:3, 1:3) )
              end do ! imesh
            end do ! mqn_m
          end do ! oqn_l
        end if

        if (afterIntegration) then
          do oqn_l = 0, 2, 2
            do mqn_m = -oqn_l, oqn_l
              lm = oqn_l * (oqn_l + 1) + 1 + mqn_m
              grGrtPsiPsiTestNat(1:3, -1:1) = matmul( conjg(Tmatrix(1:3, 1:3)), grGrtPsiPsiTestNatNat(lm, -1:1, -1:1))
              grGrtPsiPsiTest(lm, 1:3, 1:3) = matmul( grGrtPsiPsiTestNat(1:3, -1:1), conjg(Tmatrix_transposed(1:3, 1:3)) )
              grPsiGrtPsiTestNat(1:3, -1:1) = matmul( conjg(Tmatrix(1:3, 1:3)), grPsiGrtPsiTestNatNat(lm, -1:1, -1:1))
              grPsiGrtPsiTest(lm, 1:3, 1:3) = matmul( grPsiGrtPsiTestNat(1:3, -1:1), Tmatrix_transposed(1:3, 1:3) )
            end do ! mqn_m
          end do ! oqn_l
        end if
      end if

      ! Read in valence density to compare it to the output of calcscalmatelem
      allocate( r2Rho0MTlhVal( atoms%jmtd, 0:lathar%nlhd, atoms%ntype, 1) )
      r2Rho0MtlhVal(:, :, :, :) = cmplx(0., 0.)
      rewind(1040)
      read(1040) r2Rho0MtLhVal

      allocate( fReal(atoms%jmtd), fImag(atoms%jmtd) )
      fReal(:) = 0.0
      fImag(:) = 0.0
      do imesh = 1, atoms%jri(itype)
        fReal(imesh) = real( r2Rho0MTLhVal(imesh, 0, 1, 1) * clnu_atom(1, 0, iatom) )
        fImag(imesh) = aimag( r2Rho0MTLhVal(imesh, 0, 1, 1) * clnu_atom(1, 0, iatom) )
      end do ! imesh
      call intgr3LinIntp(fReal, atoms%rmsh(1,itype), atoms%dx(itype), atoms%jri(itype), intgrReal, 1)
      call intgr3LinIntp(fImag, atoms%rmsh(1,itype), atoms%dx(itype), atoms%jri(itype), intgrImag, 1)

      if (abs( 0.5 * sqrt(fpi) * cmplx(intgrReal, intgrImag) - PsiPsi(iatom) ) > 9e-7) passed = .false.

      ! Calculate the gradient of the unperturbed density.
      !   Any factors as in Fleur have been removed. In order to improve the numerical accuracy of the gradient r^2 is multiplied
      !   before and removed after the calculation of the gradient. Thereyby, subtraction of small numbers close to the core are
      !   avoided.
      call CalcGrR2FinLH( atoms, lathar, clnu_atom, nmem_atom, mlh_atom, r2Rho0MTlhVal(:, :, :, 1), r2GrRho0MTVal )

      ! Compare gradient of valence density to output of calcvecmatelem
      do idir = 1, 3
        fReal(:) = 0.0
        fImag(:) = 0.0
        do imesh = 1, atoms%jri(itype)
          fReal(imesh) = real( r2GrRho0MTVal(imesh, 1, iatom, idir) )
          fImag(imesh) = aimag( r2GrRho0MTVal(imesh, 1, iatom, idir) )
        end do ! imesh
        call intgr3LinIntp(fReal, atoms%rmsh(1,itype), atoms%dx(itype), atoms%jri(itype), intgrReal, 1)
        call intgr3LinIntp(fImag, atoms%rmsh(1,itype), atoms%dx(itype), atoms%jri(itype), intgrImag, 1)

        if ( abs( 0.25 * sqrt(fpi) * cmplx(intgrReal, intgrImag) - grPsiPsi(idir, iatom) ) > 9e-7 ) passed = .false.

        ! Compare lm-resolved integrals to the integrals of grRho lm-resolved and also the integrands which are additionaly mesh
        ! resolved to the non-integrated density gradient or density tensor geadients
        if (simpGradient) then
          oqn_l = 1
          do mqn_m = -oqn_l, oqn_l
            lm = oqn_l * (oqn_l + 1) + 1 + mqn_m
            if (beforeIntegration) then
              do imesh = 1, atoms%jri(itype)
                write(4000, '(3i8, 2f15.8)') idir, mqn_m, imesh, grPsiPsiMesh(imesh, lm, idir)
                write(4001, '(3i8, 2f15.8)') idir, mqn_m, imesh, 0.25 * r2GrRho0MTVal(imesh, lm, iatom, idir) / atoms%rmsh(imesh, itype)**2
              end do ! imesh
            end if
            if (output) then
              write(4002, '(2i8, 2f15.8)') idir, mqn_m, grPsiPsiTest(lm, idir)
            end if
            fReal(:) = 0.0
            fImag(:) = 0.0
            do imesh = 1, atoms%jri(itype)
              fReal(imesh) = 0.25 * real( r2GrRho0MtVal(imesh, lm, iatom, idir) )
              fImag(imesh) = 0.25 * aimag(r2GrRho0MtVal(imesh, lm, iatom, idir)  )
            end do ! imesh
            call intgr3LinIntp(fReal, atoms%rmsh(1,itype), atoms%dx(itype), atoms%jri(itype), intgrReal, 1)
            call intgr3LinIntp(fImag, atoms%rmsh(1,itype), atoms%dx(itype), atoms%jri(itype), intgrImag, 1)
            if (output) then
              write(4003, '(2i8, 2f15.8)') idir, mqn_m, cmplx(intgrReal, intgrImag)
            end if
            if (abs(cmplx(intgrReal, intgrImag) - grPsiPsiTest(lm, idir)) > 9e-7) passed = .false.
          end do ! mqn_m
        end if

      end do ! idir

      if (tensorGradient) then
        do idirC = 1, 3

          r2GrGrtRho0MTVal(:, :, :, :) = cmplx(0., 0.)
          call CalcGrR2FinSH( atoms, r2GrRho0MTVal(:, :, :, idirC), r2GrGrtRho0MTVal )

          if (beforeIntegration) then
            do idirR = 1, 3
              do oqn_l = 0, 2, 2
                do mqn_m = -oqn_l, oqn_l
                  lm = oqn_l * (oqn_l + 1) + 1 + mqn_m
                  do imesh = lMesh, uMesh
                    write(4010, '(5i8,2f30.8)') idirC, idirR, oqn_l, mqn_m, imesh, grGrtPsiPsiMesh(imesh, lm, idirR, idirC)
                    write(4011, '(5i8,2f30.8)') idirC, idirR, oqn_l, mqn_m, imesh, grPsiGrtPsiMesh(imesh, lm, idirR, idirC)
                    write(4012, '(5i8,2f30.8)') idirC, idirR, oqn_l, mqn_m, imesh, grGrtPsiPsiMesh(imesh, lm, idirR, idirC) + grPsiGrtPsiMesh(imesh, lm, idirR, idirC)
                    write(4013, '(5i8,2f30.8)') idirC, idirR, oqn_l, mqn_m, imesh, 0.25 * r2GrGrtRho0MTVal(imesh, lm, iatom, idirR) / atoms%rmsh(imesh, itype)**2
                  end do ! imesh
                end do ! mqn_m
              end do ! oqn_l
            end do ! idirR
          end if

          if (afterIntegration) then
            do idirR = 1, 3
              do oqn_l = 0, 2, 2
                do mqn_m = -oqn_l, oqn_l
                  lm = oqn_l * (oqn_l + 1) + 1 + mqn_m
                  if (output) then
                    write(4014, '(4i8,2f15.8)') idirC, idirR, oqn_l, mqn_m, grPsiGrtPsiTest(lm, idirR, idirC)
                    write(4015, '(4i8,2f15.8)') idirC, idirR, oqn_l, mqn_m, grGrtPsiPsiTest(lm, idirR, idirC)
                    write(4016, '(4i8,2f15.8)') idirC, idirR, oqn_l, mqn_m, grGrtPsiPsiTest(lm, idirR, idirC) + grPsiGrtPsiTest(lm, idirR, idirC)
                  end if
                  fReal(:) = 0.
                  fImag(:) = 0.
                  do imesh = 1, atoms%jri(itype)
                    fReal(imesh) = 0.25 * real( r2GrGrtRho0MtVal(imesh, lm, iatom, idirR) )
                    fImag(imesh) = 0.25 * aimag(r2GrGrtRho0MtVal(imesh, lm, iatom, idirR)  )
                  end do ! imesh
                  call intgr3LinIntp(fReal, atoms%rmsh(1, itype), atoms%dx(itype), atoms%jri(itype), intgrReal, 1)
                  call intgr3LinIntp(fImag, atoms%rmsh(1, itype), atoms%dx(itype), atoms%jri(itype), intgrImag, 1)
                  if (output) then
                    write(4017, '(4i8,2f15.8)') idirC, idirR, oqn_l, mqn_m, cmplx(intgrReal, intgrImag)
                  end if
                  if (abs(grGrtPsiPsiTest(lm, idirR, idirC) + grPsiGrtPsiTest(lm, idirR, idirC) - cmplx(intgrReal, intgrImag)) > 9e-6 ) passed = .false.
                end do ! mqn_m
              end do ! oqn_l
            end do ! idirR
          end if ! afterIntegration
        end do ! idirC
      end if
    end do ! ieqat
  end do ! itype

  if (passed) then
    write (logUnit, '(a)')'                                                      |__ passed!'
  else
    write (logUnit, '(a)')'                                                      |__ failed!'
    call juDFT_warn('Testing <grPsi|Psi>, <grGrtPsi|Psi> and <grPsi|grTPsi> failed.', &
      & calledby='TestR2orNotWfMtGradNgrNTensGrOvls', hint='Debug.')
  end if

  ! For the simple gradient it does not make much a difference whether r2 is multiplied or not, for the tensor gradient it does make
  ! a difference for the Neon test system, and one can see that without r2 multiplied it is a bit more accurate. This can come from
  ! the fact that r2 dampens down the quantitiy to be derived so much that numerical errors occur because of subtraction of small
  ! numbers.

  end subroutine TestR2orNotWfMtGradNgrNTensGrOvls

  ! Calculate the integral of gradRho gradVeff0 and compare it to the matrix element <grPsi|grVeff|Psi> in the IR.
  subroutine testIRIntegral(atoms, cell, stars, Veff0, dimens, kpts, input, results, ngdp, gdp, nv, mapGbas, gBasVec, nobd, z0, kpq2kPrVec, logUnit)

    use m_types
    use m_jpConstants, only : iu
    use m_jpPotDensHelper, only : warpIRPot, convertStar2G
    use mod_juPhonUtils, only : Calc2ArgIntIR
    use m_jpSternhHF, only : CalcMEPotIR
    use m_juDFT_NOstopNO, only : juDFT_warn

    implicit none

    type(t_atoms),     intent(in) :: atoms
    type(t_cell),      intent(in) :: cell
    type(t_stars),     intent(in) :: stars
    type(t_potential), intent(in) :: Veff0
    type(t_dimension), intent(in) :: dimens
    type(t_kpts),      intent(in) :: kpts
    type(t_input),     intent(in) :: input
    type(t_results),   intent(in) :: results

    integer,           intent(in) :: ngdp
    integer,           intent(in) :: logUnit


    integer,           intent(in) :: gdp(:, :)
    integer,           intent(in) :: nv(:, :)
    integer,           intent(in) :: mapGbas(:, :, :)
    integer,           intent(in) :: gBasVec(:, :)
    integer,           intent(in) :: nobd(:, :)
    complex,           intent(in) :: z0(:, :, :, :)
    integer,           intent(in) :: kpq2kPrVec(:, :, :)

    integer                       :: idirC
    integer                       :: idirR
    integer                       :: iG
    integer                       :: ikpt
    integer                       :: nmat
    integer                       :: iBas
    integer                       :: iband
    logical                       :: passed = .true.
    logical                       :: passedWarping

    complex                       :: integralIR(3, 3)


    complex, allocatable          :: grVeff0IR(:, :)
    complex, allocatable          :: w_grVeff0IR(:, :)
    complex, allocatable          :: rho1IRContainer(:, :, :)
    complex, allocatable          :: rho1MTContainer(:, :, :, :, :)
    complex, allocatable          :: vEff1IRContainer(:, :, :)
    complex, allocatable          :: w_vEff1IRContainer(:, :, :)
    complex, allocatable          :: vEff1MTContainer(:, :, :, :, :)
    complex, allocatable          :: vEff0Pw(:)
    complex, allocatable          :: rho0IRpw(:)
    complex, allocatable          :: grRho0IR(:, :)
    complex, allocatable          :: testIR(:, :)
    complex, allocatable          :: z1nGDummy(:, :, :)
    complex ,allocatable          :: rho0IRnCt(:, :)
    complex                       :: testIRSum(3, 3)
    real                          :: Gext(3)
    real                          :: kExt(3)
    real                          :: gBasExt(3)
    ! code

    write(logUnit, '(a)')   'Compare int dV grRho_IR grVeff_IR with <grPsi|grVeff|Psi>_IR'
    write(logUnit, '(a)')   '------------------------------------------------------------'

    allocate(rho0IRncT(stars%n3d, 1))
    rho0IRncT(:, :) = cmplx(0.,0.)
    rewind(7800)
    read(7800) rho0IRncT
    ! Convert from star expansion coefficients to plane-wave expansion coefficients
    allocate( vEff0Pw(ngdp), rho0IRpw(ngdp) )
    vEff0Pw(:) = cmplx(0., 0.)
    rho0IRpw(:) = cmplx(0., 0.)
    call convertStar2G( Veff0%vpw_uw(:, 1), vEff0Pw(:), stars, ngdp, gdp )
    call convertStar2G( rho0IRncT(:, 1), rho0IRpw(:), stars, ngdp, gdp )

    ! Perform analytical gradient of unperturbed density and effective potential in the interstitial region
    allocate( grVeff0IR(ngdp, 3), grRho0IR(ngdp, 3) )
    grVeff0IR(:, :) = cmplx(0., 0.)
    grRho0IR(:, :) = cmplx(0., 0.)
    do iG = 1, ngdp
      Gext(1:3) = matmul(cell%bmat(1:3, 1:3), gdp(1:3, iG))
      grVeff0IR(iG, 1:3) = iu * Gext(1:3) * vEff0Pw(iG)
      grRho0IR(iG, 1:3)  = iu * Gext(1:3) * rho0IRpw(iG)
    end do ! iG

    ! Warp interstitial second-order external potential
    allocate(w_grVeff0IR(ngdp, 3))
    w_grVeff0IR = cmplx(0.0, 0.0)
    do idirC = 1, 3
      call warpIRPot(stars, ngdp, idirC, gdp, grVeff0IR, w_grVeff0IR(:, idirC))
    end do ! idirC

    ! Fill container arrays with gradient coefficients so that they can be processed by EvelIntRho1Veff1
    allocate( rho1IRContainer(ngdp, 3, atoms%nat), vEff1IRContainer(ngdp, 3, atoms%nat), w_vEff1IRContainer(ngdp, 3, atoms%nat), &
      & rho1MTContainer( atoms%jmtd, (atoms%lmaxd + 1)**2, atoms%nat, 3, atoms%nat ), &
      & vEff1MTContainer( atoms%jmtd, (atoms%lmaxd + 1)**2, atoms%nat, 3, atoms%nat ) )

    rho1IRContainer(:, :, :) = cmplx(0., 0.)
    w_vEff1IRContainer(:, :, :) = cmplx(0., 0.)
    vEff1IRContainer(:, :, :) = cmplx(0., 0.)
    rho1MTContainer(:, :, :, :, :) = cmplx(0., 0.)
    vEff1MTContainer(:, :, :, :, :) = cmplx(0., 0.)

    do idirC = 1, 3
      do iG = 1, ngdp
        rho1IRContainer(iG, idirC, 1) = -grRho0IR(iG, idirC)
        w_vEff1IRContainer(iG, idirC, 1) = -w_grVeff0IR(iG, idirC)
        vEff1IRContainer(iG, idirC, 1) = -grVeff0IR(iG, idirC)
      end do ! iG
    end do ! idirC

    integralIR(:, :) = cmplx(0., 0.)
    do idirC = 1, 3
      do idirR = 1, 3
        call Calc2ArgIntIR(cell, ngdp, rho1IRContainer(:, idirR, 1), w_vEff1IRContainer(:, idirC, 1), integralIR(idirR, idirC))
      end do ! idirR
    end do ! idirC

    call testWarpingHFdynMat( atoms, cell, stars, ngdp, 1, gdp, rho1IRContainer, vEff1IRContainer, w_vEff1IRContainer, passedWarping )
    passed = passed .and. passedWarping

    if (.false.) then
      write(*, '(a)') 'Integral x 0.25'
      write(*, '(3(2(es16.8,1x),3x))') 0.25 * integralIR(1, :)
      write(*, '(3(2(es16.8,1x),3x))') 0.25 * integralIR(2, :)
      write(*, '(3(2(es16.8,1x),3x))') 0.25 * integralIR(3, :)
    end if

    allocate(testIR(4, 4))
    allocate(z1nGDummy(dimens%nbasfcn, dimens%neigd, 3 ))

    testIRSum(:, :) = cmplx(0., 0.)
    testIR(:, :) = cmplx(0., 0.)
    z1nGDummy(:, :, :) = cmplx(0., 0.)

    do ikpt = 1, kpts%nkpt
      nmat = nv(1, ikpt) + atoms%nlotot
      kExt(1:3) = matmul(cell%bmat, kpts%bk(1:3, ikpt))
      do iBas = 1, nv(1, ikpt) !+ atoms%nlotot
        gbasExt(1:3) = matmul( cell%bmat(1:3, 1:3), gBasVec(1:3, mapGbas(iBas, ikpt, 1)))
        do iband = 1, nobd(ikpt, 1)
          do idirC = 1, 3
            z1nGDummy(iBas, iband, idirC) = iu * ( kExt(idirC) + gbasExt(idirC) ) * z0(iBas, iband, ikpt, 1)
          end do ! idirC
        end do ! iband
      end do ! iBas

      do idirC = 1, 3
        do idirR = 1, 3
          testIR(:, :) = cmplx(0., 0.)
          call calcMEPotIR( stars, dimens, GbasVec(:, mapGbas(:nv(1, ikpt), ikpt, 1)), &
            & GbasVec(:, mapGbas(:nv(1, ikpt), ikpt, 1)), nv, grVeff0IR(:, idirC), z1nGDummy(:, :, idirR), z0(:, :, ikpt, 1), gdp, nmat,&
            & nobd(ikpt, 1), nobd(ikpt, 1), ikpt, 1, ikpt, ngdp, testIR, kpq2kPrVec ) !todo spin-relevant
          do iband = 1, nobd(ikpt, 1)
            testIRSum(idirR, idirC) = testIRSum(idirR, idirC) + results%w_iks(iband, ikpt, 1) * testIR(iband, iband)
          end do ! iband
        end do ! idirR
      end do ! idirC
    end do ! ikpt

    if (.false.) then
      write(*, '(a)') 'testIRSum'
      write(*, '(3(2(es16.8,1x),3x))') testIRSum(1, :)
      write(*, '(3(2(es16.8,1x),3x))') testIRSum(2, :)
      write(*, '(3(2(es16.8,1x),3x))') testIRSum(3, :)
    end if

    if (any( abs(testIRSum(:, :) - 0.25 * integralIR(:, :)) > 9e-9) ) passed = .false.

    if ( passed ) then
      write( logUnit, '(a)' ) '                                                           |__ passed!'
    else
      write( logUnit, '(a)' ) '                                                           |__ failed!'
      call JuDFT_warn('Compare <grPsi|grVeff|Psi>_IR with int dV grRho_IR grVeff_IR', calledby='testIRIntegral', hint='Check logfile for more information!')
    end if

  end subroutine testIRIntegral

  ! test the cutoffs of the warping function by performing the warping via the reciprocal representation of the step function
  subroutine testWarpingHFdynMat( atoms, cell, stars, ngdp, iDatom, gdp, rho1IR, vExt1IR, w_vExt1IR, passed )

    use m_types, only : t_atoms, t_cell, t_stars
    use m_jpConstants, only : iu, fpi, tpi
    use m_sphbes
    use m_jpPotDensHelper, only : convertStar2G

    implicit none

    type(t_atoms),              intent(in)  :: atoms
    type(t_cell),               intent(in)  :: cell
    type(t_stars),              intent(in)  :: stars

    integer,                    intent(in)  :: ngdp
    integer,                    intent(in)  :: iDatom

    integer,                    intent(in)  :: gdp(:, :)
    complex,                    intent(in)  :: rho1IR(:, :, :)
    complex,                    intent(in)  :: vExt1IR(:, :, :)
    complex,                    intent(in)  :: w_vExt1IR(:, :, :)
    logical,                    intent(out) :: passed

    real,          allocatable              :: sbes(:)
    integer                                 :: idirR
    integer                                 :: idirC
    real                                    :: pref
    integer                                 :: iG
    integer                                 :: iGp
    real                                    :: GdiffCartAbsV
    complex                                 :: stepFunc
    integer                                 :: itype
    integer                                 :: iatom
    integer                                 :: ieqat
    complex                                 :: dynMatHfFFT(3, 3)
    complex                                 :: dynMatHfSUM(3, 3)
    real                                    :: Gext(3)
    real                                    :: Gdiff(3)
    real                                    :: Gpext(3)
    real                                    :: GdiffCart(3)

    ! If we want to debug stars%ustep in stepf.F, we have to comment out the read in of wkf2 containing the stepfunction


    ! Integrate two vector quantities in the interstitial so that the Bloch character vanishes by warping one of them and
    ! calculating the dot_product
    dynMatHfFFT(:, :) = cmplx(0., 0.)
    do idirC =  1, 3
      do idirR =  1, 3
        dynMatHfFFT(idirR, idirC) = cell%omtil * dot_product( rho1IR(:, idirR, iDatom), w_vExt1IR(:, idirC, iDatom) )
      end do ! idirR
    end do ! idirC

    ! Integrate the same quantities by double-sum a product of these quantities and the step function in reciprocal space
    allocate(sbes(0:atoms%lmax(atoms%ntype)))

    sbes = cmplx(0., 0.)
    dynMatHfSUM(:, :) = cmplx(0., 0.)
    pref = fpi * atoms%rmt(1)**2
    do iG = 1, ngdp
      Gext(:) = matmul(cell%bmat, gdp(:, iG))
      do iGp = 1, ngdp
        Gdiff(1:3) = gdp(1:3, iG) - gdp(1:3, iGp)
        Gpext(1:3) = matmul(cell%bmat, gdp(1:3, iGp))
        GdiffCart(1:3) = Gext(1:3) - Gpext(1:3)
        GdiffCartAbsV = norm2(GdiffCart)
        stepFunc = cmplx(0., 0.)
        if (iG == iGp) then
          ! It's just the interstitial volume
          stepFunc = cell%omtil
          iatom = 0
          do itype = 1, atoms%ntype
            do ieqat = 1, atoms%neq(itype)
              iatom = iatom + 1
              stepFunc = stepFunc - atoms%volmts(itype)
            end do ! ieqat
          end do ! itype
        else
          iatom = 0
          do itype = 1, atoms%ntype
            call sphbes(atoms%lmax(itype), GdiffCartAbsV * atoms%rmt(itype), sbes)
            do ieqat = 1, atoms%neq(itype)
              iatom = iatom + 1
              stepFunc =  -pref * exp(-tpi * iu * dot_product(gDiff(1:3), atoms%taual(1:3, itype))) * sbes(1) / GdiffCartAbsV
            end do ! ieqat
          end do ! itype
        end if
        do idirC = 1, 3
          do idirR = 1, 3
            dynMatHfSUM(idirR, idirC) = dynMatHfSUM(idirR, idirC) + conjg(rho1IR(iG, idirR, iDatom)) * vExt1IR(iGp, idirC, iDatom) &
                                                                                                                        & * stepFunc
          end do ! idirR
        end do ! idirC
      end do ! iGp
    end do ! iG

    if (maxval(abs(dynMatHfSUM - dynMatHfFFT)) > 10e-10) then
      passed = .false.
    else
      passed = .true.
    end if
    if (.false.) then
      write(*, *) 'Test of warping'
      write(*, *) 'Max. deviation:', maxval(abs(dynMatHfSUM - dynMatHfFFT))
    end if

  end subroutine testWarpingHFdynMat

  ! This routine adds <Psi|H - eps|Psi> IR with MT and at the bottom writes out IR, MT and the sum of them, which have to vanish.
  subroutine testIR3rdMatElem( atoms, stars, lathar, dimens, usdus, sym, kpts, qpts, cell, input, logUnit, mlh_atom, clnu_atom, vEff0MtLh, rbas1, rbas2,&
      & ilo2p, nmem_atom, nRadFun, gbas, mapGbas, nv, kveclo, nobd, z, El, iloTable, eig, kpq2kPrVec, w_vEff0IR, gdp, ngdp, vEff0IR )

#include"cppmacro.h"

    use m_types, only : t_atoms, t_stars, t_sphhar, t_sym, t_cell, t_kpts, t_dimension, t_usdus, t_enpara, t_input, t_tlmplm
     
    use m_abcof3
    use m_juDFT_NOstopNO, only : juDFT_warn
    use m_jpSternhPulaySurface, only : tlmplm4H0, calcHS0MT
    use m_jpSetupDynMat, only : CalcScalBasfMatElems, CalcPsi1HepsPsi1IR
    use m_jpConstants, only : iu
    use m_jpSternhHF, only : calcMEPotIR
    use m_jpPotDensHelper, only : WarpIRPot, convertStar2G
    use m_jpSetupDynMatHelper, only : CalcHnGrV0Varphi
!
    implicit none

    ! Type parameter
    type(t_atoms),                  intent(in) :: atoms
    type(t_stars),                  intent(in) :: stars
    type(t_sphhar),                 intent(in) :: lathar
    type(t_dimension),              intent(in) :: dimens
    type(t_usdus),                  intent(in) :: usdus
    type(t_sym),                    intent(in) :: sym
    type(t_kpts),                   intent(in) :: kpts
    type(t_kpts),                   intent(in) :: qpts
    type(t_cell),                   intent(in) :: cell
    type(t_input),                  intent(in) :: input

    ! Scalar parameter
    integer,                        intent(in) :: logUnit

    ! Array parameter
    integer,                        intent(in) :: mlh_atom(:, 0:, :)
    complex,                        intent(in) :: clnu_atom(:, 0:, :)
    real,                           intent(in) :: vEff0MtLh(:, 0:, :)
    real,                           intent(in) :: rbas1(:,:,0:, :, :)
    real,                           intent(in) :: rbas2(:,:,0:, :, :)
    integer,                        intent(in) :: ilo2p(:, :)
    integer,                        intent(in) :: nmem_atom(0:, :)
    integer,                        intent(in) :: nRadFun(0:, :)
    integer,                        intent(in) :: gbas(:, :)
    integer,                        intent(in) :: mapGbas(:, :, :)
    integer,                        intent(in) :: nv(:, :)
    integer,                        intent(in) :: kveclo(:,:)
    integer,                        intent(in) :: nobd(:, :)
    MCOMPLEX,                       intent(in) :: z(:,:,:,:)
    real,                           intent(in) :: El(:, 0:, :, :)
    integer,                        intent(in) :: iloTable(:, 0:, :)
    real,                           intent(in) :: eig(:, :, :)
    integer,                        intent(in) :: kpq2kPrVec(:, :, :)
    complex,                        intent(in) :: w_vEff0IR(:,:)
    complex,                        intent(in) :: vEff0IR(:,:)
    integer,                        intent(in) :: gdp(:, :)

    ! Type variable
    type(od_inp)                               :: odi
    type(od_sym)                               :: ods

    ! Scalar variables
    integer                                    :: iatom
    integer                                    :: itype
    integer                                    :: ieqat
    integer                                    :: ptsym
    integer                                    :: ilh
    integer                                    :: oqn_l
    integer                                    :: lm_pre
    integer                                    :: imem
    integer                                    :: mqn_m
    integer                                    :: lm
    integer                                    :: imesh
    integer                                    :: lmpMax
    integer                                    :: nRadFunMax
    integer                                    :: passed
    integer                                    :: ikpt
    integer                                    :: nmat
    integer                                    :: iband
    integer                                    :: iDtype
    integer                                    :: iDatom
    integer                                    :: iDeqat
    integer                                    :: lmp
    integer                                    :: pMaxLocal
    integer                                    :: iradf
    complex                                    :: matEntryRef
    complex                                    :: wholeUnitCellIR
    integer                                    :: ngdp
    complex                                    :: overlapIRkn
    complex                                    :: overlapRevIRkn
    integer                                    :: iG

    ! Array variabels
    integer,           allocatable             :: ngoprI(:)
    real,              allocatable             :: varphi1(:, :, :)
    real,              allocatable             :: varphi2(:, :, :)
    complex,           allocatable             :: vEff0MtSpH(:, :, :)
    complex,           allocatable             :: hVarphi(:, :, :, :)
    integer,           allocatable             :: lmpT(:)
    real,              allocatable             :: r2(:)
    complex,           allocatable             :: a(:, :, :)
    complex,           allocatable             :: b(:, :, :)
    complex,           allocatable             :: bascof_lo(:, :, :, :, :)
    complex,           allocatable             :: ab0cofScl(:, :)
    complex,           allocatable             :: varphiPsi(:)
    complex,           allocatable             :: psiPsi(:, :)
    complex,           allocatable             :: varphiHPsi(:)
    complex,           allocatable             :: psiHpsi(:, :)
    integer,           allocatable             :: nlo_atom(:)
    real,          allocatable              :: varphiVarphi(:, :, :)
    complex,       allocatable              :: varphiHvarphi(:, :, :)
    complex,       allocatable              :: vEff0NsphGrVarphi(:, :, :, :, :)
    complex,       allocatable              :: r2grVeff0SphVarphi(:, :, :, :, :)
    complex,       allocatable              :: r2grVeff0SphVarphiDummy(:, :, :, :, :)
    integer,       allocatable              :: grVarphiChLout(:, :)
    integer,       allocatable              :: grVarphiChMout(:, :)
    real,          allocatable              :: grVarphiCh1(:, :, :, :)
    real,          allocatable              :: grVarphiCh2(:, :, :, :)
    complex,       allocatable             :: psiHepsPsiMT(:, :)
    complex,       allocatable             :: wholeUnitCellIRkSum(:, :)
    complex,       allocatable             :: IRtest(:, :)
    complex,       allocatable            :: IRtestSummed(:)
    complex,       allocatable            :: w_vEff0IRpw(:)
    complex,       allocatable            :: vEff0IRpw(:)
    complex,       allocatable            :: z0Dummy(:, :)
    complex,       allocatable            :: w_z0Dummy(:)
    real                                  :: kpGExt(1:3)
    write( logUnit, '(a)' ) 'PsiHepsPsi test'
    write( logUnit, '(a)' ) '---------------'

    ! We do not want the local coordinate systems to be rotated for non-representative atoms constructing the matching coefficients.
    allocate(ngoprI(atoms%nat))
    allocate( lmpT(atoms%ntype) )
    lmpT(:) = 0
    do itype = 1, atoms%ntype
      lmpT(itype) = sum( [ ( (2 * oqn_l + 1)* nRadFun(oqn_l, itype), oqn_l = 0, atoms%lmax(itype) ) ] )
    end do ! itype
    lmpMax     = maxval( lmpT(:) )
    nRadFunMax = maxval( nRadFun(:, :) )
!    ! The a are dimensioned until lmax + 2 for the tested routines to work, but we only need the abcof until lmax for this test
    allocate( varphi1(atoms%jmtd, nRadFunMax, 0:atoms%lmaxd), varphi2(atoms%jmtd, nRadFunMax, 0:atoms%lmaxd) )
    allocate( hVarphi(2, atoms%jmtd, 0:(atoms%lmaxd + 3)**2 - 1, lmpMax) )
    allocate( r2(atoms%jmtd) )
    allocate( a( dimens%nvd, 0:(atoms%lmaxd + 3)**2 - 1, atoms%nat), b(dimens%nvd, 0:(atoms%lmaxd + 3)**2 - 1, atoms%nat), &
      & bascof_lo(3, -atoms%llod:atoms%llod, 4 * atoms%llod + 2, atoms%nlod, atoms%nat) )
    allocate( varphiVarphi(lmpMax, lmpMax, atoms%ntype), varphiHvarphi(lmpMax, lmpMax, atoms%nat))
    allocate( r2grVeff0SphVarphi(2, atoms%jmtd, 0:(atoms%lmaxd + 3)**2 - 1, lmpMax, 3) )
    allocate( r2grVeff0SphVarphiDummy(2, atoms%jmtd, 0:(atoms%lmaxd + 3)**2 - 1, lmpMax, 3) )
    allocate( vEff0NsphGrVarphi(2, atoms%jmtd, 0:(atoms%lmaxd + 3)**2 - 1, lmpMax, -1:1) )
    allocate( grVarphiCh1(atoms%jmtd, 2, lmpMax, -1:1), grVarphiCh2(atoms%jmtd, 2, lmpMax, -1:1), &
            & grVarphiChLout(2, 0:atoms%lmaxd), grVarphiChMout(-atoms%lmaxd:atoms%lmaxd, -1:1) )
    allocate(ab0cofScl(maxval(nobd), lmpMax))
    allocate(varphiPsi(lmpMax))
    allocate(psiPsi(atoms%nat, maxval(nobd)))
    allocate(varphiHPsi(lmpMax))
    allocate(psiHpsi(atoms%nat, maxval(nobd)))
    allocate( psiHepsPsiMT(atoms%nat, maxval(nobd)))
    allocate( wholeUnitCellIRkSum(atoms%nat, maxval((nobd))))
        allocate(nlo_atom(atoms%nat))
        nlo_atom = 0

    allocate( z0Dummy(dimens%nvd, 1))
    z0Dummy = cmplx(0.)

    allocate(IRtest(maxval(nobd), maxval(nobd)))
    allocate(IRtestSummed(maxval(nobd)))
    IRtestSummed = cmplx(0., 0.)
    allocate( w_vEff0IRpw(ngdp) )
    allocate( vEff0IRpw(ngdp) )
    w_vEff0IRpw(:) = cmplx(0., 0.)
    vEff0IRpw(:) = cmplx(0., 0.)
    allocate(w_z0Dummy(dimens%nvd))
    w_z0Dummy(:) = cmplx(0., 0.)

    varphiVarphi(:, :, :)            = 0.
    varphiHvarphi(:, :, :)              = cmplx(0., 0.)
    r2grVeff0SphVarphi(:, :, :, :, :) = cmplx(0., 0.)
    r2grVeff0SphVarphiDummy(:, :, :, :, :) = cmplx(0., 0.)
    vEff0NsphGrVarphi(:, :, :, :, :) = cmplx(0., 0.)
    grVarphiCh1(:, :, :, :) = 0.
    grVarphiCh2(:, :, :, :) = 0.
    grVarphiChLout(:, :) = 0
    grVarphiChMout(:, :) = 0
    psiHepsPsiMT(:, :) = 0.
    wholeUnitCellIRkSum(:, :) = cmplx(0., 0.)

    ! Calculate k-independent quantities

    ! Transform unperturbed potential given in lattice harmonics into spherical harmonic representation
    allocate(vEff0MtSpH(atoms%jmtd, (atoms%lmaxd + 1)**2, atoms%nat))
    vEff0MTSpH = cmplx(0.0, 0.0)
    iatom = 0
    do itype = 1, atoms%ntype
      do ieqat = 1, atoms%neq(itype)
        iatom = iatom + 1
        ptsym = atoms%ntypsy(iatom)
        do ilh = 0, lathar%nlh(ptsym)
          oqn_l = lathar%llh(ilh, ptsym)
          lm_pre = oqn_l * (oqn_l + 1) + 1
          do imem = 1, nmem_atom(ilh, iatom)
            mqn_m = mlh_atom(imem, ilh, iatom)
            lm = lm_pre + mqn_m
            do imesh = 1, atoms%jri(itype)
              vEff0MtSpH(imesh, lm, iatom) = vEff0MtSpH(imesh, lm, iatom) + vEff0MtLh(imesh, ilh, itype) * clnu_atom(imem, ilh, iatom)
            end do ! imesh
          end do ! imem
        end do ! ilh
      end do ! ieqat
    end do ! itype

    iatom = 0
    do itype = 1, atoms%ntype
      r2(:) = 0.

      ! Precalculate radial Jacobi determinant for later integrals
      do imesh = 1, atoms%jri(itype)
        r2(imesh) = atoms%rmsh(imesh, itype) * atoms%rmsh(imesh, itype)
      end do ! imesh
      varphi1(:, :, :) = 0.
      varphi2(:, :, :) = 0.
      do oqn_l = 0, atoms%lmax(itype)
        do iradf = 1, nRadFun(oqn_l, itype)
          do imesh = 1, atoms%jri(itype)
            ! In Fleur the radial solutions u_lp are multiplied by a factor r per default to avoid an additional multiplication of the
            ! Jacobi determinant r^2 in radial integrals given spherical coordinates.
            varphi1(imesh, iradf, oqn_l) = rbas1(imesh, iradf, oqn_l, itype, 1) / atoms%rmsh(imesh, itype)
            varphi2(imesh, iradf, oqn_l) = rbas2(imesh, iradf, oqn_l, itype, 1) / atoms%rmsh(imesh, itype)
          end do ! imesh
        end do ! iradf
      end do ! oqn_l

      do ieqat = 1, atoms%neq(itype)
        iatom = iatom + 1
        hVarphi = cmplx(0.0, 0.0)
        ! Calculate the action of the hamiltonian onto an unperturbed basis function hVarphi
        call CalcHnGrV0Varphi( atoms, lathar, itype, iatom, lmpMax, El, varphi1, varphi2, nRadFun, vEff0MtSpH(:, :, iatom), vEff0MtLh, clnu_atom, &
          & nmem_atom, mlh_atom, grVarphiCh1, grVarphiCh2, grVarphiChLout, grVarphiChMout, hVarphi, vEff0NsphGrVarphi, r2grVeff0SphVarphi, r2grVeff0SphVarphiDummy )

        ! Calculate the matrix element <varphi|varphi> and <varphi|H|varphi>
        call CalcScalBasfMatElems( atoms, itype, iatom, nRadFun, r2, varphi1, varphi2, hVarphi, varphiVarphi, varphiHvarphi )
      end do ! ieqat
    end do ! itype
    deallocate(r2grVeff0SphVarphiDummy)

    ! How dangerous is interpolation at first mesh points?
!    lmp = 0
!    do oqn_l = 0, atoms%lmax(1)
!      do mqn_m = -oqn_l, oqn_l
!
!        !p = 1
!        lmp = lmp + 1
!        write(2292, '(2f15.8)') varphiVarphi(lmp, lmp, 1), 1.
!
!        ! p = 2
!        lmp = lmp + 1
!        write(2292, '(2f15.8)') varphiVarphi(lmp, lmp, 1), usdus%ddn(oqn_l, 1, 1)
!
!      end do ! mqn_m
!    end do ! oqn_l
!    NOstopNO

    ngoprI(:) = 1
    passed = .true.
    call convertStar2G(w_vEff0IR(:, 1), w_vEff0IRpw, stars, ngdp, gdp)
    call convertStar2G(vEff0IR(:, 1), vEff0IRpw, stars, ngdp, gdp)
    do ikpt = 1, kpts%nkpt

      ! Small matching coefficients
      nmat = nv(1, ikpt) + atoms%nlotot
      a(:, :, :) = cmplx(0.0, 0.0)
      b(:, :, :) = cmplx(0.0, 0.0)
      bascof_lo(:, :, :, :, :) = cmplx(0.0, 0.0)
      call abcof3( atoms%lmaxd, atoms%ntype, atoms%nat, sym%nop, dimens%nvd, dimens%jspd, 1, dimens%lmd, dimens%nbasfcn, &
        & atoms%llod, atoms%nlod, atoms%nlotot, sym%invtab, atoms%ntype, sym%mrot, ngoprI, atoms%taual, atoms%neq, atoms%lmax, &
        & atoms%rmt, cell%omtil, cell%bmat, cell%bbmat, kpts%bk(:, ikpt), gbas(1, mapGbas(:nv(1, ikpt), ikpt, 1)), &
        & gbas(2, mapGbas(:nv(1, ikpt), ikpt, 1)), gbas(3, mapGbas(:nv(1, ikpt), ikpt, 1)), nv(:, ikpt), nmat, &
        & usdus%us, usdus%dus, usdus%uds, usdus%duds, usdus%ddn, atoms%invsat, sym%invsatnr, usdus%ulos, usdus%uulon, usdus%dulon, &
        & usdus%dulos, atoms%llo, atoms%nlo, atoms%l_dulo, atoms%lapw_l, kveclo(:,ikpt), odi, ods, a, b, bascof_lo )

      psiHpsi(:, :) = cmplx(0., 0.)
      psiPsi(:, :) = cmplx(0., 0.)
      ab0cofScl(:, :) = cmplx(0., 0.)
      do iband = 1, nobd(ikpt, 1)
        iDatom = 0
        do iDtype = 1, atoms%ntype
          do iDeqat = 1, atoms%neq(iDtype)
            iDatom = iDatom + 1
            lmp = 0
            lm  = 0
            ! Large matching coefficients
            do oqn_l = 0, atoms%lmax(iDtype)
              do mqn_m = - oqn_l, oqn_l
                pMaxLocal = nRadFun(oqn_l, iDtype)
                ! p = 1
                ab0cofScl(iband, lmp + 1) = iu**oqn_l * dot_product( conjg(z(:nv(1, ikpt), iband, ikpt, 1)), a(:nv(1, ikpt), lm, iDatom) )
                ! p = 2
                ab0cofScl(iband, lmp + 2) = iu**oqn_l * dot_product( conjg(z(:nv(1, ikpt), iband, ikpt, 1)), b(:nv(1, ikpt), lm, iDatom) )
                ! Add LO contributions
                do iradf = 3, pMaxLocal
                  ! p = 1
                  ab0cofScl(iband, lmp + 1) = ab0cofScl(iband, lmp + 1) + &
                    & iu**oqn_l * dot_product( conjg(z(nv(1, ikpt) + 1:nv(1, ikpt) + atoms%nlotot, iband, ikpt, 1)), &
                                                    & bascof_lo(1, mqn_m, :atoms%nlotot, iloTable(iradf, oqn_l, iDtype), iDatom) )
                  ! p = 2
                  ab0cofScl(iband, lmp + 2) = ab0cofScl(iband, lmp + 2) + &
                    & iu**oqn_l * dot_product( conjg(z(nv(1, ikpt) + 1:nv(1, ikpt) + atoms%nlotot, iband, ikpt, 1)), &
                                                    & bascof_lo(2, mqn_m, :atoms%nlotot, iloTable(iradf, oqn_l, iDtype), iDatom) )
                  ! 2 < p < LOs for that l and that atom type
                  ab0cofScl(iband, lmp + iradf) = iu**oqn_l * dot_product( conjg(z(nv(1, ikpt) + 1:nv(1, ikpt) + atoms%nlotot, iband, ikpt, 1)), &
                                                    & bascof_lo(3, mqn_m, :atoms%nlotot, iloTable(iradf, oqn_l, iDtype), iDatom) )
                end do ! iradf
                lm = lm + 1
                lmp = lmp + pMaxLocal
              end do ! mqn_m
            end do !oqn_l

            ! Calculation of MT overlap
            varphiPsi(:) = cmplx(0., 0.)
            varphiPsi(1: lmpT(iDtype)) = &
                               & matmul( varphiVarphi(1:lmpT(iDtype), 1:lmpT(iDtype), iDtype), ab0cofScl(iband, 1:lmpT(iDtype)) )
            psiPsi(iDatom, iband) = dot_product(ab0cofScl(iband, :lmpMax), varphiPsi(:lmpMax))

            ! Setup of Hamilton matrix element with wavefunctions in bra and ket
            varphiHPsi(:) = cmplx(0.0, 0.0)
            varphiHPsi(1:lmpT(iDtype)) = matmul( varphiHvarphi(1:lmpT(iDtype), 1:lmpT(iDtype), iDatom), ab0cofScl(iband,1:lmpT(iDtype)) )
            psiHpsi(iDatom, iband) = dot_product(ab0cofScl(iband, :lmpMax), varphiHPsi(:))

            ! New
            psiHepsPsiMT(iDatom, iband) = psiHepsPsiMT(iDatom, iband) + psiHpsi(iDatom, iband) - eig(iband, ikpt, 1) * psiPsi(iDatom, iband)
!            psiHepsPsiMT(iDatom, iband) = psiHepsPsiMT(iDatom, iband) + psiPsi(iDatom, iband)

            ! vEff0IR is warped!!!!
            wholeUnitCellIR = cmplx(0., 0.)
            call calcPsi1HepsPsi1IR( kpts, qpts, stars, dimens, cell, gBas(:, mapGbas(:nv(1, ikpt), ikpt, 1)), nv, &
              & z(:nv(1, ikpt), iband, ikpt, 1), z(:nv(1, ikpt), iband, ikpt, 1), nmat, ikpt, 1, ikpt, kpq2kPrVec, &
              & w_vEff0IR, eig, iband, wholeUnitCellIR )

          !  write(3002, '(2i8,4f15.8)') iband, ikpt, wholeUnitCellIR

            wholeUnitCellIRkSum(iDatom, iband) = wholeUnitCellIRkSum(iDatom, iband) + wholeUnitCellIR

!            write(*, '(2i8,6(es15.8,2x))') ikpt, iband, wholeUnitCellIR, psiHepsPsiMT(iDatom, iband), wholeUnitCellIR + psiHepsPsiMT(iDatom, iband)

          end do ! iDeqat
        end do ! iDtype

        ! Calculate the overlap by double sum method
!        z0Dummy = cmplx(0., 0.)
!        z0Dummy(:nv(1, ikpt), 1) = z(:nv(1, ikpt), iband, ikpt, 1)
!        w_z0Dummy(:) = cmplx(0., 0.)
!
!        ! Warp the ket z0
!        call warpIRPot(stars, nv(1, ikpt), 1, gbas(:, mapGbas(:nv(1, ikpt), ikpt, 1)), z0Dummy, w_z0Dummy)
!
!        ! Calculate the interstitial integral by summing over the basis G-vectors
!        overlapIRkn = dot_product( z(:nv(1, ikpt), iband, ikpt, 1), w_z0Dummy(:nv(1, ikpt)) )
!
!        ! Sum of IR and MT overlap for every band and k-point is supposed to be 1 always.
!        write(*, '(2f15.7)') psiPsi(1, iband) + overlapIRkn
!        ! overlap works for all bands, sternheimer routine for the veff0 test somewhat creepy for p bands, not symmetric!, higher bands are not equal, s band is equal, discuss with Gustav, for the kinetic energy, we have to recheck it and compare it to the treatment in sternheimer, especially the 0.5 from the kinetic energy is added at the end of the routine and does not have to be added in the beginning. The symmetry of the IR kinetic energy should be compatible with the symmetry of the kinetic energy in the IR.
!
!        ! Kinetic energy
!        z0Dummy = cmplx(0., 0.)
!        w_z0Dummy(:) = cmplx(0., 0.)
!
!        do iG = 1, nv(1, ikpt)
!          kpgExt(1:3) = matmul(cell%bmat(1:3, 1:3), kpts%bk(1:3, ikpt) + gbas(1:3, mapGbas(iG, ikpt, 1)))
!          z0Dummy(iG, 1) = 0.5 * norm2(kpgExt(1:3))**2 * z(iG, iband, ikpt, 1)
!        end do !iG
!
!        call warpIRPot(stars, nv(1, ikpt), 1, gbas(:, mapGbas(:nv(1, ikpt), ikpt, 1)), z0Dummy, w_z0Dummy)
!
!        overlapIRkn = dot_product( z(:nv(1, ikpt), iband, ikpt, 1), w_z0Dummy(:nv(1, ikpt)) )
!
!        write(3003, '(2i8,4f15.8)') iband, ikpt, overlapIRkn
!
!        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!        z0Dummy = cmplx(0., 0.)
!        w_z0Dummy(:) = cmplx(0., 0.)
!        z0Dummy(:nv(1, ikpt), 1) = z(:nv(1, ikpt), iband, ikpt, 1)
!        call warpIRPot(stars, nv(1, ikpt), 1, gbas(:, mapGbas(:nv(1, ikpt), ikpt, 1)), z0Dummy, w_z0Dummy)
!
!        z0Dummy = cmplx(0., 0.)
!        do iG = 1, nv(1, ikpt)
!          kpgExt(1:3) = matmul(cell%bmat(1:3, 1:3), kpts%bk(1:3, ikpt) + gbas(1:3, mapGbas(iG, ikpt, 1)))
!          z0Dummy(iG, 1) = 0.5 * norm2(kpgExt(1:3))**2 * z(iG, iband, ikpt, 1)
!        end do !iG
!        overlapRevIRkn = dot_product( z0Dummy(:nv(1, ikpt), 1), w_z0Dummy(:nv(1, ikpt)) )
!        write(3004, '(2i8,4f15.8)') iband, ikpt, overlapRevIRkn
!
      end do ! iband
!      IRtest(:, :) = cmplx(0., 0.)
      call calcMEPotIR( stars, dimens, gbas(:, mapGbas(:nv(input%jspins, ikpt), ikpt, input%jspins)), &
        & gbas(:, mapGbas(:nv(input%jspins, ikpt), ikpt, input%jspins)), nv, vEff0IRpw, z(:, :, ikpt, 1), &
        & z(:, :, ikpt, 1), gdp, nmat, nobd(ikpt, 1), nobd(ikpt, 1), ikpt, 1, ikpt, ngdp, IRtest, kpq2kPrVec ) !todo spin-relevant
        do iband = 1, nobd(ikpt, 1)
          IRtestSummed(iband) = IRtestSummed(iband) + IRtest(iband, iband)
        end do ! iband
!      write(*, *) ikpt
    end do ! ikpt
    ! todo Note: p bands calculated by one method only should be the same for the Gamma point but not necessarily at the sum of all k-points
    ! in any case the two methods should deliver the same

    do iDatom = 1, atoms%nat
      do iband = 1, maxval(nobd)
        if (.false.) then
          write(*, '(2i8,6(es15.8,2x))') iband, iDatom, psiHepsPsiMT(iDatom, iband), wholeUnitCellIRkSum(iDatom, iband), psiHepsPsiMT(iDatom, iband) + wholeUnitCellIRkSum(iDatom, iband)!IRtestSummed(iband)
        end if
!        write(*, '(2i8,4(es15.8,2x))') iband, iDatom, IRtestSummed(iband), wholeUnitCellIRkSum(iDatom, iband)
        if (abs(psiHepsPsiMT(iDatom, iband) + wholeUnitCellIRkSum(iDatom, iband)) > 9e-9) passed = .false.
      end do
    end do ! iDatom

    if ( passed ) then
      write( logUnit, '(a)' ) '              |__ passed!'
    else
      write( logUnit, '(a)' ) '              |__ failed!'
      call JuDFT_warn('PsiHepsPsi test', calledby='testIR3rdMatElem', hint='Debug!')
    end if
  end subroutine testIR3rdMatElem

  ! Calculate <Psi|H - eps|Psi>_MT and compare it to the Sternheimer way of calculating it.
  subroutine testHeps1Phi( atoms, lathar, dimens, enpara, usdus, input, sym, kpts, cell, logUnit, mlh_atom, clnu_atom,      &
      & vEff0MtLh, rbas1, rbas2, uuilon, duilon, ulouilopn, ilo2p, nmem_atom, nRadFun, gbas, mapGbas, nv, kveclo, nobd, z, El,     &
      & iloTable, vEff0MTsh )

#include"cppmacro.h"

    use m_types, only : t_atoms, t_sphhar, t_sym, t_cell, t_kpts, t_dimension, t_usdus, t_enpara, t_input, t_tlmplm
     
    use m_abcof3
    use m_juDFT_NOstopNO, only : juDFT_warn
    use m_jpSternhPulaySurface, only : tlmplm4H0, calcHS0MT
    use m_jpSetupDynMat, only : CalcScalBasfMatElems
    use m_jpConstants, only : iu
    use m_jpSetupDynMatHelper, only : CalcHnGrV0Varphi
!
    implicit none

    ! Type parameter
    type(t_atoms),                  intent(in) :: atoms
    type(t_sphhar),                 intent(in) :: lathar
    type(t_dimension),              intent(in) :: dimens
    type(t_enpara),                 intent(in) :: enpara
    type(t_usdus),                  intent(in) :: usdus
    type(t_input),                  intent(in) :: input
    type(t_sym),                    intent(in) :: sym
    type(t_kpts),                   intent(in) :: kpts
    type(t_cell),                   intent(in) :: cell

    ! Scalar parameter
    integer,                        intent(in) :: logUnit

    ! Array parameter
    integer,                        intent(in) :: mlh_atom(:, 0:, :)
    complex,                        intent(in) :: clnu_atom(:, 0:, :)
    real,                           intent(in) :: vEff0MtLh(:, 0:, :)
    real,                           intent(in) :: rbas1(:,:,0:, :, :)
    real,                           intent(in) :: rbas2(:,:,0:, :, :)
    real,                           intent(in) :: uuilon(:, :)
    real,                           intent(in) :: duilon(:, :)
    real,                           intent(in) :: ulouilopn(:, :, :)
    integer,                        intent(in) :: ilo2p(:, :)
    integer,                        intent(in) :: nmem_atom(0:, :)
    integer,                        intent(in) :: nRadFun(0:, :)
    integer,                        intent(in) :: gbas(:, :)
    integer,                        intent(in) :: mapGbas(:, :, :)
    integer,                        intent(in) :: nv(:, :)
    integer,                        intent(in) :: kveclo(:,:)
    integer,                        intent(in) :: nobd(:, :)
    MCOMPLEX,                       intent(in) :: z(:,:,:,:)
    real,                           intent(in) :: El(:, 0:, :, :)
    integer,                        intent(in) :: iloTable(:, 0:, :)
    complex,                        intent(in) :: vEff0MTsh(:, :, :, :)

    ! Type variable
    type(od_inp)                               :: odi
    type(od_sym)                               :: ods
    type(t_tlmplm)                             :: td4HS0

    ! Scalar variables
    integer                                    :: iatom
    integer                                    :: itype
    integer                                    :: ieqat
    integer                                    :: ptsym
    integer                                    :: ilh
    integer                                    :: oqn_l
    integer                                    :: lm_pre
    integer                                    :: imem
    integer                                    :: mqn_m
    integer                                    :: lm
    integer                                    :: imesh
    integer                                    :: lmpMax
    integer                                    :: nRadFunMax
    integer                                    :: passed
    integer                                    :: ikpt
    integer                                    :: nmat
    integer                                    :: iband
    integer                                    :: iDtype
    integer                                    :: iDatom
    integer                                    :: iDeqat
    integer                                    :: lmp
    integer                                    :: pMaxLocal
    integer                                    :: iradf
    complex                                    :: matEntryRef

    ! Array variabels
    integer,           allocatable             :: ngoprI(:)
    real,              allocatable             :: varphi1(:, :, :)
    real,              allocatable             :: varphi2(:, :, :)
    complex,           allocatable             :: vEff0MtSpH(:, :, :)
    complex,           allocatable             :: hVarphi(:, :, :, :)
    integer,           allocatable             :: lmpT(:)
    real,              allocatable             :: r2(:)
    complex,           allocatable             :: a(:, :, :)
    complex,           allocatable             :: b(:, :, :)
    complex,           allocatable             :: bascof_lo(:, :, :, :, :)
    complex,           allocatable             :: ab0cofSH(:, :, :)
    complex,           allocatable             :: ab0cofScl(:)
    complex,           allocatable             :: varphiPsi(:)
    complex,           allocatable             :: psiPsi(:, :)
    complex,           allocatable             :: varphiHPsi(:)
    complex,           allocatable             :: psiHpsi(:, :)
    complex,           allocatable             :: h0MTRef(:, :)
    complex,           allocatable             :: s0MTRef(:, :)
    integer,           allocatable             :: nlo_atom(:)
    real,          allocatable              :: varphiVarphi(:, :, :)
    complex,       allocatable              :: varphiHvarphi(:, :, :)
    complex,       allocatable              :: vEff0NsphGrVarphi(:, :, :, :, :)
    complex,       allocatable              :: r2grVeff0SphVarphi(:, :, :, :, :)
    complex,       allocatable              :: r2grVeff0SphVarphiDummy(:, :, :, :, :)
    integer,       allocatable              :: grVarphiChLout(:, :)
    integer,       allocatable              :: grVarphiChMout(:, :)
    real,          allocatable              :: grVarphiCh1(:, :, :, :)
    real,          allocatable              :: grVarphiCh2(:, :, :, :)

    write( logUnit, '(a)' ) 'PsiHepsPsi test'
    write( logUnit, '(a)' ) '---------------'

    ! We do not want the local coordinate systems to be rotated for non-representative atoms constructing the matching coefficients.
    allocate(ngoprI(atoms%nat))
    allocate( lmpT(atoms%ntype) )
    lmpT(:) = 0
    do itype = 1, atoms%ntype
      lmpT(itype) = sum( [ ( (2 * oqn_l + 1)* nRadFun(oqn_l, itype), oqn_l = 0, atoms%lmax(itype) ) ] )
    end do ! itype
    lmpMax     = maxval( lmpT(:) )
    nRadFunMax = maxval( nRadFun(:, :) )
!    ! The a are dimensioned until lmax + 2 for the tested routines to work, but we only need the abcof until lmax for this test
    allocate( varphi1(atoms%jmtd, nRadFunMax, 0:atoms%lmaxd), varphi2(atoms%jmtd, nRadFunMax, 0:atoms%lmaxd) )
    allocate( hVarphi(2, atoms%jmtd, 0:(atoms%lmaxd + 3)**2 - 1, lmpMax) )
    allocate( r2(atoms%jmtd) )
    allocate( a( dimens%nvd, 0:(atoms%lmaxd + 3)**2 - 1, atoms%nat), b(dimens%nvd, 0:(atoms%lmaxd + 3)**2 - 1, atoms%nat), &
      & bascof_lo(3, -atoms%llod:atoms%llod, 4 * atoms%llod + 2, atoms%nlod, atoms%nat) )
    allocate( varphiVarphi(lmpMax, lmpMax, atoms%ntype), varphiHvarphi(lmpMax, lmpMax, atoms%nat))
    allocate( r2grVeff0SphVarphi(2, atoms%jmtd, 0:(atoms%lmaxd + 3)**2 - 1, lmpMax, 3) )
    allocate( r2grVeff0SphVarphiDummy(2, atoms%jmtd, 0:(atoms%lmaxd + 3)**2 - 1, lmpMax, 3) )
    allocate( vEff0NsphGrVarphi(2, atoms%jmtd, 0:(atoms%lmaxd + 3)**2 - 1, lmpMax, -1:1) )
    allocate( grVarphiCh1(atoms%jmtd, 2, lmpMax, -1:1), grVarphiCh2(atoms%jmtd, 2, lmpMax, -1:1), &
            & grVarphiChLout(2, 0:atoms%lmaxd), grVarphiChMout(-atoms%lmaxd:atoms%lmaxd, -1:1) )
    allocate(ab0cofSH(maxval(nobd), lmpMax, atoms%nat))
    allocate(ab0cofScl(lmpMax))
    allocate(varphiPsi(lmpMax))
    allocate(psiPsi(atoms%nat, maxval(nobd)))
    allocate(varphiHPsi(lmpMax))
    allocate(psiHpsi(atoms%nat, maxval(nobd)))
    allocate(h0MTRef(maxval(nobd), maxval(nobd)), s0MTRef(maxval(nobd), maxval(nobd)))
        allocate(nlo_atom(atoms%nat))
        nlo_atom = 0

    varphiVarphi(:, :, :)            = 0.
    varphiHvarphi(:, :, :)              = cmplx(0., 0.)
    r2grVeff0SphVarphi(:, :, :, :, :) = cmplx(0., 0.)
    r2grVeff0SphVarphiDummy(:, :, :, :, :) = cmplx(0., 0.)
    vEff0NsphGrVarphi(:, :, :, :, :) = cmplx(0., 0.)
    grVarphiCh1(:, :, :, :) = 0.
    grVarphiCh2(:, :, :, :) = 0.
    grVarphiChLout(:, :) = 0
    grVarphiChMout(:, :) = 0

    ! Calculate k-independent quantities

    ! Transform unperturbed potential given in lattice harmonics into spherical harmonic representation
    allocate(vEff0MtSpH(atoms%jmtd, (atoms%lmaxd + 1)**2, atoms%nat))
    vEff0MTSpH = cmplx(0.0, 0.0)
    iatom = 0
    do itype = 1, atoms%ntype
      do ieqat = 1, atoms%neq(itype)
        iatom = iatom + 1
        ptsym = atoms%ntypsy(iatom)
        do ilh = 0, lathar%nlh(ptsym)
          oqn_l = lathar%llh(ilh, ptsym)
          lm_pre = oqn_l * (oqn_l + 1) + 1
          do imem = 1, nmem_atom(ilh, iatom)
            mqn_m = mlh_atom(imem, ilh, iatom)
            lm = lm_pre + mqn_m
            do imesh = 1, atoms%jri(itype)
              vEff0MtSpH(imesh, lm, iatom) = vEff0MtSpH(imesh, lm, iatom) + vEff0MtLh(imesh, ilh, itype) * clnu_atom(imem, ilh, iatom)
            end do ! imesh
          end do ! imem
        end do ! ilh
      end do ! ieqat
    end do ! itype


    ! Matrix elements for the creation of the reference matrix element
    call tlmplm4H0( atoms, dimens, enpara, usdus, input, td4HS0, 1, logUnit, rbas1, rbas2, uuilon, duilon, ulouilopn, ilo2p, vEff0MTSpH )

    iatom = 0
    do itype = 1, atoms%ntype
      r2(:) = 0.

      ! Precalculate radial Jacobi determinant for later integrals
      do imesh = 1, atoms%jri(itype)
        r2(imesh) = atoms%rmsh(imesh, itype) * atoms%rmsh(imesh, itype)
      end do ! imesh

      do ieqat = 1, atoms%neq(itype)
        varphi1(:, :, :) = 0.
        varphi2(:, :, :) = 0.
        do oqn_l = 0, atoms%lmax(itype)
          do iradf = 1, nRadFun(oqn_l, itype)
            do imesh = 1, atoms%jri(itype)
              ! In Fleur the radial solutions u_lp are multiplied by a factor r per default to avoid an additional multiplication of the
              ! Jacobi determinant r^2 in radial integrals given spherical coordinates.
              varphi1(imesh, iradf, oqn_l) = rbas1(imesh, iradf, oqn_l, itype, 1) / atoms%rmsh(imesh, itype)
              varphi2(imesh, iradf, oqn_l) = rbas2(imesh, iradf, oqn_l, itype, 1) / atoms%rmsh(imesh, itype)
            end do ! imesh
          end do ! iradf
        end do ! oqn_l
        iatom = iatom + 1
        hVarphi = cmplx(0.0, 0.0)
        call CalcHnGrV0Varphi( atoms, lathar, itype, iatom, lmpMax, El, varphi1, varphi2, nRadFun, vEff0MtSpH(:, :, iatom), vEff0MtLh, clnu_atom, &
          & nmem_atom, mlh_atom, grVarphiCh1, grVarphiCh2, grVarphiChLout, grVarphiChMout, hVarphi, vEff0NsphGrVarphi, r2grVeff0SphVarphi, r2grVeff0SphVarphiDummy )
        call CalcScalBasfMatElems( atoms, itype, iatom, nRadFun, r2, varphi1, varphi2, hVarphi, varphiVarphi, varphiHvarphi )
      end do ! ieqat
    end do ! itype
    deallocate(r2grVeff0SphVarphiDummy)

    ! How dangerous is interpolation at first mesh points?
!    lmp = 0
!    do oqn_l = 0, atoms%lmax(1)
!      do mqn_m = -oqn_l, oqn_l
!
!        !p = 1
!        lmp = lmp + 1
!        write(2292, '(2f15.8)') varphiVarphi(lmp, lmp, 1), 1.
!
!        ! p = 2
!        lmp = lmp + 1
!        write(2292, '(2f15.8)') varphiVarphi(lmp, lmp, 1), usdus%ddn(oqn_l, 1, 1)
!
!      end do ! mqn_m
!    end do ! oqn_l
!    NOstopNO

    ngoprI(:) = 1
    passed = .true.

    do ikpt = 1, kpts%nkpt

      nmat = nv(1, ikpt) + atoms%nlotot
      a(:, :, :) = cmplx(0.0, 0.0)
      b(:, :, :) = cmplx(0.0, 0.0)
      bascof_lo(:, :, :, :, :) = cmplx(0.0, 0.0)
      call abcof3( atoms%lmaxd, atoms%ntype, atoms%nat, sym%nop, dimens%nvd, dimens%jspd, 1, dimens%lmd, dimens%nbasfcn, &
        & atoms%llod, atoms%nlod, atoms%nlotot, sym%invtab, atoms%ntype, sym%mrot, ngoprI, atoms%taual, atoms%neq, atoms%lmax, &
        & atoms%rmt, cell%omtil, cell%bmat, cell%bbmat, kpts%bk(:, ikpt), gbas(1, mapGbas(:nv(1, ikpt), ikpt, 1)), &
        & gbas(2, mapGbas(:nv(1, ikpt), ikpt, 1)), gbas(3, mapGbas(:nv(1, ikpt), ikpt, 1)), nv(:, ikpt), nmat, &
        & usdus%us, usdus%dus, usdus%uds, usdus%duds, usdus%ddn, atoms%invsat, sym%invsatnr, usdus%ulos, usdus%uulon, usdus%dulon, &
        & usdus%dulos, atoms%llo, atoms%nlo, atoms%l_dulo, atoms%lapw_l, kveclo(:,ikpt), odi, ods, a, b, bascof_lo )

      psiHpsi(:, :) = cmplx(0., 0.)
      psiPsi(:, :) = cmplx(0., 0.)
      ab0cofScl(:) = cmplx(0., 0.)
      ab0cofSH(:, :, :) = cmplx(0., 0.)


      do iband = 1, nobd(ikpt, 1)
        iDatom = 0
        do iDtype = 1, atoms%ntype
          do iDeqat = 1, atoms%neq(iDtype)
            iDatom = iDatom + 1
            ab0cofScl(:) = cmplx(0., 0.)
            lmp = 0
            lm  = 0
            do oqn_l = 0, atoms%lmax(iDtype)
              do mqn_m = - oqn_l, oqn_l
                pMaxLocal = nRadFun(oqn_l, iDtype)
                ! p = 1
                ab0cofScl(lmp + 1) = iu**oqn_l * dot_product( conjg(z(:nv(1, ikpt), iband, ikpt, 1)), a(:nv(1, ikpt), lm, iDatom) )
                ! p = 2
                ab0cofScl(lmp + 2) = iu**oqn_l * dot_product( conjg(z(:nv(1, ikpt), iband, ikpt, 1)), b(:nv(1, ikpt), lm, iDatom) )
                ! Add LO contributions
                do iradf = 3, pMaxLocal
                  ! p = 1
                  ab0cofScl(lmp + 1) = ab0cofScl(lmp + 1) + iu**oqn_l * &
                    & dot_product( conjg(z(nv(1, ikpt) + 1:nv(1, ikpt) + atoms%nlotot, iband, ikpt, 1)), &
                                                         & bascof_lo(1, mqn_m, :atoms%nlotot, iloTable(iradf, oqn_l, iDtype), iDatom) )
                  ! p = 2
                  ab0cofScl(lmp + 2) = ab0cofScl(lmp + 2) + iu**oqn_l * &
                    & dot_product( conjg(z(nv(1, ikpt) + 1:nv(1, ikpt) + atoms%nlotot, iband, ikpt, 1)), &
                                                         & bascof_lo(2, mqn_m, :atoms%nlotot, iloTable(iradf, oqn_l, iDtype), iDatom) )
                  ! 2 < p < LOs for that l and that atom type
                  ab0cofScl(lmp + iradf) = iu**oqn_l * dot_product( conjg(z(nv(1, ikpt) + 1:nv(1, ikpt) + atoms%nlotot, iband, ikpt, 1)), &
                                                         & bascof_lo(3, mqn_m, :atoms%nlotot, iloTable(iradf, oqn_l, iDtype), iDatom) )
                end do ! iradf


                ! This is a precalculation of the 1st and 3rd line of A.50. Actually the resulting lmp should have been primed or at
                ! least the p as lm1Pr = lm. But for sake of performance we place it here.
                ! sake of performance
                lm = lm + 1
                lmp = lmp + pMaxLocal
              end do ! mqn_m
            end do !oqn_l

            ! The i^l is performed within the calcHS0Mat
            lmp = 0
            lm  = 0
            do oqn_l = 0, atoms%lmax(iDtype)
              do mqn_m = - oqn_l, oqn_l
                pMaxLocal = nRadFun(oqn_l, iDtype)
                ! p = 1
                ab0cofSH(iband, lmp + 1, iDatom) = dot_product( conjg(z(:nv(1, ikpt), iband, ikpt, 1)), a(:nv(1, ikpt), lm, iDatom) )
                ! p = 2
                ab0cofSH(iband, lmp + 2, iDatom) = dot_product( conjg(z(:nv(1, ikpt), iband, ikpt, 1)), b(:nv(1, ikpt), lm, iDatom) )
                ! Add LO contributions
                do iradf = 3, pMaxLocal
                  ! p = 1
                  ab0cofSH(iband, lmp + 1, iDatom) = ab0cofSH(iband, lmp + 1, iDatom) + &
                    & dot_product( conjg(z(nv(1, ikpt) + 1:nv(1, ikpt) + atoms%nlotot, iband, ikpt, 1)), &
                                                    & bascof_lo(1, mqn_m, :atoms%nlotot, iloTable(iradf, oqn_l, iDtype), iDatom) )
                  ! p = 2
                  ab0cofSH(iband, lmp + 2, iDatom) = ab0cofSH(iband, lmp + 2, iDatom) + &
                    & dot_product( conjg(z(nv(1, ikpt) + 1:nv(1, ikpt) + atoms%nlotot, iband, ikpt, 1)), &
                                                    & bascof_lo(2, mqn_m, :atoms%nlotot, iloTable(iradf, oqn_l, iDtype), iDatom) )
                  ! 2 < p < LOs for that l and that atom type
                  ab0cofSH(iband, lmp + iradf, iDatom) = dot_product( conjg(z(nv(1, ikpt) + 1:nv(1, ikpt) + atoms%nlotot, iband, ikpt, 1)), &
                                                    & bascof_lo(3, mqn_m, :atoms%nlotot, iloTable(iradf, oqn_l, iDtype), iDatom) )
                end do ! iradf
                lm = lm + 1
                lmp = lmp + pMaxLocal
              end do ! mqn_m
            end do !oqn_l
            varphiPsi(:) = cmplx(0., 0.)
            varphiPsi(1: lmpT(iDtype)) = matmul( varphiVarphi(1:lmpT(iDtype), 1:lmpT(iDtype), iDtype), ab0cofScl(1:lmpT(iDtype)) )
            psiPsi(iDatom, iband) = dot_product(ab0cofScl(1:lmpMax), varphiPsi(1:lmpMax))
            ! Precalculation of the 1st and 3rd line in A.50 PhD thesis A. Klueppelberg.
            varphiHPsi(:) = cmplx(0.0, 0.0)
            varphiHPsi(1:lmpT(iDtype)) = matmul( varphiHvarphi(1:lmpT(iDtype), 1:lmpT(iDtype), iDatom), ab0cofScl(1:lmpT(iDtype)) )
            psiHpsi(iDatom, iband) = dot_product(ab0cofScl(1:lmpT(iDtype)), varphiHPsi(1:lmpT(iDtype)))
          end do ! iDeqat
        end do ! iDtype
      end do ! iband
      ! Reference matrix element
      h0MTRef = cmplx(0.0, 0.0)
      s0MTRef = cmplx(0.0, 0.0)
      iDatom = 0
      do iDtype = 1, atoms%ntype
        do iDeqat = 1, atoms%neq(iDtype)
          iDatom = iDatom + 1
          call calcHS0MT( atoms, usdus, td4HS0, ikpt, ikpt, iDtype, iDatom, nobd(:, 1), nobd, El, ab0cofSH(:, :, iDatom), ab0cofSH(:, :, iDatom), &
            & nRadFun, iloTable, nlo_atom, s0MTRef, h0MTRef )
          do iBand = 1, nobd(ikpt, 1)
            if ( abs( h0MTRef(iband, iband) - psiHPsi(iDatom, iband) ) > 9e-6 ) then
              passed = .false.
              write(logUnit, '(a, i3, a, i3)') '--> H matrix element failed at k-point ', ikpt, ' and band ', iband
            end  if
            if ( abs( s0MTRef(iband, iband) - PsiPsi(iDatom, iband) ) > 1e-8 ) then
              passed = .false.
              write(logUnit, '(a, i3, a, i3)') '--> S matrix element failed at k-point ', ikpt, ' and band ', iband
            end  if
          end do ! iBand
        end do ! iDeqat
      end do ! iDtype
    end do ! ikpt

    if ( passed ) then
      write( logUnit, '(a)' ) '              |__ passed!'
    else
      write( logUnit, '(a)' ) '              |__ failed!'
      call JuDFT_warn('PsiHepsPsi test faiiled', calledby='testPsiHepsPsi', hint='Check logfile for more information!')
    end if

  end subroutine testHeps1Phi

  !Compare of IR integral grRho grVeff with IR restricted <grPsi|grVeff|Psi>'
  subroutine CalcGrVeffGrtRhoInt( atoms, cell, lathar, dimens, stars, Veff0, input, ngdp, memd_atom, clnu_atom, nmem_atom, mlh_atom,      &
      & rho0IRst, gdp, rho0MT )

    use m_types
    use m_jpConstants, only : iu
    use mod_juPhonUtils, only : calcGrR2FinLH
    use m_jpPotDensHelper, only : CalcIRdVxcKern, CalcMTdVxcKern, warpIRPot, convertStar2G
    use m_jpGrVeff0, only : GenGrVeff0

    implicit none

    ! Type parameters
    type(t_atoms),                  intent(in) :: atoms
    type(t_cell),                   intent(in) :: cell
    type(t_sphhar),                 intent(in) :: lathar
    type(t_dimension),              intent(in) :: dimens
    type(t_stars),                  intent(in) :: stars
    type(t_potential),              intent(in) :: Veff0
    type(t_input),                  intent(in) :: input

    ! Scalar parameter
    integer,                        intent(in) :: ngdp
    integer,                        intent(in) :: memd_atom

    ! Array parameters
    complex,                        intent(in) :: clnu_atom(:, 0:, :)
    integer,                        intent(in) :: nmem_atom(0:, :)
    integer,                        intent(in) :: mlh_atom(:, 0:, :)
    complex,                        intent(in) :: rho0IRst(:, :)
    integer,                        intent(in) :: gdp(:, :)
    real,                           intent(in) :: rho0MT(:, 0:, :, :)

    ! Scalar variables
    integer                                    :: iG
    integer                                    :: idirC
    integer                                    :: idir
    integer                                    :: iatom
    integer                                    :: itype
    integer                                    :: ieqat
    integer                                    :: oqn_l
    integer                                    :: lm_pre
    integer                                    :: mqn_m
    integer                                    :: lm
    integer                                    :: ilh
    integer                                    :: imesh
    integer                                    :: ptsym
    integer                                    :: imem
    logical                                    :: numericalGradient
    logical                                    :: harSw =.true.
    logical                                    :: extSw =.true.
    logical                                    :: xcSw =.true.
    logical                                    :: testGoldstein
    logical                                    :: fullGrVeff0

    ! Array variables
    complex,           allocatable             :: grRho0IR(:, :)
    real,              allocatable             :: r2Veff0MT(:, :, :)
    complex,           allocatable             :: grVeff0IR(:, :)
    complex,           allocatable             :: vEff0IRpw(:)
    complex,           allocatable             :: rho0IRpw(:, :)
    real,              allocatable             :: r2Rho0MTlh(:, :, :, :)
    real,              allocatable             :: gaussWghts(:) ! gaussian weights belonging to gausPts
    complex,           allocatable             :: ylm(:, : )
    real,              allocatable             :: dKernMTGPts(:, :, :)
    complex,           allocatable             :: grVxcIRKern(:)
    complex,           allocatable             :: rho1IRContainer(:, :, :)
    complex,           allocatable             :: rho1MTContainer(:, :, :, :, :)
    complex,           allocatable             :: w_grVeff0IR(:, :)
    complex,           allocatable             :: w_vEff1IRContainer(:, :, :)
    complex,           allocatable             :: vEff1MTContainer(:, :, :, :, :)
    complex,           allocatable             :: grVeff0MT(:, :, :, :)
    complex,           allocatable             :: r2GrVeff0MT(:, :, :, :)
    complex,           allocatable             :: r2GrRho0MTsh(:, :, :, :)
    complex,           allocatable             :: grRho0MT(:, :, :, :)
    real,              allocatable             :: r2Rho0MT(:, :, :)
    complex,           allocatable             :: r2GrRho0MT(:, :, :, :)
    complex,           allocatable              :: rho0MTsh(:, :, :, :)
    complex                                    :: convIntgrl(3, 3)
    real                                       :: Gext(3)
    complex                                    :: integralIR(3, 3)
    complex                                    :: integralMT(3, 3)

    convIntgrl(:, :) = cmplx(0., 0.)


    !todo this code is copied and can be outsourced to a routine
    ! One can choose whether to calculate the numerical gradient of the unperturbed potential or to use the Weinert method for
    ! calculating it. The Weinert method has the advantage that continuity is enforced by construction at the MT boundary.
    ! It has been tested that both methods deliver the same numbers except for the just named behavior.
    numericalGradient = .false.
    if (numericalGradient) then

      ! Calculate the numerical gradient of the unperturbed effective potential.
      ! ------------------------------------------------------------------------

      ! Calculate the IR gradient.
      !   Convert the unperturbed effective potential from star expansion coefficients to plane-wave expansion coefficients
      allocate( vEff0IRpw(ngdp) )
      vEff0IRpw(:) = cmplx(0., 0.)
      ! SPIN MISSING
      call ConvertStar2G( Veff0%vpw_uw(:, 1), vEff0IRpw(:), stars, ngdp, gdp )

      ! Perform analytical gradient of effective potential in the interstitial region
      allocate( grVeff0IR(ngdp, 3) )
      grVeff0IR(:, :) = cmplx(0., 0.)
      do iG = 1, ngdp
        Gext(1:3) = matmul(cell%bmat(1:3, 1:3), gdp(1:3, iG))
        grVeff0IR(iG, 1:3) = iu * Gext(1:3) * vEff0IRpw(iG)
      end do ! iG

      ! Calculate the MT gradient.
      !   Any factors as in Fleur have been removed. In order to improve the numerical accuracy of the gradient r^2 is multiplied
      !   before and removed after the calculation of the gradient. Thereyby, subtraction of small numbers close to the core are
      !   avoided.
      allocate( r2Veff0MT( atoms%jmtd, 0:lathar%nlhd, atoms%ntype) )
      allocate( grVeff0MT( atoms%jmtd, ( atoms%lmaxd + 1 )**2, atoms%nat, 3) )
      r2Veff0MT(:, :, :) = 0.
      grVeff0MT(:, :, :, :) = cmplx(0., 0.)

      do itype = 1, atoms%ntype
        do ilh = 0, lathar%nlhd
          do imesh = 1, atoms%jri(itype)
            ! SPIN MISSING
            r2Veff0MT(imesh, ilh, itype) = Veff0%vr(imesh, ilh, itype, 1) * atoms%rmsh(imesh, itype) * atoms%rmsh(imesh, itype)
          end do ! imesh
        end do ! ilh
      end do ! itype

      call CalcGrR2FinLH( atoms, lathar, clnu_atom, nmem_atom, mlh_atom, r2Veff0MT, r2GrVeff0MT )

      do idirC = 1, 3
        iatom = 0
        do itype = 1, atoms%ntype
          do ieqat = 1, atoms%neq(itype)
            iatom = iatom + 1
            do oqn_l = 0, atoms%lmax(itype)
              lm_pre = oqn_l * (oqn_l + 1) + 1
              do mqn_m = -oqn_l, oqn_l
                lm = lm_pre + mqn_m
                do imesh = 1, atoms%jri(itype)
                  grVeff0MT(imesh, lm, iatom, idirC) = r2GrVeff0MT(imesh, lm, iatom, idirC) / atoms%rmsh(imesh, itype)**2
                end do ! imesh
              end do ! mqn_m
            end do ! oqn_l
          end do ! ieqat
        end do ! itype
      end do ! idirC

    else

      ! Use the Weinert method to calculate the gradient of the unperturbed effective potential
      ! ---------------------------------------------------------------------------------------

      ! Calculate the gradient of the unperturbed density.
      !   Any factors as in Fleur have been removed. In order to improve the numerical accuracy of the gradient r^2 is multiplied
      !   before and removed after the calculation of the gradient. Thereyby, subtraction of small numbers close to the core are
      !   avoided.
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
            lm_pre = oqn_l * (oqn_l + 1)
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
      call convertStar2G( rho0IRst(:, 1), rho0IRpw(:, 1), stars, ngdp, gdp )
      grRho0IR(:, :) = cmplx(0., 0.)
      do idir = 1, 3
        do iG = 1, ngdp
          Gext(1:3) = matmul(cell%bmat(1:3, 1:3), gdp(1:3, iG))
          grRho0IR(iG, idir)  = iu * Gext(idir) * rho0IRpw(iG, 1)
        end do ! iG
      end do ! idir

      allocate( r2Rho0MT( atoms%jmtd, 0:lathar%nlhd, atoms%ntype) )
      allocate( grRho0MT( atoms%jmtd, ( atoms%lmaxd + 1 )**2, atoms%nat, 3) )
      r2Rho0MT(:, :, :) = 0.
      grRho0MT(:, :, :, :) = cmplx(0., 0.)

      do itype = 1, atoms%ntype
        do ilh = 0, lathar%nlhd
          do imesh = 1, atoms%jri(itype)
            ! SPIN MISSING
            r2Rho0MT(imesh, ilh, itype) = rho0MT(imesh, ilh, itype, 1) * atoms%rmsh(imesh, itype) * atoms%rmsh(imesh, itype)
          end do ! imesh
        end do ! ilh
      end do ! itype

      call CalcGrR2FinLH( atoms, lathar, clnu_atom, nmem_atom, mlh_atom, r2Rho0MT, r2GrRho0MT )

      do idirC = 1, 3
        iatom = 0
        do itype = 1, atoms%ntype
          do ieqat = 1, atoms%neq(itype)
            iatom = iatom + 1
            do oqn_l = 0, atoms%lmax(itype)
              lm_pre = oqn_l * (oqn_l + 1) + 1
              do mqn_m = -oqn_l, oqn_l
                lm = lm_pre + mqn_m
                do imesh = 1, atoms%jri(itype)
                  grRho0MT(imesh, lm, iatom, idirC) = r2GrRho0MT(imesh, lm, iatom, idirC) / atoms%rmsh(imesh, itype)**2
                end do ! imesh
              end do ! mqn_m
            end do ! oqn_l
          end do ! ieqat
        end do ! itype
      end do ! idirC

      ! Precalculate quantities for the kernel of the unperturbed effective potential's gradient.
      ! SPIN MISSING
      call CalcIRdVxcKern(stars, gdp, ngdp, rho0IRst(:, 1), grVxcIRKern)
      ! SPIN MISSING
      call CalcMTdVxcKern(atoms, dimens, lathar, rho0MT(:, :, :, 1), nmem_atom, clnu_atom, mlh_atom, gaussWghts, ylm, dKernMTGPts)

      ! Generates gradient of the unperturbed effective potential by using the Weinert method.
      harSw = .true.
      extSw = .true.
      xcSw = .true.
      fullGrVeff0 = .true.
      testGoldstein = .false.
      ! SPIN MISSING
      call GenGrVeff0(atoms, cell, dimens, stars, ngdp, harSw, extSw, xcSw, gdp, rho0IRpw, rho0MTsh, grRho0IR, grRho0MT, gaussWghts, ylm, &
        & dKernMTGPts, grVxcIRKern, testGoldstein, fullGrVeff0, grVeff0IR, grVeff0MT )
      deallocate(grRho0MT)

    end if ! numericalGradient

    allocate(w_grVeff0IR(ngdp, 3))
    w_grVeff0IR = cmplx(0.0, 0.0)
    do idirC = 1, 3
      call warpIRPot(stars, ngdp, idirC, gdp, grVeff0IR, w_grVeff0IR(:, idirC))
    end do ! idirC

    ! Convert from star expansion coefficients to plane-wave expansion coefficients
 !   allocate( rho0IRpw(ngdp) )
 !   rho0IRpw(:) = cmplx(0., 0.)
 !   call convertStar2G( rho0IRst(:, 1), rho0IRpw(:), stars, ngdp, gdp )

 !   ! Perform analytical gradient of unperturbed density and effective potential in the interstitial region
 !   allocate( grRho0IR(ngdp, 3) )
 !   grRho0IR(:, :) = cmplx(0., 0.)
 !   do iG = 1, ngdp
 !     Gext(1:3) = matmul(cell%bmat(1:3, 1:3), gdp(1:3, iG))
 !     grRho0IR(iG, 1:3)  = iu * Gext(1:3) * rho0IRpw(iG)
 !   end do ! iG

    ! Gradient of density and effective potential in MT
    allocate( r2Rho0MTlh( atoms%jmtd, 0:lathar%nlhd, atoms%ntype, input%jspins) )
    allocate( grRho0MT( atoms%jmtd, ( atoms%lmaxd + 1 )**2, atoms%nat, 3) )
    r2Rho0MTlh(:, :, :, :) = 0.
    grRho0MT(:, :, :, :) = cmplx(0., 0.)

    ! Read in valence density form FLEUR multiplied with r2
    rewind(1040)
    read(1040) r2Rho0MtLh

    call calcGrR2FinLH( atoms, lathar, clnu_atom, nmem_atom, mlh_atom, r2Rho0MTlh(:, :, :, 1), r2GrRho0MTsh )

    do idirC = 1, 3
      iatom = 0
      do itype = 1, atoms%ntype
        do ieqat = 1, atoms%neq(itype)
          iatom = iatom + 1
          do oqn_l = 0, atoms%lmax(itype)
            lm_pre = oqn_l * (oqn_l + 1) + 1
            do mqn_m = -oqn_l, oqn_l
              lm = lm_pre + mqn_m
              do imesh = 1, atoms%jri(itype)
                ! For reasons of performance, this calculation is done within the idirC loop
                grRho0MT(imesh, lm, iatom, idirC) = r2GrRho0MTsh(imesh, lm, iatom, idirC) / atoms%rmsh(imesh, itype)**2
              end do ! imesh
            end do ! mqn_m
          end do ! oqn_l
        end do ! ieqat
      end do ! itype
    end do ! idirC


    ! Fill container arrays with gradient coefficients so that they can be processed by EvelIntRho1Veff1
    allocate( rho1IRContainer(ngdp, 3, atoms%nat), w_vEff1IRContainer(ngdp, 3, atoms%nat), &
      & rho1MTContainer( atoms%jmtd, (atoms%lmaxd + 1)**2, atoms%nat, 3, atoms%nat ), &
      & vEff1MTContainer( atoms%jmtd, (atoms%lmaxd + 1)**2, atoms%nat, 3, atoms%nat ) )

    w_vEff1IRContainer(:, :, :) = cmplx(0., 0.)
    rho1IRContainer(:, :, :) = cmplx(0., 0.)
    rho1MTContainer(:, :, :, :, :) = cmplx(0., 0.)
    vEff1MTContainer(:, :, :, :, :) = cmplx(0., 0.)

    ! todo think about this for more than one atom

    do idirC = 1, 3
      do iG = 1, ngdp
        rho1IRContainer(iG, idirC, 1) = grRho0IR(iG, idirC)
        write(4000, '(2i8,2f15.8)') iG, idirC, rho1IRContainer(iG, idirC, 1)
        w_vEff1IRContainer(iG, idirC, 1) = w_grVeff0IR(iG, idirC)
        write(4001, '(2i8,2f15.8)') iG, idirC, w_vEff1IRContainer(iG, idirC, 1)
      end do ! iG
    end do ! idirC

    !todo does not work for more than one atom
    do idirC = 1, 3
      ! There should be an itype
      do oqn_l = 0, atoms%lmax(1)
        do mqn_m = -oqn_l, oqn_l
          lm = oqn_l * (oqn_l + 1) + 1 + mqn_m
          do imesh = 1, atoms%jri(1)
            rho1MTContainer(imesh, lm, 1, idirC, 1) = grRho0MT(imesh, lm, 1, idirC)
            write(4002, '(3i8,2f15.8)') imesh, lm, idirC, rho1MTContainer(imesh, lm, 1, idirC, 1)
            vEff1MTContainer(imesh, lm, 1, idirC, 1) = grVeff0MT(imesh, lm, idirC, 1)
            write(4003, '(3i8,2f15.8)') imesh, lm, idirC, vEff1MTContainer(imesh, lm, 1, idirC, 1)
          end do ! imesh
        end do ! mqn_m
      end do ! oqn_l
    end do ! idirC

    ! todo should be substituted by Calc2ArgIntIR and Calc2ARgCmplxIntMT
    !call EvalIntRho1Veff1(atoms, cell, rho1IRContainer, rho1MTContainer, w_vEff1IRContainer, vEff1MTContainer, ConvIntgrl)

    write(*, '(a)') 'int rho GrVeff'
    write(*, '(3(2(es16.8,1x),3x))') convIntgrl(1, :)
    write(*, '(3(2(es16.8,1x),3x))') convIntgrl(2, :)
    write(*, '(3(2(es16.8,1x),3x))') convIntgrl(3, :)
  end subroutine CalcGrVeffGrtRhoInt

  ! There is a formula with Gaunt coefficients testing the muffin-tin gradient of the density, namely every l and m component. which gives an impression whether we have correctly set up the grVarphi which can me made to a grRho
  subroutine testGrVarphi( atoms, sym, dimens, lathar, stars, Veff0, kpts, cell, input, usdus, results, ngdp, logUnit, rho0IRst, rho0Mtlh, gdp,      &
      & clnu_atom, nmem_atom, mlh_atom, nRadFun, rbas1, rbas2, nv, ilst, GbasVec, nobd, z, kpq2kPrVec, El, kveclo, iloTable )

    use m_types, only : t_atoms, t_sym, t_dimension, t_sphhar, t_stars, t_potential, t_kpts, t_cell, t_input, t_usdus, t_results
    use m_jpConstants, only : iu, Tmatrix, fpi
    use mod_juPhonUtils, only : Derivative, CalcChannelsGrFlpNat, calcGrR2FinLH
     
    use m_abcof3
    use m_gaunt, only : gaunt1
    use m_juDFT_NOstopNO, only : juDFT_warn

    !implicit none

    ! Type parameters
    type(t_atoms),                  intent(in) :: atoms
    type(t_sym),                    intent(in) :: sym
    type(t_dimension),              intent(in) :: dimens
    type(t_sphhar),                 intent(in) :: lathar
    type(t_stars),                  intent(in) :: stars
    type(t_potential),              intent(in) :: Veff0
    type(t_kpts),                   intent(in) :: kpts
    type(t_cell),                   intent(in) :: cell
    type(t_input),                  intent(in) :: input
    type(t_usdus),                  intent(in) :: usdus
    type(t_results),                intent(in) :: results

    ! Scalar parameters
    integer,                        intent(in) :: ngdp
    integer,                        intent(in) :: logUnit

    ! Array parameters
    complex,                        intent(in) :: rho0IRst(:, :)
    real,                           intent(in) :: rho0Mtlh(:, 0:, :, :)
    integer,                        intent(in) :: gdp(:, :)
    complex,                        intent(in) :: clnu_atom(:, 0:, :)
    integer,                        intent(in) :: nmem_atom(0:, :)
    integer,                        intent(in) :: mlh_atom(:, 0:, :)
    integer,                        intent(in) :: nRadFun(0:, :)
    real,                           intent(in) :: rbas1(:,:,0:,:,:)
    real,                           intent(in) :: rbas2(:,:,0:,:,:)
    integer,                        intent(in) :: nv(:, :)
    integer,                        intent(in) :: ilst(:, :, :)
    integer,                        intent(in) :: GbasVec(:, :)
    integer,                        intent(in) :: nobd(:, :)
    MCOMPLEX,                       intent(in) :: z(:, :, :, :)
    integer,                        intent(in) :: kpq2kPrVec(:, :, :)
    real,                           intent(in) :: El(:, 0:, :, :)
    integer,                        intent(in) :: kveclo(:,:)
    integer,                        intent(in)    :: iloTable(:, 0:, :)

    ! Type variable
    type(od_inp)                                  :: odi
    type(od_sym)                                  :: ods

    ! Scalar variables
    integer                                    :: ikpt
    integer                                    :: iG
    integer                                    :: idirR
    integer                                    :: idirC
    integer                                    :: itype
    integer                                    :: ilh
    integer                                    :: imesh
    integer                                    :: oqn_l
    integer                                    :: mqn_m
    integer                                    :: lm
    integer                                    :: iatom
    integer                                    :: ieqat
    integer                                    :: lm_pre
    integer                                    :: lmpMax
    integer                                    :: nRadFunMax
    integer                                    :: iradf
    integer                                    :: iBas
    integer                                    :: iband
    integer                                    :: nmat
    integer                                    :: ispin
    integer                                    :: mqn_m2prC
    integer                                    :: lmp
    integer                                    :: pMaxLocal
    integer                                    :: idir
    integer                                    :: mqn_m2PrR
    integer                                    :: mqn_m4Pr
    integer                                    :: lmp2Pr
    integer                                    :: oqn_l2Pr
    integer                                    :: mqn_m2Pr
    integer                                    :: iradf2Pr
    integer                                    :: lmp1Pr
    integer                                    :: oqn_l1Pr
    integer                                    :: mqn_m1Pr
    integer                                    :: mqn_m3Pr
    integer                                    :: iradf1Pr
    integer                                    :: ichan
    integer                                    :: oqn_l3Pr
    integer                                    :: iradf3Pr
    integer                                    :: ptsym
    integer                                    :: imem
    logical                                    :: passed = .true.

    !! Array variables
    integer,           allocatable             :: lmpT(:)
    real,              allocatable             :: r2Rho0MTlh(:, :, :, :)
    complex,           allocatable             :: grRho0MT(:, :, :, :)
    complex,           allocatable             :: r2GrRho0MTsh(:, :, :, :)
    real,              allocatable             :: varphi1(:, :, :)
    real,              allocatable             :: varphi2(:, :, :)
    real,              allocatable             :: r2(:)
    real,              allocatable             :: delrVarphi1(:, :, :)
    real,              allocatable             :: delrVarphi2(:, :, :)
    integer,           allocatable             :: grVarphiChLout(:, :)
    integer,           allocatable             :: grVarphiChMout(:, :)
    real,              allocatable             :: grVarphiCh1(:, :, :, :)
    real,              allocatable             :: grVarphiCh2(:, :, :, :)
    real,              allocatable             :: grVarphiGauntVarphi(:, :, :, :, :, :)
    real,              allocatable             :: varphiGauntVarphi(:, :, :, :, :)
    complex,           allocatable             :: varphiGauntPsi(:)
    complex,           allocatable             :: psiGauntPsi(:, :)
    complex,           allocatable             :: grVarphiGauntPsiNat(:)
    complex,           allocatable             :: grPsiGauntPsiNat(:, :, :)
    complex,           allocatable             :: grPsiGauntPsi(:, :, :)
    complex,           allocatable             :: a(:, :, :)
    complex,           allocatable             :: b(:, :, :)
    complex,           allocatable             :: bascof_lo(:, :, :, :, :)
    integer,           allocatable             :: ngoprI(:)
    complex,           allocatable             :: ab0cofScl(:)
    complex,           allocatable             :: rho0MtSpH(:, :, :, :)
    real                                       :: Gext(3)
    real                                       :: kExt(3)

    real, allocatable :: rho0ValFleur(:, :, :, :)

      write( logUnit, '(a)' ) 'Testing the lm components of grRho0MT with Gaunt coefficient'
      write( logUnit, '(a)' ) '------------------------------------------------------------'


    ! Calculate the gradient of the density from wavefunctions
    !todo : also calculate the density from wavefunctions not only the gradient of it
    ! Quantities for initialization
    allocate( lmpT(atoms%ntype) )
    lmpT(:) = 0
    do itype = 1, atoms%ntype
      lmpT(itype) = sum( [ ( (2 * oqn_l + 1)* nRadFun(oqn_l, itype), oqn_l = 0, atoms%lmax(itype) ) ] )
    end do ! itype
    lmpMax     = maxval( lmpT(:) )
    nRadFunMax = maxval( nRadFun(:, :) )

    allocate( varphi1(atoms%jmtd, nRadFunMax, 0:atoms%lmaxd), varphi2(atoms%jmtd, nRadFunMax, 0:atoms%lmaxd) )
    allocate( grVarphiCh1(atoms%jmtd, 2, lmpMax, -1:1), grVarphiCh2(atoms%jmtd, 2, lmpMax, -1:1), &
            & grVarphiChLout(2, 0:atoms%lmaxd), grVarphiChMout(-atoms%lmaxd:atoms%lmaxd, -1:1) )
    allocate( delrVarphi1( atoms%jmtd, nRadFunMax, 0:atoms%lmaxd), delrVarphi2( atoms%jmtd, nRadFunMax, 0:atoms%lmaxd) )
    allocate( ab0cofScl(lmpMax) )
    allocate( grVarphiGauntPsiNat(lmpMax) )
    allocate( grVarphiGauntVarphi(lmpMax, lmpMax, atoms%jmtd, 2:4, -1:1, atoms%ntype))
    allocate( varphiGauntVarphi(lmpMax, lmpMax, atoms%jmtd, 1, atoms%ntype))
    allocate( varphiGauntPsi(lmpMax), psiGauntPsi(atoms%jmtd, 1) )
    allocate( rho0MtSpH( atoms%jmtd, (atoms%lmaxd + 1)**2, atoms%nat, 1 ))

    allocate(rho0ValFleur(atoms%jmtd, 0: lathar%nlhd, atoms%ntype, input%jspins))
    rho0ValFleur = 0
    rewind(1040)
    read(1040) rho0ValFleur
    iatom = 0
    rho0MtSpH = 0
    do itype = 1, atoms%ntype
      do ieqat = 1, atoms%neq(itype)
        iatom = iatom + 1
        ptsym = atoms%ntypsy(iatom)
        do ilh = 0, lathar%nlh(ptsym)
          oqn_l = lathar%llh(ilh, ptsym)
          do imem = 1, nmem_atom(ilh, iatom)
            mqn_m = mlh_atom(imem, ilh, iatom)
            lm = oqn_l * (oqn_l + 1) + mqn_m + 1
            do imesh = 1, atoms%jri(itype)
              rho0MtSpH(imesh, lm, iatom, 1) = rho0MtSpH(imesh, lm, iatom, 1) + rho0ValFleur(imesh, ilh, itype, 1) * clnu_atom(imem, ilh, iatom) / atoms%rmsh(imesh, itype) / atoms%rmsh(imesh, itype)
            end do ! imesh
            if ( lm == 1 ) then
              do imesh = 1, atoms%jri(itype)
!                rho0MtSpH(imesh, lm, iatom, 1) = rho0MtSpH(imesh, lm, iatom, 1)  * sqrt(fpi)
!todo this might be incorrect
              end do ! imesh

            end if
          end do ! imem
        end do ! ilh
      end do
    end do

    ! Gradient of density in MT
    allocate( r2Rho0MTlh( atoms%jmtd, 0:lathar%nlhd, atoms%ntype, 1) )
    allocate( grRho0MT( atoms%jmtd, ( atoms%lmaxd + 1 )**2, atoms%nat, 3) )
    r2Rho0MTlh(:, :, :, :) = 0.
    grRho0MT(:, :, :, :) = cmplx(0., 0.)


    call calcGrR2FinLH( atoms, lathar, clnu_atom, nmem_atom, mlh_atom, rho0ValFleur(:, :, :, 1), r2GrRho0MTsh )

    do idirC = 1, 3
      iatom = 0
      do itype = 1, atoms%ntype
        do ieqat = 1, atoms%neq(itype)
          iatom = iatom + 1
          do oqn_l = 0, atoms%lmax(itype)
            lm_pre = oqn_l * (oqn_l + 1) + 1
            do mqn_m = -oqn_l, oqn_l
              lm = lm_pre + mqn_m
              do imesh = 1, atoms%jri(itype)
                grRho0MT(imesh, lm, iatom, idirC) = r2GrRho0MTsh(imesh, lm, iatom, idirC) / atoms%rmsh(imesh, itype)**2
              end do ! imesh
            end do ! mqn_m
          end do ! oqn_l
        end do ! ieqat
      end do ! itype
    end do ! idirC

    grVarphiCh1(:, :, :, :) = 0.
    grVarphiCh2(:, :, :, :) = 0.
    grVarphiChLout(:, :) = 0
    grVarphiChMout(:, :) = 0
    grVarphiGauntPsiNat(:) = cmplx(0., 0.)
    ab0cofScl(:) = cmplx(0., 0.)
    grVarphiGauntVarphi(:, :, :, :, :, :) = 0.
    varphiGauntVarphi(:, :, :, :, :) = 0.
    varphiGauntPsi(:) = cmplx(0., 0.)
    psiGauntPsi(:, :) = cmplx(0., 0.)

    iatom = 0
    do itype = 1, atoms%ntype

      ! These arrays are set to zero here, because 
      varphi1(:, :, :) = 0.
      varphi2(:, :, :) = 0.
      delrVarphi1(:, :, :) = 0.
      delrVarphi2(:, :, :) = 0.
      do oqn_l = 0, atoms%lmax(itype)
        do iradf = 1, nRadFun(oqn_l, itype)
          do imesh = 1, atoms%jri(itype)
            ! In Fleur the radial solutions u_lp are multiplied by a factor r by default to avoid an additional multiplication of the
            ! Jacobi determinant r^2 in radial integrals given spherical coordinates.
            varphi1(imesh, iradf, oqn_l) = rbas1(imesh, iradf, oqn_l, itype, 1) / atoms%rmsh(imesh, itype)
            varphi2(imesh, iradf, oqn_l) = rbas2(imesh, iradf, oqn_l, itype, 1) / atoms%rmsh(imesh, itype)
          end do ! imesh
          ! Precalculate partial derivatives of varphis in r-direction since it is needed twice
          call Derivative( varphi1(1:atoms%jri(itype), iradf, oqn_l), itype, atoms, delrVarphi1(1:atoms%jri(itype), iradf, oqn_l) )
          call Derivative( varphi2(1:atoms%jri(itype), iradf, oqn_l), itype, atoms, delrVarphi2(1:atoms%jri(itype), iradf, oqn_l) )
        end do ! iradf
      end do ! oqn_l

      ! Calculate the application of the gradient and the gradient's dyadic product onto the MT basis functions (matching coefficients
      ! have no spatial dependence) and determing its scattering channels.
      grVarphiChLout(:, :) = 0
      grVarphiChMout(:, :) = 0
      grVarphiCh1(:, :, :, :) = 0.
      grVarphiCh2(:, :, :, :) = 0.
      call CalcChannelsGrFlpNat( atoms, itype, nRadFun, varphi1, varphi2, delrVarphi1, delrVarphi2, grVarphiChLout, grVarphiChMout, &
                                                                                                        & grVarphiCh1, grVarphiCh2 )

      do mqn_m4Pr = -1, 1
        do oqn_l = 1, 1!atoms%lmax(itype)
          do mqn_m = -oqn_l, oqn_l
            lm = oqn_l * (oqn_l + 1) + 1 + mqn_m
            ! natural coordinates
            lmp2Pr = 0
            do oqn_l2Pr = 0, atoms%lmax(itype)
              do mqn_m2Pr = -oqn_l2Pr, oqn_l2Pr
                do iradf2Pr = 1, nRadFun(oqn_l2Pr, itype)
                  lmp2Pr = lmp2Pr + 1
                  lmp1Pr = 0
                  do oqn_l1Pr = 0, atoms%lmax(itype)
                    do mqn_m1Pr = -oqn_l1Pr, oqn_l1Pr
                      mqn_m3Pr = grVarphiChMout(mqn_m1Pr, mqn_m4Pr)
                      do iradf1Pr = 1, nRadFun(oqn_l1Pr, itype)
                        lmp1Pr = lmp1Pr + 1
                        do ichan = 1, 2
                          oqn_l3Pr = grVarphiChLout(ichan, oqn_l1Pr)
                          if ( ( abs(mqn_m3Pr) > oqn_l3Pr )  .or. ( oqn_l3Pr < 0 ) .or. ( oqn_l3Pr > atoms%lmax(itype) ) )  cycle
                          do imesh = 1, atoms%jri(itype)
                            grVarphiGauntVarphi(lmp1Pr, lmp2Pr, imesh, lm, mqn_m4Pr, itype) = &
                              & grVarphiGauntVarphi(lmp1Pr, lmp2Pr, imesh, lm, mqn_m4Pr, itype) + &
                              & ( grVarphiCh1(imesh, ichan, lmp1Pr, mqn_m4Pr) * varphi1(imesh, iradf2Pr, oqn_l2Pr) &
                              & + grVarphiCh2(imesh, ichan, lmp1Pr, mqn_m4Pr) * varphi2(imesh, iradf2Pr, oqn_l2Pr) ) &
                                 & * gaunt1( oqn_l2Pr, oqn_l, oqn_l3Pr, mqn_m2Pr, mqn_m, mqn_m3Pr, atoms%lmax(itype) )
                          end do ! imesh
                        end do ! ichan
                      end do ! iradf1Pr
                    end do ! mqn_m1Pr
                  end do ! oqn_l1Pr
                end do ! iradf2Pr
              end do ! mqn_m2Pr
            end do ! oqn_l2Pr
          end do ! mqn_m
        end do ! oqn_l
      end do ! mqn_m4Pr

!      do oqn_l = 0, 0
!        do mqn_m = -oqn_l, oqn_l
!          lm = oqn_l * (oqn_l + 1) + 1 + mqn_m
!          lmp2Pr = 0
!          do oqn_l2Pr = 0, atoms%lmax(itype)
!            do mqn_m2Pr = -oqn_l2Pr, oqn_l2Pr
!              do iradf2Pr = 1, nRadFun(oqn_l2Pr, itype)
!                lmp2Pr = lmp2Pr + 1
!                lmp1Pr = 0
!                do oqn_l1Pr = 0, atoms%lmax(itype)
!                  do mqn_m1Pr = -oqn_l1Pr, oqn_l1Pr
!                    do iradf1Pr = 1, nRadFun(oqn_l1Pr, itype)
!                      lmp1Pr = lmp1Pr + 1
!                      do imesh = 1, atoms%jri(itype)
!                        varphiGauntVarphi(lmp1Pr, lmp2Pr, imesh, lm, itype) = varphiGauntVarphi(lmp1Pr, lmp2Pr, imesh, lm, itype) &
!                          & + ( varphi1(imesh, iradf1Pr, oqn_l1Pr) * varphi1(imesh, iradf2Pr, oqn_l2Pr) &
!                          & +   varphi2(imesh, iradf1Pr, oqn_l1Pr) * varphi2(imesh, iradf2Pr, oqn_l2Pr) ) &
!                          & * gaunt1(oqn_l2Pr, oqn_l, oqn_l1Pr, mqn_m2Pr, mqn_m, mqn_m1Pr, atoms%lmaxd)
!                      end do ! imesh
!                    end do ! iradf1Pr
!                  end do ! mqn_m1Pr
!                end do ! oqn_l1Pr
!              end do ! iradf2Pr
!            end do ! mqn_m2Pr
!          end do ! oqn_l2Pr
!        end do ! mqn_m
!      end do ! oqn_l

    end do ! itype

    ispin = 1
    allocate( a( dimens%nvd, 0:dimens%lmd, atoms%nat), b(dimens%nvd, 0:dimens%lmd, atoms%nat), &
      & bascof_lo(3, -atoms%llod:atoms%llod, 4 * atoms%llod + 2, atoms%nlod, atoms%nat) )
    allocate( ngoprI(atoms%nat) )
    a(:, :, :)               = cmplx(0., 0.)
    b(:, :, :)               = cmplx(0., 0.)
    bascof_lo(:, :, :, :, :) = cmplx(0., 0.)
    ! We do not want the local coordinate systems to be rotated for non-representative atoms constructing the matching coefficients.
    ngoprI(:) = 1
    iatom = 0
    allocate( grPsiGauntPsiNat( atoms%jmtd, (atoms%lmaxd + 1)**2, -1:1) )
    allocate( grPsiGauntPsi( atoms%jmtd, (atoms%lmaxd + 1)**2, 3) )
    grPsiGauntPsiNat = cmplx(0., 0.)
    grPsiGauntPsi = cmplx(0., 0.)
    do itype = 1, atoms%ntype
      do ieqat = 1, atoms%neq(itype)
        iatom = iatom + 1
        do ikpt = 1, kpts%nkpt
          nmat = nv(1, ikpt) + atoms%nlotot
          a(:, :, :) = cmplx(0.0, 0.0)
          b(:, :, :) = cmplx(0.0, 0.0)
          bascof_lo(:, :, :, :, :) = cmplx(0.0, 0.0)
          call abcof3( atoms%lmaxd, atoms%ntype, atoms%nat, sym%nop, dimens%nvd, dimens%jspd, 1, dimens%lmd, dimens%nbasfcn, &
            & atoms%llod, atoms%nlod, atoms%nlotot, sym%invtab, atoms%ntype, sym%mrot, ngoprI, atoms%taual, atoms%neq, atoms%lmax, &
            & atoms%rmt, cell%omtil, cell%bmat, cell%bbmat, kpts%bk(:, ikpt), gbasVec(1, ilst(:nv(1, ikpt), ikpt, 1)), &
            & gbasVec(2, ilst(:nv(1, ikpt), ikpt, 1)), gbasVec(3, ilst(:nv(1, ikpt), ikpt, 1)), nv(:, ikpt), nmat, &
            & usdus%us, usdus%dus, usdus%uds, usdus%duds, usdus%ddn, atoms%invsat, sym%invsatnr, usdus%ulos, usdus%uulon, usdus%dulon, &
            & usdus%dulos, atoms%llo, atoms%nlo, atoms%l_dulo, atoms%lapw_l, kveclo(:,ikpt), odi, ods, a, b, bascof_lo )


          do iband = 1, nobd(ikpt, 1)
            ab0cofScl(:) = cmplx(0., 0.)
            lmp = 0
            lm  = 0
            do oqn_l = 0, atoms%lmax(itype)
              do mqn_m = - oqn_l, oqn_l
                pMaxLocal = nRadFun(oqn_l, itype)
                ! p = 1
                ab0cofScl(lmp + 1) = iu**oqn_l * dot_product( conjg(z(:nv(1, ikpt), iband, ikpt, 1)), a(:nv(1, ikpt), lm, iatom) )
                ! p = 2
                ab0cofScl(lmp + 2) = iu**oqn_l * dot_product( conjg(z(:nv(1, ikpt), iband, ikpt, 1)), b(:nv(1, ikpt), lm, iatom) )
                ! Add LO contributions
                do iradf = 3, pMaxLocal
                  ! p = 1
                  ab0cofScl(lmp + 1) = ab0cofScl(lmp + 1) + iu**oqn_l * &
                    & dot_product( conjg(z(nv(1, ikpt) + 1:nv(1, ikpt) + atoms%nlotot, iband, ikpt, 1)), &
                                                         & bascof_lo(1, mqn_m, :atoms%nlotot, iloTable(iradf, oqn_l, itype), iatom) )
                  ! p = 2
                  ab0cofScl(lmp + 2) = ab0cofScl(lmp + 2) + &
                    & iu**oqn_l * dot_product( conjg(z(nv(1, ikpt) + 1:nv(1, ikpt) + atoms%nlotot, iband, ikpt, 1)), &
                                                         & bascof_lo(2, mqn_m, :atoms%nlotot, iloTable(iradf, oqn_l, itype), iatom) )
                  ! 2 < p < LOs for that l and that atom type
                  ab0cofScl(lmp + iradf) = iu**oqn_l * dot_product( conjg(z(nv(1, ikpt) + 1:nv(1, ikpt) + atoms%nlotot, iband, ikpt, 1)), &
                                                         & bascof_lo(3, mqn_m, :atoms%nlotot, iloTable(iradf, oqn_l, itype), iatom) )
                end do ! iradf


                ! This is a precalculation of the 1st and 3rd line of A.50. Actually the resulting lmp should have been primed or at
                ! least the p as lm1Pr = lm. But for sake of performance we place it here.
                ! sake of performance
                lm = lm + 1
                lmp = lmp + pMaxLocal
              end do ! mqn_m
            end do !oqn_l
            do mqn_m2PrR = -1, 1
              do oqn_l = 1, 1!atoms%lmax(itype)
                do mqn_m = -oqn_l, oqn_l
                  lm = oqn_l * (oqn_l + 1) + 1 + mqn_m
                  do imesh = 1, atoms%jri(itype)
                    grVarphiGauntPsiNat(:) = cmplx(0., 0.)
                    grVarphiGauntPsiNat(1:lmpT(itype)) = matmul(grVarphiGauntVarphi(1:lmpT(itype), 1:lmpT(itype), imesh, lm, mqn_m2PrR, itype), ab0cofScl(1:lmpT(itype)))
                    grPsiGauntPsiNat(imesh, lm, mqn_m2PrR) = grPsiGauntPsiNat(imesh, lm, mqn_m2PrR) + 2 * results%w_iks(iband, ikpt, 1) * dot_product(ab0cofScl(1:lmpT(itype)), grVarphiGauntPsiNat(1:lmpT(itype)))
                   end do ! imesh
                end do ! mqn_m
              end do ! oqn_l
            end do ! mqn_m2PrR
         !   do oqn_l = 0, 0
         !     do mqn_m = -oqn_l, oqn_l
         !       lm = oqn_l * (oqn_l + 1) + 1 + mqn_m
         !       do imesh = 1, atoms%jri(itype)
         !         varphiGauntPsi(:) = cmplx(0., 0.)
         !         varphiGauntPsi(1:lmpT(itype)) = matmul(varphiGauntVarphi(1:lmpT(itype), 1:lmpT(itype), imesh, lm, itype), ab0cofScl(1:lmpT(itype)))
         !         psiGauntPsi(imesh, lm) = psiGauntPsi(imesh, lm) + 2 * results%w_iks(iband, ikpt, 1) * dot_product(ab0cofScl(1:lmpT(itype)), varphiGauntPsi(1:lmpT(itype)))
         !       end do ! imesh
         !     end do ! mqn_m
         !   end do ! oqn_l
          end do ! iband
        end do !ikptR

        do oqn_l = 1, 1!atoms%lmax(itype)
          do mqn_m = -oqn_l, oqn_l
            lm = oqn_l * (oqn_l + 1) + 1 + mqn_m
            do imesh = 1, atoms%jri(itype)
              grPsiGauntPsi(imesh, lm, 1:3) = matmul(conjg(Tmatrix(1:3, 1:3)), grPsiGauntPsiNat(imesh, lm, -1:1))
            end do ! imesh
          end do ! mqn_m
        end do ! oqn_l

        !NOte: there are differences close to the core which are seen when we do not use the r^2 multiplication any more.
        if (.false.) then
        do idir = 1, 3
          do oqn_l = 1, 1!atoms%lmax(itype)
            do mqn_m = -oqn_l, oqn_l
              lm = oqn_l * (oqn_l + 1) + 1 + mqn_m
              do imesh = 1, atoms%jri(itype)
                write(2277, '(3i8,2f20.8)') idir, mqn_m, imesh, 2 * grPsiGauntPsi(imesh, lm, idir) * atoms%rmsh(imesh, 1)
                write(2278, '(3i8,2f20.8)') idir, mqn_m, imesh, grRho0MT(imesh, lm, iatom, idir) * atoms%rmsh(imesh, 1)
              end do ! imesh
            end do ! mqn_m
          end do ! oqn_l
        end do ! idir
      end if

        do idir = 1, 3
          do oqn_l = 1, 1!atoms%lmax(itype)
            do mqn_m = -oqn_l, oqn_l
              lm = oqn_l * (oqn_l + 1) + 1 + mqn_m
              do imesh = 1, atoms%jri(itype)
                ! Neon takes 9e-7, Al 4e-6
                if (abs(2 * grPsiGauntPsi(imesh, lm, idir) * atoms%rmsh(imesh, 1) - grRho0MT(imesh, lm, iatom, idir) * atoms%rmsh(imesh, 1)) >= 1e-5) passed = .false.
              end do ! imesh
            end do ! mqn_m
          end do ! oqn_l
        end do ! idir

       ! do oqn_l = 0, 0!atoms%lmax(itype)
       !   do mqn_m = -oqn_l, oqn_l
       !     lm = oqn_l * (oqn_l + 1) + 1 + mqn_m
       !     do imesh = 1, atoms%jri(itype)
       !       write(2279, '(2i8,2f20.8)') mqn_m, imesh, psiGauntPsi(imesh, lm)
       !       write(2280, '(2i8,2f20.8)') mqn_m, imesh, rho0MtSpH(imesh, lm, iatom, 1)
       !     end do ! imesh
       !   end do ! mqn_m
       ! end do ! oqn_l

      end do ! ieqat
    end do ! itype

    if ( passed ) then
      write( logUnit, '(a)' ) '                                                           |__ passed!'
      write( logUnit, *)
    else
      write( logUnit, '(a)' ) '                                                           |__ failed!'
      write( logUnit, *)
      call JuDFT_warn('Testing the lm components of grRho0MT with Gaunt coefficient', calledby='testGrVarphi', hint='Debug!')
    end if

  end subroutine testGrVarphi

  !Compare of MT integral grRho grVeff with MT restricted <grPsi|grVeff|Psi>'
  subroutine TestDynMatPulInt( atoms, sym, dimens, lathar, stars, input, usdus, Veff0, kpts, cell, results, ngdp, logUnit, rho0IRst, rho0MTlh, Veff0Mtlh, gdp, clnu_atom, nmem_atom, mlh_atom, nRadFun, rbas1, rbas2, nv, ilst, GbasVec, nobd, z, kpq2kPrVec, El, kveclo, iloTable)

    use m_types, only : t_atoms, t_sym, t_dimension, t_sphhar, t_stars, t_potential, t_kpts, t_cell, t_input, t_usdus, t_results
    use m_jpConstants, only : iu, Tmatrix
    use mod_juPhonUtils, only : Derivative, CalcChannelsGrFlpNat, Calc2ArgCmplxIntMT, calcGrR2FinLH
    use m_jpSternhHF, only : calcMEPotIR
    use m_jpPotDensHelper, only : warpIRPot, convertStar2G
    use m_jpSetupDynMat, only : CalcGrVarphiHepsGrtVarphiElem
     
    use m_abcof3
    use m_jpSetupDynMatHelper, only : CalcHnGrV0Varphi, CalcFnsphVarphi
    use m_juDFT_NOstopNO, only : juDFT_warn

    use m_jpTestPotential, only : checkjuPhPots

    use m_intgr

    implicit none

    ! Type parameters
    type(t_atoms),                  intent(in) :: atoms
    type(t_sym),                    intent(in) :: sym
    type(t_dimension),              intent(in) :: dimens
    type(t_sphhar),                 intent(in) :: lathar
    type(t_stars),                  intent(in) :: stars
    type(t_potential),              intent(in) :: Veff0
    type(t_kpts),                   intent(in) :: kpts
    type(t_cell),                   intent(in) :: cell
    type(t_input),                  intent(in) :: input
    type(t_usdus),                  intent(in) :: usdus
    type(t_results),                intent(in) :: results

    ! Scalar parameters
    integer,                        intent(in) :: ngdp
    integer,                        intent(in) :: logUnit


    ! Array parameters
    complex,                        intent(in) :: rho0IRst(:, :)
    real,                           intent(in) :: rho0Mtlh(:, 0:, :, :)
    real,                           intent(in) :: Veff0Mtlh(:, 0:, :, :)
    integer,                        intent(in) :: gdp(:, :)
    complex,                        intent(in) :: clnu_atom(:, 0:, :)
    integer,                        intent(in) :: nmem_atom(0:, :)
    integer,                        intent(in) :: mlh_atom(:, 0:, :)
    integer,                        intent(in) :: nRadFun(0:, :)
    real,                           intent(in) :: rbas1(:,:,0:,:,:)
    real,                           intent(in) :: rbas2(:,:,0:,:,:)
    integer,                        intent(in) :: nv(:, :)
    integer,                        intent(in) :: ilst(:, :, :)
    integer,                        intent(in) :: GbasVec(:, :)
    integer,                        intent(in) :: nobd(:, :)
    MCOMPLEX,                       intent(in) :: z(:, :, :, :)
    integer,                        intent(in) :: kpq2kPrVec(:, :, :)
    real,                           intent(in) :: El(:, 0:, :, :)
    integer,                        intent(in) :: kveclo(:,:)
    integer,                        intent(in)    :: iloTable(:, 0:, :)

    ! Type variable
    type(od_inp)                                  :: odi
    type(od_sym)                                  :: ods

    ! Scalar variables
    integer                                    :: ikpt
    integer                                    :: iG
    integer                                    :: idirR
    integer                                    :: idirC
    integer                                    :: itype
    integer                                    :: ilh
    integer                                    :: imesh
    integer                                    :: oqn_l
    integer                                    :: mqn_m
    integer                                    :: lm
    integer                                    :: iqpt
    integer                                    :: iatom
    integer                                    :: ieqat
    integer                                    :: lm_pre
    integer                                    :: lmpMax
    integer                                    :: nRadFunMax
    integer                                    :: iradf
    integer                                    :: iBas
    integer                                    :: iband
    integer                                    :: nmat
    integer                                    :: ispin
    integer                                    :: mqn_m2prC
    integer                                    :: lmp
    integer                                    :: pMaxLocal
    integer                                    :: idir
    integer                                    :: mqn_m2PrR
    logical                                    :: passed = .true.


    integer         :: lmp1Pr
    integer         :: lm1Pr
    integer         :: oqn_l1Pr
    integer         :: mqn_m1Pr
    integer         :: mqn_m3Pr
    integer         :: iradf1Pr
    integer         :: ichanPr
    integer         :: oqn_l3Pr
    integer         :: lm3Pr
    real, allocatable            :: intgrdR(:)
    real, allocatable            :: intgrdI(:)
    real :: integralR
    real :: integralI

    ! Array variables
    integer,           allocatable             :: lmpT(:)
    complex,           allocatable             :: grVeff0IR(:, :)
    complex,           allocatable             :: grRho0IR(:, :)
    complex,           allocatable             :: vEff0Pw(:)
    complex,           allocatable             :: rho0IRpw(:)
    complex,           allocatable             :: w_grVeff0IR(:, :)
    real,              allocatable             :: r2Rho0MTlh(:, :, :, :)
    real,              allocatable             :: r2Veff0MTlh(:, :, :, :)
    complex,           allocatable             :: grRho0MT(:, :, :, :)
    complex,           allocatable             :: grVeff0MT(:, :, :, :)
    complex,           allocatable             :: rho1IRContainer(:, :, :)
    complex,           allocatable             :: w_vEff1IRContainer(:, :, :)
    complex,           allocatable             :: rho1MTContainer(:, :, :, :, :)
    complex,           allocatable             :: vEff1MTContainer(:, :, :, :, :)
    complex,           allocatable             :: r2GrRho0MTsh(:, :, :, :)
    complex,           allocatable             :: r2GrVeff0Mtsh(:, :, :, :)
    real,              allocatable             :: varphi1(:, :, :)
    real,              allocatable             :: varphi2(:, :, :)
    real,              allocatable             :: r2(:)
    real,              allocatable             :: delrVarphi1(:, :, :)
    real,              allocatable             :: delrVarphi2(:, :, :)
    integer,           allocatable             :: grVarphiChLout(:, :)
    integer,           allocatable             :: grVarphiChMout(:, :)
    real,              allocatable             :: grVarphiCh1(:, :, :, :)
    real,              allocatable             :: grVarphiCh2(:, :, :, :)
    complex,           allocatable             :: ikpGz0(:, :, :, :)
    complex,           allocatable             :: vEff0MtSpH(:, :)
    complex,           allocatable             :: vEff0NsphGrVarphi(:, :, :, :, :)
    complex,           allocatable             :: r2grVeff0SphVarphi(:, :, :, :, :)
    complex,           allocatable             :: grVarphiGrtVeff0SphVarphi(:, :, :, :, :)
    real,              allocatable             :: vEff0MtLhDummy(:, :, :)
    real,              allocatable             :: grVarphiGrtVarphi(:, :, :, :, :)
    complex,           allocatable             :: grVarphiHpreGrtVarphi(:, :, :, :, :)
    complex,           allocatable             :: hVarphi(:, :, :, :)
    complex,           allocatable             :: altIntgrlIRband(:, :)
    complex,           allocatable             :: a(:, :, :)
    complex,           allocatable             :: b(:, :, :)
    complex,           allocatable             :: bascof_lo(:, :, :, :, :)
    integer,           allocatable             :: ngoprI(:)
    complex,           allocatable             :: ab0cofScl(:)
    complex,           allocatable             :: grVarphiGrtVeff0PsiNat(:)
    complex                                    :: grPsiGrtVeff0PsiNat(-1:1, 3)
    complex                                    :: grPsiGrtVeff0PsiMt(3, 3)
    real                                       :: Gext(3)
    real                                       :: kExt(3)
    complex                                    :: altIntgrlIR(3, 3)
    complex                                    :: altIntgrl(3, 3)
    integer :: ii, jj
!!!!!!!!!!!!!!!!!!!!!!
complex, allocatable :: grVarphiGrtVeff0SphVarphiMeshNat(:, :, :, :, :, :)
complex, allocatable :: grPsiGrtVeff0PsiMtMesh(:, :, :)
complex, allocatable :: grPsiGrtVeff0PsiNatMesh(:, :, :)
complex, allocatable :: grVarphiGrtVeff0PsiNatMesh(:)
complex, allocatable :: integrandInt(:, :, :)
real, allocatable :: integralRMesh(:), integralIMesh(:)
complex, allocatable :: grVeff0MTCont(:, :, :, :)
complex :: ConvIntgrl(3, 3)
complex :: integral


      write(logUnit, '(a)')   'Compare of MT integral grRho grVeff with MT restricted <grPsi|grVeff|Psi>'
      write(logUnit, '(a)')     '-------------------------------------------------------------------------'

allocate(intgrdR(atoms%jmtd), intgrdI(atoms%jmtd))
allocate(integralRMesh(atoms%jmtd), integralIMesh(atoms%jmtd))
    ! Quantities for initialization
    allocate( lmpT(atoms%ntype) )
    lmpT(:) = 0
    do itype = 1, atoms%ntype
      lmpT(itype) = sum( [ ( (2 * oqn_l + 1)* nRadFun(oqn_l, itype), oqn_l = 0, atoms%lmax(itype) ) ] )
    end do ! itype
    lmpMax     = maxval( lmpT(:) )
allocate(grPsiGrtVeff0PsiMtMesh(atoms%jmtd, 3, 3), grPsiGrtVeff0PsiNatMesh(atoms%jmtd, -1:1, 3), grVarphiGrtVeff0PsiNatMesh(lmpMax))
intgrdR(:) = 0.
intgrdI(:) = 0.
integralRMesh(:) = 0.
integralIMesh(:) = 0.
grPsiGrtVeff0PsiMtMesh(:, :, :) = cmplx(0., 0.)
grPsiGrtVeff0PsiNatMesh(:, :, :) = cmplx(0., 0.)
grVarphiGrtVeff0PsiNatMesh(:) = cmplx(0., 0.)

    altIntgrl(:, :) = cmplx(0., 0.)

!!!!!!
    ! This test is only applicable for q = 0
    iqpt = 1

    ! Convert from star expansion coefficients to plane-wave expansion coefficients
    allocate( vEff0Pw(ngdp), rho0IRpw(ngdp) )
    vEff0Pw(:) = cmplx(0., 0.)
    rho0IRpw(:) = cmplx(0., 0.)
    call convertStar2G( Veff0%vpw_uw(:, 1), vEff0Pw(:), stars, ngdp, gdp )
    call convertStar2G( rho0IRst(:, 1), rho0IRpw(:), stars, ngdp, gdp )

    ! Perform analytical gradient of unperturbed density and effective potential in the interstitial region
    allocate( grVeff0IR(ngdp, 3), grRho0IR(ngdp, 3) )
    grVeff0IR(:, :) = cmplx(0., 0.)
    grRho0IR(:, :) = cmplx(0., 0.)
    do iG = 1, ngdp
      Gext(1:3) = matmul(cell%bmat(1:3, 1:3), gdp(1:3, iG))
      grVeff0IR(iG, 1:3) = iu * Gext(1:3) * vEff0Pw(iG)
      grRho0IR(iG, 1:3)  = iu * Gext(1:3) * rho0IRpw(iG)
    end do ! iG

    ! Warp interstitial second-order external potential
    allocate(w_grVeff0IR(ngdp, 3))
    w_grVeff0IR = cmplx(0.0, 0.0)
    do idirC = 1, 3
      call warpIRPot(stars, ngdp, idirC, gdp, grVeff0IR, w_grVeff0IR(:, idirC))
    end do ! idirC

    ! Gradient of density and effective potential in MT
    allocate( r2Rho0MTlh( atoms%jmtd, 0:lathar%nlhd, atoms%ntype, input%jspins) )
    allocate( r2Veff0MTlh( atoms%jmtd, 0:lathar%nlhd, atoms%ntype, 1) )
    allocate( grRho0MT( atoms%jmtd, ( atoms%lmaxd + 1 )**2, atoms%nat, 3) )
    allocate( grVeff0MT( atoms%jmtd, ( atoms%lmaxd + 1 )**2, atoms%nat, 3) )
    r2Rho0MTlh(:, :, :, :) = 0.
    r2Veff0MTlh(:, :, :, :) = 0.
    grRho0MT(:, :, :, :) = cmplx(0., 0.)
    grVeff0MT(:, :, :, :) = cmplx(0., 0.)

    ! Read in valence density form FLEUR
    rewind(1040)
    read(1040) r2Rho0MtLh
    ! todo factor of fpi in l = 0 component?

    do itype = 1, atoms%ntype
      do ilh = 0, lathar%nlhd
        do imesh = 1, atoms%jri(itype)
          r2Veff0MTlh(imesh, ilh, itype, 1) = vEff0MTlh(imesh, ilh, itype, 1) * atoms%rmsh(imesh, itype) * atoms%rmsh(imesh, itype)
        end do
      end do
    end do

    call calcGrR2FinLH( atoms, lathar, clnu_atom, nmem_atom, mlh_atom, r2Rho0MTlh(:, :, :, 1), r2GrRho0MTsh )
    call calcGrR2FinLH( atoms, lathar, clnu_atom, nmem_atom, mlh_atom, r2Veff0MTlh(:, :, :, 1), r2GrVeff0MTsh )

    do idirC = 1, 3
      iatom = 0
      do itype = 1, atoms%ntype
        do ieqat = 1, atoms%neq(itype)
          iatom = iatom + 1
          do oqn_l = 0, atoms%lmax(itype)
            lm_pre = oqn_l * (oqn_l + 1) + 1
            do mqn_m = -oqn_l, oqn_l
              lm = lm_pre + mqn_m
              do imesh = 1, atoms%jri(itype)
                ! For reasons of performance, this calculation is done within the idirC loop
                grRho0MT(imesh, lm, iatom, idirC) = r2GrRho0MTsh(imesh, lm, iatom, idirC) / atoms%rmsh(imesh, itype)**2
                grVeff0MT(imesh, lm, iatom, idirC) = r2GrVeff0MTsh(imesh, lm, iatom, idirC) / atoms%rmsh(imesh, itype)**2
              end do ! imesh
            end do ! mqn_m
          end do ! oqn_l
        end do ! ieqat
      end do ! itype
    end do ! idirC


    ! Fill container arrays with gradient coefficients so that they can be processed by EvelIntRho1Veff1
    allocate( rho1IRContainer(ngdp, 3, atoms%nat), w_vEff1IRContainer(ngdp, 3, atoms%nat),  &
      & rho1MTContainer( atoms%jmtd, (atoms%lmaxd + 1)**2, atoms%nat, 3, atoms%nat ), &
      & vEff1MTContainer( atoms%jmtd, (atoms%lmaxd + 1)**2, atoms%nat, 3, atoms%nat ) )

    rho1IRContainer(:, :, :) = cmplx(0., 0.)
    w_vEff1IRContainer(:, :, :) = cmplx(0., 0.)
    rho1MTContainer(:, :, :, :, :) = cmplx(0., 0.)
    vEff1MTContainer(:, :, :, :, :) = cmplx(0., 0.)

!    ! todo think about this for more than one atom


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!NOTE: NO IR FOR THIS TEST TO WORK. THIS WAS ADDED LATER TO TEST SOMETHONG DIFFERENT

    do idirC = 1, 3
      do iG = 1, ngdp
        rho1IRContainer(iG, idirC, 1) = cmplx(0., 0.)!-grRho0IR(iG, idirC)
        w_vEff1IRContainer(iG, idirC, 1) = cmplx(0., 0.)!-w_grVeff0IR(iG, idirC)
      end do ! iG
    end do ! idirC

    allocate(grVeff0MTCont(atoms%jmtd, (atoms%lmaxd + 1)**2, 3, atoms%nat))
    grVeff0MTCont(:, :, :, :) = cmplx(0., 0.)
    !todo does not work for more than one atom
    do idirC = 1, 3
      ! There should be an itype
      do oqn_l = 0, atoms%lmax(1)
        do mqn_m = -oqn_l, oqn_l
          lm = oqn_l * (oqn_l + 1) + 1 + mqn_m
          do imesh = 1, atoms%jri(1)
            rho1MTContainer(imesh, lm, 1, idirC, 1) = -grRho0MT(imesh, lm, 1, idirC)
            vEff1MTContainer(imesh, lm, 1, idirC, 1) = -grVeff0MT(imesh, lm, 1, idirC)
            grVeff0MTCont(imesh, lm, idirC, 1) = grVeff0MT(imesh, lm, 1, idirC)
          end do ! imesh
        end do ! mqn_m
      end do ! oqn_l
    end do ! idirC

    !todo continious version of grVeff instead of two gradients
    ! todo when calculating continious version is there the same integral used
!    call checkjuPhPots(30, atoms, ngdp, gdp, cell, lathar, sym, grVeff0IR, grVeff0MTCont, logUnit)
!    call checkjuPhPots(30, atoms, ngdp, gdp, cell, lathar, sym, w_grVeff0IR, grVeff0MTCont, logUnit)

    ! Calculate the integral in the conventional way
    ConvIntgrl(:, :) = cmplx(0., 0.)
    do idirC = 1, 3
      do idirR = 1, 3
        do oqn_l = 0, atoms%lmax(1)
          do mqn_m = -oqn_l, oqn_l
            lm = oqn_l * (oqn_l + 1) + 1 + mqn_m
            integral = 0
            call Calc2ArgCmplxIntMT( atoms, 1, rho1MTContainer(:, lm, 1, idirR, 1), vEff1MTContainer(:, lm, 1, idirC, 1), integral)
            convIntgrl(idirR, idirC) = convIntgrl(idirR, idirC) + integral
          end do ! mqn_m
        end do ! oqn_l
      end do ! idirR
    end do ! idirC

    if (.false.) then
      write(*, '(a)') 'Matrix of conventional integrals'
      write(*, '(3(2(es16.8,1x),3x))')  convIntgrl(1, :)
      write(*, '(3(2(es16.8,1x),3x))')  convIntgrl(2, :)
      write(*, '(3(2(es16.8,1x),3x))')  convIntgrl(3, :)
    end if


    ! Setup zBra and zKet for evaluation of respective benchmark braket

    ! Quantities for initialization
    !allocate( lmpT(atoms%ntype) )
    lmpT(:) = 0
    do itype = 1, atoms%ntype
      lmpT(itype) = sum( [ ( (2 * oqn_l + 1)* nRadFun(oqn_l, itype), oqn_l = 0, atoms%lmax(itype) ) ] )
    end do ! itype
    lmpMax     = maxval( lmpT(:) )
    nRadFunMax = maxval( nRadFun(:, :) )

    allocate( varphi1(atoms%jmtd, nRadFunMax, 0:atoms%lmaxd), varphi2(atoms%jmtd, nRadFunMax, 0:atoms%lmaxd) )
    allocate( grVarphiCh1(atoms%jmtd, 2, lmpMax, -1:1), grVarphiCh2(atoms%jmtd, 2, lmpMax, -1:1), &
            & grVarphiChLout(2, 0:atoms%lmaxd), grVarphiChMout(-atoms%lmaxd:atoms%lmaxd, -1:1) )
    allocate( delrVarphi1( atoms%jmtd, nRadFunMax, 0:atoms%lmaxd), delrVarphi2( atoms%jmtd, nRadFunMax, 0:atoms%lmaxd) )
    allocate( ikpGz0(dimens%nbasfcn, dimens%neigd, 3, kpts%nkpt) )
    allocate( vEff0MtSpH( atoms%jmtd, 0:dimens%lmd), vEff0NsphGrVarphi(2, atoms%jmtd, 0:(atoms%lmaxd + 3)**2 - 1, lmpMax, -1:1), &
            & r2grVeff0SphVarphi(2, atoms%jmtd, 0:(atoms%lmaxd + 3)**2 - 1, lmpMax, 3) )
    allocate( grVarphiGrtVeff0SphVarphi(lmpMax, lmpMax, -1:1, 1:3, atoms%nat), vEff0MtLhDummy( atoms%jmtd, 0:lathar%nlhd, atoms%ntype ) )
    allocate( r2(atoms%jmtd) )
    allocate( grVarphiGrtVarphi(lmpMax, lmpMax, -1:1, -1:1, atoms%ntype), &
            & grVarphiHpreGrtVarphi(lmpMax, lmpMax, -1:1, -1:1, atoms%nat) )
    allocate( altIntgrlIRband(maxval(nobd), maxval(nobd)) )
    allocate( ab0cofScl(lmpMax) )
    allocate(grVarphiGrtVeff0PsiNat(lmpMax))

    grVarphiCh1(:, :, :, :) = 0.
    grVarphiCh2(:, :, :, :) = 0.
    grVarphiChLout(:, :) = 0
    grVarphiChMout(:, :) = 0
    ikpGz0(:, :, :, :) = cmplx(0., 0.)
    vEff0MtSpH = cmplx(0., 0.)
    vEff0NsphGrVarphi(:, :, :, :, :) = cmplx(0., 0.)
    r2grVeff0SphVarphi(:, :, :, :, :) = cmplx(0., 0.)
    grVarphiGrtVeff0SphVarphi(:, :, :, :, :) = cmplx(0., 0.)
    vEff0MtLhDummy(:, :, :) = cmplx(0., 0.)
    r2(:) = 0.
    grVarphiGrtVarphi(:, :, :, :, :) = 0.
    grVarphiHpreGrtVarphi(:, :, :, :, :) = cmplx(0., 0.)
    altIntgrlIRband(:, :) = cmplx(0., 0.)
    grVarphiGrtVeff0PsiNat(:) = cmplx(0., 0.)
    ab0cofScl(:) = cmplx(0., 0.)

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!
allocate(grVarphiGrtVeff0SphVarphiMeshNat(atoms%jmtd, lmpMax, lmpMax, -1:1, 1:3, atoms%nat))
grVarphiGrtVeff0SphVarphiMeshNat = 0.
!!!!!!!!!!!!!!!!!!!!!!!
    iatom = 0
    do itype = 1, atoms%ntype

      ! These arrays are set to zero here, because 
      varphi1(:, :, :) = 0.
      varphi2(:, :, :) = 0.
      delrVarphi1(:, :, :) = 0.
      delrVarphi2(:, :, :) = 0.
      do oqn_l = 0, atoms%lmax(itype)
        do iradf = 1, nRadFun(oqn_l, itype)
          do imesh = 1, atoms%jri(itype)
            ! In Fleur the radial solutions u_lp are multiplied by a factor r by default to avoid an additional multiplication of the
            ! Jacobi determinant r^2 in radial integrals given spherical coordinates.
            varphi1(imesh, iradf, oqn_l) = rbas1(imesh, iradf, oqn_l, itype, 1) / atoms%rmsh(imesh, itype)
            varphi2(imesh, iradf, oqn_l) = rbas2(imesh, iradf, oqn_l, itype, 1) / atoms%rmsh(imesh, itype)
          end do ! imesh
          ! Precalculate partial derivatives of varphis in r-direction since it is needed twice
          call Derivative( varphi1(1:atoms%jri(itype), iradf, oqn_l), itype, atoms, delrVarphi1(1:atoms%jri(itype), iradf, oqn_l) )
          call Derivative( varphi2(1:atoms%jri(itype), iradf, oqn_l), itype, atoms, delrVarphi2(1:atoms%jri(itype), iradf, oqn_l) )
        end do ! iradf
      end do ! oqn_l

      ! Calculate the application of the gradient and the gradient's dyadic product onto the MT basis functions (matching coefficients
      ! have no spatial dependence) and determing its scattering channels.
      grVarphiChLout(:, :) = 0
      grVarphiChMout(:, :) = 0
      grVarphiCh1(:, :, :, :) = 0.
      grVarphiCh2(:, :, :, :) = 0.
      call CalcChannelsGrFlpNat( atoms, itype, nRadFun, varphi1, varphi2, delrVarphi1, delrVarphi2, grVarphiChLout, grVarphiChMout, &
                                                                                                        & grVarphiCh1, grVarphiCh2 )

      do ieqat = 1, atoms%neq(itype)
        iatom = iatom + 1
        r2grVeff0SphVarphi(:, :, :, :, :) = cmplx(0., 0.)
        do idirC = 1, 3
          call CalcFnsphVarphi( atoms, itype, 0, nRadFun, varphi1, varphi2, r2GrVeff0MTSh(:, :, iatom, idirC), &
                                                                                          & r2grVeff0SphVarphi(:, :, :, :, idirC) )
        end do

        ! Calculate the integral
        do mqn_m2PrC = -1, 1

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!    do mqn_m2PrR = -1, 1
!      lmp = 0
!      lm = -1
!      do oqn_l = 0, atoms%lmax(itype)
!        do mqn_m = -oqn_l, oqn_l
!          lm = lm + 1
!          do iradf = 1, nRadFun(oqn_l, itype)
!            lmp = lmp + 1
!            lmp1Pr = 0
!            lm1Pr = -1
!            do oqn_l1Pr = 0, atoms%lmax(itype)
!              do mqn_m1Pr = -oqn_l1Pr, oqn_l1Pr
!                mqn_m3Pr = grVarphiChMout(mqn_m1Pr, mqn_m2PrR)
!                lm1Pr = lm1Pr  + 1
!                do iradf1Pr = 1, nRadFun(oqn_l1Pr, itype)
!                  lmp1Pr = lmp1Pr + 1
!                  ! overlap of grVarphi grVarphi
!                  do ichanPr = 1, 2
!                    oqn_l3Pr = grVarphiChLout(ichanPr, oqn_l1Pr)
!                    ! We have to catch channels where l < 0 and where m"' does not contribute, lmax + 1 is allowed due to the gradient
!                    if ( ( abs(mqn_m3Pr) > oqn_l3Pr )  .or. ( oqn_l3Pr < 0 )  )  cycle
!                    lm3Pr = oqn_l3Pr * (oqn_l3Pr + 1) + mqn_m3Pr
!
!                   ! intgrdR(:) = 0.
!                   ! intgrdI(:) = 0.
!                    do imesh = 1, atoms%jri(itype)
!      !                intgrdR(imesh) = real ( &
!      !                  &   (grVarPhiCh1(imesh, ichanPr, lmp1Pr, mqn_m2PrR) * r2grVeff0SphVarphi(1, imesh, lm3Pr, lmp, mqn_m2PrC + 2) &
!      !                  &  + grVarPhiCh2(imesh, ichanPr, lmp1Pr, mqn_m2PrR) * r2grVeff0SphVarphi(2, imesh, lm3Pr, lmp, mqn_m2PrC + 2) ) )
!      !                intgrdI(imesh) = aimag( &
!      !                  !todo why without iatom
!      !                  &   (grVarPhiCh1(imesh, ichanPr, lmp1Pr, mqn_m2PrR) * r2grVeff0SphVarphi(1, imesh, lm3Pr, lmp, mqn_m2PrC + 2) &
!      !                  &  + grVarPhiCh2(imesh, ichanPr, lmp1Pr, mqn_m2PrR) * r2grVeff0SphVarphi(2, imesh, lm3Pr, lmp, mqn_m2PrC + 2) ) )
!                      grVarphiGrtVeff0SphVarphiMeshNat(imesh, lmp1Pr, lmp, mqn_m2PrR, mqn_m2PrC + 2, iatom) =  grVarphiGrtVeff0SphVarphiMeshNat(imesh, lmp1Pr, lmp, mqn_m2PrR, mqn_m2PrC + 2, iatom) + ( &
!                        &   (grVarPhiCh1(imesh, ichanPr, lmp1Pr, mqn_m2PrR) * r2grVeff0SphVarphi(1, imesh, lm3Pr, lmp, mqn_m2PrC + 2) &
!                        &  + grVarPhiCh2(imesh, ichanPr, lmp1Pr, mqn_m2PrR) * r2grVeff0SphVarphi(2, imesh, lm3Pr, lmp, mqn_m2PrC + 2) ) )
!                    end do ! imesh
!                   ! call intgr3(intgrdR(:), atoms%rmsh(1, itype), atoms%dx(itype), atoms%jri(itype), integralR)
!                   ! call intgr3(intgrdI(:), atoms%rmsh(1, itype), atoms%dx(itype), atoms%jri(itype), integralI)
!                   ! ! < grVarphi| grTveff0Sph |varphi >
!                   ! ! idir = mqn_m2PrC + 2 due to performance reasons
!                   ! grVarphiGrtVeff0SphVarphi(lmp1Pr, lmp, mqn_m2PrR, mqn_m2PrC + 2, iatom) = &
!                   !          & grVarphiGrtVeff0SphVarphi(lmp1Pr, lmp, mqn_m2PrR, mqn_m2PrC + 2, iatom) + cmplx(integralR, integralI)
!                  end do ! ichanPr
!                end do ! iradf1Pr
!              end do ! mqn_m1Pr
!            end do ! oqn_l1Pr
!          end do ! iradf
!        end do ! mqn_m
!      end do ! oqn_l
!    end do ! mqn_m2PrR
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

          call CalcGrVarphiHepsGrtVarphiElem( atoms, itype, iatom, mqn_m2PrC, lmpMax, nRadFun, grVarphiChLout, grVarphiChMout, r2, &
            & grVarPhiCh1, grVarPhiCh2, El, lmpT, vEff0NsphGrVarphi, r2grVeff0SphVarphi, grVarphiGrtVarphi, grVarphiHpreGrtVarphi, &
            & grVarphiGrtVeff0SphVarphi)
        end do ! mqn_m2PrC
      end do ! ieqat
    end do ! itype

!    do itype = 1, atoms%ntype
!      do mqn_m2PrR  = -1, 1
!        do mqn_m2PrC = -1, 1
!          do ii = 1, lmpMax
!            do jj = 1, lmpMax
!              write(2292, '(4i8,f20.8)') ii, jj, mqn_m2PrR, mqn_m2PrC, grVarphiGrtVarphi(ii, jj, mqn_m2PrR, mqn_m2PrC, 1)
!              write(2293, '(4i8,2f20.8)') ii, jj, mqn_m2PrR, mqn_m2PrC, grVarphiHpreGrtVarphi(ii, jj, mqn_m2PrR, mqn_m2PrC, 1)
!              write(2294, '(4i8,2f20.8)') ii, jj, mqn_m2PrR, mqn_m2PrC, grVarphiGrtVeff0SphVarphi(ii, jj, mqn_m2PrR, mqn_m2PrC, 1)
!            end do ! jj
!          end do ! ii
!        end do ! mqn_m2PrC
!      end do ! mqn_m2PrR
!    end do ! itype
!
!
    ispin = 1
    allocate( a( dimens%nvd, 0:dimens%lmd, atoms%nat), b(dimens%nvd, 0:dimens%lmd, atoms%nat), &
      & bascof_lo(3, -atoms%llod:atoms%llod, 4 * atoms%llod + 2, atoms%nlod, atoms%nat) )
    allocate( ngoprI(atoms%nat) )
    a(:, :, :)               = cmplx(0., 0.)
    b(:, :, :)               = cmplx(0., 0.)
    bascof_lo(:, :, :, :, :) = cmplx(0., 0.)
    ! We do not want the local coordinate systems to be rotated for non-representative atoms constructing the matching coefficients.
    ngoprI(:) = 1
    iatom = 0
    do itype = 1, atoms%ntype
      do ieqat = 1, atoms%neq(itype)
        iatom = iatom + 1
        grPsiGrtVeff0PsiMt(:, :) = cmplx(0., 0.)
!!        altIntgrlIR(:, :) = cmplx(0., 0.)
        do ikpt = 1, kpts%nkpt
!!write(*, *) 'ikpt = ', ikpt
          nmat = nv(1, ikpt) + atoms%nlotot
          a(:, :, :) = cmplx(0.0, 0.0)
          b(:, :, :) = cmplx(0.0, 0.0)
          bascof_lo(:, :, :, :, :) = cmplx(0.0, 0.0)
          call abcof3( atoms%lmaxd, atoms%ntype, atoms%nat, sym%nop, dimens%nvd, dimens%jspd, 1, dimens%lmd, dimens%nbasfcn, &
            & atoms%llod, atoms%nlod, atoms%nlotot, sym%invtab, atoms%ntype, sym%mrot, ngoprI, atoms%taual, atoms%neq, atoms%lmax, &
            & atoms%rmt, cell%omtil, cell%bmat, cell%bbmat, kpts%bk(:, ikpt), gbasVec(1, ilst(:nv(1, ikpt), ikpt, 1)), &
            & gbasVec(2, ilst(:nv(1, ikpt), ikpt, 1)), gbasVec(3, ilst(:nv(1, ikpt), ikpt, 1)), nv(:, ikpt), nmat, &
            & usdus%us, usdus%dus, usdus%uds, usdus%duds, usdus%ddn, atoms%invsat, sym%invsatnr, usdus%ulos, usdus%uulon, usdus%dulon, &
            & usdus%dulos, atoms%llo, atoms%nlo, atoms%l_dulo, atoms%lapw_l, kveclo(:,ikpt), odi, ods, a, b, bascof_lo )

          grPsiGrtVeff0PsiNat(:, :) = cmplx(0., 0.)
          !!!!!!
!          grPsiGrtVeff0PsiNatMesh(:, :, :) = cmplx(0., 0.)
          !!!!!!!
          do iband = 1, nobd(ikpt, 1)
!write(*, *) 'iband = ', iband
            ab0cofScl(:) = cmplx(0., 0.)
            lmp = 0
            lm  = 0
            do oqn_l = 0, atoms%lmax(itype)
              do mqn_m = - oqn_l, oqn_l
                pMaxLocal = nRadFun(oqn_l, itype)
                ! p = 1
                ab0cofScl(lmp + 1) = iu**oqn_l * dot_product( conjg(z(:nv(1, ikpt), iband, ikpt, 1)), a(:nv(1, ikpt), lm, iatom) )
                ! p = 2
                ab0cofScl(lmp + 2) = iu**oqn_l * dot_product( conjg(z(:nv(1, ikpt), iband, ikpt, 1)), b(:nv(1, ikpt), lm, iatom) )
                ! Add LO contributions
                do iradf = 3, pMaxLocal
                  ! p = 1
                  ab0cofScl(lmp + 1) = ab0cofScl(lmp + 1) + &
                    & iu**oqn_l * dot_product( conjg(z(nv(1, ikpt) + 1:nv(1, ikpt) + atoms%nlotot, iband, ikpt, 1)), &
                                                         & bascof_lo(1, mqn_m, :atoms%nlotot, iloTable(iradf, oqn_l, itype), iatom) )
                  ! p = 2
                  ab0cofScl(lmp + 2) = ab0cofScl(lmp + 2) + &
                    & iu**oqn_l * dot_product( conjg(z(nv(1, ikpt) + 1:nv(1, ikpt) + atoms%nlotot, iband, ikpt, 1)), &
                                                         & bascof_lo(2, mqn_m, :atoms%nlotot, iloTable(iradf, oqn_l, itype), iatom) )
                  ! 2 < p < LOs for that l and that atom type
                  ab0cofScl(lmp + iradf) = iu**oqn_l * dot_product( conjg(z(nv(1, ikpt) + 1:nv(1, ikpt) + atoms%nlotot, iband, ikpt, 1)), &
                                                         & bascof_lo(3, mqn_m, :atoms%nlotot, iloTable(iradf, oqn_l, itype), iatom) )
                end do ! iradf


                ! This is a precalculation of the 1st and 3rd line of A.50. Actually the resulting lmp should have been primed or at
                ! least the p as lm1Pr = lm. But for sake of performance we place it here.
                ! sake of performance
                lm = lm + 1
                lmp = lmp + pMaxLocal
              end do ! mqn_m
            end do !oqn_l
            do idir = 1, 3
              do mqn_m2PrR = -1, 1
                grVarphiGrtVeff0PsiNat(1:lmpT(itype)) = cmplx(0., 0.)
                grVarphiGrtVeff0PsiNat(1:lmpT(itype)) = matmul(grVarphiGrtVeff0SphVarphi(1:lmpT(itype), 1:lmpT(itype), mqn_m2PrR, idir, iatom), ab0cofScl(1:lmpT(itype)))
                grPsiGrtVeff0PsiNat(mqn_m2PrR, idir) = grPsiGrtVeff0PsiNat(mqn_m2PrR, idir) + 2 * results%w_iks(iband, ikpt, 1) * dot_product(ab0cofScl(1:lmpT(itype)), grVarphiGrtVeff0PsiNat(1:lmpT(itype)))


                !!!!!!!!!!!!!!!!!!!!!!!!

!                do imesh = 1, atoms%jri(itype)
!                  grVarphiGrtVeff0PsiNatMesh(1:lmpT(itype)) = cmplx(0., 0.)
!                  grVarphiGrtVeff0PsiNatMesh(1:lmpT(itype)) = matmul(grVarphiGrtVeff0SphVarphiMeshNat(imesh, 1:lmpT(itype), 1:lmpT(itype), mqn_m2PrR, idir, iatom), ab0cofScl(1:lmpT(itype)))
!                  grPsiGrtVeff0PsiNatMesh(imesh, mqn_m2PrR, idir) = grPsiGrtVeff0PsiNatMesh(imesh, mqn_m2PrR, idir) + 2 * results%w_iks(iband, ikpt, 1) * dot_product(ab0cofScl(1:lmpT(itype)), grVarphiGrtVeff0PsiNatMesh(1:lmpT(itype)))
!                end do ! imesh

                !!!!!!!!
              end do ! mqn_m2PrR
            end do ! idir
          end do ! iband
          grPsiGrtVeff0PsiMt(1:3, 1:3) = grPsiGrtVeff0PsiMt(1:3, 1:3) + matmul(conjg(Tmatrix(1:3, 1:3)), grPsiGrtVeff0PsiNat(-1:1, 1:3))
          !!!!!!!!
!          do imesh = 1, atoms%jri(itype)
!            grPsiGrtVeff0PsiMtMesh(imesh, 1:3, 1:3) = grPsiGrtVeff0PsiMtMesh(imesh, 1:3, 1:3) + matmul(conjg(Tmatrix(1:3, 1:3)), grPsiGrtVeff0PsiNatMesh(imesh, -1:1, 1:3))
!          end do ! imesh
          !!!!!!!

          ikpGz0 = cmplx(0.0, 0.0)
          !kExt(1:3) = matmul(cell%bmat, kpts%bk(1:3, ikpt))
          !gExt(:) = 0.
          !do idirR = 1, 3
          !  do iBas = 1, nv(1, ikpt) !+ atoms%nlotot
          !    gExt(1:3) = matmul( cell%bmat(1:3, 1:3), gBasVec(1:3, ilst(iBas, ikpt, 1)))
          !    do iband = 1, nobd(ikpt, 1)
          !      ikpGz0(iBas, iband, idirR, ikpt) = iu * ( kExt(idirR) + gExt(idirR) ) * z(iBas, iband, ikpt, 1)
          !    end do ! iband
          !  end do ! iBas
          !end do ! idirR

          !nmat = nv(1, ikpt) + atoms%nlotot
          !do idirC = 1, 3
          !  do idirR = 1, 3
          !    altIntgrlIRband(:, :) = cmplx(0., 0.)
          !    call calcMEPotIR( input, stars, dimens, GbasVec(:, ilst(:nv(input%jspins, ikpt), ikpt, input%jspins)),                   &
          !      & GbasVec(:, ilst(:nv(input%jspins, ikpt), ikpt, input%jspins)), nv, w_vEff1IRContainer(:, idirC, 1),                  &
          !      & ikpGz0(:, :, idirR, ikpt), z(:, :, ikpt, 1), gdp, nmat, nobd(ikpt, 1), nobd(ikpt, 1), ispin, ikpt, iqpt, ikpt, ngdp, &
          !      & altIntgrlIRband(:, :), kpq2kPrVec )
          !    do iband = 1, nobd(ikpt, 1)
          !      altIntgrlIR(idirR, idirC) = altIntgrlIR(idirR, idirC) + 2 * results%w_iks(iband, ikpt, 1) * altIntgrlIRband(iband, iband)
          !    end do ! iband
          !  end do ! idirR
          !end do ! idirC
!
!
!          !do idirC = 1, 3
!          !  do idirR = 1, 3
!          !    do iband = 1, nobd(ikpt, 1)
!          !      altIntgrlIR = altIntgrlIR + altIntgrlIRband(iband, iband, idirR, idirC)
!          !    end do ! iband
!          !    altIntgrl(idirR, idirC) = altIntgrl(idirR, idirC) + altIntgrlIR(idirR, idirC) + altIntgrlMT(idirR, idirC)
!          !  end do ! idirR
!          !end do ! idirC
!

       end do !ikptR
        !altIntgrl(1:3, 1:3) = 2 * (altIntgrlIR(1:3, 1:3) + grPsiGrtVeff0PsiMt(1:3, 1:3))
        altIntgrl(1:3, 1:3) = 2 * (grPsiGrtVeff0PsiMt(1:3, 1:3))

!        allocate(integrandInt(atoms%jmtd, 3, 3))
!        integrandInt(:, :, :) = cmplx(0., 0.)
!        rewind(2284)
!        read(2284) integrandInt
!        do idirC = 1, 3
!          do idirR = 1, 3
!            do imesh = 1, atoms%jri(itype)
!!              intgrdR(imesh) = real ( integrandInt(imesh, idirR, idirC) - 2 * grPsiGrtVeff0PsiMtMesh(imesh, idirR, idirC) )
!!              intgrdI(imesh) = aimag( integrandInt(imesh, idirR, idirC) - 2 * grPsiGrtVeff0PsiMtMesh(imesh, idirR, idirC) )
!              intgrdR(imesh) = real ( 2 * grPsiGrtVeff0PsiMtMesh(imesh, idirR, idirC) )
!              intgrdI(imesh) = aimag( 2 * grPsiGrtVeff0PsiMtMesh(imesh, idirR, idirC) )
!            end do ! imesh
!            call intgr2(intgrdR, atoms%rmsh(1,itype), atoms%dx(itype), atoms%jri(itype), integralRMesh)
!            call intgr2(intgrdI, atoms%rmsh(1,itype), atoms%dx(itype), atoms%jri(itype), integralIMesh)
!            do imesh = 1, atoms%jri(itype)
!              write(2285, '(3i8,2f22.8)') idirC, idirR, imesh, integralRMesh(imesh), integralIMesh(imesh)
!            end do ! imesh
!          end do ! idirR
!        end do ! idirC

!        do idirC = 1, 3
!          do idirR = 1, 3
!            do imesh = 1, atoms%jri(itype)
!              write(2283, '(3i8,2f22.8)') idirC, idirR, imesh, 2 * grPsiGrtVeff0PsiMtMesh(imesh, idirR, idirC)
!            end do ! imesh
!          end do ! idirR
!        end do ! idirC

        if (.false.) then
          write(*, '(a)') 'Matrix of brakets'
          write(*, '(3(2(es16.8,1x),3x))') altIntgrl(1, :)
          write(*, '(3(2(es16.8,1x),3x))') altIntgrl(2, :)
          write(*, '(3(2(es16.8,1x),3x))') altIntgrl(3, :)
        end if
        if (any(abs(altIntgrl(:, :) - convIntgrl(:, :)) > 9e-6)) passed = .false.
      end do ! ieqat
    end do ! itype
    deallocate(a, b, bascof_lo)


    if ( passed ) then
      write( logUnit, '(a)' ) '                                                                        |__ passed!'
    else
      write( logUnit, '(a)' ) '                                                                        |__ failed!'
      call JuDFT_warn('Compare of MT integral grRho grVeff with MT restricted <grPsi|grVeff|Psi>', calledby='TestDynMatPulInt', hint='Check logfile for more information!')
    end if

  end subroutine TestDynMatPulInt

  ! Check cancelling within 3rd braket of Pulay contribution to the dynamical matrix when  setting a special z1
  subroutine Test3rdBraKetHeps( atoms, lathar, kpts, qpts, sym, dimens, cell, usdus, stars, Veff0, results, logUnit, nRadFun, nv, GbasVec, ilst, kveclo, mapKpq2K, &
                               & nobd, z, iloTable, eig, kpq2kPrVec, El, vEff0MtLh, nmem_atom, mlh_atom, clnu_atom, rbas1, rbas2 )

    use m_types
    use m_jpSetupDynMat, only : Add1stOrdWfPulayBraKets2DynMat, CalcScalBasfMatElems
    use m_jpSetupDynMatHelper, only : CalcHnGrV0Varphi
    use m_juDFT_NOstopNO, only : juDFT_warn
    use m_jpSetupDynMatSF, only : CalcHGrVarphi

    implicit none

    ! Type parameters
    type(t_atoms),              intent(in)  :: atoms
    type(t_sphhar),             intent(in)  :: lathar
    type(t_kpts),               intent(in)  :: kpts
    type(t_kpts),               intent(in)  :: qpts
    type(t_sym),                intent(in)  :: sym
    type(t_dimension),          intent(in)  :: dimens
    type(t_cell),               intent(in)  :: cell
    type(t_usdus),              intent(in)  :: usdus
    type(t_stars),              intent(in)  :: stars
    type(t_potential),          intent(in)  :: Veff0
    type(t_results),            intent(in)  :: results

    ! Scalar parameter
    integer,                    intent(in)  :: logUnit

    ! Array parameters
    integer,                    intent(in)  :: nRadFun(0:, :)
    integer,                    intent(in)  :: nv(:, :)
    integer,                    intent(in)  :: GbasVec(:, :)
    integer,                    intent(in)  :: ilst(:, :, :)
    integer,                    intent(in)  :: kveclo(:,:)
    integer,                    intent(in)  :: mapKpq2K(:, :)
    integer,                    intent(in)  :: nobd(:, :)
    MCOMPLEX,                   intent(in)  :: z(:,:,:,:)
    integer,                    intent(in)  :: iloTable(:, 0:, :)
    real,                       intent(in)  :: eig(:, :, :)
    integer,                    intent(in)  :: kpq2kPrVec(:, :, :)
    real,                       intent(in)  :: El(:, 0:, :, :)
    real,                       intent(in)  :: vEff0MtLh(:, 0:, :)
    integer,                    intent(in)  :: nmem_atom(0:, :)
    integer,                    intent(in)  :: mlh_atom(:, 0:, :)
    complex,                    intent(in)  :: clnu_atom(:, 0:, :)
    real,                       intent(in) :: rbas1(:, :, 0:, :, :)
    real,                       intent(in) :: rbas2(:, :, 0:, :, :)

    ! Scalar variable
    integer                                 :: ikpt
    integer                                 :: iqpt
    integer                                 :: lmpMax
    integer                                 :: itype
    integer                                 :: ieqat
    integer                                 :: iatom
    integer                                 :: iradf
    integer                                 :: imesh
    integer                                 :: oqn_l
    integer                                 :: ptsym
    integer                                 :: ilh
    integer                                 :: nRadFunMax
    integer                                 :: lm_pre
    integer                                 :: mqn_m
    integer                                 :: lm
    integer                                 :: iBas
    integer                                 :: iband
    integer                                 :: idir
    integer                                 :: imem
    integer                                 :: ii
    integer                                 :: jj
    integer                                 :: mqn_m2PrC
    logical                                 :: eps1DynMatPulTestSw
    logical                                 :: testCompTerm3rdBraKetsVarBra
    logical                                 :: testCompTerm3rdBraKetsVarKet
    logical                                 :: dynMatPu3rdBraKetHepsSw
    logical                                 :: passed = .true.

    ! Array variables
    integer,       allocatable              :: lmpT(:)
    real,          allocatable              :: grVarphiVarphi(:, :, :, :)
    real,          allocatable              :: varphiGrVarphi(:, :, :, :)
    complex,       allocatable              :: z1nG(:, :, :, :)
    complex,       allocatable              :: grVarphiHVarphi(:, :, :, :)
    real,          allocatable              :: varphiVarphi(:, :, :)
    complex,       allocatable              :: varphiHvarphi(:, :, :)
    integer,       allocatable              :: grVarphiChLout(:, :)
    integer,       allocatable              :: grVarphiChMout(:, :)
    real,          allocatable              :: grVarphiCh1(:, :, :, :)
    real,          allocatable              :: grVarphiCh2(:, :, :, :)
    real,          allocatable              :: r2(:)
    real,          allocatable              :: varphi1(:, :, :)
    real,          allocatable              :: varphi2(:, :, :)
    complex,       allocatable              :: vEff0MtSpH(:, :)
    complex,       allocatable              :: hVarphi(:, :, :, :)
    complex,       allocatable              :: vEff0NsphGrVarphi(:, :, :, :, :)
    complex,       allocatable              :: r2grVeff0SphVarphi(:, :, :, :, :)
    complex,       allocatable              :: r2grVeff0SphVarphiDummy(:, :, :, :, :)
    complex,       allocatable              :: hGrVarphi(:, :, :, :, :)
    real,          allocatable              :: varphiVarphiDummy(:, :, :)
    complex,       allocatable              :: varphiGrVeff0SphVarphi(:, :, :, :)
    complex,       allocatable              :: varphiHGrvarphi(:, :, :, :)
    complex                                 :: pulayCntrb(3, 3)

    write( logUnit, '(a)' ) 'Check cancelling within 3rd braket in Pulay contribution to the dynamical matrix.'
    write( logUnit, '(a)' ) '---------------------------------------------------------------------------------'

    iqpt = 1

    ! Quantities for initialization
    allocate( lmpT(atoms%ntype) )
    lmpT(:) = 0
    do itype = 1, atoms%ntype
      lmpT(itype) = sum( [ ( (2 * oqn_l + 1)* nRadFun(oqn_l, itype), oqn_l = 0, atoms%lmax(itype) ) ] )
    end do ! itype
    lmpMax     = maxval( lmpT(:) )
    nRadFunMax = maxval( nRadFun(:, :) )

    allocate( grVarphiVarphi(lmpMax, lmpMax, -1:1, atoms%ntype), grVarphiHVarphi(lmpMax, lmpMax, -1:1, atoms%nat ) )
    allocate( varphiVarphi(lmpMax, lmpMax, atoms%ntype), varphiHvarphi(lmpMax, lmpMax, atoms%nat))
    allocate( varphi1(atoms%jmtd, nRadFunMax, 0:atoms%lmaxd), varphi2(atoms%jmtd, nRadFunMax, 0:atoms%lmaxd) )
    allocate( z1nG(dimens%nbasfcn, 3, atoms%nat, maxval(nobd(:, :))) )
    allocate( r2(atoms%jmtd), vEff0MtSpH( atoms%jmtd, 0:dimens%lmd) )
    allocate( grVarphiCh1(atoms%jmtd, 2, lmpMax, -1:1), grVarphiCh2(atoms%jmtd, 2, lmpMax, -1:1), &
            & grVarphiChLout(2, 0:atoms%lmaxd), grVarphiChMout(-atoms%lmaxd:atoms%lmaxd, -1:1) )
    allocate( hVarphi(2, atoms%jmtd, 0:(atoms%lmaxd + 3)**2 - 1, lmpMax)  )
    allocate( vEff0NsphGrVarphi(2, atoms%jmtd, 0:(atoms%lmaxd + 3)**2 - 1, lmpMax, -1:1), &
            & r2grVeff0SphVarphi(2, atoms%jmtd, 0:(atoms%lmaxd + 3)**2 - 1, lmpMax, 3) )
    allocate( r2grVeff0SphVarphiDummy(2, atoms%jmtd, 0:(atoms%lmaxd + 3)**2 - 1, lmpMax, 3) )
    allocate( hGrVarphi(2, atoms%jmtd, (atoms%lmaxd + 1)**2, lmpMax, -1:1))
    allocate( varphiVarphiDummy(lmpMax, lmpMax, atoms%ntype) )
    allocate( varphiGrVeff0SphVarphi(lmpMax, lmpMax, atoms%nat, 3) )
    allocate( varphiHGrvarphi(lmpMax, lmpMax, atoms%nat, -1:1) )
    allocate( varphiGrvarphi(lmpMax, lmpMax, -1:1, atoms%ntype) )

    grVarphiVarphi(:, :, :, :) = 0.
    varphiGrvarphi(:, :, :, :) = 0.
    grVarphiHvarphi(:, :, :, :) = cmplx(0., 0.)
    varphiVarphi(:, :, :) = cmplx(0., 0.)
    varphiHvarphi(:, :, :) = cmplx(0., 0.)
    hVarphi(:, :, :, :) = cmplx(0., 0.)
    grVarphiChLout = 0
    grVarphiChMout = 0
    grVarphiCh1(:, :, :, :) = cmplx(0., 0.)
    grVarphiCh2(:, :, :, :) = cmplx(0., 0.)
    vEff0MtSpH(:, :) = cmplx(0., 0.)
    vEff0NsphGrVarphi(:, :, :, :, :) = cmplx(0., 0.)
    r2grVeff0SphVarphi(:, :, :, :, :) = cmplx(0., 0.)
    r2grVeff0SphVarphiDummy(:, :, :, :, :) = cmplx(0., 0.)
    hGrVarphi(:, :, :, :, :) = cmplx(0., 0.)
    varphiGrVeff0SphVarphi(:, :, :, :) = cmplx(0., 0.)
    varphiVarphiDummy(:, :, :) = 0.
    varphiHGrvarphi(:, :, :, :) = cmplx(0., 0.)


    iatom = 0
    do itype = 1, atoms%ntype
      do ieqat = 1, atoms%neq(itype)
        iatom = iatom + 1

        varphi1(:, :, :) = 0.
        varphi2(:, :, :) = 0.

        do oqn_l = 0, atoms%lmax(itype)
          do iradf = 1, nRadFun(oqn_l, itype)
            do imesh = 1, atoms%jri(itype)
              ! In Fleur the radial solutions u_lp are multiplied by a factor r by default to avoid an additional multiplication of the
              ! Jacobi determinant r^2 in radial integrals given spherical coordinates.
              varphi1(imesh, iradf, oqn_l) = rbas1(imesh, iradf, oqn_l, itype, 1) / atoms%rmsh(imesh, itype)
              varphi2(imesh, iradf, oqn_l) = rbas2(imesh, iradf, oqn_l, itype, 1) / atoms%rmsh(imesh, itype)
            end do ! imesh
          end do ! iradf
        end do ! oqn_l

        r2(:) = 0.

        ! Precalculate radial Jacobi determinant for later integrals
        do imesh = 1, atoms%jri(itype)
          r2(imesh) = atoms%rmsh(imesh, itype) * atoms%rmsh(imesh, itype)
        end do ! imesh

        vEff0MtSpH(:, :) = cmplx(0.0, 0.0)
        ptsym = atoms%ntypsy(iatom)
        do ilh = 0, lathar%nlh(ptsym)
          oqn_l = lathar%llh(ilh, ptsym)
          lm_pre = oqn_l * (oqn_l + 1)
          do imem = 1, nmem_atom(ilh, iatom)
            mqn_m = mlh_atom(imem, ilh, iatom)
            lm = lm_pre + mqn_m
            !todo one could only evaluate the vEff0MtSpH which have a contribution, i.e. the oqn_l and mqn_m which are in llh and mlh_atom
            ! maybe construct a pointer and run only over them to make it memory efficient.
            do imesh = 1, atoms%jri(itype)
              vEff0MtSpH(imesh, lm) = vEff0MtSpH(imesh, lm) + vEff0MtLh(imesh, ilh, itype) * clnu_atom(imem, ilh, iatom)
            end do ! imesh
          end do ! imem
        end do ! ilh

        hVarphi = cmplx(0.0, 0.0)
        call CalcHnGrV0Varphi( atoms, lathar, itype, iatom, lmpMax, El, varphi1, varphi2, nRadFun, vEff0MtSpH, vEff0MtLh, clnu_atom, &
          & nmem_atom, mlh_atom, grVarphiCh1, grVarphiCh2, grVarphiChLout, grVarphiChMout, hVarphi, vEff0NsphGrVarphi, r2grVeff0SphVarphi, r2grVeff0SphVarphiDummy )


        call CalcScalBasfMatElems( atoms, itype, iatom, nRadFun, r2, varphi1, varphi2, hVarphi, varphiVarphi, varphiHvarphi )
        do mqn_m2PrC = -1, 1
          call CalcHGrVarphi( atoms, itype, mqn_m2PrC, lmpMax, atoms%lmax(itype), grVarphiChMout, nRadFun, grVarPhiCh1, grVarPhiCh2,        &
                                                      & grVarphiChLout, vEff0NsphGrVarphi, El, lmpT, hGrVarphi )
        end do ! mqn_m2PrC
        do mqn_m2PrC = -1, 1
          varphiVarphiDummy(:, :, :) = 0.
          !todo attention where the hGrVarphi starts!
          call CalcScalBasfMatElems( atoms, itype, iatom, nRadFun, r2, varphi1, varphi2, hGrVarphi(:, :, :, :, mqn_m2PrC), varphiVarphiDummy, varphiHGrvarphi(:, :, :, mqn_m2PrC) )
          do jj = 1, lmpMax
            do ii = 1, lmpMax
              varphiGrVarphi(ii, jj, mqn_m2PrC, itype) = grVarphiVarphi(jj, ii, mqn_m2PrC, itype)
            end do ! ii
          end do ! jj
        end do ! mqn_m2PrC
      end do ! ieqat
    end do ! itype
    deallocate(r2grVeff0SphVarphiDummy)

    iatom = 0
    r2(:) = 1.
    iatom = 0
    do itype = 1, atoms%ntype
      do ieqat = 1, atoms%neq(itype)
        iatom = iatom + 1
        do idir = 1, 3
          call CalcScalBasfMatElems( atoms, itype, iatom, nRadFun, r2, varphi1, varphi2, r2GrVeff0SphVarphi(:, :, :, :, idir), varphiVarphiDummy, varphiGrVeff0SphVarphi(:, :, :, idir) )
        end do ! idir
      end do ! ieqat
    end do ! itype

    eps1DynMatPulTestSw = .false.
    testCompTerm3rdBraKetsVarBra = .false.
    testCompTerm3rdBraKetsVarKet = .false.
    dynMatPu3rdBraKetHepsSw = .true.
    pulayCntrb(:, :) = 0
    do ikpt = 1, kpts%nkpt
      z1nG = cmplx(0.0, 0.0)
      do iBas = 1, nv(1, ikpt) !+ atoms%nlotot
        do iband = 1, nobd(ikpt, 1)
          do idir = 1, 3
            z1nG(iBas, idir, 1, iband) = z(iBas, iband, ikpt, 1)
          end do ! iBas
        end do ! idir
      end do ! iband
      call Add1stOrdWfPulayBraKets2DynMat( atoms, kpts, qpts, sym, dimens, cell, usdus, stars, results, ikpt, iqpt, lmpMax, nRadFunMax, nv, gBasVec, &
        & ilst, kveclo, mapKpq2K, nobd, z, z1nG, iloTable, grVarphiVarphi, nRadFun, eig, kpq2kPrVec, &
       & grVarphiHvarphi, varphiVarphi, varphiHvarphi, Veff0%vpw, lmpT,&
        & eps1DynMatPulTestSw, testCompTerm3rdBraKetsVarBra, testCompTerm3rdBraKetsVarKet, dynMatPu3rdBraKetHepsSw, varphiGrVeff0SphVarphi, varphiHGrvarphi, varphiGrVarphi, pulayCntrb )
    end do ! ikpt

    if (.false.) then
      write(*, *) 'pulayCntr'
      write(*, '(3(2(es16.8,1x),3x))') pulayCntrb(1, :)
      write(*, '(3(2(es16.8,1x),3x))') pulayCntrb(2, :)
      write(*, '(3(2(es16.8,1x),3x))') pulayCntrb(3, :)
    end if

    if (any ( abs(pulayCntrb(:, :)) > 9e-7 ) ) passed = .false.

    if ( passed ) then
      write( logUnit, '(a)' ) '                                                                                |__ passed!'
      write(logUnit, *)
    else
      write( logUnit, '(a)' ) '                                                                                |__ failed!'
      write(logUnit, *)
      call JuDFT_warn('Check cancelling within 3rd braket in Pulay contribution to the dynamical matrix.', calledby='Test3rdBraKetHeps', hint='Debug!')
    end if

  end subroutine Test3rdBraKetHeps

  ! test the action of a Hamiltonian onto a gradient of phi by comparing two ways of calculating such a quantitiy
  subroutine TestActionHgrPhi(atoms, dimens, lathar, Veff0, logUnit, rbas1, rbas2, nRadFun, El, mlh_atom, nmem_atom, clnu_atom)

    use m_types
    use mod_juPhonUtils, only : Derivative, CalcChannelsGrFlpNat
    use m_jpSetupDynMat, only : CalcGrVarphiHepsGrtVarphiElem
    use m_jpSetupDynMatSF, only : CalcHGrVarphi
    use m_jpSetupDynMatHelper, only : CalcHnGrV0Varphi
    use m_juDFT_NOstopNO, only : juDFT_warn

    implicit none

    ! Type parameter
    type(t_atoms),              intent(in)  :: atoms
    type(t_dimension),          intent(in)  :: dimens
    type(t_sphhar),             intent(in)  :: lathar
    type(t_potential),          intent(in)  :: Veff0

    integer,                    intent(in)  :: logUnit

    ! Array parameters
    real,                       intent(in)  :: rbas1(:, :, 0:, :)
    real,                       intent(in)  :: rbas2(:, :, 0:, :)
    integer,                    intent(in)  :: nRadFun(0:, :)
    real,                       intent(in)  :: El(:, 0:, :, :)
    integer,                    intent(in)  :: mlh_atom(:, 0:, :)
    integer,                    intent(in)  :: nmem_atom(0:, :)
    complex,                    intent(in)  :: clnu_atom(:, 0:, :)

    ! Scalar variables
    integer                                 :: oqn_l
    integer                                 :: iradf
    integer                                 :: imesh
    integer                                 :: lmpMax
    integer                                 :: nRadFunMax
    integer                                 :: lmaxBra
    integer                                 :: iDeqat
    integer                                 :: iDtype
    integer                                 :: iDatom
    integer                                 :: ptsym
    integer                                 :: ilh
    integer                                 :: lm_pre
    integer                                 :: imem
    integer                                 :: mqn_m
    integer                                 :: lm
    integer                                 :: mqn_m2PrC
    integer                                 :: mqn_m2PrR
    integer                                 :: ii
    integer                                 :: jj
    logical                                 :: passed = .true.

    ! Array variables
    integer,       allocatable              :: lmpT(:)
    real,          allocatable              :: varphi1(:, :, :)
    real,          allocatable              :: varphi2(:, :, :)
    real,          allocatable              :: delrVarphi1(:, :, :)
    real,          allocatable              :: delrVarphi2(:, :, :)
    integer,       allocatable              :: grVarphiChLout(:, :)
    integer,       allocatable              :: grVarphiChMout(:, :)
    real,          allocatable              :: grVarphiCh1(:, :, :, :)
    real,          allocatable              :: grVarphiCh2(:, :, :, :)
    complex,       allocatable              :: vEff0MtSpH(:, :)
    real,          allocatable              :: grVarphiGrtVarphi(:, :, :, :, :)
    complex,       allocatable              :: grVarphiHpreGrtVarphi(:, :, :, :, :)
    complex,       allocatable              :: grVarphiGrtVeff0SphVarphibench(:, :, :, :, :)
    complex,       allocatable              :: grVarphiGrtVeff0SphVarphi(:, :, :, :, :)
    complex,       allocatable              :: vEff0NsphGrVarphi(:, :, :, :, :)
    complex,       allocatable              :: r2grVeff0SphVarphi(:, :, :, :, :)
    complex,       allocatable              :: r2grVeff0SphVarphiDummy(:, :, :, :, :)
    complex,       allocatable              :: hVarphi(:, :, :, :)
    complex,       allocatable              :: hGrVarphi(:, :, :, :, :)
    complex,       allocatable              :: grVarphiHGrVarphi(:, :, :, :, :)
    real,          allocatable              :: r2(:)

    write(logUnit, '(a)') 'Test consistency of Hamiltonian action onto gradVarphi.'
    write(logUnit, '(a)') '-------------------------------------------------------'

    ! Quantities for initialization
    allocate( lmpT(atoms%ntype) )
    lmpT(:) = 0
    do iDtype = 1, atoms%ntype
      lmpT(iDtype) = sum( [ ( (2 * oqn_l + 1)* nRadFun(oqn_l, iDtype), oqn_l = 0, atoms%lmax(iDtype) ) ] )
    end do ! iDtype
    lmpMax     = maxval( lmpT(:) )
    nRadFunMax = maxval( nRadFun(:, :) )

    ! Allocation of required arrays
    allocate( varphi1(atoms%jmtd, nRadFunMax, 0:atoms%lmaxd), varphi2(atoms%jmtd, nRadFunMax, 0:atoms%lmaxd) )
    allocate( delrVarphi1( atoms%jmtd, nRadFunMax, 0:atoms%lmaxd), delrVarphi2( atoms%jmtd, nRadFunMax, 0:atoms%lmaxd) )
    allocate( grVarphiCh1(atoms%jmtd, 2, lmpMax, -1:1), grVarphiCh2(atoms%jmtd, 2, lmpMax, -1:1), &
            & grVarphiChLout(2, 0:atoms%lmaxd), grVarphiChMout(-atoms%lmaxd:atoms%lmaxd, -1:1) )
    allocate( vEff0MtSpH( atoms%jmtd, 0:dimens%lmd), vEff0NsphGrVarphi(2, atoms%jmtd, 0:(atoms%lmaxd + 3)**2 - 1, lmpMax, -1:1), &
            & r2grVeff0SphVarphi(2, atoms%jmtd, 0:(atoms%lmaxd + 3)**2 - 1, lmpMax, 3) )
    allocate( r2grVeff0SphVarphiDummy(2, atoms%jmtd, 0:(atoms%lmaxd + 3)**2 - 1, lmpMax, 3) )
    allocate( hVarphi(2, atoms%jmtd, 0:(atoms%lmaxd + 3)**2 - 1, lmpMax)  )
    allocate( hGrVarphi(2, atoms%jmtd, (atoms%lmaxd + 2)**2, lmpMax, -1:1))
    allocate( grVarphiHGrVarphi(lmpMax, lmpMax, -1:1, -1:1, atoms%nat) )
    allocate( r2(atoms%jmtd) )
    allocate( grVarphiGrtVarphi(lmpMax, lmpMax, -1:1, -1:1, atoms%ntype), &
            & grVarphiHpreGrtVarphi(lmpMax, lmpMax, -1:1, -1:1, atoms%nat) )
    allocate( grVarphiGrtVeff0SphVarphibench(lmpMax, lmpMax, -1:1, 1:3, atoms%nat) )
    allocate( grVarphiGrtVeff0SphVarphi(lmpMax, lmpMax, -1:1, 1:3, atoms%nat) )

    ! Initialization
    grVarphiGrtVeff0SphVarphibench(:, :, :, :, :) = cmplx(0., 0.)
    grVarphiGrtVeff0SphVarphi(:, :, :, :, :) = cmplx(0., 0.)
    grVarphiHpreGrtVarphi(:, :, :, :, :) = cmplx(0., 0.)
    grVarphiGrtVarphi(:, :, :, :, :) = cmplx(0., 0.)
    grVarphiHGrVarphi(:, :, :, :, :) = cmplx(0., 0.)
    hGrVarphi(:, :, :, :, :) = cmplx(0., 0.)
    r2grVeff0SphVarphi(:, :, :, :, :) = cmplx(0., 0.)
    r2grVeff0SphVarphiDummy(:, :, :, :, :) = cmplx(0., 0.)
    vEff0NsphGrVarphi(:, :, :, :, :) = cmplx(0., 0.)

    iDatom = 0
    do iDtype = 1, atoms%ntype
      r2(:) = 0.
      ! Jacobi determinant of radial part
      do imesh = 1, atoms%jri(iDtype)
        r2(imesh) = atoms%rmsh(imesh, iDtype) * atoms%rmsh(imesh, iDtype)
      end do ! imesh
      do iDeqat = 1, atoms%neq(iDtype)
        iDatom = iDatom + 1

        ! Calculate the gradient of varphi
        varphi1(:, :, :) = 0.
        varphi2(:, :, :) = 0.
        delrVarphi1(:, :, :) = 0.
        delrVarphi2(:, :, :) = 0.
        do oqn_l = 0, atoms%lmax(iDtype)
          do iradf = 1, nRadFun(oqn_l, iDtype)
            do imesh = 1, atoms%jri(iDtype)
              ! In Fleur the radial solutions u_lp are multiplied by a factor r by default to avoid an additional multiplication of the
              ! Jacobi determinant r^2 in radial integrals given spherical coordinates.
              varphi1(imesh, iradf, oqn_l) = rbas1(imesh, iradf, oqn_l, iDtype) / atoms%rmsh(imesh, iDtype)
              varphi2(imesh, iradf, oqn_l) = rbas2(imesh, iradf, oqn_l, iDtype) / atoms%rmsh(imesh, iDtype)
            end do ! imesh
            ! Precalculate partial derivatives of varphis in r-direction since it is needed twice
            call Derivative( varphi1(1:atoms%jri(iDtype), iradf, oqn_l), iDtype, atoms, delrVarphi1(1:atoms%jri(iDtype), iradf, oqn_l) )
            call Derivative( varphi2(1:atoms%jri(iDtype), iradf, oqn_l), iDtype, atoms, delrVarphi2(1:atoms%jri(iDtype), iradf, oqn_l) )
          end do ! iradf
        end do ! oqn_l

        call CalcChannelsGrFlpNat( atoms, iDtype, nRadFun, varphi1, varphi2, delrVarphi1, delrVarphi2, grVarphiChLout, grVarphiChMout, &
                                                                                                        & grVarphiCh1, grVarphiCh2 )

        ! Expand the coefficients of the lattice-harmonic given potential into spherical-harmonic coefficients for the given atom.
        vEff0MtSpH(:, :) = cmplx(0.0, 0.0)
        ptsym = atoms%ntypsy(iDatom)
        do ilh = 0, lathar%nlh(ptsym)
          oqn_l = lathar%llh(ilh, ptsym)
          lm_pre = oqn_l * (oqn_l + 1)
          do imem = 1, nmem_atom(ilh, iDatom)
            mqn_m = mlh_atom(imem, ilh, iDatom)
            lm = lm_pre + mqn_m
            !todo one could only evaluate the vEff0MtSpH which have a contribution, i.e. the oqn_l and mqn_m which are in llh and mlh_atom
            ! maybe construct a pointer and run only over them to make it memory efficient.
            do imesh = 1, atoms%jri(iDtype)
              vEff0MtSpH(imesh, lm) = vEff0MtSpH(imesh, lm) + Veff0%vr(imesh, ilh, iDtype, 1) * clnu_atom(imem, ilh, iDatom)
            end do ! imesh
          end do ! imem
        end do ! ilh
        hVarphi = cmplx(0.0, 0.0)
        call CalcHnGrV0Varphi( atoms, lathar, iDtype, iDatom, lmpMax, El, varphi1, varphi2, nRadFun, vEff0MtSpH, Veff0%vr(:, :, :, 1), clnu_atom, &
          & nmem_atom, mlh_atom, grVarphiCh1, grVarphiCh2, grVarphiChLout, grVarphiChMout, hVarphi, vEff0NsphGrVarphi, r2grVeff0SphVarphi, r2grVeff0SphVarphiDummy)

        ! Call routine to calculate H grVarphi
        ! For testing reasons we want this cutoff because we project on the gradient of varphi
        lmaxBra = atoms%lmax(iDtype)! + 1
        do mqn_m2PrC = -1, 1
          call CalcHGrVarphi( atoms, iDtype, mqn_m2PrC, lmpMax, lmaxBra, grVarphiChMout, nRadFun, grVarPhiCh1, grVarPhiCh2,        &
                                                      & grVarphiChLout, vEff0NsphGrVarphi, El, lmpT, hGrVarphi )

          call CalcGrVarphiHepsGrtVarphiElem( atoms, iDtype, iDatom, mqn_m2PrC, lmpMax, nRadFun, grVarphiChLout, grVarphiChMout, r2, &
            & grVarPhiCh1, grVarPhiCh2, El, lmpT, vEff0NsphGrVarphi, r2grVeff0SphVarphi, grVarphiGrtVarphi, grVarphiHpreGrtVarphi, &
            & grVarphiGrtVeff0SphVarphibench)
        end do ! mqn_m2PrC

        ! Call routine to calculate volInt grVarphi
        call calcProjOnGrVarphi( atoms, iDtype, iDatom, nRadFun, r2, grVarphiChLout, grVarphiChMout, grVarphiCh1, grVarphiCh2, hGrVarphi, grVarphiHGrVarphi )

        ! Set Jacobi determinant to 1 because it is contained in r2grVeff0SphVarphi
        r2(:) = 1.
        call calcProjOnGrVarphi( atoms, iDtype, iDatom, nRadFun, r2, grVarphiChLout, grVarphiChMout, grVarphiCh1, grVarphiCh2, r2grVeff0SphVarphi, grVarphiGrTVeff0SphVarphi )

        if (.false.) then
          do mqn_m2PrC = -1, 1
            do mqn_m2PrR = -1, 1
              do jj = 1, lmpMax
                do ii = 1, lmpMax
                  write(4000, '(4i8,2f15.8)') mqn_m2PrC, mqn_m2PrR, ii, jj, grVarphiHgrVarphi(ii, jj, mqn_m2PrR, mqn_m2PrC, iDatom)
                  write(4001, '(4i8,2f15.8)') mqn_m2PrC, mqn_m2PrR, ii, jj, grVarphiHpreGrtVarphi(ii, jj, mqn_m2PrR, mqn_m2PrC, iDatom)
                  write(4002, '(4i8,2f15.8)') mqn_m2PrC, mqn_m2PrR, ii, jj, grVarphiGrtVeff0SphVarphibench(ii, jj, mqn_m2PrR, mqn_m2PrC + 2, iDatom)
                  write(4003, '(4i8,2f15.8)') mqn_m2PrC, mqn_m2PrR, ii, jj, grVarphiGrTVeff0SphVarphi(ii, jj, mqn_m2PrR, mqn_m2PrC + 2, iDatom)
                end do ! ii
              end do ! jj
            end do ! mqn_m2PrR
          end do ! mqn_m2PrC
        end if

        if ( any(abs(grVarphiHgrVarphi(:, :, :, :, :) - grVarphiHpreGrtVarphi(:, :, :, :, :)) > 9e-6) ) passed = .false.
        if ( any(abs(grVarphiGrtVeff0SphVarphibench(:, :, :, :, :) - grVarphiGrTVeff0SphVarphi(:, :, :, :, :)) > 9e-6) ) passed = .false.


        if (passed) then
          write (logUnit, '(a)')'                                                      |__ passed!'
          write(logUnit, *)
        else
          write (logUnit, '(a)')'                                                      |__ failed!'
          write(logUnit, *)
          call juDFT_warn('Test consistency of Hamiltonian action onto gradVarphi.', &
            & calledby='testActionHgrPhi', hint='Debug.')
        end if

      end do ! iDeqat
    end do ! iDtype

  end subroutine TestActionHgrPhi

  ! Project an array MT quantitiy which has acted on a ket with a gradvarphi in the bra
  subroutine calcProjOnGrVarphi( atoms, iDtype, iDatom, nRadFun, r2, grVarphiChLout, grVarphiChMout, grVarphiCh1, grVarphiCh2, hGrVarphi, grVarphiHGrVarphi )

    use m_types, only : t_atoms
    !use m_intgr, only : intgr3NoIntp
    use m_intgr, only : intgr3LinIntp

    implicit none

    ! Type parameters
    type(t_atoms),              intent(in)  :: atoms

    ! Scalar parameters
    integer,                    intent(in) :: iDtype
    integer,                    intent(in) :: iDatom

    ! Array parameters
    integer,                    intent(in)  :: nRadFun(0:, :)
    real,                       intent(in)  :: r2(:)
    integer,                    intent(in)  :: grVarphiChLout(:, 0:)
    integer,                    intent(in)  :: grVarphiChMout(-atoms%lmaxd:, -1:)
    real,                       intent(in)  :: grVarphiCh1(:, :, :, -1:)
    real,                       intent(in)  :: grVarphiCh2(:, :, :, -1:)
    complex,                    intent(in)  :: hGrVarphi(:, :, :, :, -1:)
    complex,                    intent(out) :: grVarphiHGrVarphi(:, :, -1:, -1:, :)

    ! Scalar variables
    integer                                 :: mqn_m2PrC
    integer                                 :: mqn_m2PrR
    integer                                 :: lmp
    integer                                 :: oqn_l
    integer                                 :: mqn_m
    integer                                 :: iradf
    integer                                 :: lmp1Pr
    integer                                 :: oqn_l1Pr
    integer                                 :: mqn_m1Pr
    integer                                 :: iradf1Pr
    integer                                 :: ichan1Pr
    integer                                 :: oqn_l3Pr
    integer                                 :: mqn_m3Pr
    integer                                 :: imesh
    integer                                 :: lm3Pr
    real                                    :: integralR
    real                                    :: integralI

    real,          allocatable              :: intgrdR(:)
    real,          allocatable              :: intgrdI(:)

    allocate( intgrdR(atoms%jri(iDtype)), intgrdI(atoms%jri(iDtype)) )

    grVarphiHGrVarphi(:, :, :, :, :) = cmplx(0., 0.)

    do mqn_m2PrC = -1, 1
      do mqn_m2PrR = -1, 1
        lmp = 0
        do oqn_l = 0, atoms%lmax(iDtype)
          do mqn_m = -oqn_l, oqn_l
            do iradf = 1, nRadFun(oqn_l, iDtype)
              lmp = lmp + 1
              lmp1Pr = 0
              do oqn_l1Pr = 0, atoms%lmax(iDtype)
                do mqn_m1Pr = -oqn_l1Pr, oqn_l1Pr
                  mqn_m3Pr = grVarphiChMout(mqn_m1Pr, mqn_m2PrR)
                  do iradf1Pr = 1, nRadFun(oqn_l1Pr, iDtype)
                    lmp1Pr = lmp1Pr + 1
                    do ichan1Pr = 1, 2
                      oqn_l3Pr = grVarphiChLout(ichan1Pr, oqn_l1Pr)
                      if ( ( abs(mqn_m3Pr) > oqn_l3Pr )  .or. ( oqn_l3Pr < 0 ) )  cycle
                      lm3Pr = oqn_l3Pr * (oqn_l3Pr + 1) + mqn_m3Pr + 1
                      intgrdR = 0.
                      intgrdI = 0.
                      do imesh = 1, atoms%jri(iDtype)
                        intgrdR(imesh) = r2(imesh) * &
                              & real( grVarphiCh1(imesh, ichan1Pr, lmp1Pr, mqn_m2PrR) * hGrVarphi(1, imesh, lm3Pr, lmp, mqn_m2PrC) &
                                  & + grVarphiCh2(imesh, ichan1Pr, lmp1Pr, mqn_m2PrR) * hGrVarphi(2, imesh, lm3Pr, lmp, mqn_m2PrC) )
                        intgrdI(imesh) = r2(imesh) * &
                             & aimag( grVarphiCh1(imesh, ichan1Pr, lmp1Pr, mqn_m2PrR) * hGrVarphi(1, imesh, lm3Pr, lmp, mqn_m2PrC) &
                                  & + grVarphiCh2(imesh, ichan1Pr, lmp1Pr, mqn_m2PrR) * hGrVarphi(2, imesh, lm3Pr, lmp, mqn_m2PrC) )
                      end do ! imesh
                      call intgr3LinIntp(intgrdR(:), atoms%rmsh(1, iDtype), atoms%dx(iDtype), atoms%jri(iDtype), integralR, 1)
                      call intgr3LinIntp(intgrdI(:), atoms%rmsh(1, iDtype), atoms%dx(iDtype), atoms%jri(iDtype), integralI, 1)
                      grVarphiHGrVarphi(lmp1Pr, lmp, mqn_m2PrR, mqn_m2PrC, iDatom) = &
                                        & grVarphiHGrVarphi(lmp1Pr, lmp, mqn_m2PrR, mqn_m2PrC, iDatom) + cmplx(integralR, integralI)
                    end do ! ichan1Pr
                  end do ! iradf1Pr
                end do ! mqn_m1Pr
              end do ! oqn_l1Pr
            end do ! iradf
          end do ! mqn_m
        end do ! oqn_l
      end do ! mqn_m2PrR
    end do ! mqn_m2PrC

  end subroutine calcProjOnGrVarphi

end module m_jpTestDynMatPulay
