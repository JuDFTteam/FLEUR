module m_jpTestGoldsteinCond

  implicit none

  contains

  ! Calls every test, that checks the Goldstone condition partly or in total
  ! todo Everywhere write whether the rho is given in stars or lattice harmonics or not, the same for the potentials
  subroutine TestGoldsteinCond( atoms, cell, lathar, stars, dimens, kpts, qpts, results, Veff0, sym, usdus, enpara, input, ngdp, &
      & memd_atom, logUnit, testGoldsteinRemainingSw, test1st2ndPulDynMatEps1, test1st2ndPulDynMatCancel, testGoldsteinSurfSw, &
      & rho0IR, rho0MT, clnu_atom, nmem_atom, mlh_atom, gdp, gBas, mapGbas, kveclo, nv, kpq2kPrVec, eig, nobd, mapKpq2K, nRadFun, &
      & rbas1, rbas2, El, iloTable, z0, uuilon, duilon, ulouilopn, ilo2p )

#include "cppmacro.h"

    use m_types

    implicit none

    type(t_atoms),                  intent(in) :: atoms
    type(t_cell),                   intent(in) :: cell
    type(t_sphhar),                 intent(in) :: lathar
    type(t_stars),                  intent(in) :: stars
    type(t_dimension),              intent(in) :: dimens
    type(t_kpts),                   intent(in) :: kpts
    type(t_kpts),                   intent(in) :: qpts
    type(t_results),                intent(in) :: results
    type(t_potential),              intent(in) :: Veff0
    type(t_sym),                    intent(in) :: sym
    type(t_usdus),                  intent(in) :: usdus
    type(t_enpara),                 intent(in) :: enpara
    type(t_input),                  intent(in) :: input

    ! Scalar parameters
    integer,                        intent(in) :: ngdp
    integer,                        intent(in) :: memd_atom
    integer,                        intent(in) :: logUnit
    logical,                        intent(in) :: testGoldsteinRemainingSw
    logical,                        intent(in) :: test1st2ndPulDynMatEps1
    logical,                        intent(in) :: test1st2ndPulDynMatCancel
    logical,                        intent(in) :: testGoldsteinSurfSw

    ! Array parameters
    complex,                        intent(in) :: rho0IR(:, :)
    real,                           intent(in) :: rho0MT(:, 0:, :, :)
    complex,                        intent(in) :: clnu_atom(:, 0:, :)
    integer,                        intent(in) :: nmem_atom(0:, :)
    integer,                        intent(in) :: mlh_atom(:, 0:, :)
    integer,                        intent(in) :: gdp(:, :)
    integer,                        intent(in) :: gBas(:, :)
    integer,                        intent(in) :: mapGbas(:, :, :)
    integer,                        intent(in) :: kveclo(:,:)
    integer,                        intent(in) :: nv(:, :)
    integer,                        intent(in) :: kpq2kPrVec(:, :, :)
    real,                           intent(in) :: eig(:, :, :)
    integer,                        intent(in) :: nobd(:, :)
    integer,                        intent(in) :: mapKpq2K(:, :)
    integer,                        intent(in) :: nRadFun(0:, :)
    real,                           intent(in) :: rbas1(:, :, 0:, :, :)
    real,                           intent(in) :: rbas2(:, :, 0:, :, :)
    real,                           intent(in) :: El(:, 0:, :, :)
    integer,                        intent(in) :: iloTable(:, 0:, :)
    MCOMPLEX,                       intent(in) :: z0(:, :, :, :)
    real,                           intent(in) :: uuilon(:, :)
    real,                           intent(in) :: duilon(:, :)
    real,                           intent(in) :: ulouilopn(:, :, :)
    integer,                        intent(in) :: ilo2p(:, :)
    ! code

    ! Scalar local variables
    logical                                    :: passedDMPuMEWithEps1

    if ( test1st2ndPulDynMatEps1 .or. test1st2ndPulDynMatCancel .or. testGoldsteinSurfSw .or. testGoldsteinRemainingSw) then
      write(*, *)
      write(*, '(a)') 'Initiating Goldstein condition test(s)...'
      write(*, '(a)') '-----------------------------------------'
      write (logUnit, * )
      write (logUnit, '(a)') 'Dynamical matrix Goldstein condition test(s).'
      write (logUnit, '(a)') '*********************************************'
      write (logUnit, * )
    else
      write(*, '(a)') '------------------------------------------------------'
      write(*, '(a)') 'DISABLED dynamical matrix Goldstein condition test(s)!'
      write(*, '(a)') '------------------------------------------------------'
      write(logUnit, *)
      write(logUnit, '(a)') 'DISABLED dynamical matrix Goldstein condition test(s)!'
      write(logUnit, '(a)') '******************************************************'
      return
    end if



    if (test1st2ndPulDynMatEps1) then
      write(*, '(2x,a)') 'Performing TestDynMatPuWithEps1...'
      write (logUnit, '(a)')   'Compare eps(1) with parts of the Pulay contribution matrix elements to the dynamical matrix.'
      write ( logUnit, '(a)' ) '--------------------------------------------------------------------------------------------'
      call TestDynMatPuWithEps1( atoms, lathar, dimens, input, stars, kpts, cell, sym, usdus, enpara, results, ngdp, qpts, &
        & mlh_atom, nmem_atom, clnu_atom, nobd, nv, gBas, mapGbas, z0, gdp, kpq2kPrVec, kveclo, nRadFun, rho0MT, rho0IR, memd_atom, iloTable, El,  &
        & rbas1, rbas2, uuilon, duilon, ulouilopn, ilo2p, Veff0%vr, mapKpq2K, eig, Veff0%vpw_uw, Veff0%vpw, logUnit, passedDMPuMEWithEps1 )
      if ( passedDMPuMEWithEps1 ) then
        write(logUnit, '(a)') '                                                                                           |__ passed!'
        write(logUnit, *)
      else
        write(logUnit, '(a)') '                                                                                           |__ failed!'
        write(logUnit, *)
      end if
    else
      write(*, '(2x,a)') 'DISABLED TestDynMatPuWithEps1!'

      write (logUnit, '(a)')   'Compare eps(1) with parts of the Pulay contribution matrix elements to the dynamical matrix.'
      write ( logUnit, '(a)' ) '--------------------------------------------------------------------------------------------'
      write (logUnit, '(a)')   '                                                                                           |__ DISABLED!'
      write (logUnit, *)
    end if

    if (test1st2ndPulDynMatCancel) then
      write(*, '(2x,a)') 'Performing TestComp2ndN1stOrdBasFunc...'
      call TestComp2ndN1stOrdBasFunc(atoms, lathar, dimens, kpts, cell, sym, qpts, usdus, stars, results, .true., .true., logUnit, nRadFun, mlh_atom, &
        & nmem_atom, clnu_atom, Veff0%vr(:, :, :, 1), rbas1, rbas2, El, nobd, nv, z0, mapKpq2K, eig, gBas, mapGbas, kveclo, iloTable,&
                             & kpq2kPrVec, Veff0%vpw )
    else
      write(*, '(2x,a)') 'DISABLED TestComp2ndN1stOrdBasFunc!'

      write (logUnit, '(a)')   'Test Goldstone condition for the Pulay terms.'
      write (logUnit, '(a)' ) '---------------------------------------------'
      write(logUnit, '(a)') '                                            |__ DISABLED!'
      write (logUnit, *)
    end if

    if (testGoldsteinSurfSw) then
      write(*, '(2x,a)') 'Performing TestGoldsteinSurf...'
      call TestGoldsteinSurf( atoms, dimens, stars, cell, results, Veff0, kpts, qpts, lathar, sym, usdus, ngdp, logUnit, nRadFun, &
        & nobd, gdp, mapKpq2K, kpq2kPrVec, nv, mapGbas, gBas, z0, eig, rbas1(:, :, :, :, 1), rbas2(:, :, :, :, 1), nmem_atom, &
        & mlh_atom, clnu_atom, kveclo, iloTable, El )
    else
      write(*, '(2x,a)') 'DISABLED TestGoldsteinSurf!'

      write (logUnit, '(a)')   'Test Goldstone condition for the surface integrals.'
      write (logUnit, '(a)' ) '---------------------------------------------------'
      write (logUnit, '(a)')  '                                                  |__ DISABLED!'
      write (logUnit, *)
    end if

    if (testGoldsteinRemainingSw) then
      write(*, '(2x,a)') 'Performing TestGoldsteinRemaining...'
      call testGoldsteinRemaining( atoms, cell, lathar, stars, dimens, kpts, qpts, results, Veff0, sym, usdus, ngdp, memd_atom, rho0IR, rho0MT, clnu_atom, nmem_atom, mlh_atom, gdp, gBas, mapGbas, kveclo, nv, kpq2kPrVec, eig, nobd, z0, Veff0%vpw, mapKpq2K, nRadFun, rbas1(:, :, :, :, :), rbas2(:, :, :, :, :), El, iloTable, logUnit )
    else
      write(*, '(2x,a)') 'DISABLED TestGoldsteinRemaining!'

      write(logUnit, '(a)') 'Testing sum of terms that should vanish due to Goldstone condition meeting scaled Schroedinger equation'
      write(logUnit, '(a)') '-------------------------------------------------------------------------------------------------------'
      write (logUnit, '(a)')'                                                                                                       |__ DISABLED!'
      write (logUnit, *)
    end if

    write(*, *)

  end subroutine TestGoldsteinCond

  !>-------------------------------------------------------------------------------------------------------------------------------------
  !> @author
  !> Christian-Roman Gerhorst
  !>
  !> @brief
  !> Tests canceling of remaining terms for Goldstone condition at the Pulay matrix elements.
  !>
  !> @details
  !> For q = 0 the phonon energies of the system have to be zero. Therefore, the dynamical matrix should be zero for q = 0. Also,
  !> z1 = -i (k + G) z0 and first variations of the gradient are the gradients of this quantities. After exploiting this, a lot
  !> of terms cancel each other in the matrix element contributions of the Pulay terms for the dynamical matrix. But, some still
  !> remain and should in principle cancel because the Schroedinger equation holds. In this case, it is something like a scaled
  !> Schroedinger equation and the Hamiltonian is not programmed in a hermitian way to ensure consistency with IR and MT and the
  !> this scaled Schroedinger equation projected onto another wave function has to cancel away which is tested here for the
  !> Goldstone condition to be fulfilled.
  !>
  !> @param[in] atoms        : Atoms type, see types.f90.
  !> @param[in] dimens       : Dimension type, see types.f90.
  !> @param[in] kpts         : K-points type, see types.f90.
  !> @param[out] eig         : Contains Kohn\--Sham eigenvalues.
  !> @param[out] kveclo      : Basis G-vectors of local orbitals.
  !> @param[out] nv          : Number of LAPW G-basis vectors for given k-point.
  !> @param[out] mapGbas     : For various k-points G-basis vectors occur more than once, thus they are only stored once in juPhon. This
  !>                           pointer array contains the right index for GbasVec array to "unfold" G-basis vectors again.
  !> @param[out] gBas        : G-basis vectors
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
  subroutine testGoldsteinRemaining( atoms, cell, lathar, stars, dimens, kpts, qpts, results, Veff0, sym, usdus, &
      & ngdp, memd_atom, rho0IRst, rho0MTlh, clnu_atom, nmem_atom, mlh_atom, gdp, gBas, gBasUnwrap, kveclo, nv, kpq2kPrVec, eig,   &
      & nobd, z, w_vEff0IRst, mapKpq2K, nRadFun, rbas1, rbas2, El, iloTable, logUnit )

#include "cppmacro.h"

    use m_types, only : t_atoms, t_cell, t_sphhar, t_stars, t_dimension, t_kpts, t_results, t_potential, t_sym, t_usdus
    use m_jpSetupDynMat, only : calcPsi1HepsPsi1IR
    use m_jpSetupDynMatSF, only : CalcSFintIRPsi1HepsPsi, CalcSFintIRPsiHepsPsi1
    use m_jpPotDensHelper, only : WarpIRPot, CalcIRdVxcKern, CalcMTdVxcKern, convertStar2G
    use m_jpGrVeff0, only : GenGrVeff0
    use m_jpConstants, only : iu
    use mod_juPhonUtils, only : Calc2ArgIntIR, CalcGrR2FinLh
    use m_jpSetupDynMatHelper, only : CalcHnGrV0Varphi
    use m_jpSetupDynMat, only : CalcScalBasfMatElems, CalcSelfAdjCorrPhi
    use m_abcof3, only : abcof3
    use m_od_types, only : od_inp, od_sym
    use mod_juPhonUtils, only : Derivative
    use m_juDFT_NOstopNO, only : juDFT_warn

    implicit none

    ! Type parameters
    type(t_atoms),                  intent(in) :: atoms
    type(t_cell),                   intent(in) :: cell
    type(t_sphhar),                 intent(in) :: lathar
    type(t_stars),                  intent(in) :: stars
    type(t_dimension),              intent(in) :: dimens
    type(t_kpts),                   intent(in) :: kpts
    type(t_kpts),                   intent(in) :: qpts
    type(t_results),                intent(in) :: results
    type(t_potential),              intent(in) :: Veff0
    type(t_sym),                    intent(in) :: sym
    type(t_usdus),                  intent(in) :: usdus

    ! Scalar parameters
    integer,                        intent(in) :: ngdp
    integer,                        intent(in) :: memd_atom
    integer,                       intent(in)  :: logUnit

    ! Array parameters
    complex,                        intent(in) :: rho0IRst(:, :)
    real,                           intent(in) :: rho0MTlh(:, 0:, :, :)
    complex,                        intent(in) :: clnu_atom(:, 0:, :)
    integer,                        intent(in) :: nmem_atom(0:, :)
    integer,                        intent(in) :: mlh_atom(:, 0:, :)
    integer,                        intent(in) :: gdp(:, :)
    integer,                        intent(in) :: gBas(:, :)
    integer,                        intent(in) :: gBasUnwrap(:, :, :)
    integer,                        intent(in) :: kveclo(:,:)
    integer,                        intent(in) :: nv(:, :)
    integer,                        intent(in) :: kpq2kPrVec(:, :, :)
    real,                           intent(in) :: eig(:, :, :)
    integer,                        intent(in) :: nobd(:, :)
    complex,                        intent(in) :: w_vEff0IRst(:, :)
    integer,                        intent(in) :: mapKpq2K(:, :)
    integer,                        intent(in) :: nRadFun(0:, :)
    real,                           intent(in) :: rbas1(:, :, 0:, :, :)
    real,                           intent(in) :: rbas2(:, :, 0:, :, :)
    real,                           intent(in) :: El(:, 0:, :, :)
    integer,                        intent(in) :: iloTable(:, 0:, :)
    MCOMPLEX,                       intent(in) :: z(:, :, :, :)

    ! Type variable
    type(od_inp)                               :: odi
    type(od_sym)                               :: ods

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
    logical                                    :: fullGrVeff0
    logical                                    :: testGoldstein
    logical                                    :: harSw
    logical                                    :: extSw
    logical                                    :: xcSw
    integer                                    :: iG
    integer                                    :: iband
    integer                                    :: nmat
    integer                                    :: ikpt
    integer                                    :: iqpt
    complex                                    :: wholeUnitCellIR
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
    logical                                    :: passed = .true.
    logical                                    :: coretail


    ! rho1IR     : interstitial first variation of density
    ! rho1MT     : muffin-tin first variation of density
    ! w_vEff1IR  : warped first variation of effective potential in interstitial
    ! vEff1MT    : first variation of effective potential in muffin-tin
    ! Array variables
    complex,           allocatable             :: rho1IRContainer(:, :, :)
    complex,           allocatable             :: rho1MTContainer(:, :, :, :, :)
    complex,           allocatable             :: w_vEff1IRContainer(:, :, :)
    complex,           allocatable             :: vEff1MTContainer(:, :, :, :, :)
    real,              allocatable             :: gaussWghts(:) ! gaussian weights belonging to gausPts
    complex,           allocatable             :: ylm(:, : )
    real,              allocatable             :: dKernMTGPts(:, :, :)
    real,              allocatable             :: r2Rho0MT(:, :, :)
    complex,           allocatable             :: grRho0MT(:, :, :, :)
    complex,           allocatable             :: r2GrRho0MT(:, :, :, :)
    complex,           allocatable             :: grVxcIRKern(:)
    complex,           allocatable             :: grVeff0IR(:, :)
    complex,           allocatable             :: grVeff0MT(:, :, :, :)
    complex,           allocatable             :: grRho0IR(:, :)
    complex,           allocatable             :: w_grVeff0IR(:, :)
    complex,           allocatable             :: rho0IRPw(:, :)
    complex,           allocatable             :: z1nG(:, :, :, :)
    complex,           allocatable             :: z1nGMat(:, :, :, :, :)
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
    complex,           allocatable             :: r2grVeff0Varphi(:, :, :, :, :)
    complex,           allocatable             :: aKpq(:, :, :)
    complex,           allocatable             :: bKpq(:, :, :)
    complex,           allocatable             :: bascof_loKpq(:, :, :, :, :)
    complex,           allocatable             :: z1Gext(:)
    complex,           allocatable             :: varphiPsi(:)
    complex,           allocatable             :: varphiHPsi(:)
    complex,           allocatable             :: varphiPsiCjg(:)
    complex,           allocatable             :: varphiHPsiCjg(:)
    complex,           allocatable             :: dynMatPu(:, :)
    complex,           allocatable             :: rho0IRnCt(:, :)
    complex,           allocatable              :: rho0MTsh(:, :, :, :)
    complex                                    :: psiHepsPsi(3, 3)
    complex                                    :: psiHepsPsiCjg(3, 3)
    real                                       :: Gext(3)
    real                                       :: kext(3)
    real                                       :: kpqExt(3)
    complex                                    :: surfInt(3, 3)
    complex                                    :: surfInts(3, 3)
    complex                                    :: dynMatPuMt(3, 3)
    complex                                    :: dynMatMT2(3, 3)

    ! The coretail has to be false, otherwise the test does not work. As the gradient of the coretails in the interstitial density
    ! is not compensated in the other terms. Within the production code, however, the z1 are converged such that they compensate for
    ! the coretail correction that is added to the gradient of the density.
    coretail = .false.
    write(logUnit, *)
    write(logUnit, *)
    write(logUnit, '(a)') 'Testing sum of terms that should vanish due to Goldstone condition meeting scaled Schroedinger equation'
    write(logUnit, '(a)') '-------------------------------------------------------------------------------------------------------'

    ! Define cumulative variable which should in sum be zero
    allocate (dynMatPu(3 * atoms%nat, 3* atoms%nat))
    dynMatPu(:, :) = cmplx(0., 0.)

    allocate( rho1IRContainer(ngdp, 3, atoms%nat), w_vEff1IRContainer(ngdp, 3, atoms%nat) )
    allocate( rho1MTContainer(atoms%jmtd, (atoms%lmaxd + 1)**2, atoms%nat, 3, atoms%nat), &
            & vEff1MTContainer(atoms%jmtd, (atoms%lmaxd + 1)**2, atoms%nat, 3, atoms%nat) )
    rho1IRContainer(:, :, :) = cmplx(0., 0.)
    w_vEff1IRContainer(:, :, :) = cmplx(0., 0.)
    vEff1MTContainer(:, :, :, :, :) = cmplx(0., 0.)
    rho1MTContainer(:, :, :, :, :) = cmplx(0., 0.)

    allocate( rho0IRpw(ngdp, 1) )
    rho0IRpw(:, 1) = cmplx(0., 0.)
    ! SPIN MISSING
    call ConvertStar2G( rho0IRst(:, 1), rho0IRpw(:, 1), stars, ngdp, gdp )

    ! Perform analytical gradient of effective potential in the interstitial region
    allocate( grRho0IR(ngdp, 3) )
    grRho0IR(:, :) = cmplx(0., 0.)
    do idir = 1, 3
      do iG = 1, ngdp
        Gext(1:3) = matmul(cell%bmat(1:3, 1:3), gdp(1:3, iG))
        grRho0IR(iG, 1:3) = iu * Gext(1:3) * rho0IRpw(iG, 1)
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
              rho0MTsh(imesh, lm, iatom, 1) = rho0MTsh(imesh, lm, iatom, 1) + rho0MTlh(imesh, ilh, itype, 1) &
                                                                                                    & * clnu_atom(imem, ilh, iatom)
            end do ! imesh
          end do ! imem
        end do ! ilh
      end do ! ieqat
    end do ! itype

    ! Calculate the gradient of the muffin-tin density for the Weinert calculation of the potential's gradient
    allocate( r2Rho0MT( atoms%jmtd, 0:lathar%nlhd, atoms%ntype) )
    allocate( grRho0MT( atoms%jmtd, ( atoms%lmaxd + 1 )**2, atoms%nat, 3) )
    r2Rho0MT(:, :, :) = 0.
    grRho0MT(:, :, :, :) = cmplx(0., 0.)

    do itype = 1, atoms%ntype
      do ilh = 0, lathar%nlhd
        do imesh = 1, atoms%jri(itype)
          ! SPIN MISSING
          r2Rho0MT(imesh, ilh, itype) = rho0MTlh(imesh, ilh, itype, 1) * atoms%rmsh(imesh, itype) * atoms%rmsh(imesh, itype)
        end do ! imesh
      end do ! ilh
    end do ! itype

    call CalcGrR2FinLH( atoms, lathar, clnu_atom, nmem_atom, mlh_atom, r2Rho0MT, r2GrRho0MT )

    do idir = 1, 3
      iatom = 0
      do itype = 1, atoms%ntype
        do ieqat = 1, atoms%neq(itype)
          iatom = iatom + 1
          do oqn_l = 0, atoms%lmax(itype)! + 1
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

    ! Precalculate quantities for the kernel of the unperturbed effective potential's gradient.
    ! SPIN MISSING
    call CalcIRdVxcKern(stars, gdp, ngdp, rho0IRst(:, 1), grVxcIRKern)
    ! SPIN MISSING
    call CalcMTdVxcKern(atoms, dimens, lathar, rho0MTlh(:, :, :, 1), nmem_atom, clnu_atom, mlh_atom, gaussWghts, ylm, dKernMTGPts)

    ! Generates gradient of the unperturbed effective potential by using the Weinert method.
    harSw = .true.
    extSw = .true.
    xcSw = .true.
    fullGrVeff0 = .true.
    testGoldstein = .false.
    ! SPIN MISSING
    call GenGrVeff0(atoms, cell, dimens, stars, ngdp, harSw, extSw, xcSw, gdp, rho0IRpw, rho0MTsh, grRho0IR, grRho0MT, &
      & gaussWghts, ylm, dKernMTGPts, grVxcIRKern, testGoldstein, fullGrVeff0, grVeff0IR, grVeff0MT )
    !todo do not forget to warp and to fill the containers

    rho0IRpw(:, 1) = cmplx(0., 0.)
    if (coretail) then
      ! SPIN MISSING
      call ConvertStar2G( rho0IRst(:, 1), rho0IRpw(:, 1), stars, ngdp, gdp )
    else
      allocate(rho0IRncT(stars%n3d, 1))
      rho0IRncT(:, :) = cmplx(0.,0.)
      rewind(7800)
      read(7800) rho0IRncT
      call ConvertStar2G( rho0IRncT(:, 1), rho0IRpw(:, 1), stars, ngdp, gdp )
    end if

    ! Perform analytical gradient of effective potential in the interstitial region
    grRho0IR(:, :) = cmplx(0., 0.)
    do idir = 1, 3
      do iG = 1, ngdp
        Gext(1:3) = matmul(cell%bmat(1:3, 1:3), gdp(1:3, iG))
        grRho0IR(iG, 1:3) = iu * Gext(1:3) * rho0IRpw(iG, 1)
      end do ! iG
    end do ! idir

    allocate(w_grVeff0IR(ngdp, 3))
    w_grVeff0IR = cmplx(0.0, 0.0)
    do idir = 1, 3
      call WarpIRPot(stars, ngdp, idir, gdp, grVeff0IR, w_grVeff0IR(:, idir))
    end do ! idir

    ! Fill up IR container variables
    do idir = 1, 3
      do iG = 1, ngdp
        rho1IRContainer(iG, idir, 1) = -grRho0IR(iG, idir)
        w_vEff1IRContainer(iG, idir, 1) = -w_grVeff0IR(iG, idir)
      end do ! iG
    end do ! idir

    ! Calculate unit cell integral for IR
    do idirC = 1, 3
      do idirR = 1, 3
        call Calc2ArgIntIR(cell, ngdp, rho1IRContainer(:, idirR, 1), w_vEff1IRContainer(:, idirC, 1), dynMatPu(idirR, idirC))
      end do ! idirR
    end do ! idirC

    if (.false.) then
      write(*, '(a)') 'surfInt grVeff grRho in IR'
      write(*, '(3(2(es16.8,1x),3x))') dynMatPu(1, :)
      write(*, '(3(2(es16.8,1x),3x))') dynMatPu(2, :)
      write(*, '(3(2(es16.8,1x),3x))') dynMatPu(3, :)
    end if

    ! Calculate first contribution of 3rd braket for the IR
    psiHepsPsi(:, :) = cmplx(0., 0.)
    allocate( z1nG(dimens%nbasfcn, 3, atoms%nat, maxval(nobd(:, :))) )
    allocate( z1nGMat(dimens%nbasfcn, 3, 3, atoms%nat, maxval(nobd(:, :))) )
    do ikpt = 1, kpts%nkpt

      nmat = nv(1, ikpt) + atoms%nlotot
      z1nGMat(:, :, :, :, :) = cmplx(0., 0.)
      kExt(1:3) = matmul(cell%bmat, kpts%bk(1:3, ikpt))
      gExt(:) = 0.
      do iBas = 1, nv(1, ikpt) !+ atoms%nlotot
        Gext(1:3) = matmul( cell%bmat(1:3, 1:3), gBas(1:3, gBasUnwrap(iBas, ikpt, 1)))
        do iband = 1, nobd(ikpt, 1)
          do idirC = 1, 3
            do idirR = 1, 3
              z1nGMat(iBas, idirR, idirC, 1, iband) = ( kExt(idirR) + Gext(idirR) ) * ( kExt(idirC) + Gext(idirC) ) * z(iBas, iband, ikpt, 1)
            end do ! idirR
          end do ! idirC
        end do ! iband
      end do ! iBas
      do iband = 1, nobd(ikpt, 1)
        do idirC = 1, 3
          do idirR = 1, 3
            wholeUnitCellIR = cmplx(0., 0.)
            call calcPsi1HepsPsi1IR( kpts, qpts, stars, dimens, cell, gBas(:, gBasUnwrap(:nv(1, ikpt), ikpt, 1)), nv, &
              & z(:, iband, ikpt, 1), z1nGMat(:, idirR, idirC, 1, iband), nmat, ikpt, iqpt, ikpt, kpq2kPrVec, &
              & w_vEff0IRst, eig, iband, wholeUnitCellIR )
            psiHepsPsi(idirR, idirC) = psiHepsPsi(idirR, idirC) + 2 * results%w_iks(iband, ikpt, 1) * wholeUnitCellIR

            wholeUnitCellIR = cmplx(0., 0.)
            call calcPsi1HepsPsi1IR( kpts, qpts, stars, dimens, cell, gBas(:, gBasUnwrap(:nv(1, ikpt), ikpt, 1)), nv, &
              & z1nGMat(:, idirR, idirC, 1, iband), z(:, iband, ikpt, 1), nmat, ikpt, iqpt, ikpt, kpq2kPrVec, &
              & w_vEff0IRst, eig, iband, wholeUnitCellIR )
            psiHepsPsi(idirR, idirC) = psiHepsPsi(idirR, idirC) + 2 * results%w_iks(iband, ikpt, 1) * wholeUnitCellIR
          end do ! idirR
        end do ! idirC
      end do ! iband
    end do ! ikpt

    if (.false.) then
      write(*, '(a)') 'IR contribution with grad in bra + contribution with grad in ket'
      write(*, '(3(2(es16.8,1x),3x))') psiHepsPsi(1, :)
      write(*, '(3(2(es16.8,1x),3x))') psiHepsPsi(2, :)
      write(*, '(3(2(es16.8,1x),3x))') psiHepsPsi(3, :)
    end if
    psiHepsPsi(:, :) = cmplx(0., 0.)
    do ikpt = 1, kpts%nkpt

      nmat = nv(1, ikpt) + atoms%nlotot
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
      do iband = 1, nobd(ikpt, 1)
        do idirC = 1, 3
          do idirR = 1, 3
            call calcPsi1HepsPsi1IR( kpts, qpts, stars, dimens, cell, gBas(:, gBasUnwrap(:nv(1, ikpt), ikpt, 1)), nv, &
              & z1nG(:, idirR, 1, iband), z1nG(:, idirC, 1, iband), nmat, ikpt, iqpt, ikpt, kpq2kPrVec, &
              & w_vEff0IRst, eig, iband, wholeUnitCellIR )
            psiHepsPsi(idirR, idirC) = psiHepsPsi(idirR, idirC) + 4 * results%w_iks(iband, ikpt, 1) * wholeUnitCellIR
          end do ! idirR
        end do ! idirC
      end do ! iband
    end do ! ikpt

    if (.false.) then
      write(*, '(a)') 'IR contribution with grad in bra and ket'
      write(*, '(3(2(es16.8,1x),3x))') psiHepsPsi(1, :)
      write(*, '(3(2(es16.8,1x),3x))') psiHepsPsi(2, :)
      write(*, '(3(2(es16.8,1x),3x))') psiHepsPsi(3, :)
    end if

    dynMatPu(1:3, 1:3) = dynMatPu(1:3, 1:3) + psiHepsPsi(1:3, 1:3)


    ! Calculate surface IR integrals containing a Psi1 with z1 = -i (k + G) z0
    iqpt = 1
    surfInts(:, :) = cmplx(0., 0.)
    surfInt(:, :) = cmplx(0., 0.)
    call CalcSFintIRPsi1HepsPsi( atoms, dimens, stars, cell, results, Veff0, kpts, qpts, ngdp, iqpt, 1, 1, 1, nobd,&
      & gdp, mapKpq2K, kpq2kPrVec, nv, gBasUnwrap, gBas, z, eig, surfInt,  .true. )
    if (.false.) then
      write(*, '(a)') 'surfInt'
      write(*, '(3(2(es16.8,1x),3x))') surfInt(1, :)
      write(*, '(3(2(es16.8,1x),3x))') surfInt(2, :)
      write(*, '(3(2(es16.8,1x),3x))') surfInt(3, :)
    end if

    surfInts(1:3, 1:3) = surfInts(1:3, 1:3) + surfInt(1:3, 1:3)

    surfInt(:, :) = cmplx(0., 0.)
    call CalcSFintIRPsiHepsPsi1( atoms, dimens, stars, cell, results, Veff0, kpts, qpts, ngdp, iqpt, 1, 1, 1,&
      & nobd, gdp, mapKpq2K, kpq2kPrVec, nv, gBasUnwrap, gBas, z, eig, surfInt, .true. )
    if (.false.) then
      write(*, '(a)') 'surfIntCjg'
      write(*, '(3(2(es16.8,1x),3x))') surfInt(1, :)
      write(*, '(3(2(es16.8,1x),3x))') surfInt(2, :)
      write(*, '(3(2(es16.8,1x),3x))') surfInt(3, :)
    end if
    surfInts(1:3, 1:3) = surfInts(1:3, 1:3) + surfInt(1:3, 1:3)

    dynMatPu(1:3, 1:3) = dynMatPu(1:3, 1:3) + surfInts(1:3, 1:3)
    if (.false.) then
      write(*, '(a)') 'sum of everything without MT'
      write(*, '(3(2(es16.8,1x),3x))') dynMatPu(1, :)
      write(*, '(3(2(es16.8,1x),3x))') dynMatPu(2, :)
      write(*, '(3(2(es16.8,1x),3x))') dynMatPu(3, :)
    end if

    ! Calculate contributions of the 1st and the 2nd braket for the MT with a z1 and i (k + G), z1 is set to the q=0 solution again

    ! We do not want the local coordinate systems to be rotated for non-representative atoms constructing the matching coefficients.
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
            & r2grVeff0SphVarphi(2, atoms%jmtd, 0:(atoms%lmaxd + 3)**2 - 1, lmpMax, 3), &
            & r2grVeff0Varphi(2, atoms%jmtd, 0:(atoms%lmaxd + 3)**2 - 1, lmpMax, 3) )
    allocate( hVarphi(2, atoms%jmtd, 0:(atoms%lmaxd + 3)**2 - 1, lmpMax)  )
    allocate( grVarphiCh1(atoms%jmtd, 2, lmpMax, -1:1), grVarphiCh2(atoms%jmtd, 2, lmpMax, -1:1), &
            & grVarphiChLout(2, 0:atoms%lmaxd), grVarphiChMout(-atoms%lmaxd:atoms%lmaxd, -1:1) )

    varphiVarphi(:, :, :) = cmplx(0., 0.)
    varphiHvarphi(:, :, :) = cmplx(0., 0.)
    vEff0NsphGrVarphi(:, :, :, :, :) = cmplx(0., 0.)
    r2grVeff0SphVarphi(:, :, :, :, :) = cmplx(0., 0.)
    r2grVeff0Varphi(:, :, :, :, :) = cmplx(0., 0.)
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

        hVarphi = cmplx(0.0, 0.0)
        call CalcHnGrV0Varphi( atoms, lathar, itype, iatom, lmpMax, El, varphi1, varphi2, nRadFun, vEff0MtSpH, Veff0%vr(:, :, :, 1), clnu_atom, &
          & nmem_atom, mlh_atom, grVarphiCh1, grVarphiCh2, grVarphiChLout, grVarphiChMout, hVarphi, vEff0NsphGrVarphi, r2grVeff0SphVarphi, r2GrVeff0Varphi )

        deallocate( vEff0MtSpH )

        ! Calculate all radial integrals with no gradients.
        ! Calculate scalar basis function matrix elements
        call CalcScalBasfMatElems( atoms, itype, iatom, nRadFun, r2, varphi1, varphi2, hVarphi, varphiVarphi, varphiHvarphi )
      end do ! ieqat
    end do ! itype

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


    dynMatPuMt(:, :) = cmplx(0., 0.)

    do ikpt = 1, kpts%nkpt
      ikpq = ikpt

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
               !     psiHepsPsi(idirR, idirC) = dot_product( abcofSumMat(:lmpT(iDtypeB), iDatomB, idirR, idirC), varphiPsi(:lmpT(iDtypeB)) )
               !     dynMatPuMT(idirR + (iDatomB - 1) * 3, idirC + (iDatomA - 1) * 3) =     &
               !       &   dynMatPuMT(idirR + (iDatomB - 1) * 3, idirC + (iDatomA - 1) * 3) &
               !todo this does not work for more than one atom
                      dynMatPuMt(idirR, idirC) = dynMatPuMt(idirR, idirC)  &
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
      write(*, '(a)') 'PsiHepsPsi + psiHepsPsiCjg in muffin-tin'
      write(*, '(3(2(es16.8,1x),3x))') dynMatPuMt(1, :)
      write(*, '(3(2(es16.8,1x),3x))') dynMatPuMt(2, :)
      write(*, '(3(2(es16.8,1x),3x))') dynMatPuMt(3, :)
    end if

    ! The Hamiltonian is programmed non-hermitianly, to be consistent with the IR. Therefore, if we want to use a correction term
    ! we need actually the correction between the spherical Hamiltonian applied to the left or applied to the right. As the basis
    ! functions are basically the same in the bra and the ket, it is sufficient to twist the matching coefficients, so that one
    ! can avoid the correction terms. Using the correction term would lead to the calculation of both the application of the
    ! spherical Hamiltonian to the left and the right. So we avoid overhead by twisting the matching coefficients.
    if ( any(abs(dynMatPuMt(:, :) + dynMatPu(:, :)) > 9e-7) ) passed = .false.

    if (passed) then
      write (logUnit, '(a)')'                                                                                                       |__ passed!'
    else
      write (logUnit, '(a)')'                                                                                                       |__ failed!'
      call juDFT_warn('Testing sum of terms that should vanish due to Goldstone condition meeting scaled Schroedinger equation failed', &
        & calledby='testGoldsteinRemaining', hint='Debug.')
    end if

  end subroutine testGoldsteinRemaining

  !>-------------------------------------------------------------------------------------------------------------------------------------
  !> @author
  !> Christian-Roman Gerhorst and Markus Betzinger
  !>
  !> @brief
  !> We test terms in 1st, 2nd and 3rd braket of Pulay matrix elements to hold condition in Equation 7.18.
  !>
  !> @details
  !. We identify terms which are similiar in the 1st, 2nd and 3rd braket of the Pulay matrix elements.
  !> In the Goldstone mode, they have to fulfill the same condition as eps1 for q = 0, i.e. they have to vanish.
  !?
  !> @param[out] atoms       : Atoms type, see types.f90.
  !> @param[out] cell        : Unit cell type, see types.f90.
  !> @param[out] sym         : Symmetries type, see types.f90.
  !> @param[out] stars       : Stars type, see types.f90.
  !> @param[out] lathar      : Lattice harmonics type, see types.f90.
  !> @param[out] input       : Input type, see types.f90.
  !> @param[out] dimens      : Dimension type, see types.f90.
  !> @param[out] enpara      : Energy parameter type, see types.f90.
  !> @param[out] kpts        : K-points type, see types.f90.
  !> @param[out] qpts        : Q-points type, see types.f90.
  !> @param[out] usdus       : Type containing quantities consisting of the radial solutions, see types.f90.
  !> @param[out] results     : Results type, see types.f90.
  !> @param[in]  logUnit     : Unit number for juPhon.log.
  !> @param[out] ngdp        : Number of G-vectors for potentials and densities.
  !> @param[out] memd_atom   : Maximal number of members in all lattice harmonics for all atoms.
  !> @param[out] kveclo      : Basis G-vectors of local orbitals.
  !> @param[out] nv          : Number of LAPW G-basis vectors for given k-point.
  !> @param[out] nRadFun     : Number of radial functions per orbital quantum number l and atom type.
  !> @param[out] iloTable    : Number of local orbital if orbital quantum number l, atom type and index p > 2 from nRadFun is given.
  !> @param[out] gBas        : G-basis vectors
  !> @param[out] mapGbas     : For various k-points G-basis vectors occur more than once, thus they are only stored once in juPhon. This
  !>                           pointer array contains the right index for GbasVec array to "unfold" G-basis vectors again.
  !> @param[out] mapKpq2K    : For a given k-point index and q-point index, this array gives the k-point set index of the result k + q
  !>                           (already mapped back to Brillouin zone).
  !> @param[out] eig         : Contains Kohn\--Sham eigenvalues.
  !> @param[out] rbas1       : Large components of radial solution, its energy derivative and u_LO
  !> @param[out] rbas2       : Small components of radial solution, its energy derivative and u_LO
  !> @param[out] El          : Contains LAPW and LO energy parameters.
  !> @param[out] z0           : Kohn-Sham eigenvectors.
  !> @param[out] rho0IR      : Star coefficients of the unperturbed and converged interstitial density parsed from Fleur.
  !> @param[out] rho0MT      : Radial coefficients of the unperturbed and converged muffin-tin densities parsed from Fleur.
  !> @param[out] gdp         : G-vectors of potentials and densities.
  !> @param[out] mlh_atom    : Magnetic quantum number m of lattice harmonic members for every atom.
  !> @param[out] nmem_atom   : Number of lattice harmonic members for every atom.
  !> @param[out] clnu_atom   : Phase mediating between stars and plane waves.
  !>
  !>--------------------------------------------------------------------------------------------------------------------------------------
  subroutine TestDynMatPuWithEps1( atoms, lathar, dimens, input, stars, kpts, cell, sym, usdus, enpara, results, ngdp, qpts,       &
      & mlh_atom, nmem_atom, clnu_atom, nobd, nv, gbas, mapGbas, z0, gdp, kpq2kPrVec, kveclo, nRadFun, rho0MT, rho0IR, memd_atom,     &
      & iloTable, El, rbas1, rbas2, uuilon, duilon, ulouilopn, ilo2p, vEff0MTLh, mapKpq2K, eig, vEff0IRst, w_Veff0IRst, logUnit, passed )

    use m_types, only : t_atoms, t_sphhar, t_dimension, t_input, t_stars, t_kpts, t_cell, t_sym, t_usdus, t_enpara, t_noco,        &
                                                                                                               & t_tlmplm, t_results
    use m_od_types, only : od_inp, od_sym

    use m_jpSternhHF, only : CalcMEPotIR, Tlmplm4V, CalcVsumMT
    use m_jpGrVeff0, only : GenGrVeff0
    use m_jpPotDensHelper, only : WarpIRPot, CalcIRdVxcKern, CalcMTdVxcKern, ConvertStar2G
    use m_abcof, only : Abcof
    use mod_juPhonUtils, only : Derivative, CalcChannelsGrFlpNat, CalcGrR2FinLh
    use m_jpSetupDynMat, only : CalcSelfAdjCorrection, CalcScalBasfMatElems, CalcVecBasfMatElems,                &
                                                & Add2ndOrdWfPulayBraKets2DynMat, CalcSelfAdjCorrPhi, Add1stOrdWfPulayBraKets2DynMat
    use m_jpSetupDynMatHelper, only : CalcHnGrV0Varphi
    use m_jpSetupDynMatSF, only : CalcHGrVarphi
    use m_jpConstants, only : iu
    use m_juDFT_NOstopNO, only : juDFT_warn

    implicit none

#include "cppmacro.h" ! Using MCOMPLEX

    ! Type parameters
    type(t_atoms),                  intent(in)  :: atoms
    type(t_sphhar),                 intent(in)  :: lathar
    type(t_dimension),              intent(in)  :: dimens
    type(t_input),                  intent(in)  :: input
    type(t_stars),                  intent(in)  :: stars
    type(t_kpts),                   intent(in)  :: kpts
    type(t_cell),                   intent(in)  :: cell
    type(t_sym),                    intent(in)  :: sym
    type(t_usdus),                  intent(in)  :: usdus
    type(t_enpara),                 intent(in)  :: enpara
    type(t_kpts),                   intent(in)  :: qpts
    type(t_results),                intent(in)  :: results

    ! Scalar parameters
    integer,                        intent(in)  :: ngdp
    integer,                        intent(in)  :: logUnit

    ! Array parameters
    integer,                        intent(in)  :: mlh_atom(:,0:,:)
    integer,                        intent(in)  :: nmem_atom(0:, :)
    complex,                        intent(in)  :: clnu_atom(:,0:,:)
    integer,                        intent(in)  :: nobd(:, :)
    integer,                        intent(in)  :: nv(:, :)
    integer,                        intent(in)  :: gbas(:, :)
    integer,                        intent(in)  :: mapGbas(:, :, :)
    MCOMPLEX,                       intent(in)  :: z0(:, :, :, :) ! Attention: see zBar
    integer,                        intent(in)  :: gdp(:, :)
    integer,                        intent(in)  :: kpq2kPrVec(:, :, :)
    integer,                        intent(in)  :: kveclo(:,:)
    integer,                        intent(in)  :: nRadFun(0:, :)
    integer,                        intent(in)  :: iloTable(:, 0:, :)
    real,                           intent(in)  :: El(:, 0:, :, :)
    real,                           intent(in)  :: rbas1(:, :, 0:, :, :)
    real,                           intent(in)  :: rbas2(:, :, 0:, :, :)
    real,                           intent(in)  :: uuilon(:, :)
    real,                           intent(in)  :: duilon(:, :)
    real,                           intent(in)  :: ulouilopn(:, :, :)
    integer,                        intent(in)  :: ilo2p(:, :)
    real,                           intent(in)  :: vEff0MTLh(:, 0:, :, :)
    integer,                        intent(in)  :: mapKpq2K(:, :)
    real,                           intent(in)  :: eig(:, :, :)
    complex,                        intent(in)  :: vEff0IRst(:,:)
    complex,                        intent(in)  :: w_Veff0IRst(:,:)
    complex,                        intent(in)  :: rho0IR(:, :)
    real,                           intent(in)  :: rho0MT(:, 0:, :, :)
    integer,                        intent(in)  :: memd_atom
    logical,                        intent(out) :: passed

    ! Type variables
    type(t_noco)                                :: noco
    type(od_inp)                                :: odi
    type(od_sym)                                :: ods
    type(t_tlmplm)                              :: tdV(3)
    type(t_tlmplm)                              :: tdV2(3)

    ! Scalar variables
    integer                                     :: itype
    integer                                     :: ilh
    integer                                     :: imesh
    integer                                     :: idir
    integer                                     :: iatom
    integer                                     :: ieqat
    integer                                     :: oqn_l
    integer                                     :: lm_pre
    integer                                     :: mqn_m
    integer                                     :: lm
    integer                                     :: ikpt
    integer                                     :: nmat
    integer                                     :: iDatom
    integer                                     :: iDtype
    integer                                     :: iDeqat
    integer                                     :: iqpt
    integer                                     :: nrBndBra
    integer                                     :: nrBndKet
    integer                                     :: maxlmp
    integer                                     :: lm2
    integer                                     :: lmp
    integer                                     :: ilo
    integer                                     :: iband
    integer                                     :: lmpMax
    integer                                     :: nRadFunMax
    integer                                     :: iradf
    integer                                     :: mqn_m2PrC
    integer                                     :: ptsym
    integer                                     :: imem
    integer                                     :: mqn_m2PrR !todo do we really need variables for R and C???
    integer                                     :: ichan
    logical                                     :: eps1DynMatPulTestSw
    integer                                     :: iBas
    logical                                     :: testComp2ndN1stOrdBasFuncSw
    logical                                     :: testCompTerm3rdBraKetsVarBra
    logical                                     :: testCompTerm3rdBraKetsVarKet
    logical                                     :: dynMatPu3rdBraKetHepsSw
    logical                                     :: numericalGradient
    integer                                     :: iG
    integer                                     :: ii
    integer                                     :: jj
    logical                                     :: harSw =.true.
    logical                                     :: extSw =.true.
    logical                                     :: xcSw =.true.
    logical                                     :: testGoldstein
    logical                                     :: fullGrVeff0
    logical                                     :: passed1st2nd
    logical                                     :: passed3rd

    ! Array variables
    real,              allocatable              :: r2Rho0MT(:, :, :)
    complex,           allocatable              :: grRho0MT(:, :, :, :) ! Dummy quantity at the moment
    complex,           allocatable              :: r2GrRho0MT(:, :, :, :)
    complex,           allocatable              :: acof(:, :, :)
    complex,           allocatable              :: bcof(:, :, :)
    complex,           allocatable              :: ccof(:, :, :, :)
    integer,           allocatable              :: ngoprI(:)
    complex,           allocatable              :: mCoef(:, :, :)
    integer,           allocatable              :: nlo_atom(:)
    complex,           allocatable              :: sumVMTsContainer(:, :, :, :)
    complex,           allocatable              :: sumVMTs2Container(:, :, :, :)
    complex,           allocatable              :: psiGrVeffPsiIR(:, :)
    complex,           allocatable              :: psiGrVeffPsiMT(:, :)
    integer,           allocatable              :: lmpT(:)
    real,              allocatable              :: varphi1(:, :, :)
    real,              allocatable              :: varphi2(:, :, :)
    real,              allocatable              :: r2(:)
    real,              allocatable              :: delrVarphi1(:, :, :)
    real,              allocatable              :: delrVarphi2(:, :, :)
    real,              allocatable              :: delrGrVarphiCh1(:, :, :, :)
    real,              allocatable              :: delrGrVarphiCh2(:, :, :, :)
    real,              allocatable              :: grVarphiCh1(:, :, :, :)
    real,              allocatable              :: grVarphiCh2(:, :, :, :)
    integer,           allocatable              :: grVarphiChLout(:, :)
    integer,           allocatable              :: grVarphiChMout(:, :)
    complex,           allocatable              :: hVarphi(:, :, :, :)
    complex,           allocatable              :: z1nG(:, :, :, :)
    real,              allocatable              :: varphiVarphi(:, :, :)
    complex,           allocatable              :: varphiHvarphi(:, :, :)
    real,              allocatable              :: grVarphiVarphi(:, :, :, :)
    complex,           allocatable              :: grVarphiHVarphi(:, :, :, :)
    complex,           allocatable              :: vEff0MtSpH(:, :)
    complex,           allocatable              :: vEff0NsphGrVarphi(:, :, :, :, :)
    complex,           allocatable              :: r2grVeff0SphVarphi(:, :, :, :, :)
    complex,           allocatable              :: r2grVeff0Varphi(:, :, :, :, :)
    real,              allocatable              :: r2Veff0MT(:, :, :)
    complex,           allocatable              :: grVeff0MT(:, :, :, :)
    complex,           allocatable              :: r2GrVeff0MT(:, :, :, :)
    complex,           allocatable              :: vEff0IRpw(:)
    complex,           allocatable              :: grVeff0IR(:, :)
    complex,           allocatable              :: w_grVeff0IR(:, :)
    complex,           allocatable              :: w_vEff1IRContainer(:, :, :)
    real,              allocatable              :: gaussWghts(:) ! gaussian weights belonging to gausPts
    complex,           allocatable              :: ylm(:, : )
    real,              allocatable              :: dKernMTGPts(:, :, :)
    complex,           allocatable              :: grVxcIRKern(:)
    complex,           allocatable              :: hGrVarphi(:, :, :, :, :)
    real,              allocatable              :: varphiVarphiDummy(:, :, :)
    complex,           allocatable              :: varphiGrVeff0SphVarphi(:, :, :, :)
    complex,           allocatable              :: varphiGrVeff0Varphi(:, :, :, :)
    complex,           allocatable              :: varphiHGrvarphi(:, :, :, :)
    real,              allocatable              :: varphiGrVarphi(:, :, :, :)
    complex,           allocatable              :: rho0MTsh(:, :, :, :)
    complex,           allocatable              :: grRho0IR(:, :)
    complex,           allocatable              :: rho0IRpw(:, :)
    complex                                     :: psiGrVeffPsiUC(3)
    real                                        :: Gext(3)
    complex                                     :: cntrbMatElem12(3, 3)
    complex                                     :: cntrbMatElem3(3, 3)


    ! This is a test for q = 0
    iqpt = 1

    ! Calculation of the Equation ?, left-hand side.
    ! -----------------------------------------------

    ! One can choose whether to calculate the numerical gradient of the unperturbed potential or to use the Weinert method for
    ! calculating it. The Weinert method has the advantage that continuity is enforced by construction at the MT boundary.
    ! It has been tested that both methods deliver the same numbers except for the just named behavior.
    numericalGradient = .true.
    if (numericalGradient) then

      ! Calculate the numerical gradient of the unperturbed effective potential.
      ! ------------------------------------------------------------------------

      ! Calculate the IR gradient.
      !   Convert the unperturbed effective potential from star expansion coefficients to plane-wave expansion coefficients
      allocate( vEff0IRpw(ngdp) )
      vEff0IRpw(:) = cmplx(0., 0.)
      ! SPIN MISSING
      call ConvertStar2G( vEff0IRst(:, 1), vEff0IRpw(:), stars, ngdp, gdp )

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
            r2Veff0MT(imesh, ilh, itype) = vEff0MTLh(imesh, ilh, itype, 1) * atoms%rmsh(imesh, itype) * atoms%rmsh(imesh, itype)
          end do ! imesh
        end do ! ilh
      end do ! itype

      call CalcGrR2FinLH( atoms, lathar, clnu_atom, nmem_atom, mlh_atom, r2Veff0MT, r2GrVeff0MT )

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
                  grVeff0MT(imesh, lm, iatom, idir) = r2GrVeff0MT(imesh, lm, iatom, idir) / atoms%rmsh(imesh, itype)**2
                end do ! imesh
              end do ! mqn_m
            end do ! oqn_l
          end do ! ieqat
        end do ! itype
      end do ! idir

    else

      ! Use the Weinert method to calculate the gradient of the unperturbed effective potential
      ! ---------------------------------------------------------------------------------------

      ! Calculate the gradient of the unperturbed density.
      !   Any factors as in Fleur have been removed. In order to improve the numerical accuracy of the gradient r^2 is multiplied
      !   before and removed after the calculation of the gradient. Thereyby, subtraction of small numbers close to the core are
      !   avoided.
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

      ! Precalculate quantities for the kernel of the unperturbed effective potential's gradient.
      ! SPIN MISSING
      call CalcIRdVxcKern(stars, gdp, ngdp, rho0IR(:, 1), grVxcIRKern)
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

    end if ! numericalGradient

    ! Warp the interstitial gradient of the unperturbed effective potential
    allocate(w_grVeff0IR(ngdp, 3))
    w_grVeff0IR = cmplx(0.0, 0.0)
    do idir = 1, 3
      call WarpIRPot(stars, ngdp, idir, gdp, grVeff0IR, w_grVeff0IR(:, idir))
    end do ! idir

    allocate( w_vEff1IRContainer(ngdp, 3, atoms%nat) )
    w_vEff1IRContainer(:, :, :) = cmplx(0., 0.)

    allocate( sumVMTsContainer(atoms%jmtd, (atoms%lmaxd + 1 )**2, 3, atoms%nat) )
    allocate( sumVMTs2Container(atoms%jmtd, (atoms%lmaxd + 1 )**2, 3, atoms%nat) )
    sumVMTsContainer(:, :, :, :) = cmplx(0., 0.)
    sumVMTs2Container(:, :, :, :) = cmplx(0., 0.)

    ! Deactivate symmetry operations to rotate MTs of non-representative atoms.
    allocate(ngoprI(atoms%nat))
    ngoprI(:) = 1

    psiGrVeffPsiUC = cmplx(0., 0.)

    ! Allocate and initialize noco arrays and deactivate noco calculations and spin spiral calculations
    allocate(noco%alph(atoms%ntype), noco%beta(atoms%ntype))
    noco%alph(:) = 1.
    noco%beta(:) = 1.
    noco%qss(:)  = 0.
    noco%l_noco  = .false.
    noco%l_ss    = .false.


    iDatom = 0
    do iDtype = 1, atoms%ntype
      do iDeqat = 1, atoms%neq(iDtype)
        iDatom = iDatom + 1


        ! We use routines expecting the first variation of the effective potential so as well a selective displacement of single
        ! atoms iDatom. Therefore, we use containers and put the negative gradient of the unperturbed effective potential into these
        ! because this should be the case for q = 0 from an analytical point of view.
        w_vEff1IRContainer(1:ngdp, 1:3, iDatom) = -w_grVeff0IR(1:ngdp, 1:3)

        sumVMTsContainer = cmplx(0.0, 0.0)
        sumVMTs2Container = cmplx(0.0, 0.0)
        iatom = 0
        ! We choose here between two version because the order of indices for the potential and gradient routines is still
        ! inconsistent.
        ! TODO ISSUE NO 7: INCONSISTENCY OF ARRAY INDEX ORDER, IF-THEN CAN THEN BE REMOVED
        if (numericalGradient) then
          do itype = 1, atoms%ntype
            do ieqat = 1, atoms%neq(itype)
              iatom = iatom + 1
              do idir = 1, 3
                do oqn_l = 0, atoms%lmax(itype)
                  lm_pre = oqn_l * (oqn_l + 1) + 1
                  do mqn_m = -oqn_l, oqn_l
                    lm = lm_pre + mqn_m
                    lm2 = lm_pre - mqn_m
                    do imesh = 1, atoms%jri(itype)
                      sumVMTsContainer(imesh, lm, idir, iatom) = -grVeff0MT(imesh, lm, iatom, idir)
                      sumVMTs2Container(imesh, lm2, idir, iatom) = (-1)**(mqn_m) * conjg(sumVMTsContainer(imesh, lm, idir, iatom))
                    end do ! imesh
                  end do ! mqn_m
                end do ! oqn_l
              end do ! idir
            end do ! ieqat
          end do ! itype
        else
          do itype = 1, atoms%ntype
            do ieqat = 1, atoms%neq(itype)
              iatom = iatom + 1
              do idir = 1, 3
                do oqn_l = 0, atoms%lmax(itype)
                  lm_pre = oqn_l * (oqn_l + 1) + 1
                  do mqn_m = -oqn_l, oqn_l
                    lm = lm_pre + mqn_m
                    lm2 = lm_pre - mqn_m
                    do imesh = 1, atoms%jri(itype)
                      sumVMTsContainer(imesh, lm, idir, iatom) = -grVeff0MT(imesh, lm, idir, iatom)
                      sumVMTs2Container(imesh, lm2, idir, iatom) = (-1)**(mqn_m) * conjg(sumVMTsContainer(imesh, lm, idir, iatom))
                    end do ! imesh
                  end do ! mqn_m
                end do ! oqn_l
              end do ! idir
            end do ! ieqat
          end do ! itype
        end if

        ! Calculate the k-independent part of the MT potential matrix element <BasFunc| -grVeff0 |BasFunc>  As the potential, has a
        ! complex part and we want, in a first step, obey to the Fleur routine structure of tlmplm and hssr_wu, we need a seperate
        ! treatment of the upper and lower half of these tlmplm matrices.
        do idir = 1, 3
          ! SPIN MISSING
          call Tlmplm4V( atoms, lathar, dimens, enpara, usdus, input, tdV(idir), 1, 1, idir, sumVMTsContainer, rbas1, rbas2,      &
                                                                                      & uuilon, duilon, ulouilopn, ilo2p, nlo_atom )

          ! SPIN MISSING
          call Tlmplm4V( atoms, lathar, dimens, enpara, usdus, input, tdV2(idir), 1, 1, idir, sumVMTs2Container, rbas1, rbas2,     &
                                                                                      & uuilon, duilon, ulouilopn, ilo2p, nlo_atom )
        end do


        do ikpt = 1, kpts%nkpt

          ! SPIN MISSING
          nrBndBra = nobd(ikpt, 1)
          nrBndKet = nobd(ikpt, 1)

          allocate( psiGrVeffPsiIR(nrBndBra, nrBndKet) )
          allocate( psiGrVeffPsiMT(nrBndBra, nrBndKet ) )

          do idir = 1, 3

            ! Calculate the intersitital part of <Basfun|-grVeff0|Basfun>
            ! -----------------------------------------------------------
            ! nmat has to be the dimension of the bras due to the concept of CalcMEPotIR
            nmat = nv(1, ikpt) + atoms%nlotot
            psiGrVeffPsiIR = cmplx(0., 0.)
            ! SPIN MISSING
            call CalcMEPotIR( stars, dimens, gbas(:, mapGbas(:nv(1, ikpt), ikpt, 1)), &
              & gbas(:, mapGbas(:nv(1, ikpt), ikpt, 1)), nv, w_vEff1IRContainer(:, idir, iDatom), z0(:, :, ikpt, 1),        &
              & z0(:, :, ikpt, 1), gdp, nmat, nrBndBra, nrBndKet, ikpt, iqpt, ikpt, ngdp, psiGrVeffPsiIR, kpq2kPrVec )


            ! Calculate the muffin-tin part of <Basfun|-grVeff0|Basfun>
            ! Call of abcof for standard kets of MT potentials matrix element
            allocate(acof(nrBndKet, 0:dimens%lmd, atoms%nat), bcof(nrBndKet, 0:dimens%lmd, atoms%nat), &
              &ccof(-atoms%llod:atoms%llod, nrBndKet, atoms%nlod, atoms%nat))
            acof(:, :, :) = cmplx(0.0, 0.0)
            bcof(:, :, :) = cmplx(0.0, 0.0)
            ccof(:, :, :, :) = cmplx(0.0, 0.0)
            ! SPIN MISSING
            call Abcof ( atoms%lmaxd, atoms%ntype, dimens%neigd, nrBndKet, atoms%nat, sym%nop, dimens%nvd, dimens%jspd, dimens%lmd, &
              & dimens%nbasfcn, atoms%llod, atoms%nlod, atoms%nlotot, sym%invtab, atoms%ntype, sym%mrot, ngoprI, atoms%taual, atoms%neq, &
              & atoms%lmax, atoms%rmt, cell%omtil, cell%bmat, cell%bbmat, kpts%bk(:, ikpt), &
              & gbas(1, mapGbas(:nv(1, ikpt), ikpt, 1)), gbas(2, mapGbas(:nv(1, ikpt), ikpt, 1)), &
              & gbas(3, mapGbas(:nv(1, ikpt), ikpt, 1)), nv(:, ikpt), nmat, nrBndKet, z0(:, :, ikpt, 1), usdus%us(:, :, 1), &
              & usdus%dus(:, :, 1), usdus%uds, usdus%duds(:, :, 1), usdus%ddn(:, :, 1), atoms%invsat, sym%invsatnr, usdus%ulos(:, :, 1), &
              & usdus%uulon(:, :, 1), usdus%dulon(:, :, 1),  usdus%dulos(:, :, 1), atoms%llo, atoms%nlo, atoms%l_dulo, &
              & atoms%lapw_l, noco%l_noco, noco%l_ss, 1, noco%alph, noco%beta, noco%qss, kveclo(:, ikpt), odi, ods, &
              & acof, bcof, ccof)


            ! Rearrange potential acofs, bcofs, ccofs ensuring an consistent handling of LO terms within the loop the structure compared to LAPWs
            maxlmp = maxval( (/ (sum( (/ ((2*oqn_l+1)* nRadFun(oqn_l, iDtype), oqn_l = 0, atoms%lmax(iDtype)) /) ), iDtype=1, atoms%ntype) /) )
            allocate(mCoef(nrBndKet, maxlmp, atoms%nat))
            iatom  = 0
            mCoef =  0
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
                    mCoef(:nrBndKet, lmp, iatom) = acof(:nrBndKet, lm, iatom)
                    !p = 2
                    lmp = lmp + 1
                    mCoef(:nrBndKet, lmp, iatom) = bcof(:nrBndKet, lm, iatom)
                    !LOs
                    do iradf = 3, nRadFun(oqn_l, itype)
                      ilo = iloTable(iradf, oqn_l, itype)
                      lmp = lmp + 1
                      mCoef(:nrBndKet, lmp, iatom) = ccof(mqn_m, :nrBndKet, ilo, iatom)
                    end do
                  end do
                end do
              end do
            end do
            deallocate(acof, bcof, ccof)

            ! Calculate the matrix element of the MT potentials
            psiGrVeffPsiMT = cmplx(0.0, 0.0)
            ! SPIN MISSING
            call CalcVsumMT( atoms, tdV(idir), tdV2(idir), ikpt, ikpt, nobd(:, 1), nobd, mCoef, mCoef, nRadFun, &
              & iloTable, nlo_atom, psiGrVeffPsiMT )

            deallocate( mCoef )

            ! Sum up the k-dependent contributions and only take the entries diagonal in bands
            do iband = 1, nrBndKet
              psiGrVeffPsiUC(idir) = psiGrVeffPsiUC(idir) + psiGrVeffPsiIR(iband, iband) + psiGrVeffPsiMT(iband, iband)
            end do ! iband

          end do ! idir

          deallocate( psiGrVeffPsiIR, psiGrVeffPsiMT )

        end do ! ikpt

      end do ! iDeqat
    end do ! iDtype

    ! Calculation of the Equation ?, right-hand side.
    ! -----------------------------
    ! Quantities for initialization
    allocate( lmpT(atoms%ntype) )
    lmpT(:) = 0
    do iDtype = 1, atoms%ntype
      lmpT(iDtype) = sum( [ ( (2 * oqn_l + 1)* nRadFun(oqn_l, iDtype), oqn_l = 0, atoms%lmax(iDtype) ) ] )
    end do ! iDtype
    lmpMax     = maxval( lmpT(:) )
    nRadFunMax = maxval( nRadFun(:, :) )

    allocate( varphi1(atoms%jmtd, nRadFunMax, 0:atoms%lmaxd), varphi2(atoms%jmtd, nRadFunMax, 0:atoms%lmaxd) )
    allocate( r2(atoms%jmtd) )
    allocate( delrVarphi1( atoms%jmtd, nRadFunMax, 0:atoms%lmaxd), delrVarphi2( atoms%jmtd, nRadFunMax, 0:atoms%lmaxd) )
    allocate( delrGrVarphiCh1(atoms%jmtd, 2, lmpMax, -1:1), delrGrVarphiCh2(atoms%jmtd, 2, lmpMax, -1:1) )
    allocate( grVarphiCh1(atoms%jmtd, 2, lmpMax, -1:1), grVarphiCh2(atoms%jmtd, 2, lmpMax, -1:1), &
            & grVarphiChLout(2, 0:atoms%lmaxd), grVarphiChMout(-atoms%lmaxd:atoms%lmaxd, -1:1) )
    allocate( hVarphi(2, atoms%jmtd, 0:(atoms%lmaxd + 3)**2 - 1, lmpMax))
    allocate( varphiVarphi(lmpMax, lmpMax, atoms%ntype), varphiHvarphi(lmpMax, lmpMax, atoms%nat))
    allocate( varphiGrVarphi(lmpMax, lmpMax, -1:1, atoms%ntype), grVarphiVarphi(lmpMax, lmpMax, -1:1, atoms%ntype), grVarphiHVarphi(lmpMax, lmpMax, -1:1, atoms%nat ) )
    allocate( vEff0MtSpH( atoms%jmtd, 0:dimens%lmd))
    allocate( vEff0NsphGrVarphi(2, atoms%jmtd, 0:(atoms%lmaxd + 3)**2 - 1, lmpMax, -1:1), &
                                                     & r2grVeff0SphVarphi(2, atoms%jmtd, 0:(atoms%lmaxd + 3)**2 - 1, lmpMax, 3), &
                                                     & r2grVeff0Varphi(2, atoms%jmtd, 0:(atoms%lmaxd + 3)**2 - 1, lmpMax, 3))
    allocate( hGrVarphi(2, atoms%jmtd, (atoms%lmaxd + 1)**2, lmpMax, -1:1))
    allocate( varphiVarphiDummy(lmpMax, lmpMax, atoms%ntype) )
    allocate( varphiGrVeff0SphVarphi(lmpMax, lmpMax, atoms%nat, 3) )
    allocate( varphiGrVeff0Varphi(lmpMax, lmpMax, atoms%nat, 3) )
    allocate( varphiHGrvarphi(lmpMax, lmpMax, atoms%nat, -1:1) )

    hGrVarphi(:, :, :, :, :) = cmplx(0., 0.)
    varphiGrVeff0SphVarphi(:, :, :, :) = cmplx(0., 0.)
    varphiGrVeff0Varphi(:, :, :, :) = cmplx(0., 0.)
    varphiVarphiDummy(:, :, :) = 0.
    varphiHGrvarphi(:, :, :, :) = cmplx(0., 0.)
    varphiVarphi(:, :, :)                         = 0.
    varphiHvarphi(:, :, :)                        = cmplx(0., 0.)
    grVarphiVarphi(:, :, :, :)                    = 0.
    varphiGrVarphi(:, :, :, :)                    = 0.
    grVarphiHVarphi(:, :, :, :)                   = cmplx(0., 0.)
    vEff0NsphGrVarphi(:, :, :, :, :)         = cmplx(0., 0.)
    r2grVeff0SphVarphi(:, :, :, :, :)        = cmplx(0., 0.)
    r2grVeff0Varphi(:, :, :, :, :)        = cmplx(0., 0.)

    iDatom = 0
    do iDtype = 1, atoms%ntype

      ! Precalculate radial Jacobi determinant for later integrals
      r2(:) = 0.
      do imesh = 1, atoms%jri(iDtype)
        r2(imesh) = atoms%rmsh(imesh, iDtype) * atoms%rmsh(imesh, iDtype)
      end do ! imesh

      ! Prepare the radial solutions u and their derivatives
      varphi1(:, :, :) = 0.
      varphi2(:, :, :) = 0.
      delrVarphi1(:, :, :) = 0.
      delrVarphi2(:, :, :) = 0.
      do oqn_l = 0, atoms%lmax(iDtype)
        do iradf = 1, nRadFun(oqn_l, iDtype)
          do imesh = 1, atoms%jri(iDtype)
            ! In Fleur the radial solutions u_lp are multiplied by a factor r by default to avoid an additional multiplication of the
            ! Jacobi determinant r^2 in radial integrals given spherical coordinates.
            varphi1(imesh, iradf, oqn_l) = rbas1(imesh, iradf, oqn_l, iDtype, 1) / atoms%rmsh(imesh, iDtype)
            varphi2(imesh, iradf, oqn_l) = rbas2(imesh, iradf, oqn_l, iDtype, 1) / atoms%rmsh(imesh, iDtype)
          end do ! imesh
          ! Precalculate partial derivatives of varphis in r-direction since it is needed twice
          call Derivative( varphi1(1:atoms%jri(iDtype), iradf, oqn_l), iDtype, atoms, delrVarphi1(1:atoms%jri(iDtype), iradf, oqn_l) )
          call Derivative( varphi2(1:atoms%jri(iDtype), iradf, oqn_l), iDtype, atoms, delrVarphi2(1:atoms%jri(iDtype), iradf, oqn_l) )
        end do ! iradf
      end do ! oqn_l

      ! Calculate the application of the gradient and the gradient's dyadic product onto the MT basis functions (matching coefficients
      ! have no spatial dependence) and determine its scattering channels.
      grVarphiChLout(:, :) = 0
      grVarphiChMout(:, :) = 0
      grVarphiCh1(:, :, :, :) = 0.
      grVarphiCh2(:, :, :, :) = 0.
      call CalcChannelsGrFlpNat( atoms, iDtype, nRadFun, varphi1, varphi2, delrVarphi1, delrVarphi2, grVarphiChLout, grVarphiChMout, &
                                                                                                        & grVarphiCh1, grVarphiCh2 )

      ! Precalculate partial derivatives of the varphis' gradients in r-direction since it is required two times. The grGrt version
      ! is only set zero and allocated because it is needed as an API parameter later.
      delrGrVarphiCh1(:, :, :, :) = 0.
      delrGrVarphiCh2(:, :, :, :) = 0.
      do mqn_m2PrR = -1, 1
        lmp = 0
        do oqn_l = 0, atoms%lmax(iDtype)
          do mqn_m = -oqn_l, oqn_l
            do iradf = 1, nRadFun(oqn_l, iDtype)
              lmp = lmp + 1
              do ichan = 1, 2
                call Derivative( grVarphiCh1(1:atoms%jri(iDtype), ichan, lmp, mqn_m2PrR), iDtype, atoms, delrGrVarphiCh1(1:atoms%jri(iDtype), ichan, lmp, mqn_m2PrR) )
                call Derivative( grVarphiCh2(1:atoms%jri(iDtype), ichan, lmp, mqn_m2PrR), iDtype, atoms, delrGrVarphiCh2(1:atoms%jri(iDtype), ichan, lmp, mqn_m2PrR) )
              end do ! ichan
            end do !iradf
          end do ! mqn_m
        end do ! oqn_l
      end do ! mqn_m2PrR


      ! Calculate H |varphi>
      do iDeqat = 1, atoms%neq(iDtype)
        iDatom = iDatom + 1

        ! Expand the coefficients of the lattice-harmonic given potential into spherical-harmonic coefficients for the given atom.
        vEff0MtSpH(:, :) = cmplx(0.0, 0.0)
        ptsym = atoms%ntypsy(iDatom)
        do ilh = 0, lathar%nlh(ptsym)
          oqn_l = lathar%llh(ilh, ptsym)
          lm_pre = oqn_l * (oqn_l + 1)
          do imem = 1, nmem_atom(ilh, iDatom)
            mqn_m = mlh_atom(imem, ilh, iDatom)
            lm = lm_pre + mqn_m
            do imesh = 1, atoms%jri(iDtype)
              ! SPIN MISSING
              vEff0MtSpH(imesh, lm) = vEff0MtSpH(imesh, lm) + vEff0MtLh(imesh, ilh, iDtype, 1) * clnu_atom(imem, ilh, iDatom)
            end do ! imesh
          end do ! imem
        end do ! ilh

        ! Calculate the action of the Hamiltonian onto a basis function.
        ! For the integrals we have to go until lmax + 2, but later, we only have to multiply the abcof until lmax
        hVarphi = cmplx(0.0, 0.0)
        call CalcHnGrV0Varphi( atoms, lathar, iDtype, iDatom, lmpMax, El, varphi1, varphi2, nRadFun, vEff0MtSpH, vEff0MtLh(:, :, :, 1), clnu_atom, &
          & nmem_atom, mlh_atom, grVarphiCh1, grVarphiCh2, grVarphiChLout, grVarphiChMout, hVarphi, vEff0NsphGrVarphi,             &
          & r2grVeff0SphVarphi, r2GrVeff0Varphi )

        ! Calculate all radial integrals with no gradients.
        call CalcScalBasfMatElems( atoms, iDtype, iDatom, nRadFun, r2, varphi1, varphi2, hVarphi, varphiVarphi, varphiHvarphi )

        ! Calculate all radial integrals with one gradient.
        call CalcVecBasfMatElems( atoms, iDtype, 2, nRadFun, r2, grVarphiChLout, grVarphiChMout, varphi1, varphi2, &
          & grVarPhiCh1, grVarphiCh2, hVarphi, grVarphiVarphi(:, :, :, iDtype),  grVarphiHVarphi(:, :, :, iDatom) )

        do mqn_m2PrC = -1, 1
          call CalcHGrVarphi( atoms, iDtype, mqn_m2PrC, lmpMax, atoms%lmax(iDtype), grVarphiChMout, nRadFun, grVarPhiCh1, grVarPhiCh2,        &
                                                      & grVarphiChLout, vEff0NsphGrVarphi, El, lmpT, hGrVarphi )
        end do ! mqn_m2PrC
        do mqn_m2PrC = -1, 1
          varphiVarphiDummy(:, :, :) = 0.
          !todo attention where the hGrVarphi starts!
          call CalcScalBasfMatElems( atoms, iDtype, iDatom, nRadFun, r2, varphi1, varphi2, hGrVarphi(:, :, :, :, mqn_m2PrC), varphiVarphiDummy, varphiHGrvarphi(:, :, :, mqn_m2PrC) )
          do jj = 1, lmpMax
            do ii = 1, lmpMax
              varphiGrVarphi(ii, jj, mqn_m2PrC, iDtype) = grVarphiVarphi(jj, ii, mqn_m2PrC, iDtype)
            end do ! ii
          end do ! jj
        end do ! mqn_m2PrC

      end do ! iDeqat
    end do !iDtype

    iatom = 0
    r2(:) = 1.
    iatom = 0
    do itype = 1, atoms%ntype
      do ieqat = 1, atoms%neq(itype)
        iatom = iatom + 1
        do idir = 1, 3
          call CalcScalBasfMatElems( atoms, itype, iatom, nRadFun, r2, varphi1, varphi2, r2GrVeff0SphVarphi(:, :, :, :, idir), varphiVarphiDummy, varphiGrVeff0SphVarphi(:, :, :, idir) )
          call CalcScalBasfMatElems( atoms, itype, iatom, nRadFun, r2, varphi1, varphi2, r2GrVeff0Varphi(:, :, :, :, idir), varphiVarphiDummy, varphiGrVeff0Varphi(:, :, :, idir) )
        end do ! idir
      end do ! ieqat
    end do ! itype

    ! get rid of unrequired arrays
    deallocate( vEff0MtSpH, r2, hVarphi, r2grVeff0SphVarphi, vEff0NsphGrVarphi )
    deallocate( varphi1, varphi2, grVarphiCh1, grVarphiCh2, grVarphiChLout, grVarphiChMout )
    deallocate( delrVarphi1, delrVarphi2, delrGrVarphiCh1, delrGrVarphiCh2 )

    eps1DynMatPulTestSw = .true.
    testComp2ndN1stOrdBasFuncSw = .false.
    dynMatPu3rdBraKetHepsSw = .false.
    cntrbMatElem12(:, :) = cmplx(0., 0.)
    cntrbMatElem3(:, :) = cmplx(0., 0.)

    ! This loop can be parallelized later
    allocate( z1nG(dimens%nbasfcn, 3, atoms%nat, maxval(nobd(:, :))) )
    do ikpt = 1, kpts%nkpt

      z1nG = cmplx(0.0, 0.0)
      do iband = 1, nobd(ikpt, 1)
        do idir = 1, 3
          do iBas = 1, nv(1, ikpt) + atoms%nlotot
            z1nG(iBas, idir, 1, iband) = 0.5 * z0(iBas, iband, ikpt, 1)
          end do ! iBas
        end do ! idir
      end do ! iband

      call Add2ndOrdWfPulayBraKets2DynMat( atoms, cell, dimens, sym, kpts, qpts, usdus, results, iqpt, ikpt, nRadFunMax, lmpMax, &
        & eps1DynMatPulTestSw, testComp2ndN1stOrdBasFuncSw, mapKpq2K, eig, nobd, nv, &
        & gbas, mapGbas, kveclo, z0, z1nG, nRadFun, iloTable, varphiVarphi, varphiHvarphi, grVarphiVarphi, grVarphiHVarphi, &
        & lmpT, varphiGrVeff0SphVarphi, varphiHGrvarphi, varphiGrVarphi, varphiGrVeff0Varphi, cntrbMatElem12 )

      call Add1stOrdWfPulayBraKets2DynMat( atoms, kpts, qpts, sym, dimens, cell, usdus, stars, results, ikpt, iqpt, lmpMax, nRadFunMax, nv, gbas, mapGbas, kveclo, &
        & mapKpq2K, nobd, z0, z1nG, iloTable, grVarphiVarphi, nRadFun, eig, kpq2kPrVec, grVarphiHvarphi, &
        & varphiVarphi, varphiHvarphi, w_Veff0IRst, lmpT, &
        & eps1DynMatPulTestSw, testCompTerm3rdBraKetsVarBra, testCompTerm3rdBraKetsVarKet, dynMatPu3rdBraKetHepsSw, varphiGrVeff0SphVarphi, varphiHGrvarphi, varphiGrVarphi, cntrbMatElem3 )

    end do ! ikpt

    ! Does the alogorithm pass the test?
    passed1st2nd = .true.
    passed3rd = .true.

    do idir = 1, 3

      if ( abs( psiGrVeffPsiUC(idir) - cntrbMatElem12(idir, idir) ) >= 1e-6 ) then
        write (*, '(a)')   'Test comparing eps(1) with parts of the Pulay contribution matrix elements to the dynamical matrix failed!'
        write(*, '(a)') 'Direction and left-hand side contribution 1 and 2'
        write(*, '(i8,2f15.8)') idir, psiGrVeffPsiUC(idir)
        write(*, '(a)') 'Right-hand side. Contribution 1 and 2, only diagonal relevant!'
        write(*, '(3(2(es16.8,1x),3x))') cntrbMatElem12(1, :)
        write(*, '(3(2(es16.8,1x),3x))') cntrbMatElem12(2, :)
        write(*, '(3(2(es16.8,1x),3x))') cntrbMatElem12(3, :)
        passed1st2nd = .false.
      end if

      if ( abs( psiGrVeffPsiUC(idir) - cntrbMatElem3(idir, idir) ) >= 1e-6 ) then
        write (*, '(a)')   'Test comparing eps(1) with parts of the Pulay contribution matrix elements to the dynamical matrix failed!'
        write(*, '(a)') 'Direction and left-hand side contribution 3'
        write(*, '(i8,2f15.8)') idir, psiGrVeffPsiUC(idir)
        write(*, '(a)') 'Right-hand side. Contribution 3, only diagonal relevant!'
        write(*, '(3(2(es16.8,1x),3x))') cntrbMatElem3(1, :)
        write(*, '(3(2(es16.8,1x),3x))') cntrbMatElem3(2, :)
        write(*, '(3(2(es16.8,1x),3x))') cntrbMatElem3(3, :)
        passed3rd = .false.
      end if

    end do ! idir

    passed = passed1st2nd .and. passed3rd
    if (.not.passed) then
      call juDFT_warn('Compare eps(1) with parts of the Pulay contribution matrix elements to the dynamical matrix failed', &
        & calledby='TestDynMatPuWithEps1', hint='Debug.')
    end if
  end subroutine TestDynMatPuWithEps1

  !>-------------------------------------------------------------------------------------------------------------------------------------
  !> @author
  !> Christian-Roman Gerhorst and Markus Betzinger
  !>
  !> @brief
  !> Set z1, differently, such that within the 1st, 2nd and 3rd Pulay braket of the dynamical matrix, terms cancel each other.
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
  !> @param[out] dimens      : Dimension type, see types.f90.
  !> @param[out] qpts        : Q-points type, see types.f90.
  !> @param[out] usdus       : Type containing quantities consisting of the radial solutions, see types.f90.
  !> @param[out] results     : Results type, see types.f90.
  !> @param[in]  logUnit     : Unit number for juPhon.log.
  !> @param[out] kveclo      : Basis G-vectors of local orbitals.
  !> @param[out] nv          : Number of LAPW G-basis vectors for given k-point.
  !> @param[out] nRadFun     : Number of radial functions per orbital quantum number l and atom type.
  !> @param[out] iloTable    : Number of local orbital if orbital quantum number l, atom type and index p > 2 from nRadFun is given.
  !> @param[out] GbasVec     : G-basis vectors
  !> @param[out] ilst        : For various k-points G-basis vectors occur more than once, thus they are only stored once in juPhon. This
  !>                           pointer array contains the right index for GbasVec array to "unfold" G-basis vectors again.
  !> @param[out] mapKpq2K    : For a given k-point index and q-point index, this array gives the k-point set index of the result k + q
  !>                           (already mapped back to Brillouin zone).
  !> @param[out] eig         : Contains Kohn\--Sham eigenvalues.
  !> @param[out] rbas1       : Large components of radial solution, its energy derivative and u_LO
  !> @param[out] rbas2       : Small components of radial solution, its energy derivative and u_LO
  !> @param[out] El          : Contains LAPW and LO energy parameters.
  !> @param[out] z           : Kohn-Sham eigenvectors.
  !>                           potential parsed from Fleur.
  !> @param[out] mlh_atom    : Magnetic quantum number m of lattice harmonic members for every atom.
  !> @param[out] nmem_atom   : Number of lattice harmonic members for every atom.
  !> @param[out] clnu_atom   : Phase mediating between stars and plane waves.
  !> @param[out] nobd        : Number of occupied bands per k-point and spin
  !>
  !> @note : This test for the 3rd braket only works if the matrix elements with one basis variation in ket or bra are zero by
  !> default (for example for neon)
  !> Tests whether for certain choice of z1 certain contributions in the 1st, 2nd and 3rd raket of the dynamical matrix cancel
  !>--------------------------------------------------------------------------------------------------------------------------------------
  subroutine TestComp2ndN1stOrdBasFunc( atoms, lathar, dimens, kpts, cell, sym, qpts, usdus, stars, results, testComp2ndN1stOrdBasFuncSw,   &
      & testCompTerm3rdBraKets, logUnit, nRadFun, mlh_atom, nmem_atom, clnu_atom, vEff0MTLh, rbas1, rbas2, El, nobd, nv, z, mapKpq2K, eig,  &
      & GbasVec, ilst, kveclo, iloTable, kpq2kPrVec, vEff0IR )

    use m_types
    use m_jpSetupDynMat, only : CalcScalBasfMatElems, CalcVecBasfMatElems, Add2ndOrdWfPulayBraKets2DynMat, Add1stOrdWfPulayBraKets2DynMat
    use mod_juPhonUtils, only : Derivative, CalcChannelsGrFlpNat
    use m_jpConstants, only : iu
    use m_jpSetupDynMatHelper, only : CalcHnGrV0Varphi
    use m_juDFT_NOstopNO, only : juDFT_warn
    use m_jpSetupDynMatSF, only : CalcHGrVarphi

    implicit none

    ! Type parameters
    type(t_atoms),                  intent(in) :: atoms
    type(t_sphhar),                 intent(in) :: lathar
    type(t_dimension),              intent(in) :: dimens
    type(t_kpts),                   intent(in) :: kpts
    type(t_cell),                   intent(in) :: cell
    type(t_sym),                    intent(in) :: sym
    type(t_kpts),                   intent(in) :: qpts
    type(t_usdus),                  intent(in) :: usdus
    type(t_stars),                  intent(in) :: stars
    type(t_results),                intent(in) :: results

    ! Scalar parameters
    logical,                        intent(in) :: testComp2ndN1stOrdBasFuncSw
    logical,                        intent(in) :: testCompTerm3rdBraKets
    integer,                        intent(in) :: logUnit

    ! Array parameters
    integer,                        intent(in) :: nRadFun(0:, :)
    integer,                        intent(in) :: mlh_atom(:,0:,:)
    integer,                        intent(in) :: nmem_atom(0:, :)
    complex,                        intent(in) :: clnu_atom(:,0:,:)
    real,                           intent(in) :: vEff0MTLh(:, 0:, :)
    real,                           intent(in) :: rbas1(:, :, 0:, :, :)
    real,                           intent(in) :: rbas2(:, :, 0:, :, :)
    real,                           intent(in) :: El(:, 0:, :, :)
    integer,                        intent(in) :: iloTable(:, 0:, :)
    integer,                        intent(in) :: nobd(:, :)
    integer,                        intent(in) :: nv(:, :)
    MCOMPLEX,                       intent(in) :: z(:, :, :, :) ! Attention: see zBar
    integer,                        intent(in) :: mapKpq2K(:, :)
    real,                           intent(in) :: eig(:, :, :)
    integer,                        intent(in) :: GbasVec(:, :)
    integer,                        intent(in) :: ilst(:, :, :)
    integer,                        intent(in) :: kveclo(:,:)
    integer,                        intent(in) :: kpq2kPrVec(:, :, :)
    complex,                        intent(in) :: vEff0IR(:,:)

    ! Scalar variables
    integer                                    :: itype
    integer                                    :: oqn_l
    integer                                    :: lmpMax
    integer                                    :: nRadFunMax
    integer                                    :: iatom
    integer                                    :: imesh
    integer                                    :: iradf
    integer                                    :: mqn_m2PrR !todo do we really need variables for R and C???
    integer                                    :: lmp
    integer                                    :: mqn_m
    integer                                    :: ichan
    integer                                    :: ieqat
    integer                                    :: ptsym
    integer                                    :: ilh
    integer                                    :: lm_pre
    integer                                    :: imem
    integer                                    :: mqn_m2PrC
    integer                                    :: lm
    integer                                    :: iband
    integer                                    :: idir
    integer                                    :: iBas
    logical                                    :: eps1DynMatPulTestSw
    integer                                    :: ikpt
    integer                                    :: iqpt
    integer                                     :: ii
    integer                                     :: jj
    logical                                    :: testCompTerm3rdBraKetsVarBra
    logical                                    :: testCompTerm3rdBraKetsVarKet
    logical                                    :: dynMatPu3rdBraKetHepsSw
    logical                                    :: passed = .true.

    ! Array variables
    integer,           allocatable             :: lmpT(:)
    real,              allocatable             :: varphi1(:, :, :)
    real,              allocatable             :: varphi2(:, :, :)
    real,              allocatable             :: r2(:)
    real,              allocatable             :: delrVarphi1(:, :, :)
    real,              allocatable             :: delrVarphi2(:, :, :)
    real,              allocatable             :: delrGrVarphiCh1(:, :, :, :)
    real,              allocatable             :: delrGrVarphiCh2(:, :, :, :)
    real,              allocatable             :: grVarphiCh1(:, :, :, :)
    real,              allocatable             :: grVarphiCh2(:, :, :, :)
    integer,           allocatable             :: grVarphiChLout(:, :)
    integer,           allocatable             :: grVarphiChMout(:, :)
    complex,           allocatable             :: hVarphi(:, :, :, :)
    complex,           allocatable             :: vnsphEff0Varphi(:, :, :, :)
    complex,           allocatable             :: grVnsphEff0Varphi(:, :, :, :, :)
    real,              allocatable             :: varphiVarphi(:, :, :)
    complex,           allocatable             :: varphiHvarphi(:, :, :)
    real,              allocatable             :: grVarphiVarphi(:, :, :, :)
    real,              allocatable             :: varphiGrVarphi(:, :, :, :)
    complex,           allocatable             :: grVarphiHVarphi(:, :, :, :)
    complex,           allocatable             :: grVarphiVnsphVarphiSurf(:, :)
    complex,           allocatable             :: vEff0NsphGrVarphi(:, :, :, :, :)
    complex,           allocatable             :: r2grVeff0SphVarphi(:, :, :, :, :)
    complex,           allocatable             :: r2grVeff0Varphi(:, :, :, :, :)
    complex,           allocatable             :: vEff0MtSpH(:, :)
    complex,           allocatable             :: z1nG(:, :, :, :)
    complex,       allocatable              :: hGrVarphi(:, :, :, :, :)
    real,          allocatable              :: varphiVarphiDummy(:, :, :)
    complex,       allocatable              :: varphiGrVeff0SphVarphi(:, :, :, :)
    complex,       allocatable              :: varphiGrVeff0Varphi(:, :, :, :)
    complex,       allocatable              :: varphiHGrvarphi(:, :, :, :)
    complex                                    :: pulayCntrb(3, 3)
    real                                       :: kExt(3)
    real                                       :: gbasExt(3)

    write (logUnit, '(a)')   'Test Goldstone condition for the Pulay terms.'
    write (logUnit, '(a)' ) '---------------------------------------------'

    ! Quantities for initialization
    allocate( lmpT(atoms%ntype) )
    lmpT(:) = 0
    do itype = 1, atoms%ntype
      lmpT(itype) = sum( [ ( (2 * oqn_l + 1)* nRadFun(oqn_l, itype), oqn_l = 0, atoms%lmax(itype) ) ] )
    end do ! itype
    lmpMax     = maxval( lmpT(:) )
    nRadFunMax = maxval( nRadFun(:, :) )

    allocate( varphi1(atoms%jmtd, nRadFunMax, 0:atoms%lmaxd), varphi2(atoms%jmtd, nRadFunMax, 0:atoms%lmaxd) )
    allocate( r2(atoms%jmtd) )
    allocate( delrVarphi1( atoms%jmtd, nRadFunMax, 0:atoms%lmaxd), delrVarphi2( atoms%jmtd, nRadFunMax, 0:atoms%lmaxd) )
    allocate( delrGrVarphiCh1(atoms%jmtd, 2, lmpMax, -1:1), delrGrVarphiCh2(atoms%jmtd, 2, lmpMax, -1:1) )
    allocate( grVarphiCh1(atoms%jmtd, 2, lmpMax, -1:1), grVarphiCh2(atoms%jmtd, 2, lmpMax, -1:1), &
            & grVarphiChLout(2, 0:atoms%lmaxd), grVarphiChMout(-atoms%lmaxd:atoms%lmaxd, -1:1) )
    allocate( hVarphi(2, atoms%jmtd, 0:(atoms%lmaxd + 3)**2 - 1, lmpMax), vnsphEff0Varphi(2, atoms%jmtd, 0:(atoms%lmaxd + 3)**2 - 1, lmpMax), &
      & grVnsphEff0Varphi(2, atoms%jmtd, 0:(atoms%lmaxd + 3)**2 - 1, lmpMax, -1:1 ) )
    allocate( varphiVarphi(lmpMax, lmpMax, atoms%ntype), varphiHvarphi(lmpMax, lmpMax, atoms%nat))
    allocate( grVarphiVarphi(lmpMax, lmpMax, -1:1, atoms%ntype), varphiGrVarphi(lmpMax, lmpMax, -1:1, atoms%ntype), grVarphiHVarphi(lmpMax, lmpMax, -1:1, atoms%nat ) )
    allocate( grVarphiVnsphVarphiSurf(lmpMax, lmpMax) )
    allocate( vEff0MtSpH( atoms%jmtd, 0:dimens%lmd))
    allocate( vEff0NsphGrVarphi(2, atoms%jmtd, 0:(atoms%lmaxd + 3)**2, lmpMax, -1:1), r2grVeff0Varphi(2, atoms%jmtd, 0:(atoms%lmaxd + 3)**2 - 1, lmpMax, 3), &
    &r2grVeff0SphVarphi(2, atoms%jmtd, 0:(atoms%lmaxd + 3)**2 - 1, lmpMax, 3))
    allocate( hGrVarphi(2, atoms%jmtd, (atoms%lmaxd + 1)**2, lmpMax, -1:1))
    allocate( varphiVarphiDummy(lmpMax, lmpMax, atoms%ntype) )
    allocate( varphiGrVeff0SphVarphi(lmpMax, lmpMax, atoms%nat, 3) )
    allocate( varphiGrVeff0Varphi(lmpMax, lmpMax, atoms%nat, 3) )
    allocate( varphiHGrvarphi(lmpMax, lmpMax, atoms%nat, -1:1) )

    hGrVarphi(:, :, :, :, :) = cmplx(0., 0.)
    varphiGrVeff0SphVarphi(:, :, :, :) = cmplx(0., 0.)
    varphiGrVeff0Varphi(:, :, :, :) = cmplx(0., 0.)
    varphiVarphiDummy(:, :, :) = 0.
    varphiHGrvarphi(:, :, :, :) = cmplx(0., 0.)
    varphiVarphi(:, :, :)            = 0.
    varphiHvarphi(:, :, :)              = cmplx(0., 0.)
    grVarphiVarphi(:, :, :, :)       = 0.
    varphiGrVarphi(:, :, :, :)       = 0.
    grVarphiHVarphi(:, :, :, :)         = cmplx(0., 0.)
    vEff0NsphGrVarphi(:, :, :, :, :)    = cmplx(0., 0.)
    r2grVeff0SphVarphi(:, :, :, :, :)   = cmplx(0., 0.)
    r2grVeff0Varphi(:, :, :, :, :)   = cmplx(0., 0.)

    iatom = 0
    do itype = 1, atoms%ntype

      varphi1(:, :, :) = 0.
      varphi2(:, :, :) = 0.
      r2(:) = 0.

      ! Precalculate radial Jacobi determinant for later integrals
      do imesh = 1, atoms%jri(itype)
        r2(imesh) = atoms%rmsh(imesh, itype) * atoms%rmsh(imesh, itype)
      end do ! imesh

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
      ! todo Test: Compare to conventional routine by summing up the scattering channels
      call CalcChannelsGrFlpNat( atoms, itype, nRadFun, varphi1, varphi2, delrVarphi1, delrVarphi2, grVarphiChLout, grVarphiChMout, &
                                                                                                        & grVarphiCh1, grVarphiCh2 )
      ! Precalculate partial derivatives of the varphis' gradients in r-direction since it is required two times.
      delrGrVarphiCh1(:, :, :, :) = 0.
      delrGrVarphiCh2(:, :, :, :) = 0.
      do mqn_m2PrR = -1, 1
        lmp = 0
        do oqn_l = 0, atoms%lmax(itype)
          do mqn_m = -oqn_l, oqn_l
            do iradf = 1, nRadFun(oqn_l, itype)
              lmp = lmp + 1
              do ichan = 1, 2
                call Derivative( grVarphiCh1(1:atoms%jri(itype), ichan, lmp, mqn_m2PrR), itype, atoms, delrGrVarphiCh1(1:atoms%jri(itype), ichan, lmp, mqn_m2PrR) )
                call Derivative( grVarphiCh2(1:atoms%jri(itype), ichan, lmp, mqn_m2PrR), itype, atoms, delrGrVarphiCh2(1:atoms%jri(itype), ichan, lmp, mqn_m2PrR) )
              end do ! ichan
            end do !iradf
          end do ! mqn_m
        end do ! oqn_l
      end do ! mqn_m2PrR

      deallocate( delrVarphi1, delrVarphi2, delrGrVarphiCh1, delrGrVarphiCh2 )

      ! Calculate H |varphi>
      do ieqat = 1, atoms%neq(itype)
        iatom = iatom + 1
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
              vEff0MtSpH(imesh, lm) = vEff0MtSpH(imesh, lm) + vEff0MtLh(imesh, ilh, itype) * clnu_atom(imem, ilh, iatom)
            end do ! imesh
          end do ! imem
        end do ! ilh

        hVarphi = cmplx(0.0, 0.0)
        vnsphEff0Varphi = cmplx(0.0, 0.0)
        grVnsphEff0Varphi(:, :, :, :, :) = cmplx(0.0, 0.0)
        call CalcHnGrV0Varphi( atoms, lathar, itype, iatom, lmpMax, El, varphi1, varphi2, nRadFun, vEff0MtSpH, vEff0MtLh, clnu_atom, &
          & nmem_atom, mlh_atom, grVarphiCh1, grVarphiCh2, grVarphiChLout, grVarphiChMout, hVarphi, vEff0NsphGrVarphi, r2grVeff0SphVarphi, r2GrVeff0Varphi )

        deallocate( vEff0MtSpH, grVnsphEff0Varphi )
        ! todo for the integrals we have to go until lmax + 2 but later we only have to multiply the abcof until lmax, so revert it in jpInit!

                        ! we sink interests of modularity here for sake of performance so that not every integral a new loop structure is
                        ! set up. Therefore we cannot have one routine for one integral
        ! Calculate all radial integrals with no gradients.
        ! The overlap can be taken from fleur.
        ! Calculate scalar basis function matrix elements
        call CalcScalBasfMatElems( atoms, itype, iatom, nRadFun, r2, varphi1, varphi2, hVarphi, varphiVarphi, varphiHvarphi )
        ! not optimal with the first dimension
        ! This should be R for row but we want to avoid overhead setting up a new loop.
        call CalcVecBasfMatElems( atoms, itype, 2, nRadFun, r2, grVarphiChLout, grVarphiChMout, varphi1, varphi2, &
          & grVarPhiCh1, grVarphiCh2, hVarphi, grVarphiVarphi(:, :, :, itype), grVarphiHVarphi(:, :, :, iatom) )

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
   end do !itype
   iatom = 0
   r2(:) = 1.
   iatom = 0
   do itype = 1, atoms%ntype
     do ieqat = 1, atoms%neq(itype)
       iatom = iatom + 1
       do idir = 1, 3
         call CalcScalBasfMatElems( atoms, itype, iatom, nRadFun, r2, varphi1, varphi2, r2GrVeff0SphVarphi(:, :, :, :, idir), varphiVarphiDummy, varphiGrVeff0SphVarphi(:, :, :, idir) )
         call CalcScalBasfMatElems( atoms, itype, iatom, nRadFun, r2, varphi1, varphi2, r2GrVeff0Varphi(:, :, :, :, idir), varphiVarphiDummy, varphiGrVeff0Varphi(:, :, :, idir) )
       end do ! idir
     end do ! ieqat
   end do ! itype

    ! get rid of unrequired arrays
    deallocate( varphi1, varphi2, grVarphiCh1, grVarphiCh2, grVarphiChLout, grVarphiChMout )
    deallocate( r2 )
    deallocate( hVarphi, vnsphEff0Varphi )
    eps1DynMatPulTestSw = .false.
    pulayCntrb(:, :) = cmplx(0., 0.)

    ! This loop can be parallelized later
    iqpt = 1
    allocate( z1nG(dimens%nbasfcn, 3, atoms%nat, maxval(nobd(:, :))) )
    do ikpt = 1, kpts%nkpt
      z1nG = cmplx(0.0, 0.0)
      kExt(1:3) = matmul(cell%bmat, kpts%bk(1:3, ikpt))
      gbasExt(:) = 0.
      do iBas = 1, nv(1, ikpt) !+ atoms%nlotot
        gbasExt(1:3) = matmul( cell%bmat(1:3, 1:3), gBasVec(1:3, ilst(iBas, ikpt, 1)))
        do iband = 1, nobd(ikpt, 1)
          do idir = 1, 3
            z1nG(iBas, idir, 1, iband) = -0.5 * iu * ( kExt(idir) + gbasExt(idir) ) * z(iBas, iband, ikpt, 1)
          end do ! iBas
        end do ! idir
      end do ! iband

      ! Terms are activated and deactivated such that for the above Goldstein condition choice of z1nG they cancel each other.
      call Add2ndOrdWfPulayBraKets2DynMat( atoms, cell, dimens, sym, kpts, qpts, usdus, results, iqpt, ikpt, nRadFunMax, lmpMax, eps1DynMatPulTestSw, testComp2ndN1stOrdBasFuncSw, mapKpq2K, eig, nobd, nv, &
        & gBasVec, ilst, kveclo, z, z1nG, nRadFun, iloTable, varphiVarphi, varphiHvarphi, grVarphiVarphi, grVarphiHVarphi, &
        & lmpT, varphiGrVeff0SphVarphi, varphiHGrvarphi, varphiGrVarphi, varphiGrVeff0Varphi, pulayCntrb )

    end do ! ikpt

    if (.false.) then
      write(*, *) 'pulayCntr1'
      write(*, '(3(2(es16.8,1x),3x))') pulayCntrb(1, :)
      write(*, '(3(2(es16.8,1x),3x))') pulayCntrb(2, :)
      write(*, '(3(2(es16.8,1x),3x))') pulayCntrb(3, :)
    end if
    if (any(abs(pulayCntrb(:, :)) > 9e-9)) passed = .false.

    if (testCompTerm3rdBraKets) then
      dynMatPu3rdBraKetHepsSw = .false.

      testCompTerm3rdBraKetsVarKet = .true.
      testCompTerm3rdBraKetsVarBra = .false.
      pulayCntrb(:, :) = cmplx(0., 0.)
      do ikpt = 1, kpts%nkpt
        z1nG = cmplx(0.0, 0.0)
        kExt(1:3) = matmul(cell%bmat, kpts%bk(1:3, ikpt))
        gbasExt(:) = 0.
        do iBas = 1, nv(1, ikpt) !+ atoms%nlotot
          gbasExt(1:3) = matmul( cell%bmat(1:3, 1:3), gBasVec(1:3, ilst(iBas, ikpt, 1)))
          do iband = 1, nobd(ikpt, 1)
            do idir = 1, 3
              z1nG(iBas, idir, 1, iband) = -iu * ( kExt(idir) + gbasExt(idir) ) * z(iBas, iband, ikpt, 1)
            end do ! iBas
          end do ! idir
        end do ! iband

        ! Terms are activated and deactivated such that for the above Goldstein condition choice of z1nG they cancel each other.
        call Add1stOrdWfPulayBraKets2DynMat( atoms, kpts, qpts, sym, dimens, cell, usdus, stars, results, ikpt, iqpt, lmpMax, nRadFunMax, nv, GbasVec, &
        & ilst, kveclo, mapKpq2K, nobd, z, z1nG, iloTable, grVarphiVarphi, nRadFun, eig, kpq2kPrVec, &
        & grVarphiHvarphi, varphiVarphi, varphiHvarphi, vEff0IR, lmpT,&
        & eps1DynMatPulTestSw, testCompTerm3rdBraKetsVarBra, testCompTerm3rdBraKetsVarKet, dynMatPu3rdBraKetHepsSw, varphiGrVeff0SphVarphi, varphiHGrvarphi, varphiGrVarphi, pulayCntrb )

      end do ! ikpt
      if (.false.) then
        write(*, *) 'pulayCntr2'
        write(*, '(3(2(es16.8,1x),3x))') pulayCntrb(1, :)
        write(*, '(3(2(es16.8,1x),3x))') pulayCntrb(2, :)
        write(*, '(3(2(es16.8,1x),3x))') pulayCntrb(3, :)
      end if
      if (any(abs(pulayCntrb(:, :)) > 9e-9)) passed = .false.

      testCompTerm3rdBraKetsVarKet = .false.
      testCompTerm3rdBraKetsVarBra = .true.
      pulayCntrb(:, :) = cmplx(0., 0.)
      do ikpt = 1, kpts%nkpt
        z1nG = cmplx(0.0, 0.0)
        kExt(1:3) = matmul(cell%bmat, kpts%bk(1:3, ikpt))
        gbasExt(:) = 0.
        do iBas = 1, nv(1, ikpt) !+ atoms%nlotot
          gbasExt(1:3) = matmul( cell%bmat(1:3, 1:3), gBasVec(1:3, ilst(iBas, ikpt, 1)))
          do iband = 1, nobd(ikpt, 1)
            do idir = 1, 3
              z1nG(iBas, idir, 1, iband) = -iu * ( kExt(idir) + gbasExt(idir) ) * z(iBas, iband, ikpt, 1)
            end do ! iBas
          end do ! idir
        end do ! iband

        ! Terms are activated and deactivated such that for the above Goldstein condition choice of z1nG they cancel each other.
        call Add1stOrdWfPulayBraKets2DynMat( atoms, kpts, qpts, sym, dimens, cell, usdus, stars, results, ikpt, iqpt, lmpMax, nRadFunMax, nv, gBasVec, &
          & ilst, kveclo, mapKpq2K, nobd, z, z1nG, iloTable, grVarphiVarphi, nRadFun, eig, kpq2kPrVec, &
          & grVarphiHvarphi, varphiVarphi, varphiHvarphi, vEff0IR, lmpT,&
          & eps1DynMatPulTestSw, testCompTerm3rdBraKetsVarBra, testCompTerm3rdBraKetsVarKet, dynMatPu3rdBraKetHepsSw, varphiGrVeff0SphVarphi, varphiHGrvarphi, varphiGrVarphi, pulayCntrb )

      end do ! ikpt

      if (any(abs(pulayCntrb(:, :)) > 9e-9)) passed = .false.

      if (.false.) then
        write(*, *) 'pulayCntr3'
        write(*, '(3(2(es16.8,1x),3x))') pulayCntrb(1, :)
        write(*, '(3(2(es16.8,1x),3x))') pulayCntrb(2, :)
        write(*, '(3(2(es16.8,1x),3x))') pulayCntrb(3, :)
      end if
    end if

    if ( passed ) then
      write(logUnit, '(a)') '                                            |__ passed!'
    else
      write(logUnit, '(a)') '                                            |__ failed!'
      call juDFT_warn('Test Goldstone condition for the Pulay terms failed', &
        & calledby='TestComp2ndN1stOrdBasFunc', hint='Debug.')
    end if

    write(logUnit, * )
    write(logUnit, * )

  end subroutine TestComp2ndN1stOrdBasFunc

  !>-------------------------------------------------------------------------------------------------------------------------------------
  !> @author
  !> Christian-Roman Gerhorst and Markus Betzinger
  !>
  !> @brief
  !> Test Goldstein condition in surface integrals of dynamical matrix (Equation 7.119 Phd Klüppelberg).
  !>
  !> @details
  !> Neglecting the prefactor 2 in the 1st and 2nd line of 7.119 it should cancel for the IR with the terms of the third line
  !> containing a gradient of a wavefunction. In the muffin-tin the first-order wave function should reduce to the gradient of
  !> it. The latter is fullfilled provided q = 0. We check this part of the Goldstone condition within this test.
  !>
  !>
  !>
  !> @param[out] atoms       : Atoms type, see types.f90.
  !> @param[out] cell        : Unit cell type, see types.f90.
  !> @param[out] sym         : Symmetries type, see types.f90.
  !> @param[out] stars       : Stars type, see types.f90.
  !> @param[out] lathar      : Lattice harmonics type, see types.f90.
  !> @param[out] kpts        : K-points type, see types.f90.
  !> @param[out] qpts        : Q-points type, see types.f90.
  !> @param[out] Veff0       : Type containing unperturbed output potentials of Fleur, see types.f90.
  !> @param[out] usdus       : Type containing quantities consisting of the radial solutions, see types.f90.
  !> @param[out] results     : Results type, see types.f90.
  !> @param[in]  logUnit     : Unit number for juPhon.log.
  !> @param[out] kveclo      : Basis G-vectors of local orbitals.
  !> @param[out] nv          : Number of LAPW G-basis vectors for given k-point.
  !> @param[out] nRadFun     : Number of radial functions per orbital quantum number l and atom type.
  !> @param[out] iloTable    : Number of local orbital if orbital quantum number l, atom type and index p > 2 from nRadFun is given.
  !> @param[out] gBas        : G-basis vectors
  !> @param[out] mapGbas     : For various k-points G-basis vectors occur more than once, thus they are only stored once in juPhon. This
  !>                           pointer array contains the right index for gBas array to "unfold" G-basis vectors again.
  !> @param[out] mapKpq2K    : For a given k-point index and q-point index, this array gives the k-point set index of the result k + q
  !>                           (already mapped back to Brillouin zone).
  !> @param[out] eig         : Contains Kohn\--Sham eigenvalues.
  !> @param[out] rbas1       : Large components of radial solution, its energy derivative and u_LO
  !> @param[out] rbas2       : Small components of radial solution, its energy derivative and u_LO
  !> @param[out] El          : Contains LAPW and LO energy parameters.
  !> @param[out] z0           : Kohn-Sham eigenvectors.
  !> @param[out] gdp         : G-vectors of potentials and densities.
  !> @param[out] mlh_atom    : Magnetic quantum number m of lattice harmonic members for every atom.
  !> @param[out] nmem_atom   : Number of lattice harmonic members for every atom.
  !> @param[out] clnu_atom   : Phase mediating between stars and plane waves.
  !> @param[out] nobd        : Number of occupied bands per k-point and spin
  !>--------------------------------------------------------------------------------------------------------------------------------
  subroutine testGoldsteinSurf( atoms, dimens, stars, cell, results, Veff0, kpts, qpts, lathar, sym, usdus, ngdp, logUnit, nRadFun, nobd, gdp, mapKpq2K, kpq2kPrVec, nv, mapGbas, gBas, z0, eig, rbas1, rbas2, nmem_atom, mlh_atom, clnu_atom, kveclo, iloTable, El )

    use m_types
    use mod_juPhonUtils, only : Derivative, CalcChannelsGrFlpNat
    use m_jpSetupDynMatHelper, only : CalcHnGrV0Varphi
    use m_jpSetupDynMatSF, only : CalcSFintIRPsi1HepsPsi, CalcSFintIRPsiHepsPsi1, CalcSFintMTPsi1HepsPsiAndPsiHepsPsi1ExpCoeffVar, CalcSFintMTPsi1HepsPsiAndPsiHepsPsi1BasVarikpG, CalcSFintIRgradPsiHepsPsi, CalcSFintIRPsiHepsGradPsi
    use m_juDFT_NOstopNO, only : juDFT_warn

    implicit none

    type(t_atoms),                 intent(in)  :: atoms
    type(t_dimension),             intent(in)  :: dimens
    type(t_stars),                 intent(in)  :: stars
    type(t_cell),                  intent(in)  :: cell
    type(t_results),               intent(in)  :: results
    type(t_potential),             intent(in)  :: Veff0
    type(t_kpts),                  intent(in)  :: kpts
    type(t_kpts),                  intent(in)  :: qpts
    type(t_sphhar),                intent(in)  :: lathar
    type(t_sym),                   intent(in)  :: sym
    type(t_usdus),                 intent(in)  :: usdus

    ! Scalar parameters
    integer,                       intent(in)  :: ngdp
    integer,                       intent(in)  :: logUnit

    ! Array parameters
    integer,                       intent(in)  :: nRadFun(0:, :)
    integer,                       intent(in)  :: nobd(:, :)
    integer,                       intent(in)  :: gdp(:, :)
    integer,                       intent(in)  :: mapKpq2K(:, :)
    integer,                       intent(in)  :: kpq2kPrVec(:, :, :)
    integer,                       intent(in)  :: nv(:, :)
    integer,                       intent(in)  :: mapGbas(:, :, :)
    integer,                       intent(in)  :: gBas(:, :)
    MCOMPLEX,                      intent(in)  :: z0(:, :, :, :)
    real,                          intent(in)  :: eig(:, :, :)
    real,                          intent(in)  :: rbas1(:, :, 0:, :)
    real,                          intent(in)  :: rbas2(:, :, 0:, :)
    integer,                       intent(in)  :: nmem_atom(0:, :)
    integer,                       intent(in)  :: mlh_atom(:,0:,:)
    complex,                       intent(in)  :: clnu_atom(:,0:,:)
    integer,                       intent(in)  :: kveclo(:,:)
    integer,                       intent(in)  :: iloTable(:, 0:, :)
    real,                          intent(in)  :: El(:, 0:, :, :)

    ! Scalar variables
    integer                                    :: iDatomA
    integer                                    :: iDtypeA
    integer                                    :: iDatomB
    integer                                    :: iDtypeB
    integer                                    :: iDeqatA
    integer                                    :: iDeqatB
    integer                                    :: oqn_l
    integer                                    :: iradf
    integer                                    :: imesh
    integer                                    :: lmpMax
    integer                                    :: ptsym
    integer                                    :: ilh
    integer                                    :: lm_pre
    integer                                    :: imem
    integer                                    :: lm
    integer                                    :: mqn_m
    integer                                    :: nRadFunMax
    integer                                    :: idirR
    integer                                    :: idirC
    logical                                    :: testGoldstein
    integer                                    :: iqpt
    logical                                    :: passed = .true.

    ! Array variables
    integer,           allocatable             :: lmpT(:)
    real,              allocatable             :: varphi1(:, :, :)
    real,              allocatable             :: varphi2(:, :, :)
    real,              allocatable             :: delrVarphi1(:, :, :)
    real,              allocatable             :: delrVarphi2(:, :, :)
    integer,           allocatable             :: grVarphiChLout(:, :)
    integer,           allocatable             :: grVarphiChMout(:, :)
    real,              allocatable             :: grVarphiCh1(:, :, :, :)
    real,              allocatable             :: grVarphiCh2(:, :, :, :)
    complex,           allocatable             :: vEff0MtSpH(:, :)
    complex,           allocatable             :: hFullVarphi(:, :, :, :)
    complex,           allocatable             :: r2grVeff0SphVarphi(:, :, :, :, :)
    complex,           allocatable             :: r2grVeff0Varphi(:, :, :, :, :)
    complex,           allocatable             :: vEff0NsphGrVarphi(:, :, :, :, :)
    complex,           allocatable             :: surfInts(:, :)
    complex                                    :: surfInt(3, 3)

    write (logUnit, '(a)')   'Test Goldstone condition for the surface integrals.'
    write ( logUnit, '(a)' ) '---------------------------------------------------'

    testGoldstein = .true.
    iqpt = 1

    ! Required for initialization
    allocate( lmpT(atoms%ntype) )
    lmpT(:) = 0
    do iDtypeB = 1, atoms%ntype
      lmpT(iDtypeB) = sum( [ ( (2 * oqn_l + 1)* nRadFun(oqn_l, iDtypeB), oqn_l = 0, atoms%lmax(iDtypeB) ) ] )
    end do ! itype
    lmpMax     = maxval( lmpT(:) )
    nRadFunMax = maxval( nRadFun(:, :) )

    ! Allocate and initialize arrays within the code
    allocate( varphi1(atoms%jmtd, nRadFunMax, 0:atoms%lmaxd), varphi2(atoms%jmtd, nRadFunMax, 0:atoms%lmaxd) )
    allocate( delrVarphi1( atoms%jmtd, nRadFunMax, 0:atoms%lmaxd), delrVarphi2( atoms%jmtd, nRadFunMax, 0:atoms%lmaxd) )
    allocate( grVarphiCh1(atoms%jmtd, 2, lmpMax, -1:1), grVarphiCh2(atoms%jmtd, 2, lmpMax, -1:1), &
            & grVarphiChLout(2, 0:atoms%lmaxd), grVarphiChMout(-atoms%lmaxd:atoms%lmaxd, -1:1) )
    allocate( vEff0MtSpH( atoms%jmtd, 0:dimens%lmd), vEff0NsphGrVarphi(2, atoms%jmtd, 0:(atoms%lmaxd + 3)**2 - 1, lmpMax, -1:1), &
            & r2grVeff0SphVarphi(2, atoms%jmtd, 0:(atoms%lmaxd + 3)**2 - 1, lmpMax, 3), &
            & r2grVeff0Varphi(2, atoms%jmtd, 0:(atoms%lmaxd + 3)**2 - 1, lmpMax, 3) )
    allocate( hFullVarphi(2, atoms%jmtd, 0:(atoms%lmaxd + 3)**2 - 1, lmpMax)  )
    allocate(surfInts(3 * atoms%nat, 3 * atoms%nat))


    surfInts(:, :) = cmplx(0., 0.)
    iDatomA = 0
    do iDtypeA = 1, atoms%ntype
      do iDeqatA = 1, atoms%neq(iDtypeA)
        iDatomA = iDatomA + 1

        iDatomB = 0
        do iDtypeB = 1, atoms%ntype
          do iDeqatB = 1, atoms%neq(iDtypeB)
            iDatomB = iDatomB + 1

            ! Divide out factor r from FLEUR generated basis function
            varphi1(:, :, :) = 0.
            varphi2(:, :, :) = 0.
            delrVarphi1(:, :, :) = cmplx(0., 0.)
            delrVarphi2(:, :, :) = cmplx(0., 0.)
            do oqn_l = 0, atoms%lmax(iDtypeB)
              do iradf = 1, nRadFun(oqn_l, iDtypeB)
                do imesh = 1, atoms%jri(iDtypeB)
                  ! In Fleur the radial solutions u_lp are multiplied by a factor r by default to avoid an additional multiplication of the
                  ! Jacobi determinant r^2 in radial integrals given spherical coordinates.
                  varphi1(imesh, iradf, oqn_l) = rbas1(imesh, iradf, oqn_l, iDtypeB) / atoms%rmsh(imesh, iDtypeB)
                  varphi2(imesh, iradf, oqn_l) = rbas2(imesh, iradf, oqn_l, iDtypeB) / atoms%rmsh(imesh, iDtypeB)
                end do ! imesh
                call Derivative( varphi1(1:atoms%jri(iDtypeB), iradf, oqn_l), iDtypeB, atoms, delrVarphi1(1:atoms%jri(iDtypeB), iradf, oqn_l) )
                call Derivative( varphi2(1:atoms%jri(iDtypeB), iradf, oqn_l), iDtypeB, atoms, delrVarphi2(1:atoms%jri(iDtypeB), iradf, oqn_l) )
              end do ! iradf
            end do ! oqn_l

            ! Calculate the application of the gradient  onto the MT basis functions (matching coefficients
            ! have no spatial dependence) and determing its scattering channels.
            grVarphiChLout(:, :) = 0
            grVarphiChMout(:, :) = 0
            grVarphiCh1(:, :, :, :) = 0.
            grVarphiCh2(:, :, :, :) = 0.
            call CalcChannelsGrFlpNat( atoms, iDtypeB, nRadFun, varphi1, varphi2, delrVarphi1, delrVarphi2, grVarphiChLout, grVarphiChMout, &
                                                                                                            & grVarphiCh1, grVarphiCh2 )

            vEff0MtSpH(:, :) = cmplx(0.0, 0.0)
            ptsym = atoms%ntypsy(iDatomB)
            do ilh = 0, lathar%nlh(ptsym)
              oqn_l = lathar%llh(ilh, ptsym)
              lm_pre = oqn_l * (oqn_l + 1)
              do imem = 1, nmem_atom(ilh, iDatomB)
                mqn_m = mlh_atom(imem, ilh, iDatomB)
                lm = lm_pre + mqn_m
                !todo one could only evaluate the vEff0MtSpH which have a contribution, i.e. the oqn_l and mqn_m which are in llh and mlh_atom
                ! maybe construct a pointer and run only over them to make it memory efficient.
                do imesh = 1, atoms%jri(iDtypeB)
                  vEff0MtSpH(imesh, lm) = vEff0MtSpH(imesh, lm) + Veff0%vr(imesh, ilh, iDtypeB, 1) * clnu_atom(imem, ilh, iDatomB)
                end do ! imesh
              end do ! imem
            end do ! ilh

            hFullVarphi = cmplx(0.0, 0.0)
            call CalcHnGrV0Varphi( atoms, lathar, iDtypeB, iDatomB, lmpMax, El, varphi1, varphi2, nRadFun, vEff0MtSpH, Veff0%vr(:, :, :, 1), clnu_atom, &
              & nmem_atom, mlh_atom, grVarphiCh1, grVarphiCh2, grVarphiChLout, grVarphiChMout, hFullVarphi, vEff0NsphGrVarphi, r2grVeff0SphVarphi, r2grVeff0Varphi )
            vEff0NsphGrVarphi(:, :, :, :, :) = cmplx(0., 0.)
            r2grVeff0SphVarphi(:, :, :, :, :) = cmplx(0., 0.)
            r2grVeff0Varphi(:, :, :, :, :) = cmplx(0., 0.)

            ! Calculates IR surface integral with integrand Psi1 (H - eps) Psi0
            surfInt(:, :) = cmplx(0., 0.)
            call CalcSFintIRPsi1HepsPsi( atoms, dimens, stars, cell, results, Veff0, kpts, qpts, ngdp, iqpt, iDtypeB, iDatomB, iDatomA, nobd,&
                                                  & gdp, mapKpq2K, kpq2kPrVec, nv, mapGbas, gBas, z0, eig, surfInt, .true. )
            do idirC = 1, 3
              do idirR = 1, 3
                surfInts(idirR + (iDatomB - 1) * 3, idirC + (iDatomA - 1)) = surfInts(idirR + (iDatomB - 1) * 3, idirC + (iDatomA - 1) * 3) + surfInt(idirR, idirC)
              end do
            end do

            ! Calculates IR surface integral with integrand Psi0 ( H - eps) Psi1
            surfInt(:, :) = cmplx(0., 0.)
            call CalcSFintIRPsiHepsPsi1( atoms, dimens, stars, cell, results, Veff0, kpts, qpts, ngdp, iqpt, iDtypeB, iDatomB, iDatomA,&
                                                  & nobd, gdp, mapKpq2K, kpq2kPrVec, nv, mapGbas, gBas, z0, eig, surfInt, .true. )
            do idirC = 1, 3
              do idirR = 1, 3
                surfInts(idirR + (iDatomB - 1) * 3, idirC + (iDatomA - 1)) = surfInts(idirR + (iDatomB - 1) * 3, idirC + (iDatomA - 1) * 3) + surfInt(idirR, idirC)
              end do
            end do

            ! Calculates the MT surface integrals sum with the integrands Psi1* (H - eps) Psi0 and Psi0* (H - eps) Psi1 for q = 0
            surfInt(:, :) = cmplx(0., 0.)
            call CalcSFintMTPsi1HepsPsiAndPsiHepsPsi1ExpCoeffVar( atoms, sym, dimens, usdus, kpts, cell, results, iqpt, iDtypeB, &
              & iDatomB, iDatomA, lmpMax, nRadFun, eig, varphi1, varphi2, hFullVarphi, mapKpq2K,  gBas, mapGbas, nv,&
              & kveclo, z0, lmpT, nobd, iloTable, surfInt, .true. )
            do idirC = 1, 3
              do idirR = 1, 3
                surfInts(idirR + (iDatomB - 1) * 3, idirC + (iDatomA - 1)) = surfInts(idirR + (iDatomB - 1) * 3, idirC + (iDatomA - 1) * 3) + surfInt(idirR, idirC)
              end do
            end do

          if (iDatomA == iDatomB) then

            ! Calculates the MT surface integrals sum with the integrands (i(k + G)Psi0))* (H - eps) Psi0 and Psi0* (H - eps) (i(k + G)Psi0))
            surfInt(:, :) = cmplx(0., 0.)
            call CalcSFintMTPsi1HepsPsiAndPsiHepsPsi1BasVarikpG( atoms, sym, dimens, usdus, kpts, cell, results, lmpMax, iDtypeA,      &
              & iDatomA, nRadFun, eig, hFullVarphi, gBas, mapGbas, nv, kveclo, z0, nobd, lmpT, iloTable, varphi1, varphi2, surfInt )
            do idirC = 1, 3
              do idirR = 1, 3
                surfInts(idirR + (iDatomB - 1) * 3, idirC + (iDatomA - 1)) = surfInts(idirR + (iDatomB - 1) * 3, idirC + (iDatomA - 1) * 3) + surfInt(idirR, idirC)
              end do
            end do

            ! Calcualtes the IR surface integral with the integrand grad Psi (H - eps) Psi
            surfInt(:, :) = cmplx(0., 0.)
            call CalcSFintIRgradPsiHepsPsi( atoms, stars, cell, dimens, kpts, results, Veff0, ngdp, iDtypeB, iDatomB, nobd, eig, gBas, &
                                                                                                 & mapGbas, nv, gdp, z0, surfInt )
            do idirC = 1, 3
              do idirR = 1, 3
                surfInts(idirR + (iDatomB - 1) * 3, idirC + (iDatomA - 1)) = surfInts(idirR + (iDatomB - 1) * 3, idirC + (iDatomA - 1) * 3) + surfInt(idirR, idirC)
              end do
            end do

            ! Calcualtes the IR surface integral with the integrand Psi (H - eps) grad Psi
            surfInt(:, :) = cmplx(0., 0.)
            call CalcSFintIRPsiHepsGradPsi( atoms, stars, cell, dimens, kpts, results, Veff0, ngdp, iDtypeB, iDatomB, nobd, eig, gBas, &
            &  mapGbas, nv, gdp, z0, surfInt )
            do idirC = 1, 3
              do idirR = 1, 3
                surfInts(idirR + (iDatomB - 1) * 3, idirC + (iDatomA - 1)) = surfInts(idirR + (iDatomB - 1) * 3, idirC + (iDatomA - 1) * 3) + surfInt(idirR, idirC)
              end do
            end do
          end if
          end do ! iDeqatA
        end do ! iDtypeA
      end do ! iDeqatB
    end do ! iDtypeB

    if (any(abs(surfInts(:, :)) > 9e-9)) passed = .false.
    if ( passed ) then
      write(logUnit, '(a)') '                                                  |__ passed!'
    else
      write(logUnit, '(a)') '                                                  |__ failed!'
      call juDFT_warn('Test Goldstone condition for the surface integrals. failed', &
        & calledby='testGoldsteinSurf', hint='Debug.')
    end if
  end subroutine testGoldsteinSurf

end module m_jpTestGoldsteinCond
