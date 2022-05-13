module m_jpTestDynMatSurf

#include "cppmacro.h"
  implicit none

  contains

  subroutine TestDynMatSurf( atoms, dimens, stars, sym, cell, kpts, qpts, input, lathar, usdus, Veff0, results, testComp3ArgSFIntsSw, testComp2ArgSFIntsSw, testComp2ArgGrSFIntsSw, nRadFun, El, eig,  &
                          & gdp, rbas1, rbas2, nmem_atom, mlh_atom, clnu_atom, z0, nv, gbas, mapGbas, nobd, kveclo, iloTable, ngdp,&
                          & kpq2kPrVec, mapKpq2K, logUnit, memd_atom, rho0IR, rho0MT )

    use m_types

    implicit none

    type(t_atoms),                  intent(in)  :: atoms
    type(t_dimension),              intent(in)  :: dimens
    type(t_stars),                  intent(in)  :: stars
    type(t_sym),                    intent(in)  :: sym
    type(t_cell),                   intent(in)  :: cell
    type(t_kpts),                   intent(in)  :: kpts
    type(t_kpts),                   intent(in)  :: qpts
    type(t_input),                  intent(in)  :: input
    type(t_sphhar),                 intent(in)  :: lathar
    type(t_usdus),                  intent(in)  :: usdus
    type(t_potential),              intent(in)  :: Veff0
    type(t_results),                intent(in)  :: results

    ! Scalar parameter
    integer,                        intent(in)  :: ngdp
    integer,                        intent(in)  :: logUnit
    integer,                        intent(in)  :: memd_atom
    logical,                        intent(in)  :: testComp3ArgSFIntsSw
    logical,                        intent(in)  :: testComp2ArgSFIntsSw
    logical,                        intent(in)  :: testComp2ArgGrSFIntsSw

    ! Array parameter
    integer,                        intent(in)  :: nRadFun(0:, :)
    real,                           intent(in)  :: El(:, 0:, :, :)
    real,                           intent(in)  :: rbas1(:, :, 0:, :, :)
    real,                           intent(in)  :: rbas2(:, :, 0:, :, :)
    integer,                        intent(in)  :: nmem_atom(0:, :)
    integer,                        intent(in)  :: mlh_atom(:,0:,:)
    complex,                        intent(in)  :: clnu_atom(:,0:,:)
    MCOMPLEX,                       intent(in)  :: z0(:,:,:,:)
    integer,                        intent(in)  :: nv(:, :)
    integer,                        intent(in)  :: gbas(:, :)
    integer,                        intent(in)  :: mapGbas(:, :, :)
    integer,                        intent(in)  :: nobd(:, :)
    integer,                        intent(in)  :: kveclo(:,:)
    integer,                        intent(in)  :: iloTable(:, 0:, :)
    integer,                        intent(in)  :: gdp(:, :)
    real,                           intent(in)  :: eig(:, :, :)
    integer,                        intent(in)  :: kpq2kPrVec(:, :, :)
    integer,                        intent(in)  :: mapKpq2K(:, :)
    complex,                        intent(in)  :: rho0IR(:,:)
    real,                           intent(in)  :: rho0MT(:, 0:, :, :)


    if ( testComp3ArgSFIntsSw .or. testComp2ArgSFIntsSw .or. testComp2ArgGrSFIntsSw ) then
      write(*, *)
      write(*, '(a)') 'Initiating dynamical matrix surface integral test(s)...'
      write(*, '(a)') '-------------------------------------------------------'
      write(logUnit, *)
      write(logUnit, '(a)') 'Dynamical matrix surface integral test(s)'
      write(logUnit, '(a)') '*****************************************'
      write(logUnit, *)
    else
      write(*, '(a)') '---------------------------------------------------'
      write(*, '(a)') 'DISABLED dynamical matrix surface integral test(s)!'
      write(*, '(a)') '---------------------------------------------------'
      write(logUnit, '(a)')
      write(logUnit, '(a)') 'DISABLED dynamical matrix surface integral test(s)!'
      write(logUnit, '(a)') '***************************************************'
      return
    end if

    if (testComp3ArgSFIntsSw) then
      write(*, '(2x,a)') 'Performing TestComp3ArgSFInts...'
      ! Compares surface integrals of the form Psi* (H - eps) Psi in the IR and MT with different routines.
      call TestComp3ArgSFInts( atoms, dimens, stars, sym, cell, kpts, qpts, input, lathar, usdus, Veff0, results, nRadFun, eig, kpq2kPrVec, El,   &
                          & gdp, rbas1(:, :, :, :, 1), rbas2(:, :, :, :, 1), nmem_atom, mlh_atom, clnu_atom, z0, nv, gbas, mapGbas, nobd, kveclo, iloTable, ngdp, mapKpq2K, logUnit )
      else
      write(*, '(2x,a)') 'DISABLED TestComp3ArgSFInts!'
      write( logUnit, '(a)' ) 'Test of general surface integral for dynamical matrix'
      write( logUnit, '(a)' ) '-----------------------------------------------------'
      write( logUnit, '(a)' ) '                                                    |__ DISABLED!'
      write (logUnit, *)

    end if

    if (testComp2ArgSFIntsSw) then
      write(*, '(2x,a)') 'Performing TestComp2ArgSFInts...'
      ! Compares the tested 3-argument surface integral sum_kn dS Psi*_kn grVeff Psi_kn to surface integral dS rho grVeff
      call TestComp2ArgSFInts( atoms, sym, cell, lathar, stars, dimens, kpts, results, usdus, ngdp, logUnit, memd_atom,     &
        & clnu_atom, nmem_atom, rbas1(:, :, :, :, 1), rbas2(:, :, :, :, 1), mlh_atom, rho0IR, rho0MT, gdp, nobd, gbas, mapgbas, nv,&
        & z0, nRadFun, kveclo, iloTable , kpq2kPrVec)
    else
      write(*, '(2x,a)') 'DISABLED TestComp2ArgSFInts!'
      write( logUnit, '(a)' ) 'Test of dynamical matrix 2-argument surface integral'
      write( logUnit, '(a)' ) '----------------------------------------------------'
      write( logUnit, '(a)' ) '                                                    |__ DISABLED!'
      write (logUnit, *)
    end if

    if (testComp2ArgGrSFIntsSw) then
      write(*, '(2x,a)') 'Performing TestComp2ArgGrSFInts...'
      call TestComp2ArgGrSFInts( atoms, dimens, lathar, cell, kpts, qpts, stars, Veff0, results, sym, usdus, ngdp, logUnit, gdp, nv, gbas, mapGbas, nobd, z0, nmem_atom, mlh_atom, clnu_atom, nRadFun, rbas1(:, :, :, :, 1), rbas2(:, :, :, :, 1), kveclo, iloTable, mapKpq2K, kpq2kPrVec, eig )
    else
      write(*, '(2x,a)') 'DISABLED TestComp2ArgGrSFInts!'
      write( logUnit, * )
      write( logUnit, '(a)' ) 'Test of surface integrals containing wave function gradients'
      write( logUnit, '(a)' ) '------------------------------------------------------------'
      write( logUnit, '(a)' ) '                                                           |__ DISABLED!'
      write (logUnit, *)
    end if

  end subroutine TestDynMatSurf


  !>-------------------------------------------------------------------------------------------------------------------------------------
  !> @author
  !> Christian-Roman Gerhorst
  !>
  !> @brief
  !> Compare surfInt dS rho grVeff and surfInt sum_kn surfInt dS Psi*_kn grVeff Psi_kn in the IR and the MT
  !>
  !> @details
  !> The k-point and band sum over the wavefunctions is equal to the density. Therefore, we can compare the tested 3-argument
  !> integral with the 2-argument integral. The second argument is the gradient of the potential.
  !>
  !> @note
  !> The surface integrals in the IR and the MT are the same in the order of magnitude of the integrand factor continuity, i.e. the
  !> density and the gradient of the potential are continious up to 1e-6 approx., therefore also the surface integrals should be
  !> mirroring this continuity by being equal up to 1e-6.
  !>
  !> @param[out] atoms       : Atoms type, see types.f90.
  !> @param[out] cell        : Unit cell type, see types.f90.
  !> @param[out] sym         : Symmetries type, see types.f90.
  !> @param[out] stars       : Stars type, see types.f90.
  !> @param[out] lathar      : Lattice harmonics type, see types.f90.
  !> @param[out] dimens      : Dimension type, see types.f90.
  !> @param[out] kpts        : K-points type, see types.f90.
  !> @param[out] qpts        : Q-points type, see types.f90.
  !> @param[out] Veff0       : Type containing unperturbed output potentials of Fleur, see types.f90.
  !> @param[out] usdus       : Type containing quantities consisting of the radial solutions, see types.f90.
  !> @param[out] results     : Results type, see types.f90.
  !> @param[in]  logUnit     : Unit number for juPhon.log.
  !> @param[out] ngdp        : Number of G-vectors for potentials and densities.
  !> @param[out] kveclo      : Basis G-vectors of local orbitals.
  !> @param[out] nv          : Number of LAPW G-basis vectors for given k-point.
  !> @param[out] nRadFun     : Number of radial functions per orbital quantum number l and atom type.
  !> @param[out] iloTable    : Number of local orbital if orbital quantum number l, atom type and index p > 2 from nRadFun is given.
  !> @param[out] rbas1       : Large components of radial solution, its energy derivative and u_LO
  !> @param[out] rbas2       : Small components of radial solution, its energy derivative and u_LO
  !> @param[out] El          : Contains LAPW and LO energy parameters.
  !> @param[out] z0          : Kohn-Sham eigenvectors.
  !> @param[out] gdp         : G-vectors of potentials and densities.
  !> @param[out] mlh_atom    : Magnetic quantum number m of lattice harmonic members for every atom.
  !> @param[out] nmem_atom   : Number of lattice harmonic members for every atom.
  !> @param[out] clnu_atom   : Phase mediating between stars and plane waves.
  !> @param[out] nobd        : Number of occupied bands per k-point and spin
  !>
  !> @ todo complete documentation
  !>--------------------------------------------------------------------------------------------------------------------------------------
  subroutine TestComp2ArgSFInts( atoms, sym, cell, lathar, stars, dimens, kpts, results, usdus, ngdp, logUnit, memd_atom, clnu_atom, nmem_atom,  &
                            & rbas1, rbas2, mlh_atom, rho0IR, rho0MT, gdp, nobd, gbas, mapGbas,&
                            & nv, z0, nRadFun, kveclo, iloTable , kpq2kPrVec)
    use m_types, only : t_atoms, t_sym, t_cell, t_sphhar, t_stars, t_dimension, t_potential, t_kpts, t_results, t_usdus
     
    use m_jpPotDensHelper, only : CalcIRdVxcKern, CalcMTdVxcKern
    use m_jpGrVeff0, only : GenGrVeff0
    use m_jpSetupDynMatSF, only : IRcoeffVeffUv, calcSurfVeff, PrepareMTSurfIntDM, CalcSIntMT
    use mod_juPhonUtils, only : calcGrR2FinLH
    use m_jpPotDensHelper, only : convertStar2G
    use m_jpConstants, only : iu, fpi
    use m_jpSetupDynMatHelper, only : CalcSurfIntIRDynMat, CalcSurfIntMTDynMat, CalcFnsphVarphi
    use m_abcof3
    use m_juDFT_NOstopNO, only : juDFT_warn

    implicit none

    ! Type parameter
    type(t_atoms),                  intent(in) :: atoms
    type(t_sym),                    intent(in) :: sym
    type(t_cell),                   intent(in) :: cell
    type(t_sphhar),                 intent(in) :: lathar
    type(t_stars),                  intent(in) :: stars
    type(t_dimension),              intent(in) :: dimens
    type(t_kpts),                   intent(in) :: kpts
    type(t_results),                intent(in) :: results
    type(t_usdus),                  intent(in) :: usdus

    ! Scalar parameter
    integer,                        intent(in) :: ngdp
    integer,                        intent(in) :: logUnit
    integer,                        intent(in) :: memd_atom

    ! Array parameters
    complex,                        intent(in) :: clnu_atom(:, 0:, :)
    integer,                        intent(in)  :: nRadFun(0:, :)
    integer,                        intent(in) :: nmem_atom(0:, :)
    integer,                        intent(in) :: mlh_atom(:, 0:, :)
    real,                           intent(in) :: rho0MT(:, 0:, :, :)
    integer,                        intent(in) :: gdp(:, :)
    complex,                        intent(in) :: rho0IR(:,:)
    integer,                        intent(in)  :: nobd(:, :)
    integer,                        intent(in)  :: gbas(:, :)
    integer,                        intent(in)  :: mapGbas(:, :, :)
    integer,                        intent(in)  :: nv(:, :)
    real,                           intent(in)  :: rbas1(:, :, 0:, :)
    real,                           intent(in)  :: rbas2(:, :, 0:, :)
    MCOMPLEX,                       intent(in)  :: z0(:,:,:,:)
    integer,                        intent(in)  :: kveclo(:,:)
    integer,                        intent(in)  :: iloTable(:, 0:, :)
    integer,                        intent(in)  :: kpq2kPrVec(:, :, :)

    ! Type variables
    type(od_inp)                                :: odi
    type(od_sym)                                :: ods

    ! Scalar variables
    integer                                    :: idir
    integer                                    :: iG
    integer                                    :: iatom
    integer                                    :: ieqat
    integer                                    :: iDtype
    integer                                    :: itype
    integer                                    :: ilh
    integer                                    :: idirR
    integer                                    :: idirC
    integer                                    :: iDatom
    integer                                    :: iDeqat
    integer                                    :: oqn_l
    integer                                    :: lm_pre
    integer                                    :: mqn_m
    integer                                    :: lm
    integer                                    :: imesh
    logical                                    :: testGoldstein
    integer                                    :: coScale
    integer                                    :: ikpt
    integer                                    :: iband
    integer                                    :: lmpMax
    integer                                    :: nRadFunMax
    integer                                    :: iradf
    integer                                    :: chanMaxBra
    integer                                    :: chanMaxKet
    integer                                    :: lmaxBra
    integer                                    :: ptsym
    integer                                    :: mqn_m2PrBra
    integer                                    :: mqn_m2PrKet
    integer                                    :: lmp
    integer                                    :: nmat
    integer                                    :: pMaxLocal
    integer                                    :: imem
    logical                                    :: harSw = .true.
    logical                                    :: extSw = .true.
    logical                                    :: xcSw = .true.
    logical                                    :: fullGrVeff0 = .true.
    logical                                    :: passed = .true.

    ! Array variables
    complex,           allocatable             :: rho0IRpw(:, :)
    complex,           allocatable             :: grVxcIRKern(:)
    real,              allocatable             :: gaussWghts(:) ! gaussian weights belonging to gausPts
    complex,           allocatable             :: ylm(:, :)
    real,              allocatable             :: dKernMTGPts(:, :, :)
    real,              allocatable             :: r2Rho0MTVal(:, :, :, :)
    real,              allocatable             :: rho0MTVal(:, :, :, :)
    complex,           allocatable             :: grRho0MT(:, :, :, :)
    complex,           allocatable             :: r2GrRho0MT(:, :, :, :)
    complex,           allocatable             :: grVeff0MT(:, :, :, :)
    complex,           allocatable             :: grVeff0IR(:, :)
    complex,           allocatable             :: veffUvIR(:, :)
    real,              allocatable             :: varphi1(:, :, :)
    real,              allocatable             :: varphi2(:, :, :)
    complex,           allocatable             :: surfIntVFast(:, :)
    complex,           allocatable             :: grVeff0Varphi(:, :, :, :)
    integer,           allocatable             :: lmpT(:)
    complex,           allocatable             :: hFullNoAbcofBK(:, :, :, :)
    complex,           allocatable             :: overlapNoAbcofBK(:, :, :, :)
    integer,           allocatable             :: muKet(:, :)
    integer,           allocatable             :: muBra(:, :)
    real,              allocatable             :: varphiBra1(:, :, :, :)
    real,              allocatable             :: varphiKet1(:, :, :, :)
    real,              allocatable             :: varphiBra2(:, :, :, :)
    real,              allocatable             :: varphiKet2(:, :, :, :)
    integer,           allocatable             :: lambdaKet(:, :)
    integer,           allocatable             :: lambdaBra(:, :)
    complex,           allocatable             :: a(:, :, :)
    complex,           allocatable             :: b(:, :, :)
    complex,           allocatable             :: bascof_lo(:, :, :, :, :)
    integer,           allocatable             :: ngoprI(:)
    complex,           allocatable             :: ab0cofKet(:, :)
    complex,           allocatable             :: overlap(:)
    complex,           allocatable             :: surfIntMT(:)
    real,              allocatable             :: r2Rho0MT(:, :, :)
    complex,           allocatable             :: rho0IRnCt(:, :)
    complex,           allocatable              :: rho0MTsh(:, :, :, :)
    complex,           allocatable              :: grRho0IR(:, :)
    complex                                    :: surfIntIR2Arg(3, 3)
    complex                                    :: surfIntMT2Arg(3, 3)
    complex                                    :: surfIntSum(3, 3)
    complex                                    :: surfIntSumMT(3, 3)
    real                                       :: Gext(3)

    write( logUnit, '(a)' ) 'Test of dynamical matrix 2-argument surface integral '
    write( logUnit, '(a)' ) '-----------------------------------------------------'

    surfIntIR2Arg(:, :) = cmplx(0., 0.)
    surfIntMT2Arg(:, :) = cmplx(0., 0.)

    ! NOTE: It is relevant whether we take the Weinert routine or the numerical gradients of the section-wise defined effective
    ! potential's gradient, because we are dealing with surface integrals and do not ensure continuity when not using the Weinert
    ! method.
    
    !NOTE: Probably the upper finding was due to a bug and placing the conjg wrongly in the MT surface integral. We can test this out
    ! when we use the surface integral in the deprecated module for the dynamical matrix that uses the numerical gradients of the
    ! quantities to calculate. For reasons of consistency we still use the Weinert generation of the potential.

     ! Use the Weinert method to calculate the gradient of the unperturbed effective potential
     ! ---------------------------------------------------------------------------------------

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
     ! Calculate the gradient of the unperturbed density.
     !   Any factors as in Fleur have been removed. In order to improve the numerical accuracy of the gradient r^2 is multiplied
     !   before and removed after the calculation of the gradient. Thereyby, subtraction of small numbers close to the core are
     !   avoided.
     allocate( r2Rho0MT( atoms%jmtd, 0:lathar%nlhd, atoms%ntype) )
     allocate( grRho0MT( atoms%jmtd, ( atoms%lmaxd + 1 )**2, atoms%nat, 3) )
     r2Rho0MT(:, :, :) = 0.
     grRho0MT(:, :, :, :) = cmplx(0., 0.)

     do iDtype = 1, atoms%ntype
       do ilh = 0, lathar%nlhd
         do imesh = 1, atoms%jri(iDtype)
           ! SPIN MISSING
           r2Rho0MT(imesh, ilh, iDtype) = rho0MT(imesh, ilh, iDtype, 1) * atoms%rmsh(imesh, iDtype) * atoms%rmsh(imesh, iDtype)
         end do ! imesh
       end do ! ilh
     end do ! iDtype

     call CalcGrR2FinLH( atoms, lathar, clnu_atom, nmem_atom, mlh_atom, r2Rho0MT, r2GrRho0MT )

     do idirR = 1, 3
       iDatom = 0
       do iDtype = 1, atoms%ntype
         do iDeqat = 1, atoms%neq(iDtype)
           iDatom = iDatom + 1
           do oqn_l = 0, atoms%lmax(iDtype)
             lm_pre = oqn_l * (oqn_l + 1) + 1
             do mqn_m = -oqn_l, oqn_l
               lm = lm_pre + mqn_m
               do imesh = 1, atoms%jri(iDtype)
                 grRho0MT(imesh, lm, iDatom, idirR) = r2GrRho0MT(imesh, lm, iDatom, idirR) / atoms%rmsh(imesh, iDtype)**2
               end do ! imesh
             end do ! mqn_m
           end do ! oqn_l
         end do ! iDeqat
       end do ! iDtype
     end do ! idirR

     ! Precalculate quantities for the kernel of the unperturbed effective potential's gradient.
     ! SPIN MISSING
     call CalcIRdVxcKern(stars, gdp, ngdp, rho0IR(:, 1), grVxcIRKern)
     ! SPIN MISSING
     call CalcMTdVxcKern(atoms, dimens, lathar, rho0MT(:, :, :, 1), nmem_atom, clnu_atom, mlh_atom, gaussWghts, ylm, dKernMTGPts)

    allocate(rho0IRpw(ngdp, 1))
    rho0IRpw(:, :) = cmplx(0., 0.)
    call convertStar2G( rho0IR(:, 1), rho0IRpw(:, 1), stars, ngdp, gdp )
    allocate(grRho0IR(ngdp, 3))
    grRho0IR(:, :) = cmplx(0., 0.)
    do idir = 1, 3
      do iG = 1, ngdp
        Gext(1:3) = matmul(cell%bmat(1:3, 1:3), gdp(1:3, iG))
        grRho0IR(iG, idir)  = iu * Gext(idir) * rho0IRpw(iG, 1)
      end do ! iG
    end do ! idir

     ! Generates gradient of the unperturbed effective potential by using the Weinert method.
     harSw = .true.
     extSw = .true.
     xcSw = .true.
     fullGrVeff0 = .true.
     testGoldstein = .false.
     ! SPIN MISSING
     call GenGrVeff0(atoms, cell, dimens, stars, ngdp, harSw, extSw, xcSw, gdp, rho0IRpw, rho0MTsh, grRho0IR, grRho0MT, gaussWghts, ylm, &
       & dKernMTGPts, grVxcIRKern, testGoldstein, fullGrVeff0, grVeff0IR, grVeff0MT )

    rho0IRpw(:, :) = cmplx(0., 0.)
    allocate(rho0IRncT(stars%n3d, 1))
    rho0IRncT(:, :) = cmplx(0.,0.)
    rewind(7800)
    read(7800) rho0IRncT
    call convertStar2G( rho0IRncT(:, 1), rho0IRpw(:, 1), stars, ngdp, gdp )

     ! Read in valence density, because we only evaluate the valence density by a multiplication of the wave functions.
    allocate( r2Rho0MTVal( atoms%jmtd, 0:lathar%nlhd, atoms%ntype, 1) )
    r2Rho0MTVal = 0.
    rewind(1040)
    read(1040) r2Rho0MTVal
    allocate( rho0MTVal( atoms%jmtd, 0:lathar%nlhd, atoms%ntype, 1) )
    rho0MTVal(:, :, :, :) = 0.
    do iDtype = 1, atoms%ntype
      do ilh = 0, lathar%nlhd
        do imesh = 1, atoms%jri(iDtype)
          rho0MTVal(imesh, ilh, iDtype, 1) = r2Rho0MTVal(imesh, ilh, iDtype, 1) / atoms%rmsh(imesh, iDtype) / atoms%rmsh(imesh, iDtype)
        end do
      end do
    end do

    ! Calculate the 2 argument surface integrals of rho0Valence and grVeff0
    call CalcSurfIntIRDynMat( atoms, cell, ngdp, ngdp, gdp, gdp, rho0IRpw(:, 1), grVeff0IR, [0., 0., 0.], surfIntIR2Arg )
    call CalcSurfIntMTDynMat( atoms, lathar, clnu_atom, nmem_atom, mlh_atom, rho0MTVal(:, :, :, 1), grVeff0MT, surfIntMT2Arg)

    if (.false.) then
      write(*, '(a)') 'IR Surface integral'
      write(*, '(3(2(es16.8,1x),3x))') surfIntIR2Arg(1, :)
      write(*, '(3(2(es16.8,1x),3x))') surfIntIR2Arg(2, :)
      write(*, '(3(2(es16.8,1x),3x))') surfIntIR2Arg(3, :)

      write(*, '(a)') 'MT Surface integral'
      write(*, '(3(2(es16.8,1x),3x))') surfIntMT2Arg(1, :)
      write(*, '(3(2(es16.8,1x),3x))') surfIntMT2Arg(2, :)
      write(*, '(3(2(es16.8,1x),3x))') surfIntMT2Arg(3, :)
    end if

    ! Quantities for initialization
    allocate( lmpT(atoms%ntype) )
    lmpT(:) = 0
    do iDtype = 1, atoms%ntype
      lmpT(iDtype) = sum( [ ( (2 * oqn_l + 1)* nRadFun(oqn_l, iDtype), oqn_l = 0, atoms%lmax(iDtype) ) ] )
    end do ! itype
    lmpMax     = maxval( lmpT(:) )
    nRadFunMax = maxval( nRadFun(:, :) )
    coScale = 1
    allocate( veffUvIR(3, (coScale * atoms%lmaxd + 1)**2) )
    veffUvIR(:, :) = cmplx(0., 0.)

    allocate( surfIntVFast(maxval(nobd), 3))
    surfIntVFast(:, :) = cmplx(0., 0.)

    surfIntSum(:, :) = cmplx(0., 0.)
    surfIntSumMT(:, :) = cmplx(0., 0.)


    allocate( varphi1(atoms%jmtd, nRadFunMax, 0:atoms%lmaxd), varphi2(atoms%jmtd, nRadFunMax, 0:atoms%lmaxd) )
    allocate( ngoprI(atoms%nat) )
    allocate( grVeff0Varphi(2, atoms%jmtd, 0:(atoms%lmaxd + 3)**2 - 1, lmpMax) )
    allocate( hFullNoAbcofBK(lmpMax, lmpMax, 3, atoms%nat), overlapNoAbcofBK(lmpMax, lmpMax, 3, atoms%ntype) )
    allocate( lambdaKet(2, 0:atoms%lmaxd), lambdaBra(2, 0:atoms%lmaxd) )
    allocate( muKet(-atoms%lmaxd:atoms%lmaxd, -1:1), muBra(-atoms%lmaxd:atoms%lmaxd, -1:1) )
    allocate( varphiBra1(atoms%jmtd, 2, lmpMax, -1:1), varphiBra2(atoms%jmtd, 2, lmpMax, -1:1), &
            & varphiKet1(atoms%jmtd, 2, lmpMax, -1:1), varphiKet2(atoms%jmtd, 2, lmpMax, -1:1) )
    allocate( a( dimens%nvd, 0:dimens%lmd, atoms%nat), b(dimens%nvd, 0:dimens%lmd, atoms%nat), &
      & bascof_lo(3, -atoms%llod:atoms%llod, 4 * atoms%llod + 2, atoms%nlod, atoms%nat) )
    allocate( ab0cofKet(lmpMax, maxval(nobd)) )
    allocate(overlap(maxval(nobd)))
    allocate(surfIntMT(maxval(nobd)))
    iDatom = 0
    do iDtype = 1, atoms%ntype
      do iDeqat = 1, atoms%neq(iDtype)
        iDatom = iDatom + 1

        varphi1(:, :, :) = 0.
        varphi2(:, :, :) = 0.
        do oqn_l = 0, atoms%lmax(iDtype)
          do iradf = 1, nRadFun(oqn_l, iDtype)
            do imesh = 1, atoms%jri(iDtype)
              ! In Fleur the radial solutions u_lp are multiplied by a factor r by default to avoid an additional multiplication of the
              ! Jacobi determinant r^2 in radial integrals given spherical coordinates.
              varphi1(imesh, iradf, oqn_l) = rbas1(imesh, iradf, oqn_l, iDtype) / atoms%rmsh(imesh, iDtype)
              varphi2(imesh, iradf, oqn_l) = rbas2(imesh, iradf, oqn_l, iDtype) / atoms%rmsh(imesh, iDtype)
            end do ! imesh
          end do ! iradf
        end do ! oqn_l

        ngoprI(:) = 1

        do idirR = 1, 3
          ! Calculate the IR surf integral of Psi* grVeff0 Psi
          veffUvIR(:, :) = cmplx(0., 0.)
          call IRcoeffVeffUv( atoms, stars, cell, iDtype, iDatom, ngdp, coScale, gdp, grVeff0IR(:, idirR), veffUvIR )

          ! Action of gradient of potential onto MT basis functions without basis matching coefficients
          grVeff0Varphi = cmplx(0., 0.)
          call CalcFnsphVarphi( atoms, iDtype, 0, nRadFun, varphi1, varphi2, grVeff0MT(:, :, idirR, iDatom), grVeff0Varphi )

          ! Use the general 3 argument routine to check the result of the 2-argument surface integral
          hFullNoAbcofBK(:, :, :, :) = cmplx(0., 0.)
          overlapNoAbcofBK(:, :, :, :) = cmplx(0., 0.)
          chanMaxBra = 1
          chanMaxKet = 1
          lmaxBra = atoms%lmax(iDtype)
          muBra(:, :) = 0
          muKet(:, :) = 0
          mqn_m2PrBra = -1
          mqn_m2PrKet = -1
          varphiBra1 = 0.
          varphiKet1 = 0.
          varphiBra2 = 0.
          varphiKet2 = 0.
          do mqn_m = -atoms%lmax(iDtype), atoms%lmax(iDtype)
            muKet(mqn_m, -1) = mqn_m
            muBra(mqn_m, -1) = mqn_m
          end do ! mqn_m
          do oqn_l = 0, atoms%lmax(iDtype)
            lambdaKet(1, oqn_l) = oqn_l
            lambdaBra(1, oqn_l) = oqn_l
          end do ! oqn_l
          lmp = 0
          do oqn_l = 0, atoms%lmax(iDtype)
            do mqn_m = -oqn_l, oqn_l
              do iradf = 1, nRadFun(oqn_l, iDtype)
                lmp = lmp + 1
                varphiBra1(atoms%jri(iDtype), 1, lmp, -1) = varphi1(atoms%jri(iDtype), iradf, oqn_l )
                varphiKet1(atoms%jri(iDtype), 1, lmp, -1) = varphi1(atoms%jri(iDtype), iradf, oqn_l )
                varphiBra2(atoms%jri(iDtype), 1, lmp, -1) = varphi2(atoms%jri(iDtype), iradf, oqn_l )
                varphiKet2(atoms%jri(iDtype), 1, lmp, -1) = varphi2(atoms%jri(iDtype), iradf, oqn_l )
              end do ! iradf
            end do ! mqn_m
          end do ! oqn_l

          call PrepareMTSurfIntDM( atoms, iDtype, chanMaxBra, chanMaxKet, lmaxBra, lmpMax, mqn_m2PrBra, mqn_m2PrKet, muKet, muBra,     &
            & lambdaKet, lambdaBra, nRadFun, varphiKet1, varphiKet2, varphiBra1, varphiBra2, grVeff0Varphi, atoms%rmt(iDtype)**2, hFullNoAbcofBK, overlapNoAbcofBK )

          do ikpt = 1, kpts%nkpt
            surfIntVFast = cmplx(0., 0.)
            call calcSurfVeff( atoms, kpts, cell, dimens, ikpt, ikpt, coScale, gbas, veffUvIR, nv, nobd(ikpt, 1), &
              & mapGbas, z0(:, :, ikpt, 1), z0(:, :, ikpt, 1), iDtype, iDatom, surfIntVFast )

            ! Calculate the small basis matching coefficients at k.
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

            ! Calculate the large basis matching coefficients at k
            ab0cofKet(:, :) = cmplx(0., 0.)
            do iband = 1, nobd(ikpt, 1)
              lmp = 0
              lm  = 0
              do oqn_l = 0, atoms%lmax(iDtype)
                do mqn_m = - oqn_l, oqn_l
                  pMaxLocal = nRadFun(oqn_l, iDtype)
                  ! p = 1
                  ab0cofKet(lmp + 1, iband) = iu**oqn_l * dot_product( conjg(z0(:nv(1, ikpt), iband, ikpt, 1)), a(:nv(1, ikpt), lm, iDatom) )
                  ! p = 2
                  ab0cofKet(lmp + 2, iband) = iu**oqn_l * dot_product( conjg(z0(:nv(1, ikpt), iband, ikpt, 1)), b(:nv(1, ikpt), lm, iDatom) )
                  ! Add LO contributions
                  do iradf = 3, pMaxLocal
                    ! p = 1
                    ab0cofKet(lmp + 1, iband) = ab0cofKet(lmp + 1, iband) + iu**oqn_l * &
                      & dot_product( conjg(z0(nv(1, ikpt) + 1:nv(1, ikpt) + atoms%nlotot, iband, ikpt, 1)), &
                                                           & bascof_lo(1, mqn_m, :atoms%nlotot, iloTable(iradf, oqn_l, iDtype), iDatom) )
                    ! p = 2
                    ab0cofKet(lmp + 2, iband) = ab0cofKet(lmp + 2, iband) + iu**oqn_l * &
                      & dot_product( conjg(z0(nv(1, ikpt) + 1:nv(1, ikpt) + atoms%nlotot, iband, ikpt, 1)), &
                                                           & bascof_lo(2, mqn_m, :atoms%nlotot, iloTable(iradf, oqn_l, iDtype), iDatom) )
                    ! 2 < p < LOs for that l and that atom type
                    ab0cofKet(lmp + iradf, iband) = iu**oqn_l * dot_product( conjg(z0(nv(1, ikpt) + 1:nv(1, ikpt) + atoms%nlotot, iband, ikpt, 1)), &
                                                           & bascof_lo(3, mqn_m, :atoms%nlotot, iloTable(iradf, oqn_l, iDtype), iDatom) )
                  end do ! iradf

                  lm = lm + 1
                  lmp = lmp + pMaxLocal
                end do ! mqn_m
              end do !oqn_l
            end do ! iband

            do idirC = 1, 3

              ! Calculate the 3 -argument surface integral which has been prepared before
              overlap(:) = cmplx(0., 0.)
              surfIntMT(:) = cmplx(0., 0.)
              call CalcSintMT( ikpt, ikpt, idirC, iDatom, iDtype, lmpMax, nobd, hFullNoAbcofBK, overlapNoAbcofBK, ab0cofKet, ab0cofKet,   &
                                                                                                              & overlap, surfIntMT )
              do iband = 1, nobd(ikpt, 1)
                surfIntSum(idirR, idirC) = surfIntSum(idirR, idirC) + results%w_iks(iband, ikpt, 1) * surfIntVFast(iband, idirC)
                surfIntSumMT(idirR, idirC) = surfIntSumMT(idirR, idirC) + results%w_iks(iband, ikpt, 1) * surfIntMT(iband)
              end do ! iband
            end do ! idirC
          end do ! ikpt
        end do ! idirR

      end do ! iDeqat
    end do ! iDtype

    if (.false.) then
      write(*, '(a)') 'IR Surface integral benchmark'
      write(*, '(3(2(es16.8,1x),3x))') 2 * surfIntSum(1, :)
      write(*, '(3(2(es16.8,1x),3x))') 2 * surfIntSum(2, :)
      write(*, '(3(2(es16.8,1x),3x))') 2 * surfIntSum(3, :)

      write(*, '(a)') 'MT Surface integral benchmark'
      write(*, '(3(2(es16.8,1x),3x))') 2 * surfIntSumMT(1, :)
      write(*, '(3(2(es16.8,1x),3x))') 2 * surfIntSumMT(2, :)
      write(*, '(3(2(es16.8,1x),3x))') 2 * surfIntSumMT(3, :)
    end if

    ! Compare the results of the different modes of calculation
    do idirC = 1, 3
      do idirR = 1, 3
        if ( ( abs(surfIntIR2Arg(idirR, idirC) - 2 * surfIntSum(idirR, idirC)) >= 5e-5 ) &!.or. &
        !    & ( abs(surfIntMT2Arg(idirR, idirC) - 2 * surfIntSumMT(idirR, idirC)) >= 5e-7 ) .or. &
        !    & (abs(surfIntIR2Arg(idirR, idirC) + surfIntMT2Arg(idirR, idirC)) >=5e-5) .or. &
        !    & (abs( surfIntSum(idirR, idirC) + surfIntSumMT(idirR, idirC) ) >=5e-5) &
            &) passed = .false.
      end do ! idirR
    end do ! idirC

    if ( passed ) then
      write( logUnit, '(a)' ) '                                                     |__ passed!'
      write( logUnit, * )
    else
      write( logUnit, '(a)' ) '                                                     |__ failed!'
      write( logUnit, * )
      call JuDFT_warn('Test of dynamical matrix 2-argument surface integral', calledby='TestComp2ArgSFInts', &
        & hint='Activate and check debug output', &
        & file='m_jpTestDynMatSurf_mod.F90')
    end if
    !todo
    ! call with grVeff0Varphi and change z1 to z0 in there with testMode variable then at least on the diagonal it should be the same
    !call CalcSFintMTPsi1HepsPsiAndPsiHepsPsi1ExpCoeffVar( atoms, sym, dimens, usdus, kpts, cell, results, iqpt, iDtypeB, iDatomB, iDatomA, lmpMax, nRadFun, eig, varphi1, varphi2, hVarphi, &
    !& mapKpq2K, gBas, mapGbas, nv, kveclo, z0, lmpT, nobd, iloTable, surfInt )
  end subroutine TestComp2ArgSFInts

  !>-------------------------------------------------------------------------------------------------------------------------------------
  !> @author
  !> Christian-Roman Gerhorst
  !>
  !> @brief
  !> Compares output of different routines to calculate surface integrals dS Psi*_k'n (H - eps) Psi_k'n in the lR and the MT.
  !>
  !> @details
  !> For the Sternheimer SCC, surface integrals in the IR and the MT for different k-points and bands were programmed. For the
  !> dynamical matrix, we only need surface integrals that have equal bands but might have different k-points in the wave functions.
  !> Therefore, the routines can be simplified and we compare both versions in the IR and the MT, respectively.
  !>
  !> @note
  !> If the k-points are equal, the sum over all k-points and bands is zero in the IR and the MT. This does not hold for different
  !> k-points. Moreover, the IR and the MT surface integral are only the same, when we neglect the kinetic energy, because this
  !> implies the second derivative of the LAPW basis function which is not continious and therefore the difference of the surface
  !> integrals cannot vanish
  !>
  !> @param[out] atoms       : Atoms type, see types.f90.
  !> @param[out] cell        : Unit cell type, see types.f90.
  !> @param[out] sym         : Symmetries type, see types.f90.
  !> @param[out] stars       : Stars type, see types.f90.
  !> @param[out] lathar      : Lattice harmonics type, see types.f90.
  !> @param[out] input       : Input type, see types.f90.
  !> @param[out] dimens      : Dimension type, see types.f90.
  !> @param[out] kpts        : K-points type, see types.f90.
  !> @param[out] qpts        : Q-points type, see types.f90.
  !> @param[out] Veff0       : Type containing unperturbed output potentials of Fleur, see types.f90.
  !> @param[out] usdus       : Type containing quantities consisting of the radial solutions, see types.f90.
  !> @param[out] results     : Results type, see types.f90.
  !> @param[in]  logUnit     : Unit number for juPhon.log.
  !> @param[out] ngdp        : Number of G-vectors for potentials and densities.
  !> @param[out] kveclo      : Basis G-vectors of local orbitals.
  !> @param[out] nv          : Number of LAPW G-basis vectors for given k-point.
  !> @param[out] nRadFun     : Number of radial functions per orbital quantum number l and atom type.
  !> @param[out] iloTable    : Number of local orbital if orbital quantum number l, atom type and index p > 2 from nRadFun is given.
  !> @param[out] mapKpq2K    : For a given k-point index and q-point index, this array gives the k-point set index of the result k + q
  !>                           (already mapped back to Brillouin zone).
  !> @param[out] eig         : Contains Kohn\--Sham eigenvalues.
  !> @param[out] rbas1       : Large components of radial solution, its energy derivative and u_LO
  !> @param[out] rbas2       : Small components of radial solution, its energy derivative and u_LO
  !> @param[out] El          : Contains LAPW and LO energy parameters.
  !> @param[out] z0          : Kohn-Sham eigenvectors.
  !> @param[out] gdp         : G-vectors of potentials and densities.
  !> @param[out] mlh_atom    : Magnetic quantum number m of lattice harmonic members for every atom.
  !> @param[out] nmem_atom   : Number of lattice harmonic members for every atom.
  !> @param[out] clnu_atom   : Phase mediating between stars and plane waves.
  !> @param[out] nobd        : Number of occupied bands per k-point and spin
  !>
  !> @ todo complete documentation
  !>--------------------------------------------------------------------------------------------------------------------------------------
  subroutine TestComp3ArgSFInts( atoms, dimens, stars, sym, cell, kpts, qpts, input, lathar, usdus, Veff0, results, nRadFun, eig,  &
      & kpq2kPrVec, El, gdp, rbas1, rbas2, nmem_atom, mlh_atom, clnu_atom, z0, nv, gbas, mapGbas, nobd, kveclo, iloTable, ngdp,    &
      & mapKpq2K, logUnit )

    use m_types, only : t_atoms, t_dimension, t_stars, t_sym, t_cell, t_kpts, t_input, t_sphhar, t_usdus, t_potential, t_results, t_noco
     
    use m_jpSetupDynMatSF, only : PrepareMTSurfIntDM, CalcSIntMT, calcSurfVeff, CalcSintKinEnergOvl, IRcoeffVeffUv
    use m_jpSternhPulaySurface, only : CalcSintKinEps, CalcSfVeffFast
    use m_jpTestSternheimer, only : PrepareMTSurfInt, CalcSurfIntMT
    use m_juDFT_NOstopNO, only : juDFT_warn
    use m_jpConstants
    use m_abcof3
    use mod_juPhonUtils, only : Derivative, CalcChannelsGrFlpNat
    use m_jpPotDensHelper, only : ConvertStar2G
    use m_jpSetupDynMatHelper, only : CalcHnGrV0Varphi

    implicit none

    ! Type parameter
    type(t_atoms),                  intent(in)  :: atoms
    type(t_dimension),              intent(in)  :: dimens
    type(t_stars),                  intent(in)  :: stars
    type(t_sym),                    intent(in)  :: sym
    type(t_cell),                   intent(in)  :: cell
    type(t_kpts),                   intent(in)  :: kpts
    type(t_kpts),                   intent(in)  :: qpts
    type(t_input),                  intent(in)  :: input
    type(t_sphhar),                 intent(in)  :: lathar
    type(t_usdus),                  intent(in)  :: usdus
    type(t_potential),              intent(in)  :: Veff0
    type(t_results),                intent(in)  :: results

    ! Scalar parameter
    integer,                        intent(in)  :: ngdp
    integer,                        intent(in)  :: logUnit

    ! Array parameter
    integer,                        intent(in)  :: nRadFun(0:, :)
    real,                           intent(in)  :: El(:, 0:, :, :)
    real,                           intent(in)  :: rbas1(:, :, 0:, :)
    real,                           intent(in)  :: rbas2(:, :, 0:, :)
    integer,                        intent(in)  :: nmem_atom(0:, :)
    integer,                        intent(in)  :: mlh_atom(:,0:,:)
    complex,                        intent(in)  :: clnu_atom(:,0:,:)
    MCOMPLEX,                       intent(in)  :: z0(:,:,:,:)
    integer,                        intent(in)  :: nv(:, :)
    integer,                        intent(in)  :: gbas(:, :)
    integer,                        intent(in)  :: mapGbas(:, :, :)
    integer,                        intent(in)  :: nobd(:, :)
    integer,                        intent(in)  :: kveclo(:,:)
    integer,                        intent(in)  :: iloTable(:, 0:, :)
    integer,                        intent(in)  :: gdp(:, :)
    real,                           intent(in)  :: eig(:,:,:)
    integer,                        intent(in)  :: kpq2kPrVec(:, :, :)
    integer,                        intent(in)  :: mapKpq2K(:, :)

    ! Type variables
    type(od_inp)                                :: odi
    type(od_sym)                                :: ods
    type(t_noco)                                :: noco

    ! Scalar variables
    integer                                     :: iDtype
    integer                                     :: iDatom
    integer                                     :: idir
    integer                                     :: itype
    integer                                     :: oqn_l
    integer                                     :: lmpMax
    integer                                     :: ikpt
    integer                                     :: coScale
    integer                                     :: iband
    integer                                     :: iqpt
    integer                                     :: iBas
    integer                                     :: ikpq
    integer                                     :: iDeqat
    logical                                     :: passed
    integer                                     :: nmat
    integer                                     :: lmp
    integer                                     :: lm
    integer                                     :: mqn_m
    integer                                     :: pMaxLocal
    integer                                     :: iradf
    integer                                     :: imesh
    integer                                     :: nRadFunMax
    integer                                     :: ptsym
    integer                                     :: ilh
    integer                                     :: lm_pre
    integer                                     :: imem
    integer                                     :: chanMaxBra
    integer                                     :: chanMaxKet
    integer                                     :: lmaxBra
    integer                                     :: mqn_m2PrBra
    integer                                     :: mqn_m2PrKet

    ! Array variables
    complex,           allocatable              :: hFullNoAbcofBK(:, :, :, :)
    complex,           allocatable              :: overlapNoAbcofBK(:, :, :, :)
    complex,           allocatable              :: hFullNoAbcofBKbench(:, :, :, :)
    complex,           allocatable              :: overlapNoAbcofBKbench(:, :, :, :)
    integer,           allocatable              :: lmpT(:)
    complex,           allocatable              :: overlap(:)
    complex,           allocatable              :: surfIntMT(:)
    complex,           allocatable              :: overlapbench(:, :)
    complex,           allocatable              :: surfIntMTbench(:, :)
    complex,           allocatable              :: veffUvIR(:, :)
    complex,           allocatable              :: surfIntVFast(:, :)
    complex,           allocatable              :: surfIntVeffSternh(:, :, :)
    complex,           allocatable              :: surfIntTeps(:, :)
    complex,           allocatable              :: surfIntKinSternh(:, :)
    complex,           allocatable              :: surfInt(:, :)
    complex,           allocatable              :: surfIntOvlSternh(:, :)
    real,              allocatable              :: gBasBra(:, :)
    real,              allocatable              :: gBasKet(:, :)
    integer,           allocatable              :: ngoprI(:)
    complex,           allocatable              :: a(:, :, :)
    complex,           allocatable              :: b(:, :, :)
    complex,           allocatable              :: bascof_lo(:, :, :, :, :)
    complex,           allocatable              :: aKpq(:, :, :)
    complex,           allocatable              :: bKpq(:, :, :)
    complex,           allocatable              :: bascof_loKpq(:, :, :, :, :)
    complex,           allocatable              :: ab0cofBra(:, :)
    complex,           allocatable              :: ab0cofKet(:, :)
    real,          allocatable                  :: varphi1(:, :, :)
    real,          allocatable                  :: varphi2(:, :, :)
    real,          allocatable                  :: delrVarphi1(:, :, :)
    real,          allocatable                  :: delrVarphi2(:, :, :)
    integer,       allocatable                  :: grVarphiChLout(:, :)
    integer,       allocatable                  :: grVarphiChMout(:, :)
    real,          allocatable                  :: grVarphiCh1(:, :, :, :)
    real,          allocatable                  :: grVarphiCh2(:, :, :, :)
    complex,       allocatable                  :: vEff0MtSpH(:, :)
    complex,       allocatable                  :: hVarphi(:, :, :, :)
    real,          allocatable                  :: varphiBra1(:, :, :, :)
    real,          allocatable                  :: varphiKet1(:, :, :, :)
    real,          allocatable                  :: varphiBra2(:, :, :, :)
    real,          allocatable                  :: varphiKet2(:, :, :, :)
    integer,       allocatable                  :: lambdaKet(:, :)
    integer,       allocatable                  :: lambdaBra(:, :)
    integer,       allocatable                  :: muKet(:, :)
    integer,       allocatable                  :: muBra(:, :)
    complex,       allocatable                  :: vEff0NsphGrVarphi(:, :, :, :, :)
    complex,       allocatable                  :: r2grVeff0SphVarphi(:, :, :, :, :)
    complex,       allocatable                  :: r2grVeff0SphVarphiDummy(:, :, :, :, :)
    complex,       allocatable                  :: vpw_eff_uw(:)
    complex                                     :: intIrSf(3)
    complex                                     :: intMtSf(3)

    write( logUnit, '(a)' ) 'Test of general surface integral for dynamical matrix'
    write( logUnit, '(a)' ) '-----------------------------------------------------'

    ! Quantities for initialization
    allocate( lmpT(atoms%ntype) )
    lmpT(:) = 0
    do itype = 1, atoms%ntype
      lmpT(itype) = sum( [ ( (2 * oqn_l + 1)* nRadFun(oqn_l, itype), oqn_l = 0, atoms%lmax(itype) ) ] )
    end do ! itype
    lmpMax     = maxval( lmpT(:) )
    nRadFunMax = maxval( nRadFun(:, :) )
    coScale = 1

    ! Allocation of required in-loop arrays
    allocate( veffUvIR(3, (coScale * atoms%lmaxd + 1)**2) )
    allocate( hFullNoAbcofBK(lmpMax, lmpMax, 3, atoms%nat), overlapNoAbcofBK(lmpMax, lmpMax, 3, atoms%ntype) )
    allocate( hFullNoAbcofBKbench(lmpMax, lmpMax, 3, atoms%nat), overlapNoAbcofBKbench(lmpMax, lmpMax, 3, atoms%ntype) )
    allocate( gBasBra(3, dimens%nvd) )
    allocate( gBasKet(3, dimens%nvd) )
    allocate( surfIntKinSternh(maxval(nobd(:, 1)), maxval(nobd(:, 1))) )
    allocate( surfIntOvlSternh(maxval(nobd(:, 1)), maxval(nobd(:, 1))) )
    allocate( surfIntTeps(maxval(nobd), maxval(nobd)) )
    allocate( surfInt(maxval(nobd), maxval(nobd)) )
    allocate( surfIntVFast(maxval(nobd), 3))
    allocate( surfIntVeffSternh(maxval(nobd), maxval(nobd), 3) )
    allocate(overlap(maxval(nobd)))
    allocate(surfIntMT(maxval(nobd)))
    allocate(overlapbench(maxval(nobd), maxval(nobd)))
    allocate(surfIntMTbench(maxval(nobd), maxval(nobd)))
    allocate( ngoprI(atoms%nat) )
    allocate( a( dimens%nvd, 0:dimens%lmd, atoms%nat), b(dimens%nvd, 0:dimens%lmd, atoms%nat), &
      & bascof_lo(3, -atoms%llod:atoms%llod, 4 * atoms%llod + 2, atoms%nlod, atoms%nat) )
    allocate( aKpq( dimens%nvd, 0:dimens%lmd, atoms%nat), bKpq(dimens%nvd, 0:dimens%lmd, atoms%nat), &
      & bascof_loKpq(3, -atoms%llod:atoms%llod, 4 * atoms%llod + 2, atoms%nlod, atoms%nat) )
    allocate( ab0cofBra(lmpMax, maxval(nobd)), ab0cofKet(lmpMax, maxval(nobd)) )
    allocate( noco%alph(atoms%ntype), noco%beta(atoms%ntype) ) !Up to now those variables are only of dummy character
    allocate( varphi1(atoms%jmtd, nRadFunMax, 0:atoms%lmaxd), varphi2(atoms%jmtd, nRadFunMax, 0:atoms%lmaxd) )
    allocate( delrVarphi1( atoms%jmtd, nRadFunMax, 0:atoms%lmaxd), delrVarphi2( atoms%jmtd, nRadFunMax, 0:atoms%lmaxd) )
    allocate( vEff0MtSpH( atoms%jmtd, 0:dimens%lmd), vEff0NsphGrVarphi(2, atoms%jmtd, 0:(atoms%lmaxd + 3)**2 - 1, lmpMax, -1:1), &
            & r2grVeff0SphVarphi(2, atoms%jmtd, 0:(atoms%lmaxd + 3)**2 - 1, lmpMax, 3) )
    allocate( r2grVeff0SphVarphiDummy(2, atoms%jmtd, 0:(atoms%lmaxd + 3)**2 - 1, lmpMax, 3) )
    allocate( hVarphi(2, atoms%jmtd, 0:(atoms%lmaxd + 3)**2 - 1, lmpMax)  )
    allocate( grVarphiCh1(atoms%jmtd, 2, lmpMax, -1:1), grVarphiCh2(atoms%jmtd, 2, lmpMax, -1:1), &
            & grVarphiChLout(2, 0:atoms%lmaxd), grVarphiChMout(-atoms%lmaxd:atoms%lmaxd, -1:1) )
    allocate( varphiBra1(atoms%jmtd, 2, lmpMax, -1:1), varphiBra2(atoms%jmtd, 2, lmpMax, -1:1), &
            & varphiKet1(atoms%jmtd, 2, lmpMax, -1:1), varphiKet2(atoms%jmtd, 2, lmpMax, -1:1) )
    allocate( lambdaKet(2, 0:atoms%lmaxd), lambdaBra(2, 0:atoms%lmaxd) )
    allocate( muKet(-atoms%lmaxd:atoms%lmaxd, -1:1), muBra(-atoms%lmaxd:atoms%lmaxd, -1:1) )
    allocate( vpw_eff_uw( ngdp ) )

    lambdaKet = 0
    lambdaBra = 0
    muKet = 0
    muBra = 0

    ! Initialization
    passed = .true.
    iDatom = 0
    vEff0NsphGrVarphi(:, :, :, :, :)    = cmplx(0., 0.)
    r2grVeff0SphVarphi(:, :, :, :, :)   = cmplx(0., 0.)
    r2grVeff0SphVarphiDummy(:, :, :, :, :)   = cmplx(0., 0.)

    vpw_eff_uw = 0
    call convertStar2G( Veff0%vpw_uw(:, 1), vpw_eff_uw, stars, ngdp, gdp )

    do iDtype = 1, atoms%ntype
      do iDeqat = 1, atoms%neq(iDtype)
        iDatom = iDatom + 1

        veffUvIR(:, :) = cmplx(0., 0.)
        call IRcoeffVeffUv( atoms, stars, cell, iDtype, iDatom, ngdp, coScale, gdp, vpw_eff_uw, veffUvIR )

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

        ! Calculate the application of the gradient  onto the MT basis functions (matching coefficients have no spatial dependence) and determing its scattering channels.
        grVarphiChLout(:, :) = 0
        grVarphiChMout(:, :) = 0
        grVarphiCh1(:, :, :, :) = 0.
        grVarphiCh2(:, :, :, :) = 0.
        call CalcChannelsGrFlpNat( atoms, iDtype, nRadFun, varphi1, varphi2, delrVarphi1, delrVarphi2, grVarphiChLout, grVarphiChMout, &
                                                                                                        & grVarphiCh1, grVarphiCh2 )
        ! Calculate the application of the full Hamiltonian onto a basis function, interwine the non-spherical effective potential
        ! with the ket and do the same for the gradient of the spherical potential and the ket.
        hVarphi = cmplx(0.0, 0.0)
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
        call CalcHnGrV0Varphi( atoms, lathar, iDtype, iDatom, lmpMax, El, varphi1, varphi2, nRadFun, vEff0MtSpH, Veff0%vr(:, :, :, 1), clnu_atom, &
          & nmem_atom, mlh_atom, grVarphiCh1, grVarphiCh2, grVarphiChLout, grVarphiChMout, hVarphi, vEff0NsphGrVarphi, r2grVeff0SphVarphi, r2grVeff0SphVarphiDummy )

        ! Use the general MT surface integral routine to prepare the calculation of  into dS Psi* H - eps Psi
        hFullNoAbcofBK(:, :, :, :) = cmplx(0., 0.)
        overlapNoAbcofBK(:, :, :, :) = cmplx(0., 0.)
        hFullNoAbcofBKbench(:, :, :, :) = cmplx(0., 0.)
        overlapNoAbcofBKbench(:, :, :, :) = cmplx(0., 0.)
        chanMaxBra = 1
        chanMaxKet = 1
        lmaxBra = atoms%lmax(iDtype)
        muBra(:, :) = 0
        muKet(:, :) = 0
        mqn_m2PrBra = -1
        mqn_m2PrKet = -1
        varphiBra1 = 0.
        varphiKet1 = 0.
        varphiBra2 = 0.
        varphiKet2 = 0.
        do mqn_m = -atoms%lmax(iDtype), atoms%lmax(iDtype)
          muKet(mqn_m, -1) = mqn_m
          muBra(mqn_m, -1) = mqn_m
        end do ! mqn_m
        do oqn_l = 0, atoms%lmax(iDtype)
          lambdaKet(1, oqn_l) = oqn_l
          lambdaBra(1, oqn_l) = oqn_l
        end do ! oqn_l
        lmp = 0
        do oqn_l = 0, atoms%lmax(iDtype)
          do mqn_m = -oqn_l, oqn_l
            do iradf = 1, nRadFun(oqn_l, iDtype)
              lmp = lmp + 1
              varphiBra1(atoms%jri(iDtype), 1, lmp, -1) = varphi1(atoms%jri(iDtype), iradf, oqn_l )
              varphiKet1(atoms%jri(iDtype), 1, lmp, -1) = varphi1(atoms%jri(iDtype), iradf, oqn_l )
              varphiBra2(atoms%jri(iDtype), 1, lmp, -1) = varphi2(atoms%jri(iDtype), iradf, oqn_l )
              varphiKet2(atoms%jri(iDtype), 1, lmp, -1) = varphi2(atoms%jri(iDtype), iradf, oqn_l )
            end do ! iradf
          end do ! mqn_m
        end do ! oqn_l
        call PrepareMTSurfIntDM( atoms, iDtype, chanMaxBra, chanMaxKet, lmaxBra, lmpMax, mqn_m2PrBra, mqn_m2PrKet, muKet, muBra,     &
          & lambdaKet, lambdaBra, nRadFun, varphiKet1, varphiKet2, varphiBra1, varphiBra2, hVarphi, atoms%rmt(iDtype)**2, hFullNoAbcofBK, overlapNoAbcofBK )

        ! Here, the more specific routine, originating from the Sternheimer SCC is used for the MT calculation of into dS Psi* H - eps Psi
        hFullNoAbcofBKbench(:, :, :, :) = cmplx(0., 0.)
        overlapNoAbcofBKbench(:, :, :, :) = cmplx(0., 0.)
        call PrepareMTSurfInt( atoms, dimens, kpts, lathar, Veff0, iDtype, iDatom, nRadFun, El, rbas1, rbas2, nmem_atom, mlh_atom,&
          & clnu_atom, hFullNoAbcofBKbench, overlapNoAbcofBKbench, z0)

        write(*, *) 'Fix q after nobd is cleared before commit'
        do iqpt = 1, 1!4

          surfIntTeps = cmplx(0., 0.)
          surfInt = cmplx(0., 0.)
          surfIntKinSternh = cmplx(0., 0.)
          surfIntOvlSternh = cmplx(0., 0.)
          surfIntVFast = cmplx(0., 0.)
          surfIntVeffSternh = cmplx(0., 0.)
          intIrSf(:) = cmplx(0., 0.)
          intMtSf(:) = cmplx(0., 0.)
          overlap = cmplx(0., 0.)
          surfIntMT = cmplx(0., 0.)
          ! We do not want the local coordinate systems to be rotated for non-representative atoms constructing the matching coefficients.
          ngoprI(:) = 1

          do ikpt = 1, kpts%nkpt

            ikpq = mapKpq2K(ikpt, iqpt)

            ! Prepare different G-sets for the ket, the bra and for the kinetic energy term.
            gBasBra(:, :) = 0.
            gBasKet(:, :) = 0.
            ! NOTE: One has to care which k is here
            do iBas = 1, nv(1, ikpt)
              gBasKet(1:3, iBas) = real( gBas(1:3, mapGbas(iBas, ikpt, 1)) )
            end do ! iBas
            ! NOTE: One has to care which k is here
            do iBas = 1, nv(1, ikpq)
              gBasBra(1:3, iBas) = real( gBas(1:3, mapGbas(iBas, ikpq, 1)) + qpts%bk(1:3, iqpt) + kpq2kPrVec(1:3, ikpt, iqpt) )
            end do ! iBas

            ! Generalized routine to calcualte oint dS Psi* V Psi
            surfIntVFast = cmplx(0., 0.)
            call calcSurfVeff( atoms, kpts, cell, dimens, ikpq, ikpt, coScale, gbas, veffUvIR, nv, nobd(ikpt, 1), &
              & mapGbas, z0(:, :, ikpq, 1), z0(:, :, ikpt, 1), iDtype, iDatom, surfIntVFast )

            ! Sternheimer routine  to calcualte oint dS Psi* V Psi
            surfIntVeffSternh = cmplx(0., 0.)
            call calcSfVeffFast( atoms, kpts, qpts, cell, dimens, ikpt, ikpq, ikpt,  kpq2kPrVec, coScale, gBas, veffUvIR, nv, nobd(:, 1), nobd, mapGbas, z0,      &
              & iDtype, iDatom, surfIntVeffSternh)

            ! Calculate the basis matching coefficients at k + q = k' and the matching coefficients at k.
            nmat = nv(1, ikpq) + atoms%nlotot
            aKpq(:, :, :) = cmplx(0.0, 0.0)
            bKpq(:, :, :) = cmplx(0.0, 0.0)
            bascof_loKpq(:, :, :, :, :) = cmplx(0.0, 0.0)
            call abcof3( atoms%lmaxd, atoms%ntype, atoms%nat, sym%nop, dimens%nvd, dimens%jspd, 1, dimens%lmd, dimens%nbasfcn, &
              & atoms%llod, atoms%nlod, atoms%nlotot, sym%invtab, atoms%ntype, sym%mrot, ngoprI, atoms%taual, atoms%neq, atoms%lmax, &
              & atoms%rmt, cell%omtil, cell%bmat, cell%bbmat, kpts%bk(:, ikpq), gbas(1, mapGbas(:nv(1, ikpq), ikpq, 1)), &
              & gbas(2, mapGbas(:nv(1, ikpq), ikpq, 1)), gbas(3, mapGbas(:nv(1, ikpq), ikpq, 1)), nv(:, ikpq), nmat, &
              & usdus%us, usdus%dus, usdus%uds, usdus%duds, usdus%ddn, atoms%invsat, sym%invsatnr, usdus%ulos, usdus%uulon, usdus%dulon, &
              & usdus%dulos, atoms%llo, atoms%nlo, atoms%l_dulo, atoms%lapw_l, kveclo(:,ikpq), odi, ods, aKpq, bKpq, bascof_loKpq )

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

            ! Calculate the large matching coefficients as known from the Fleur abcof routine for the ket at k
            ab0cofKet(:, :) = cmplx(0., 0.)
            do iband = 1, nobd(ikpt, 1)
              lmp = 0
              lm  = 0
              do oqn_l = 0, atoms%lmax(iDtype)
                do mqn_m = - oqn_l, oqn_l
                  pMaxLocal = nRadFun(oqn_l, iDtype)
                  ! p = 1
                  ab0cofKet(lmp + 1, iband) = iu**oqn_l * dot_product( conjg(z0(:nv(1, ikpt), iband, ikpt, 1)), a(:nv(1, ikpt), lm, iDatom) )
                  ! p = 2
                  ab0cofKet(lmp + 2, iband) = iu**oqn_l * dot_product( conjg(z0(:nv(1, ikpt), iband, ikpt, 1)), b(:nv(1, ikpt), lm, iDatom) )
                  ! Add LO contributions
                  do iradf = 3, pMaxLocal
                    ! p = 1
                    ab0cofKet(lmp + 1, iband) = ab0cofKet(lmp + 1, iband) + iu**oqn_l * &
                      & dot_product( conjg(z0(nv(1, ikpt) + 1:nv(1, ikpt) + atoms%nlotot, iband, ikpt, 1)), &
                                                           & bascof_lo(1, mqn_m, :atoms%nlotot, iloTable(iradf, oqn_l, iDtype), iDatom) )
                    ! p = 2
                    ab0cofKet(lmp + 2, iband) = ab0cofKet(lmp + 2, iband) + iu**oqn_l * &
                      & dot_product( conjg(z0(nv(1, ikpt) + 1:nv(1, ikpt) + atoms%nlotot, iband, ikpt, 1)), &
                                                           & bascof_lo(2, mqn_m, :atoms%nlotot, iloTable(iradf, oqn_l, iDtype), iDatom) )
                    ! 2 < p < LOs for that l and that atom type
                    ab0cofKet(lmp + iradf, iband) = iu**oqn_l * dot_product( conjg(z0(nv(1, ikpt) + 1:nv(1, ikpt) + atoms%nlotot, iband, ikpt, 1)), &
                                                           & bascof_lo(3, mqn_m, :atoms%nlotot, iloTable(iradf, oqn_l, iDtype), iDatom) )
                  end do ! iradf

                  lm = lm + 1
                  lmp = lmp + pMaxLocal
                end do ! mqn_m
              end do !oqn_l
            end do ! iband

            ! Calculate the large matching coefficients as known from the Fleur abcof routine for the ket at k
            ab0cofBra(:, :) = cmplx(0., 0.)
            do iband = 1, nobd(ikpt, 1)
              lmp = 0
              lm  = 0
              do oqn_l = 0, atoms%lmax(iDtype)
                do mqn_m = - oqn_l, oqn_l
                  pMaxLocal = nRadFun(oqn_l, iDtype)
                  ! p = 1
                  ab0cofBra(lmp + 1, iband) = iu**oqn_l * dot_product( conjg(z0(:nv(1, ikpq), iband, ikpq, 1)), aKpq(:nv(1, ikpq), lm, iDatom) )
                  ! p = 2
                  ab0cofBra(lmp + 2, iband) = iu**oqn_l * dot_product( conjg(z0(:nv(1, ikpq), iband, ikpq, 1)), bKpq(:nv(1, ikpq), lm, iDatom) )
                  ! Add LO contributions
                  do iradf = 3, pMaxLocal
                    ! p = 1
                    ab0cofBra(lmp + 1, iband) = ab0cofBra(lmp + 1, iband) + iu**oqn_l * &
                      & dot_product( conjg(z0(nv(1, ikpq) + 1:nv(1, ikpq) + atoms%nlotot, iband, ikpq, 1)), &
                                                           & bascof_loKpq(1, mqn_m, :atoms%nlotot, iloTable(iradf, oqn_l, iDtype), iDatom) )
                    ! p = 2
                    ab0cofBra(lmp + 2, iband) = ab0cofBra(lmp + 2, iband) + iu**oqn_l * &
                      & dot_product( conjg(z0(nv(1, ikpq) + 1:nv(1, ikpq) + atoms%nlotot, iband, ikpq, 1)), &
                                                           & bascof_loKpq(2, mqn_m, :atoms%nlotot, iloTable(iradf, oqn_l, iDtype), iDatom) )
                    ! 2 < p < LOs for that l and that atom type
                    ab0cofBra(lmp + iradf, iband) = iu**oqn_l * dot_product( conjg(z0(nv(1, ikpq) + 1:nv(1, ikpq) + atoms%nlotot, iband, ikpq, 1)), &
                                                           & bascof_loKpq(3, mqn_m, :atoms%nlotot, iloTable(iradf, oqn_l, iDtype), iDatom) )
                  end do ! iradf

                  lm = lm + 1
                  lmp = lmp + pMaxLocal
                end do ! mqn_m
              end do !oqn_l
            end do ! iband

            do idir = 1, 3
              ! Generalized routine to calculate the overlap in the IR of the wave function and the kinetic energy term
              surfIntTeps = cmplx(0.0, 0.0)
              surfInt = cmplx(0.0, 0.0)
              call calcSintKinEnergOvl( atoms, cell, kpts, results, iDtype, iDatom, ikpq, ikpt, idir, &
                & nv(1, :), nobd(ikpt, 1), gBasBra, gBasKet, gBasKet, gBasBra, z0(:, :, ikpq, 1), z0(:, :, ikpt, 1), surfIntTeps, surfInt )

              ! Sternheimer routine to calculate the overlap in the IR of the wave function and the kinetic energy term
              surfIntKinSternh = cmplx(0., 0.)
              surfIntOvlSternh = cmplx(0., 0.)
              call calcSintKinEps( atoms, cell, kpts, qpts, iDtype, iDatom, nobd(ikpt, 1), nobd(ikpt, 1), ikpt, ikpq, iqpt, idir, &
                & nv(1, :), gbas, mapGbas, z0(:, :, ikpq, 1), z0(:, :, ikpt, 1), surfIntKinSternh,          &
                & surfIntOvlSternh, kpq2kPrVec )


              ! Generalized routine to calculate the surface integral oint dS Psi* (H - eps) Psi in the MT
              overlap(:) = cmplx(0., 0.)
              surfIntMT(:) = cmplx(0., 0.)
              call CalcSintMT( ikpq, ikpt, idir, iDatom, iDtype, lmpMax, nobd, hFullNoAbcofBK, overlapNoAbcofBK, ab0cofBra, ab0cofKet,   &
                                                                                                              & overlap, surfIntMT )

              ! Sternheimer routine to calculate the surface integral oint dS Psi* (H - eps) Psi in the MT
              overlapbench(:, :) = cmplx(0., 0.)
              surfIntMTbench(:, :) = cmplx(0., 0.)
              call CalcSurfIntMT(atoms, dimens, sym, cell, kpts, input, usdus, ikpq, ikpt, idir, iDatom, iDtype, nv, gbas, mapGbas, z0(:, :, ikpq, 1), z0(:, :, ikpt, 1), nobd, &
                & nobd(:, 1), kveclo, nRadFun, iloTable, hFullNoAbcofBKbench, overlapNoAbcofBKbench, overlapbench, surfIntMTbench)

              do iband = 1, nobd(ikpt, 1)

                ! we only compare the IR routines or the MT routines, respectively. The reasons is given in the above note within
                ! the documentation of this routine
             !   if ( abs( (surfIntTeps(iband, iband) + surfIntVFast(iband, idir) - eig(iband, ikpt, 1) * surfInt(iband, iband)) &
             !          & - (surfIntKinSternh(iband, iband) + surfIntVeffSternh(iband, iband, idir) &
             !          &   - eig(iband, ikpt, 1) * surfIntOvlSternh(iband, iband)) ) > 1e-7 ) passed = .false.
             ! cannot be compared because of different kinetic energy
                if ( abs( (surfIntVFast(iband, idir) - eig(iband, ikpt, 1) * surfInt(iband, iband)) &
                       & - (surfIntVeffSternh(iband, iband, idir) &
                       &   - eig(iband, ikpt, 1) * surfIntOvlSternh(iband, iband)) ) > 1e-5 ) passed = .false.

                if ( abs((surfIntMT(iband) - eig(iband, ikpt, 1) * overlap(iband)) - (surfIntMTbench(iband, iband) &
                  & - eig(iband, ikpt, 1) * overlapbench(iband, iband))) > 1e-5 ) passed = .false.

                ! Output for debugging
                if (.false.) then
                  write(3010, '(3i8, 2f15.8)') idir, ikpt, iband, &
                                              &  surfIntVFast(iband, idir)- eig(iband, ikpt, 1) * surfInt(iband, iband)
                  write(3011, '(3i8, 2f15.8)') idir, ikpt, iband,  &
                                    &  surfIntVeffSternh(iband, iband, idir) - eig(iband, ikpt, 1) * surfIntOvlSternh(iband, iband)
                  write(3012, '(3i8, 2f15.8)') idir, ikpt, iband, -(surfIntMT(iband) - eig(iband, ikpt, 1) * overlap(iband))
                  write(3013, '(3i8, 2f15.8)') idir, ikpt, iband, &
                                                & -(surfIntMTbench(iband, iband) - eig(iband, ikpt, 1) * overlapbench(iband, iband))
                end if

                ! Summing up the bands and k-points for all displacement directions
                intIrSf(idir) = intIrSf(idir) + results%w_iks(iband, ikpt, 1) * (surfIntTeps(iband, iband) + surfIntVFast(iband, idir) - eig(iband, ikpt, 1) * surfInt(iband, iband))
                intMtSf(idir) = intMtSf(idir) + results%w_iks(iband, ikpt, 1) * (surfIntMT(iband) - eig(iband, ikpt, 1) * overlap(iband))
              end do ! iband
            end do ! idir
          end do ! ikpt

          if (iqpt == 1) then
            ! Should be quite accurate because we do not integrate over the radial mesh
            if ( any(abs(intIRSf(:)) > 1e-9) .or. any(abs(intMTSf(:)) > 1e-9) ) passed = .false.
          end if

          ! Output the displacement directions
          if (.false.) then
            write(*, *) intIRSf(1:3)
            write(*, *) intMtSf(1:3)
          end if

        end do ! iDeqat
      end do ! iDtype
    end do ! iqpt

    if ( passed ) then
      write( logUnit, '(a)' ) '                                                    |__ passed!'
      write( logUnit, * )
    else
      write( logUnit, '(a)' ) '                                                    |__ failed!'
      write( logUnit, * )
      call JuDFT_warn('Test of general surface integral for dynamical matrix failed.', calledby='TestComp3ArgSFInts', hint='Uncomment debug output to analyze error!', &
        & file='m_jpTestDynMatSurf_mod.F90')
    end if

  end subroutine TestComp3ArgSFInts

  ! Tests genereal 3argument surface routine in which one of the arguments is a gradient of a wavefunction. They should be the same
  ! as two argument surface integrals with gradrho and veff
  subroutine TestComp2ArgGrSFInts( atoms, dimens, lathar, cell, kpts, qpts, stars, Veff0, results, sym, usdus, ngdp, logUnit, gdp, nv, gbas, mapGbas, nobd, z0, nmem_atom, mlh_atom, clnu_atom, nRadFun, rbas1, rbas2, kveclo, iloTable, mapKpq2K, kpq2kPrVec, eig )

    use m_types, only : t_atoms, t_dimension, t_sphhar, t_cell, t_kpts, t_stars, t_potential, t_results, t_sym, t_usdus
    use m_jpConstants, only : iu, Tmatrix
    use mod_juPhonUtils, only : Derivative, CalcChannelsGrFlpNat, CalcGrR2FinLH
    use m_jpSetupDynMatSF, only :  IRcoeffVeffUv, CalcSurfVeff, PrepareMTSurfIntDM, CalcSIntMT, CalcSintKinEnergOvl, CalcSFintIRPsi1HepsPsi, CalcSFintIRPsiHepsPsi1
     
    use m_jpSetupDynMatHelper, only : CalcFnsphGrVarphi, CalcSurfIntIRDynMat, CalcSurfIntMTDynMat
    use m_jpSetupDynMatHelper, only : CalcFnsphVarphi
    use m_abcof3
    use m_juDFT_NOstopNO, only : juDFT_warn
    use m_jpPotDensHelper, only : convertStar2G

    implicit none

    ! Type parameter
    type(t_atoms),                  intent(in) :: atoms
    type(t_dimension),              intent(in) :: dimens
    type(t_sphhar),                 intent(in) :: lathar
    type(t_cell),                   intent(in) :: cell
    type(t_kpts),                   intent(in) :: kpts
    type(t_kpts),                   intent(in) :: qpts
    type(t_stars),                  intent(in) :: stars
    type(t_potential),              intent(in) :: Veff0
    type(t_results),                intent(in) :: results
    type(t_sym),                    intent(in) :: sym
    type(t_usdus),                  intent(in) :: usdus

    ! Scalar variable
    integer,                        intent(in) :: ngdp
    integer,                       intent(in)  :: logUnit

    ! Array parameter
    real,                           intent(in)  :: rbas1(:, :, 0:, :)
    real,                           intent(in)  :: rbas2(:, :, 0:, :)
    integer,                        intent(in) :: gdp(:, :)
    integer,                        intent(in)  :: nv(:, :)
    integer,                        intent(in)  :: nobd(:, :)
    integer,                        intent(in)  :: gbas(:, :)
    integer,                        intent(in)  :: mapGbas(:, :, :)
    integer,                        intent(in)  :: nmem_atom(0:, :)
    integer,                        intent(in)  :: mlh_atom(:,0:,:)
    complex,                        intent(in)  :: clnu_atom(:,0:,:)
    integer,                        intent(in)  :: nRadFun(0:, :)
    MCOMPLEX,                       intent(in)  :: z0(:,:,:,:)
    integer,                        intent(in)  :: kveclo(:,:)
    integer,                        intent(in)  :: iloTable(:, 0:, :)
    integer,                        intent(in)  :: mapKpq2K(:, :)
    integer,                        intent(in)  :: kpq2kPrVec(:, :, :)
    real,                          intent(in)  :: eig(:, :, :)

    ! Type variables
    type(od_inp)                                :: odi
    type(od_sym)                                :: ods

    ! Scalar parameter
    integer                                    :: iG
    integer                                    :: iDatom
    integer                                    :: iDeqat
    integer                                    :: iDtype
    integer                                    :: coScale
    integer                                    :: ikpt
    integer                                    :: iBas
    integer                                    :: iband
    integer                                    :: idirR
    integer                                    :: idirC
    integer                                    :: oqn_l
    integer                                    :: lm_pre
    integer                                    :: mqn_m
    integer                                    :: lm
    integer                                    :: imesh
    integer                                    :: lmpMax
    integer                                    :: nRadFunMax
    integer                                    :: iradf
    integer                                     :: ptsym
    integer                                     :: ilh
    integer                                     :: chanMaxBra
    integer                                     :: chanMaxKet
    integer                                     :: lmaxBra
    integer                                     :: mqn_m2PrKet
    integer                                     :: lmp
    integer                                     :: mqn_m2PrBra
    integer                                     :: imem
    integer                                     :: nmat
    integer                                     :: pMaxLocal
    logical                                     :: passed = .true.
    logical                                     :: testMode = .true.

    ! Array variable
    complex,           allocatable             :: rho0IRpw(:)
    complex,           allocatable             :: Veff0IRpw(:)
    complex,           allocatable             :: surfIntVFast(:, :)
    complex,           allocatable             :: grRho0IR(:, :)
    complex,           allocatable             :: veffUvIR(:, :)
    complex,           allocatable             :: gradz0(:, :, :)
    complex,           allocatable             :: grRho0MT(:, :, :, :)
    real,              allocatable             :: r2Rho0MTVal(:, :, :, :)
    complex,           allocatable             :: r2GrRho0MT(:, :, :, :)
    integer,           allocatable             :: lmpT(:)
    complex,           allocatable             :: surfIntMT(:)
    complex,           allocatable             :: hFullNoAbcofBK(:, :, :, :)
    complex,           allocatable             :: overlapNoAbcofBK(:, :, :, :)
    complex,           allocatable             :: ab0cofKet(:, :)
    complex,           allocatable             :: a(:, :, :)
    complex,           allocatable             :: b(:, :, :)
    complex,           allocatable             :: bascof_lo(:, :, :, :, :)
    complex,           allocatable             :: veff0Varphi(:, :, :, :)
    integer,           allocatable             :: muKet(:, :)
    real,              allocatable             :: grVarphiCh1(:, :, :, :)
    real,              allocatable             :: grVarphiCh2(:, :, :, :)
    integer,           allocatable             :: lambdaKet(:, :)
    real,              allocatable             :: varphi1(:, :, :)
    real,              allocatable             :: varphi2(:, :, :)
    real,              allocatable             :: varphiKet1(:, :, :, :)
    real,              allocatable             :: varphiKet2(:, :, :, :)
    integer,           allocatable             :: grVarphiChLout(:, :)
    integer,           allocatable             :: grVarphiChMout(:, :)
    real,              allocatable             :: delrVarphi1(:, :, :)
    real,              allocatable             :: delrVarphi2(:, :, :)
    complex,           allocatable             :: vEff0MtSpH(:, :)
    complex,           allocatable             :: overlap(:)
    integer,           allocatable             :: ngoprI(:)
    complex,           allocatable             :: vEff0GrVarphi(:, :, :, :, :)
    complex,           allocatable             :: surfIntTeps(:, :)
    complex,           allocatable             :: surfInt(:, :)
    real,              allocatable             :: gBasMapped(:, :)
    complex,           allocatable             :: rho0IRnCt(:, :)
    complex                                    :: intMtSfNat(-1:1, 3)
    complex                                    :: intMtSf(3, 3)
    real                                       :: Gext(3)
    real                                       :: kExt(3)
    complex                                    :: surfIntIR2Arg(3, 3)
    complex                                    :: surfIntMT2Arg(3, 3)
    complex                                    :: surfIntSum(3, 3)
    complex                                    :: surfIntSumO(3, 3)
    complex                                    :: overlapSumNat(-1:1, 3)
    complex                                    :: overlapSum(3, 3)

    write( logUnit, '(a)' ) 'Test of surface integrals containing wave function gradients'
    write( logUnit, '(a)' ) '------------------------------------------------------------'

    ! This switch is for the surface integrals containing psi1 to make them use gradPsi0 instead of psi1
    testMode = .true.

    !-------------------------------------------------------------------------------------------------------------------------------
    ! Use 2-argument surface integral to calculate grRho vEff in the IR and the MT
    !-------------------------------------------------------------------------------------------------------------------------------

    ! Initialization
    surfIntIR2Arg(:, :) = cmplx(0., 0.)
    surfIntMT2Arg(:, :) = cmplx(0., 0.)
    surfIntSum(:, :) = cmplx(0., 0.)
    surfIntSumO(:, :) = cmplx(0., 0.)

    allocate(rho0IRpw(ngdp), Veff0IRpw(ngdp))
    rho0IRpw = cmplx(0., 0.)
    Veff0IRpw = cmplx(0., 0.)

    allocate(rho0IRncT(stars%n3d, 1))
    rho0IRncT(:, :) = cmplx(0.,0.)
    rewind(7800)
    read(7800) rho0IRncT
    call convertStar2G( rho0IRncT(:, 1), rho0IRpw, stars, ngdp, gdp )
    call convertStar2G( Veff0%vpw_uw(:, 1), vEff0IRpw , stars, ngdp, gdp )

    ! Perform analytical gradient of density in the interstitial region
    allocate( grRho0IR(ngdp, 3) )
    grRho0IR(:, :) = cmplx(0., 0.)
    do iG = 1, ngdp
      Gext(1:3) = matmul(cell%bmat(1:3, 1:3), gdp(1:3, iG))
      grRho0IR(iG, 1:3) = iu * Gext(1:3) * rho0IRpw(iG)
    end do ! iG


    ! Calculate IR surface integral of grRho veff0
    call CalcSurfIntIRDynMat( atoms, cell, ngdp, ngdp, gdp, gdp, Veff0IRpw, grRho0IR, [0., 0., 0.], surfIntIR2Arg )

    ! Initialization
    allocate( grRho0MT( atoms%jmtd, ( atoms%lmaxd + 2 )**2, 3, atoms%nat) )
    allocate( r2Rho0MTVal( atoms%jmtd, 0:lathar%nlhd, atoms%ntype, 1) )
    grRho0MT(:, :, :, :) = cmplx(0., 0.)
    r2Rho0MTVal(:, :, :, :) = cmplx(0., 0.)

    ! Read in valence density which is implicitely multiplied with r^2 already in FLEUR.
    rewind(1040)
    read(1040) r2Rho0MTVal

    ! Calculate numerical muffin-tin gradient of valence density multiplied with r^2
    call calcGrR2FinLH( atoms, lathar, clnu_atom, nmem_atom, mlh_atom, r2Rho0MTVal(:, :, :, 1), r2GrRho0MT )

    ! The factor r^2 is divided out to result in the pure gradient of the density
    do idirR = 1, 3
      iDatom = 0
      do iDtype = 1, atoms%ntype
        do iDeqat = 1, atoms%neq(iDtype)
          iDatom = iDatom + 1
          do oqn_l = 0, atoms%lmax(iDtype)
            lm_pre = oqn_l * (oqn_l + 1) + 1
            do mqn_m = -oqn_l, oqn_l
              lm = lm_pre + mqn_m
              do imesh = 1, atoms%jri(iDtype)
                grRho0MT(imesh, lm, idirR, iDatom) = r2GrRho0MT(imesh, lm, iDatom, idirR) / atoms%rmsh(imesh, iDtype) &
                                                                                                       & / atoms%rmsh(imesh, iDtype)
              end do ! imesh
            end do ! mqn_m
          end do ! oqn_l
        end do ! ieqat
      end do ! itype
    end do ! idirR

    ! Calculate the MT surface integral of the gradient of rho and the effective potential.
    call CalcSurfIntMTDynMat( atoms, lathar, clnu_atom, nmem_atom, mlh_atom, Veff0%vr(:, :, :, 1), grRho0MT, surfIntMT2Arg)

    if (.false.) then
      write(*, '(a)') 'IR Surface integral'
      write(*, '(3(2(es16.8,1x),3x))') surfIntIR2Arg(1, :)
      write(*, '(3(2(es16.8,1x),3x))') surfIntIR2Arg(2, :)
      write(*, '(3(2(es16.8,1x),3x))') surfIntIR2Arg(3, :)

      write(*, '(a)') 'MT Surface integral'
      write(*, '(3(2(es16.8,1x),3x))') surfIntMT2Arg(1, :)
      write(*, '(3(2(es16.8,1x),3x))') surfIntMT2Arg(2, :)
      write(*, '(3(2(es16.8,1x),3x))') surfIntMT2Arg(3, :)
    end if

    ! Initialization
    coScale = 1
    allocate( veffUvIR(3, (coScale * atoms%lmaxd + 1)**2) )
    allocate( gradz0(dimens%nvd, maxval(nobd), 3) )
    allocate( surfIntVFast(maxval(nobd), 3))

    ! Quantities for initialization
    allocate( lmpT(atoms%ntype) )
    lmpT(:) = 0
    do iDtype = 1, atoms%ntype
      lmpT(iDtype) = sum( [ ( (2 * oqn_l + 1)* nRadFun(oqn_l, iDtype), oqn_l = 0, atoms%lmax(iDtype) ) ] )
    end do ! itype
    lmpMax     = maxval( lmpT(:) )
    nRadFunMax = maxval( nRadFun(:, :) )
    allocate(surfIntMT(maxval(nobd)), overlap(maxval(nobd)) )
    allocate( hFullNoAbcofBK(lmpMax, lmpMax, 3, atoms%nat), overlapNoAbcofBK(lmpMax, lmpMax, 3, atoms%ntype) )
    allocate( ab0cofKet(lmpMax, maxval(nobd)) )
    allocate( a( dimens%nvd, 0:dimens%lmd, atoms%nat), b(dimens%nvd, 0:dimens%lmd, atoms%nat), &
      & bascof_lo(3, -atoms%llod:atoms%llod, 4 * atoms%llod + 2, atoms%nlod, atoms%nat) )
    allocate( veff0Varphi(2, atoms%jmtd, 0:(atoms%lmaxd + 3)**2 - 1, lmpMax) )
    allocate( varphiKet1(atoms%jmtd, 2, lmpMax, -1:1), varphiKet2(atoms%jmtd, 2, lmpMax, -1:1) )
    allocate( lambdaKet(2, 0:atoms%lmaxd) )
    allocate( muKet(-atoms%lmaxd:atoms%lmaxd, -1:1) )
    allocate( grVarphiCh1(atoms%jmtd, 2, lmpMax, -1:1), grVarphiCh2(atoms%jmtd, 2, lmpMax, -1:1), &
            & grVarphiChLout(2, 0:atoms%lmaxd), grVarphiChMout(-atoms%lmaxd:atoms%lmaxd, -1:1) )
    allocate( delrVarphi1( atoms%jmtd, nRadFunMax, 0:atoms%lmaxd), delrVarphi2( atoms%jmtd, nRadFunMax, 0:atoms%lmaxd) )
    allocate( ngoprI(atoms%nat) )
    allocate( varphi1(atoms%jmtd, nRadFunMax, 0:atoms%lmaxd), varphi2(atoms%jmtd, nRadFunMax, 0:atoms%lmaxd) )
    allocate( vEff0MtSpH( atoms%jmtd, 0:dimens%lmd))
    allocate( surfIntTeps(maxval(nobd), maxval(nobd)) )
    allocate( surfInt(maxval(nobd), maxval(nobd)) )
    allocate( gBasMapped(3, dimens%nvd) )
    ngoprI(:) = 1
    surfIntTeps(:, :) = cmplx(0., 0.)
    surfInt(:, :) = cmplx(0., 0.)
    gBasMapped(:, :) = 0.


    iDatom = 0
    do iDtype = 1, atoms%ntype
      do iDeqat = 1, atoms%neq(iDtype)
        iDatom = iDatom + 1

        ! Prepare calculation of surface integral with Psi* Psi and Veff
        veffUvIR(:, :) = cmplx(0., 0.)
        call IRcoeffVeffUv( atoms, stars, cell, iDtype, iDatom, ngdp, coScale, gdp, vEff0IRpw, veffUvIR )

        varphi1(:, :, :) = 0.
        varphi2(:, :, :) = 0.
        delrVarphi1(:, :, :) = cmplx(0., 0.)
        delrVarphi2(:, :, :) = cmplx(0., 0.)
        do oqn_l = 0, atoms%lmax(iDtype)
          do iradf = 1, nRadFun(oqn_l, iDtype)
            do imesh = 1, atoms%jri(iDtype)
              ! In Fleur the radial solutions u_lp are multiplied by a factor r by default to avoid an additional multiplication of the
              ! Jacobi determinant r^2 in radial integrals given spherical coordinates.
              varphi1(imesh, iradf, oqn_l) = rbas1(imesh, iradf, oqn_l, iDtype) / atoms%rmsh(imesh, iDtype)
              varphi2(imesh, iradf, oqn_l) = rbas2(imesh, iradf, oqn_l, iDtype) / atoms%rmsh(imesh, iDtype)
            end do ! imesh
            call Derivative( varphi1(1:atoms%jri(iDtype), iradf, oqn_l), iDtype, atoms, delrVarphi1(1:atoms%jri(iDtype), iradf, oqn_l) )
            call Derivative( varphi2(1:atoms%jri(iDtype), iradf, oqn_l), iDtype, atoms, delrVarphi2(1:atoms%jri(iDtype), iradf, oqn_l) )
          end do ! iradf
        end do ! oqn_l

        ! Calculate the application of the gradient  onto the MT basis functions (matching coefficients
        ! have no spatial dependence) and determing its scattering channels.
        grVarphiChLout(:, :) = 0
        grVarphiChMout(:, :) = 0
        grVarphiCh1(:, :, :, :) = 0.
        grVarphiCh2(:, :, :, :) = 0.
        call CalcChannelsGrFlpNat( atoms, iDtype, nRadFun, varphi1, varphi2, delrVarphi1, delrVarphi2, grVarphiChLout, grVarphiChMout, &
                                                                                                        & grVarphiCh1, grVarphiCh2 )

        ! Calculate the interwining of the effective potential with the wave function in the ket which is of the array form as the
        ! applicaiton of H onto varphi.
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
        veff0Varphi(:, :, :, :) = cmplx(0., 0.)
        call CalcFnsphVarphi( atoms, iDtype, 0, nRadFun, varphi1, varphi2, vEff0MtSpH(:, :), veff0Varphi )

        ! Set the parameters for the genreral muffin-tin surface integral in which the gradient of the wavefunction is made in the
        ! bra and in the ket there is the normal wave function
        chanMaxBra = 2
        chanMaxKet = 1
        lmaxBra = atoms%lmax(iDtype)
        muKet(:, :) = 0
        mqn_m2PrKet = -1
        varphiKet1 = 0.
        varphiKet2 = 0.
        do mqn_m = -atoms%lmax(iDtype), atoms%lmax(iDtype)
          muKet(mqn_m, -1) = mqn_m
        end do ! mqn_m
        do oqn_l = 0, atoms%lmax(iDtype)
          lambdaKet(1, oqn_l) = oqn_l
        end do ! oqn_l
        lmp = 0
        do oqn_l = 0, atoms%lmax(iDtype)
          do mqn_m = -oqn_l, oqn_l
            do iradf = 1, nRadFun(oqn_l, iDtype)
              lmp = lmp + 1
              varphiKet1(atoms%jri(iDtype), 1, lmp, -1) = varphi1(atoms%jri(iDtype), iradf, oqn_l )
              varphiKet2(atoms%jri(iDtype), 1, lmp, -1) = varphi2(atoms%jri(iDtype), iradf, oqn_l )
            end do ! iradf
          end do ! mqn_m
        end do ! oqn_l

        intMtSfNat(:, :) = cmplx(0., 0.)
        overlapSumNat(:, :) = cmplx(0., 0.)
        do mqn_m2PrBra = -1, 1
        !todo this is not optimal because interstitial integrals are executed threee times
          hFullNoAbcofBK(:, :, :, :) = cmplx(0., 0.)
          overlapNoAbcofBK(:, :, :, :) = cmplx(0., 0.)
          ! Calculate MT surface integral of basis functions grVarphi Veff Varphi
          call PrepareMTSurfIntDM( atoms, iDtype, chanMaxBra, chanMaxKet, lmaxBra, lmpMax, mqn_m2PrBra, mqn_m2PrKet, muKet, grVarphiChMout,     &
            & lambdaKet, grVarphiChLout, nRadFun, varphiKet1, varphiKet2, grVarphiCh1, grVarphiCh2, veff0Varphi, atoms%rmt(iDtype)**2, hFullNoAbcofBK, overlapNoAbcofBK )
          overlap = cmplx(0., 0.)
          surfIntMT = cmplx(0., 0.)
          do ikpt = 1, kpts%nkpt

            ! Generate the analytical derivative of the wave functions
            surfIntVFast = cmplx(0., 0.)
            gradz0(:, :, :) = cmplx(0., 0.)
            kExt(1:3) = matmul(cell%bmat, kpts%bk(1:3, ikpt))
            gExt(:) = 0.
            gBasMapped(:, :) = 0.
            do iBas = 1, nv(1, ikpt) !+ atoms%nlotot
              gExt(1:3) = matmul( cell%bmat(1:3, 1:3), gBas(1:3, mapGbas(iBas, ikpt, 1)))
              do iband = 1, nobd(ikpt, 1)
                do idirR = 1, 3
                  gradz0(iBas, iband, idirR) = iu * ( kExt(idirR) + gExt(idirR) ) * z0(iBas, iband, ikpt, 1)
                end do ! idirR
              end do ! iband
              gBasMapped(1:3, iBas) = real( gBas(1:3, mapGbas(iBas, ikpt, 1)) )
            end do ! iBas

           !todo directions correct?
            do idirR = 1, 3
            ! The surface integral of gradPsi Veff Psi with the old method
              call calcSurfVeff( atoms, kpts, cell, dimens, ikpt, ikpt, coScale, gbas, veffUvIR, nv, nobd(ikpt, 1), &
                & mapGbas, gradz0(:, :, idirR), z0(:, :, ikpt, 1), iDtype, iDatom, surfIntVFast )
              do idirC = 1, 3
              ! Old benchmark routine from the Sternheimer equation calculating gradPsi Veff Psi, and gradPsi Psi
              surfIntTeps = cmplx(0.0, 0.0)
              surfInt = cmplx(0.0, 0.0)
              call calcSintKinEnergOvl( atoms, cell, kpts, results, iDtype, iDatom, ikpt, ikpt, idirC, &
                & nv(1, :), nobd(ikpt, 1), gBasMapped, gBasMapped, gBasMapped, gBasMapped, gradz0(:, :, idirR), z0(:, :, ikpt, 1), surfIntTeps, surfInt )
                do iband = 1, nobd(ikpt, 1)
                  surfIntSum(idirR, idirC) = surfIntSum(idirR, idirC) + results%w_iks(iband, ikpt, 1) * surfIntVFast(iband, idirC)
                  surfIntSumO(idirR, idirC) = surfIntSumO(idirR, idirC) + results%w_iks(iband, ikpt, 1) * surfInt(iband, iband)
                end do ! iband
              end do ! idirC
            end do ! idirR

            ! Generate the small matching coefficients
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

            ! Generate the large matching coefficients by multiplication of z0
            ab0cofKet(:, :) = cmplx(0., 0.)
            do iband = 1, nobd(ikpt, 1)
              lmp = 0
              lm  = 0
              do oqn_l = 0, atoms%lmax(iDtype)
                do mqn_m = - oqn_l, oqn_l
                  pMaxLocal = nRadFun(oqn_l, iDtype)
                  ! p = 1
                  ab0cofKet(lmp + 1, iband) = iu**oqn_l * dot_product( conjg(z0(:nv(1, ikpt), iband, ikpt, 1)), a(:nv(1, ikpt), lm, iDatom) )
                  ! p = 2
                  ab0cofKet(lmp + 2, iband) = iu**oqn_l * dot_product( conjg(z0(:nv(1, ikpt), iband, ikpt, 1)), b(:nv(1, ikpt), lm, iDatom) )
                  ! Add LO contributions
                  do iradf = 3, pMaxLocal
                    ! p = 1
                    ab0cofKet(lmp + 1, iband) = ab0cofKet(lmp + 1, iband) + iu**oqn_l * &
                      & dot_product( conjg(z0(nv(1, ikpt) + 1:nv(1, ikpt) + atoms%nlotot, iband, ikpt, 1)), &
                                                           & bascof_lo(1, mqn_m, :atoms%nlotot, iloTable(iradf, oqn_l, iDtype), iDatom) )
                    ! p = 2
                    ab0cofKet(lmp + 2, iband) = ab0cofKet(lmp + 2, iband) + iu**oqn_l * &
                      & dot_product( conjg(z0(nv(1, ikpt) + 1:nv(1, ikpt) + atoms%nlotot, iband, ikpt, 1)), &
                                                           & bascof_lo(2, mqn_m, :atoms%nlotot, iloTable(iradf, oqn_l, iDtype), iDatom) )
                    ! 2 < p < LOs for that l and that atom type
                    ab0cofKet(lmp + iradf, iband) = iu**oqn_l * dot_product( conjg(z0(nv(1, ikpt) + 1:nv(1, ikpt) + atoms%nlotot, iband, ikpt, 1)), &
                                                           & bascof_lo(3, mqn_m, :atoms%nlotot, iloTable(iradf, oqn_l, iDtype), iDatom) )
                  end do ! iradf


                  lm = lm + 1
                  lmp = lmp + pMaxLocal
                end do ! mqn_m
              end do !oqn_l
            end do ! iband

            do idirC = 1, 3

                ! Complete the general surface integral multiplying the large matching coefficients
                overlap(:) = cmplx(0., 0.)
                surfIntMT(:) = cmplx(0., 0.)
                call CalcSintMT( ikpt, ikpt, idirC, iDatom, iDtype, lmpMax, nobd, hFullNoAbcofBK, overlapNoAbcofBK, ab0cofKet, ab0cofKet,   &
                                                                                                                & overlap, surfIntMT )
              do iband = 1, nobd(ikpt, 1)
                intMtSfNat(mqn_m2PrBra, idirC) = intMtSfNat(mqn_m2PrBra, idirC) + results%w_iks(iband, ikpt, 1) * (surfIntMT(iband))
                overlapSumNat(mqn_m2PrBra, idirC) = overlapSumNat(mqn_m2PrBra, idirC) + results%w_iks(iband, ikpt, 1) * (overlap(iband))
              end do ! iband
            end do ! idirC

          end do ! ikpt
        end do ! mqn_m2PrBra

        ! We have gradients of wave functions in natural coordinates to be transformed into cartesian coordinates
        intMtSf(1:3, 1:3) = matmul(conjg(Tmatrix(1:3, 1:3)), intMtSfNat(-1:1, 1:3))
        overlapSum(1:3, 1:3) = matmul(conjg(Tmatrix(1:3, 1:3)), overlapSumNat(-1:1, 1:3))

    if (.false.) then
      write(*, '(a)') 'IR Surface integral benchmark'
      write(*, '(3(2(es16.8,1x),3x))') 4 * surfIntSum(1, :) / 3
      write(*, '(3(2(es16.8,1x),3x))') 4 * surfIntSum(2, :) / 3
      write(*, '(3(2(es16.8,1x),3x))') 4 * surfIntSum(3, :) / 3

      write(*, '(a)') 'IR Surface integral benchmarko'
      write(*, '(3(2(es16.8,1x),3x))') 4 * surfIntSumO(1, :) / 3
      write(*, '(3(2(es16.8,1x),3x))') 4 * surfIntSumO(2, :) / 3
      write(*, '(3(2(es16.8,1x),3x))') 4 * surfIntSumO(3, :) / 3

      write(*, '(a)') 'MT Surface integral benchmark'
      write(*, '(3(2(es16.8,1x),3x))') 4 * intMtSf(1, :)
      write(*, '(3(2(es16.8,1x),3x))') 4 * intMtSf(2, :)
      write(*, '(3(2(es16.8,1x),3x))') 4 * intMtSf(3, :)

      write(*, '(a)') 'MT Surface integral benchmarko'
      write(*, '(3(2(es16.8,1x),3x))') 4 * overlapSum(1, :)
      write(*, '(3(2(es16.8,1x),3x))') 4 * overlapSum(2, :)
      write(*, '(3(2(es16.8,1x),3x))') 4 * overlapSum(3, :)
    end if

    do idirC = 1, 3
      do idirR = 1, 3
        if ( abs(surfIntIR2Arg(idirR, idirC) - 4 * surfIntSum(idirR, idirC) / 3 ) > 5e-7 ) passed = .false.
        if ( abs(surfIntMT2Arg(idirR, idirC) - 4 * intMtSf(idirR, idirC)) > 5e-5 ) passed = .false.
        if ( abs(4 * surfIntSum(idirR, idirC) / 3 + 4 * intMtSf(idirR, idirC)) > 5e-5) passed = .false.
        if ( abs(4 * surfIntSumO(idirR, idirC) / 3 + 4 * overlapSum(idirR, idirC)) > 5e-5) passed = .false.
      end do ! idirR
    end do ! idirC

    ! We repeat the whole test but now the gradient is in the ket
    allocate( vEff0GrVarphi(2, atoms%jmtd, 0:(atoms%lmaxd + 3)**2 - 1, lmpMax, -1:1) )
    vEff0GrVarphi = cmplx(0., 0.)
    do mqn_m2PrKet = -1, 1
      call calcFnsphGrVarphi( atoms, iDtype, 0, mqn_m2PrKet, nRadFun, grVarphiCh1, grVarphiCh2, grVarphiChLout, grVarphiChMout, vEff0MtSpH, &
                                                                                      & vEff0GrVarphi(:, :, :, :, mqn_m2PrKet) )
    end do ! mqn_m2PrKet

    chanMaxBra = 1
    chanMaxKet = 2
    lmaxBra = atoms%lmax(iDtype)
    muKet(:, :) = 0
    mqn_m2PrBra = -1
    varphiKet1 = 0.
    varphiKet2 = 0.
    do mqn_m = -atoms%lmax(iDtype), atoms%lmax(iDtype)
      muKet(mqn_m, -1) = mqn_m
    end do ! mqn_m
    do oqn_l = 0, atoms%lmax(iDtype)
      lambdaKet(1, oqn_l) = oqn_l
    end do ! oqn_l
    lmp = 0
    do oqn_l = 0, atoms%lmax(iDtype)
      do mqn_m = -oqn_l, oqn_l
        do iradf = 1, nRadFun(oqn_l, iDtype)
          lmp = lmp + 1
          !todo this should be varphiBra
          varphiKet1(atoms%jri(iDtype), 1, lmp, -1) = varphi1(atoms%jri(iDtype), iradf, oqn_l )
          varphiKet2(atoms%jri(iDtype), 1, lmp, -1) = varphi2(atoms%jri(iDtype), iradf, oqn_l )
        end do ! iradf
      end do ! mqn_m
    end do ! oqn_l

    overlapSumNat(:, :) = cmplx(0., 0.)
    intMtSfNat(:, :) = cmplx(0., 0.)
    surfIntSum(:, :) = cmplx(0., 0.)
    surfIntSumO(:, :) = cmplx(0., 0.)
    do mqn_m2PrKet = -1, 1
    !todo this is not optimal because interstitial integrals are executed threee times
      hFullNoAbcofBK(:, :, :, :) = cmplx(0., 0.)
      overlapNoAbcofBK(:, :, :, :) = cmplx(0., 0.)
      ! Calculate overlap varphi GradVarphi and varphi Veff0 gradVarphi
      call PrepareMTSurfIntDM( atoms, iDtype, chanMaxBra, chanMaxKet, lmaxBra, lmpMax, mqn_m2PrBra, mqn_m2PrKet, grVarphiChMout,  muKet,    &
        & grVarphiChLout, lambdaKet, nRadFun, grVarphiCh1, grVarphiCh2, varphiKet1, varphiKet2, vEff0GrVarphi(:, :, :, :, mqn_m2PrKet), atoms%rmt(iDtype)**2, hFullNoAbcofBK, overlapNoAbcofBK )
      ! Note for q = 0 the matrices should be diagonal in the end for the dynamical matrix

      overlap = cmplx(0., 0.)
      surfIntMT = cmplx(0., 0.)
      do ikpt = 1, kpts%nkpt

        ! As we do not store the grad z0
        surfIntVFast = cmplx(0., 0.)
        kExt(1:3) = matmul(cell%bmat, kpts%bk(1:3, ikpt))
        gExt(:) = 0.
        gradz0(:, :, :) = cmplx(0., 0.)
        gBasMapped(:, :) = 0.
        do iBas = 1, nv(1, ikpt) !+ atoms%nlotot
          gExt(1:3) = matmul( cell%bmat(1:3, 1:3), gBas(1:3, mapGbas(iBas, ikpt, 1)))
          do iband = 1, nobd(ikpt, 1)
            do idirR = 1, 3
              gradz0(iBas, iband, idirR) = iu * ( kExt(idirR) + gExt(idirR) ) * z0(iBas, iband, ikpt, 1)
            end do ! idirR
          end do ! iband
          gBasMapped(1:3, iBas) = real( gBas(1:3, mapGbas(iBas, ikpt, 1)) )
        end do ! iBas

       !todo directions correct?
        do idirR = 1, 3
          ! Calculate the IR integral of Psi veff0 grPsi with the Sternheimer routine
          call calcSurfVeff( atoms, kpts, cell, dimens, ikpt, ikpt, coScale, gbas, veffUvIR, nv, nobd(ikpt, 1), &
            & mapGbas, z0(:, :, ikpt, 1), gradz0(:, :, idirR), iDtype, iDatom, surfIntVFast )
          do idirC = 1, 3
            surfIntTeps = cmplx(0.0, 0.0)
            surfInt = cmplx(0.0, 0.0)
            ! Calculate overlap of Psi grad Psi with Sternheimer routine
            call calcSintKinEnergOvl( atoms, cell, kpts, results, iDtype, iDatom, ikpt, ikpt, idirC, &
              & nv(1, :), nobd(ikpt, 1), gbasMapped, gBasMapped, gBasMapped, gBasMapped, z0(:, :, ikpt, 1), gradz0(:, :, idirR), surfIntTeps, surfInt )
            do iband = 1, nobd(ikpt, 1)
              surfIntSum(idirR, idirC) = surfIntSum(idirR, idirC) + results%w_iks(iband, ikpt, 1) * surfIntVFast(iband, idirC)
              surfIntSumO(idirR, idirC) = surfIntSumO(idirR, idirC) + results%w_iks(iband, ikpt, 1) * surfInt(iband, iband)
            end do ! iband
          end do ! idirC
        end do ! idirR

        ! Small coefficients have to bei calculated again because we do not store k-dependent quantities
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

        ! Large matching coefficients are calculated multiplying the z0
        ab0cofKet(:, :) = cmplx(0., 0.)
        do iband = 1, nobd(ikpt, 1)
          lmp = 0
          lm  = 0
          do oqn_l = 0, atoms%lmax(iDtype)
            do mqn_m = - oqn_l, oqn_l
              pMaxLocal = nRadFun(oqn_l, iDtype)
              ! p = 1
              ab0cofKet(lmp + 1, iband) = iu**oqn_l * dot_product( conjg(z0(:nv(1, ikpt), iband, ikpt, 1)), a(:nv(1, ikpt), lm, iDatom) )
              ! p = 2
              ab0cofKet(lmp + 2, iband) = iu**oqn_l * dot_product( conjg(z0(:nv(1, ikpt), iband, ikpt, 1)), b(:nv(1, ikpt), lm, iDatom) )
              ! Add LO contributions
              do iradf = 3, pMaxLocal
                ! p = 1
                ab0cofKet(lmp + 1, iband) = ab0cofKet(lmp + 1, iband) + iu**oqn_l * &
                  & dot_product( conjg(z0(nv(1, ikpt) + 1:nv(1, ikpt) + atoms%nlotot, iband, ikpt, 1)), &
                                                       & bascof_lo(1, mqn_m, :atoms%nlotot, iloTable(iradf, oqn_l, iDtype), iDatom) )
                ! p = 2
                ab0cofKet(lmp + 2, iband) = ab0cofKet(lmp + 2, iband) + iu**oqn_l * &
                  & dot_product( conjg(z0(nv(1, ikpt) + 1:nv(1, ikpt) + atoms%nlotot, iband, ikpt, 1)), &
                                                       & bascof_lo(2, mqn_m, :atoms%nlotot, iloTable(iradf, oqn_l, iDtype), iDatom) )
                ! 2 < p < LOs for that l and that atom type
                ab0cofKet(lmp + iradf, iband) = iu**oqn_l * dot_product( conjg(z0(nv(1, ikpt) + 1:nv(1, ikpt) + atoms%nlotot, iband, ikpt, 1)), &
                                                       & bascof_lo(3, mqn_m, :atoms%nlotot, iloTable(iradf, oqn_l, iDtype), iDatom) )
              end do ! iradf

              lm = lm + 1
              lmp = lmp + pMaxLocal
            end do ! mqn_m
          end do !oqn_l
        end do ! iband

        do idirC = 1, 3

            overlap(:) = cmplx(0., 0.)
            surfIntMT(:) = cmplx(0., 0.)
            call CalcSintMT( ikpt, ikpt, idirC, iDatom, iDtype, lmpMax, nobd, hFullNoAbcofBK, overlapNoAbcofBK, ab0cofKet, ab0cofKet,   &
                                                                                                            & overlap, surfIntMT )
          do iband = 1, nobd(ikpt, 1)
            intMtSfNat(mqn_m2PrKet, idirC) = intMtSfNat(mqn_m2PrKet, idirC) + results%w_iks(iband, ikpt, 1) * (surfIntMT(iband))
            overlapSumNat(mqn_m2PrKet, idirC) = overlapSumNat(mqn_m2PrKet, idirC) + results%w_iks(iband, ikpt, 1) * (overlap(iband))
          end do ! iband
        end do ! idirC

      end do ! ikpt
    end do ! mqn_m2PrKet
    intMtSf(1:3, 1:3) = matmul(Tmatrix(1:3, 1:3), intMtSfNat(-1:1, 1:3))
    overlapSum(1:3, 1:3) = matmul(Tmatrix(1:3, 1:3), overlapSumNat(-1:1, 1:3))

    if (.false.) then
      write(*, '(a)') 'IR Surface integral benchmark'
      write(*, '(3(2(es16.8,1x),3x))') 4 * surfIntSum(1, :) / 3
      write(*, '(3(2(es16.8,1x),3x))') 4 * surfIntSum(2, :) / 3
      write(*, '(3(2(es16.8,1x),3x))') 4 * surfIntSum(3, :) / 3

      write(*, '(a)') 'IR Surface integral benchmarko'
      write(*, '(3(2(es16.8,1x),3x))') 4 * surfIntSumO(1, :) / 3
      write(*, '(3(2(es16.8,1x),3x))') 4 * surfIntSumO(2, :) / 3
      write(*, '(3(2(es16.8,1x),3x))') 4 * surfIntSumO(3, :) / 3

      write(*, '(a)') 'MT Surface integral benchmark'
      write(*, '(3(2(es16.8,1x),3x))') 4 * intMtSf(1, :)
      write(*, '(3(2(es16.8,1x),3x))') 4 * intMtSf(2, :)
      write(*, '(3(2(es16.8,1x),3x))') 4 * intMtSf(3, :)

      write(*, '(a)') 'MT Surface integral benchmarko'
      write(*, '(3(2(es16.8,1x),3x))') 4 * overlapSum(1, :)
      write(*, '(3(2(es16.8,1x),3x))') 4 * overlapSum(2, :)
      write(*, '(3(2(es16.8,1x),3x))') 4 * overlapSum(3, :)
    end if

    do idirC = 1, 3
      do idirR = 1, 3
        if ( abs(surfIntIR2Arg(idirR, idirC) - 4 * surfIntSum(idirR, idirC) / 3) > 5e-7 ) passed = .false.
        if ( abs(surfIntMT2Arg(idirR, idirC) - 4 * intMtSf(idirR, idirC)) > 5e-5 ) passed = .false.
        if ( abs(4 * surfIntSum(idirR, idirC) / 3 + 4 * intMtSf(idirR, idirC)) > 5e-5 ) passed = .false.
        if ( abs(4 * surfIntSumO(idirR, idirC) / 3 + 4 * overlapSum(idirR, idirC)) > 5e-5 ) passed = .false.
      end do ! idirR
    end do ! idirC

      end do ! iDeqat
    end do ! iDtype

    if ( passed ) then
      write( logUnit, '(a)' ) '                                                           |__ passed!'
      write( logUnit, * )
    else
      write( logUnit, '(a)' ) '                                                           |__ failed!'
      write( logUnit, * )
      call JuDFT_warn('Test of surface integrals containing wave function gradients', calledby='TestComp2ArgGrSFInts', &
        & hint='Activate and check debug output', &
        & file='m_jpTestDynMatSurf_mod.F90')
    end if

  end subroutine TestComp2ArgGrSFInts

end module m_jpTestDynMatSurf
