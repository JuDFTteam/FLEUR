module m_jpSetupDynMat
    USE m_constants

!NOTE: Coretail correction is implicetly in the dynamical matrix in the density and their gradients and variations of as well as
! in the potential
  implicit none

  contains

  subroutine SetupDynamicMatrix(fmpi, noco, nococonv, oneD, atoms, input, sym, cell, lathar, stars, kpts, qpts, usdus, results, Veff0, iqpt, ngdp, ngpqdp, gdp, mlh_atom, nmem_atom, clnu_atom, &
      & rho0IRst, rho1IR, rho1MT, vExt1MT, vEff1IR, vEff1MT, vEff0IR, vEff0MT, rho0MT, E2ndOrdII, El, eig, rbas1, rbas2, &
      & iloTable, nv, nobd, ilst, GbasVec, z, kveclo, nRadFun, mapKpq2K, kpq2kPrVec, gpqdp, memd_atom, logUnit, vXC0IR, eXCIR, vXC0MT, eXCMT, vExt1IR_final, vHar1IR_final, vHar1MT_final, grRho0IR, grRho0MT, grVext0IR, grVext0MT, grVeff0IR, grVeff0MT, dynMat, rho1MTDelta, vExt1MTDelta, vExt1MTq0, vHar1MTDelta, vHar1MTq0, vXc1MTDelta, vXc1MTq0, vExt1MTnoVol, grVeff0MThxc, vEff1MTnoVol, vH1MTnoVol, vExt1MTnoVolnoq, vExt1noqIR_final, rho1MTz0, grVCoul0IR_DM_SF, grVCoul0MT_DM_SF, vCoul1IRtempNoVol, vCoul1MTtempNoVol )

    use m_dfpt_init, only : convertStar2G
    use m_types
    use m_jpSetupDynMatSF, only : SetupDynMatSF
    use m_juDFT_stop, only : juDFT_error

    implicit none

    ! Type parameters
    type(t_mpi),                  intent(in)  :: fmpi
    type(t_nococonv),                intent(in)  :: nococonv
    type(t_noco),                intent(in)  :: noco
    type(t_oneD),                intent(in)  :: oneD
    type(t_atoms),                  intent(in)  :: atoms
    type(t_input),                  intent(in)  :: input
    type(t_sym),                    intent(in)  :: sym
    type(t_stars),                  intent(in)  :: stars
    type(t_cell),                   intent(in)  :: cell
    type(t_sphhar),                 intent(in)  :: lathar
    type(t_kpts),                   intent(in)  :: kpts
    type(t_kpts),                   intent(in)  :: qpts
    type(t_usdus),                  intent(in)  :: usdus
    type(t_results),                intent(in)  :: results
    type(t_potden),              intent(in)  :: Veff0

    ! Scalar parameters
    integer,                        intent(in)  :: iqpt
    integer,                        intent(in)  :: ngdp
    integer,                        intent(in)  :: ngpqdp
    integer,                       intent(in)  :: logUnit
    integer,                       intent(in)  :: memd_atom

    ! Array parameters
    integer,                        intent(in)  :: mlh_atom(:, 0:, :)
    integer,                        intent(in)  :: nmem_atom(0:, :)
    complex,                        intent(in)  :: clnu_atom(:, 0:, :)
    complex,                        intent(in)  :: rho0IRst(:, :)
    complex,                        intent(in)  :: rho1IR(:, :, :)
    complex,                        intent(in)  :: rho1MT(:, :, :, :, :)
    complex,                        intent(in)  :: rho1MTz0(:, :, :, :)
    complex,                        intent(in)  :: vExt1MT(:, :, :, :, :)
    complex,                        intent(in)  :: vExt1MTnoVol(:, :, :, :, :)
    complex,                        intent(in)  :: vEff1MTnoVol(:, :, :, :, :)
    complex,                        intent(in)  :: vExt1MTnoVolnoq(:, :, :, :, :)
    complex,                        intent(in)  :: grVCoul0IR_DM_SF(:, :)
    complex,                        intent(in)  :: grVCoul0MT_DM_SF(:, :, :, :)
    complex,                        intent(in)  :: vH1MTnoVol(:, :, :, :, :)
    complex,                        intent(in)  :: vEff1IR(:, :, :)
    complex,                        intent(in)  :: vEff1MT(:, :, :, :, :)
    real,                           intent(in)  :: rho0MT(:, :, :, :)
    complex,                        intent(in)  :: E2ndOrdII(:, :)
    integer,                        intent(in)  :: gdp(:, :)
    integer,                        intent(in)  :: gpqdp(:, :)
    real,                           intent(in)  :: El(:, 0:, :, :)
    real,                           intent(in)  :: eig(:, :, :)
    real,                           intent(in)  :: rbas1(:,:,0:,:,:)
    real,                           intent(in)  :: rbas2(:,:,0:,:,:)
    integer,                        intent(in)  :: iloTable(:, 0:, :)
    integer,                        intent(in)  :: nv(:, :)
    integer,                        intent(in)  :: nobd(:, :)
    integer,                        intent(in)  :: ilst(:, :, :)
    integer,                        intent(in)  :: GbasVec(:, :)
    complex,                       intent(in)  :: z(:, :, :, :)
    integer,                        intent(in)  :: kveclo(:,:)
    integer,                        intent(in)  :: nRadFun(:, :)
    complex,                        intent(in)  :: vEff0IR(:,:)
    real,                           intent(in)  :: vEff0MT(:, 0:, :)
    integer,                        intent(in)  :: mapKpq2K(:, :)
    integer,                        intent(in)  :: kpq2kPrVec(:, :, :)
    complex,                        intent(in)  :: vXC0IR(:, :)
    complex,                        intent(in)  :: eXCIR(:)
    real,                           intent(in)  :: vXC0MT(:, 0:, :, :)
    real,                           intent(in)  :: eXCMT(:, 0:, :)
    complex,                        intent(in)  :: vExt1IR_final(:, :, :)
    complex,                        intent(in)  :: vExt1noqIR_final(:, :, :)
    complex,                        intent(in)  :: vHar1IR_final(:, :, :)
    complex,                        intent(in)  :: vHar1MT_final(:, :, :, :, :)
    complex,                        intent(in)  :: grVext0IR(:, :)
    complex,                        intent(in)  :: grVext0MT(:, :, :, :)
    complex,                        intent(in)  :: grVeff0IR(:, :)
    complex,                        intent(in)  :: grVeff0MT(:, :, :, :)
    complex,                        intent(in)  :: grVeff0MThxc(:, :, :, :)
    complex,                        intent(in)  :: grRho0IR(:, :)
    complex,                        intent(in)  :: grRho0MT(:, :, :, :)
    complex,                        intent(in)  :: rho1MTDelta(:, :, :, :, :)
    complex,                        intent(in)  :: vExt1MTDelta(:, :, :, :, :)
    complex,                        intent(in)  :: vExt1MTq0(:, :, :, :, :)
    complex,                        intent(in)  :: vHar1MTDelta(:, :, :, :, :)
    complex,                        intent(in)  :: vHar1MTq0(:, :, :, :, :)
    complex,                        intent(in)  :: vXc1MTDelta(:, :, :, :, :)
    complex,                        intent(in)  :: vXc1MTq0(:, :, :, :, :)
    complex,                        intent(in)  :: vCoul1IRtempNoVol(:, :)
    complex,                        intent(in)  :: vCoul1MTtempNoVol(:, :, :, :)
    complex,           allocatable, intent(out) :: dynMat(:, :)

    ! Scalar variables
    integer                                     :: idirR
    integer                                     :: idirC
    integer                                     :: iDtype
    integer                                     :: iDatom
    integer                                     :: iDeqat
    integer                                     :: idir
    integer                                     :: iatom
    integer                                     :: itype
    integer                                     :: ieqat
    integer                                     :: oqn_l
    integer                                     :: mqn_m
    integer                                     :: lm
    integer                                     :: imesh


    ! Array variables
    complex, allocatable                        :: dynMatHF(:, :)
    complex, allocatable                        :: dynMatPu(:, :)
    complex, allocatable                        :: dynMatSf(:, :)
    complex, allocatable                        :: rho0IRpw(:)


    if (.false.) then
      ! Attention! 2 lines of vext1Delta and vhar1Delta must be uncommented
      iDatom = 0
      do iDtype = 1, atoms%ntype
        do iDeqat = 1, atoms%neq(iDtype)
          iDatom = iDatom + 1
          do idir = 1, 3
            iatom = 0
            do itype = 1, atoms%ntype
              do ieqat = 1, atoms%neq(itype)
                iatom = iatom + 1
                do oqn_l = 0, atoms%lmax(itype)
                  do mqn_m = -oqn_l, oqn_l
                    lm = oqn_l * (oqn_l + 1) + 1 + mqn_m
                    do imesh = 1, atoms%jri(itype)
                      write(7501, '(5i8,2es15.8)') iatom, idir, oqn_l, mqn_m, imesh, rho1MTDelta(imesh, lm, iatom, idir, iDatom)
                      write(7502, '(5i8,2es15.8)') iatom, idir, oqn_l, mqn_m, imesh, rho1MT(imesh, lm, iatom, idir, iDatom)
                      write(7503, '(5i8,2es15.8)') iatom, idir, oqn_l, mqn_m, imesh, grRho0MT(imesh, lm, iDatom, idir)
                      write(7504, '(5i8,2es15.8)') iatom, idir, oqn_l, mqn_m, imesh, vExt1MTDelta(imesh, lm, iatom, idir, iDatom)
                      write(7505, '(5i8,2es15.8)') iatom, idir, oqn_l, mqn_m, imesh, vExt1MTq0(imesh, lm, iatom, idir, iDatom)
                      write(7506, '(5i8,2es15.8)') iatom, idir, oqn_l, mqn_m, imesh, vHar1MTDelta(imesh, lm, iatom, idir, iDatom)
                      write(7507, '(5i8,2es15.8)') iatom, idir, oqn_l, mqn_m, imesh, vHar1MTq0(imesh, lm, iatom, idir, iDatom)
                      write(7508, '(5i8,2es15.8)') iatom, idir, oqn_l, mqn_m, imesh, vXc1MTDelta(imesh, lm, iatom, idir, iDatom)
                      write(7509, '(5i8,2es15.8)') iatom, idir, oqn_l, mqn_m, imesh, vXc1MTq0(imesh, lm, iatom, idir, iDatom)
                    end do ! imesh
                  end do ! mqn_m
                end do ! oqn-l
              end do ! ieqat
            end do ! itype
          end do ! idir
        end do ! iDeqat
      end do ! iDtype
      call juDFT_error('Old juPhon stopcall.', calledby='SetupDynamicMatrix')
    end if

    allocate (dynMat(3 * atoms%nat, 3 * atoms%nat))

    dynMat = cmplx(0.0, 0.0)

    allocate( rho0IRpw(ngdp) )
    rho0IRpw(:) = cmplx(0., 0.)

    ! Convert star representation of unperturbed density to plane-wave representation
    call convertStar2G( rho0IRst(:, 1), rho0IRpw, stars, ngdp, gdp )

    ! 23.09.21 - Discussion notes:
    ! - Check the interstitial surface integrals on whether it's always \vec{e} * V^{(1)T} as it should be. Also: correct sign?
    ! - HF 176,2,MT: zeroing is not correct.
    ! - HF Surface IR: ensure grad only of V_{ext,\alpha}; full Vext is not correct.
    ! - HF Surface MT: zeroing is not correct.
    ! - v1Deltas/=v1goods ---> Calculate and pass! ---> Put integrals after Pulay integration of rho1*VH

    ! Calculate the Hellmann-Feynman contribution to the dynamical matrix
    call SetupDynMatHF(atoms, sym, cell, lathar, stars, ngdp, ngpqdp, gdp, mlh_atom, nmem_atom, clnu_atom, rho0IRpw, rho0MT, grRho0IR, grRho0MT,  &
      & rho1IR, rho1MT, grVext0IR, grVext0MT, vExt1IR_final, vExt1MT, E2ndOrdII, dynMatHF, rho1MTDelta, vExt1MTDelta, vExt1MTq0, iqpt, vExt1MTnoVol, vExt1MTnoVolnoq, vExt1noqIR_final)

    ! Calculate the Pulay contribution to the dynamical matrix
    call SetupDynMatPu( fmpi, noco, nococonv, oneD, atoms, stars, lathar, input, sym, kpts, qpts, cell, usdus, results, iqpt, ngdp, ngpqdp, gdp, mapKpq2K, rho1IR, rho1MT, &
      & vEff1IR, vEff1MT, grRho0IR, grRho0MT, grVeff0IR, grVeff0MT, El, eig, rbas1, rbas2, iloTable, nv, nobd, ilst, GbasVec, z, kveclo, nRadFun, clnu_atom, nmem_atom,    &
      & mlh_atom, vEff0IR, vEff0MT, kpq2kPrVec, dynMatPu, rho1MTDelta, vExt1MTDelta, vExt1MTq0, vHar1MTDelta, vHar1MTq0, vXc1MTDelta, vXc1MTq0, grVeff0MThxc, vEff1MTnoVol, &
      & vExt1MTnoVol, vH1MTnoVol, rho1MTz0 )

    ! Calculate the surface contribution to the dynamical matrix
    !!!latest:
    call SetupDynMatSF( fmpi, noco, nococonv, oneD, atoms, input, stars, cell, results, Veff0, kpts, qpts, lathar, sym, usdus, ngdp, iqpt, logUnit, &
      & memd_atom, nobd, gdp, mapKpq2K, rbas1(:, :, :, :, 1), rbas2(:, :, :, :, 1), nmem_atom, mlh_atom, clnu_atom, kveclo, iloTable, kpq2kPrVec, nv, ilst, &
      & gBasVec, nRadFun, z, eig, El, rho0IRpw, rho0MT, ngpqdp, gpqdp, rho1IR, rho1MTDelta, vXC0IR, eXCIR, vXC0MT, eXCMT, vExt1IR_final, &
      & vExt1MT, vHar1IR_final, vHar1MT_final, grRho0IR, grRho0MT, grVeff0IR, grVeff0MT, vEff0MT, grVCoul0IR_DM_SF, grVCoul0MT_DM_SF, vCoul1IRtempNoVol, vCoul1MTtempNoVol, dynMatSf )

    ! Add up all three contributions
    dynMat(:, :) = &
      &   dynMatHF(:, :) &
      & + dynMatPu(:, :) &
      & + dynMatSf(:, :)

    !if (compPhon) then
      !open(114,file='000_dynmat_grob',form='FORMATTED',position='append',action='WRITE',status='UNKNOWN')
        !write(114,*) 'HF:'
        !write(114,*) dynMatHF
        !write(114,*) 'Pulay:'
        !write(114,*) dynMatPu
        !write(114,*) 'Surface:'
        !write(114,*) dynMatSf
        !write(114,*) 'Full:'
        !write(114,*) dynMat
      !close(114)
    !end if

    do idirC=1,3
      do idirR=1,3
        if (idirR.ne.idirC) then
          write(565,*) 'HF Offdiag full'
          write(565,*) idirR, idirC, dynMatHF(idirR,idirC)
        end if
      end do
    end do

    !if (compPhon) then
    !dynMat(:, :) = dynMatHF(:, :)
    !end if

  end subroutine SetupDynamicMatrix

  ! We orientated at int_nv from fleurnewatbroyden version
  ! We have to use variables here which are not q-dependent. The q-dependency is outsourced.
  !Note: this is deprecated but inteteresting for the thesis
!                      ! Add delta distribution on diagonal
!                      ! The loop runs through 3 times for each diagonal element therefore we have to divide by 9 instead of 3
                       ! the division by rmsh is necessary if the routine gets r^2 rho!
                       ! Due to the upper text, we can do a linear interpolation here. We choose the value of rho at the first
                       ! mesh point. Still, we do not know, whether this is a good interpolation. Therefore, we can either proceed
                       ! on the linear appendix of rho heading the zero or we just renormalize the MT to the surface integral to
                       ! fulfill the Goldstone condition. Going further on the mesh up to zero is a bit problematic because actually
                       ! we do also not know the potential and its variation in mesh points smaller than the first mesh point.
                       ! Therefore, it is okay to renormalize here to fulfill the Goldstone condition. What is interesting is the
                       ! difference to zero not how much zero is. Still we have to keep that in mind when evaluating other integrals
                       ! how accurate we can evaluate them.
                       ! If we move a bit from the first mesh point the result is really sensitive, the question is whether this
                       ! makes sense or to move more than the mesh point to zero or just take the first mesh point and renormalize
                       ! to the surface integral in the MT.
  subroutine SetupDynMatHF(atoms, sym, cell, lathar, stars, ngdp, ngpqdp, gdp, mlh_atom, nmem_atom, clnu_atom, rho0IRpw, rho0MT, grRho0IR, grRho0MT,  &
      & rho1IR, rho1MT, grVext0IR, grVext0MT, vExt1IR, vExt1MT, E2ndOrdII, dynMatHF, rho1MTDelta, vExt1MTDelta, vExt1MTq0, iqpt, vExt1MTnoVol, vExt1MTnoVolnoq, vExt1IRnoq )

    use m_types, only : t_atoms, t_sym, t_cell, t_sphhar, t_stars
    use m_juDFT_stop, only : juDFT_error

    implicit none

    ! Type parameters
    type(t_atoms),               intent(in)  :: atoms
    type(t_sym),               intent(in)  :: sym
    type(t_cell),                intent(in)  :: cell
    type(t_sphhar),              intent(in)  :: lathar
    type(t_stars),               intent(in)  :: stars

    ! Scalar parameters
    integer,                     intent(in)  :: ngdp
    integer,                     intent(in)  :: ngpqdp
    integer,                        intent(in)  :: iqpt

    ! Array parameters
    integer,                     intent(in)  :: gdp(:, :)
    integer,                     intent(in)  :: mlh_atom(:, 0:, :)
    integer,                     intent(in)  :: nmem_atom(0:, :)
    complex,                     intent(in)  :: clnu_atom(:, 0:, :)
    complex,                     intent(in)  :: rho0IRpw(:)
    real,                        intent(in)  :: rho0MT(:, 0:, :, :)
    complex,                     intent(in)  :: grRho0IR(:, :)
    complex,                     intent(in)  :: grRho0MT(:, :, :, :)
    complex,                     intent(in)  :: rho1IR(:, :, :)
    complex,                     intent(in)  :: rho1MT(:, :, :, :, :)
    !todo warping of vExt1IR should be done here
    complex,                     intent(in)  :: grVext0IR(:, :)
    complex,                     intent(in)  :: grVext0MT(:, :, :, :)
    complex,                     intent(in)  :: vExt1IR(:, :, :)
    complex,                     intent(in)  :: vExt1IRnoq(:, :, :)
    complex,                     intent(in)  :: vExt1MT(:, :, :, :, :)
    complex,                     intent(in)  :: vExt1MTnoVol(:, :, :, :, :)
    complex,                     intent(in)  :: vExt1MTnoVolnoq(:, :, :, :, :)
    complex,                     intent(in)  :: E2ndOrdII(:, :)
    complex,                     intent(in)  :: rho1MTDelta(:, :, :, :, :)
    complex,                        intent(in)  :: vExt1MTDelta(:, :, :, :, :)
    complex,                        intent(in)  :: vExt1MTq0(:, :, :, :, :)
    complex,        allocatable, intent(out) :: dynMatHF(:, :)

    ! Scalar variables
    integer                                  :: iAatom
    integer                                  :: iBatom
    integer                                  :: iCatom
    integer                                  :: iAtype
    integer                                  :: iBtype
    integer                                  :: iCtype
    integer                                  :: iAeqat
    integer                                  :: iBeqat
    integer                                  :: iCeqat
    integer                                  :: idirC
    integer                                  :: idirR
    integer                                  :: oqn_l
    integer                                  :: mqn_m
    integer                                  :: lm
    integer                                  :: lm_pre, imesh
    complex                                  :: integral, integralsum
    logical                                  :: finiteQoptimization
    logical                                  :: testoptimization
    integer                                  :: idir
    integer                                  :: iG

    ! Array variables
    complex,        allocatable              :: w_grVext0IR(:, :), rho1MTdummy(:,:,:,:,:), rho1IRdummy(:,:,:)
    complex,        allocatable              :: w_vExt1IR(:, :, :), w_vExt1IRnoq(:, :, :)
    complex,        allocatable              :: vExt1MTContainer(:, :, :, :)
    complex                                  :: integral3x3(3, 3)
    complex                                  :: dynMatHFTest(3, 3)
    complex                                  :: dynMatCompInt1(3, 3)
    complex                                  :: dynMatCompInt1MgradInt(3, 3)

    ! For finite q, the integral rho1 Vext1 in the MT can be split up so that the part of it corresponding to gradRho gradVext1
    ! can be canceled analytically.
    ! For q = = nothing is changed
    finiteQoptimization = .false.
    !if (compPhon) then
      finiteQoptimization = .false.
    !end if
    testoptimization = .false.

    dynMatHFTest(:, :) = cmplx(0., 0.)
    dynMatCompInt1(:, :) = cmplx(0., 0.)
    dynMatCompInt1MgradInt(:, :) = cmplx(0., 0.)

    allocate( dynMatHF( 3 * atoms%nat, 3 * atoms%nat) )
    allocate( w_grVext0IR(ngdp, 3) )
    dynMatHF(:, :) = cmplx(0., 0.)
    w_grVext0IR(:, :) = cmplx(0., 0.)

    !allocate(rho1MTdummy, mold=rho1MT)
    !allocate(rho1IRdummy, mold=rho1IR)
    !rho1MTdummy=cmplx(1.0,0.0)
    !rho1IRdummy=cmplx(0.0,0.0)

    !do iG=1, ngdp
    !  if ((gdp(1, iG).eq.0).and.(gdp(2, iG).eq.0).and.(gdp(3, iG).eq.0)) then
    !    rho1IRdummy(iG, :, 1)=cmplx(1.0,0.0)
    !  end if
    !end do

    allocate( vExt1MTContainer(atoms%jmtd, (atoms%lmaxd + 2)**2, 3, atoms%nat) )

    !todo beware to use shifted G-set here
    allocate(w_vExt1IR(ngdp, 3, atoms%nat))
    allocate(w_vExt1IRnoq(ngdp, 3, atoms%nat))
    w_vExt1IR(:, :, :) = cmplx(0., 0.)
    w_vExt1IRnoq(:, :, :) = cmplx(0., 0.)
    iAatom = 0
    do iAtype = 1, atoms%ntype
      do iAeqat = 1, atoms%neq(iAtype)
        iAatom = iAatom + 1
        do idir = 1, 3
      !      call warpIRPot(stars, ngpqdp, idir, gpqdp, vExt1IRtemp, w_vExt1IRtemp(:, idir))
          call warpIRPot(stars, ngdp, idir, gdp, vExt1IR(:, :, iAatom), w_vExt1IR(:, idir, iAatom))
          call warpIRPot(stars, ngdp, idir, gdp, vExt1IRnoq(:, :, iAatom), w_vExt1IRnoq(:, idir, iAatom))
        end do ! idir
      end do ! iAeqat
    end do ! iAtype


    ! Warp interstitial gradient of the external unperturbed potential
    do idirC = 1, 3
      call warpIRPot( stars, ngdp, idirC, gdp, grVext0IR(:, :), w_grVext0IR(:, idirC) )
    end do

    !if (compPhon) then
      !open(109,file='000_dynmat',form='FORMATTED',position='append',action='WRITE',status='UNKNOWN')
      !write(109,*) 'iqph', iqpt
      !write(109,*) 'Term 1 IR, then MT'
    !end if
    iAatom = 0
    do iAtype = 1, atoms%ntype
      do iAeqat = 1, atoms%neq(iAtype)
        iAatom = iAatom + 1
        do idirC = 1, 3
          iBatom = 0
          do iBtype = 1, atoms%ntype
            do iBeqat = 1, atoms%neq(iBtype)
              iBatom = iBatom + 1
              do idirR = 1, 3

                ! Interstitial part of 7.112b PhDAK
                ! Calculate IR integral of rho1 and vExt1IR (the latter is warped with the step function)
                ! Altogether with IR up to 1e-7
                ! For the IR quantities, we take all G-vectors up to Gmax as there might be contributions for G-vectors between
                ! 2kmax and Gmax whenever coretails are expanded in the IR
                ! Note: Not using the shifted G-set for the variated quantities leads to a non-hermitian and complex dynamical matrix
                integral = cmplx(0., 0.)
                call Calc2ArgIntIR( cell, ngpqdp, rho1IR(:ngpqdp, idirR, iBatom), w_vExt1IR(:ngpqdp, idirC, iAatom), integral )
                !call Calc2ArgIntIR( cell, ngpqdp, rho1IRdummy(:ngpqdp, idirR, iBatom), w_vExt1IR(:ngpqdp, idirC, iAatom), integral )

                !(5.3.176), 1st integral, IR
                dynMatHF(idirR + (iBatom - 1) * 3, idirC + (iAatom - 1) * 3) = &
                  & dynMatHF(idirR + (iBatom - 1) * 3, idirC + (iAatom - 1) * 3) + integral
                write(465,*) 'HF rho1 V1 IR'
                write(465,*) idirR, idirC, integral
                if ((iBatom.eq.iAatom).and.(idirR.eq.idirC)) then
                  write(565,*) 'HF rho1 V1 IR'
                  write(565,*) idirR, idirC, integral
                end if
                !if (compPhon) then
                  !write(109,*) idirR, idirC, iBatom, iAatom
                  !write(109,*) 'IR1', real(integral), aimag(integral)
                !end if

                !!!Note: 000_dynmat here

                if ((iqpt == 1) .or. (.not.finiteQoptimization)) then
                  if ((idirR == 1) .and. (idirC == 1)) then
                    write(*, *) "Original method rho1 Vext1 1/2"
                  end if
                  integralsum = cmplx(0.0,0.0)
                  iCatom = 0
                  do iCtype = 1, atoms%ntype
                    do iCeqat = 1, atoms%neq(iCtype)
                      iCatom = iCatom + 1
                      ! TODO rho1MT up to l + 1, at vext1MT up to l, what is correct?
                      do oqn_l = 0, atoms%lmax(iCtype)
                        do mqn_m = -oqn_l, oqn_l
                          lm = oqn_l * (oqn_l + 1) + 1 + mqn_m

                          ! Muffin-tin part of 7.112b PhDAK
                          ! This integral is significantly dependent on the interpolation starting from the first mesh point to zero.
                          ! This is due to the fact that the integrand itself is about 1e6 multiplied by about 1e-5 of the mesh so the
                          ! product although still to be integrated is really significant and a powerful setting to influence the
                          ! result of this integral. It can be even influenced such that the trace in Vext2 and Eii2 can be subtracted
                          ! again on the diagonal and still the Goldstone condition can be fulfilled. As the integrand returns to zero
                          ! again and would be really zero at r = 0, we decided for Gustavs idea for a triangle interpolation.
                          ! This seems to be a simple but rather effective interpolation. The question of having an interpolation
                          ! or not influences this integral on the order of 1e1, so it is quiet important.
                          ! Moreover, we decided to evaluate the integral up to the sixth mesh point and then to interpolate. This is
                          ! due to the fact that the numerical gradient of the density is due to some boundary effects quiet turbulent
                          ! at the core on a logarithmic x-axis. Still, we do not know the curve of the integrand exactly and probably
                          ! have no good interpolation. Therefore, we are not able to show the mathematical relation derived by Fabian
                          ! in an exact way.
                          ! The product of rho1 and vext1 for q = 0 the l=1 component is 1e16
                          integral = cmplx(0., 0.)
                          !if (compPhon) then
                        !    call Calc2ArgCmplxIntMT( atoms, iCtype, rho1MT(:, lm, iCatom, idirR, iBatom), &
                              !& vExt1MT(:, lm, iCatom, idirC, iAatom), integral )
                          !call Calc2ArgCmplxIntMT( atoms, iCtype, rho1MTdummy(:, lm, iCatom, idirR, iBatom), &
                          !  & vExt1MT(:, lm, iCatom, idirC, iAatom), integral )
                         ! else
                            !call Calc2ArgCmplxIntMT( atoms, iCtype, rho1MT(:, lm, iCatom, idirR, iBatom), &
                            !    & vExt1MTnoVol(:, lm, iCatom, idirC, iAatom), integral )
                            call Calc2ArgCmplxIntMT( atoms, iCtype, rho1MT(:, lm, iCatom, idirR, iBatom), &
                                & vExt1MT(:, lm, iCatom, idirC, iAatom), integral )
                            !!!CRGfix
                            !call Calc2ArgCmplxIntMT( atoms, iCtype, rho1MT(:, lm, iCatom, idirR, iBatom), &
                            !    & vExt1MTnoVol(:, lm, iCatom, idirC, iAatom), integral )
                          !end if
                          !write(466,*) rho1MTDelta(:, lm, iCatom, idirR, iBatom)
                          !write(469,*) vExt1MTnoVol(:, lm, iCatom, idirC, iAatom), integral
                          !(5.3.176), 1st integral, MT
                          dynMatHF(idirR + (iBatom - 1) * 3, idirC + (iAatom - 1) * 3) = &
                            &dynMatHF(idirR + (iBatom - 1) * 3, idirC + (iAatom - 1) * 3) + integral

                          !if (compPhon) then
                            integralsum=integralsum+integral
                            !write(109,*) 'MT1', real(integral), aimag(integral)
                          !end if

                          !!!Note: 000_dynmat here

                        end do ! mqn_m
                      end do ! oqn_l
                    end do ! iCeqat
                  end do ! iCtype
                  !if (compPhon) then
                  write(465,*) 'HF rho1 V1 MT'
                  write(465,*) idirR, idirC, integralsum
                  if ((iBatom.eq.iAatom).and.(idirR.eq.idirC)) then
                    write(565,*) 'HF rho1 V1 MT'
                    write(565,*) idirR, idirC, integralsum
                  end if
                  !write(109,*) 'MT1', real(integralsum), aimag(integralsum)
                  !end if


                else if ((iqpt /= 1) .and. finiteQoptimization) then
                  if ((idirR == 1) .and. (idirC == 1)) then
                    write(*, *) "Optimized method rho1 Vext1"
                  end if
                  integralsum=cmplx(0.0,0.0)
                  dynMatHFTest(:, :) = cmplx(0., 0.)
                  iCatom = 0
                  do iCtype = 1, atoms%ntype
                    do iCeqat = 1, atoms%neq(iCtype)
                      iCatom = iCatom + 1
                      ! TODO rho1MT up to l + 1, at vext1MT up to l, what is correct?
                      do oqn_l = 0, atoms%lmax(iCtype)
                        do mqn_m = -oqn_l, oqn_l
                          lm = oqn_l * (oqn_l + 1) + 1 + mqn_m


                          ! q-dependant part of rho, q-dependent part of Vext
                          integral = cmplx(0., 0.)
                          call Calc2ArgCmplxIntMT( atoms, iCtype, rho1MTDelta(:, lm, iCatom, idirR, iBatom), &
                            & vExt1MTDelta(:, lm, iCatom, idirC, iAatom), integral )

                          dynMatHFTest(idirR + (iBatom - 1) * 3, idirC + (iAatom - 1) * 3) = &
                            &dynMatHFTest(idirR + (iBatom - 1) * 3, idirC + (iAatom - 1) * 3) + integral

                          !!!Note: 000_dynmat here

                          if (testoptimization) then
                            dynMatCompInt1(idirR + (iBatom - 1) * 3, idirC + (iAatom - 1) * 3) = &
                              &dynMatCompInt1(idirR + (iBatom - 1) * 3, idirC + (iAatom - 1) * 3) + integral
                          end if

                          !!!Note: 000_dynmat here

                          ! q=0 part of rho, q-dependant part of Vext
                          integral = cmplx(0., 0.)
                          call Calc2ArgCmplxIntMT( atoms, iCtype, -grRho0MT(:, lm, iCatom, idirR), &
                            & vExt1MTDelta(:, lm, iCatom, idirC, iAatom), integral )

                          dynMatHFTest(idirR + (iBatom - 1) * 3, idirC + (iAatom - 1) * 3) = &
                            &dynMatHFTest(idirR + (iBatom - 1) * 3, idirC + (iAatom - 1) * 3) + integral

                          !!!Note: 000_dynmat here

                          if (testoptimization) then
                            dynMatCompInt1(idirR + (iBatom - 1) * 3, idirC + (iAatom - 1) * 3) = &
                              &dynMatCompInt1(idirR + (iBatom - 1) * 3, idirC + (iAatom - 1) * 3) + integral
                          end if

                          !!!Note: 000_dynmat here

                          ! q-dependent part of rho, q = 0 part of Vext
                          integral = cmplx(0., 0.)
                          call Calc2ArgCmplxIntMT( atoms, iCtype, rho1MTDelta(:, lm, iCatom, idirR, iBatom), &
                            & vExt1MTq0(:, lm, iCatom, idirC, iAatom), integral )

                          dynMatHFTest(idirR + (iBatom - 1) * 3, idirC + (iAatom - 1) * 3) = &
                            &dynMatHFTest(idirR + (iBatom - 1) * 3, idirC + (iAatom - 1) * 3) + integral

                          if (testoptimization) then
                            dynMatCompInt1(idirR + (iBatom - 1) * 3, idirC + (iAatom - 1) * 3) = &
                              &dynMatCompInt1(idirR + (iBatom - 1) * 3, idirC + (iAatom - 1) * 3) + integral
                          end if


                          if (testoptimization) then
                            ! Make substitution complete with q=0 parts of rho1 and Vext1
                            integral = cmplx(0., 0.)
                            call Calc2ArgCmplxIntMT( atoms, iCtype, -grRho0MT(:, lm, iCatom, idirR), &
                              & vExt1MTq0(:, lm, iCatom, idirC, iAatom), integral )

                            dynMatCompInt1(idirR + (iBatom - 1) * 3, idirC + (iAatom - 1) * 3) = &
                              &dynMatCompInt1(idirR + (iBatom - 1) * 3, idirC + (iAatom - 1) * 3) + integral

                            ! Subtract integral to compare
                            integral = cmplx(0., 0.)
                            call Calc2ArgCmplxIntMT( atoms, iCtype, rho1MT(:, lm, iCatom, idirR, iBatom), &
                              & vExt1MT(:, lm, iCatom, idirC, iAatom), integral )

                            dynMatCompInt1(idirR + (iBatom - 1) * 3, idirC + (iAatom - 1) * 3) = &
                              &dynMatCompInt1(idirR + (iBatom - 1) * 3, idirC + (iAatom - 1) * 3) - integral

                          end if ! test Optimization

                        end do ! mqn_m
                      end do ! oqn_l
                    end do ! iCeqat
                  end do ! iCtype
                  !write(109,*) 'MT', real(dynMatHFTest(idirR + (iBatom - 1) * 3, idirC + (iAatom - 1) * 3)), aimag(dynMatHFTest(idirR + (iBatom - 1) * 3, idirC + (iAatom - 1) * 3))
                  dynMatHF(idirR + (iBatom - 1) * 3, idirC + (iAatom - 1) * 3) = &
                    &dynMatHF(idirR + (iBatom - 1) * 3, idirC + (iAatom - 1) * 3) + dynMatHFTest(idirR + (iBatom - 1) * 3, idirC + (iAatom - 1) * 3)

                  if (testoptimization) then
                    if (iAatom == iBatom) then

                      dynMatCompInt1MgradInt(idirR + (iBatom - 1) * 3, idirC + (iAatom - 1) * 3) = &
                                                                    & dynMatHFTest(idirR + (iBatom - 1) * 3, idirC + (iAatom - 1) * 3)

                      iCatom = 0
                      do iCtype = 1, atoms%ntype
                        do iCeqat = 1, atoms%neq(iCtype)
                          iCatom = iCatom + 1

                          do oqn_l = 0, atoms%lmax(iCatom)
                            do mqn_m = -oqn_l, oqn_l
                              lm = oqn_l * ( oqn_l + 1 ) + 1 + mqn_m

                              integral = cmplx(0., 0.)
                              call Calc2ArgCmplxIntMT( atoms, iCtype, rho1MT(:, lm, iCatom, idirR, iBatom), &
                                & vExt1MT(:, lm, iCatom, idirC, iAatom), integral )

                              dynMatCompInt1MgradInt(idirR + (iBatom - 1) * 3, idirC + (iAatom - 1) * 3) = &
                                &dynMatCompInt1MgradInt(idirR + (iBatom - 1) * 3, idirC + (iAatom - 1) * 3) - integral

                              ! MT vanishes up to 9e-8
                              ! MT volume integral of grRho and grVext0
                              integral = cmplx(0., 0.)
                              call Calc2ArgCmplxIntMT( atoms, iCtype, grRho0MT(:, lm, iCatom, idirR), grVext0MT(:, lm,  idirC, iCatom), integral )

                              dynMatCompInt1MgradInt(idirR + (iBatom - 1) * 3, idirC + (iAatom - 1) * 3) = &
                                & dynMatCompInt1MgradInt(idirR + (iBatom - 1) * 3, idirC + (iAatom - 1) * 3) + integral

                            end do ! oqn_l
                          end do ! mqn_m

                        end do ! iCeqat
                      end do ! iCtype

                    end if  ! iAtom == iBatom
                  end if ! testoptimization
                end if ! Optimization enabled for finite q

                if ( iAatom == iBatom ) then

                  ! It should only be the gradient of Vext from MT iAatom; is that given? Same for the latter MT
                  ! IR volume integral of grRho0IR and w_grVext0IR
                  integral = cmplx(0., 0.)
                  call Calc2ArgIntIR( cell, ngdp, grRho0IR(:ngdp, idirR), w_vExt1IRnoq(:ngdp, idirC, iAatom), integral )
                  !call Calc2ArgIntIR( cell, ngdp, rho1IRdummy(:ngdp, idirR, iBatom), w_grVext0IR(:ngdp, idirC), integral )

                  !(5.3.176), 2nd integral, IR
                  dynMatHF(idirR + (iBatom - 1) * 3, idirC + (iAatom - 1) * 3) = &
                     & dynMatHF(idirR + (iBatom - 1) * 3, idirC + (iAatom - 1) * 3) + integral
                  write(465,*) 'HF grRho V1_q0 IR'
                  write(465,*) idirR, idirC, integral
                  if ((iBatom.eq.iAatom).and.(idirR.eq.idirC)) then
                    write(565,*) 'HF grRho V1_q0 IR'
                    write(565,*) idirR, idirC, integral
                  end if
                  !if (compPhon) then
                    !write(109,*) 'IR2', real(-integral), aimag(-integral)
                  !end if

                  if ((iqpt == 1) .or. (.not.finiteQoptimization)) then

                    if ((idirR==1) .and. (idirC==1)) then
                      write(*, *) "Original method rho1 Vext1 2/2"
                    end if
                    integralsum=cmplx(0.0,0.0)
                    iCatom = 0
                    do iCtype = 1, atoms%ntype
                      do iCeqat = 1, atoms%neq(iCtype)
                        iCatom = iCatom + 1

                        do oqn_l = 0, atoms%lmax(iCatom)
                          do mqn_m = -oqn_l, oqn_l
                            lm = oqn_l * ( oqn_l + 1 ) + 1 + mqn_m

                            ! MT vanishes up to 9e-8
                            ! MT volume integral of grRho and grVext0
                            integral = cmplx(0., 0.)
                            !if (compPhon) then
                             ! call Calc2ArgCmplxIntMT( atoms, iCtype, grRho0MT(:, lm, iCatom, idirR), grVext0MT(:, lm,  idirC, iCatom), integral )
                              !call Calc2ArgCmplxIntMT( atoms, iCtype, rho1MTdummy(:, lm, iCatom, idirR, iBatom), grVext0MT(:, lm,  idirC, iCatom), integral )
                            !dynMatHF(idirR + (iBatom - 1) * 3, idirC + (iAatom - 1) * 3) = &
                            !  & dynMatHF(idirR + (iBatom - 1) * 3, idirC + (iAatom - 1) * 3) - integral
                            !else
                              !call Calc2ArgCmplxIntMT( atoms, iCtype, grRho0MT(:, lm, iCatom, idirR), vExt1MTnoVolnoq(:, lm, iCatom, idirC, iAatom), integral )
                              call Calc2ArgCmplxIntMT( atoms, iCtype, grRho0MT(:, lm, iCatom, idirR), grVext0MT(:, lm,  idirC, iCatom), integral )
                              !!!CRGfix
                              !!call Calc2ArgCmplxIntMT( atoms, iCtype, grRho0MT(:, lm, iCatom, idirR), grVext0MT(:, lm,  idirC, iCatom), integral )
                            !dynMatHF(idirR + (iBatom - 1) * 3, idirC + (iAatom - 1) * 3) = &
                            !  & dynMatHF(idirR + (iBatom - 1) * 3, idirC + (iAatom - 1) * 3) - integral
                            !write(467,*) grRho0MT(:, lm, iCatom, idirR), vExt1MTnoVolnoq(:, lm, iCatom, idirC, iAatom), integral
                            !end if

                            !(5.3.176), 2nd integral, MT
                            integralsum=integralsum-integral

                          end do ! oqn_l
                        end do ! mqn_m
                      end do ! iCeqat
                    end do ! iCtype
                    write(465,*) 'HF grRho V1_q0 MT'
                    write(465,*) idirR, idirC, integralsum
                    if ((iBatom.eq.iAatom).and.(idirR.eq.idirC)) then
                      write(565,*) 'HF grRho V1_q0 MT'
                      write(565,*) idirR, idirC, integralsum
                    end if

                    !if (compPhon) then
                    !  write(109,*) 'MT2', real(integralsum), aimag(integralsum)
                    !end if
                  end if ! if q = 0 or for finite q if no optimization activated

                end if ! iAatom == iBatom

                ! Evaluation of 7.112d
                !Note: This minus is not part of the theory. However, it should be there, because we have charge neutratility, coming
                ! in the system, so the integrals rho vext (e- Z) and rho Vhar/Vxc (e- e-) and Eii(Z Z) should compensate eacht other
                ! With the minus in here and the trace activated, the Eii2 compensates up to 1e-4 with the other parts of the HF
                ! contribution. The charge neutrality should also hold after a spatial derivative of rho0 Vext yielding rho1 Vext1
                ! rho0 Vext2 and Eii2. We leave the minus now, although we cannot explain analytically for the moment where it comes
                ! from.

                !if (compPhon) then
                  !if (idirC.eq.idirR.and.iAatom.eq.iBatom) then
                    !(5.3.177), 2nd term
                    dynMatHF(idirR + (iBatom - 1) * 3, idirC + (iAatom - 1) * 3) =                                                    &
                      & dynMatHF(idirR + (iBatom - 1) * 3, idirC + (iAatom - 1) * 3) +                                                &
                      & E2ndOrdII(idirR + (iBatom - 1) * 3, idirC + (iAatom - 1) * 3)
                    !write(109,*) '3', real(E2ndOrdII(idirR + (iBatom - 1) * 3, idirC + (iAatom - 1) * 3)), aimag(E2ndOrdII(idirR + (iBatom - 1) * 3, idirC + (iAatom - 1) * 3))
                    write(465,*) 'HF Eii2'
                    write(465,*) idirR, idirC, E2ndOrdII(idirR + (iBatom - 1) * 3, idirC + (iAatom - 1) * 3)
                    if ((iBatom.eq.iAatom).and.(idirR.eq.idirC)) then
                      write(565,*) 'HF Eii2'
                      write(565,*) idirR, idirC, E2ndOrdII(idirR + (iBatom - 1) * 3, idirC + (iAatom - 1) * 3)
                    end if
                  !end if
                !else
                !  dynMatHF(idirR + (iBatom - 1) * 3, idirC + (iAatom - 1) * 3) =                                                    &
                !    & dynMatHF(idirR + (iBatom - 1) * 3, idirC + (iAatom - 1) * 3) +                                                &
                !    & E2ndOrdII(idirR + (iBatom - 1) * 3, idirC + (iAatom - 1) * 3)
                !end if

              end do ! idirR
            end do ! iBeqat
          end do ! iBtype
        end do ! idirC

        do idirC = 1, 3
          iCatom = 0
          do iCtype = 1, atoms%ntype
            do iCeqat = 1, atoms%neq(iCtype)
              iCatom = iCatom + 1
              do oqn_l = 0, atoms%lmax(iCtype)
                lm_pre = oqn_l * (oqn_l + 1) + 1
                do mqn_m = -oqn_l, oqn_l
                  lm = lm_pre + mqn_m
                  do imesh = 1, atoms%jri(iCtype)
                    vExt1MTContainer(imesh, lm, idirC, iCatom) = -grVext0MT(imesh, lm,  idirC, iCatom)!vExt1MTnoVolnoq(imesh, lm, iCatom, idirC, iAatom)
                    !write(467,*) vExt1MTContainer(imesh, lm, idirC, iCatom)
                  end do ! imesh
                end do ! mqn_m
              end do ! oqn_l
            end do ! iDeqat
          end do ! iDtypeB
        end do ! idirC

        integral3x3(:, :) = cmplx(0., 0.)

        call CalcSurfIntIRDynMat( atoms, cell, ngdp, ngdp, gdp, gdp, rho0IRpw(:), vExt1IRnoq(:, :, iAatom), [0., 0., 0.], integral3x3 )

        do idirC = 1, 3
          do idirR = 1, 3
          !(5.3.177), integral, IR
          dynMatHF(idirR + (iAatom - 1) * 3, idirC + (iAatom - 1) * 3) = &
            & dynMatHF(idirR + (iAatom - 1) * 3, idirC + (iAatom - 1) * 3) - integral3x3(idirR, idirC)
          write(465,*) 'HF SF IR'
          write(465,*) idirR, idirC, -integral3x3(idirR, idirC)
          end do ! idirR
        end do ! idirC

        iCatom = 0
        do iCtype = 1, atoms%ntype
          do iCeqat = 1, atoms%neq(iCtype)
            iCatom = iCatom + 1

            integral3x3 = cmplx(0., 0.)

            call CalcSurfIntMTDynMat( atoms, sym, lathar, clnu_atom, nmem_atom, mlh_atom, rho0MT(:, :, :, 1), vExt1MTContainer, integral3x3)
            do idirC = 1, 3
              do idirR = 1, 3
                !(5.3.177), integral, MT
                dynMatHF(idirR + (iAatom - 1) * 3, idirC + (iAatom - 1) * 3) = &
                  & dynMatHF(idirR + (iAatom - 1) * 3, idirC + (iAatom -1) * 3) - integral3x3(idirR, idirC)
                write(465,*) 'HF SF MT'
                write(465,*) idirR, idirC, -integral3x3(idirR, idirC)
              end do ! idirR
            end do ! idirC
          end do ! iCeqat
        end do ! iCtype


!        ! IR surface integral of rho0 and grVext0
!        ! TODO grRho up to l + 1, at grVext up to l, what is correct?
!        integral3x3(:, :) = cmplx(0., 0.)
!        ! Within this routine grVext is transposed
!        ! Surface integrals vanish up to 2e-9 in sum
!        call CalcSurfIntIRDynMat( atoms, cell, ngdp, ngdp, gdp, gdp, rho0IRpw(:), grVext0IR, [0., 0., 0.], integral3x3 )
!
!        do idirC = 1, 3
!          do idirR = 1, 3
!            if (.not.compPhon) then
!            !(5.3.177), integral, IR
!            dynMatHF(idirR + (iAatom - 1) * 3, idirC + (iAatom - 1) * 3) = &
!              & dynMatHF(idirR + (iAatom - 1) * 3, idirC + (iAatom - 1) * 3) + integral3x3(idirR, idirC)
!            end if
!          end do ! idirR
!        end do ! idirC
!
!        ! MT surface integral of rho0 and grVext
!        iCatom = 0
!        do iCtype = 1, atoms%ntype
!          do iCeqat = 1, atoms%neq(iCtype)
!            iCatom = iCatom + 1
!
!            integral3x3 = cmplx(0., 0.)
!            call CalcSurfIntMTDynMat( atoms, lathar, clnu_atom, nmem_atom, mlh_atom, rho0MT(:, :, :, 1), grVext0MT, integral3x3)
!
!            do idirC = 1, 3
!              do idirR = 1, 3
!                if (.not.compPhon) then
!                  !(5.3.177), integral, MT
!                  dynMatHF(idirR + (iAatom - 1) * 3, idirR + (iAatom - 1) * 3) = &
!                    & dynMatHF(idirR + (iAatom - 1) * 3, idirR + (iAatom -1) * 3) + integral3x3(idirR, idirC)
!                end if
!              end do ! idirR
!            end do ! idirC
!
!          end do ! iCeqat
!        end do ! iCtype

      end do ! iAeqat
    end do ! iAtype

    !if (compPhon) then
      !close(109)
    !end if

    ! Check of comparison between optimization and original method
    if (any(abs(dynMatCompInt1(:, :)) > 1e-8)) call juDFT_error('rho1 Vext1 integral not the same as optimization', calledby='SetupDynMatHF')
    if (any(abs(dynMatCompInt1MgradInt(:, :)) > 1e-8)) call juDFT_error('rho1 Vext1 integral -gradRho Veff1 integral not the same as optimization', calledby='SetupDynMatHF')

    if (.false.) then
      write(*, '(a)') 'Test Matrix'
      write(*, '(3(2(es16.8,1x),3x))') dynMatCompInt1(1, :)
      write(*, '(3(2(es16.8,1x),3x))') dynMatCompInt1(2, :)
      write(*, '(3(2(es16.8,1x),3x))') dynMatCompInt1(3, :)

      write(*, '(a)') 'Test Matrix'
      write(*, '(3(2(es16.8,1x),3x))') dynMatCompInt1MgradInt(1, :)
      write(*, '(3(2(es16.8,1x),3x))') dynMatCompInt1MgradInt(2, :)
      write(*, '(3(2(es16.8,1x),3x))') dynMatCompInt1MgradInt(3, :)
    end if

  end subroutine SetupDynMatHF


  !>--------------------------------------------------------------------------------------------------------------------------------
  !> @author
  !> Christian-Roman Gerhorst
  !>
  !> @brief
  !> Operates the evaluation of 6.52 PhD thesis A. Klueppelberg, i.e. the Pulay contributions to the setup of the dynamical matrix.
  !>
  !> @details
  !> Equation 6.52 is subdivided in 7.114, 7.115 and 7.118 (all equations from PhD thesis A. Klueppelberg). For every of the latter
  !> equations a respective subroutine is called.
  !>--------------------------------------------------------------------------------------------------------------------------------
  subroutine SetupDynMatPu(fmpi, noco, nococonv, oneD, atoms, stars, lathar, input, sym, kpts, qpts, cell, usdus, results, iqpt, ngdp, ngpqdp, gdp, mapKpq2K, rho1IR, rho1MT,      &
      & vEff1IR, vEff1MT, grRho0IR, grRho0MT, grVeff0IR, grVeff0MT, El, eig, rbas1, rbas2, iloTable, nv, nobd, ilst, GbasVec, z, kveclo, nRadFun, clnu_atom, nmem_atom,    &
      & mlh_atom, vEff0IR, vEff0MT, kpq2kPrVec, dynMatPu, rho1MTDelta, vExt1MTDelta, vExt1MTq0, vHar1MTDelta, vHar1MTq0, vXc1MTDelta, vXc1MTq0, grVeff0MThxc, vEff1MTnoVol, &
      & vExt1MTnoVol, vH1MTnoVol, rho1MTz0 )

    use m_types

    implicit none

    ! Type parameters
    type(t_mpi),                  intent(in)  :: fmpi
    type(t_nococonv),                intent(in)  :: nococonv
    type(t_noco),                intent(in)  :: noco
    type(t_oneD),                intent(in)  :: oneD
    type(t_atoms),     intent(in)  :: atoms
    type(t_stars),     intent(in)  :: stars
    type(t_sphhar),    intent(in)  :: lathar
    type(t_input),     intent(in)  :: input
    type(t_sym),       intent(in)  :: sym
    type(t_kpts),      intent(in)  :: kpts
    type(t_kpts),      intent(in)  :: qpts
    type(t_cell),      intent(in)  :: cell
    type(t_usdus),     intent(in)  :: usdus
    type(t_results),   intent(in)  :: results

    integer,           intent(in)  :: iqpt
    integer,           intent(in)  :: ngdp
    integer,           intent(in)  :: ngpqdp

    ! rho1IR     : interstitial first variation of density
    ! rho1MT     : muffin-tin first variation of density
    ! vEff1IR  : warped first variation of effective potential in interstitial
    ! vEff1MT    : first variation of effective potential in muffin-tin
    ! El         : LAPW energy parameters
    ! eig        : Kohn-Sham eigen energies
    ! rbas1      : large component of radial basis functions
    ! rbas2      : small component of radial basis functions
    ! iloTable   : index array to reorganize LOs
    ! nv         : number of all eigenvalues
    ! nobd       : number of occupied eigenvalues
    ! ilst       : index array to avoid redundant G-basis vectors
    ! GbasVec    : G-basis vector
    ! z          : Kohn-Sham wavefunction expansion coefficients
    ! kveclo     : index array for selecting correct G-vector for respective LO
    ! nRadFun    : number of radial functions for certain l and k-point
    ! clnu_atom  : clnu given that the non-representative MTs are not rotated
    ! nmem_atom  : nmem given that the non-representative MTs are not rotated
    ! mlh_atom   : mlh given that the non-representative MTs are not rotated
    ! vEff0IR    : unperturbed effective potential in interstitial
    ! vEff0MT    : unperturbed effective potential in muffin-tin
    ! mapKpq2K   : index of k' = k + q, i.e. which k-point set index has the (maybe backfolded) sum of k and q
    ! kpq2kPrVec : contains backfolding vector, if k' = k + q required a backfolding in the unit cell.
    integer,           intent(in)  :: gdp(:, :)
    complex,           intent(in)  :: rho1IR(:, :, :)
    complex,           intent(in)  :: rho1MT(:, :, :, :, :)
    complex,           intent(in)  :: rho1MTz0(:, :, :, :)
    complex,           intent(in)  :: vEff1IR(:, :, :)
    complex,           intent(in)  :: vEff1MT(:, :, :, :, :)
    complex,           intent(in)  :: vEff1MTnoVol(:, :, :, :, :)
    complex,           intent(in)  :: vExt1MTnoVol(:, :, :, :, :)
    complex,           intent(in)  :: vH1MTnoVol(:, :, :, :, :)
    complex,                    intent(in)  :: grRho0IR(:, :)
    complex,                    intent(in)  :: grRho0MT(:, :, :, :)
    complex,                    intent(in)  :: grVeff0IR(:, :)
    complex,                    intent(in)  :: grVeff0MT(:, :, :, :)
    complex,                    intent(in)  :: grVeff0MThxc(:, :, :, :)
    real,              intent(in)  :: El(:, 0:, :, :)
    real,              intent(in)  :: eig(:, :, :)
    real,              intent(in)  :: rbas1(:,:,0:,:,:)
    real,              intent(in)  :: rbas2(:,:,0:,:,:)
    integer,           intent(in)  :: iloTable(:, 0:, :)
    integer,           intent(in)  :: nv(:, :)
    integer,           intent(in)  :: nobd(:, :)
    integer,           intent(in)  :: ilst(:, :, :)
    integer,           intent(in)  :: GbasVec(:, :)
    complex,          intent(in)  :: z(:, :, :, :)
    integer,           intent(in)  :: kveclo(:,:)
    integer,           intent(in)  :: nRadFun(:, :)
    complex,           intent(in)  :: clnu_atom(:, 0:, :)
    integer,           intent(in)  :: nmem_atom(0:, :)
    integer,           intent(in)  :: mlh_atom(:, 0:, :)
    complex,           intent(in)  :: vEff0IR(:,:)
    real,              intent(in)  :: vEff0MT(:, 0:, :)
    integer,           intent(in)  :: mapKpq2K(:, :)
    integer,           intent(in)  :: kpq2kPrVec(:, :, :)
    complex,           intent(in)  :: vExt1MTDelta(:, :, :, :, :)
    complex,           intent(in)  :: vExt1MTq0(:, :, :, :, :)
    complex,           intent(in)  :: vHar1MTDelta(:, :, :, :, :)
    complex,           intent(in)  :: vHar1MTq0(:, :, :, :, :)
    complex,           intent(in)  :: vXc1MTDelta(:, :, :, :, :)
    complex,           intent(in)  :: vXc1MTq0(:, :, :, :, :)
    complex,           intent(in)  :: rho1MTDelta(:, :, :, :, :)
    complex, allocatable,         intent(out) :: dynMatPu(:, :)

    complex, allocatable           :: dynMatPuInt(:, :)
    complex, allocatable           :: dynMatPuME(:, :)
    integer :: idirC, idirR

    allocate( dynMatPu( 3 * atoms%nat, 3 * atoms%nat) )
    dynMatPu = cmplx(0., 0.)

    ! Evaluate 7.114 PhD thesis A. Klueppelberg with recasted core contributions and contributions from recasted Pulay matrix elements.
    call EvalIntRho1Veff1(atoms, cell, stars, ngpqdp, gdp, rho1IR, rho1MT, vEff1IR, vEff1MT, grRho0MT, grVeff0MT, dynMatPuInt, rho1MTDelta, vExt1MTDelta, vExt1MTq0, vHar1MTDelta, vHar1MTq0, vXc1MTDelta, vXc1MTq0, iqpt, ngdp, grVeff0MThxc, vEff1MTnoVol, vExt1MTnoVol, vH1MTnoVol, rho1MTz0 )

    !write(470,*) dynMatPuInt

    ! Evaluate valence contribution of H - eps brakets
    call EvalPuHepsBraKetsVal(fmpi, noco, nococonv, oneD, atoms, input, cell, sym, lathar, stars, kpts, qpts, usdus, results, iqpt, nRadFun, vEff0IR, vEff0MT, clnu_atom,&
    & nmem_atom, mlh_atom, rbas1(:, :, :, :, 1), rbas2(:, :, :, :, 1), El, mapKpq2K, nobd, nv, gBasVec, ilst, kveclo, z, iloTable, &
    & eig, kpq2kPrVec, dynMatPuME )

    do idirC = 1, 3
      do idirR = 1, 3
        write(465,*) 'Pu Brakets'
        write(465,*) idirR, idirC, dynMatPuME(idirR, idirC)
        write(565,*) 'Pu Brakets'
        write(565,*) idirR, idirC, dynMatPuME(idirR, idirC)
      end do
    end do


    ! Add up contribution from Pulay volume integral and the matrix elements containing 1st and 2nd order variation of the WFs.
    dynMatPu(:, :) =  &
     &   dynMatPuInt  &
     & + dynMatPuME


  end subroutine SetupDynMatPu

  !>-------------------------------------------------------------------------------------------------------------------------------------
  !> @author
  !> Christian-Roman Gerhorst
  !>
  !> @brief
  !> Evaluates the integral in 7.114 PhD thesis A. Klueppelberg (PhDAK)
  !>
  !> @details
  !> The integral is split into an interstitial part and a muffin-tin part and evaluated according to 7.113 PhDAK. We can do this
  !> because we have only mixing terms with vanishing Bloch characters and due to time inversion symmetry we can express +q- by
  !> -q-quantities
  !>
  !> @note
  !> Core terms are included as the integrated quantities are already full-electron quantities. This routine does not work anymore
  !> in this form if time inversion symmetry is broken!
  !>
  !> @todo
  !> Account for spin
  !>
  !>-------------------------------------------------------------------------------------------------------------------------------------
  subroutine EvalIntRho1Veff1( atoms, cell, stars, ngpqdp, gdp, rho1IR, rho1MT, vEff1IR, vEff1MT, grRho0MT, grVeff0MT, dynMatPuInt, rho1MTDelta, vExt1MTDelta, vExt1MTq0, vHar1MTDelta, vHar1MTq0, vXc1MTDelta, vXc1MTq0, iqpt, ngdp, grVeff0MThxc, vEff1MTnoVol, vExt1MTnoVol, vH1MTnoVol, rho1MTz0 )

    use m_types

    implicit none

    ! Type parameters
    type(t_atoms),              intent(in)  :: atoms
    type(t_cell),               intent(in)  :: cell
    type(t_stars),              intent(in)  :: stars


    ! Scalar parameters
    integer,                    intent(in)  :: ngpqdp
    integer,                    intent(in)  :: ngdp
    integer,                    intent(in)  :: iqpt
    ! rho1IR     : interstitial first variation of density
    ! rho1MT     : muffin-tin first variation of density
    ! w_vEff1IR  : warped first variation of effective potential in interstitial
    ! vEff1MT    : first variation of effective potential in muffin-tin
    ! Array parameter
    integer,                    intent(in)  :: gdp(:, :)
    complex,                    intent(in)  :: rho1IR(:, :, :)
    complex,                    intent(in)  :: rho1MT(:, :, :, :, :)
    complex,                    intent(in)  :: rho1MTz0(:, :, :, :)
    complex,                    intent(in)  :: vEff1IR(:, :, :)
    complex,                    intent(in)  :: vEff1MT(:, :, :, :, :)
    complex,                    intent(in)  :: vEff1MTnoVol(:, :, :, :, :)
    complex,                    intent(in)  :: vExt1MTnoVol(:, :, :, :, :)
    complex,                    intent(in)  :: vH1MTnoVol(:, :, :, :, :)
    complex,                    intent(in)  :: grRho0MT(:, :, :, :)
    complex,                    intent(in)  :: grVeff0MT(:, :, :, :)
    complex,                    intent(in)  :: grVeff0MThxc(:, :, :, :)
    complex,                    intent(in)  :: vExt1MTDelta(:, :, :, :, :)
    complex,                    intent(in)  :: vExt1MTq0(:, :, :, :, :)
    complex,                    intent(in)  :: vHar1MTDelta(:, :, :, :, :)
    complex,                    intent(in)  :: vHar1MTq0(:, :, :, :, :)
    complex,                    intent(in)  :: vXc1MTDelta(:, :, :, :, :)
    complex,                    intent(in)  :: vXc1MTq0(:, :, :, :, :)
    complex,                    intent(in)  :: rho1MTDelta(:, :, :, :, :)
    complex,     allocatable,   intent(out) :: dynMatPuInt(:, :)

    ! variables beginning with i are loop variables. A stands for alpha, B for beta, G for gamma, D for displaced.
    ! Local variables
    integer                                 :: iAtype
    integer                                 :: iAeqat
    integer                                 :: idir
    integer                                 :: iAatom
    integer                                 :: iDtypeA
    integer                                 :: iDatomA
    integer                                 :: idirC
    integer                                 :: iBatom
    integer                                 :: iDtypeB
    integer                                 :: iDeqatB
    integer                                 :: idirR
    integer                                 :: iGatom
    integer                                 :: itypeG
    integer                                 :: ieqatG
    integer                                 :: oqn_l
    integer                                 :: mqn_m
    integer                                 :: lm
    complex                                 :: integral, integralsum, integralsum2
    logical                                 :: finiteQoptimization
    logical                                 :: testoptimization
    complex                                 :: dynMatHFTest(3, 3)
    complex                                 :: dynMatCompInt1(3, 3)
    complex                                 :: dynMatCompInt1MgradInt(3, 3)
    integer :: imesh

    complex,                   allocatable  :: w_vEff1IR(:, :, :)

    ! For finite q, the integral rho1 Veff1 in the MT can be split up so that the part of it corresponding to gradRho gradVeff1
    ! can be canceled analytically.
    ! Also, Veff is split up into Vext, Vhar and Vxc
    ! For q = = nothing is changed
    finiteQoptimization = .false.
    testoptimization = .false.

    allocate( dynMatPuInt( 3 * atoms%nat, 3 * atoms%nat ) )
    dynMatPuInt = cmplx(0., 0.)

    dynMatHFTest(:, :) = cmplx(0., 0.)
    dynMatCompInt1(:, :) = cmplx(0., 0.)
    dynMatCompInt1MgradInt(:, :) = cmplx(0., 0.)


    !todo beware to use shifted G-set here
    allocate(w_vEff1IR(ngpqdp, 3, atoms%nat))
    w_vEff1IR(:, :, :) = cmplx(0., 0.)
    iAatom = 0
    do iAtype = 1, atoms%ntype
      do iAeqat = 1, atoms%neq(iAtype)
        iAatom = iAatom + 1
        do idir = 1, 3
      !      call warpIRPot(stars, ngpqdp, idir, gpqdp,  vEff1IR(:, :, iAatom), w_vEff1IR(:, idir, iAatom))
          call warpIRPot(stars, ngdp, idir, gdp, vEff1IR(:, :, iAatom), w_vEff1IR(:, idir, iAatom))
        end do ! idir
      end do ! iAeqat
    end do ! iAtype

    iAatom = 0
    do iDtypeA = 1, atoms%ntype
      do iDatomA = 1, atoms%neq(iDtypeA)
        iAatom = iAatom + 1
        do idirC = 1, 3
          iBatom = 0
          do iDtypeB = 1, atoms%ntype
            do iDeqatB = 1, atoms%neq(iDtypeB)
              iBatom = iBatom + 1
              do idirR = 1, 3

                ! Evaluate IR part of rho1 w_vEff1IR. Note that the first argument of the dot_product is complex conjugated. The outer product is
                ! performed by choosing the correct vector component idirR or idirC of the first-order quantities.
                integral = cmplx(0., 0.)
                call Calc2ArgIntIR( cell, ngpqdp, rho1IR(:ngpqdp, idirR, iBatom), w_vEff1IR(:ngpqdp, idirC, iAatom), integral)

                !(5.3.178), IR
                !dynMatPuInt(idirR + (iBatom - 1) * 3, idirC + (iAatom - 1) * 3) = &
                !  & dynMatPuInt(idirR + (iBatom - 1) * 3, idirC + (iAatom - 1) * 3) + integral
                !write(465,*) '(5.3.178), IR'
                !write(465,*) idirR, idirC, integral
                if ((iqpt == 1) .or. (.not.finiteQoptimization)) then
                  if ((idirR == 1) .and. (idirC == 1)) then
                    write(*, *) "Original method rho1 Veff1 1/2"
                  end if

                  integralsum=cmplx(0.0,0.0)
                  iGatom = 0
                  do itypeG = 1, atoms%ntype
                    do ieqatG = 1, atoms%neq(itypeG)
                      iGatom = iGatom + 1
                      do oqn_l = 0, atoms%lmax(itypeG)
                        do mqn_m = -oqn_l, oqn_l
                          lm = oqn_l * (oqn_l + 1) + 1 + mqn_m

                          ! Evaluate the MT part of rho1 vEff1IR
                          ! Muffin-tin part of 7.112b PhDAK
                          ! This integral is significantly dependent on the interpolation starting from the first mesh point to zero.
                          ! This is due to the fact that the integrand itself is about 1e6 multiplied by about 1e-5 of the mesh so the
                          ! product although still to be integrated is really significant and a powerful setting to influence the
                          ! result of this integral. It can be even influenced such that the trace in Vext2 and Eii2 can be subtracted
                          ! again on the diagonal and still the Goldstone condition can be fulfilled. As the integrand returns to zero
                          ! again and would be really zero at r = 0, we decided for Gustavs idea for a triangle interpolation.
                          ! This seems to be a simple but rather effective interpolation. The question of having an interpolation
                          ! or not influences this integral on the order of 1e1, so it is quiet important.
                          ! Moreover, we decided to evaluate the integral up to the sixth mesh point and then to interpolate. This is
                          ! due to the fact that the numerical gradient of the density is due to some boundary effects quiet turbulent
                          ! at the core on a logarithmic x-axis. Still, we do not know the curve of the integrand exactly and probably
                          ! have no good interpolation. Therefore, we are not able to show the mathematical relation derived by Fabian
                          ! in an exact way.
                          ! The product of rho1 and vext1 for q = 0 the l component is 1e16
                      ! TODO rho1MT up to l + 1, at vext1MT up to l, what is correct?
                       !   integral = cmplx(0., 0.)
                       !   call Calc2ArgCmplxIntMT( atoms, itypeG, rho1MT(:, lm, iGatom, idirR, iBatom), vEff1MT(:, lm, iGatom, idirC, iAatom), integral)

                       !   dynMatPuInt(idirR + (iBatom - 1) * 3, idirC + (iAatom - 1) * 3) = &
                       !     & dynMatPuInt(idirR + (iBatom - 1) * 3, idirC + (iAatom - 1) * 3) + integral

                          integral = cmplx(0., 0.)
                          if (iGatom.eq.iBatom) then
                            call Calc2ArgCmplxIntMT( atoms, itypeG, rho1MTz0(:, lm, idirR, iBatom), vEff1MTnoVol(:, lm, iGatom, idirC, iAatom), integral)
                          end if

                          !(5.3.178), MT
                          dynMatPuInt(idirR + (iBatom - 1) * 3, idirC + (iAatom - 1) * 3) = &
                            & dynMatPuInt(idirR + (iBatom - 1) * 3, idirC + (iAatom - 1) * 3) + integral
                          integralsum=integralsum+integral
                        end do ! mqn_m
                      end do ! oqn_l

                    end do ! ieqatG
                  end do ! itypeG
                  write(465,*) 'Pu Int rho1z0 V1good'
                  write(465,*) idirR, idirC, integralsum
                  write(565,*) 'Pu Int rho1z0 V1good'
                  write(565,*) idirR, idirC, integralsum

                  integralsum2=cmplx(0.0,0.0)
                  do oqn_l = 0, atoms%lmax(iDtypeA)
                    do mqn_m = -oqn_l, oqn_l
                      lm = oqn_l * (oqn_l + 1) + 1 + mqn_m
                      integral = cmplx(0., 0.)
                      call Calc2ArgCmplxIntMT( atoms, iDtypeA, rho1MTDelta(:, lm, iAatom, idirR, iBatom), grVeff0MT(:, lm, idirC, iAatom), integral)

                      !(5.3.185)
                      dynMatPuInt(idirR + (iBatom - 1) * 3, idirC + (iAatom - 1) * 3) = &
                        & dynMatPuInt(idirR + (iBatom - 1) * 3, idirC + (iAatom - 1) * 3) + integral
                      integralsum2=integralsum2+integral
                    end do ! mqn_m
                  end do ! oqn_l

                  integralsum=cmplx(0.0,0.0)
                  !integralsum2=cmplx(0.0,0.0)
                  do oqn_l = 0, atoms%lmax(iDtypeB)
                    do mqn_m = -oqn_l, oqn_l
                      lm = oqn_l * (oqn_l + 1) + 1 + mqn_m
                      integral = cmplx(0., 0.)
                      call Calc2ArgCmplxIntMT( atoms, iDtypeB, grRho0MT(:, lm, idirR, iBatom), vEff1MTnoVol(:, lm, iBatom, idirC, iAatom), integral)

                  !    !(5.3.186) [2 components]
                  !    ! 1/2
                  !!!latest: The + is a hotfix!!!
                      dynMatPuInt(idirR + (iBatom - 1) * 3, idirC + (iAatom - 1) * 3) = &
                        & dynMatPuInt(idirR + (iBatom - 1) * 3, idirC + (iAatom - 1) * 3) - integral
                  !
                      integralsum=integralsum-integral
                      integral = cmplx(0., 0.)
                      !call Calc2ArgCmplxIntMT( atoms, iDtypeB, grVeff0MT(:, lm, idirR, iAatom), rho1MTDelta(:, lm, iBatom, idirC, iAatom), integral)
                      !!!CRGfix
                      !call Calc2ArgCmplxIntMT( atoms, iDtypeB, grVeff0MThxc(:, lm, idirR, iAatom), rho1MTDelta(:, lm, iBatom, idirC, iAatom), integral)
                  !
                  !    integral = cmplx(0., 0.)
                  !    call Calc2ArgCmplxIntMT( atoms, iDtypeB, grRho0MT(:, lm, idirR, iBatom), vH1MTnoVol(:, lm, iBatom, idirC, iAatom), integral)
                  !
                  !    ! 2/2
                      !dynMatPuInt(idirR + (iBatom - 1) * 3, idirC + (iAatom - 1) * 3) = &
                      !  & dynMatPuInt(idirR + (iBatom - 1) * 3, idirC + (iAatom - 1) * 3) + integral
                      !integralsum2=integralsum2+integral
                    end do ! mqn_m
                  end do ! oqn_l
                  write(465,*) 'Pu Int grRho V1good'
                  write(465,*) idirR, idirC, integralsum
                  write(465,*) 'Pu Int rho1good grV'
                  write(465,*) idirR, idirC, integralsum2
                  write(565,*) 'Pu Int grRho V1good'
                  write(565,*) idirR, idirC, integralsum
                  write(565,*) 'Pu Int rho1good grV'
                  write(565,*) idirR, idirC, integralsum2

                else if ((iqpt /= 1) .and. finiteQoptimization) then
                  if ((idirR == 1) .and. (idirC == 1)) then
                    write(*, *) "Optimized method rho1 Veff1"
                  end if

                  iGatom = 0
                  do itypeG = 1, atoms%ntype
                    do ieqatG = 1, atoms%neq(itypeG)
                      iGatom = iGatom + 1
                      do oqn_l = 0, atoms%lmax(itypeG)
                        do mqn_m = -oqn_l, oqn_l
                          lm = oqn_l * (oqn_l + 1) + 1 + mqn_m

                          ! q-dependant part of rho, q-dependant part of Vhar
                          integral = cmplx(0., 0.)
                          call Calc2ArgCmplxIntMT( atoms, itypeG, rho1MTDelta(:, lm, iGatom, idirR, iBatom), &
                            & vHar1MTDelta(:, lm, iGatom, idirC, iAatom), integral )

                          dynMatHFTest(idirR + (iBatom - 1) * 3, idirC + (iAatom - 1) * 3) = &
                            &dynMatHFTest(idirR + (iBatom - 1) * 3, idirC + (iAatom - 1) * 3) + integral
                          if (testoptimization) then
                            dynMatCompInt1(idirR + (iBatom - 1) * 3, idirC + (iAatom - 1) * 3) = &
                              &dynMatCompInt1(idirR + (iBatom - 1) * 3, idirC + (iAatom - 1) * 3) + integral
                          end if

                          ! q=0 part of rho, q-dependant part of Vhar
                          integral = cmplx(0., 0.)
                          call Calc2ArgCmplxIntMT( atoms, itypeG, -grRho0MT(:, lm, iGatom, idirR), &
                            & vHar1MTDelta(:, lm, iGatom, idirC, iAatom), integral )

                          !Is this basically the 2nd integral of (5.3.186)?
                          dynMatHFTest(idirR + (iBatom - 1) * 3, idirC + (iAatom - 1) * 3) = &
                            &dynMatHFTest(idirR + (iBatom - 1) * 3, idirC + (iAatom - 1) * 3) + integral
                          if (testoptimization) then
                            dynMatCompInt1(idirR + (iBatom - 1) * 3, idirC + (iAatom - 1) * 3) = &
                              &dynMatCompInt1(idirR + (iBatom - 1) * 3, idirC + (iAatom - 1) * 3) + integral
                          end if

                          ! q-dependant part of rho, q=0 part of Vhar
                          integral = cmplx(0., 0.)
                          call Calc2ArgCmplxIntMT( atoms, itypeG, rho1MTDelta(:, lm, iGatom, idirR, iBatom), &
                            & vHar1MTq0(:, lm, iGatom, idirC, iAatom), integral )

                          dynMatHFTest(idirR + (iBatom - 1) * 3, idirC + (iAatom - 1) * 3) = &
                            &dynMatHFTest(idirR + (iBatom - 1) * 3, idirC + (iAatom - 1) * 3) + integral
                          if (testoptimization) then
                            dynMatCompInt1(idirR + (iBatom - 1) * 3, idirC + (iAatom - 1) * 3) = &
                              &dynMatCompInt1(idirR + (iBatom - 1) * 3, idirC + (iAatom - 1) * 3) + integral
                          end if

                          if (testoptimization) then
                            ! q=0 parts of rho and Vhar
                            integral = cmplx(0., 0.)
                            call Calc2ArgCmplxIntMT( atoms, itypeG, -grRho0MT(:, lm, iGatom, idirR), &
                              & vHar1MTq0(:, lm, iGatom, idirC, iAatom), integral )

                            dynMatCompInt1(idirR + (iBatom - 1) * 3, idirC + (iAatom - 1) * 3) = &
                              &dynMatCompInt1(idirR + (iBatom - 1) * 3, idirC + (iAatom - 1) * 3) + integral
                          end if



                          ! q-dependant part of rho and q = 0 part of vext1
                          integral = cmplx(0., 0.)
                          call Calc2ArgCmplxIntMT( atoms, itypeG, rho1MTDelta(:, lm, iGatom, idirR, iBatom), &
                            & vExt1MTDelta(:, lm, iGatom, idirC, iAatom), integral )

                          dynMatHFTest(idirR + (iBatom - 1) * 3, idirC + (iAatom - 1) * 3) = &
                            &dynMatHFTest(idirR + (iBatom - 1) * 3, idirC + (iAatom - 1) * 3) + integral
                          if (testoptimization) then
                            dynMatCompInt1(idirR + (iBatom - 1) * 3, idirC + (iAatom - 1) * 3) = &
                              &dynMatCompInt1(idirR + (iBatom - 1) * 3, idirC + (iAatom - 1) * 3) + integral
                          end if

                          ! q = 0 part of rho and q-dependant part of Vext1
                          integral = cmplx(0., 0.)
                          call Calc2ArgCmplxIntMT( atoms, itypeG, -grRho0MT(:, lm, iGatom, idirR), &
                            & vExt1MTDelta(:, lm, iGatom, idirC, iAatom), integral )

                          !Is this basically the 1st integral of (5.3.186)?
                          dynMatHFTest(idirR + (iBatom - 1) * 3, idirC + (iAatom - 1) * 3) = &
                            &dynMatHFTest(idirR + (iBatom - 1) * 3, idirC + (iAatom - 1) * 3) + integral
                          if (testoptimization) then
                            dynMatCompInt1(idirR + (iBatom - 1) * 3, idirC + (iAatom - 1) * 3) = &
                              &dynMatCompInt1(idirR + (iBatom - 1) * 3, idirC + (iAatom - 1) * 3) + integral
                          end if

                          ! q-dependant part of rho and q=0 part of vExt1
                          integral = cmplx(0., 0.)
                          call Calc2ArgCmplxIntMT( atoms, itypeG, rho1MTDelta(:, lm, iGatom, idirR, iBatom), &
                            & vExt1MTq0(:, lm, iGatom, idirC, iAatom), integral )

                          dynMatHFTest(idirR + (iBatom - 1) * 3, idirC + (iAatom - 1) * 3) = &
                            &dynMatHFTest(idirR + (iBatom - 1) * 3, idirC + (iAatom - 1) * 3) + integral
                          if (testoptimization) then
                            dynMatCompInt1(idirR + (iBatom - 1) * 3, idirC + (iAatom - 1) * 3) = &
                              &dynMatCompInt1(idirR + (iBatom - 1) * 3, idirC + (iAatom - 1) * 3) + integral
                          end if

                          if (testoptimization) then
                            ! q=0 parts of rho and vExt1
                            integral = cmplx(0., 0.)
                            call Calc2ArgCmplxIntMT( atoms, itypeG, -grRho0MT(:, lm, iGatom, idirR), &
                              & vExt1MTq0(:, lm, iGatom, idirC, iAatom), integral )

                            dynMatCompInt1(idirR + (iBatom - 1) * 3, idirC + (iAatom - 1) * 3) = &
                              &dynMatCompInt1(idirR + (iBatom - 1) * 3, idirC + (iAatom - 1) * 3) + integral
                          end if


                          ! q-dependant parts of rho and vxc1
                          integral = cmplx(0., 0.)
                          call Calc2ArgCmplxIntMT( atoms, itypeG, rho1MTDelta(:, lm, iGatom, idirR, iBatom), &
                            & vXc1MTDelta(:, lm, iGatom, idirC, iAatom), integral )

                          dynMatHFTest(idirR + (iBatom - 1) * 3, idirC + (iAatom - 1) * 3) = &
                            &dynMatHFTest(idirR + (iBatom - 1) * 3, idirC + (iAatom - 1) * 3) + integral
                          if (testoptimization) then
                            dynMatCompInt1(idirR + (iBatom - 1) * 3, idirC + (iAatom - 1) * 3) = &
                              &dynMatCompInt1(idirR + (iBatom - 1) * 3, idirC + (iAatom - 1) * 3) + integral
                          end if

                          ! q = 0 part of rho and q-dependant part of Vxc1
                          integral = cmplx(0., 0.)
                          call Calc2ArgCmplxIntMT( atoms, itypeG, -grRho0MT(:, lm, iGatom, idirR), &
                            & vXc1MTDelta(:, lm, iGatom, idirC, iAatom), integral )

                          dynMatHFTest(idirR + (iBatom - 1) * 3, idirC + (iAatom - 1) * 3) = &
                            &dynMatHFTest(idirR + (iBatom - 1) * 3, idirC + (iAatom - 1) * 3) + integral
                          if (testoptimization) then
                            dynMatCompInt1(idirR + (iBatom - 1) * 3, idirC + (iAatom - 1) * 3) = &
                              &dynMatCompInt1(idirR + (iBatom - 1) * 3, idirC + (iAatom - 1) * 3) + integral
                          end if

                          ! q-dependant part of rho1 and q = 0 part of vxc1
                          integral = cmplx(0., 0.)
                          call Calc2ArgCmplxIntMT( atoms, itypeG, rho1MTDelta(:, lm, iGatom, idirR, iBatom), &
                            & vXc1MTq0(:, lm, iGatom, idirC, iAatom), integral )

                          dynMatHFTest(idirR + (iBatom - 1) * 3, idirC + (iAatom - 1) * 3) = &
                            &dynMatHFTest(idirR + (iBatom - 1) * 3, idirC + (iAatom - 1) * 3) + integral
                          if (testoptimization) then
                            dynMatCompInt1(idirR + (iBatom - 1) * 3, idirC + (iAatom - 1) * 3) = &
                              &dynMatCompInt1(idirR + (iBatom - 1) * 3, idirC + (iAatom - 1) * 3) + integral
                          end if

                          if (testoptimization) then
                            ! q=0 parts of vxc1 and of rho
                            integral = cmplx(0., 0.)
                            call Calc2ArgCmplxIntMT( atoms, itypeG, -grRho0MT(:, lm, iGatom, idirR), &
                              & vXc1MTq0(:, lm, iGatom, idirC, iAatom), integral )

                            dynMatCompInt1(idirR + (iBatom - 1) * 3, idirC + (iAatom - 1) * 3) = &
                              &dynMatCompInt1(idirR + (iBatom - 1) * 3, idirC + (iAatom - 1) * 3) + integral

                            integral = cmplx(0., 0.)
                            call Calc2ArgCmplxIntMT( atoms, itypeG, rho1MT(:, lm, iGatom, idirR, iBatom), vEff1MT(:, lm, iGatom, idirC, iAatom), integral)

                            dynMatCompInt1(idirR + (iBatom - 1) * 3, idirC + (iAatom - 1) * 3) = &
                              & dynMatCompInt1(idirR + (iBatom - 1) * 3, idirC + (iAatom - 1) * 3) - integral
                          end if

                        end do ! mqn_m
                      end do ! oqn_l
                    end do ! ieqatG
                  end do ! itypeG
                  ! add the splitted integrals that substitute the muffin integral of rho1 Veff1
                  dynMatPuInt(idirR + (iBatom - 1) * 3, idirC + (iAatom - 1) * 3) = &
                    & dynMatPuInt(idirR + (iBatom - 1) * 3, idirC + (iAatom - 1) * 3) + dynMatHFTest(idirR + (iBatom - 1) * 3, idirC + (iAatom - 1) * 3)

                  if (testoptimization) then
                    if (iAatom == iBatom) then
                      dynMatCompInt1MgradInt(idirR + (iBatom - 1) * 3, idirC + (iAatom - 1) * 3) = &
                                                                   & dynMatHFTest(idirR + (iBatom - 1) * 3, idirC + (iAatom - 1) * 3)

                      iGatom = 0
                      do itypeG = 1, atoms%ntype
                        do ieqatG = 1, atoms%neq(itypeG)
                          iGatom = iGatom + 1

                          do oqn_l = 0, atoms%lmax(iGatom)
                            do mqn_m = -oqn_l, oqn_l
                              lm = oqn_l * ( oqn_l + 1 ) + 1 + mqn_m

                              integral = cmplx(0., 0.)
                              call Calc2ArgCmplxIntMT( atoms, itypeG, rho1MT(:, lm, iGatom, idirR, iBatom), vEff1MT(:, lm, iGatom, idirC, iAatom), integral)
                              dynMatCompInt1MgradInt(idirR + (iBatom - 1) * 3, idirC + (iAatom - 1) * 3) = &
                                &dynMatCompInt1MgradInt(idirR + (iBatom - 1) * 3, idirC + (iAatom - 1) * 3) - integral

                              ! MT vanishes up to 9e-8
                              ! MT volume integral of grRho and grVext0
                              integral = cmplx(0., 0.)
                              call Calc2ArgCmplxIntMT( atoms, itypeG, grRho0MT(:, lm, iGatom, idirR), grVeff0MT(:, lm, idirC, iGatom), integral)
                              dynMatCompInt1MgradInt(idirR + (iBatom - 1) * 3, idirC + (iAatom - 1) * 3) = &
                                & dynMatCompInt1MgradInt(idirR + (iBatom - 1) * 3, idirC + (iAatom - 1) * 3) + integral

                            end do ! oqn_l
                          end do ! mqn_m

                        end do ! ieqatG
                      end do ! itypeG

                    end if ! iAatom == iBatom
                  end if ! test optimization

                end if ! for finite q optimization added

                ! Note: We combine here the valence and the core contribution of the integral grRho grVeff. When recasting the equations, such that we can implement it like this, a surface integral gradPsi* Veff Psi occurs in the valence part (cancelling a problematic term in the surface integrals and as well the complex conjugate Psi* Veff gradPsi occurs, also cancelling the respective counter term stemming from the Psi1 basis correction. As far as the core terms are concerned, the same terms occur.
                ! WE SUPPOSE SUCH TERMS TO BE ZERO AT THE MOMENT BECAUSE THE CORE WAVE FUNCTION IS SUPPOSED TO BE ZERO AT THE BOUNDARY OF THE MT BALL (CORE TAILS ARE NOT RELEVANT AT THE MOMENT)
                if ((iqpt == 1) .or. (.not.finiteQoptimization)) then
                  if (iAatom == iBatom) then

                    if ((idirR==1) .and. (idirC==1)) then
                      write(*, *) "Original method rho1 Veff1 2/2"
                    end if

                    ! Evaluate the MT part of rho1 vEff1IR
                    iGatom = 0
                    do itypeG = 1, atoms%ntype
                      do ieqatG = 1, atoms%neq(itypeG)
                        iGatom = iGatom + 1
                        do oqn_l = 0, atoms%lmax(itypeG)
                          do mqn_m = -oqn_l, oqn_l
                            lm = oqn_l * (oqn_l + 1) + 1 + mqn_m

                !            integral = cmplx(0., 0.)
                !            call Calc2ArgCmplxIntMT( atoms, itypeG, grRho0MT(:, lm, iGatom, idirR), grVeff0MT(:, lm, idirC, iGatom), integral)

                !            dynMatPuInt(idirR + (iBatom - 1) * 3, idirC + (iAatom - 1) * 3) = &
                !              & dynMatPuInt(idirR + (iBatom - 1) * 3, idirC + (iAatom - 1) * 3) - integral

                          end do ! mqn_m
                        end do ! oqn_l
                      end do ! ieqatG
                    end do ! itypeG
                  end if ! iAtom = iBatom
                end if ! q = 0 or no optimization enabled

              end do ! idirR
            end do ! iDeqatB
          end do ! iDtypeB
        end do ! idirC
      end do ! iDatomA
    end do ! iDtypeA

    if (any(abs(dynMatCompInt1(:, :)) > 1e-8)) call juDFT_error('rho1 Veff1 integral not the same as optimization', calledby='EvalIntRho1Veff1')
    if (any(abs(dynMatCompInt1MgradInt(:, :)) > 1e-8)) call juDFT_error('rho1 Veff1 integral -gradRho Veff1 integral not the same as optimization', calledby='EvalIntRho1Veff1')

    if (.false.) then
      write(*, '(a)') 'Test Matrix'
      write(*, '(3(2(es16.8,1x),3x))') dynMatCompInt1(1, :)
      write(*, '(3(2(es16.8,1x),3x))') dynMatCompInt1(2, :)
      write(*, '(3(2(es16.8,1x),3x))') dynMatCompInt1(3, :)

      write(*, '(a)') 'Test Matrix'
      write(*, '(3(2(es16.8,1x),3x))') dynMatCompInt1MgradInt(1, :)
      write(*, '(3(2(es16.8,1x),3x))') dynMatCompInt1MgradInt(2, :)
      write(*, '(3(2(es16.8,1x),3x))') dynMatCompInt1MgradInt(3, :)
    end if

  end subroutine EvalIntRho1Veff1


  ! Evaluate Pulay (H - eps) braket
  !>-------------------------------------------------------------------------------------------------------------------------------------
  !> @author
  !> Christian-Roman Gerhorst
  !>
  !> @brief
  !>
  !>
  !> @details
  !>
  !> @note
  !> Note: One could make this routine faster by pulling out the Gaunt coefficients from the gradient routine and using it after the > calculation of the radial integrals. In this way one would only have to calculate the radial integrals for every l and then multiply
  !> the Gaunt coefficients. This would be more efficient but we choose the traditional way to minimize the error probability. After this
  !> programm is running one has a good benchmark.
  !>
  !> @todo
  !> Account for the spin.
  !>-------------------------------------------------------------------------------------------------------------------------------------
  subroutine EvalPuHepsBraKetsVal(fmpi, noco, nococonv, oneD, atoms, input, cell, sym, lathar, stars, kpts, qpts, usdus, results, iqpt, nRadFun, vEff0IR, vEff0MtLh, clnu_atom,&
      & nmem_atom, mlh_atom, rbas1, rbas2, El, mapKpq2K, nobd, nv, gBasVec, ilst, kveclo, z, iloTable, eig, kpq2kPrVec, dynMatPu )

    use m_types
    use m_jpSetupDynMatSF, only : CalcChannelsGrFlpNat, CalcChannelsGrGrtFlpNat, readInz1, CalcHnGrV0Varphi, CalcHGrVarphi
    use m_dfpt_init, only : Derivative

    implicit none

    ! Type parameter
    type(t_mpi),                  intent(in)  :: fmpi
    type(t_nococonv),                intent(in)  :: nococonv
    type(t_noco),                intent(in)  :: noco
    type(t_oneD),                intent(in)  :: oneD
    type(t_atoms),              intent(in)  :: atoms
    type(t_input),              intent(in)  :: input
    type(t_sphhar),             intent(in)  :: lathar
    type(t_stars),              intent(in)  :: stars
    type(t_cell),               intent(in)  :: cell
    type(t_sym),                intent(in)  :: sym
    type(t_kpts),               intent(in)  :: kpts
    type(t_kpts),               intent(in)  :: qpts
    type(t_usdus),              intent(in)  :: usdus
    type(t_results),            intent(in)  :: results

    ! Scalar parameter
    integer,                    intent(in)  :: iqpt

    ! Array parameter
    ! rbas1      : large component of radial basis functions
    ! rbas2      : small component of radial basis functions
    ! nRadFun    : number of radial functions for certain l and k-point
    ! vEff0IR    : unperturbed effective potential in interstitial
    ! vEff0MTLh  : unperturbed effective potential in muffin-tin
    ! nmem_atom  : nmem given that the non-representative MTs are not rotated
    ! mlh_atom   : mlh given that the non-representative MTs are not rotated
    ! clnu_atom  : clnu given that the non-representative MTs are not rotated
    ! El         : LAPW energy parameters
    ! mapKpq2K   : index of k' = k + q, i.e. which k-point set index has the (maybe backfolded) sum of k and q
    ! nobd       : number of occupied eigenvalues
    ! nv         : number of all eigenvalues
    ! GbasVec    : G-basis vector
    ! ilst       : index array to avoid redundant G-basis vectors
    ! kveclo     : index array for selecting correct G-vector for respective LO
    ! z          : Kohn-Sham wavefunction expansion coefficients
    ! iloTable   : index array to reorganize LOs
    ! eig        : Kohn-Sham eigen energies
    ! kpq2kPrVec : contains backfolding vector, if k' = k + q required a backfolding in the unit cell.

    real,                       intent(in)  :: rbas1(:, :, 0:, :)
    real,                       intent(in)  :: rbas2(:, :, 0:, :)
    integer,                    intent(in)  :: nRadFun(0:, :)
    complex,                    intent(in)  :: vEff0IR(:,:)
    real,                       intent(in)  :: vEff0MtLh(:, 0:, :)
    integer,                    intent(in)  :: nmem_atom(0:, :)
    integer,                    intent(in)  :: mlh_atom(:, 0:, :)
    complex,                    intent(in)  :: clnu_atom(:, 0:, :)
    real,                       intent(in)  :: El(:, 0:, :, :)
    integer,                    intent(in)  :: mapKpq2K(:, :)
    integer,                    intent(in)  :: nobd(:, :)
    integer,                    intent(in)  :: nv(:, :)
    integer,                    intent(in)  :: GbasVec(:, :)
    integer,                    intent(in)  :: ilst(:, :, :)
    integer,                    intent(in)  :: kveclo(:,:)
    complex,                   intent(in)  :: z(:,:,:,:)
    integer,                    intent(in)  :: iloTable(:, 0:, :)
    real,                       intent(in)  :: eig(:, :, :)
    integer,                    intent(in)  :: kpq2kPrVec(:, :, :)
    complex, allocatable,       intent(out) :: dynMatPu(:, :)

    ! Scalar variables
    integer                                 :: lmpMax
    integer                                 :: nRadFunMax
    integer                                 :: iatom
    integer                                 :: itype
    integer                                 :: imesh
    integer                                 :: oqn_l
    integer                                 :: iradf
    integer                                 :: mqn_m2PrR
    integer                                 :: mqn_m
    integer                                 :: ichan
    integer                                 :: lmp
    integer                                 :: mqn_m2PrC
    integer                                 :: ieqat
    integer                                 :: ptsym
    integer                                 :: ilh
    integer                                 :: lm_pre
    integer                                 :: imem
    integer                                 :: lm
    integer                                 :: ikpt
    integer                                 :: ii
    integer                                 :: jj
    logical                                 :: eps1DynMatPulTestSw
    logical                                 :: testComp2ndN1stOrdBasFuncSw
    logical                                 :: testCompTerm3rdBraKetsVarKet
    logical                                 :: testCompTerm3rdBraKetsVarBra
    logical                                 :: dynMatPu3rdBraKetHepsSw

    ! Array variables
    integer,       allocatable              :: lmpT(:)
    complex,       allocatable              :: z1nG(:, :, :, :)
    real,          allocatable              :: varphi1(:, :, :)
    real,          allocatable              :: varphi2(:, :, :)
    real,          allocatable              :: r2(:)
    real,          allocatable              :: delrVarphi1(:, :, :)
    real,          allocatable              :: delrVarphi2(:, :, :)
    integer,       allocatable              :: grVarphiChLout(:, :)
    integer,       allocatable              :: grVarphiChMout(:, :)
    real,          allocatable              :: grVarphiCh1(:, :, :, :)
    real,          allocatable              :: grVarphiCh2(:, :, :, :)
    real,          allocatable              :: delrGrVarphiCh1(:, :, :, :)
    real,          allocatable              :: delrGrVarphiCh2(:, :, :, :)
    complex,       allocatable              :: vEff0MtSpH(:, :)
    complex,       allocatable              :: hVarphi(:, :, :, :)
    real,          allocatable              :: varphiVarphi(:, :, :)
    complex,       allocatable              :: varphiHvarphi(:, :, :)
    real,          allocatable              :: grVarphiVarphi(:, :, :, :)
    complex,       allocatable              :: grVarphiHpreGrtVarphi(:, :, :, :, :)
    complex,       allocatable              :: grVarphiHVarphi(:, :, :, :)
    complex,       allocatable              :: vEff0NsphGrVarphi(:, :, :, :, :)
    complex,       allocatable              :: r2grVeff0SphVarphi(:, :, :, :, :)
    complex,       allocatable              :: r2grVeff0Varphi(:, :, :, :, :)
    complex,       allocatable              :: hGrVarphi(:, :, :, :, :)
    real,          allocatable              :: varphiVarphiDummy(:, :, :)
    complex,       allocatable              :: varphiGrVeff0SphVarphi(:, :, :, :)
    complex,       allocatable              :: varphiGrVeff0Varphi(:, :, :, :)
    complex,       allocatable              :: varphiHGrvarphi(:, :, :, :)
    real,          allocatable              :: varphiGrVarphi(:, :, :, :)
    real                                    :: kExt(1:3)
    real                                    :: gbasExt(1:3)
    integer :: iBas
    integer :: iband
    integer :: idir
    integer :: lmaxBra

    allocate(dynMatPu(3 * atoms%nat, 3 * atoms%nat))
    dynMatPu(:, :) = cmplx(0., 0.)


    ! Quantities for initialization
    allocate( lmpT(atoms%ntype) )
    lmpT(:) = 0
    do itype = 1, atoms%ntype
      lmpT(itype) = sum( [ ( (2 * oqn_l + 1)* nRadFun(oqn_l, itype), oqn_l = 0, atoms%lmax(itype) ) ] )
    end do ! itype
    lmpMax     = maxval( lmpT(:) )
    nRadFunMax = maxval( nRadFun(:, :) )

    allocate( z1nG(SIZE(z(:,1,1,1)), 3, atoms%nat, maxval(nobd(:, :))) )
    allocate( varphi1(atoms%jmtd, nRadFunMax, 0:atoms%lmaxd), varphi2(atoms%jmtd, nRadFunMax, 0:atoms%lmaxd) )
    allocate( r2(atoms%jmtd) )
    allocate( delrVarphi1( atoms%jmtd, nRadFunMax, 0:atoms%lmaxd), delrVarphi2( atoms%jmtd, nRadFunMax, 0:atoms%lmaxd) )
    allocate( grVarphiCh1(atoms%jmtd, 2, lmpMax, -1:1), grVarphiCh2(atoms%jmtd, 2, lmpMax, -1:1), &
            & grVarphiChLout(2, 0:atoms%lmaxd), grVarphiChMout(-atoms%lmaxd:atoms%lmaxd, -1:1) )
    allocate( hVarphi(2, atoms%jmtd, 0:(atoms%lmaxd + 3)**2 - 1, lmpMax)  )
    allocate( varphiVarphi(lmpMax, lmpMax, atoms%ntype), varphiHvarphi(lmpMax, lmpMax, atoms%nat))
    allocate( grVarphiVarphi(lmpMax, lmpMax, -1:1, atoms%ntype), grVarphiHVarphi(lmpMax, lmpMax, -1:1, atoms%nat ) )
    allocate( vEff0MtSpH( atoms%jmtd, 0:atoms%lmaxd*(atoms%lmaxd+2)), vEff0NsphGrVarphi(2, atoms%jmtd, 0:(atoms%lmaxd + 3)**2 - 1, lmpMax, -1:1), &
            & r2grVeff0SphVarphi(2, atoms%jmtd, 0:(atoms%lmaxd + 3)**2 - 1, lmpMax, 3) )
    allocate( r2grVeff0Varphi(2, atoms%jmtd, 0:(atoms%lmaxd + 3)**2 - 1, lmpMax, 3))
    allocate( hGrVarphi(2, atoms%jmtd, (atoms%lmaxd + 2)**2, lmpMax, -1:1))
    allocate( varphiVarphiDummy(lmpMax, lmpMax, atoms%ntype) )
    allocate( varphiGrVeff0SphVarphi(lmpMax, lmpMax, atoms%nat, 3) )
    allocate( varphiGrVeff0Varphi(lmpMax, lmpMax, atoms%nat, 3) )
    allocate( varphiHGrvarphi(lmpMax, lmpMax, atoms%nat, -1:1) )
    allocate( varphiGrVarphi(lmpMax, lmpMax, -1:1, atoms%ntype))

    hGrVarphi(:, :, :, :, :) = cmplx(0., 0.)
    varphiGrVeff0SphVarphi(:, :, :, :) = cmplx(0., 0.)
    varphiGrVeff0Varphi(:, :, :, :) = cmplx(0., 0.)
    varphiVarphiDummy(:, :, :) = 0.
    varphiHGrvarphi(:, :, :, :) = cmplx(0., 0.)
    varphiVarphi(:, :, :)            = 0.
    varphiHvarphi(:, :, :)              = cmplx(0., 0.)
    grVarphiVarphi(:, :, :,  :)         = 0.
    grVarphiHVarphi(:, :, :, :)         = cmplx(0., 0.)
    vEff0NsphGrVarphi(:, :, :, :, :)    = cmplx(0., 0.)
    r2grVeff0SphVarphi(:, :, :, :, :)   = cmplx(0., 0.)
    r2grVeff0Varphi(:, :, :, :, :)   = cmplx(0., 0.)
    varphiGrVarphi(:, :, :, :) = 0.

    iatom = 0
    do itype = 1, atoms%ntype

    lmaxBra = atoms%lmax(itype)

      ! Precalculate radial Jacobi determinant for later integrals
      r2(:) = 0.
      do imesh = 1, atoms%jri(itype)
        r2(imesh) = atoms%rmsh(imesh, itype) * atoms%rmsh(imesh, itype)
      end do ! imesh

      varphi1(:, :, :) = 0.
      varphi2(:, :, :) = 0.
      delrVarphi1(:, :, :) = 0.
      delrVarphi2(:, :, :) = 0.
      do oqn_l = 0, atoms%lmax(itype)
        do iradf = 1, nRadFun(oqn_l, itype)
          do imesh = 1, atoms%jri(itype)
            ! In Fleur the radial solutions u_lp are multiplied by a factor r by default to avoid an additional multiplication of the
            ! Jacobi determinant r^2 in radial integrals given spherical coordinates.
            varphi1(imesh, iradf, oqn_l) = rbas1(imesh, iradf, oqn_l, itype) / atoms%rmsh(imesh, itype)
            varphi2(imesh, iradf, oqn_l) = rbas2(imesh, iradf, oqn_l, itype) / atoms%rmsh(imesh, itype)
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

      deallocate( delrVarphi1, delrVarphi2 )

      ! Calculate H |varphi>
      do ieqat = 1, atoms%neq(itype)
        iatom = iatom + 1
        !todo block begin                             place this block into calchngrv0varphi
        ! Expand the coefficients of the lattice-harmonic given potential into spherical-harmonic coefficients for the given atom.
        vEff0MtSpH(:, :) = cmplx(0.0, 0.0)
        ptsym = sym%ntypsy(iatom)
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

        !todo block end.............................................

        hVarphi = cmplx(0.0, 0.0)
        call CalcHnGrV0Varphi( atoms, sym, lathar, itype, iatom, lmpMax, El, varphi1, varphi2, nRadFun, vEff0MtSpH, vEff0MtLh, clnu_atom, &
          & nmem_atom, mlh_atom, grVarphiCh1, grVarphiCh2, grVarphiChLout, grVarphiChMout, hVarphi, vEff0NsphGrVarphi, r2grVeff0SphVarphi, r2grVeff0Varphi )

        do mqn_m2PrC = -1, 1
          call CalcHGrVarphi( atoms, itype, mqn_m2PrC, lmpMax, lmaxBra, grVarphiChMout, nRadFun, grVarPhiCh1, grVarPhiCh2,        &
                                                      & grVarphiChLout, vEff0NsphGrVarphi, El, lmpT, hGrVarphi )
        end do ! mqn_m2PrC

        deallocate( vEff0MtSpH )

        do mqn_m2PrC = -1, 1
          varphiVarphiDummy(:, :, :) = 0.
          !todo attention where the hGrVarphi starts!
          call CalcScalBasfMatElems( atoms, itype, iatom, nRadFun, r2, varphi1, varphi2, hGrVarphi(:, :, :, :, mqn_m2PrC), varphiVarphiDummy, varphiHGrvarphi(:, :, :, mqn_m2PrC) )
        end do ! mqn_m2PrC


        ! Calculate all radial integrals with no gradients.jjjj
        ! Calculate scalar basis function matrix elements

        !call CalcScalBasfMatElems( atoms, itype, iatom, nRadFun, r2, varphi1, varphi2, hVarphi, varphiVarphi, varphiHvarphi )
        !!! Symmetrized kinetic energy:
        call CalcScalBasfMatElems( atoms, itype, iatom, nRadFun, r2, varphi1, varphi2, hVarphi, varphiVarphi, varphiHvarphi, usdus )

        ! not optimal with the first dimension
        ! This should be R for row but we want to avoid overhead setting up a new loop.
        !todo also have a look at calcproj on grVarphi in deprecated tests
        call CalcVecBasfMatElems( atoms, itype, 2, nRadFun, r2, grVarphiChLout, grVarphiChMout, varphi1, varphi2, &
          & grVarPhiCh1, grVarphiCh2, hVarphi, grVarphiVarphi(:, :, :, itype), grVarphiHVarphi(:, :, :, iatom) )

        do mqn_m2PrC = -1, 1
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
    deallocate( varphi1, varphi2, grVarphiCh1, grVarphiCh2, grVarphiChLout, grVarphiChMout, &
      & r2grVeff0SphVarphi, vEff0NsphGrVarphi, r2, hVarphi, varphiVarphiDummy, hGrVarphi )

    write(*, *) 'All quantities for DynMatPu have been prepared'
    ! This is what Aaron Klueppelberg calls first and second braket.
    eps1DynMatPulTestSw = .false.
    testComp2ndN1stOrdBasFuncSw = .false.
    testCompTerm3rdBraKetsVarKet = .false.
    testCompTerm3rdBraKetsVarBra = .false.
    dynMatPu3rdBraKetHepsSw = .false.
! todo we still have to make it band dependent for metals later
    write(*, *) 'bug for other systems than neon with more than one atom, also some lines below'
    write(*, *) 'why not z1 at nobd(ikpq)?'
    do ikpt = 1, kpts%nkpt
      z1nG = cmplx(0.0, 0.0)
      call ReadInz1( atoms, ikpt, iqpt, mapKpq2K(ikpt, iqpt), nobd, nv, z1nG )
!      kExt(1:3) = matmul(cell%bmat, kpts%bk(1:3, ikpt))
!      gbasExt(:) = 0.
!      do iBas = 1, nv(1, ikpt) !+ atoms%nlotot
!        gbasExt(1:3) = matmul( cell%bmat(1:3, 1:3), gBasVec(1:3, ilst(iBas, ikpt, 1)))
!        do iband = 1, nobd(ikpt, 1)
!          do idir = 1, 3
!            z1nG(iBas, idir, 1, iband) = -iu * ( kExt(idir) + gbasExt(idir) ) * z(iBas, iband, ikpt, 1)
!          end do ! iBas
!        end do ! idir
!      end do ! iband
      !(5.3.179/180)
      !call Add2ndOrdWfPulayBraKets2DynMat( atoms, cell, dimens, sym, kpts, qpts, usdus, results, iqpt, ikpt, nRadFunMax, lmpMax, eps1DynMatPulTestSw, testComp2ndN1stOrdBasFuncSw, mapKpq2K, eig, nobd, nv, &
      !  & gBasVec, ilst, kveclo, z, z1nG, nRadFun, iloTable, varphiVarphi, varphiHvarphi, grVarphiVarphi, grVarphiHVarphi, &
      !  &  lmpT, varphiGrVeff0SphVarphi, varphiHGrvarphi, varphiGrVarphi, varphiGrVeff0Varphi, dynMatPu )
      !write(470,*) dynMatPu
      !(5.3.181)
      !call Add1stOrdWfPulayBraKets2DynMat( atoms, kpts, qpts, sym, dimens, cell, usdus, stars, results, ikpt, iqpt, lmpMax, nRadFunMax, nv, GbasVec, ilst, kveclo, &
      !& mapKpq2K, nobd, z, z1nG, iloTable, grVarphiVarphi, nRadFun, eig, kpq2kPrVec, grVarphiHvarphi, &
      !& varphiVarphi, varphiHvarphi, vEff0IR, lmpT, eps1DynMatPulTestSw, testCompTerm3rdBraKetsVarBra, testCompTerm3rdBraKetsVarKet, dynMatPu3rdBraKetHepsSw, varphiGrVeff0SphVarphi, varphiHGrvarphi, varphiGrVarphi, dynMatPu )
      !write(470,*) dynMatPu
      call AddAlexPulayBraKets2DynMat( fmpi, noco, nococonv, oneD, atoms, input, kpts, qpts, sym, cell, usdus, stars, results, ikpt, iqpt, lmpMax, nRadFunMax, nv, GbasVec, ilst, kveclo, &
      & mapKpq2K, nobd, z, z1nG, iloTable, grVarphiVarphi, nRadFun, eig, kpq2kPrVec, grVarphiHvarphi, &
      & varphiVarphi, varphiHvarphi, vEff0IR, lmpT, eps1DynMatPulTestSw, testCompTerm3rdBraKetsVarBra, testCompTerm3rdBraKetsVarKet, dynMatPu3rdBraKetHepsSw, &
      & varphiGrVeff0SphVarphi, varphiHGrvarphi, varphiGrVarphi, dynMatPu )
    end do ! ikpt
    !write(*, '(a)') '1st, 2nd and 3rd Pulay Braket term added to dynamical matrix.'
    write(*, '(a)') 'Alex Pulay Braket terms added to dynamical matrix.'

  end subroutine EvalPuHepsBraKetsVal

  subroutine CalcGrVarphiHepsGrtVarphiElem( atoms, itype, iatom, mqn_m2PrC, lmpMax, nRadFun, grVarphiChLout, grVarphiChMout, r2, grVarPhiCh1, &
      & grVarPhiCh2, El, lmpT, vEff0NsphGrVarphi, r2grVeff0SphVarphi, grVarphiGrtVarphi, grVarphiHpreGrtVarphi,&
      & grVarphiGrtVeff0SphVarphi)

    use m_types, only : t_atoms
    !use m_intgr, only : intgr3NoIntp
    use m_intgr, only : intgr3!LinIntp ! TODO: Is this ok?
    !use m_intgr, only : intgr3

    implicit none

    type(t_atoms),              intent(in)  :: atoms

    integer,                    intent(in)  :: itype
    integer,                    intent(in)  :: iatom
    integer,                    intent(in)  :: mqn_m2PrC
    integer,                    intent(in)  :: lmpMax

    integer,                    intent(in)  :: nRadFun(0:, :)
    integer,                    intent(in)  :: grVarphiChLout(:, 0:)
    integer,                    intent(in)  :: grVarphiChMout(-atoms%lmaxd:, -1:)
    real,                       intent(in)  :: r2(:)
    real,                       intent(in)  :: grVarphiCh1(:, :, :, -1:)
    real,                       intent(in)  :: grVarphiCh2(:, :, :, -1:)
    real,                       intent(in)  :: El(:, 0:, :, :)
    integer,                    intent(in)  :: lmpT(:)
    complex,                    intent(in)  :: vEff0NsphGrVarphi(:, :, 0:, :, -1:)
    complex,                    intent(in)  :: r2grVeff0SphVarphi(:, :, 0:, :, :)
    real,                     intent(inout) :: grVarphiGrtVarphi(:, :, -1:, -1:, :)
    complex,                    intent(inout) :: grVarphiHpreGrtVarphi(:, :, -1:, -1:, :)
    complex,                    intent(inout) :: grVarphiGrtVeff0SphVarphi(:, :, -1:, :, :)

    integer                                 :: lmp
    integer                                 :: oqn_l
    integer                                 :: lm_pre
    integer                                 :: mqn_m
    integer                                 :: lm
    integer                                 :: iradf
    integer                                 :: oqn_l3Pr
    integer                                 :: iradfKet
    integer                                 :: imesh
    real                                    :: integralR
    real                                    :: integralI
    integer                                 :: lmp1Pr
    integer                                 :: oqn_l1Pr
    integer                                 :: lm1Pr_pre
    integer                                 :: mqn_m1Pr
    integer                                 :: lm1Pr
    integer                                 :: iradf1Pr
    integer                                 :: ichanPr
    integer                                 :: ichan
    integer                                 :: lm3Pr
    integer                                 :: mqn_m2PrR
    integer                                 :: mqn_m3Pr

    real,          allocatable              :: intgrdR(:)
    real,          allocatable              :: intgrdI(:)
    complex,       allocatable              :: grVarphiVeff0NsphGrvarphi(:, :)
    complex,       allocatable              :: hsphActedOnGrVarphi(:, :)

    allocate( intgrdR(atoms%jri(itype)), intgrdI(atoms%jri(itype)) )
    allocate( grVarphiVeff0NsphGrvarphi(lmpMax, lmpMax) )
    allocate( hsphActedOnGrVarphi(lmpMax, lmpMax) )

    do mqn_m2PrR = -1, 1
      grVarphiVeff0NsphGrvarphi = cmplx(0., 0.)
      hsphActedOnGrVarphi(:, :) = 0
      lmp = 0
      lm = -1
      do oqn_l = 0, atoms%lmax(itype)
        do mqn_m = -oqn_l, oqn_l
          lm = lm + 1
          do iradf = 1, nRadFun(oqn_l, itype)
            lmp = lmp + 1
            lmp1Pr = 0
            lm1Pr = -1
            do oqn_l1Pr = 0, atoms%lmax(itype)
              do mqn_m1Pr = -oqn_l1Pr, oqn_l1Pr
                mqn_m3Pr = grVarphiChMout(mqn_m1Pr, mqn_m2PrR)
                lm1Pr = lm1Pr  + 1
                do iradf1Pr = 1, nRadFun(oqn_l1Pr, itype)
                  lmp1Pr = lmp1Pr + 1
                  ! overlap of grVarphi grVarphi
                  do ichanPr = 1, 2
                    oqn_l3Pr = grVarphiChLout(ichanPr, oqn_l1Pr)
                    ! We have to catch channels where l < 0 and where m"' does not contribute, lmax + 1 is allowed due to the gradient
                    if ( ( abs(mqn_m3Pr) > oqn_l3Pr )  .or. ( oqn_l3Pr < 0 )  )  cycle
                    lm3Pr = oqn_l3Pr * (oqn_l3Pr + 1) + mqn_m3Pr
                    do ichan = 1, 2
                      ! < grVarphi | grtVarphi >
                      if ( ( abs(grVarphiChMout(mqn_m, mqn_m2PrC)) > grVarphiChLout(ichan, oqn_l) ) .or. &
                                                                                       & (grVarphiChLout(ichan, oqn_l) < 0) ) cycle
                      if ( lm3Pr /= grVarphiChLout(ichan, oqn_l) * (grVarphiChLout(ichan, oqn_l) + 1) &
                                                                                        & + grVarphiChMout(mqn_m, mqn_m2PrC) ) cycle
                      intgrdR(:) = 0.
                      do imesh = 1, atoms%jri(itype)
                        intgrdR(imesh) = r2(imesh) * &
                        & ( grVarphiCh1(imesh, ichanPr, lmp1Pr, mqn_m2PrR) * grVarphiCh1(imesh, ichan, lmp, mqn_m2PrC) &
                        & + grVarphiCh2(imesh, ichanPr, lmp1Pr, mqn_m2PrR) * grVarphiCh2(imesh, ichan, lmp, mqn_m2PrC) )
                      end do ! imesh
                      call intgr3(intgrdR(:), atoms%rmsh(1, itype), atoms%dx(itype), atoms%jri(itype), integralR)! TODO: Is this ok?
                      !call intgr3(intgrdR(:), atoms%rmsh(1, itype), atoms%dx(itype), atoms%jri(itype), integralR)
!
                      grVarphiGrtVarphi(lmp1Pr, lmp, mqn_m2prR, mqn_m2prC, itype) = &
                                                           & grVarphiGrtVarphi(lmp1Pr, lmp, mqn_m2prR, mqn_m2prC, itype) + integralR
                    end do ! ichan

                    intgrdR(:) = 0.
                    intgrdI(:) = 0.
                    do imesh = 1, atoms%jri(itype)
                      intgrdR(imesh) = r2(imesh) * real ( &
                        &   (grVarPhiCh1(imesh, ichanPr, lmp1Pr, mqn_m2PrR) * vEff0NsphGrVarphi(1, imesh, lm3Pr, lmp, mqn_m2PrC) &
                        &  + grVarPhiCh2(imesh, ichanPr, lmp1Pr, mqn_m2PrR) * vEff0NsphGrVarphi(2, imesh, lm3Pr, lmp, mqn_m2PrC) ) )
                      intgrdI(imesh) = r2(imesh) * aimag( &
                        &   (grVarPhiCh1(imesh, ichanPr, lmp1Pr, mqn_m2PrR) * vEff0NsphGrVarphi(1, imesh, lm3Pr, lmp, mqn_m2PrC) &
                        &  + grVarPhiCh2(imesh, ichanPr, lmp1Pr, mqn_m2PrR) * vEff0NsphGrVarphi(2, imesh, lm3Pr, lmp, mqn_m2PrC) ) )
                    end do ! imesh
                    call intgr3(intgrdR(:), atoms%rmsh(1, itype), atoms%dx(itype), atoms%jri(itype), integralR)! TODO: Is this ok?
                    call intgr3(intgrdI(:), atoms%rmsh(1, itype), atoms%dx(itype), atoms%jri(itype), integralI )! TODO: Is this ok?
                    ! < grVarphi | Veff0nsph | grVarphi >
                    grVarphiVeff0NsphGrvarphi(lmp1Pr, lmp) = grVarphiVeff0NsphGrvarphi(lmp1Pr, lmp) + cmplx(integralR, integralI)

                    intgrdR(:) = 0.
                    intgrdI(:) = 0.
                    do imesh = 1, atoms%jri(itype)
                      intgrdR(imesh) = real ( &
                        &   (grVarPhiCh1(imesh, ichanPr, lmp1Pr, mqn_m2PrR) * r2grVeff0SphVarphi(1, imesh, lm3Pr, lmp, mqn_m2PrC + 2) &
                        &  + grVarPhiCh2(imesh, ichanPr, lmp1Pr, mqn_m2PrR) * r2grVeff0SphVarphi(2, imesh, lm3Pr, lmp, mqn_m2PrC + 2) ) )
                      intgrdI(imesh) = aimag( &
                        !todo why without iatom
                        &   (grVarPhiCh1(imesh, ichanPr, lmp1Pr, mqn_m2PrR) * r2grVeff0SphVarphi(1, imesh, lm3Pr, lmp, mqn_m2PrC + 2) &
                        &  + grVarPhiCh2(imesh, ichanPr, lmp1Pr, mqn_m2PrR) * r2grVeff0SphVarphi(2, imesh, lm3Pr, lmp, mqn_m2PrC + 2) ) )
                    end do ! imesh
                    call intgr3(intgrdR(:), atoms%rmsh(1, itype), atoms%dx(itype), atoms%jri(itype), integralR)! TODO: Is this ok?
                    call intgr3(intgrdI(:), atoms%rmsh(1, itype), atoms%dx(itype), atoms%jri(itype), integralI)! TODO: Is this ok?
                    ! < grVarphi| grTveff0Sph |varphi >
                    ! idir = mqn_m2PrC + 2 due to performance reasons
                    grVarphiGrtVeff0SphVarphi(lmp1Pr, lmp, mqn_m2PrR, mqn_m2PrC + 2, iatom) = &
                             & grVarphiGrtVeff0SphVarphi(lmp1Pr, lmp, mqn_m2PrR, mqn_m2PrC + 2, iatom) + cmplx(integralR, integralI)
                  end do ! ichanPr
                  if (iradf < 3) then
                    hsphActedOnGrVarphi(lmp1Pr, lmp) = hsphActedOnGrVarphi(lmp1Pr, lmp) &
                                            & + El(1, oqn_l, itype, 1) * grVarphiGrtVarphi(lmp1Pr, lmp, mqn_m2prR, mqn_m2prC, itype)
                    if (iradf == 2) then
                      hsphActedOnGrVarphi(lmp1Pr, lmp) = hsphActedOnGrVarphi(lmp1Pr, lmp) &
                                                                     & + grVarphiGrtVarphi(lmp1Pr, lmp, mqn_m2prR, mqn_m2prC, itype)
                    end if ! (iradf == 2)
                  !else
                  !  hsphActedOnGrVarphi(lmp1Pr, lmp) = hsphActedOnGrVarphi(lmp1Pr, lmp) &
                  !                      & + El(iradf, oqn_l, itype, 1) * grVarphiGrtVarphi(lmp1Pr, lmp, mqn_m2prR, mqn_m2prC, itype)
                  end if ! (iradf < 3)
                end do ! iradf1Pr
              end do ! mqn_m1Pr
            end do ! oqn_l1Pr
          end do ! iradf
        end do ! mqn_m
      end do ! oqn_l
      grVarphiHpreGrtVarphi(1:lmpT(itype), 1:lmpT(itype), mqn_m2PrR, mqn_m2PrC, iatom) = &
        &   hsphActedOnGrVarphi(1:lmpT(itype), 1:lmpT(itype)) + grVarphiVeff0NsphGrvarphi(1:lmpT(itype), 1:lmpT(itype))
!      grVarphiHpreGrtVarphi(1:lmpT(itype), 1:lmpT(itype), mqn_m2PrR, mqn_m2PrC, iatom) = &
!        &   grVarphiVeff0NsphGrvarphi(1:lmpT(itype), 1:lmpT(itype))
    end do ! mqn_m2PrR

  end subroutine CalcGrVarphiHepsGrtVarphiElem

  subroutine CalcVecBasfMatElems( atoms, itype, chanMax, nRadFun, r2, grVarphiChLout, grVarphiChMout, varphi1, varphi2, grVarPhiCh1,&
      & grVarphiCh2, hVarphi, grVarphiVarphi, grVarphiHVarphi )

    use m_types, only : t_atoms
!    use m_intgr, only : intgr3NoIntp
    use m_intgr, only : intgr3!LinIntp! TODO: Is this ok?
!    use m_intgr, only : intgr3

    implicit none

    type(t_atoms),             intent(in)  :: atoms

    !todo check the starting values of the intent ins
    integer,                   intent(in)  :: itype
    integer,                   intent(in)  :: chanMax

    integer,                   intent(in)  :: nRadFun(0:, :)
    real,                      intent(in)  :: r2(:)
    integer,                   intent(in)  :: grVarphiChLout(:, 0:)
    integer,                   intent(in)  :: grVarphiChMout(-atoms%lmaxd:, -1:)
    real,                      intent(in)  :: varphi1(:, :, 0:)
    real,                      intent(in)  :: varphi2(:, :, 0:)
    real,                      intent(in)  :: grVarphiCh1(:, :, :, -1:)
    real,                      intent(in)  :: grVarphiCh2(:, :, :, -1:)
    complex,                   intent(in)  :: hVarphi(:, :, 0:, :)
    real,                      intent(out) :: grVarphiVarphi(:, :, -1:)
    complex,                   intent(out) :: grVarphiHVarphi(:, :, -1:)

    integer                                :: lmp
    integer                                :: oqn_l
    integer                                :: mqn_m
    integer                                :: lm
    integer                                :: ichan
    integer                                :: oqn_l3Pr
    integer                                :: mqn_m3Pr
    integer                                :: iradf
    integer                                :: imesh
    real                                   :: integralR
    real                                   :: integralI
    integer                                :: lmp1Pr
    integer                                :: oqn_l1Pr
    integer                                :: mqn_m1Pr
    integer                                :: lm1Pr
    integer                                :: iradf1Pr
    integer                                :: lm3Pr
    integer                                :: mqn_m2Pr

    real,          allocatable             :: intgrdR(:)
    real,          allocatable             :: intgrdI(:)

    ! oqn_l3Pr is only a temporary variable

    allocate( intgrdR(atoms%jri(itype)), intgrdI(atoms%jri(itype)) )
    !todo check whether this setting of zero is okay everywhere and wheter we should intent out or intent inout
    grVarphiVarphi = cmplx(0., 0.)
    grVarphiHVarphi(:, :, :) = cmplx(0., 0.)
    do mqn_m2Pr = -1, 1
      lmp1Pr = 0
      do oqn_l1Pr = 0, atoms%lmax(itype)
        do mqn_m1Pr = -oqn_l1Pr, oqn_l1Pr
          mqn_m = grVarphiChMout(mqn_m1Pr, mqn_m2Pr)
          do iradf1Pr = 1, nRadFun(oqn_l1Pr, itype)
            lmp1Pr = lmp1Pr + 1

            ! Actually we should loop over lmp1Pr but for better performance we place the loop here
            do ichan = 1, chanMax
              oqn_l = grVarphiChLout(ichan, oqn_l1Pr)
              if ( ( abs(mqn_m) > oqn_l ) .or. ( oqn_l < 0 ) .or. ( oqn_l > atoms%lmax(itype) ) ) cycle
              lmp = 0
              ! Determine lmp1Pr by counting until scattering channel
              ! Approach oqn_l1Pr
              do oqn_l3Pr = 0, oqn_l - 1
                do mqn_m3Pr = -oqn_l3Pr, oqn_l3Pr
                  do iradf = 1, nRadFun(oqn_l3Pr, itype)
                    lmp = lmp + 1
                  end do ! iradf1
                end do ! mqn_m3Pr
              end do ! oqn_l3Pr
              ! Approach mqn_m1Pr
              do mqn_m3Pr = -oqn_l, mqn_m - 1
                do iradf = 1, nRadFun(oqn_l, itype)
                  lmp = lmp + 1
                end do ! iradf
              end do ! mqn_m3Pr
              do iradf = 1, nRadFun(oqn_l, itype)
                lmp = lmp + 1
                intgrdR(:) = 0.
                do imesh = 1, atoms%jri(itype)
                  intgrdR(imesh) = r2(imesh) * ( grVarPhiCh1(imesh, ichan, lmp1Pr, mqn_m2Pr) * varphi1(imesh, iradf, oqn_l)       &
                    & + grVarPhiCh2(imesh, ichan, lmp1Pr, mqn_m2Pr) * varphi2(imesh, iradf, oqn_l) )
                end do ! imesh
                !call intgr3NoIntp(intgrdR(:), atoms%rmsh(1, itype), atoms%dx(itype), atoms%jri(itype), integralR)
                call intgr3(intgrdR(:), atoms%rmsh(1, itype), atoms%dx(itype), atoms%jri(itype), integralR)! TODO: Is this ok?
                !call intgr3(intgrdR(:), atoms%rmsh(1, itype), atoms%dx(itype), atoms%jri(itype), integralR)
                ! < grVarphi | varphi >
                grVarphiVarphi(lmp1Pr, lmp, mqn_m2Pr) = grVarphiVarphi(lmp1Pr, lmp, mqn_m2Pr) + integralR
              end do ! iradf
            end do ! ichan

            lmp = 0
            lm = -1
            mqn_m3Pr = grVarphiChMout(mqn_m1Pr, mqn_m2Pr)
            do oqn_l = 0, atoms%lmax(itype)
              do mqn_m = -oqn_l, oqn_l
                lm = lm + 1
                do iradf = 1, nRadFun(oqn_l1Pr, itype)
                  lmp = lmp + 1
                  do ichan = 1, chanMax
                    oqn_l3Pr = grVarphiChLout(ichan, oqn_l1Pr)
                    ! We have to catch channels where l < 0 and where m"' does not contribute, lmax + 1 is allowed due to the gradient
                    if ( ( abs(mqn_m3Pr) > oqn_l3Pr )  .or. ( oqn_l3Pr < 0 )  )  cycle
                    lm3Pr = oqn_l3Pr * (oqn_l3Pr + 1) + mqn_m3Pr
                    intgrdR(:) = 0.
                    intgrdI(:) = 0.
                    do imesh = 1, atoms%jri(itype)
                      intgrdR(imesh) = r2(imesh) * real ( &
                        &   (grVarPhiCh1(imesh, ichan, lmp1Pr, mqn_m2Pr) * hVarphi(1, imesh, lm3Pr, lmp) &
                        &  + grVarPhiCh2(imesh, ichan, lmp1Pr, mqn_m2Pr) * hVarphi(2, imesh, lm3Pr, lmp) ) )
                      intgrdI(imesh) = r2(imesh) * aimag( &
                        &   (grVarPhiCh1(imesh, ichan, lmp1Pr, mqn_m2Pr) * hVarphi(1, imesh, lm3Pr, lmp) &
                        &  + grVarPhiCh2(imesh, ichan, lmp1Pr, mqn_m2Pr) * hVarphi(2, imesh, lm3Pr, lmp) ) )
                    end do ! imesh
            !        call intgr3NoIntp(intgrdR(:), atoms%rmsh(1, itype), atoms%dx(itype), atoms%jri(itype), integralR)
                call intgr3(intgrdR(:), atoms%rmsh(1, itype), atoms%dx(itype), atoms%jri(itype), integralR)! TODO: Is this ok?
            !        call intgr3NoIntp(intgrdI(:), atoms%rmsh(1, itype), atoms%dx(itype), atoms%jri(itype), integralI)
                call intgr3(intgrdI(:), atoms%rmsh(1, itype), atoms%dx(itype), atoms%jri(itype), integralI)! TODO: Is this ok?
                    ! < grVarphi| H |varphi >
                    grVarphiHvarphi(lmp1Pr, lmp, mqn_m2Pr) = grVarphiHvarphi(lmp1Pr, lmp, mqn_m2Pr) &
                      & + cmplx(integralR, integralI)
                  end do ! ichan
                end do ! iradf1Pr
              end do ! mqn_m1Pr
            end do ! oqn_l1Pr
          end do ! p
        end do ! mqn_m
      end do ! oqn_l
    end do ! mqn_m1Pr

  end subroutine CalcVecBasfMatElems

  subroutine CalcScalBasfMatElems( atoms, itype, iatom, nRadFun, r2, varphi1, varphi2, hVarphi, varphiVarphi, varphiHvarphi, usdus )

    use m_types, only : t_atoms, t_usdus
 !   use m_intgr, only : intgr3NoIntp
    use m_intgr, only : intgr3!LinIntp! TODO: Is this ok?
 !   use m_intgr, only : intgr3

    implicit none

    ! Type parameters
    type(t_atoms),              intent(in)  :: atoms
    type(t_usdus), optional,    intent(in) :: usdus

    ! Scalar parameters
    integer,                    intent(in)  :: itype
    integer,                    intent(in)  :: iatom

    ! Array parameters
    integer,                    intent(in)  :: nRadFun(0:, :)
    real,                       intent(in)  :: r2(:)
    real,                       intent(in)  :: varphi1(:, :, 0:)
    real,                       intent(in)  :: varphi2(:, :, 0:)
    complex,                    intent(in)  :: hVarphi(:, :, 0:, :)
    real,                       intent(inout) :: varphiVarphi(:, :, :)
    complex,                    intent(inout) :: varphiHvarphi(:, :, :)

    ! Scalar variabels
    integer                                 :: lmp
    integer                                 :: oqn_l
    integer                                 :: mqn_m
    integer                                 :: lm
    integer                                 :: iradf
    integer                                 :: iradf1Pr
    integer                                 :: imesh
    real                                    :: integralR
    real                                    :: integralI
    integer                                 :: lmp1Pr
    integer                                 :: oqn_l1Pr
    integer                                 :: mqn_m1Pr
    integer                                 :: lm1Pr

    ! Array variabels
    real,          allocatable              :: intgrdR(:)
    real,          allocatable              :: intgrdI(:)

    allocate( intgrdR(atoms%jri(itype)), intgrdI(atoms%jri(itype)) )
    lmp = 0
    lm = -1
    do oqn_l = 0, atoms%lmax(itype)
      do mqn_m = -oqn_l, oqn_l
        lm = lm + 1
        ! todo we have to check whether this is consistent with the integrations routines
        ! We place here the filling of the overlap matrix because it is diagonal in lm. For performance reasons we might reduce
        ! readability of the code.
        !todo beware not iatom but itype in varphivarphi!!!!!!
        ! todo delete this block if radial integration is not the same with intgr0 and intgr3
        !varphiVarphi(1, 1, lm, iatom) = 1
        !do p1Pr = 3, nRadFun(oqn_l, itype)
        !  varphiVarphi(p1Pr, 1, lm, iatom) = usdus%uulon(oqn_l, itype, 1)
        !end do ! p1Pr
        !varphiVarphi(2, 2, lm, iatom) =  usdus%ddn(oqn_l, itype, 1)
        !do p1Pr = 3, nRadFun(oqn_l, itype)
        !  varphiVarphi(p1Pr,2, lm, iatom) = usdus%dulon(oqn_l, itype, 1)
        !end do ! p1Pr
        !do p = 3, nRadFun(oqn_l, itype)
        !  do p1Pr = 3, nRadFun(oqn_l, itype)
        !    varphiVarphi(p1Pr, p, lm, iatom) = usdus%uloulopn(iloTable(p1Pr, oqn_l, itype), iloTable(p, oqn_l, itype), oqn_l, itype)
        !  end do ! p1Pr
        !end do ! p

        ! Calculation of <varphi|varphi>
        do iradf = 1, nRadFun(oqn_l, itype)
          do iradf1Pr = 1, nRadFun(oqn_l, itype)
            intgrdR(:) = 0.
            do imesh = 1, atoms%jri(itype)
              intgrdR(imesh) = r2(imesh) * (varphi1(imesh, iradf1Pr, oqn_l) * varphi1(imesh, iradf, oqn_l) &
                & + varphi2(imesh, iradf1Pr, oqn_l) * varphi2(imesh, iradf, oqn_l))
            end do ! imesh
            call intgr3(intgrdR(:), atoms%rmsh(1, itype), atoms%dx(itype), atoms%jri(itype), integralR )! TODO: Is this ok?
            varphiVarphi(lmp + iradf1Pr, lmp + iradf, itype) = integralR
          end do ! iradf1Pr
        end do ! iradf
        ! Calculation of <varphi|H|varphi>
        do iradf = 1, nRadFun(oqn_l, itype)
          lmp = lmp + 1
          lmp1Pr = 0
          lm1Pr = -1
          do oqn_l1Pr = 0, atoms%lmax(itype)
            do mqn_m1Pr = -oqn_l1Pr, oqn_l1Pr
              lm1Pr = lm1Pr + 1
              do iradf1Pr = 1, nRadFun(oqn_l1Pr, itype)
                lmp1Pr = lmp1Pr + 1
                intgrdR(:) = 0.
                intgrdI(:) = 0.
                do imesh = 1, atoms%jri(itype)
                  ! The hVarphi arrays contain already the Gaunt coefficients and are dependent on the lm1Pr the bra varphi can be
                  ! matched to.
                  intgrdR(imesh) = r2(imesh) * real((varphi1(imesh, iradf1Pr, oqn_l1Pr) * hVarphi(1, imesh, lm1Pr, lmp) &
                                  & + varphi2(imesh, iradf1Pr, oqn_l1Pr) * hVarphi(2, imesh, lm1Pr, lmp)))
                  intgrdI(imesh) = r2(imesh) * aimag((varphi1(imesh, iradf1Pr, oqn_l1Pr) * hVarphi(1, imesh, lm1Pr, lmp) &
                                  & + varphi2(imesh, iradf1Pr, oqn_l1Pr) * hVarphi(2, imesh, lm1Pr, lmp)))
                end do ! imesh
!                call intgr3(intgrdR(:), atoms%rmsh(1, itype), atoms%dx(itype), atoms%jri(itype), integralR)
!                call intgr3(intgrdI(:), atoms%rmsh(1, itype), atoms%dx(itype), atoms%jri(itype), integralI)
!                call intgr3NoIntp(intgrdR(:), atoms%rmsh(1, itype), atoms%dx(itype), atoms%jri(itype), integralR)
!                call intgr3NoIntp(intgrdI(:), atoms%rmsh(1, itype), atoms%dx(itype), atoms%jri(itype), integralI)
                call intgr3(intgrdR(:), atoms%rmsh(1, itype), atoms%dx(itype), atoms%jri(itype), integralR )! TODO: Is this ok?
                call intgr3(intgrdI(:), atoms%rmsh(1, itype), atoms%dx(itype), atoms%jri(itype), integralI )! TODO: Is this ok?
                varphiHvarphi(lmp1Pr, lmp, iatom) = cmplx(integralR, integralI)

                !!!For symmetrized kinetic energy:
                if (present(usdus)) then
                  if ((oqn_l.eq.oqn_l1Pr).and.(mqn_m.eq.mqn_m1Pr)) then
                    if (iradf1Pr.eq.1) then
                      if (iradf.eq.1) then
                        varphiHvarphi(lmp1Pr, lmp, iatom) = varphiHvarphi(lmp1Pr, lmp, iatom) &
                                                        & + 0.5*usdus%us(oqn_l, itype, 1)*usdus%dus(oqn_l, itype, 1)*atoms%rmt(itype)**2
                      end if
                      if (iradf.eq.2) then
                        varphiHvarphi(lmp1Pr, lmp, iatom) = varphiHvarphi(lmp1Pr, lmp, iatom) &
                                                        & + 0.5*usdus%us(oqn_l, itype, 1)*usdus%duds(oqn_l, itype, 1)*atoms%rmt(itype)**2
                      end if
                    end if
                    if (iradf1Pr.eq.2) then
                      if (iradf.eq.1) then
                        varphiHvarphi(lmp1Pr, lmp, iatom) = varphiHvarphi(lmp1Pr, lmp, iatom) &
                                                        & + 0.5*usdus%uds(oqn_l, itype, 1)*usdus%dus(oqn_l, itype, 1)*atoms%rmt(itype)**2
                      end if
                      if (iradf.eq.2) then
                        varphiHvarphi(lmp1Pr, lmp, iatom) = varphiHvarphi(lmp1Pr, lmp, iatom) &
                                                        & + 0.5*usdus%uds(oqn_l, itype, 1)*usdus%duds(oqn_l, itype, 1)*atoms%rmt(itype)**2
                      end if
                    end if
                  end if
                end if
              end do ! iradf1Pr
            end do ! mqn_m1Pr
          end do ! oqn_l1Pr
        end do ! iradf
      end do ! mqn_m
    end do ! oqn_l

  end subroutine CalcScalBasfMatElems

  !>-------------------------------------------------------------------------------------------------------------------------------------
  !> @author
  !> Christian-Roman Gerhorst
  !>
  !> @brief
  !> Implementation of A.52 PhD thesis Aaron Klueppelberg (PhDAK). There lacks an index p at the radial solutions u.
  !>
  !> @details
  !> If we want to express the second braket by the first braket in 7.718 PhDAK and use the self-adjointness of the Hamiltonian, we
  !> have to correct the kinetic energy by A.52 PhDAK. This routine implements the correction starting with 0.5 in the last line of
  !> A.52. PhDAK.
  !>
  !> @note
  !> This routine can be used for the correction of matrix elements containing one or two gradients if the maximal number of
  !> channels chanMax is correctly set. Remember that the input functions f have to be of the structure which comes out of the
  !> routine CalcChannelsGrFlpNat
  !>
  !> @todo
  !> Account for the spin.
  !>-------------------------------------------------------------------------------------------------------------------------------------
  subroutine CalcSelfAdjCorrection(atoms, chanMax, nRadFunMax, itype, mqn_m2Pr, nRadFun, grFlpChLout, grFlpChMout,   &
      & varphi1, varphi2, delrVarphi1, delrVarphi2, braFuncCh1, braFuncCh2, delrBraFuncCh1, delrBraFuncCh2, selfAdjCorrection)

      use m_types, only : t_atoms
      use m_juDFT_stop, only : juDFT_error

      implicit none

      ! Type paramter
      type(t_atoms),                  intent(in)  :: atoms

      ! Scalar parameter
      ! chanMax: maximal number of scattering channels
      integer,                        intent(in)  :: chanMax
      integer,                        intent(in)  :: nRadFunMax
      integer,                        intent(in)  :: itype
      integer,                        intent(in)  :: mqn_m2Pr


      ! Array parameter
      ! nradfun : number of radial solutions
      ! grVarphiChLout : resulting l channels within scattering channels of input functions
      ! varphi1 : radial basis functions without r
      ! varphi2 : radial basis function without r
      ! delrVarphi1 : radial derivative of large radial basis function divided by r
      ! delrVarphi2 : radial derivative of small radial basis function divided by r
      ! braFunc1 : large component of function f in A.52 PhDAK
      ! braFunc2 : small component of function f in A.52 PhDAK
      ! delrBraFunc1 : large component of radial derivative of f in A.52 PhDAK
      ! delrBraFunc2 : small component of radial derivative of f in A.52 PhDAK
      ! selfAdjCorrection : result of the routine
      integer,                        intent(in)  :: nRadFun(0:, :)
      real,                           intent(in)  :: varphi1(:, :, 0:)
      real,                           intent(in)  :: varphi2(:, :, 0:)
      real,                           intent(in)  :: delrVarphi1(:, :, 0:)
      real,                           intent(in)  :: delrVarphi2(:, :, 0:)
      real,                           intent(in)  :: braFuncCh1(:, :, :, -1:)
      real,                           intent(in)  :: braFuncCh2(:, :, :, -1:)
      real,                           intent(in)  :: delrBraFuncCh1(:, :, :, -1:)
      real,                           intent(in)  :: delrBraFuncCh2(:, :, :, -1:)
      integer,                        intent(in)  :: grFlpChLout(:, 0:)
      integer,                        intent(in)  :: grFlpChMout(-atoms%lmaxd:, -1:)
      real,                         intent(inout) :: selfAdjCorrection(:, :)

      ! Scalar variable
      integer                                     :: oqn_l1Pr
      integer                                     :: mqn_m1Pr
      integer                                     :: iradf1Pr
      integer                                     :: oqn_l
      integer                                     :: mqn_m
      integer                                     :: oqn_l3Pr
      integer                                     :: mqn_m3Pr
      integer                                     :: ichan
      integer                                     :: iradf
      integer                                     :: lmp
      integer                                     :: lmp1Pr
      integer                                     :: imeshBound

      imeshBound = atoms%jri(itype)


      ! Add up scattering channel of the gradient to enable calculation of the correction term
      lmp = 0
      do oqn_l = 0, atoms%lmax(itype)
        do mqn_m = -oqn_l, oqn_l
          mqn_m1Pr = grFlpChMout(mqn_m, mqn_m2Pr)
          do iradf = 1, nRadFun(oqn_l, itype)
            lmp = lmp + 1
            do ichan = 1, chanMax
              oqn_l1Pr = grFlpChLout(ichan, oqn_l)
              if ( ( oqn_l1Pr < 0 ) .or. ( oqn_l1Pr > atoms%lmax(itype) ) .or. abs(mqn_m1Pr) > oqn_l1Pr ) cycle
              lmp1Pr = 0
              ! Approach oqn_l1Pr
              do oqn_l3Pr = 0, oqn_l1Pr - 1
                do mqn_m3Pr = -oqn_l3Pr, oqn_l3Pr
                  do iradf1Pr = 1, nRadFun(oqn_l3Pr, itype)
                    lmp1Pr = lmp1Pr + 1
                  end do ! iradf1Pr
                end do ! mqn_m3Pr
              end do ! oqn_l3Pr
              ! Approach mqn_m1Pr
              do mqn_m3Pr = -oqn_l1Pr, mqn_m1Pr - 1
                do iradf1Pr = 1, nRadFun(oqn_l1Pr, itype)
                  lmp1Pr = lmp1Pr + 1
                end do ! iradf1Pr
              end do ! mqn_m3Pr
              ! Calculate at oqn_l1Pr and mqn_m1Pr for all radial functions
              do iradf1Pr = 1, nRadFun(oqn_l1Pr, itype)
                lmp1Pr = lmp1Pr + 1
                selfAdjCorrection(lmp1Pr, lmp) = selfAdjCorrection(lmp1Pr, lmp) + 0.5 *        &
                  & ( ( delrVarphi1(imeshBound, iradf1Pr, oqn_l1Pr) * braFuncCh1(imeshBound, ichan, lmp, mqn_m2Pr)                 &
                  &   + delrVarphi2(imeshBound, iradf1Pr, oqn_l1Pr) * braFuncCh2(imeshBound, ichan, lmp, mqn_m2Pr) )               &
                  & - ( varphi1(imeshBound, iradf1Pr, oqn_l1Pr)     * delrBraFuncCh1( imeshBound, ichan, lmp, mqn_m2Pr)            &
                      + varphi2(imeshBound, iradf1Pr, oqn_l1Pr)     * delrBraFuncCh2( imeshBound, ichan, lmp, mqn_m2Pr) ) )
              end do ! iradf1Pr
            end do ! ichan
          end do ! iradf
        end do ! mqn_m
      end do ! oqn_l

  end subroutine CalcSelfAdjCorrection

  !>-------------------------------------------------------------------------------------------------------------------------------------
  !> @author
  !> Christian-Roman Gerhorst
  !>
  !> @brief
  !> Implementation of A.52 PhD thesis Aaron Klueppelberg (PhDAK). There lacks an index p at the radial solutions u.
  !>
  !> @details
  !> If we want to express the second braket by the first braket in 7.718 PhDAK and use the self-adjointness of the Hamiltonian, we
  !> have to correct the kinetic energy by A.52 PhDAK. This routine implements the correction starting with 0.5 in the last line of
  !> A.52. PhDAK.
  !>
  !> @note
  !> This routine can be used for the correction of matrix elements containing one or two gradients if the maximal number of
  !> channels chanMax is correctly set. Remember that the input functions f have to be of the structure which comes out of the
  !> routine CalcChannelsGrFlpNat
  !>
  !> @todo
  !> Account for the spin.
  !>-------------------------------------------------------------------------------------------------------------------------------------
  subroutine CalcSelfAdjCorrPhi( atoms, itype, nRadFun, varphi1, varphi2, delrVarphi1, delrVarphi2, selfAdjCorVarphi )

      use m_types, only : t_atoms

      implicit none

      ! Type parameters
      type(t_atoms), intent(in)  :: atoms

      ! Scalar parameters
      integer,       intent(in)  :: itype

      ! Array parameters
      ! nradfun : number of radial solutions
      ! grVarphiChLout : resulting l channels within scattering channels of input functions
      ! varphi1 : radial basis functions without r
      ! varphi2 : radial basis function without r
      ! delrVarphi1 : radial derivative of large radial basis function divided by r
      ! delrVarphi2 : radial derivative of small radial basis function divided by r
      ! braFunc1 : large component of function f in A.52 PhDAK
      ! braFunc2 : small component of function f in A.52 PhDAK
      ! delrBraFunc1 : large component of radial derivative of f in A.52 PhDAK
      ! delrBraFunc2 : small component of radial derivative of f in A.52 PhDAK
      ! selfAdjCorrection : result of the routine
      integer,       intent(in)  :: nRadFun(0:, :)
      real,          intent(in)  :: varphi1(:, :, 0:)
      real,          intent(in)  :: varphi2(:, :, 0:)
      real,          intent(in)  :: delrVarphi1(:, :, 0:)
      real,          intent(in)  :: delrVarphi2(:, :, 0:)
      real,          intent(out) :: selfAdjCorVarphi(:, :)

      integer                    :: oqn_l1Pr
      integer                    :: lm_pre1Pr
      integer                    :: mqn_m1Pr
      integer                    :: iradf
      integer                    :: iradf1Pr
      integer                    :: lm1Pr
      integer                    :: lmp1Pr

      selfAdjCorVarphi(:, :) = cmplx(0., 0.)
      lmp1Pr = 0
      do oqn_l1Pr = 0, atoms%lmax(itype)
        lm_pre1Pr = oqn_l1Pr * (oqn_l1Pr + 1) + 1
        do mqn_m1Pr = -oqn_l1Pr, oqn_l1Pr
          lm1Pr = oqn_l1Pr * (oqn_l1Pr + 1) + mqn_m1Pr
          do iradf = 1, nRadFun(oqn_l1Pr, itype)
            do iradf1Pr = 1, nRadFun(oqn_l1Pr, itype)
              selfAdjCorVarphi(lmp1Pr + iradf, lmp1Pr + iradf1Pr) = selfAdjCorVarphi(lmp1Pr + iradf, lmp1Pr + iradf1Pr) + 2 *    &
                & ( delrVarphi1(atoms%jri(itype), iradf1Pr, oqn_l1Pr) * varphi1(atoms%jri(itype), iradf, oqn_l1Pr)                 &
                & + delrVarphi2(atoms%jri(itype), iradf1Pr, oqn_l1Pr) * varphi2(atoms%jri(itype), iradf, oqn_l1Pr)                 &
                & - varphi1(atoms%jri(itype), iradf1Pr, oqn_l1Pr)     * delrVarphi1(atoms%jri(itype), iradf, oqn_l1Pr)             &
                & - varphi2(atoms%jri(itype), iradf1Pr, oqn_l1Pr)     * delrVarphi2(atoms%jri(itype), iradf, oqn_l1Pr) )
            end do ! iradf
          end do ! iradf1Pr
          lmp1Pr = lmp1Pr + nRadFun(oqn_l1Pr, itype)
        end do ! mqn_m1Pr
      end do ! oqn_l1Pr

  end subroutine CalcSelfAdjCorrPhi

  !>-------------------------------------------------------------------------------------------------------------------------------------
  !> @author
  !> Christian-Roman Gerhorst
  !>
  !> @brief
  !> Implementation of A.50; A.51 is implemented recyling A.50 and using the correction A.52 in PhD thesis A. Klueppelberg
  !>
  !> @details
  !> Within this routine the various scalar-, vector- and matrix-like large matching coefficients are determined and multiplied to the
  !> matrix elements of the basis functions given as input parameters. To evaluate A.51 we complex conjugate A.50 but have to correct the
  !> error that the kinetic energy of the Hamiltonian is not self-adjoint in LAPW if we deal with gradient or two-fold gradients onto the
  !> basis functions. Both equations, A.50 and A.51, are added to the dynMatPu array. We only deal with the MT here so we do not have to
  !> account for the back-folding vector.
  !>
  !> @note
  !> The filling for the z1 is not optimal in memory access, but minimizes array sizes, and ensures a linear run through the array later.
  !>
  !> @todo
  !> Account for the spin.
  !> Some arrays may be saved in the Heps calculation. we do not need one array for H and one for the overlap we can write into the
  !> heps array at once!
  !>-------------------------------------------------------------------------------------------------------------------------------------
  subroutine Add2ndOrdWfPulayBraKets2DynMat( fmpi, noco, nococonv, oneD, atoms, input, cell, sym, kpts, qpts, usdus, results, iqpt, ikpt, pMax, lmpMax, eps1DynMatPulTestSw, testComp2ndN1stOrdBasFuncSw, mapKpq2K, eig, nobd, nv, &
      & gbas, ilst, kveclo, z, z1nG, nRadFun, iloTable, varphiVarphi, varphiHvarphi, grVarphiVarphiNat, grVarphiHVarphiNat, &
      & lmpT, varphiGrVeff0SphVarphi, varphiHGrvarphiNat, varphiGrVarphiNat, varphiGrVeff0Varphi, dynMatPu )

    use m_abcof3
    use m_types
    use m_jp2ndOrdQuant, only : outerProduct, outerProductME
    use m_juDFT_stop, only : juDFT_error

    implicit none

    ! Type parameters
    type(t_mpi),                  intent(in)  :: fmpi
    type(t_nococonv),                intent(in)  :: nococonv
    type(t_noco),                intent(in)  :: noco
    type(t_oneD),                intent(in)  :: oneD
    type(t_atoms),                  intent(in)    :: atoms
    type(t_input),                  intent(in)    :: input
    type(t_cell),                   intent(in)    :: cell
    type(t_sym),                    intent(in)    :: sym
    type(t_kpts),                   intent(in)    :: kpts
    type(t_kpts),                   intent(in)    :: qpts
    type(t_usdus),                  intent(in)    :: usdus
    type(t_results),                intent(in)    :: results

    ! Scalar parameters
    integer,                        intent(in)    :: iqpt
    integer,                        intent(in)    :: ikpt
    integer,                        intent(in)    :: pMax
    integer,                        intent(in)    :: lmpMax
    logical,                        intent(in)    :: eps1DynMatPulTestSw
    logical,                        intent(in)    :: testComp2ndN1stOrdBasFuncSw

    ! Array parameters
    integer,                        intent(in)    :: lmpT(:)
    integer,                        intent(in)    :: mapKpq2K(:, :)
    real,                           intent(in)    :: eig(:, :, :)
    integer,                        intent(in)    :: nobd(:, :)
    integer,                        intent(in)    :: nv(:, :)
    integer,                        intent(in)    :: gbas(:, :)
    integer,                        intent(in)    :: ilst(:, :, :)
    integer,                        intent(in)    :: kveclo(:,:)
    complex,                       intent(in)    :: z(:,:,:,:)
    complex,                        intent(in)    :: z1nG(:, :, :, :)
    integer,                        intent(in)    :: nRadFun(0:, :)
    integer,                        intent(in)    :: iloTable(:, 0:, :)
    real,                           intent(in)    :: varphiVarphi(:, :, :)
    complex,                        intent(in)    :: varphiHvarphi(:, :, :)
    real,                           intent(in)    :: grVarphiVarphiNat(:, :, :, :)
    complex,                        intent(in)    :: grVarphiHVarphiNat(:, :, :, :)
    complex,                        intent(in)    :: varphiGrVeff0SphVarphi(:, :, :, :)
    complex,                        intent(in)    :: varphiGrVeff0Varphi(:, :, :, :)
    complex,                        intent(in)    :: varphiHGrvarphiNat(:, :, :, :)
    real,                           intent(in)    :: varphiGrVarphiNat(:, :, :, :)
    complex,                        intent(inout) :: dynMatPu(:, :)


    ! Type variable
    type(od_inp)                                  :: odi
    type(od_sym)                                  :: ods
    type(t_lapw) :: lapw

    ! Scalar variables
    integer                                       :: iBas
    integer                                       :: nmat
    integer                                       :: iDtypeA
    integer                                       :: iDeqatA
    integer                                       :: iDatomA
    integer                                       :: iband
    integer                                       :: lmp
    integer                                       :: lm
    integer                                       :: oqn_l
    integer                                       :: mqn_m
    integer                                       :: pMaxLocal
    integer                                       :: iradf
    integer                                       :: iDtypeB
    integer                                       :: iDeqatB
    integer                                       :: iDatomB
    integer                                       :: idirR
    integer                                       :: idirC
    integer                                       :: ikpq
    complex                                       :: matElemHeps
    integer :: nk



    ! Array variables
    integer,           allocatable                :: ngoprI(:)
    complex,           allocatable                :: a(:, :, :)
    complex,           allocatable                :: b(:, :, :)
    complex,           allocatable                :: bascof_lo(:, :, :, :, :)
    complex,           allocatable                :: aKpq(:, :, :)
    complex,           allocatable                :: bKpq(:, :, :)
    complex,           allocatable                :: bascof_loKpq(:, :, :, :, :)
    complex,           allocatable                :: z1Gext(:)
    real,              allocatable                :: gbasExt(:, :)
    real,              allocatable                :: gbasExtKpq(:, :)
    real,              allocatable                :: gBasMatExt(:, :, :)
    complex,           allocatable                :: ab0cofScl(:)
    complex,           allocatable                :: varphiPsi(:)
    complex,           allocatable                :: psiVarphi(:)
    complex,           allocatable                :: psiHVarphi(:)
    complex,           allocatable                :: ab0cofVec(:, :)
    complex,           allocatable                :: abcofSumVec(:, :)
    complex,           allocatable                :: grVarphiPsiNat(:, :)
    complex,           allocatable                :: varphiGrPsiNat(:)
    complex,           allocatable                :: ab1cofVec(:, :, :)
    complex,           allocatable                :: abcofMat(:)
    complex,           allocatable                :: abcofSumMat(:, :, :, :)
    complex,           allocatable                :: abcofSumMatAlpha(:)
    complex,           allocatable                :: varphiHPsi(:)
    complex,           allocatable                :: grVarphiHPsiNat(:, :)
    complex,           allocatable                :: varphiHGrPsiNat(:)
    complex,           allocatable                :: varphiGrVeff0SphPsi(:)
    complex,           allocatable                :: varphiGrVeff0Psi(:)
    complex                                       :: psiHepsPsi(3, 3)
    complex                                       :: psiGrPsiNat(3, 3)
    complex                                       :: grPsiPsiNat(3, 3)
    complex                                       :: grPsiHPsiNat(3, 3)
    complex                                       :: psiHGrPsiNat(3, 3)
    complex                                       :: psiGrVeff0sphPsi(3, 3)
    complex                                       :: psiGrVeff0Psi(3, 3)
    complex                                       :: grPsiHepsPsiNat(3, 3)
    complex                                       :: psiHepsGrPsiNat(3, 3)
    complex                                       :: grPsiHepsPsi(3, 3)
    complex                                       :: psiHepsGrPsi(3, 3)
    complex                                       :: grPsiHPsiNatT(3, 3)
    complex                                       :: psiHGrPsiNatT(3, 3)
    complex                                       :: psiGrVeff0sphPsiT(3, 3)
    complex                                       :: psiGrVeff0PsiT(3, 3)
    complex                                       :: grPsiPsiNatT(3, 3)
    complex                                       :: psiGrPsiNatT(3, 3)
    complex                                       :: grPsiHepsPsiNatT(3, 3)
    complex                                       :: psiHepsGrPsiNatT(3, 3)
    complex                                       :: grPsiHepsPsiT(3, 3)
    complex                                       :: psiHepsGrPsiT(3, 3)
    real                                          :: kExt(3)
    real                                          :: kExtkExtT(3, 3)
    real                                          :: kpqExt(3)

    ! We do not want the local coordinate systems to be rotated for non-representative atoms constructing the matching coefficients.
    allocate( ngoprI(atoms%nat) )
    ngoprI(:) = 1

    ! Allocate matching coefficients of MT basis functions
    allocate( a( MAXVAL(nv), 0:atoms%lmaxd*(atoms%lmaxd+2), atoms%nat), b(MAXVAL(nv), 0:atoms%lmaxd*(atoms%lmaxd+2), atoms%nat), &
      & bascof_lo(3, -atoms%llod:atoms%llod, 4 * atoms%llod + 2, atoms%nlod, atoms%nat) )
    allocate( aKpq( MAXVAL(nv), 0:atoms%lmaxd*(atoms%lmaxd+2), atoms%nat), bKpq(MAXVAL(nv), 0:atoms%lmaxd*(atoms%lmaxd+2), atoms%nat), &
      & bascof_loKpq(3, -atoms%llod:atoms%llod, 4 * atoms%llod + 2, atoms%nlod, atoms%nat) )
    allocate( z1Gext(MAXVAL(nv)) )
    allocate( gbasExt( SIZE(z(:,1,1,1)), 3), gbasExtKpq( SIZE(z(:,1,1,1)), 3), gBasMatExt( SIZE(z(:,1,1,1)), 3, 3 ) )
    allocate( ab0cofScl(lmpMax), ab0cofVec(lmpMax, 3), abcofSumVec(lmpMax, 3), ab1cofVec( lmpMax, atoms%nat, 3), abcofMat(pMax), &
      & abcofSumMat(lmpMax, atoms%nat, 3, 3), abcofSumMatAlpha(lmpMax) )
    allocate( varphiPsi(lmpMax), grVarphiPsiNat(lmpMax, 3), varphiHPsi(lmpMax), &
      & grVarphiHPsiNat(lmpMax, 3) )
    allocate( varphiGrPsiNat(lmpMax), varphiHGrPsiNat(lmpMax), varphiGrVeff0sphPsi(lmpMax) )
    allocate( varphiGrVeff0Psi(lmpMax) )
    allocate( psiVarphi(lmpMax), psiHVarphi(lmpMax) )

    ! This loop can be parallelized later

      ikpq = mapKpq2K(ikpt, iqpt)

      ! We have terms with factors of k, k + q or G, where a transformation into cartesian coordinates is required.
      kExt(1:3) = matmul(cell%bmat, kpts%bk(1:3, ikpt))
      kExtkExtT(1:3, 1:3) = outerProduct(kExt(1:3), kExt(1:3))
      kpqExt(1:3) = matmul(cell%bmat, kpts%bk(1:3, ikpq))

      gbasExt(:, :) = 0.
      ! We assign the G-basis vectors once so later we can just loop without jumps over iBas.
      do iBas = 1, nv(1, ikpt)
        gbasExt(iBas, 1:3) = matmul( cell%bmat(1:3, 1:3), gbas(1:3, ilst(iBas, ikpt, 1)))
      end do ! iBas
      do iBas = nv(1, ikpt) + 1, nv(1, ikpt) + atoms%nlotot
        gBasExt(iBas, 1:3) = matmul( cell%bmat(1:3, 1:3), gbas(1:3, ilst(kveclo(iBas, ikpt), ikpt, 1)) )
      end do ! iBas

      gbasExtKpq(:, :) = 0.
      ! The nv has to be taken at ikpq. We assign the G-basis vectors once so later we can just loop without jumps over iBas.
      do iBas = 1, nv(1, ikpq)
        gbasExtKpq(iBas, 1:3) = matmul( cell%bmat, gbas(1:3, ilst(iBas, ikpq, 1)) )
      end do ! iBas
      do iBas = nv(1, ikpq) + 1, nv(1, ikpq) + atoms%nlotot
        gbasExtKpq(iBas, 1:3) = matmul( cell%bmat, gbas(1:3, ilst(kveclo(iBas, ikpq), ikpq, 1)) )
      end do ! iBas

      gBasMatExt(:, :, :) = 0.
      do idirC = 1, 3
        do idirR = 1, 3
          do iBas = 1, nv(1, ikpt) + atoms%nlotot
            ! The dyadic product is rewritten by a proper choice of idir loop variables
            gbasMatExt(iBas, idirR, idirC) = gbasExt(iBas, idirR) * gbasExt(iBas, idirC)
          end do ! iBas
        end do ! idirR
      end do ! idirC

      ! Calculate the basis matching coefficients at k + q = k' and the matching coefficients at k.
      nmat = nv(1, ikpq) + atoms%nlotot
      aKpq(:, :, :) = cmplx(0.0, 0.0)
      bKpq(:, :, :) = cmplx(0.0, 0.0)
      bascof_loKpq(:, :, :, :, :) = cmplx(0.0, 0.0)
      !call abcof3( atoms%lmaxd, atoms%ntype, atoms%nat, sym%nop, MAXVAL(nv), input%jspins, 1, atoms%lmaxd*(atoms%lmaxd+2), SIZE(z(:,1,1,1)), &
        !& atoms%llod, atoms%nlod, atoms%nlotot, sym%invtab, atoms%ntype, sym%mrot, ngoprI, atoms%taual, atoms%neq, atoms%lmax, &
        !& atoms%rmt, cell%omtil, cell%bmat, cell%bbmat, kpts%bk(:, ikpq), gbas(1, ilst(:nv(1, ikpq), ikpq, 1)), &
        !& gbas(2, ilst(:nv(1, ikpq), ikpq, 1)), gbas(3, ilst(:nv(1, ikpq), ikpq, 1)), nv(:, ikpq), nmat, &
        !& usdus%us, usdus%dus, usdus%uds, usdus%duds, usdus%ddn, sym%invsat, sym%invsatnr, usdus%ulos, usdus%uulon, usdus%dulon, &
        !& usdus%dulos, atoms%llo, atoms%nlo, atoms%l_dulo, atoms%lapw_l, kveclo(:,ikpq), odi, ods, aKpq, bKpq, bascof_loKpq )
        nk=fmpi%k_list(ikpq)
        CALL lapw%init(input, noco, nococonv, kpts, atoms, sym, nk, cell, .FALSE., fmpi)
        CALL abcof3(input, atoms, sym, 1, cell, kpts%bk(:, ikpq), lapw, &
                            usdus, oneD, 1, lapw%dim_nvd(), aKpq, bKpq, bascof_loKpq)

      nmat = nv(1, ikpt) + atoms%nlotot
      a(:, :, :) = cmplx(0.0, 0.0)
      b(:, :, :) = cmplx(0.0, 0.0)
      bascof_lo(:, :, :, :, :) = cmplx(0.0, 0.0)
     ! call abcof3( atoms%lmaxd, atoms%ntype, atoms%nat, sym%nop, MAXVAL(nv), input%jspins, 1, atoms%lmaxd*(atoms%lmaxd+2), SIZE(z(:,1,1,1)), &
    !    & atoms%llod, atoms%nlod, atoms%nlotot, sym%invtab, atoms%ntype, sym%mrot, ngoprI, atoms%taual, atoms%neq, atoms%lmax, &
    !    & atoms%rmt, cell%omtil, cell%bmat, cell%bbmat, kpts%bk(:, ikpt), gbas(1, ilst(:nv(1, ikpt), ikpt, 1)), &
    !    & gbas(2, ilst(:nv(1, ikpt), ikpt, 1)), gbas(3, ilst(:nv(1, ikpt), ikpt, 1)), nv(:, ikpt), nmat, &
    !    & usdus%us, usdus%dus, usdus%uds, usdus%duds, usdus%ddn, sym%invsat, sym%invsatnr, usdus%ulos, usdus%uulon, usdus%dulon, &
    !    & usdus%dulos, atoms%llo, atoms%nlo, atoms%l_dulo, atoms%lapw_l, kveclo(:,ikpt), odi, ods, a, b, bascof_lo )
    nk=fmpi%k_list(ikpt)
    CALL lapw%init(input, noco, nococonv, kpts, atoms, sym, nk, cell, .FALSE., fmpi)
    CALL abcof3(input, atoms, sym, 1, cell, kpts%bk(:, ikpt), lapw, &
                        usdus, oneD, 1, lapw%dim_nvd(), a, b, bascof_lo)

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
                ab0cofScl(lmp + 1) = ImagUnit**oqn_l * dot_product( conjg(z(:nv(1, ikpt), iband, ikpt, 1)), a(:nv(1, ikpt), lm, iDatomA) )
                ! p = 2
                ab0cofScl(lmp + 2) = ImagUnit**oqn_l * dot_product( conjg(z(:nv(1, ikpt), iband, ikpt, 1)), b(:nv(1, ikpt), lm, iDatomA) )
                ! Add LO contributions
                do iradf = 3, pMaxLocal
                  ! p = 1
                  ab0cofScl(lmp + 1) = ab0cofScl(lmp + 1) + ImagUnit**oqn_l * &
                    & dot_product( conjg(z(nv(1, ikpt) + 1:nv(1, ikpt) + atoms%nlotot, iband, ikpt, 1)), &
                                                         & bascof_lo(1, mqn_m, :atoms%nlotot, iloTable(iradf, oqn_l, iDtypeA), iDatomA) )
                  ! p = 2
                  ab0cofScl(lmp + 2) = ab0cofScl(lmp + 2) + ImagUnit**oqn_l * &
                    & dot_product( conjg(z(nv(1, ikpt) + 1:nv(1, ikpt) + atoms%nlotot, iband, ikpt, 1)), &
                                                         & bascof_lo(2, mqn_m, :atoms%nlotot, iloTable(iradf, oqn_l, iDtypeA), iDatomA) )
                  ! 2 < p < LOs for that l and that atom type
                  ab0cofScl(lmp + iradf) = ImagUnit**oqn_l * dot_product( conjg(z(nv(1, ikpt) + 1:nv(1, ikpt) + atoms%nlotot, iband, ikpt, 1)), &
                                                         & bascof_lo(3, mqn_m, :atoms%nlotot, iloTable(iradf, oqn_l, iDtypeA), iDatomA) )
                end do ! iradf


                ! This is a precalculation of the 1st and 3rd line of A.50. Actually the resulting lmp should have been primed or at
                ! least the p as lm1Pr = lm. But for sake of performance we place it here.
                ! sake of performance
                lm = lm + 1
                lmp = lmp + pMaxLocal
              end do ! mqn_m
            end do !oqn_l

            ! Calculate vector-like large matching coefficients at k-point k with the unperturbed expansion coefficients z.
            ! This vector has to be fully available so we choose idirR
            ab0cofVec(:, :) = cmplx(0.0, 0.0)
            abcofSumVec(:, :) = cmplx(0.0, 0.0)
            do idirR = 1, 3
              lmp = 0
              lm  = 0
              do oqn_l = 0, atoms%lmax(iDtypeA)
                do mqn_m = -oqn_l, oqn_l
                  pMaxLocal = nRadFun(oqn_l, iDtypeA)
                  ! p = 1
                  ab0cofVec(lmp + 1, idirR) = ImagUnit**oqn_l * dot_product( conjg( z(:nv(1, ikpt), iband, ikpt, 1)), &
                                                                      & gbasExt(:nv(1, ikpt), idirR) * a(:nv(1, ikpt), lm, iDatomA) )
                  ! p = 2
                  ab0cofVec(lmp + 2, idirR) = ImagUnit**oqn_l * dot_product( conjg( z(:nv(1, ikpt), iband, ikpt, 1)), &
                                                                      & gbasExt(:nv(1, ikpt), idirR) * b(:nv(1, ikpt), lm, iDatomA) )
                  do iradf = 3, pMaxLocal
                    ! p = 1
                    ab0cofVec(lmp + 1, idirR) = ab0cofVec(lmp + 1, idirR) + ImagUnit**oqn_l * &
                      & dot_product( conjg( z(nv(1, ikpt) + 1:nv(1, ikpt) + atoms%nlotot, iband, ikpt, 1) ), &
                                                       &   gbasExt(nv(1, ikpt) + 1:nv(1, ikpt) + atoms%nlotot, idirR) &
                                                       & * bascof_lo(1, mqn_m, :atoms%nlotot, iloTable(iradf, oqn_l, iDtypeA), iDatomA) )
                    ! p = 2
                    ab0cofVec(lmp + 2, idirR) = ab0cofVec(lmp + 2, idirR) + ImagUnit**oqn_l * &
                      & dot_product( conjg( z(nv(1, ikpt) + 1:nv(1, ikpt) + atoms%nlotot, iband, ikpt, 1) ), &
                                                       &   gbasExt(nv(1, ikpt) + 1:nv(1, ikpt) + atoms%nlotot, idirR) &
                                                       & * bascof_lo(2, mqn_m, :atoms%nlotot, iloTable(iradf, oqn_l, iDtypeA), iDatomA) )
                    ! 2 < p < LOs for that l and that atom type
                    ab0cofVec(lmp + iradf, idirR) = ImagUnit**oqn_l * dot_product( conjg( z(nv(1, ikpt) + 1:nv(1, ikpt) + atoms%nlotot, iband, ikpt, 1) ),&
                                                       &   gbasExt(nv(1, ikpt) + 1:nv(1, ikpt) + atoms%nlotot, idirR) &
                                                       & * bascof_lo(3, mqn_m, :atoms%nlotot, iloTable(iradf, oqn_l, iDtypeA), iDatomA) )
                  end do ! iradf

                  ! This help array has to be calulated in a seperate loop here because earlier ab0cofVec would not be available and the
                  ! next command needs abcofSumVec already.
                  do iradf = 1, pMaxLocal
                    abcofSumVec(lmp + iradf, idirR) = ImagUnit * (ab0cofScl(lmp + iradf) * kExt(idirR) + ab0cofVec(lmp + iradf, idirR))
                  end do ! iradf

                  ! Precalculation of the 2nd and the 4th line in A.51. Due to performance, the lmp in the resulting quantity are not
                  ! primed
                  lm = lm + 1
                  lmp = lmp + pMaxLocal
                end do ! mqn_m
              end do ! oqn_l
            end do ! idirR

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
                      ab1cofVec(lmp + 1, iDatomB, idirR) = ImagUnit**oqn_l * &
                                     & dot_product( conjg( z1nG(:nv(1, ikpq), idirR, iDatomA, iband) ), aKpq(:nv(1, ikpq), lm, iDatomB) )
                      ! p = 2
                      ab1cofVec(lmp + 2, iDatomB, idirR) = ImagUnit**oqn_l * &
                                     & dot_product( conjg( z1nG(:nv(1, ikpq), idirR, iDatomA, iband) ), bKpq(:nv(1, ikpq), lm, iDatomB) )
                      do iradf = 3, pMaxLocal
                        ! p = 1
                        ab1cofVec(lmp + 1, iDatomB, idirR) = ab1cofVec(lmp + 1, iDatomB, idirR) &
                          & + ImagUnit**oqn_l * dot_product( conjg( z1nG(nv(1, ikpq) + 1:nv(1, ikpq) + atoms%nlotot, idirR, iDatomA, iband) ), &
                                                      & bascof_loKpq(1, mqn_m, :atoms%nlotot, iloTable(iradf, oqn_l, iDtypeB), iDatomB) )
                        ! p = 2
                        ab1cofVec(lmp + 2, iDatomB, idirR) = ab1cofVec(lmp + 2, iDatomB, idirR) &
                          & + ImagUnit**oqn_l * dot_product(conjg(z1nG(nv(1, ikpq) + 1:nv(1, ikpq) + atoms%nlotot, idirR, iDatomA, iband)), &
                                                      & bascof_loKpq(2, mqn_m, :atoms%nlotot, iloTable(iradf, oqn_l, iDtypeB), iDatomB) )
                        ! 2 < p < LOs for that l and that atom type
                        ab1cofVec(lmp + iradf, iDatomB, idirR) = &
                          & ImagUnit**oqn_l * dot_product( conjg( z1nG(nv(1, ikpq) + 1:nv(1, ikpq) + atoms%nlotot, idirR, iDatomA, iband) ),&
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
            grVarphiHPsiNat(:, :) = cmplx(0.0, 0.0)
            grVarphiPsiNat(:, :) = cmplx(0., 0.)
            varphiHGrPsiNat(:) = cmplx(0.0, 0.0)
            varphiGrPsiNat(:) = cmplx(0., 0.)
            do idirC = 1, 3
              do idirR = 1, 3
                lm  = 0
                lmp = 0
                abcofSumMatAlpha(:) = cmplx(0.0, 0.0)
                do oqn_l = 0, atoms%lmax(iDtypeA)
                  do mqn_m = -oqn_l, oqn_l
                    pMaxLocal = nRadFun(oqn_l, iDtypeA)
                    abcofMat(:) = cmplx(0.0, 0.0)
                    ! p = 1
                    abcofMat(1) = ImagUnit**oqn_l * dot_product( conjg( z(:nv(1, ikpt), iband, ikpt, 1) ), &
                                                                & GbasMatExt(:nv(1, ikpt), idirR, idirC) * a(:nv(1, ikpt), lm, iDatomA) )
                    ! p = 2
                    abcofMat(2) = ImagUnit**oqn_l * dot_product( conjg( z(:nv(1, ikpt), iband, ikpt, 1) ), &
                                                                & GbasMatExt(:nv(1, ikpt), idirR, idirC) * b(:nv(1, ikpt), lm, iDatomA) )
                    ! Add LOs
                    do iradf = 3, pMaxLocal
                      ! p = 1
                      abcofMat(1) = abcofMat(1)  + ImagUnit**oqn_l * dot_product( conjg( z(nv(1, ikpt) + 1:nv(1, ikpt) + atoms%nlotot, iband, ikpt, 1) ), &
                                                       &   GbasMatExt(nv(1, ikpt) + 1:nv(1, ikpt) + atoms%nlotot, idirR, idirC) &
                                                       & * bascof_lo(1, mqn_m, :atoms%nlotot, iloTable(iradf, oqn_l, iDtypeA), iDatomA) )
                      ! p = 2
                      abcofMat(2) = abcofMat(2)  + ImagUnit**oqn_l * dot_product( conjg( z(nv(1, ikpt) + 1:nv(1, ikpt) + atoms%nlotot, iband, ikpt, 1) ), &
                                                       &   GbasMatExt(nv(1, ikpt) + 1:nv(1, ikpt) + atoms%nlotot, idirR, idirC) &
                                                       & * bascof_lo(2, mqn_m, :atoms%nlotot, iloTable(iradf, oqn_l, iDtypeA), iDatomA) )
                      ! 2 < p < LOs for that l and that atom type
                      abcofMat(iradf) = ImagUnit**oqn_l * dot_product(conjg(z(nv(1, ikpt) + 1:nv(1, ikpt) + atoms%nlotot, iband, ikpt, 1)), &
                                                       &   GbasMatExt(nv(1, ikpt) + 1:nv(1, ikpt) + atoms%nlotot, idirR, idirC) &
                                                       & * bascof_lo(3, mqn_m, :atoms%nlotot, iloTable(iradf, oqn_l, iDtypeA), iDatomA) )
                    end do ! iradf

                    do iradf = 1, pMaxLocal
                      ! Sum up the large matching coefficients with matrix character at k-point ikpt required for the 3rd line in A.50
                      ! PhD thesis of A. Klueppelberg. The dyadic product is expressed by a proper choice of the idir loop variables.
                      ! Here, the lmp should have been primed actually, but for sake of performance we place it here. We do not conjugate
                      ! this quantity here but when the dot_product is performed.
                      abcofSumMatAlpha(lmp + iradf) = ab0cofScl(lmp + iradf) * kExtkExtT(idirR, idirC) + &
                        & abcofMat(iradf) + ab0cofVec(lmp + iradf, idirR) * kExt(idirC) + kExt(idirR) * ab0cofVec(lmp + iradf, idirC)
                    end do ! iradf

                    lm = lm + 1
                    lmp = lmp + pMaxLocal
                  end do ! mqn_m
                end do ! oqn_l
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
                        abcofMat(1) = ImagUnit**oqn_l * dot_product( conjg(z1Gext(:nv(1, ikpq)) ), aKpq(:nv(1, ikpq), lm, iDatomB) )
                        ! p = 2
                        abcofMat(2) = ImagUnit**oqn_l * dot_product( conjg(z1Gext(:nv(1, ikpq)) ), bKpq(:nv(1, ikpq), lm, iDatomB) )
                        do iradf = 3, pMaxLocal
                          ! p = 1
                          abcofMat(1) = abcofMat(1) &
                            & + ImagUnit**oqn_l * dot_product( conjg( z1Gext(nv(1, ikpq) + 1:nv(1, ikpq) + atoms%nlotot)), &
                                                      & bascof_loKpq(1, mqn_m, :atoms%nlotot, iloTable(iradf, oqn_l, iDtypeB), iDatomB) )
                          ! p = 2
                          abcofMat(2) = abcofMat(2) &
                            & + ImagUnit**oqn_l * dot_product( conjg( z1Gext(nv(1, ikpq) + 1:nv(1, ikpq) + atoms%nlotot)), &
                                                     & bascof_loKpq(2, mqn_m, :atoms%nlotot, iloTable(iradf, oqn_l, iDtypeB), iDatomB) )
                          ! 2 < p < LOs for that l and that atom type
                          abcofMat(iradf) = ImagUnit**oqn_l * dot_product( conjg(z1Gext(nv(1, ikpq) + 1:nv(1, ikpq) + atoms%nlotot) ) &
                                                    &, bascof_loKpq(3, mqn_m, :atoms%nlotot, iloTable(iradf, oqn_l, iDtypeB), iDatomB) )
                        end do ! p

                        do iradf = 1, nRadFun(oqn_l, iDtypeB)
                          ! Precalculation of the 1st line in A.50 PhD thesis A. Klueppelberg. Here a plus sign stand before the 2
                          ! because abcofSumMat is complex conjugated later!
                          ! todo B should be A here in ab1cofVec, because A is displaced, one iterates over MTs B, in ab1cofVec,
                          !      abcofMat also assumes A to be displaced because of z1Gext, although B should be displaced, therefore
                          !      abcofSumMat should be dependent on iDatomB(displaced) and iDatomA(iterate over)
                          abcofSumMat(lmp + iradf, iDatomB, idirR, idirC) = &
                            & ImagUnit * ( ab1cofVec(lmp + iradf, iDatomB, idirR) * kpqExt(idirC) + abcofMat(iradf) )
                        ! This line should be activated if we want to test whehter it cancels away without looking at the that are remaining for the Goldstone condition.
                        !    & iu * ( ab1cofVec(lmp + iradf, iDatomB, idirR) * kpqExt(idirC) + abcofMat(iradf) )
                        end do ! iradf
                        lm = lm + 1
                        lmp = lmp + pMaxLocal
                      end do ! mqn_m
                    end do ! oqn_l
                    if ( (iDatomB == iDatomA) .and. .not.eps1DynMatPulTestSw ) then
                      ! The 3rd line in A.50 PhD thesis A. Klueppelberg is added sothat we have all large matching coefficients of the
                      ! bra for <varphi|H - eps|varphi>.
                      ! todo the displaced atom is actually A, instead it should be B, which is irrele
                      abcofSumMat(:lmpT(iDtypeB), iDatomB, idirR, idirC) = abcofSumMat(:lmpT(iDtypeB), iDatomB, idirR, idirC) &
                                                                                                       & - abcofSumMatAlpha(:lmpT(iDtypeA))
                    end if
                  end do ! iDeqatB
                end do ! iDtypeB

              end do ! idirR
              ! Precalculation of the 2nd and the 4th line in A.50 PhD thesis A. Klueppelberg
              grVarphiHPsiNat(:lmpT(iDtypeA), idirC) = matmul( grVarphiHvarphiNat(:lmpT(iDtypeA), :lmpT(iDtypeA), idirC, iDatomA), ab0cofScl(:lmpT(iDtypeA)) )
              ! Precalculation of the 2nd and the 4th line in A.50. todo really placing it here?
              grVarphiPsiNat(1:lmpT(iDtypeA), idirC) = matmul( grVarphiVarphiNat(1:lmpT(iDtypeA), 1:lmpT(iDtypeA), idirC, iDtypeA), ab0cofScl(1:lmpT(iDtypeA)) )
            end do ! idirC


            ! Precalculation of the 1st and 3rd line in A.50 PhD thesis A. Klueppelberg.
            varphiPsi(:) = cmplx(0., 0.)
            varphiPsi(1:lmpT(iDtypeA)) = matmul( varphiVarphi(1:lmpT(iDtypeA), 1:lmpT(iDtypeA), iDtypeA), ab0cofScl(1:lmpT(iDtypeA)) )
            varphiHPsi(:) = cmplx(0.0, 0.0)
            varphiHPsi(1:lmpT(iDtypeA)) = matmul( varphiHvarphi(1:lmpT(iDtypeA), 1:lmpT(iDtypeA), iDatomA), ab0cofScl(1:lmpT(iDtypeA)) )


            iDatomB = 0
            do iDtypeB = 1, atoms%ntype
              do iDeqatB = 1, atoms%neq(iDtypeB)
                iDatomB = iDatomB + 1

                ! todo The resetting of these arrays might be wrong in this place, it was before within the if-statement iDatomB=iDatomA
                grPsiHPsiNat(:, :) = cmplx(0.0, 0.0)
                grPsiPsiNat(:, :) = cmplx(0., 0.)
                grPsiHPsiNatT(:, :) = cmplx(0.0, 0.0)
                grPsiPsiNatT(:, :) = cmplx(0., 0.)
                psiHGrPsiNat(:, :) = cmplx(0., 0.)
                psiHGrPsiNatT(:, :) = cmplx(0., 0.)
                psiGrVeff0sphPsi(:, :) = cmplx(0., 0.)
                psiGrVeff0sphPsiT(:, :) = cmplx(0., 0.)
                psiGrVeff0Psi(:, :) = cmplx(0., 0.)
                psiGrVeff0PsiT(:, :) = cmplx(0., 0.)
                psiGrPsiNat(:, :) = cmplx(0., 0.)
                psiGrPsiNatT(:, :) = cmplx(0., 0.)
                grPsiHepsPsiNatT(:, :) = cmplx(0., 0.)
                psiHepsGrPsiNatT(:, :) = cmplx(0., 0.)
                if ( (iDatomB == iDatomA) .and. .not.eps1DynMatPulTestSw ) then
                  ! Calculation of the 4th line in A.50 where we have the transpose. Therefore this has to be in a seperate loop structure
                  do idirC = 1, 3
                    do idirR = 1, 3
                      ! Calculation of the 4th line in A.50 PhD thesis A. Klueppelberg. This is added to the 2nd line in A.50 PhD thesis
                      ! A. Klueppelberg later.
                      grPsiHPsiNat(idirR, idirC)  = &
                                       & - dot_product( abcofSumVec(:lmpT(iDtypeB), idirR), grVarphiHPsiNat(:lmpT(iDtypeB), idirC) )
                      ! Transposed version of upper line
                      grPsiHPsiNatT(idirR, idirC) = &
                                       & - dot_product( abcofSumVec(:lmpT(iDtypeB), idirC), grVarphiHPsiNat(:lmpT(iDtypeB), idirR) )

                      ! Calculation of the 4th line in A.50. ! This has to stay in a separate loop construction of idirR, idirC, because
                      grPsiPsiNat(idirR, idirC)  = &
                                        & - dot_product( abcofSumVec(:lmpT(iDtypeB), idirR), grVarphiPsiNat(:lmpT(iDtypeB), idirC) )
                      ! Transposed version of upper line
                      grPsiPsiNatT(idirR, idirC) = &
                        !todo maybe sign incorrect
                                        & - dot_product( abcofSumVec(:lmpT(iDtypeB), idirC), grVarphiPsiNat(:lmpT(iDtypeB), idirR) )


                      varphiHGrPsiNat(:) = cmplx(0., 0.)
                      varphiHGrPsiNat(1:lmpT(iDtypeB)) = matmul(varphiHGrvarphiNat(1:lmpT(iDtypeB), 1:lmpT(iDtypeB), iDatomB, idirC), abcofSumVec(1:lmpT(iDtypeB), idirR))
                      psiHGrPsiNat(idirR, idirC) = -dot_product(ab0cofScl(1:lmpT(iDtypeB)), varphiHGrPsiNat(1:lmpT(iDtypeB)))

                      varphiHGrPsiNat(:) = cmplx(0., 0.)
                      varphiHGrPsiNat(1:lmpT(iDtypeB)) = matmul(varphiHGrvarphiNat(1:lmpT(iDtypeB), 1:lmpT(iDtypeB), iDatomB, idirR), abcofSumVec(1:lmpT(iDtypeB), idirC))
                      psiHGrPsiNatT(idirR, idirC) = -dot_product(ab0cofScl(1:lmpT(iDtypeB)), varphiHGrPsiNat(1:lmpT(iDtypeB)))

                      varphiGrVeff0sphPsi(:) = cmplx(0., 0.)
                      varphiGrVeff0sphPsi(1:lmpT(iDtypeB)) = matmul(varphiGrVeff0SphVarphi(1:lmpT(iDtypeB), 1:lmpT(iDtypeB), iDatomB, idirC), abcofSumVec(1:lmpT(iDtypeB), idirR))
                      psiGrVeff0sphPsi(idirR, idirC) = -dot_product(ab0cofScl(1:lmpT(iDtypeB)), varphiGrVeff0sphPsi(1:lmpT(iDtypeB)))

                      varphiGrVeff0Psi(:) = cmplx(0., 0.)
                      varphiGrVeff0Psi(1:lmpT(iDtypeB)) = matmul(varphiGrVeff0Varphi(1:lmpT(iDtypeB), 1:lmpT(iDtypeB), iDatomB, idirR), abcofSumVec(1:lmpT(iDtypeB), idirC))
                      psiGrVeff0Psi(idirR, idirC) = 2 * dot_product(ab0cofScl(1:lmpT(iDtypeB)), varphiGrVeff0Psi(1:lmpT(iDtypeB)))

                      varphiGrVeff0sphPsi(:) = cmplx(0., 0.)
                      varphiGrVeff0sphPsi(1:lmpT(iDtypeB)) = matmul(varphiGrVeff0SphVarphi(1:lmpT(iDtypeB), 1:lmpT(iDtypeB), iDatomB, idirR), abcofSumVec(1:lmpT(iDtypeB), idirC))
                      psiGrVeff0sphPsiT(idirR, idirC) = -dot_product(ab0cofScl(1:lmpT(iDtypeB)), varphiGrVeff0sphPsi(1:lmpT(iDtypeB)))

                      varphiGrVeff0Psi(:) = cmplx(0., 0.)
                      varphiGrVeff0Psi(1:lmpT(iDtypeB)) = matmul(varphiGrVeff0Varphi(1:lmpT(iDtypeB), 1:lmpT(iDtypeB), iDatomB, idirC), ab0cofScl(1:lmpT(iDtypeB)) )
                      psiGrVeff0PsiT(idirR, idirC) = 2 * dot_product(abcofSumVec(1:lmpT(iDtypeB), idirR), varphiGrVeff0Psi(1:lmpT(iDtypeB)))

                      ! todo should be alpha here instead of beta
                      varphiGrPsiNat(:) = cmplx(0., 0.)
                      varphiGrPsiNat(1:lmpT(iDtypeB)) = matmul(varphiGrVarphiNat(1:lmpT(iDtypeB), 1:lmpT(iDtypeB), idirC, iDtypeB), abcofSumVec(1:lmpT(iDtypeB), idirR))
                      psiGrPsiNat(idirR, idirC) = -dot_product(ab0cofScl(1:lmpT(iDtypeB)), varphiGrPsiNat(1:lmpT(iDtypeB)))

                      varphiGrPsiNat(:) = cmplx(0., 0.)
                      varphiGrPsiNat(1:lmpT(iDtypeB)) = matmul(varphiGrVarphiNat(1:lmpT(iDtypeB), 1:lmpT(iDtypeB), idirR, iDtypeB), abcofSumVec(1:lmpT(iDtypeB), idirC))
                      psiGrPsiNatT(idirR, idirC) = -dot_product(ab0cofScl(1:lmpT(iDtypeB)), varphiGrPsiNat(1:lmpT(iDtypeB)))

                    end do ! idirR
                  end do ! idirC
                  if ( testComp2ndN1stOrdBasFuncSw ) then
                    do idirC = 1, 3
                      do idirR = 1, 3
                        if ( abs( grPsiHPsiNat(idirR, idirC) - grPsiHPsiNatT(idirC, idirR) ) > 1e-10 ) then
                          call juDFT_error( 'Test activated by testComp2ndN1stOrdBasFuncSw failed!', &
                          & calledby='Add2ndOrdWfPulayBraKets2DynMat', hint='grPsiHPsiNat is not transposed to grPsiHPsiNatT' )
                        end if
                        if ( abs( grPsiPsiNat(idirR, idirC) - grPsiPsiNatT(idirC, idirR) ) > 1e-10 ) then
                          call juDFT_error( 'Test activated by testComp2ndN1stOrdBasFuncSw failed!', &
                          & calledby='Add2ndOrdWfPulayBraKets2DynMat', hint='grPsiPsiNat is not transposed to grPsiPsiNatT' )
                        end if
                        if ( abs( psiHGrPsiNat(idirR, idirC) - psiHGrPsiNatT(idirC, idirR) ) > 1e-10 ) then
                          call juDFT_error( 'Test activated by testComp2ndN1stOrdBasFuncSw failed!', &
                          & calledby='Add2ndOrdWfPulayBraKets2DynMat', hint='grPsiHPsiNat is not transposed to grPsiHPsiNatT' )
                        end if
                        if ( abs( psiGrPsiNat(idirR, idirC) - psiGrPsiNatT(idirC, idirR) ) > 1e-10 ) then
                          call juDFT_error( 'Test activated by testComp2ndN1stOrdBasFuncSw failed!', &
                          & calledby='Add2ndOrdWfPulayBraKets2DynMat', hint='grPsiPsiNat is not transposed to grPsiPsiNatT' )
                        end if
                        if ( abs( psiGrVeff0sphPsi(idirR, idirC) - psiGrVeff0sphPsiT(idirC, idirR) ) > 1e-10 ) then
                          call juDFT_error( 'Test activated by testComp2ndN1stOrdBasFuncSw failed!', &
                          & calledby='Add2ndOrdWfPulayBraKets2DynMat', hint='grPsiHPsiNat is not transposed to grPsiHPsiNatT' )
                        end if
                      end do ! idirR
                    end do ! idirC
                    psiGrVeff0sphPsiT(:, :) = cmplx(0., 0.)
                  else
                    grPsiHepsPsiNatT(1:3, 1:3) = grPsiHPsiNatT(1:3, 1:3) - eig(iband, ikpt, 1) * grPsiPsiNatT(1:3, 1:3)
                    psiHepsGrPsiNatT(1:3, 1:3) = psiHGrPsiNatT(1:3, 1:3) - eig(iband, ikpt, 1) * psiGrPsiNatT(1:3, 1:3)
                  end if
                end if ! iDatomA = iDatomB

                psiHepsPsi(:, :) = cmplx(0., 0.)
                grPsiHepsPsiNat(:, :) = cmplx(0., 0.)
                psiHepsGrPsiNat(:, :) = cmplx(0., 0.)
                do idirC = 1, 3
                  do idirR = 1, 3

                    ! todo Here for multiple atoms there is a bug in the next 10 lines, because iDatomB and iDatomA are not addressed correctly
                    psiVarphi(:) = cmplx(0., 0.)
                    psiVarphi(1:lmpT(iDtypeA)) = matmul( varphiVarphi(1:lmpT(iDtypeA), 1:lmpT(iDtypeA), iDtypeA), abcofSumMat(1:lmpT(iDtypeB), iDatomB, idirR, idirC))
                    psiHVarphi(:) = cmplx(0., 0.)
                    psiHVarphi(1:lmpT(iDtypeA)) = matmul( varphiHVarphi(1:lmpT(iDtypeA), 1:lmpT(iDtypeA), iDtypeA), abcofSumMat(1:lmpT(iDtypeB), iDatomB, idirR, idirC) )
                    ! Calculation of the 2nd and 4th line of A.50 PhD thesis A. Klueppelberg. The distinction between atom alpha and
                    ! beta in abcofSumMat has already be performed upwards.
                    psiHepsPsi(idirR, idirC) = dot_product( abcofSumMat(:lmpT(iDtypeB), iDatomB, idirR, idirC), varphiHPsi(:lmpT(iDtypeB)) ) &
                                               & - eig(iband, ikpt, 1) * dot_product( abcofSumMat(:lmpT(iDtypeA), iDatomB, idirR, idirC), varphiPsi(:lmpT(iDtypeA)) ) &
                                              &+dot_product( ab0cofScl(1:lmpT(iDtypeA)), psiHVarphi(1:lmpT(iDtypeA)) ) &
                                             - eig(iband, ikpt, 1) * dot_product( ab0cofScl(1:lmpT(iDtypeA)), psiVarphi(1:lmpT(iDtypeA)) )

                    !! Calculation of A.52 for f = varphi
                    !!todo check the assignment to the atoms also why does the loop over alpha not end up there

                    ! The 4th line of A.50 PhD thesis A. Klueppelberg is calculated and added to its 2nd line
                    ! Here, the ab1cofVec should be actually displaced B and iterating of A, it is vice versa here.
                    grPsiHPsiNat(idirR, idirC) = grPsiHPsiNat(idirR, idirC) &
                                                           & - 2 * dot_product( ab1cofVec(:lmpT(iDtypeB), iDatomB, idirR), grVarphiHPsiNat(:lmpT(iDtypeB), idirC) )
                    grPsiPsiNat(idirR, idirC) = grPsiPsiNat(idirR, idirC) &
                                                           & - 2 * dot_product(ab1cofVec(:lmpT(iDtypeB), iDatomB, idirR), grVarphiPsiNat(:lmpT(iDtypeB), idirC) )
                    grPsiHepsPsiNat(idirR, idirC) = grPsiHPsiNat(idirR, idirC) - eig(iband, ikpt, 1) * grPsiPsiNat(idirR, idirC)

                    varphiHGrPsiNat(:) = cmplx(0., 0.)
                    varphiHGrPsiNat(1:lmpT(iDtypeB)) = matmul(varphiHGrvarphiNat(1:lmpT(iDtypeB), 1:lmpT(iDtypeB), iDatomB, idirC), ab1cofVec(1:lmpT(iDtypeB), iDatomB, idirR))
                    psiHGrPsiNat(idirR, idirC) = psiHGrPsiNat(idirR, idirC) &
                                                     & - 2 * dot_product( ab0cofScl(1:lmpT(iDtypeB)), varphiHGrPsiNat(1:lmpT(iDtypeB)))


                    varphiGrVeff0sphPsi(:) = cmplx(0., 0.)
                    varphiGrVeff0sphPsi(1:lmpT(iDtypeB)) = matmul(varphiGrVeff0SphVarphi(1:lmpT(iDtypeB), 1:lmpT(iDtypeB), iDatomB, idirC), ab1cofVec(1:lmpT(iDtypeB), iDatomB, idirR))
                    psiGrVeff0sphPsi(idirR, idirC) = psiGrVeff0sphPsi(idirR, idirC) - 2 * dot_product(ab0cofScl(1:lmpT(iDtypeB)), varphiGrVeff0sphPsi(1:lmpT(iDtypeB)))

                    varphiGrVeff0Psi(:) = cmplx(0., 0.)
                    varphiGrVeff0Psi(1:lmpT(iDtypeB)) = matmul(varphiGrVeff0Varphi(1:lmpT(iDtypeB), 1:lmpT(iDtypeB), iDatomB, idirR), ab1cofVec(1:lmpT(iDtypeB), iDatomB, idirC))
                    psiGrVeff0Psi(idirR, idirC) = psiGrVeff0Psi(idirR, idirC) + 2 * dot_product(ab0cofScl(1:lmpT(iDtypeB)), varphiGrVeff0Psi(1:lmpT(iDtypeB)))

                    varphiGrVeff0Psi(:) = cmplx(0., 0.)
                    varphiGrVeff0Psi(1:lmpT(iDtypeB)) = matmul(varphiGrVeff0Varphi(1:lmpT(iDtypeB), 1:lmpT(iDtypeB), iDatomB, idirC), ab0cofScl(1:lmpT(iDtypeB)))
                    psiGrVeff0PsiT(idirR, idirC) = psiGrVeff0PsiT(idirC, idirR) + 2 * dot_product(ab1cofVec(1:lmpT(iDtypeB), iDatomB, idirR), varphiGrVeff0Psi(1:lmpT(iDtypeB)))

                    ! todo actually be alpha here instead of beta
                    varphiGrPsiNat(:) = cmplx(0., 0.)
                    varphiGrPsiNat(1:lmpT(iDtypeB)) = matmul(varphiGrVarphiNat(1:lmpT(iDtypeB), 1:lmpT(iDtypeB), idirC, iDtypeB), ab1cofVec(1:lmpT(iDtypeB), iDatomB, idirR))
                    psiGrPsiNat(idirR, idirC) = psiGrPsiNat(idirR, idirC) &
                      & - 2 * dot_product( ab0cofScl(1:lmpT(iDtypeB)), varphiGrPsiNat(1:lmpT(iDtypeB)))
                    psiHepsGrPsiNat(idirR, idirC) = psiHGrPsiNat(idirR, idirC) - eig(iband, ikpt, 1) * psiGrPsiNat(idirR, idirC)

                  end do ! idirR
                end do ! idirC

                ! Transformations of natural coordinates into cartesian coordinates. Therefore, we have to close the loops over idirC,
                ! and idirR
                grPsiHepsPsi(:, :) = cmplx(0.0, 0.0)
                grPsiHepsPsiT(:, :) = cmplx(0.0, 0.0)
                grPsiHepsPsi(1:3, 1:3) = matmul(grPsiHepsPsiNat(1:3, 1:3), transpose(conjg(Tmatrix0(1:3, 1:3))))
                grPsiHepsPsiT(1:3, 1:3) = matmul(conjg(Tmatrix0(1:3, 1:3)), grPsiHepsPsiNatT(1:3, 1:3))
                grPsiHepsPsi(1:3, 1:3) = grPsiHepsPsi(1:3, 1:3) + grPsiHepsPsiT(1:3, 1:3)

                psiHepsGrPsi(:, :) = cmplx(0., 0.)
                psiHepsGrPsiT(:, :) = cmplx(0., 0.)
                psiHepsGrPsi(1:3, 1:3) = matmul(psiHepsGrPsiNat(1:3, 1:3), transpose(conjg(Tmatrix0(1:3, 1:3))))
                psiHepsGrPsiT(1:3, 1:3) = matmul(Tmatrix0(1:3, 1:3), psiHepsGrPsiNatT(1:3, 1:3))
                psiHepsGrPsi(1:3, 1:3) = psiHepsGrPsi(1:3, 1:3) + psiHepsGrPsiT(1:3, 1:3) - psiGrVeff0sphPsi(1:3, 1:3) - psiGrVeff0sphPsiT(1:3, 1:3)
                do idirC = 1, 3
                  do idirR = 1, 3

                    ! Setup of complete A.50 in PhD thesis A. Klueppelberg.
                    matElemHeps = psiHepsPsi(idirR, idirC) !+ psiGrVeff0Psi(idirR, idirC) + psiGrVeff0PsiT(idirR, idirC)!+ grPsiHepsPsi(idirR, idirC) + psiHepsGrPsi(idirR, idirC)

                    ! Add the equation A.50 and A.51 with its corrections to the Dynamical Matrix.
                    dynMatPu(idirR + (iDatomB - 1) * 3, idirC + (iDatomA - 1) * 3) =     &
                      &   dynMatPu(idirR + (iDatomB - 1) * 3, idirC + (iDatomA - 1) * 3) &
                      & + 2 * results%w_iks(iband, ikpt, 1) * matElemHeps
                    !write(465,*) 'Pu 2      ', idirR, idirC, iband, ikpt
                    !write(465,*) 2 * results%w_iks(iband, ikpt, 1) * matElemHeps
                  end do ! idirR
                end do ! idirC
              end do ! iDeqatB
            end do ! iDtypeB
          end do ! iDeqatA
        end do ! iDtypeA
      end do ! iband

  end subroutine Add2ndOrdWfPulayBraKets2DynMat


  !>-------------------------------------------------------------------------------------------------------------------------------------
  !> @author
  !> Christian-Roman Gerhorst
  !>
  !> @brief
  !> Implementation of A.53 PhD thesis A. Klueppelberg behind the second equals sign.
  !>
  !> @details
  !> Third braket within the Pulay contributions to the Dynamical matrix (see 7.118 PhD thesis A. Klueppelberg). The matrix elements of
  !> basis functions are given as input parameter. This subroutine creates the adequate versions of the large matching coefficients and
  !> multiplies them to the bais function matrix elements to gain A.53 PhD thesis A. Klueppelberg. Finally, it is added to the dynMatPu
  !> array.
  !>
  !> @note
  !> The filling for the z1 is not optimal in memory access, but minimizes array sizes, and ensures a linear run through the array later.
  !>
  !> @todo
  !> Account for the spin.
  !>-------------------------------------------------------------------------------------------------------------------------------------
  subroutine Add1stOrdWfPulayBraKets2DynMat( fmpi, noco, nococonv, oneD, atoms, input, kpts, qpts, sym, cell, usdus, stars, results, ikpt, iqpt, lmpMax, nRadFunMax, nv, gBas, &
      & gBasUnwrap, kveclo, mapKpq2K, nobd, z, z1nG, iloTable, grVarphiVarphiNat, nRadFun, eig, kpq2kPrVec, &
      & grVarphiHvarphiNat, varphiVarphi, varphiHvarphi, vEff0IR, lmpT, eps1DynMatPulTestSw, &
      & testCompTerm3rdBraKetsVarBra, testCompTerm3rdBraKetsVarKet, dynMatPu3rdBraKetHepsSw, varphiGrVeff0SphVarphi, varphiHGrvarphiNat, varphiGrVarphiNat, dynMatPu )

    use m_types
    use m_abcof3

    implicit none

    ! Type parameters
    type(t_mpi),                  intent(in)  :: fmpi
    type(t_nococonv),                intent(in)  :: nococonv
    type(t_noco),                intent(in)  :: noco
    type(t_oneD),                intent(in)  :: oneD
    type(t_atoms),                  intent(in)    :: atoms
    type(t_input),                  intent(in)    :: input
    type(t_kpts),                   intent(in)    :: kpts
    type(t_kpts),                   intent(in)    :: qpts
    type(t_sym),                    intent(in)    :: sym
    type(t_cell),                   intent(in)    :: cell
    type(t_usdus),                  intent(in)    :: usdus
    type(t_stars),                  intent(in)    :: stars
    type(t_results),                intent(in)    :: results

    ! Scalar parameters
    integer,                        intent(in)    :: ikpt
    integer,                        intent(in)    :: iqpt
    integer,                        intent(in)    :: lmpMax
    integer,                        intent(in)    :: nRadFunMax
    logical,                        intent(in)    :: eps1DynMatPulTestSw
    logical,                        intent(in)    :: testCompTerm3rdBraKetsVarBra
    logical,                        intent(in)    :: testCompTerm3rdBraKetsVarKet
    logical,                        intent(in)    :: dynMatPu3rdBraKetHepsSw

    ! Array parameters
    integer,                        intent(in)    :: nv(:, :)
    integer,                        intent(in)    :: lmpT(:)
    integer,                        intent(in)    :: gBas(:, :)
    integer,                        intent(in)    :: gBasUnwrap(:, :, :)
    integer,                        intent(in)    :: kveclo(:,:)
    integer,                        intent(in)    :: mapKpq2K(:, :)
    integer,                        intent(in)    :: nobd(:, :)
    complex,                       intent(in)    :: z(:,:,:,:)
    complex,                        intent(in)    :: z1nG(:, :, :, :)
    integer,                        intent(in)    :: iloTable(:, 0:, :)
    real,                           intent(in)    :: grVarphiVarphiNat(:, :, :, :)
    integer,                        intent(in)    :: nRadFun(0:, :)
    real,                           intent(in)    :: eig(:, :, :)
    integer,                        intent(in)    :: kpq2kPrVec(:, :, :)
    complex,                        intent(in)    :: grVarphiHvarphiNat(:, :, :, :)
    real,                           intent(in)    :: varphiVarphi(:, :, :)
    complex,                        intent(in)    :: varphiHvarphi(:, :, :)
    complex,                        intent(in)    :: varphiGrVeff0SphVarphi(:, :, :, :)
    complex,                        intent(in)    :: varphiHGrvarphiNat(:, :, :, :)
    complex,                        intent(in)    :: vEff0IR(:,:)
    real,                           intent(in)    :: varphiGrVarphiNat(:, :, :, :)
    complex,                        intent(inout) :: dynMatPu(:, :)

    ! Type variable
    type(od_inp)                                  :: odi
    type(od_sym)                                  :: ods
    type(t_lapw) :: lapw

    ! Scalar variables
    integer                                       :: nmat
    integer                                       :: ikpq
    integer                                       :: iBas
    integer                                       :: iband
    integer                                       :: iDatomA
    integer                                       :: iDtypeA
    integer                                       :: iDeqatA
    integer                                       :: lm
    integer                                       :: lmp
    integer                                       :: oqn_l
    integer                                       :: mqn_m
    integer                                       :: iradf
    integer                                       :: idirC
    integer                                       :: lmp1Pr
    integer                                       :: iatomG
    integer                                       :: itypeG
    integer                                       :: ieqatG
    integer                                       :: idirR
    integer                                       :: oqn_l1Pr
    integer                                       :: mqn_m1Pr
    integer                                       :: pMaxLocal
    complex                                       :: wholeUnitCellIR
    complex                                       :: wholeUnitCellMT
    integer                                       :: iDatomB
    integer                                       :: iDtypeB
    integer                                       :: iDeqatB
    integer                                       :: iradf1Pr
    complex                                       :: sAdjCorrPsiGrPsiNat
    integer :: nk

    ! Array variabels
    integer,           allocatable                :: ngoprI(:)
    complex,           allocatable                :: a(:, :, :)
    complex,           allocatable                :: b(:, :, :)
    complex,           allocatable                :: bascof_lo(:, :, :, :, :)
    complex,           allocatable                :: aKpq(:, :, :)
    complex,           allocatable                :: bKpq(:, :, :)
    complex,           allocatable                :: bascof_loKpq(:, :, :, :, :)
    real,              allocatable                :: gBasExt(:, :)
    complex,           allocatable                :: ab0cofScl(:, :)
    complex,           allocatable                :: ab0cofVec(:)
    !complex,           allocatable                :: ab0KpGcof(:, :, :) seems obsolete
    complex,           allocatable                :: grPsiHepsVarphiOffDiagNat(:, :, :)
    complex,           allocatable                :: varphiHepsGrPsiOffDiagNat(:, :, :)
    complex,           allocatable                :: allMTvarphiHepsPsi(:, :)
    complex,           allocatable                :: ab1cofVec(:, :, :, :)
    complex,           allocatable                :: grVarphiHepsPsiOffDiagNat(:)
    complex,           allocatable                :: psiHepsgrVarphiNat(:)
    complex,           allocatable                :: grVarphiHepsPsiNat(:)
    complex,           allocatable                :: varphiHepsGrPsiOffDiag(:, :)
    !complex,           allocatable                :: psiHepsVarphi(:, :, :) ! seems obsolete
    complex,           allocatable                :: abGrPsiSumCof(:, :, :, :)
    complex,           allocatable                :: psiHepsgrPsiNat(:, :)
    complex,           allocatable                :: varphiHepsGrPsiNat(:)
    complex,           allocatable                :: varphiGrVeff0sphPsi(:)
    complex,           allocatable                :: varphiHepsPsi(:, :)
    complex,           allocatable                :: abGrPsiSumCofDiag(:, :, :)
    complex,           allocatable                :: varphiHepsPsiOffDiag(:)
    complex                                       :: grPsiHepsPsiNat(3, 3)
    complex                                       :: grPsiHepsPsi(3, 3)
    complex                                       :: psiHepsgrPsi(3, 3)
    complex                                       :: psiHepsPsi(3, 3)
    complex                                       :: psiGrVeff0sphPsi(3, 3)
    real                                          :: kExt(3)

    ! We do not want the local coordinate systems to be rotated for non-representative atoms constructing the matching coefficients.
    allocate(ngoprI(atoms%nat))
    ngoprI(:) = 1

    allocate( gBasExt(MAXVAL(nv), 3) )
    ! Small matching coefficients
    allocate( a( MAXVAL(nv), 0:atoms%lmaxd*(atoms%lmaxd+2), atoms%nat), b(MAXVAL(nv), 0:atoms%lmaxd*(atoms%lmaxd+2), atoms%nat), &
            & bascof_lo(3, -atoms%llod:atoms%llod, 4 * atoms%llod + 2, atoms%nlod, atoms%nat) )
    allocate( aKpq( MAXVAL(nv), 0:atoms%lmaxd*(atoms%lmaxd+2), atoms%nat), bKpq(MAXVAL(nv), 0:atoms%lmaxd*(atoms%lmaxd+2), atoms%nat), &
            & bascof_loKpq(3, -atoms%llod:atoms%llod, 4 * atoms%llod + 2, atoms%nlod, atoms%nat) )
    ! Large matching coefficients
    allocate( ab0cofScl(lmpMax, atoms%nat), ab0cofVec(lmpMax), & !ab0KpGcof(lmpMax, 3, atoms%nat), &
            & ab1cofVec(lmpMax, 3, atoms%nat, atoms%nat), abGrPsiSumCof(lmpMax, 3, atoms%nat, atoms%nat) )
    ! Temporary help arrays
    allocate( grPsiHepsVarphiOffDiagNat(lmpMax, 3, atoms%nat), varphiHepsGrPsiOffDiagNat(lmpMax, 3, atoms%nat), &
            & allMTvarphiHepsPsi(lmpMax, atoms%nat), grVarphiHepsPsiOffDiagNat(lmpMax), &
            & psiHepsgrVarphiNat(lmpMax), &
            & grVarphiHepsPsiNat(lmpMax), varphiHepsGrPsiOffDiag(lmpMax, 3), &!psiHepsVarphi(lmpMax, 3, atoms%nat), &
            & varphiHepsPsi(lmpMax, 3), abGrPsiSumCofDiag(lmpMax, 3, atoms%nat) )
    allocate( psiHepsgrPsiNat(3, 3), varphiHepsGrPsiNat(lmpMax), varphiGrVeff0sphPsi(lmpMax), varphiHepsPsiOffDiag(lmpMax))

      ikpq = mapKpq2K(ikpt, iqpt)

      ! Matching coefficients of MT basis functions at k
      nmat = nv(1, ikpt) + atoms%nlotot
      a(:, :, :) = cmplx(0.0, 0.0)
      b(:, :, :) = cmplx(0.0, 0.0)
      bascof_lo(:, :, :, :, :) = cmplx(0.0, 0.0)
      !call abcof3( atoms%lmaxd, atoms%ntype, atoms%nat, sym%nop, MAXVAL(nv), input%jspins, 1, atoms%lmaxd*(atoms%lmaxd+2), SIZE(z(:,1,1,1)), &
        !& atoms%llod, atoms%nlod, atoms%nlotot, sym%invtab, atoms%ntype, sym%mrot, ngoprI, atoms%taual, atoms%neq, atoms%lmax, &
        !& atoms%rmt, cell%omtil, cell%bmat, cell%bbmat, kpts%bk(:, ikpt), gBas(1, gBasUnwrap(:nv(1, ikpt), ikpt, 1)), &
        !& gBas(2, gBasUnwrap(:nv(1, ikpt), ikpt, 1)), gBas(3, gBasUnwrap(:nv(1, ikpt), ikpt, 1)), nv(:, ikpt), nmat, &
        !& usdus%us, usdus%dus, usdus%uds, usdus%duds, usdus%ddn, sym%invsat, sym%invsatnr, usdus%ulos, usdus%uulon, usdus%dulon, &
        !& usdus%dulos, atoms%llo, atoms%nlo, atoms%l_dulo, atoms%lapw_l, kveclo(:,ikpt), odi, ods, a, b, bascof_lo )

        nk=fmpi%k_list(ikpt)
        CALL lapw%init(input, noco, nococonv, kpts, atoms, sym, nk, cell, .FALSE., fmpi)
        CALL abcof3(input, atoms, sym, 1, cell, kpts%bk(:, ikpt), lapw, &
                            usdus, oneD, 1, lapw%dim_nvd(), a, b, bascof_lo)

      ! Matching coefficients of MT basis functions at k + q = k' (ikpq). In the MT we must not account for the backfolding vector.
      nmat = nv(1, ikpq) + atoms%nlotot
      aKpq(:, :, :) = cmplx(0.0, 0.0)
      bKpq(:, :, :) = cmplx(0.0, 0.0)
      bascof_loKpq(:, :, :, :, :) = cmplx(0.0, 0.0)
      !call abcof3( atoms%lmaxd, atoms%ntype, atoms%nat, sym%nop, MAXVAL(nv), input%jspins, 1, atoms%lmaxd*(atoms%lmaxd+2), SIZE(z(:,1,1,1)), &
        !& atoms%llod, atoms%nlod, atoms%nlotot, sym%invtab, atoms%ntype, sym%mrot, ngoprI, atoms%taual, atoms%neq, atoms%lmax, &
        !& atoms%rmt, cell%omtil, cell%bmat, cell%bbmat, kpts%bk(:, ikpq), gBas(1, gBasUnwrap(:nv(1, ikpq), ikpq, 1)), &
        !& gBas(2, gBasUnwrap(:nv(1, ikpq), ikpq, 1)), gBas(3, gBasUnwrap(:nv(1, ikpq), ikpq, 1)), nv(:, ikpq), nmat, &
        !& usdus%us, usdus%dus, usdus%uds, usdus%duds, usdus%ddn, sym%invsat, sym%invsatnr, usdus%ulos, usdus%uulon, usdus%dulon, &
        !& usdus%dulos, atoms%llo, atoms%nlo, atoms%l_dulo, atoms%lapw_l, kveclo(:,ikpq), odi, ods, aKpq, bKpq, bascof_loKpq )

        nk=fmpi%k_list(ikpq)
        CALL lapw%init(input, noco, nococonv, kpts, atoms, sym, nk, cell, .FALSE., fmpi)
        CALL abcof3(input, atoms, sym, 1, cell, kpts%bk(:, ikpq), lapw, &
                            usdus, oneD, 1, lapw%dim_nvd(), aKpq, bKpq, bascof_loKpq)

      ! Quantities we need later as decoration for the matching coefficients
      kExt(1:3) = matmul(cell%bmat(1:3, 1:3), kpts%bk(1:3, ikpt))
      gBasExt(:, :) = 0.
      do iBas = 1, nv(1, ikpt)
        ! We unwrap the gBas vectors once, so that we have a linear memory run within the loops later.
        gBasExt(iBas, 1:3) = matmul( cell%bmat(1:3, 1:3), gBas(1:3, gBasUnwrap(iBas, ikpt, 1)) )
      end do ! iBas


      ! Note: We should review switching the iband and the atom A loops during parallelization. For serial execution this way is optimal.
      do iband = 1, nobd(ikpt, 1)

        ab0cofScl(:, :) = cmplx(0.0, 0.0)
!        ab0KpGcof(:, :, :) = cmplx(0.0, 0.0)
        ! todo seems to be obsolete
!        psiHepsVarphi(:, :, :) = cmplx(0.0, 0.0)
        ab1cofVec(:, :, :, :) = cmplx(0.0, 0.0)
        allMTvarphiHepsPsi(:, :) = cmplx(0.0, 0.0)
        abGrPsiSumCof(:, :, :, :) = cmplx(0., 0.)
        iDatomA = 0
        do iDtypeA = 1, atoms%ntype
          do iDeqatA = 1, atoms%neq(iDtypeA)
            iDatomA = iDatomA + 1
            ! Calculate the large matching coefficients for all atoms or atom combinations before setting up the Pulay matrix element
            ! with its complicated atom combinations preventing a linear run through all atoms by one loop.
            lm  = 0
            lmp = 0
            do oqn_l = 0, atoms%lmax(iDtypeA)
              do mqn_m = - oqn_l, oqn_l
                pMaxLocal = nRadFun(oqn_l, iDtypeA)
                ! p = 1
                ab0cofScl(lmp + 1, iDatomA) = ImagUnit**oqn_l * dot_product(conjg(z(:nv(1, ikpt), iband, ikpt, 1)), a(:nv(1, ikpt), lm, iDatomA))
                ! p = 2
                ab0cofScl(lmp + 2, iDatomA) = ImagUnit**oqn_l * dot_product(conjg(z(:nv(1, ikpt), iband, ikpt, 1)), b(:nv(1, ikpt), lm, iDatomA))
                ! Add LO contributions
                do iradf = 3, pMaxLocal
                  ! p = 1
                  ab0cofScl(lmp + 1, iDatomA) = ab0cofScl(lmp + 1, iDatomA) &
                    & + ImagUnit**oqn_l * dot_product( conjg(z(nv(1, ikpt) + 1 :nv(1, ikpt) + atoms%nlotot, iband, ikpt, 1)), &
                                                        & bascof_lo(1, mqn_m, :atoms%nlotot, iloTable(iradf, oqn_l, iDtypeA), iDatomA) )
                  ! p = 2
                  ab0cofScl(lmp + 2, iDatomA) = ab0cofScl(lmp + 2, iDatomA) &
                    & + ImagUnit**oqn_l * dot_product( conjg(z(nv(1, ikpt) + 1 :nv(1, ikpt) + atoms%nlotot, iband, ikpt, 1)), &
                                                        & bascof_lo(2, mqn_m, :atoms%nlotot, iloTable(iradf, oqn_l, iDtypeA), iDatomA) )
                  ! 2 < p < LOs for that l and that atom type
                  ab0cofScl(lmp + iradf, iDatomA) = ImagUnit**oqn_l * dot_product( conjg(z(nv(1, ikpt) + 1 :nv(1, ikpt) + atoms%nlotot, iband, ikpt, 1)), &
                                                        & bascof_lo(3, mqn_m, :atoms%nlotot, iloTable(iradf, oqn_l, iDtypeA), iDatomA) )
                end do ! iradf
                lm = lm + 1
                lmp = lmp + pMaxLocal
              end do ! mqn_m
            end do ! oqn_l

            do idirR = 1, 3
              lm =  0
              lmp = 0
              ab0cofVec(:) = cmplx(0., 0.)
              do oqn_l = 0, atoms%lmax(iDtypeA)
                do mqn_m = -oqn_l, oqn_l
                  pMaxLocal = nRadFun(oqn_l, iDtypeA)

                  ! p = 1
                  ab0cofVec(lmp + 1) = &
                    &ImagUnit**oqn_l * dot_product(conjg(z(:nv(1, ikpt), iband, ikpt, 1)), gBasExt(:nv(1, ikpt), idirR) * a(:nv(1, ikpt), lm, iDatomA))
                  ! p = 2
                  ab0cofVec(lmp + 2) = &
                    &ImagUnit**oqn_l * dot_product(conjg(z(:nv(1, ikpt), iband, ikpt, 1)), gBasExt(:nv(1, ikpt), idirR) * b(:nv(1, ikpt), lm, iDatomA))
                  do iradf = 3, pMaxLocal
                    ! p = 1
                    ab0cofVec(lmp + 1) = ab0cofVec(lmp + 1) + ImagUnit**oqn_l * dot_product(conjg(z(nv(1, ikpt) + 1 :nv(1, ikpt) + atoms%nlotot, iband, ikpt, 1)), &
                                                      &   gBasExt(kveclo(:atoms%nlotot, ikpt), idirR) &
                                                      & * bascof_lo(1, mqn_m, :atoms%nlotot, iloTable(iradf, oqn_l, iDtypeA), iDatomA) )
                    ! p = 2
                    ab0cofVec(lmp + 2) = ab0cofVec(lmp + 2) + ImagUnit**oqn_l * dot_product(conjg(z(nv(1, ikpt) + 1 :nv(1, ikpt) + atoms%nlotot, iband, ikpt, 1)), &
                                                      &   gBasExt(kveclo(:atoms%nlotot, ikpt), idirR) &
                                                      & * bascof_lo(2, mqn_m, :atoms%nlotot, iloTable(iradf, oqn_l, iDtypeA), iDatomA) )
                    ! 2 < p < LOs for that l and that atom type
                    ab0cofVec(lmp + iradf) = ImagUnit**oqn_l * dot_product(conjg(z(nv(1, ikpt) + 1 :nv(1, ikpt) + atoms%nlotot, iband, ikpt, 1)), &
                                                      &   gBasExt(kveclo(:atoms%nlotot, ikpt), idirR) &
                                                      & * bascof_lo(3, mqn_m, :atoms%nlotot, iloTable(iradf, oqn_l, iDtypeA), iDatomA) )
                  end do ! iradf

                  lm = lm + 1
                  lmp = lmp + nRadFun(oqn_l, iDtypeA)
                end do ! mqn_m
              end do ! oqn_l
              ! This is supposed to be the sum of abcofs which belong to
              abGrPsiSumCof(1:lmpT(iDtypeA), idirR, iDatomA, iDatomA) = &
                                           & ImagUnit * ( kExt(idirR) * ab0cofScl(1:lmpT(iDtypeA), iDatomA) + ab0cofVec(1:lmpT(iDtypeA)) )
            end do ! idirR
            ! Diagonal in atoms
            abGrPsiSumCofDiag(:, :, :) = cmplx(0., 0.)
            abGrPsiSumCofDiag(1:lmpT(iDtypeA), 1:3, iDatomA) = abGrPsiSumCof(1:lmpT(iDtypeA), 1:3, iDatomA, iDatomA)

            iDatomB = 0
            do iDtypeB = 1, atoms%ntype
              do iDeqatB = 1, atoms%neq(iDtypeB)
                iDatomB = iDatomB + 1
                do idirR = 1, 3
                  lm  = 0
                  lmp = 0
                  do oqn_l = 0, atoms%lmax(iDtypeB)
                    do mqn_m = -oqn_l, oqn_l
                      ! This has to be iDtypeB, as the MT basis (small) matching coefficients are dependent on iDtypeB
                      pMaxLocal = nRadFun(oqn_l, iDtypeB)
                      ! This quantity has also to be calculated for all atom combinations before proceeding as we do here by closing the
                      ! atom alpha loops.
                      ab1cofVec(lmp + 1, idirR, iDatomB, iDatomA) = &
                                         & ImagUnit**oqn_l * dot_product(conjg(z1nG(:nv(1, ikpq), idirR, iDatomA, iband)), aKpq(:nv(1, ikpq), lm, iDatomB))
                      ab1cofVec(lmp + 2, idirR, iDatomB, iDatomA) = &
                                         & ImagUnit**oqn_l * dot_product(conjg(z1nG(:nv(1, ikpq), idirR, iDatomA, iband)), bKpq(:nv(1, ikpq), lm, iDatomB))
                      do iradf = 3, pMaxLocal
                        ab1cofVec(lmp + 1, idirR, iDatomB, iDatomA) = ab1cofVec(lmp + 1, idirR, iDatomB, iDatomA) &
                                     & + ImagUnit**oqn_l * dot_product( conjg( z1nG(nv(1, ikpq) + 1 :nv(1, ikpq) + atoms%nlotot, idirR, iDatomA, iband) ), &
                                                     & bascof_loKpq(1, mqn_m, :atoms%nlotot, iloTable(iradf, oqn_l, iDtypeB), iDatomB) )
                        ab1cofVec(lmp + 2, idirR, iDatomB, iDatomA) = ab1cofVec(lmp + 2, idirR, iDatomB, iDatomA) &
                                     & + ImagUnit**oqn_l * dot_product( conjg( z1nG(nv(1, ikpq) + 1 :nv(1, ikpq) + atoms%nlotot, idirR, iDatomA, iband) ), &
                                                     & bascof_loKpq(2, mqn_m, :atoms%nlotot, iloTable(iradf, oqn_l, iDtypeB), iDatomB) )
                        ab1cofVec(lmp + iradf, idirR, iDatomB, iDatomA) = &
                                     &   ImagUnit**oqn_l * dot_product( conjg( z1nG(nv(1, ikpq) + 1 :nv(1, ikpq) + atoms%nlotot, idirR, iDatomA, iband) ), &
                                                     & bascof_loKpq(3, mqn_m, :atoms%nlotot, iloTable(iradf, oqn_l, iDtypeB), iDatomB) )
                      end do ! iradf
                      ! This is a precalculation for the quantity evaluated in the MTs of the whole unit cell (1st term of A.53 in
                      ! PhD thesis A. Klüppelberg. Actually we would loop over a seperate type variable not being iDtypeB. But here,
                      ! for sake of performance, we use the given loop structure.
                      lm = lm + 1
                      lmp = lmp + nRadFun(oqn_l, iDtypeB)
                    end do ! mqn_m
                  end do ! oqn_l
                  if (eps1DynMatPulTestSw) then
                    abGrPsiSumCof(1:lmpT(iDtypeA), idirR, iDatomB, iDatomA) = ab1cofVec(1:lmpT(iDtypeA), idirR, iDatomB, iDatomA)
                  else
                    abGrPsiSumCof(1:lmpT(iDtypeA), idirR, iDatomB, iDatomA) = abGrPsiSumCof(1:lmpT(iDtypeA), idirR, iDatomB, iDatomA) &
                                                                             & + ab1cofVec(1:lmpT(iDtypeA), idirR, iDatomB, iDatomA)
                  end if
                end do ! idirR
              end do ! iDeqatB
            end do ! iDtypeB
            ! ATTENTION: ab0KpGcof, ab0cofScl, grPsiHepsVarphiOffDiagNat, psiHepsVarphi, ab1cofVec and allMTvarphiHepsPsi require the two
            !            following A loops to be closed to be available for all atoms later. This way might not be optimal for systems
            !            with a larger number of atoms than three and should be reviewed when parallelizing the code.
          end do ! iDeqatA
        end do ! iDtypeA

        iDatomA = 0
        varphiHepsPsi(:, :) = cmplx(0., 0.)
        do iDtypeA = 1, atoms%ntype
          do iDeqatA = 1, atoms%neq(iDtypeA)
            iDatomA = iDatomA + 1
            do idirC = 1, 3
            !todo better docs!      ! Calculate -eps <gradVarphi|varphi>* + HkinSelfAdjointCorrection
              ! ( = <varphi|H - eps|gradvarphi> - <grad varphi|H|varphi>* ) in the 3rd term of PhD A. Klueppelberg.

              ! Zero if eps1DynMatPulTestSw because of abGrPsiSumCofDiag
              varphiHepsPsi(1:lmpT(iDtypeA), idirC) = matmul(varphiHvarphi(1:lmpT(iDtypeA), 1:lmpT(iDtypeA), iDatomA),             &
                                                                       & abGrPsiSumCofDiag( 1:lmpT(iDtypeA), idirC, iDatomA ) )    &
                              - eig(iband, ikpt, 1) * matmul(varphiVarphi(1:lmpT(iDtypeA), 1:lmpT(iDtypeA), iDtypeA),              &
                                                                            & abGrPsiSumCofDiag( 1:lmpT(iDtypeA), idirC, iDatomA ) )
            end do ! idirC

            ! Loop over all MTs gamma, not only the displaced ones.
            iDatomB = 0
            do iDtypeB = 1, atoms%ntype
              do iDeqatB = 1, atoms%neq(iDtypeB)
                iDatomB = iDatomB + 1
                psiHepsPsi(:, :) = cmplx(0., 0.)
                varphiHepsPsiOffDiag(:) = cmplx(0., 0.)
                grPsiHepsPsiNat(:, :) = cmplx(0., 0.)
                psiHepsGrPsiNat(:, :) = cmplx(0., 0.)
                psiGrVeff0sphPsi(:, :) = cmplx(0., 0.)

                ! This is a precalculation of the off-diagonal term in the 5th line of A.53 in the PhD thesis A. Klueppelberg. The dyadic
                ! product is expressed by an adequate choice of the idir loop variables.
                do idirC = 1, 3
                  varphiHepsPsiOffDiag(1:lmpT(iDtypeB)) = matmul(varphiHvarphi(1:lmpT(iDtypeB), 1:lmpT(iDtypeB), iDtypeB),         &
                                                                             & ab1cofVec(1:lmpT(iDtypeB), idirC, iDatomB, iDatomA) )&
                                  - eig(iband, ikpt, 1) * matmul(varphiVarphi(1:lmpT(iDtypeB), 1:lmpT(iDtypeB), iDtypeB),          &
                                                                             & ab1cofVec(1:lmpT(iDtypeB), idirC, iDatomB, iDatomA) )
                  do idirR = 1, 3
                    grVarphiHepsPsiNat(:) = cmplx(0., 0.)
!                    varphiHepsGrPsiNat(:) = cmplx(0., 0.)
                    ! This term counts if eps1DynMatPulTestSw
                    if (.not.testCompTerm3rdBraKetsVarBra ) then
                      psiHepsPsi(idirR, idirC) = dot_product(ab1cofVec(1:lmpT(iDtypeA), idirR, iDatomA, iDatomB),             &
                                                                                           & varphiHepsPsi(1:lmpT(iDtypeA), idirC) )
                    end if
                    ! This term counts if eps1DynMatPulTestSw
                    !if ( .not.testCompTerm3rdBraKetsVarKet ) then
                    !psiHepsPsi(idirR, idirC) = psiHepsPsi(idirR, idirC) +                                        &
                    !                                  & dot_product(abGrPsiSumCofDiag(1:lmpT(iDtypeB), idirR, iDatomB),            &
                    !                                                                       & varphiHepsPsiOffDiag(1:lmpT(iDtypeB)))
                    !end if
                    if (iDatomA == iDatomB .and. .not.eps1DynMatPulTestSw ) then
                      psiHepsPsi(idirR, idirC) = psiHepsPsi(idirR, idirC) +                                      &
                           & dot_product(abGrPsiSumCofDiag(1:lmpT(iDtypeA), idirR, iDatomA), varphiHepsPsi(1:lmpT(iDtypeA), idirC))
                    end if

                    ! This can even be made more performant by adding the overlap with the Hamiltonian before multiplying with the
                    ! ab1cofVec
                    varphiHepsGrPsiNat(:) = cmplx(0., 0.)
                    varphiHepsGrPsiNat(1:lmpT(iDtypeA)) = matmul(varphiHGrvarphiNat(1:lmpT(iDtypeA), 1:lmpT(iDtypeA), iDatomA, idirC), ab0cofScl(:lmpT(iDtypeA), iDatomA)) &
                    & - eig(iband, ikpt, 1) * matmul(varphiGrVarphiNat(1:lmpT(iDtypeA), 1:lmpT(iDtypeA), idirC, iDtypeA), ab0cofScl(:lmpT(iDtypeA), iDatomA))


                    psiHepsGrPsiNat(idirR, idirC) = dot_product( abGrPsiSumCof(1:lmpT(iDtypeA), idirR, iDatomA, iDatomB), varphiHepsGrPsiNat(1:lmpT(iDtypeA)) )

                    varphiGrVeff0sphPsi(:) = cmplx(0., 0.)
                    varphiGrVeff0sphPsi(1:lmpT(iDtypeB)) = matmul(varphiGrVeff0SphVarphi(1:lmpT(iDtypeA), 1:lmpT(iDtypeA), iDatomA, idirC), ab0cofScl(1:lmpT(iDtypeA), iDatomA))
                    psiGrVeff0sphPsi(idirR, idirC) = dot_product(abGrPsiSumCof(1:lmpT(iDtypeA), idirR, iDatomA, iDatomB), varphiGrVeff0sphPsi(1:lmpT(iDtypeA)))



                    ! Calculation of the diagonal terms in the 5th line of the result in A.53  in the PhD thesis of A. Klueppelberg.
                    grVarphiHepsPsiNat(1:lmpT(iDtypeB)) = &
                      & matmul( grVarphiHvarphiNat(1:lmpT(iDtypeB), 1:lmpT(iDtypeB), idirR, iDatomB), &
                                                                      & abGrPsiSumCof(1:lmpT(iDtypeB), idirC, iDatomB, iDatomA) )  &
                      & - eig(iband, ikpt, 1) * matmul( grVarphiVarphiNat(1:lmpT(iDtypeB), 1:lmpT(iDtypeB), idirR, iDtypeB),    &
                      & abGrPsiSumCof(1:lmpT(iDtypeB), idirC, iDatomB, iDatomA) )

                    grPsiHepsPsiNat(idirR, idirC) = dot_product(ab0cofScl(:lmpT(iDtypeB), iDatomB), &
                                                                                              & grVarphiHepsPsiNat(:lmpT(iDtypeB)) )
                  end do ! idirR
                end do ! idirC

                psiHepsgrPsi(:, :) = cmplx(0.0, 0.0)
                psiHepsgrPsi(1:3, 1:3) = matmul(psiHepsgrPsiNat(1:3, 1:3), transpose(Tmatrix0(1:3, 1:3)))
                psiHepsgrPsi(1:3, 1:3) = psiHepsgrPsi(1:3, 1:3) - psiGrVeff0sphPsi(1:3, 1:3)

                grPsiHepsPsi(:, :) = 0.
                grPsiHepsPsi(1:3, 1:3) = matmul(conjg(Tmatrix0(1:3, 1:3)), grPsiHepsPsiNat(1:3, 1:3))

                do idirC = 1, 3
                  iatomG = 0
                  do itypeG = 1, atoms%ntype
                    do ieqatG = 1, atoms%neq(itypeG)
                      iatomG = iatomG + 1
                      allMTvarphiHepsPsi(:lmpT(itypeG), iatomG) = &
                        & matmul( varphiHvarphi(1:lmpT(itypeG), 1:lmpT(itypeG), iatomG), &
                                                                              & ab1cofVec(:lmpT(itypeG), idirC, iatomG, iDatomA) ) &
                        & - eig(iband, ikpt, 1) * matmul( varphivarphi(1:lmpT(itypeG), 1:lmpT(itypeG), itypeG), &
                                                                              & ab1cofVec(:lmpT(itypeG), idirC, iatomG, iDatomA) )
                    end do ! ieqatG
                  end do ! itypeG

                  do idirR = 1, 3
                    wholeUnitCellIR = cmplx(0.0, 0.0)
                    wholeUnitCellMT = cmplx(0.0, 0.0)

                    if( .not.eps1DynMatPulTestSw .and. .not.testCompTerm3rdBraKetsVarKet .and. .not.testCompTerm3rdBraKetsVarBra ) then
                      ! todo solve the kimax problem down under does it start form 0 or 1 we see it at the comment 0 or 1
                      ! todo The shift of the Gsets might be wrong kpq2kPrVec is for both basis funcrtions so cancels away actually
                      ! therefore one can also leave it away
                      ! nmat is at k + q, as given before the last abcof3 call in this routine
                      call calcPsi1HepsPsi1IR( kpts, qpts, stars, cell, gBas(:, gBasUnwrap(:nv(1, ikpq), ikpq, 1)), nv, &
                        & z1nG(:, idirR, iDatomB, iband), z1nG(:, idirC, iDatomA, iband), nmat, ikpt, iqpt, ikpq, kpq2kPrVec, &
                        & vEff0IR, eig, iband, wholeUnitCellIR )

                      iatomG = 0
                      do itypeG = 1, atoms%ntype
                        do ieqatG = 1, atoms%neq(itypeG)
                          iatomG = iatomG + 1
                          wholeUnitCellMT = wholeUnitCellMT &
                            & + dot_product( ab1cofVec(:lmpT(itypeG), idirR, iatomG, iDatomB), &
                                                                         & allMTvarphiHepsPsi(:lmpT(itypeG), iatomG) )
                        end do ! ieqatG
                      end do ! itypeG
                    end if
                    if (testCompTerm3rdBraKetsVarKet) then
                      grPsiHepsPsi(:, :) = cmplx(0., 0.)
                    else if (testCompTerm3rdBraKetsVarBra) then
                      psiHepsgrPsi(:, :) = cmplx(0., 0.)
                    end if

                    if (dynMatPu3rdBraKetHepsSw) then
                      psiHepsPsi(:, :) = cmplx(0., 0.)
                      psiHepsgrPsi(:, :) = cmplx(0., 0.)
                      grPsiHepsPsi(:, :) = cmplx(0., 0.)
                    end if

                    ! 4 = 2 * 2, degenerated spin and mixed term before <Psi1 |H - eps| Psi1>
                    dynMatPu(idirR + (iDatomB - 1) * 3, idirC + (iDatomA - 1) * 3) = &
                      & dynMatPu(idirR + (iDatomB - 1) * 3, idirC + (iDatomA - 1) * 3) + 4 * results%w_iks(iband, ikpt, 1) * ( &

                      !&  wholeUnitCellIR &
                      !& + wholeUnitCellMT  &
                      & + psiHepsPsi(idirR, idirC) &
                 !     & - psiHepsgrPsi(idirR, idirC) &
                 !     &- grPsiHepsPsi(idirR, idirC)&
                      &)
                    !write(465,*) 'Pu 1      ', idirR, idirC, iband, ikpt
                    !write(465,*) 4 * results%w_iks(iband, ikpt, 1) * ( wholeUnitCellIR + wholeUnitCellMT  + psiHepsPsi(idirR, idirC) )
                  end do ! idirR
                end do ! idirC
              end do ! iDeqatB
            end do ! iDtypeB
          end do ! iDeqatA
        end do ! iDtypeA
      end do ! iband

  end subroutine Add1stOrdWfPulayBraKets2DynMat

  subroutine AddAlexPulayBraKets2DynMat( fmpi, noco, nococonv, oneD, atoms, input, kpts, qpts, sym, cell, usdus, stars, results, ikpt, iqpt, lmpMax, nRadFunMax, nv, gBas, &
      & gBasUnwrap, kveclo, mapKpq2K, nobd, z, z1nG, iloTable, grVarphiVarphiNat, nRadFun, eig, kpq2kPrVec, &
      & grVarphiHvarphiNat, varphiVarphi, varphiHvarphi, vEff0IR, lmpT, eps1DynMatPulTestSw, &
      & testCompTerm3rdBraKetsVarBra, testCompTerm3rdBraKetsVarKet, dynMatPu3rdBraKetHepsSw, varphiGrVeff0SphVarphi, varphiHGrvarphiNat, varphiGrVarphiNat, dynMatPu )

    use m_types
    use m_abcof3

    implicit none

    ! Type parameters
    type(t_mpi),                  intent(in)  :: fmpi
    type(t_nococonv),                intent(in)  :: nococonv
    type(t_noco),                intent(in)  :: noco
    type(t_oneD),                intent(in)  :: oneD
    type(t_atoms),                  intent(in)    :: atoms
    type(t_input),                  intent(in)    :: input
    type(t_kpts),                   intent(in)    :: kpts
    type(t_kpts),                   intent(in)    :: qpts
    type(t_sym),                    intent(in)    :: sym
    type(t_cell),                   intent(in)    :: cell
    type(t_usdus),                  intent(in)    :: usdus
    type(t_stars),                  intent(in)    :: stars
    type(t_results),                intent(in)    :: results

    ! Scalar parameters
    integer,                        intent(in)    :: ikpt
    integer,                        intent(in)    :: iqpt
    integer,                        intent(in)    :: lmpMax
    integer,                        intent(in)    :: nRadFunMax
    logical,                        intent(in)    :: eps1DynMatPulTestSw
    logical,                        intent(in)    :: testCompTerm3rdBraKetsVarBra
    logical,                        intent(in)    :: testCompTerm3rdBraKetsVarKet
    logical,                        intent(in)    :: dynMatPu3rdBraKetHepsSw

    ! Array parameters
    integer,                        intent(in)    :: nv(:, :)
    integer,                        intent(in)    :: lmpT(:)
    integer,                        intent(in)    :: gBas(:, :)
    integer,                        intent(in)    :: gBasUnwrap(:, :, :)
    integer,                        intent(in)    :: kveclo(:,:)
    integer,                        intent(in)    :: mapKpq2K(:, :)
    integer,                        intent(in)    :: nobd(:, :)
    complex,                       intent(in)    :: z(:,:,:,:)
    complex,                        intent(in)    :: z1nG(:, :, :, :)
    integer,                        intent(in)    :: iloTable(:, 0:, :)
    real,                           intent(in)    :: grVarphiVarphiNat(:, :, :, :)
    integer,                        intent(in)    :: nRadFun(0:, :)
    real,                           intent(in)    :: eig(:, :, :)
    integer,                        intent(in)    :: kpq2kPrVec(:, :, :)
    complex,                        intent(in)    :: grVarphiHvarphiNat(:, :, :, :)
    real,                           intent(in)    :: varphiVarphi(:, :, :)
    complex,                        intent(in)    :: varphiHvarphi(:, :, :)
    complex,                        intent(in)    :: varphiGrVeff0SphVarphi(:, :, :, :)
    complex,                        intent(in)    :: varphiHGrvarphiNat(:, :, :, :)
    complex,                        intent(in)    :: vEff0IR(:,:)
    real,                           intent(in)    :: varphiGrVarphiNat(:, :, :, :)
    complex,                        intent(inout) :: dynMatPu(:, :)

    ! Type variable
    type(od_inp)                                  :: odi
    type(od_sym)                                  :: ods
    type(t_lapw) :: lapw

    ! Scalar variables
    integer                                       :: nmat
    integer                                       :: ikpq
    integer                                       :: iBas
    integer                                       :: iband
    integer                                       :: iDatomA
    integer                                       :: iDtypeA
    integer                                       :: iDeqatA
    integer                                       :: lm
    integer                                       :: lmp
    integer                                       :: oqn_l
    integer                                       :: mqn_m
    integer                                       :: iradf
    integer                                       :: idirC
    integer                                       :: lmp1Pr
    integer                                       :: iatomG
    integer                                       :: itypeG
    integer                                       :: ieqatG
    integer                                       :: idirR
    integer                                       :: oqn_l1Pr
    integer                                       :: mqn_m1Pr
    integer                                       :: pMaxLocal
    complex                                       :: wholeUnitCellIR
    complex                                       :: wholeUnitCellMT
    integer                                       :: iDatomB
    integer                                       :: iDtypeB
    integer                                       :: iDeqatB
    integer                                       :: iradf1Pr
    complex                                       :: sAdjCorrPsiGrPsiNat
    integer :: nk

    ! Array variabels
    integer,           allocatable                :: ngoprI(:)
    complex,           allocatable                :: a(:, :, :)
    complex,           allocatable                :: b(:, :, :)
    complex,           allocatable                :: bascof_lo(:, :, :, :, :)
    complex,           allocatable                :: aKpq(:, :, :)
    complex,           allocatable                :: bKpq(:, :, :)
    complex,           allocatable                :: bascof_loKpq(:, :, :, :, :)
    real,              allocatable                :: gBasExt(:, :)
    real,              allocatable                :: gqBasExt(:, :)
    complex,           allocatable                :: ab0cofScl(:, :)
    complex,           allocatable                :: ab0cofVec(:)
    !complex,           allocatable                :: ab0KpGcof(:, :, :) seems obsolete
    complex,           allocatable                :: grPsiHepsVarphiOffDiagNat(:, :, :)
    complex,           allocatable                :: varphiHepsGrPsiOffDiagNat(:, :, :)
    complex,           allocatable                :: allMTvarphiHepsPsi(:, :)
    complex,           allocatable                :: ab1cofVec(:, :, :, :)
    complex,           allocatable                :: grVarphiHepsPsiOffDiagNat(:)
    complex,           allocatable                :: psiHepsgrVarphiNat(:)
    complex,           allocatable                :: grVarphiHepsPsiNat(:)
    complex,           allocatable                :: varphiHepsGrPsiOffDiag(:, :)
    !complex,           allocatable                :: psiHepsVarphi(:, :, :) ! seems obsolete
    complex,           allocatable                :: abGrPsiSumCof(:, :, :, :)
    complex,           allocatable                :: psiHepsgrPsiNat(:, :)
    complex,           allocatable                :: varphiHepsGrPsiNat(:)
    complex,           allocatable                :: varphiGrVeff0sphPsi(:)
    complex,           allocatable                :: varphiHepsPsi(:, :)
    complex,           allocatable                :: abGrPsiSumCofDiag(:, :, :)
    complex,           allocatable                :: varphiHepsPsiOffDiag(:)
    complex,           allocatable                :: bigAscal(:,:), bigAvec(:,:,:), bigAmat(:,:,:,:), bigAz1vec(:,:,:,:), bigAz1mat(:,:,:,:,:)
    complex,           allocatable                :: helpscal(:), helpvec(:,:), helpmat(:,:,:)
    complex                                       :: grPsiHepsPsiNat(3, 3)
    complex                                       :: grPsiHepsPsi(3, 3)
    complex                                       :: psiHepsgrPsi(3, 3)
    complex                                       :: psiHepsPsi(3, 3)
    complex                                       :: psiGrVeff0sphPsi(3, 3)
    real                                          :: kExt(3)

    ! We do not want the local coordinate systems to be rotated for non-representative atoms constructing the matching coefficients.
    allocate(ngoprI(atoms%nat))
    ngoprI(:) = 1

    allocate( gBasExt(MAXVAL(nv), 3) )
    allocate( gqBasExt(MAXVAL(nv), 3) )
    ! Small matching coefficients
    allocate( a( MAXVAL(nv), 0:atoms%lmaxd*(atoms%lmaxd+2), atoms%nat), b(MAXVAL(nv), 0:atoms%lmaxd*(atoms%lmaxd+2), atoms%nat), &
            & bascof_lo(3, -atoms%llod:atoms%llod, 4 * atoms%llod + 2, atoms%nlod, atoms%nat) )
    allocate( aKpq( MAXVAL(nv), 0:atoms%lmaxd*(atoms%lmaxd+2), atoms%nat), bKpq(MAXVAL(nv), 0:atoms%lmaxd*(atoms%lmaxd+2), atoms%nat), &
            & bascof_loKpq(3, -atoms%llod:atoms%llod, 4 * atoms%llod + 2, atoms%nlod, atoms%nat) )
    ! Large matching coefficients
    allocate( ab0cofScl(lmpMax, atoms%nat), ab0cofVec(lmpMax), & !ab0KpGcof(lmpMax, 3, atoms%nat), &
            & ab1cofVec(lmpMax, 3, atoms%nat, atoms%nat), abGrPsiSumCof(lmpMax, 3, atoms%nat, atoms%nat) )
    ! Large matchigs for Alex
    allocate( bigAscal(lmpMax, atoms%nat), bigAvec(lmpMax, 3, atoms%nat), bigAmat(lmpMax, 3, 3, atoms%nat), &
            & bigAz1vec(lmpMax, 3, atoms%nat, atoms%nat), bigAz1mat(lmpMax, 3, 3, atoms%nat, atoms%nat) )
    ! Help Alex
    allocate( helpscal(lmpMax) )!, helpvec(lmpMax, 3), helpmat(lmpMax, 3, 3) )
    ! Temporary help arrays
    allocate( grPsiHepsVarphiOffDiagNat(lmpMax, 3, atoms%nat), varphiHepsGrPsiOffDiagNat(lmpMax, 3, atoms%nat), &
            & allMTvarphiHepsPsi(lmpMax, atoms%nat), grVarphiHepsPsiOffDiagNat(lmpMax), &
            & psiHepsgrVarphiNat(lmpMax), &
            & grVarphiHepsPsiNat(lmpMax), varphiHepsGrPsiOffDiag(lmpMax, 3), &!psiHepsVarphi(lmpMax, 3, atoms%nat), &
            & varphiHepsPsi(lmpMax, 3), abGrPsiSumCofDiag(lmpMax, 3, atoms%nat) )
    allocate( psiHepsgrPsiNat(3, 3), varphiHepsGrPsiNat(lmpMax), varphiGrVeff0sphPsi(lmpMax), varphiHepsPsiOffDiag(lmpMax))

    ikpq = mapKpq2K(ikpt, iqpt)

    ! Matching coefficients of MT basis functions at k
    nmat = nv(1, ikpt) + atoms%nlotot
    a(:, :, :) = cmplx(0.0, 0.0)
    b(:, :, :) = cmplx(0.0, 0.0)
    bascof_lo(:, :, :, :, :) = cmplx(0.0, 0.0)
    !call abcof3( atoms%lmaxd, atoms%ntype, atoms%nat, sym%nop, MAXVAL(nv), input%jspins, 1, atoms%lmaxd*(atoms%lmaxd+2), SIZE(z(:,1,1,1)), &
     ! & atoms%llod, atoms%nlod, atoms%nlotot, sym%invtab, atoms%ntype, sym%mrot, ngoprI, atoms%taual, atoms%neq, atoms%lmax, &
      !& atoms%rmt, cell%omtil, cell%bmat, cell%bbmat, kpts%bk(:, ikpt), gBas(1, gBasUnwrap(:nv(1, ikpt), ikpt, 1)), &
      !& gBas(2, gBasUnwrap(:nv(1, ikpt), ikpt, 1)), gBas(3, gBasUnwrap(:nv(1, ikpt), ikpt, 1)), nv(:, ikpt), nmat, &
      !& usdus%us, usdus%dus, usdus%uds, usdus%duds, usdus%ddn, sym%invsat, sym%invsatnr, usdus%ulos, usdus%uulon, usdus%dulon, &
      !& usdus%dulos, atoms%llo, atoms%nlo, atoms%l_dulo, atoms%lapw_l, kveclo(:,ikpt), odi, ods, a, b, bascof_lo )
      nk=fmpi%k_list(ikpt)
      CALL lapw%init(input, noco, nococonv, kpts, atoms, sym, nk, cell, .FALSE., fmpi)
      CALL abcof3(input, atoms, sym, 1, cell, kpts%bk(:, ikpt), lapw, &
                          usdus, oneD, 1, lapw%dim_nvd(), a, b, bascof_lo)

    ! Matching coefficients of MT basis functions at k + q = k' (ikpq). In the MT we must not account for the backfolding vector.
    nmat = nv(1, ikpq) + atoms%nlotot
    aKpq(:, :, :) = cmplx(0.0, 0.0)
    bKpq(:, :, :) = cmplx(0.0, 0.0)
    bascof_loKpq(:, :, :, :, :) = cmplx(0.0, 0.0)
    !call abcof3( atoms%lmaxd, atoms%ntype, atoms%nat, sym%nop, MAXVAL(nv), input%jspins, 1, atoms%lmaxd*(atoms%lmaxd+2), SIZE(z(:,1,1,1)), &
     ! & atoms%llod, atoms%nlod, atoms%nlotot, sym%invtab, atoms%ntype, sym%mrot, ngoprI, atoms%taual, atoms%neq, atoms%lmax, &
     ! & atoms%rmt, cell%omtil, cell%bmat, cell%bbmat, kpts%bk(:, ikpq), gBas(1, gBasUnwrap(:nv(1, ikpq), ikpq, 1)), &
      !& gBas(2, gBasUnwrap(:nv(1, ikpq), ikpq, 1)), gBas(3, gBasUnwrap(:nv(1, ikpq), ikpq, 1)), nv(:, ikpq), nmat, &
      !& usdus%us, usdus%dus, usdus%uds, usdus%duds, usdus%ddn, sym%invsat, sym%invsatnr, usdus%ulos, usdus%uulon, usdus%dulon, &
      !& usdus%dulos, atoms%llo, atoms%nlo, atoms%l_dulo, atoms%lapw_l, kveclo(:,ikpq), odi, ods, aKpq, bKpq, bascof_loKpq )
      nk=fmpi%k_list(ikpq)
      CALL lapw%init(input, noco, nococonv, kpts, atoms, sym, nk, cell, .FALSE., fmpi)
      CALL abcof3(input, atoms, sym, 1, cell, kpts%bk(:, ikpq), lapw, &
                          usdus, oneD, 1, lapw%dim_nvd(), aKpq, bKpq, bascof_loKpq)


    ! Quantities we need later as decoration for the matching coefficients
    !kExt(1:3) = matmul(cell%bmat(1:3, 1:3), kpts%bk(1:3, ikpt))
    gBasExt(:, :) = 0.
    do iBas = 1, nv(1, ikpt)
      gBasExt(iBas, 1:3) = matmul( cell%bmat(1:3, 1:3), gBas(1:3, gBasUnwrap(iBas, ikpt, 1)) )
    end do
    gqBasExt(:, :) = 0.
    do iBas = 1, nv(1, ikpq)
      gqBasExt(iBas, 1:3) = matmul( cell%bmat(1:3, 1:3), gBas(1:3, gBasUnwrap(iBas, ikpq, 1)) + qpts%bk(1:3, iqpt) + kpq2kPrVec(1:3, ikpt, iqpt) )
    end do

    do iband = 1, nobd(ikpt, 1)
      bigAscal(:,:) = cmplx(0.0, 0.0)
      bigAvec(:,:,:) = cmplx(0.0, 0.0)
      bigAmat(:,:,:,:) = cmplx(0.0, 0.0)
      bigAz1vec(:,:,:,:) = cmplx(0.0, 0.0)
      bigAz1mat(:,:,:,:,:) = cmplx(0.0, 0.0)

      iDatomA = 0
      do iDtypeA = 1, atoms%ntype
        do iDeqatA = 1, atoms%neq(iDtypeA)
          iDatomA = iDatomA + 1
          ! Calculate the large matching coefficients for all atoms or atom combinations before setting up the Pulay matrix element
          ! with its complicated atom combinations preventing a linear run through all atoms by one loop.
          lm  = 0
          lmp = 0
          do oqn_l = 0, atoms%lmax(iDtypeA)
            do mqn_m = - oqn_l, oqn_l
              pMaxLocal = nRadFun(oqn_l, iDtypeA)
              bigAscal(lmp + 1, iDatomA) = ImagUnit**oqn_l * dot_product(conjg(z(:nv(1, ikpt), iband, ikpt, 1)), a(:nv(1, ikpt), lm, iDatomA))
              bigAscal(lmp + 2, iDatomA) = ImagUnit**oqn_l * dot_product(conjg(z(:nv(1, ikpt), iband, ikpt, 1)), b(:nv(1, ikpt), lm, iDatomA))
              lm = lm + 1
              lmp = lmp + pMaxLocal
            end do ! mqn_m
          end do ! oqn_l

          do idirR = 1, 3
            lm =  0
            lmp = 0
            do oqn_l = 0, atoms%lmax(iDtypeA)
              do mqn_m = -oqn_l, oqn_l
                pMaxLocal = nRadFun(oqn_l, iDtypeA)
                bigAvec(lmp + 1, idirR, iDatomA) = ImagUnit**oqn_l * ImagUnit * dot_product(conjg(z(:nv(1, ikpt), iband, ikpt, 1)), gBasExt(:nv(1, ikpt), idirR) * a(:nv(1, ikpt), lm, iDatomA))
                bigAvec(lmp + 2, idirR, iDatomA) = ImagUnit**oqn_l * ImagUnit * dot_product(conjg(z(:nv(1, ikpt), iband, ikpt, 1)), gBasExt(:nv(1, ikpt), idirR) * b(:nv(1, ikpt), lm, iDatomA))
                lm = lm + 1
                lmp = lmp + nRadFun(oqn_l, iDtypeA)
              end do ! mqn_m
            end do ! oqn_l
          end do ! idirR

          do idirC = 1, 3
            do idirR = 1, 3
              lm =  0
              lmp = 0
              do oqn_l = 0, atoms%lmax(iDtypeA)
                do mqn_m = -oqn_l, oqn_l
                  pMaxLocal = nRadFun(oqn_l, iDtypeA)
                  bigAmat(lmp + 1, idirR, idirC, iDatomA) = &
                         & ImagUnit**oqn_l *  dot_product(conjg(z(:nv(1, ikpt), iband, ikpt, 1)), gBasExt(:nv(1, ikpt), idirR) &
                                                                                            & * gBasExt(:nv(1, ikpt), idirC) * a(:nv(1, ikpt), lm, iDatomA))
                  bigAmat(lmp + 2, idirR, idirC, iDatomA) = &
                         & ImagUnit**oqn_l * dot_product(conjg(z(:nv(1, ikpt), iband, ikpt, 1)), gBasExt(:nv(1, ikpt), idirR) &
                                                                                            & * gBasExt(:nv(1, ikpt), idirC) * b(:nv(1, ikpt), lm, iDatomA))
                  lm = lm + 1
                  lmp = lmp + nRadFun(oqn_l, iDtypeA)
                end do ! mqn_m
              end do ! oqn_l
            end do ! idirR
          end do !idirC

          iDatomB = 0
          do iDtypeB = 1, atoms%ntype
            do iDeqatB = 1, atoms%neq(iDtypeB)
              iDatomB = iDatomB + 1
              do idirR = 1, 3
                lm  = 0
                lmp = 0
                do oqn_l = 0, atoms%lmax(iDtypeA)
                  do mqn_m = -oqn_l, oqn_l
                    pMaxLocal = nRadFun(oqn_l, iDtypeA)
                    bigAz1vec(lmp + 1, idirR, iDatomA, iDatomB) = &
                                       & ImagUnit**oqn_l * dot_product(conjg(z1nG(:nv(1, ikpq), idirR, iDatomB, iband)), aKpq(:nv(1, ikpq), lm, iDatomA))
                    bigAz1vec(lmp + 2, idirR, iDatomA, iDatomB) = &
                                       & ImagUnit**oqn_l * dot_product(conjg(z1nG(:nv(1, ikpq), idirR, iDatomB, iband)), bKpq(:nv(1, ikpq), lm, iDatomA))
                    lm = lm + 1
                    lmp = lmp + nRadFun(oqn_l, iDtypeA)
                  end do ! mqn_m
                end do ! oqn_l
              end do ! idirR
              do idirC = 1, 3
                do idirR = 1, 3
                  lm  = 0
                  lmp = 0
                  do oqn_l = 0, atoms%lmax(iDtypeA)
                    do mqn_m = -oqn_l, oqn_l
                      pMaxLocal = nRadFun(oqn_l, iDtypeA)
                      bigAz1mat(lmp + 1, idirR, idirC, iDatomA, iDatomB) = &
                                         & ImagUnit**oqn_l * ImagUnit * dot_product(conjg(z1nG(:nv(1, ikpq), idirR, iDatomB, iband)), &
                                                                      & gqBasExt(:nv(1, ikpq), idirC) * aKpq(:nv(1, ikpq), lm, iDatomA))
                      bigAz1mat(lmp + 2, idirR, idirC, iDatomA, iDatomB) = &
                                         & ImagUnit**oqn_l * ImagUnit * dot_product(conjg(z1nG(:nv(1, ikpq), idirR, iDatomB, iband)), &
                                                                      & gqBasExt(:nv(1, ikpq), idirC) * bKpq(:nv(1, ikpq), lm, iDatomA))
                      lm = lm + 1
                      lmp = lmp + nRadFun(oqn_l, iDtypeA)
                    end do ! mqn_m
                  end do ! oqn_l
                end do ! idirR
              end do ! idirC
            end do ! iDeqatB
          end do ! iDtypeB
        end do ! iDeqatA
      end do ! iDtypeA

      iDatomA = 0
      do iDtypeA = 1, atoms%ntype
        do iDeqatA = 1, atoms%neq(iDtypeA)
          iDatomA = iDatomA + 1
          iDatomB = 0
          do iDtypeB = 1, atoms%ntype
            do iDeqatB = 1, atoms%neq(iDtypeB)
              iDatomB = iDatomB + 1
              psiHepsPsi(:, :) = cmplx(0., 0.)

              do idirC = 1, 3
                do idirR = 1, 3

                  helpScal(1:lmpT(iDtypeA)) = matmul( varphiHvarphi(1:lmpT(iDtypeA), 1:lmpT(iDtypeA), iDtypeA), bigAvec(1:lmpT(iDtypeA), idirC, iDatomA) ) &
                    & - eig(iband, ikpt, 1) * matmul( varphiVarphi(1:lmpT(iDtypeA), 1:lmpT(iDtypeA), iDtypeA), bigAvec(1:lmpT(iDtypeA), idirC, iDatomA) )
                  psiHepsPsi(idirR, idirC) = dot_product( bigAz1vec(1:lmpT(iDtypeA), idirR, iDatomA, iDatomB), helpScal(1:lmpT(iDtypeA)) )

                  write(475,*) 'z1 ikGz0', idirR, idirC
                  write(475,*) dot_product( bigAz1vec(1:lmpT(iDtypeA), idirR, iDatomA, iDatomB), helpScal(1:lmpT(iDtypeA)) )

                  helpScal(1:lmpT(iDtypeA)) = matmul( varphiHvarphi(1:lmpT(iDtypeA), 1:lmpT(iDtypeA), iDtypeA), bigAscal(1:lmpT(iDtypeA), iDatomA) ) &
                    & - eig(iband, ikpt, 1) * matmul( varphiVarphi(1:lmpT(iDtypeA), 1:lmpT(iDtypeA), iDtypeA), bigAscal(1:lmpT(iDtypeA), iDatomA) )
                  psiHepsPsi(idirR, idirC) = psiHepsPsi(idirR, idirC) + &
                                           & dot_product( bigAz1mat(1:lmpT(iDtypeA), idirR, idirC, iDatomA, iDatomB), helpScal(1:lmpT(iDtypeA)) )

                  write(475,*) 'ikGqz1 z0', idirR, idirC
                  write(475,*) dot_product( bigAz1mat(1:lmpT(iDtypeA), idirR, idirC, iDatomA, iDatomB), helpScal(1:lmpT(iDtypeA)) )

                  helpScal(1:lmpT(iDtypeA)) = matmul( varphiHvarphi(1:lmpT(iDtypeA), 1:lmpT(iDtypeA), iDtypeA), bigAz1vec(1:lmpT(iDtypeA), idirR, iDatomA, iDatomB) ) &
                    & - eig(iband, ikpt, 1) * matmul( varphiVarphi(1:lmpT(iDtypeA), 1:lmpT(iDtypeA), iDtypeA), bigAz1vec(1:lmpT(iDtypeA), idirR, iDatomA, iDatomB) )
                  psiHepsPsi(idirR, idirC) = psiHepsPsi(idirR, idirC) + &
                                           & dot_product( bigAvec(1:lmpT(iDtypeA), idirC, iDatomA), helpScal(1:lmpT(iDtypeA)) )

                  write(475,*) 'ikGz0 z1', idirR, idirC
                  write(475,*) dot_product( bigAvec(1:lmpT(iDtypeA), idirC, iDatomA), helpScal(1:lmpT(iDtypeA)) )

                  helpScal(1:lmpT(iDtypeA)) = matmul( varphiHvarphi(1:lmpT(iDtypeA), 1:lmpT(iDtypeA), iDtypeA), bigAz1mat(1:lmpT(iDtypeA), idirR, idirC, iDatomA, iDatomB) ) &
                    & - eig(iband, ikpt, 1) * matmul( varphiVarphi(1:lmpT(iDtypeA), 1:lmpT(iDtypeA), iDtypeA), bigAz1mat(1:lmpT(iDtypeA), idirR, idirC, iDatomA, iDatomB) )
                  psiHepsPsi(idirR, idirC) = psiHepsPsi(idirR, idirC) + &
                                           & dot_product( bigAscal(1:lmpT(iDtypeA), iDatomA), helpScal(1:lmpT(iDtypeA)) )

                  write(475,*) 'z0 ikGqz1', idirR, idirC
                  write(475,*) dot_product( bigAscal(1:lmpT(iDtypeA), iDatomA), helpScal(1:lmpT(iDtypeA)) )


                  if (iDatomA == iDatomB) then
                    helpScal(1:lmpT(iDtypeA)) = matmul( varphiHvarphi(1:lmpT(iDtypeA), 1:lmpT(iDtypeA), iDtypeA), bigAmat(1:lmpT(iDtypeA), idirR, idirC, iDatomA) ) &
                      & - eig(iband, ikpt, 1) * matmul( varphiVarphi(1:lmpT(iDtypeA), 1:lmpT(iDtypeA), iDtypeA), bigAmat(1:lmpT(iDtypeA), idirR, idirC, iDatomA) )

                    psiHepsPsi(idirR, idirC) = psiHepsPsi(idirR, idirC) - &
                                           &  dot_product( bigAscal(1:lmpT(iDtypeA), iDatomA), helpScal(1:lmpT(iDtypeA)) )

                    write(475,*) 'z0 g2z0', idirR, idirC
                    write(475,*) -dot_product( bigAscal(1:lmpT(iDtypeA), iDatomA), helpScal(1:lmpT(iDtypeA)) )

                    helpScal(1:lmpT(iDtypeA)) = matmul( varphiHvarphi(1:lmpT(iDtypeA), 1:lmpT(iDtypeA), iDtypeA), bigAvec(1:lmpT(iDtypeA), idirC, iDatomA) ) &
                      & - eig(iband, ikpt, 1) * matmul( varphiVarphi(1:lmpT(iDtypeA), 1:lmpT(iDtypeA), iDtypeA), bigAvec(1:lmpT(iDtypeA), idirC, iDatomA) )

                    psiHepsPsi(idirR, idirC) = psiHepsPsi(idirR, idirC) + &
                                           &  dot_product( bigAvec(1:lmpT(iDtypeA), idirR, iDatomA), helpScal(1:lmpT(iDtypeA)) )

                    write(475,*) 'ikGz0 ikGz0 1', idirR, idirC
                    write(475,*) dot_product( bigAvec(1:lmpT(iDtypeA), idirR, iDatomA), helpScal(1:lmpT(iDtypeA)) )

                    helpScal(1:lmpT(iDtypeA)) = matmul( varphiHvarphi(1:lmpT(iDtypeA), 1:lmpT(iDtypeA), iDtypeA), bigAvec(1:lmpT(iDtypeA), idirR, iDatomA) ) &
                      & - eig(iband, ikpt, 1) * matmul( varphiVarphi(1:lmpT(iDtypeA), 1:lmpT(iDtypeA), iDtypeA), bigAvec(1:lmpT(iDtypeA), idirR, iDatomA) )

                    psiHepsPsi(idirR, idirC) = psiHepsPsi(idirR, idirC) + &
                                           &  dot_product( bigAvec(1:lmpT(iDtypeA), idirC, iDatomA), helpScal(1:lmpT(iDtypeA)) )

                    write(475,*) 'ikGz0 ikGz0 2', idirR, idirC
                    write(475,*) dot_product( bigAvec(1:lmpT(iDtypeA), idirC, iDatomA), helpScal(1:lmpT(iDtypeA)) )

                    helpScal(1:lmpT(iDtypeA)) = matmul( varphiHvarphi(1:lmpT(iDtypeA), 1:lmpT(iDtypeA), iDtypeA), bigAscal(1:lmpT(iDtypeA), iDatomA) ) &
                      & - eig(iband, ikpt, 1) * matmul( varphiVarphi(1:lmpT(iDtypeA), 1:lmpT(iDtypeA), iDtypeA), bigAscal(1:lmpT(iDtypeA), iDatomA) )

                    psiHepsPsi(idirR, idirC) = psiHepsPsi(idirR, idirC) - &
                                           &  dot_product( bigAmat(1:lmpT(iDtypeA), idirR, idirC, iDatomA), helpScal(1:lmpT(iDtypeA)) )

                    write(475,*) 'g2z0 z0', idirR, idirC
                    write(475,*) -dot_product( bigAmat(1:lmpT(iDtypeA), idirR, idirC, iDatomA), helpScal(1:lmpT(iDtypeA)) )

                  end if

                  dynMatPu(idirR + (iDatomB - 1) * 3, idirC + (iDatomA - 1) * 3) = &
                & dynMatPu(idirR + (iDatomB - 1) * 3, idirC + (iDatomA - 1) * 3) + 2 * results%w_iks(iband, ikpt, 1) * psiHepsPsi(idirR, idirC)
                end do ! idirR
              end do ! idirC
            end do ! iDeqatB
          end do ! iDtypeB
        end do ! iDeqatA
      end do ! iDtypeA
    end do ! iband

  end subroutine AddAlexPulayBraKets2DynMat

  !--------------------------------------------------------------------------------------------------------------------------------------
  !> @author
  !> Christian-Roman Gerhorst, Gregor Michaliczek
  !>
  !> @brief
  !> Sets up interstitial part of matrix element with first variation effective potential within the Sternheimer equation
  !>
  !> @details
  !> Adds Coulomb- & XC-Potential to the interstitial Hamiltonian h.
  !> It is assumed that the input Hamiltonian already incorporates the
  !> kinetic energy or that this contribution is added afterwards. In
  !> the latter case the Hamiltonian is assumed to be initialized.
  !> It is also assumed that vpw contains thetaV(G), where theta is the
  !> step function and V the potential.
  !>
  !> @todo modernize method!!!
  !> @todo this routine is similiar to calcMEPotIR within jpSternhHF maybe unite it
  !--------------------------------------------------------------------------------------------------------------------------------------
  subroutine calcPsi1HepsPsi1IR( kpts, qpts, stars, cell, Gbas, nv, z1Bra, z1Ket, nmat, ikpt, iqpt, ikpq, kpq2kPrVec, vEff0IR, eig, iband,&
      & hepsIR )

      USE m_constants, ONLY : pimach
      use m_types

#ifdef CPP_FFTW
               use :: FFTW3
#endif

      IMPLICIT NONE

! scalar arguments

      ! Type parameter
      type(t_kpts),      intent(in)  :: kpts
      type(t_kpts),      intent(in)  :: qpts
      type(t_stars),     intent(in)  :: stars
      type(t_cell),      intent(in)  :: cell

      ! Array parameter
      integer,           intent(in)  :: Gbas(:, :)
      integer,           intent(in)  :: nv(:, :)
      complex,           intent(in)  :: z1Bra(:) ! basis functions in PW expansion
      complex,           intent(in)  :: z1Ket(:) ! basis functions in PW expansion
      integer,           intent(in)  :: kpq2kPrVec(:, :, :)
      complex,           intent(in)  :: vEff0IR(:,:)
      real,              intent(in)  :: eig(:, :, :)
      complex,           intent(out) :: hepsIR ! interstitial matrix element with H - eps 1 inside for k summed over bands

      ! Scalar parameter
      integer,           intent(in) :: nmat
      integer,           intent(in) :: ikpq
      integer,           intent(in) :: ikpt
      integer,           intent(in) :: iqpt
      integer,           intent(in) :: iband

      ! local scalar variables
      integer                       :: i, j, istar, m, n, iv, il, im, in, gMirr(3)
      real                          :: fftTime, matElemTime
      integer                       :: maxG(3), nfft(3)
      integer                       :: ifftds, ifft1ds, ifft2ds, ifft3ds
      integer                       :: smx1,smx2,smx3
      integer                       :: xis,yis,zis,xil,yil,zil,smallIndex,largeIndex
      integer                       :: xIndex,yIndex,zIndex
      real                          :: testG, time2, time1
      integer                       :: iBas
      complex                       :: h, s
      integer*8 backwardPlan
      integer*8 forwardPlan
      integer*8 backwardPlanKin
      integer*8 forwardPlanKin
      integer*8 backwardPlanOvl
      integer*8 forwardPlanOvl

! local arrays

      integer                       :: iv1d(MAXVAL(nv),1)
      integer, parameter            :: boxSizesMaxIndex = 200
      integer                       :: fftBoxSizes(1:boxSizesMaxIndex) !array with optimal FFT box sizes larger or equal to array index.
      complex, allocatable          :: theta(:)
      complex, allocatable          :: thetaV(:)

      complex, allocatable          :: tempGrid(:)
      complex, allocatable          :: tempGridOvl(:)
      complex, allocatable          :: zFFTBox(:)
      complex, allocatable          :: zFFTBoxKin(:)
      complex, allocatable          :: zFFTBoxOvl(:)
      complex, allocatable          :: kpGz1(:)
      real                          :: kpqGExt(3)
      complex                       :: vThetaZ(SIZE(z1Bra(:)))
      complex                       :: thetaZ(SIZE(z1Bra(:)))
      complex                       :: thetaZ2(SIZE(z1Bra(:)))
      complex  CPP_BLAS_cdotc
      external CPP_BLAS_cdotc

      intrinsic isign,real,cmplx,aimag,conjg

      external dfftw_plan_dft_3d
      external dfftw_execute_dft
      external dfftw_destroy_plan

! initialization of variables

     fftBoxSizes = [  1,  2,  3,  4,  5,  6,  7,  8,  9, 10, 11, 12, 13, 14, 15, 16, 20, 20, 20, 20, 21, 25, 25, 25, 25, 26, 27, 28, &
                   & 30, 30, 32, 32, 36, 36, 36, 36, 40, 40, 40, 40, 42, 42, 44, 44, 45, 48, 48, 48, 50, 50, 54, 54, 54, 54, 56, 56, &
                   & 60, 60, 60, 60, 64, 64, 64, 64, 66, 66, 72, 72, 72, 72, 72, 72, 75, 75, 75, 80, 80, 80, 80, 80, 84, 84, 84, 84, &
                   & 90, 90, 90, 90, 90, 90, 96, 96, 96, 96, 96, 96, 98, 98,100,100,104,104,104,104,105,108,108,108,112,112, 112,112,&
                   &120,120,120,120,120,120,120,120, 128,128,128,128,128,128,128,128,130,130, 132,132,135,135,135,144,144,144,144,144, &
                   &144,144,144,144,150,150,150,150,150,150, 160,160,160,160,160,160,160,160,160,160, 162,162,168,168,168,168,168,168, &
                   & 180,180, 180,180,180,180,180,180,180,180,180,180, 192,192,192,192,192,192,192,192,192,192, 192,192,200,200,200,200, &
                   & 200,200,200,200]

    hepsIR = cmplx(0.0, 0.0)

    ! todo insert q!
    ! Action of kinetic energy
    allocate( kpGz1(nv(1, ikpq)) )
    kpGz1 = cmplx(0., 0.)
    do iBas = 1, nv(1, ikpq)
      kpqGExt(1:3) = matmul(cell%bmat(1:3, 1:3), kpts%bk(1:3, ikpq) + Gbas(1:3, iBas))
      kpGz1(iBas) = norm2(kpqGext(1:3))**2 * z1Ket(iBas)
    end do ! iBas

! This subroutine performs the following steps:
! 1. Generation of thetaV(r) by performing a FFT on w_vEff1IR.
! 2. Generate basis functions on real-space grid by performing a FFT
!    on z. +++

!***********************
!* Initialization step *
!***********************

      fftTime = 0.0
      matElemTime = 0.0

      ! todo ifft1ds wirklich richtig auch für bra?
      ! determine maxG
      maxG = 0
      DO iv = 1, nv(1,ikpq)
         ! ilst has been applied outside
         il = Gbas(1, iv)
         im = Gbas(2, iv)
         in = Gbas(3, iv)
         testG = abs(il)
         IF (testG > maxG(1)) THEN
            maxG(1) = CEILING(testG)
         END IF
         testG = abs(im)
         IF (testG > maxG(2)) THEN
            maxG(2) = CEILING(testG)
         END IF
         testG = abs(in)
         IF (testG > maxG(3)) THEN
            maxG(3) = CEILING(testG)
         END IF
      END DO
      ! produkt darstellen aus Potential und Stufenfunktion, daher 2 * maxG (satz von parseval)
      smx1 = 3*maxG(1) ! bis 2 kmax(=maxG) + Gmax
      smx2 = 3*maxG(2)
      smx3 = 3*maxG(3)

      ifft1ds=(2*smx1+1)
      ifft2ds=(2*smx2+1)
      ifft3ds=(2*smx3+1)

      IF (ifft1ds <= boxSizesMaxIndex) THEN ! kann ich drin lassen
         ifft1ds = fftBoxSizes(ifft1ds)
      END IF
      IF (ifft2ds <= boxSizesMaxIndex) THEN
         ifft2ds = fftBoxSizes(ifft2ds)
      END IF
      IF (ifft3ds <= boxSizesMaxIndex) THEN
         ifft3ds = fftBoxSizes(ifft3ds)
      END IF

      ifftds=ifft1ds*ifft2ds*ifft3ds

!***************************************************************
!* Step 1: Generation of thetaV(r) by performing a FFT on w_vEff1IR: *
!***************************************************************

!  ---> Initialization

      ALLOCATE (tempGrid(0:27*stars%mx1*stars%mx2*stars%mx3-1))
      ALLOCATE (tempGridOvl(0:27*stars%mx1*stars%mx2*stars%mx3-1))

      ALLOCATE (thetaV(0:ifftds-1))
      allocate (theta(0:ifftds -1))
      thetaV = (0.0,0.0)
      theta = (0.0,0.0)
      ! todo talk with gregor ifft1ds and ifft3ds are interchanged
      CALL dfftw_plan_dft_3d(backwardPlan, ifft1ds,ifft2ds,ifft3ds, thetaV, thetaV, FFTW_BACKWARD, FFTW_ESTIMATE)
      !todo can'we use the same plan here as above?
      CALL dfftw_plan_dft_3d(backwardPlanOvl, ifft1ds,ifft2ds,ifft3ds, theta, theta, FFTW_BACKWARD, FFTW_ESTIMATE)
      thetaV = (0.0,0.0)
      theta = (0.0,0.0)
      tempGrid = (0.0,0.0)
      tempGridOvl= (0.0,0.0)

!  ---> put stars onto the large fft-grid "tempGrid"

      DO istar=0, stars%kimax
      ! vEff0IR is already warped and only has to be decorated with the phase for non-symorphic systems
         tempGrid(stars%igfft(istar, 2)) = vEff0IR(stars%igfft(istar, 1), 1) * stars%pgfft(istar) ! changed!
         tempGridOvl(stars%igfft(istar, 2)) = stars%ustep(stars%igfft(istar, 1)) * stars%pgfft(istar) ! changed!
      ENDDO

!  ---> reduce large fft grid "tempGrid" to small fft grid "thetaV" ("realFFTBox")
      DO xIndex = -smx1,smx1 ! aufpassen bzgl. indices
         xis = xIndex
         xil = xIndex
         IF (xIndex < 0) THEN
            xis = xis+ifft1ds
            xil = xil+3*stars%mx1
         END IF
         DO yIndex = -smx2,smx2
            yis = yIndex
            yil = yIndex
            IF (yIndex < 0) THEN
               yis = yis+ifft2ds
               yil = yil+3*stars%mx2
            END IF
            DO zIndex = -smx3,smx3
               zis = zIndex
               zil = zIndex
               IF (zIndex < 0) THEN
                  zis = zis+ifft3ds
                  zil = zil+3*stars%mx3
               END IF
               smallIndex = (xis+ifft1ds*yis+ifft1ds*ifft2ds*zis)
               largeIndex = (xil+3*stars%mx1*yil+9*stars%mx1*stars%mx2*zil)
               thetaV(smallIndex) = tempGrid(largeIndex)
               theta(smallIndex) = tempGridOvl(largeIndex)
            END DO
         END DO
      END DO

      DEALLOCATE (tempGrid)
      DEALLOCATE (tempGridOvl)

!  ---> perform 3D FFT on "thetaV" from reciprocal space to real space

      CALL dfftw_execute_dft(backwardPlan, thetaV, thetaV)
      CALL dfftw_destroy_plan(backwardPlan)

      CALL dfftw_execute_dft(backwardPlanOvl, theta, theta)
      CALL dfftw_destroy_plan(backwardPlanOvl)


!*********************************************************************
!* Step 2: Generate basis functions on real-space grid by performing *
!*         a FFT on z. Then calculate matrix elements H_ij.          *
!*********************************************************************

   !   !todo beware of the minus sign in front of k_i, nv from ket
      !changed begin
      DO iv = 1, nv(1, ikpq)
         il = Gbas(1, iv)! + kpq2kPrVec(1, ikpt, iqpt)
         im = Gbas(2, iv)! + kpq2kPrVec(2, ikpt, iqpt)
         in = Gbas(3, iv)! + kpq2kPrVec(3, ikpt, iqpt)
         IF (il < 0) THEN
            il = il + ifft1ds
         END IF
         IF (im < 0) THEN
            im = im + ifft2ds
         END IF
         IF (in < 0) THEN
            in = in + ifft3ds
         END IF
         iv1d(iv,1) = il+ifft1ds*im+ifft1ds*ifft2ds*in
      END DO
      !changed end


      ALLOCATE (zFFTBox(0:ifftds-1))
      ALLOCATE (zFFTBoxKin(0:ifftds-1))
      ALLOCATE (zFFTBoxOvl(0:ifftds-1))
      zFFTBox = (0.0,0.0)
      zFFTBoxKin = (0.0,0.0)
      zFFTBoxOvl = (0.0,0.0)

      ! todo discuss with gregor order of ifft1ds!
      ! todo discuss with gregor whether we have to make so much plans
      CALL dfftw_plan_dft_3d(backwardPlan, ifft1ds,ifft2ds,ifft3ds, zFFTBox, zFFTBox, FFTW_BACKWARD, FFTW_MEASURE)
      CALL dfftw_plan_dft_3d(forwardPlan, ifft1ds,ifft2ds,ifft3ds, zFFTBox, zFFTBox, FFTW_FORWARD, FFTW_MEASURE)
      CALL dfftw_plan_dft_3d(backwardPlanKin, ifft1ds,ifft2ds,ifft3ds, zFFTBoxKin, zFFTBoxKin, FFTW_BACKWARD, FFTW_MEASURE)
      CALL dfftw_plan_dft_3d(forwardPlanKin, ifft1ds,ifft2ds,ifft3ds, zFFTBoxKin, zFFTBoxKin, FFTW_FORWARD, FFTW_MEASURE)
      CALL dfftw_plan_dft_3d(backwardPlanOvl, ifft1ds,ifft2ds,ifft3ds, zFFTBoxOvl, zFFTBoxOvl, FFTW_BACKWARD, FFTW_MEASURE)
      CALL dfftw_plan_dft_3d(forwardPlanOvl, ifft1ds,ifft2ds,ifft3ds, zFFTBoxOvl, zFFTBoxOvl, FFTW_FORWARD, FFTW_MEASURE)

      CALL cpu_time(time1)
      zFFTBox = (0.0,0.0)
      zFFTBoxKin = (0.0,0.0)
      zFFTBoxOvl = (0.0,0.0)
      DO iv = 1, nv(1, ikpq)
         zFFTBox(iv1d(iv,1)) = z1Ket(iv)
         zFFTBoxOvl(iv1d(iv,1)) = z1Ket(iv)
         zFFTBoxKin(iv1d(iv,1)) = kpGz1(iv)
      END DO

      CALL dfftw_execute_dft(backwardPlan, zFFTBox, zFFTBox)
      CALL dfftw_execute_dft(backwardPlanOvl, zFFTBoxOvl, zFFTBoxOvl)
      CALL dfftw_execute_dft(backwardPlanKin, zFFTBoxKin, zFFTBoxKin)


      DO j = 0, ifftds-1
         zFFTBox(j) = zFFTBox(j) * thetaV(j)
         zFFTBoxOvl(j) = zFFTBoxOvl(j) * theta(j)
         zFFTBoxKin(j) = zFFTBoxKin(j) * theta(j)
      END DO

      CALL dfftw_execute_dft(forwardPlan, zFFTBox, zFFTBox)
      CALL dfftw_execute_dft(forwardPlanOvl, zFFTBoxOvl, zFFTBoxOvl)
      CALL dfftw_execute_dft(forwardPlanKin, zFFTBoxKin, zFFTBoxKin)

      CALL cpu_time(time2)
      fftTime = fftTime + time2 - time1
      CALL cpu_time(time1)

      vThetaZ = cmplx(0.0,0.0)
      thetaZ = cmplx(0.0, 0.0)
      thetaZ2 = cmplx(0.0, 0.0)
      DO iv = 1, nv(1, ikpq)
         vThetaZ(iv) = zFFTBox(iv1d(iv,1))
         vThetaZ(iv) = vThetaZ(iv) / ifftds
         thetaZ(iv) = zFFTBoxOvl(iv1d(iv,1))
         thetaZ(iv) = thetaZ(iv) / ifftds
         thetaZ2(iv) = zFFTBoxKin(iv1d(iv,1))
         thetaZ2(iv) = thetaZ2(iv) / ifftds
      END DO
      !changed begin
      !conjugation of cdotc z1Bra is being conjugated automatically

      ! vkin
      ! todo what is with omtil?
      ! todo symmetric or non-symmetric version?
      s = cmplx(0.0, 0.0)
      h = cmplx(0.0, 0.0)
      s = CPP_BLAS_cdotc(nmat,z1Bra(1), 1,thetaZ,1)
      ! factor 0.5 from kinetic energy
      h = 0.5 * CPP_BLAS_cdotc(nmat,z1Bra(1), 1,thetaZ2,1)
      ! veffIR
      h = h + CPP_BLAS_cdotc(nmat,z1Bra(1),1,vThetaZ,1)
      hepsIR = h - eig(iband, ikpt, 1) * s!nmat von bra!
      !changed end

      CALL dfftw_destroy_plan(forwardPlan)
      CALL dfftw_destroy_plan(backwardPlan)
      CALL dfftw_destroy_plan(forwardPlanOvl)
      CALL dfftw_destroy_plan(backwardPlanOvl)
      CALL dfftw_destroy_plan(forwardPlanKin)
      CALL dfftw_destroy_plan(backwardPlanKin)

  end subroutine calcPsi1HepsPsi1IR

  subroutine WarpIRPot( stars, ngpqdp, idir, gpqdp, pot, pot_warped )

    use m_types
    use m_cfft

    implicit none

    ! Type parameters
    type(t_stars),                                intent(in)  :: stars

    ! Scalar parameters
    integer,                                      intent(in)  :: ngpqdp
    integer,                                      intent(in)  :: idir

    ! Array parameters
    integer,                                      intent(in)  :: gpqdp(:, :)
    complex,                                      intent(in)  :: pot(:, :)
    complex,                                      intent(out) :: pot_warped(:)

    ! Array variables
    integer,                         allocatable              :: igfft(:)
    real, allocatable :: mygfft(:, :)
    real,                            allocatable              :: gfft(:, :)
    integer                                                  :: nfft(3)
    integer                                                   :: gabs(3)


    ! Scalar variables
    integer                                                   :: ifftd
    integer                                                   :: iG
    real                                                      :: scaling
    integer                                                   :: imesh
    integer :: idir1

    REAL                :: pgfftF(0:(2*stars%mx1+1)*(2*stars%mx2+1)*(2*stars%mx3+1)-1)
    INTEGER             :: igfftF(0:(2*stars%mx1+1)*(2*stars%mx2+1)*(2*stars%mx3+1)-1,2)

    nfft = [3 * stars%mx1, 3 * stars%mx2, 3 * stars%mx3]
    ifftd = product(nfft)

    allocate( igfft(ngpqdp), gfft(0: ifftd - 1, 2) )
    allocate( mygfft(0: ifftd - 1, 2) )
    gabs = 0
    igfft = 0
    do iG = 1, ngpqdp
      do idir1 = 1, 3
        if ( gpqdp(idir1, iG) >= 0 ) then
          gabs(idir1) = gpqdp(idir1, iG)
        else
          gabs(idir1) = gpqdp(idir1, iG) + nfft(idir1)
        end if
      end do
      igfft(iG) = gabs(1) + gabs(2) * nfft(1) + gabs(3) * nfft(1) * nfft(2)
    end do

!    open(1000, file='igfftjp', form='formatted')
!    do iG = 1, ngdp
!    !  write (1000, '(i10,2x,i3,2x,i3,2x,i3,2x,i8)') ( gdp(1, iG) + stars%mx1 ) + 100 * ( gdp(2, iG) + stars%mx2 ) + 10000 * ( gdp(3, iG) + stars%mx3 ), &
!    !    & gdp(1, iG), gdp(2, iG), gdp(3, iG), igfft(iG)
!      write (1000, '(i10,2x,i3,2x,i3,2x,i3,2x,i8)') iG, &
!        & gdp(1, iG), gdp(2, iG), gdp(3, iG), igfft(iG)
!    end do
!    close(1000)

    gfft = 0
    pot_warped=0

      do iG = 1, ngpqdp
        gfft(igfft(iG), 1) = real(pot(iG, idir))
        gfft(igfft(iG), 2) = aimag(pot(iG, idir))
      end do
   !   mygfft = 0
   !   mygfft = gfft

!      open(1000, file='foo', form='formatted')
!      do imesh = 0, ubound( gfft, dim=1)
!        write (1000, '(i8,f20.13)') imesh, gfft(imesh, 1)
!        write (1000, '(i8,f20.13)') imesh, gfft(imesh, 2)
!      end do
!      close(1000)


      call cfft(gfft(0, 1), gfft(0, 2), ifftd, nfft(1), nfft(1), 1)
      call cfft(gfft(0, 1), gfft(0, 2), ifftd, nfft(2), nfft(1) * nfft(2), 1)
      call cfft(gfft(0, 1), gfft(0, 2), ifftd, nfft(3), ifftd, 1)

!      open(1000, file='foo1', form='formatted')
!      do imesh = 0, ubound( gfft, dim=1)
!        write (1000, '(i8,f20.13)') imesh, gfft(imesh, 1)
!        write (1000, '(i8,f20.13)') imesh, gfft(imesh, 2)
!      end do
!      close(1000)

   !   call cfft(mygfft(:, 1), mygfft(:, 2), ifftd, nfft(1), nfft(1), 1)
   !   call cfft(mygfft(:, 1), mygfft(:, 2), ifftd, nfft(2), nfft(1) * nfft(2), 1)
   !   call cfft(mygfft(:, 1), mygfft(:, 2), ifftd, nfft(3), ifftd, 1)


!      open(1000, file='foo2', form='formatted')
!      do imesh = 0, ubound( mygfft, dim=1)
!        write (1000, '(i8,f20.13)') imesh, mygfft(imesh, 1)
!        write (1000, '(i8,f20.13)') imesh, mygfft(imesh, 2)
!      end do
!      close(1000)

!    write (29999, *) stars%ufft
      !gfft = 0
      do imesh = 0, ifftd-1
        gfft(imesh, :) = gfft(imesh, :) * stars%ufft(imesh) ! todo is ufft correctlz initialized
       ! gfft(imesh, :) = stars%ufft(imesh) ! todo is ufft correctlz initialized
      end do

      call cfft(gfft(:, 1), gfft(:, 2), ifftd, nfft(1), nfft(1), -1)
      call cfft(gfft(:, 1), gfft(:, 2), ifftd, nfft(2), nfft(1) * nfft(2), -1)
      call cfft(gfft(:, 1), gfft(:, 2), ifftd, nfft(3), ifftd, -1)

!      open(1000, file='foo5', form='formatted')
!      do imesh = 0, ubound( mygfft, dim=1)
!        write (1000, '(i8,f25.15)') imesh, gfft(imesh, 1)
!        write (1000, '(i8,f25.15)') imesh, gfft(imesh, 2)
!      end do
!      close(1000)

!      igfftF = 0
!      pgfftF = 0
!      open(1000, file='mapArray', form='unformatted')
!      read(1000) igfftF, pgfftF
!      close(1000)

      !write (*, *) 'kimax und ngdp', stars%kimax, ngdp

      scaling = 1. / real(ifftd)
      do iG=1, ngpqdp
        pot_warped(iG) = pot_warped(iG) +  cmplx( gfft(igfft(iG), 1), gfft(igfft(iG), 2) ) * scaling
      end do
  end subroutine warpIRPot

  subroutine Calc2ArgIntIR(cell, ngdp, f, w_g, integral)

    use m_types

    implicit none

    ! Type parameter
    type(t_cell), intent(in)  :: cell

    ! Scalar parameter
    integer,      intent(in)  :: ngdp
    complex,      intent(out) :: integral

    ! Array parameter
    complex,      intent(in)  :: f(:)
    complex,      intent(in)  :: w_g(:)


    ! The complex conjugation is done implicetly for the function passed as first argument of dot_product
    integral = cell%omtil * dot_product( f(:ngdp), w_g(:ngdp) )

  end subroutine Calc2ArgIntIR

  ! Calculates a 2 argument voluyme integral in the MT
  subroutine Calc2ArgCmplxIntMT( atoms, itype, f, g, integral )

    use m_types
    !use m_intgr, only : intgr3
    use m_intgr, only : intgr3!LinIntp! TODO: Is this ok?

    implicit none

    ! Type parameters
    type(t_atoms),             intent(in)  :: atoms

    ! Scalar parameters
    integer,                   intent(in)  :: itype
    complex,                   intent(out) :: integral

    ! Array parameters
    complex,                   intent(in)  :: f(:)
    complex,                   intent(in)  :: g(:)

    ! Scalar variables
    integer                                :: imesh
    real                                   :: integralReal
    real                                   :: integralImag

    ! Array variables
    real,          allocatable             :: intgrdR(:)
    real,          allocatable             :: intgrdI(:)


    allocate( intgrdR(atoms%jmtd), intgrdI(atoms%jmtd) )

    intgrdR(:) = 0.
    intgrdI(:) = 0.
    integral = cmplx(0., 0.)

    ! The intgr3LinIntp subroutine only accepts real quantities, so we split the integral into real and imaginary part.
    do imesh = 1, atoms%jri(itype)
      intgrdR(imesh) = real( atoms%rmsh(imesh, itype)**2 * conjg(f(imesh)) * g(imesh) )
      intgrdI(imesh) = aimag( atoms%rmsh(imesh, itype)**2 * conjg(f(imesh)) * g(imesh) )
    end do ! imesh
    !call intgr3(intgrdR, atoms%rmsh(1,itype), atoms%dx(itype), atoms%jri(itype), integralReal)
    !call intgr3(intgrdI, atoms%rmsh(1,itype), atoms%dx(itype), atoms%jri(itype), integralImag)
    call intgr3( intgrdR, atoms%rmsh(1, itype), atoms%dx(itype), atoms%jri(itype), integralReal )! TODO: Is this ok?
    call intgr3( intgrdI, atoms%rmsh(1, itype), atoms%dx(itype), atoms%jri(itype), integralImag )! TODO: Is this ok?

    integral = integral + cmplx( integralReal, integralImag )

  end subroutine Calc2ArgCmplxIntMT

  subroutine CalcSurfIntIRDynMat( atoms, cell, ngdp1, ngdp2, gdp1, gdp2, rho0IRpw, grVext0IR, qpoint, surfInt )

    use m_types
    use m_ylm
    use m_sphbes

    implicit none

    ! Type parameter
    type(t_atoms),        intent(in)  :: atoms
    type(t_cell),         intent(in)  :: cell

    ! Scalar parameter
    integer,              intent(in)  :: ngdp1
    integer,              intent(in)  :: ngdp2

    ! Array parameter
    integer,              intent(in)  :: gdp1(:, :)
    integer,              intent(in)  :: gdp2(:, :)
    complex,              intent(in)  :: rho0IRpw(:)
    complex,              intent(in)  :: grVext0IR(:, :)
    real,                 intent(in)  :: qpoint(:)
    complex,              intent(out) :: surfInt(3, 3)

    ! Scalar variables
    integer                           :: idirC
    integer                           :: idirR
    integer                           :: iatom
    integer                           :: itype
    integer                           :: ieqat
    integer                           :: iG
    integer                           :: it
    integer                           :: iGp
    complex                           :: phaseFac
    complex                           :: tSummedCitY1t
    complex                           :: pref
    complex                           :: surfIntNonMat

    ! Array variables
    real                              :: gSum(3)
    real                              :: gSumCart(3)
    complex                           :: ylm(4)
    real                              :: sbes(0:1)

    surfInt(:, :) = cmplx(0., 0.)
    iatom = 0
    do itype = 1, atoms%ntype
      pref = fpi_const * ImagUnit * atoms%rmt(itype) * atoms%rmt(itype)
      do ieqat = 1, atoms%neq(itype)
        iatom = iatom + 1
        do iG = 1, ngdp1
          do iGp = 1, ngdp2

            gSum(1:3) = gdp1(1:3, iG) + gdp2(1:3, iGp) + qpoint(1:3)
            gSumCart(1:3) = matmul( cell%bmat(1:3, 1:3), gSum(1:3) )

            ylm(:) = cmplx(0., 0.)
            call ylm4( 1, gSumCart, ylm )

            sbes(:) = 0
            call sphbes(1, norm2(gSumCart) * atoms%rmt(itype), sbes)

            phaseFac = exp( ImagUnit * tpi_const * dot_product(gSum(1:3), atoms%taual(1:3, iatom)))

            surfIntNonMat = pref * phaseFac * sbes(1) * rho0IRpw(iG)

            do idirC = 1, 3
              do idirR = 1, 3
                ! Corresponds to the magnetic quantum number m for the orbital quantum number l = 1
                tSummedCitY1t = cmplx(0., 0.)
                do it = -1, 1
                  tSummedCitY1t = tSummedCitY1t + c_im(idirR, it + 2) * ylm(it + 3)
                end do ! it

                surfInt(idirR, idirC) = surfInt(idirR, idirC) - surfIntNonMat * tSummedCitY1t * grVext0IR(iGp, idirC)
              end do ! idirR
            end do ! idirC
          end do ! iGp
        end do ! iG
      end do ! ieqat
    end do ! itype

  end subroutine CalcSurfIntIRDynMat

  subroutine CalcSurfIntMTDynMat(atoms, sym, lathar, clnu_atom, nmem_atom, mlh_atom, rho0MT, grVext0MT, surfInt)

    use m_types, only : t_atoms, t_sym, t_sphhar
    use m_gaunt, only : Gaunt1

    implicit none

    ! Type parameters
    type(t_atoms),                     intent(in)  :: atoms
    type(t_sym),                     intent(in)  :: sym
    type(t_sphhar),                    intent(in)  :: lathar

    ! Array parameters
    complex,                           intent(in)  :: clnu_atom(:, 0:, :)
    integer,                           intent(in)  :: nmem_atom(0:, :)
    integer,                           intent(in)  :: mlh_atom(:, 0:, :)
    real,                              intent(in)  :: rho0MT(:, 0:, :)
    complex,                           intent(in)  :: grVext0MT(:, :, :, :)
    complex,                           intent(out) :: surfInt(3, 3)

    ! Scalar variables
    integer                                        :: iatom
    integer                                        :: itype
    integer                                        :: ieqat
    integer                                        :: ptsym
    integer                                        :: ilh
    integer                                        :: oqn_l
    integer                                        :: imem
    integer                                        :: mqn_m
    integer                                        :: mqn_mp
    integer                                        :: oqn_lpp
    integer                                        :: mqn_mpp
    integer                                        :: lmpp
    real                                           :: gauntFactor
    integer                                        :: idirC
    integer                                        :: idirR

    surfInt(:, :) = 0
    iatom = 0
    do itype = 1, atoms%ntype
      do ieqat = 1, atoms%neq(itype)
        iatom = iatom + 1
        ptsym = sym%ntypsy(iatom)
        do ilh = 0, lathar%nlh(ptsym)
          oqn_l = lathar%llh(ilh, ptsym)
          do imem = 1, nmem_atom(ilh, iatom)
            mqn_m = mlh_atom(imem, ilh, iatom)
            do mqn_mp = -1, 1
              ! lmax + 1 is okay for grVext1
             ! do oqn_lpp = abs(oqn_l - 1), oqn_l + 1
              do oqn_lpp = 0, atoms%lmax(itype)
                do mqn_mpp = -oqn_lpp, oqn_lpp
             !   mqn_mpp = mqn_m + mqn_mp
                lmpp = oqn_lpp * (oqn_lpp + 1) + 1 + mqn_mpp
!                gauntFactor = Gaunt1( oqn_lpp, oqn_l, 1, mqn_mpp, mqn_m, mqn_mp, atoms%lmax(itype) + 1)
                gauntFactor = Gaunt1( oqn_lpp, oqn_l, 1, mqn_mpp, mqn_m, mqn_mp, atoms%lmax(itype))
                do idirC = 1, 3
                  do idirR = 1, 3
                  !todo check whether rho0MT or and c_im should be conjugated or grVext0MT should e conjugatd. This current solution gives the best results.
                    surfInt(idirR, idirC) = surfInt(idirR, idirC) + c_im(idirR, mqn_mp + 2)  * atoms%rmt(itype)**2         &
                      & * rho0MT(atoms%jri(itype), ilh, itype) * clnu_atom(imem, ilh, iatom) *                           &
                      & conjg(grVext0MT(atoms%jri(itype), lmpp, idirC, iatom)) * gauntFactor
                  end do ! idirC
                end do ! idirR
                end do
              end do ! oqn_lpp
            end do ! mqn_mp
          end do ! imem
        end do ! ilh
      end do ! ieqat
    end do ! itype

  end subroutine CalcSurfIntMTDynMat

end module m_jpSetupDynMat
