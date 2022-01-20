module m_jpSetupDynMatSF

    USE m_constants

  implicit none

contains

  subroutine SetupDynMatSF( fmpi, noco, nococonv, oneD, atoms, input, stars, cell, results, Veff0, kpts, qpts, lathar, sym, usdus, ngdp, iqpt, logUnit, &
      & memd_atom, nobd, gdp, mapKpq2K, rbas1, rbas2, nmem_atom, mlh_atom, clnu_atom, kveclo, iloTable, kpq2kPrVec, nv, mapGbas, &
      & gBas, nRadFun, z0, eig, El, rho0IRpw, rho0MT, ngpqdp, gpqdp, rho1IRPW, rho1MT, vXC0IRst, eXCIRst, vXC0MTlh, eXCMTlh, &
      & vExt1IR, vExt1MT, vHar1IR, vHar1MT, grRho0IR, grRho0MT, grVeff0IR, grVeff0MT, vEff0MT, grVCoul0IR_DM_SF, grVCoul0MT_DM_SF, vCoul1IRtempNoVol, vCoul1MTtempNoVol, dynMatSf )

    use m_types
    use m_dfpt_init, only : Derivative, convertStar2G, mt_gradient_old

    implicit none

    type(t_mpi),                  intent(in)  :: fmpi
    type(t_nococonv),                intent(in)  :: nococonv
    type(t_noco),                intent(in)  :: noco
    type(t_oneD),                intent(in)  :: oneD
    type(t_atoms),                 intent(in)  :: atoms
    type(t_input),                 intent(in)  :: input
    type(t_stars),                 intent(in)  :: stars
    type(t_cell),                  intent(in)  :: cell
    type(t_results),               intent(in)  :: results
    type(t_potden),             intent(in)  :: Veff0
    type(t_kpts),                  intent(in)  :: kpts
    type(t_kpts),                  intent(in)  :: qpts
    type(t_sphhar),                intent(in)  :: lathar
    type(t_sym),                   intent(in)  :: sym
    type(t_usdus),                 intent(in)  :: usdus

    ! Scalar parameters
    integer,                       intent(in)  :: ngdp
    integer,                       intent(in)  :: ngpqdp
    integer,                       intent(in)  :: iqpt
    integer,                       intent(in)  :: logUnit
    integer,                       intent(in)  :: memd_atom

    ! Array parameters
    integer,                       intent(in)  :: nRadFun(0:, :)
    integer,                       intent(in)  :: nobd(:, :)
    integer,                       intent(in)  :: gdp(:, :)
    integer,                       intent(in)  :: gpqdp(:, :)
    integer,                       intent(in)  :: mapKpq2K(:, :)
    integer,                       intent(in)  :: kpq2kPrVec(:, :, :)
    integer,                       intent(in)  :: nv(:, :)
    integer,                       intent(in)  :: mapGbas(:, :, :)
    integer,                       intent(in)  :: gBas(:, :)
    COMPLEX,                       intent(in)  :: z0(:, :, :, :)
    real,                          intent(in)  :: eig(:, :, :)
    real,                          intent(in)  :: rbas1(:, :, 0:, :)
    real,                          intent(in)  :: rbas2(:, :, 0:, :)
    integer,                       intent(in)  :: nmem_atom(0:, :)
    integer,                       intent(in)  :: mlh_atom(:,0:,:)
    complex,                       intent(in)  :: clnu_atom(:,0:,:)
    integer,                       intent(in)  :: kveclo(:,:)
    integer,                       intent(in)  :: iloTable(:, 0:, :)
    real,                          intent(in)  :: El(:, 0:, :, :)
    real,                          intent(in)  :: rho0MT(:, 0:, :, :)
    complex,                       intent(in)  :: rho0IRpw(:)
    complex,                       intent(in)  :: rho1IRPW(:, :, :)
    complex,                       intent(in)  :: rho1MT(:, :, :, :, :)
    complex,                       intent(in)  :: vXC0IRst(:, :)
    complex,                       intent(in)  :: eXCIRst(:)
    real,                          intent(in)  :: vXC0MTlh(:, 0:, :, :)
    real,                          intent(in)  :: eXCMTlh(:, 0:, :)
    complex,                       intent(in)  :: vExt1IR(:, :, :)
    complex,                       intent(in)  :: vHar1IR(:, :, :)
    complex,                       intent(in)  :: vHar1MT(:, :, :, :, :)
    complex,                       intent(in)  :: vExt1MT(:, :, :, :, :)
    complex,                       intent(in)  :: vCoul1IRtempNoVol(:, :)
    complex,                       intent(in)  :: vCoul1MTtempNoVol(:, :, :, :)
    complex,                       intent(in)  :: grVCoul0IR_DM_SF(:, :)
    complex,                       intent(in)  :: grVCoul0MT_DM_SF(:, :, :, :)
    complex,                       intent(in)  :: grRho0IR(:, :)
    complex,                       intent(in)  :: grRho0MT(:, :, :, :)
    complex,                       intent(in)  :: grVeff0IR(:, :)
    complex,                       intent(in)  :: grVeff0MT(:, :, :, :)
    real,                          intent(in)  :: vEff0MT(:, 0:, :)
    complex,  allocatable  ,       intent(out) :: dynMatSf(:, :)

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
    integer                                    :: idir
    integer                                    :: iDeqat
    integer                                    :: iG
    integer                                    :: nRadFunMax
    logical                                    :: testGoldstein
    logical                                    :: fullGrVeff0
    logical                                    :: harSw = .true.
    logical                                    :: extSw = .true.
    logical                                    :: xcSw = .true.
    logical                                    :: vextFull = .true.
    integer                                    :: idirR
    integer                                    :: idirC
    logical                                    :: testMode

    ! Array variables
    integer,           allocatable             :: lmpT(:)
    real,              allocatable             :: varphi1(:, :, :)
    complex,           allocatable             :: grVxcIRKern(:)
    real,              allocatable             :: dKernMTGPts(:, :, :)
    complex,           allocatable             :: ylm2(:, :)
    real,              allocatable             :: gaussWghts(:) ! gaussian weights belonging to gausPts
    complex,           allocatable             :: grVxc0MT(:, :, :, :)
    complex,           allocatable             :: grExcMT(:, :, :, :)
    complex,           allocatable             :: ylm1(:, : )
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
    complex,           allocatable             :: r2grVeff0SphVarphiDummy(:, :, :, :, :)
    complex,           allocatable             :: vEff0NsphGrVarphi(:, :, :, :, :)
    real,              allocatable             :: r2Rho0MT(:, :, :)
    complex,           allocatable             :: r2GrRho0MT(:, :, :, :)
    real,              allocatable             :: r2Vxc0MT(:, :, :)
    complex,           allocatable             :: r2GrVxc0MT(:, :, :, :)
    real,              allocatable             :: r2ExcMT(:, :, :)
    complex,           allocatable             :: r2GrExcMT(:, :, :, :)
    complex,           allocatable             :: rho1PWzero(:, :)
    complex,           allocatable             :: rho1MTzero(:, :, :, :)
    complex,           allocatable             :: grRho0MTContainer(:, :, :, :)
    complex,           allocatable             :: vExt1MTContainer(:, :, :, :)
    complex,           allocatable             :: vHar1MTContainer(:, :, :, :)
    complex,           allocatable             :: vCoulExt1IRContainer(:, :)
    complex,           allocatable             :: vCoulExt1MTContainer(:, :, :, :)
    complex,           allocatable             :: rho1MTgoodContainer(:, :, :, :)
    real,              allocatable             :: v0MTcontainer(:, :, :, :)
    complex,           allocatable             :: vXC0IRpw(:)
    complex,           allocatable             :: eXCIRpw(:)
    complex,           allocatable             :: grExcIR(:, :)
    complex,           allocatable             :: grVxc0IR(:, :)
    complex,           allocatable             :: surfInts(:, :, :, :)
    real                                       :: Gext(3)
    complex                                    :: surfInt(3, 3)
    complex                                    :: surfIntTest(3, 3)
    integer :: lmd

    surfIntTest(:, :) = cmplx(.0, .0)
    allocate( dynMatSF( 3 * atoms%nat, 3 * atoms%nat) )
    dynMatSF = cmplx(0., 0.)

    ! Some routines have a flag for setting the z1 = -i(k + G)z0 so that they cancel
    testMode = .false.

    ! Required to initialize
    allocate( lmpT(atoms%ntype) )
    lmpT(:) = 0
    do iDtypeB = 1, atoms%ntype
      lmpT(iDtypeB) = sum( [ ( (2 * oqn_l + 1)* nRadFun(oqn_l, iDtypeB), oqn_l = 0, atoms%lmax(iDtypeB) ) ] )
    end do ! itype
    lmpMax     = maxval( lmpT(:) )
    nRadFunMax = maxval( nRadFun(:, :) )

    lmd=atoms%lmaxd*(atoms%lmaxd+2)
    allocate( varphi1(atoms%jmtd, nRadFunMax, 0:atoms%lmaxd), varphi2(atoms%jmtd, nRadFunMax, 0:atoms%lmaxd) )
    allocate( delrVarphi1( atoms%jmtd, nRadFunMax, 0:atoms%lmaxd), delrVarphi2( atoms%jmtd, nRadFunMax, 0:atoms%lmaxd) )
    allocate( grVarphiCh1(atoms%jmtd, 2, lmpMax, -1:1), grVarphiCh2(atoms%jmtd, 2, lmpMax, -1:1), &
            & grVarphiChLout(2, 0:atoms%lmaxd), grVarphiChMout(-atoms%lmaxd:atoms%lmaxd, -1:1) )
    allocate( vEff0MtSpH( atoms%jmtd, 0:lmd), vEff0NsphGrVarphi(2, atoms%jmtd, 0:(atoms%lmaxd + 3)**2 - 1, lmpMax, -1:1), &
            & r2grVeff0SphVarphi(2, atoms%jmtd, 0:(atoms%lmaxd + 3)**2 - 1, lmpMax, 3) )
    allocate( r2grVeff0SphVarphiDummy(2, atoms%jmtd, 0:(atoms%lmaxd + 3)**2 - 1, lmpMax, 3) )
    allocate( hFullVarphi(2, atoms%jmtd, 0:(atoms%lmaxd + 3)**2 - 1, lmpMax)  )
    allocate( eXCIRpw(ngdp), vXC0IRpw(ngdp) )
    allocate( grExcIR(ngdp, 3), grVxc0IR(ngdp, 3), vCoulExt1IRContainer(ngdp, 3) )
    allocate( grRho0MTContainer(atoms%jmtd, (atoms%lmaxd + 2)**2, 3, atoms%nat), &
            & vExt1MTContainer(atoms%jmtd, (atoms%lmaxd + 2)**2, 3, atoms%nat), &
            & vHar1MTContainer(atoms%jmtd, (atoms%lmaxd + 2)**2, 3, atoms%nat), &
            & rho1MTgoodContainer(atoms%jmtd, (atoms%lmaxd + 2)**2, 3, atoms%nat), &
            & vCoulExt1MTContainer(atoms%jmtd, (atoms%lmaxd + 2)**2, 3, atoms%nat), &
            & v0MTContainer(atoms%jmtd, 0:lathar%nlhd, atoms%ntype, 1) )
    allocate(surfInts(3, 3, atoms%nat, atoms%nat))

    ! Calculating the numerical gradients of Exc0 and Vxc0 read in from FLEUR
    r2grVeff0SphVarphi(:, :, :, :, :) = cmplx(0., 0.)
    r2grVeff0SphVarphiDummy(:, :, :, :, :) = cmplx(0., 0.)
    vEff0NsphGrVarphi(:, :, :, :, :) = cmplx(0., 0.)
    eXCIRpw(:) = cmplx(0., 0.)
    vXC0IRpw(:) = cmplx(0., 0.)
    grRho0MTContainer = cmplx(0., 0.)
    vExt1MTContainer = cmplx(0., 0.)
    vHar1MTContainer = cmplx(0., 0.)
    rho1MTgoodContainer = cmplx(0., 0.)
    vCoulExt1MTContainer = cmplx(0., 0.)
    v0MTContainer = 0.0
    vCoulExt1IRContainer(:, :) = cmplx(0., 0.)

    call convertStar2G( vXC0IRst(:, 1), vXC0IRpw, stars, ngdp, gdp )
    call convertStar2G( eXCIRst(:), eXCIRpw, stars, ngdp, gdp )

    ! Perform analytical gradient of xc-quantities
    grExcIR(:, :) = cmplx(0., 0.)
    grVxc0IR(:, :) = cmplx(0., 0.)
    do idir = 1, 3
      do iG = 1, ngdp
        Gext(1:3) = matmul(cell%bmat(1:3, 1:3), gdp(1:3, iG))
        grExcIR(iG, 1:3) = ImagUnit * Gext(1:3) * excIRpw(iG)
        grVxc0IR(iG, 1:3) = ImagUnit * Gext(1:3) * vxc0IRpw(iG)
      end do ! iG
    end do ! idir

    do idir = 1, 3
      do iG = 1, ngdp
        vCoulExt1IRContainer(iG, idir) = vCoul1IRtempNoVol(iG, idir)  + vExt1IR(iG, idir, 1)
      end do ! iG
    end do ! idir

    ! Perform mt gradient of xc-quantities
    allocate( r2Vxc0MT( atoms%jmtd, 0:lathar%nlhd, atoms%ntype) )
    allocate( grVxc0MT( atoms%jmtd, ( atoms%lmaxd + 2 )**2, 3, atoms%nat  ) )
    allocate( r2ExcMT( atoms%jmtd, 0:lathar%nlhd, atoms%ntype) )
    allocate( grExcMT( atoms%jmtd, ( atoms%lmaxd + 2 )**2, 3, atoms%nat  ) )
    r2Vxc0MT(:, :, :) = 0.
    grVxc0MT(:, :, :, :) = cmplx(0., 0.)
    r2ExcMT(:, :, :) = 0.
    grExcMT(:, :, :, :) = cmplx(0., 0.)

    do iDtypeB = 1, atoms%ntype
      do ilh = 0, lathar%nlhd
        do imesh = 1, atoms%jri(iDtypeB)
          ! SPIN MISSING
          r2Vxc0MT(imesh, ilh, iDtypeB) = vXC0MTlh(imesh, ilh, iDtypeB, 1) * atoms%rmsh(imesh, iDtypeB) * atoms%rmsh(imesh, iDtypeB)
          r2ExcMT(imesh, ilh, iDtypeB) = eXCMTlh(imesh, ilh, iDtypeB) * atoms%rmsh(imesh, iDtypeB) * atoms%rmsh(imesh, iDtypeB)
        end do ! imesh
      end do ! ilh
    end do ! iDtypeB

    call mt_gradient_old( atoms, lathar, sym, clnu_atom, nmem_atom, mlh_atom, r2Vxc0MT, r2GrVxc0MT )
    call mt_gradient_old( atoms, lathar, sym, clnu_atom, nmem_atom, mlh_atom, r2ExcMT, r2GrExcMT )

    !todo we should implement a consistent order of indices
    do idir = 1, 3
      iDatomB = 0
      do iDtypeB = 1, atoms%ntype
        do iDeqat = 1, atoms%neq(iDtypeB)
          iDatomB = iDatomB + 1
          do oqn_l = 0, atoms%lmax(iDtypeB)
            lm_pre = oqn_l * (oqn_l + 1) + 1
            do mqn_m = -oqn_l, oqn_l
              lm = lm_pre + mqn_m
              do imesh = 1, atoms%jri(iDtypeB)
                grVxc0MT(imesh, lm, idir, iDatomB) = r2GrVxc0MT(imesh, lm, iDatomB, idir)! / atoms%rmsh(imesh, iDtypeB)**2
                grExcMT(imesh, lm, idir, iDatomB ) = r2GrExcMT(imesh, lm, iDatomB, idir)! / atoms%rmsh(imesh, iDtypeB)**2
              end do ! imesh
            end do ! mqn_m
          end do ! oqn_l
        end do ! iDeqat
      end do ! iDtypeB
    end do ! idir

    ! todo We need to implement this container because the order of indices between the gradient of the potential and the first
    ! variation of the potential is not consistent
    do idir = 1, 3
      iDatomB = 0
      do iDtypeB = 1, atoms%ntype
        do iDeqat = 1, atoms%neq(iDtypeB)
          iDatomB = iDatomB + 1
          do oqn_l = 0, atoms%lmax(iDtypeB)
            lm_pre = oqn_l * (oqn_l + 1) + 1
            do mqn_m = -oqn_l, oqn_l
              lm = lm_pre + mqn_m
              do imesh = 1, atoms%jri(iDtypeB)
                grRho0MTContainer(imesh, lm, idir, iDatomB) = grRho0MT(imesh, lm, iDatomB, idir)
                ! todo: wrong for polyatomic systems
                vCoulExt1MTContainer(imesh, lm, idir, iDatomB) = vCoul1MTtempNoVol(imesh, lm, iDatomB, idir) + grVCoul0MT_DM_SF(imesh, lm, idir, iDatomB) + vExt1MT(imesh, lm, iDatomB, idir, iDatomB)
              end do ! imesh
            end do ! mqn_m
          end do ! oqn_l
        end do ! iDeqat
      end do ! iDtypeB
    end do ! idir
    !do idir = 1, 3
    !  iDatomA = 0
    !  do iDtypeA = 1, atoms%ntype
    !    do iDeqat = 1, atoms%neq(iDtypeA)
    !      iDatomA = iDatomA + 1
    !      do oqn_l = 0, atoms%lmax(iDtypeA)
    !        lm_pre = oqn_l * (oqn_l + 1) + 1
    !        do mqn_m = -oqn_l, oqn_l
    !          lm = lm_pre + mqn_m
    !          do imesh = 1, atoms%jri(iDtypeA)
    !            rho1MTgoodContainer(imesh, lm, idir, iDatomA) = rho1MT(imesh, lm, iDatomA, idir, iDatomB)
    !          end do ! imesh
    !        end do ! mqn_m
    !      end do ! oqn_l
    !    end do ! iDeqat
    !  end do ! iDtypeB
    !end do ! idir

    !iDatomA = 0
    !do iDtypeA = 1, atoms%ntype
    !  do ilh = 0, lathar%nlhd
    !    do imesh = 1, atoms%jri(iDtypeA)
    !      v0MTContainer(imesh, ilh, iDtypeA, 1) = vEff0MT(imesh, ilh, iDtypeA)
    !    end do ! imesh
    !  end do ! ilh
    !end do ! iDtypeB

    surfInts(:, :, :, :) = cmplx(0., 0.)
    iDatomA = 0
    do iDtypeA = 1, atoms%ntype
      do iDeqatA = 1, atoms%neq(iDtypeA)
        iDatomA = iDatomA + 1

    ! todo We need to implement this container because the order of indices between the gradient of the potential and the first
    ! variation of the potential is not consistent
        do idir = 1, 3
          iDatomB = 0
          do iDtypeB = 1, atoms%ntype
            do iDeqat = 1, atoms%neq(iDtypeB)
              iDatomB = iDatomB + 1
              do oqn_l = 0, atoms%lmax(iDtypeB)
                lm_pre = oqn_l * (oqn_l + 1) + 1
                do mqn_m = -oqn_l, oqn_l
                  lm = lm_pre + mqn_m
                  do imesh = 1, atoms%jri(iDtypeB)
                    vHar1MTContainer(imesh, lm, idir, iDatomB) = vHar1MT(imesh, lm, iDatomB, idir, iDatomA)
                    rho1MTgoodContainer(imesh, lm, idir, iDatomB) = rho1MT(imesh, lm, iDatomB, idir, iDatomA)
                  end do ! imesh
                end do ! mqn_m
              end do ! oqn_l
            end do ! iDeqat
          end do ! iDtypeB
        end do ! idir

        do idir = 1, 3
          iDatomB = 0
          do iDtypeB = 1, atoms%ntype
            do iDeqat = 1, atoms%neq(iDtypeB)
              iDatomB = iDatomB + 1
              do oqn_l = 0, atoms%lmax(iDtypeB)
                lm_pre = oqn_l * (oqn_l + 1) + 1
                do mqn_m = -oqn_l, oqn_l
                  lm = lm_pre + mqn_m
                  do imesh = 1, atoms%jri(iDtypeB)
                    vExt1MTContainer(imesh, lm, idir, iDatomB) = vExt1MT(imesh, lm, iDatomB, idir, iDatomA)
                  end do ! imesh
                end do ! mqn_m
              end do ! oqn_l
            end do ! iDeqat
          end do ! iDtypeB
        end do ! idir

        iDatomB = 0
        do iDtypeB = 1, atoms%ntype
          do iDeqatB = 1, atoms%neq(iDtypeB)
            iDatomB = iDatomB + 1

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
                if (.false.) then
                  call Derivative( varphi1(1:atoms%jri(iDtypeB), iradf, oqn_l), iDtypeB, atoms, delrVarphi1(1:atoms%jri(iDtypeB), iradf, oqn_l) )
                  call Derivative( varphi2(1:atoms%jri(iDtypeB), iradf, oqn_l), iDtypeB, atoms, delrVarphi2(1:atoms%jri(iDtypeB), iradf, oqn_l) )
                end if
              end do ! iradf
            end do ! oqn_l

            ! Calculate the application of the gradient and the gradient's dyadic product onto the MT basis functions (matching coefficients
            ! have no spatial dependence) and determing its scattering channels.
            grVarphiChLout(:, :) = 0
            grVarphiChMout(:, :) = 0
            grVarphiCh1(:, :, :, :) = 0.
            grVarphiCh2(:, :, :, :) = 0.
            if (.false.) then
              call CalcChannelsGrFlpNat( atoms, iDtypeB, nRadFun, varphi1, varphi2, delrVarphi1, delrVarphi2, grVarphiChLout, grVarphiChMout, &
                                                                                                              & grVarphiCh1, grVarphiCh2 )
            end if

            vEff0MtSpH(:, :) = cmplx(0.0, 0.0)
            ptsym = sym%ntypsy(iDatomB)
            do ilh = 0, lathar%nlh(ptsym)
              oqn_l = lathar%llh(ilh, ptsym)
              lm_pre = oqn_l * (oqn_l + 1)
              do imem = 1, nmem_atom(ilh, iDatomB)
                mqn_m = mlh_atom(imem, ilh, iDatomB)
                lm = lm_pre + mqn_m
                !todo one could only evaluate the vEff0MtSpH which have a contribution, i.e. the oqn_l and mqn_m which are in llh and mlh_atom
                ! maybe construct a pointer and run only over them to make it memory efficient.
                do imesh = 1, atoms%jri(iDtypeB)
                  vEff0MtSpH(imesh, lm) = vEff0MtSpH(imesh, lm) + Veff0%mt(imesh, ilh, iDtypeB, 1) * clnu_atom(imem, ilh, iDatomB)
                end do ! imesh
              end do ! imem
            end do ! ilh

            hFullVarphi = cmplx(0.0, 0.0)
            call CalcHnGrV0Varphi( atoms, sym, lathar, iDtypeB, iDatomB, lmpMax, El, varphi1, varphi2, nRadFun, vEff0MtSpH, Veff0%mt(:, :, :, 1), clnu_atom, &
              & nmem_atom, mlh_atom, grVarphiCh1, grVarphiCh2, grVarphiChLout, grVarphiChMout, hFullVarphi, vEff0NsphGrVarphi, r2grVeff0SphVarphi, r2grVeff0SphVarphiDummy )
            deallocate(r2grVeff0SphVarphiDummy)


            ! Note: We need always a pair of matrix elements to hermitize out the phase that hinders us from the z1 cancelling with
            !       the i(k + G) z0. That is Psi1 heps Psi does only cancle with gradPsi Heps Psi up to 9e-4, but the Psi H - eps Psi1
            !       does cancel with Psi H - eps gradPsi up to -9e-4 so we have error cancellation. Therefore, we need also a pair
            !       of matrix elements that cancel each other. The same holds for the MT. Only then, we get the same effect as
            !       in the real density.
            ! IR Psi1 Heps Psi
            surfInt(:, :) = cmplx(0., 0.)
            call CalcSFintIRPsi1HepsPsi( atoms, input, stars, cell, results, Veff0, kpts, qpts, ngdp, iqpt, iDtypeB, iDatomB, iDatomA, nobd,&
                                                  & gdp, mapKpq2K, kpq2kPrVec, nv, mapGbas, gBas, z0, eig, surfInt, testMode )
            if (testMode) then
              surfInts(1:3, 1:3, iDatomB, iDatomA) = surfInts(1:3, 1:3, iDatomB, iDatomA) + surfInt(1:3, 1:3)
            else
              !(5.3.184), 1st integral
              surfInts(1:3, 1:3, iDatomB, iDatomA) = surfInts(1:3, 1:3, iDatomB, iDatomA) + 2 * transpose(surfInt(1:3, 1:3)) / 2
            end if

            do idirC=1,3
              do idirR=1,3
                write(465,*) 'PuSF Psi1*'
                write(465,*) idirR, idirC, surfInt(idirR, idirC)
                write(565,*) 'PuSF Psi1*'
                write(565,*) idirR, idirC, surfInt(idirR, idirC)
              end do
            end do

            ! IR Psi Heps Psi1
            surfInt(:, :) = cmplx(0., 0.)
            call CalcSFintIRPsiHepsPsi1( atoms, input, stars, cell, results, Veff0, kpts, qpts, ngdp, iqpt, iDtypeB, iDatomB, iDatomA,&
                                                  & nobd, gdp, mapKpq2K, kpq2kPrVec, nv, mapGbas, gBas, z0, eig, surfInt, testMode )
            if (testMode) then
              surfInts(1:3, 1:3, iDatomB, iDatomA) = surfInts(1:3, 1:3, iDatomB, iDatomA) + surfInt(1:3, 1:3)
            else
              !(5.3.184), 2nd integral
              surfInts(1:3, 1:3, iDatomB, iDatomA) = surfInts(1:3, 1:3, iDatomB, iDatomA) + 2 * surfInt(1:3, 1:3) / 2
            end if

            do idirC=1,3
              do idirR=1,3
                write(465,*) 'PuSF Psi1'
                write(465,*) idirR, idirC, surfInt(idirR, idirC)
                write(565,*) 'PuSF Psi1'
                write(565,*) idirR, idirC, surfInt(idirR, idirC)
              end do
            end do

            ! MT Psi1 Heps Psi and Psi Heps Psi1 with z1 varied
            !surfInt(:, :) = cmplx(0., 0.)
            !call CalcSFintMTPsi1HepsPsiAndPsiHepsPsi1ExpCoeffVar( atoms, sym, dimens, usdus, kpts, cell, results, iqpt, iDtypeB, &
            !  & iDatomB, iDatomA, lmpMax, nRadFun, eig, varphi1, varphi2, hFullVarphi, mapKpq2K,  gBas, mapGbas, nv,&
            !  & kveclo, z0, lmpT, nobd, iloTable, surfInt, testMode )
!            surfInts(1:3, 1:3, iDatomB, iDatomA) = surfInts(1:3, 1:3, iDatomB, iDatomA) + 2 * surfInt(1:3, 1:3)

            !(5.3.183), 1st integral [2 components]
            ! 1/2
            ! 2 Vext1 rho0 IR
            !surfInt(:, :) = cmplx(0., 0.)
            !call CalcSurfIntIRDynMat( atoms, cell, ngdp, ngpqdp, gdp, gpqdp, rho0IRpw, 2 * vExt1IR(:, :, iDatomA), qpts%bk(:, iqpt), surfInt )
            !surfInts(1:3, 1:3, iDatomB, iDatomA) = surfInts(1:3, 1:3, iDatomB, iDatomA) + surfInt(1:3, 1:3)

            !do idirC=1,3
              !do idirR=1,3
                !write(465,*) '(5.3.183), 1st integral, 2Vext component'
                !write(465,*) idirR, idirC, surfInt(idirR, idirC)
              !end do
            !end do

            ! Vhar1 rho0IR
            ! 2/2
            !surfInt(:, :) = cmplx(0., 0.)
            !call CalcSurfIntIRDynMat( atoms, cell, ngdp, ngpqdp, gdp, gpqdp, rho0IRpw, vHar1IR(:, :, iDatomA), qpts%bk(:, iqpt), surfInt )
            !surfInts(1:3, 1:3, iDatomB, iDatomA) = surfInts(1:3, 1:3, iDatomB, iDatomA) + surfInt(1:3, 1:3)

            !do idirC=1,3
            !  do idirR=1,3
            !    write(465,*) '(5.3.183), 1st integral, VH component'
            !    write(465,*) idirR, idirC, surfInt(idirR, idirC)
            !  end do
            !end do

            !!!latest:
            surfInt(:, :) = cmplx(0., 0.)
            call CalcSurfIntMTDynMat( atoms, sym, lathar, clnu_atom, nmem_atom, mlh_atom, Veff0%mt(:, :, :, 1), rho1MTgoodContainer, surfInt)
            surfInts(1:3, 1:3, iDatomB, iDatomA) = surfInts(1:3, 1:3, iDatomB, iDatomA) - transpose(surfInt(1:3, 1:3))

            do idirC=1,3
              do idirR=1,3
                write(465,*) 'PuSF rho1good V'
                write(465,*) idirR, idirC, -surfInt(idirC, idirR)
                write(565,*) 'PuSF rho1good V'
                write(565,*) idirR, idirC, -surfInt(idirC, idirR)
              end do
            end do


            surfInt(:, :) = cmplx(0., 0.)
            call CalcSurfIntIRDynMat( atoms, cell, ngdp, ngdp, gdp, gdp, rho0IRpw, vCoulExt1IRContainer, qpts%bk(:, iqpt), surfInt )
            surfInts(1:3, 1:3, iDatomB, iDatomA) = surfInts(1:3, 1:3, iDatomB, iDatomA) + surfInt(1:3, 1:3)
            surfIntTest(1:3, 1:3) = surfIntTest(1:3, 1:3) + surfInt(1:3, 1:3)

            surfInt(:, :) = cmplx(0., 0.)
            call CalcSurfIntIRDynMat( atoms, cell, ngdp, ngdp, gdp, gdp, rho0IRpw, grVCoul0IR_DM_SF, qpts%bk(:, 1), surfInt )
            surfInts(1:3, 1:3, iDatomB, iDatomA) = surfInts(1:3, 1:3, iDatomB, iDatomA) + surfInt(1:3, 1:3)
            surfIntTest(1:3, 1:3) = surfIntTest(1:3, 1:3) + surfInt(1:3, 1:3)

            surfInt(:, :) = cmplx(0., 0.)
            call CalcSurfIntMTDynMat( atoms, sym, lathar, clnu_atom, nmem_atom, mlh_atom, rho0MT(:, :, :, 1), vCoulExt1MTContainer, surfInt)
            surfInts(1:3, 1:3, iDatomB, iDatomA) = surfInts(1:3, 1:3, iDatomB, iDatomA) + surfInt(1:3, 1:3)
            surfIntTest(1:3, 1:3) = surfIntTest(1:3, 1:3) + surfInt(1:3, 1:3)
            write(*, *) 'surfintTest'
            write(*, '(3(2(es16.8,1x),3x))') surfIntTest(1, :)
            write(*, '(3(2(es16.8,1x),3x))') surfIntTest(2, :)
            write(*, '(3(2(es16.8,1x),3x))') surfIntTest(3, :)
            !todo check whether rho0MT and the potential is correctly run through
            !(5.3.182), 1st integral [2 components]; provide "good" quantities as input!!!
            ! 1/2
            ! 2 Vext1 rho0MT
            !surfInt(:, :) = cmplx(0., 0.)
            !call CalcSurfIntMTDynMat( atoms, lathar, clnu_atom, nmem_atom, mlh_atom, rho0MT(:, :, :, 1), 2 * vExt1MTContainer, surfInt)
            !surfInts(1:3, 1:3, iDatomB, iDatomA) = surfInts(1:3, 1:3, iDatomB, iDatomA) + surfInt(1:3, 1:3)

            !do idirC=1,3
            !  do idirR=1,3
            !    write(465,*) '(5.3.182), 1st integral, 2Vext component'
            !    write(465,*) idirR, idirC, surfInt(idirR, idirC)
            !  end do
            !end do

            !todo check whether rho0MT and the potential is correctly run through
            ! 2/2
            ! Vhar rho0MT
            !surfInt(:, :) = cmplx(0., 0.)
            !call CalcSurfIntMTDynMat( atoms, lathar, clnu_atom, nmem_atom, mlh_atom, rho0MT(:, :, :, 1), vHar1MTContainer, surfInt)
            !surfInts(1:3, 1:3, iDatomB, iDatomA) = surfInts(1:3, 1:3, iDatomB, iDatomA) + surfInt(1:3, 1:3)

            !do idirC=1,3
            !  do idirR=1,3
            !    write(465,*) '(5.3.182), 1st integral, VH component'
            !    write(465,*) idirR, idirC, surfInt(idirR, idirC)
            !  end do
            !end do

            if (iDatomA == iDatomB) then

              ! Basis set correction part i(k + G) of MT Psi1 Heps Psi and Psi Heps Psi1
              !surfInt(:, :) = cmplx(0., 0.)
              !call CalcSFintMTPsi1HepsPsiAndPsiHepsPsi1BasVarikpG( atoms, sym, dimens, usdus, kpts, cell, results, lmpMax, iDtypeA,      &
              !  & iDatomA, nRadFun, eig, hFullVarphi, gBas, mapGbas, nv, kveclo, z0, nobd, lmpT, iloTable, varphi1, varphi2, surfInt )
!              surfInts(1:3, 1:3, iDatomB, iDatomA) = surfInts(1:3, 1:3, iDatomB, iDatomA) + 2 * surfInt(1:3, 1:3)

              if (.false.) then

                ! gradPsi Heps Psi in MT
                surfInt(:, :) = cmplx(0., 0.)
                call CalcSFintMTgradPsiHepsPsi( fmpi, noco, nococonv, oneD, atoms, input, cell, kpts, sym, results, usdus, iDtypeB, iDatomB, nRadFun, &
                  & lmpMax, nobd, nv, gBas, mapGbas, kveclo, z0, eig, hFullVarphi, iloTable, grVarphiChLout, &
                  & grVarphiChMout, grVarphiCh1, grVarphiCh2, varphi1, varphi2, surfInt )

                ! Basis set correction part grad in Bra
                !surfInts(1:3, 1:3, iDatomB, iDatomA) = surfInts(1:3, 1:3, iDatomB, iDatomA) - 2 * transpose(surfInt(1:3, 1:3))

                ! Surface integral with gradient in Bra
                !surfInts(1:3, 1:3, iDatomB, iDatomA) = surfInts(1:3, 1:3, iDatomB, iDatomA) + 2 * surfInt(1:3, 1:3)

                surfInt(:, :) = cmplx(0., 0.)
                ! Psi Heps gradPsi  MT
                call CalcSFintMTPsiHepsGradPsi( fmpi, noco, nococonv, oneD, atoms, input, kpts, sym, cell, usdus, results, iDtypeB, iDatomB, varphi1, varphi2, &
                  & nv, El, gBas, eig, lmpMax, mapGbas, nRadFun, kveclo, nobd, z0, iloTable, grVarphiChLout, grVarphiChMout, &
                  & r2grVeff0SphVarphi, vEff0NsphGrVarphi, lmpT, grVarphiCh1, grVarphiCh2, surfInt )

                ! Basis set correction part grad in Ket
                !surfInts(1:3, 1:3, iDatomB, iDatomA) = surfInts(1:3, 1:3, iDatomB, iDatomA) - 2 * transpose(surfInt(1:3, 1:3))

                ! Surface integral with gradient in Ket
                !surfInts(1:3, 1:3, iDatomB, iDatomA) = surfInts(1:3, 1:3, iDatomB, iDatomA) + 2 * surfInt(1:3, 1:3)

              end if

              !(5.3.184), 3rd integral [3 components]
              ! 1/3
              ! gradPsi Teps Psi IR
              surfInt(:, :) = cmplx(0., 0.)
              call CalcSFintIRgradPsiHepsPsi( atoms, input, stars, cell, kpts, results, Veff0, ngdp, iDtypeB, iDatomB, nobd, eig, gBas, &
                                                                                                 & mapGbas, nv, gdp, z0, surfInt )
              surfInts(1:3, 1:3, iDatomB, iDatomA) = surfInts(1:3, 1:3, iDatomB, iDatomA) + surfInt(1:3, 1:3)

              do idirC=1,3
                do idirR=1,3
                  write(465,*) 'PuSF grad Psi*'
                  write(465,*) idirR, idirC, surfInt(idirR, idirC)
                  write(565,*) 'PuSF grad Psi*'
                  write(565,*) idirR, idirC, surfInt(idirR, idirC)
                end do
              end do

              ! 2/3
              ! psi Teps grPsi IR
              surfInt(:, :) = cmplx(0., 0.)
              call CalcSFintIRPsiHepsGradPsi( atoms, input, stars, cell, kpts, results, Veff0, ngdp, iDtypeB, iDatomB, nobd, eig, gBas, &
                mapGbas, nv, gdp, z0, surfInt )
              surfInts(1:3, 1:3, iDatomB, iDatomA) = surfInts(1:3, 1:3, iDatomB, iDatomA) + transpose(surfInt(1:3, 1:3))

              do idirC=1,3
                do idirR=1,3
                  write(465,*) 'PuSF grad Psi'
                  write(465,*) idirR, idirC, surfInt(idirR, idirC)
                  write(565,*) 'PuSF grad Psi'
                  write(565,*) idirR, idirC, surfInt(idirR, idirC)
                end do
              end do

              ! 3/3
              ! psi GrVeff Psi IR
              !surfInt(:, :) = cmplx(0., 0.)
              !call CalcSFintIRPsigrVeffPsi( atoms, stars, cell, dimens, kpts, results, ngdp, iDtypeB, iDatomB, nobd, eig, gBas, &
              !  mapGbas, nv, gdp, z0, grVeff0IR, surfInt )
              !surfInts(1:3, 1:3, iDatomB, iDatomA) = surfInts(1:3, 1:3, iDatomB, iDatomA) + surfInt(1:3, 1:3)

              !do idirC=1,3
              !  do idirR=1,3
              !    write(465,*) '(5.3.184), 3rd integral, grad V component'
              !    write(465,*) idirR, idirC, surfInt(idirR, idirC)
              !  end do
              !end do

              ! psi grVeff Psi MT
              !surfInt(:, :) = cmplx(0., 0.)
              !call CalcSurfintMTPsigrVeff0Psi( atoms, input, sym, dimens, cell, kpts, usdus, results, lmpMax, iDtypeB, iDatomB, nobd, gbas, &
              !  & mapGbas, nRadFun, nv, z0, kveclo, varphi1, varphi2, grVeff0MT, iloTable, surfInt )
              !surfInts(1:3, 1:3, iDatomB, iDatomA) = surfInts(1:3, 1:3, iDatomB, iDatomA) + surfInt(1:3, 1:3)

              !(5.3.183), 2nd integral [4 components]
              ! 1/4
              ! rho grVxc0IR IR
              !surfInt(:, :) = cmplx(0., 0.)
              !call CalcSurfIntIRDynMat( atoms, cell, ngdp, ngdp, gdp, gdp, rho0IRpw, grVxc0IR, [0., 0., 0.], surfInt )
              !surfInts(1:3, 1:3, iDatomB, iDatomA) = surfInts(1:3, 1:3, iDatomB, iDatomA) - transpose(surfInt(1:3, 1:3))

              !do idirC=1,3
              !  do idirR=1,3
              !    write(465,*) '(5.3.183), 2nd integral, rho (grad Vxc) component'
              !    write(465,*) idirR, idirC, -surfInt(idirC, idirR)
              !  end do
              !end do

              ! 2/4
              ! vxc grRho IR
              !surfInt(:, :) = cmplx(0., 0.)
              !call CalcSurfIntIRDynMat( atoms, cell, ngdp, ngdp, gdp, gdp, vXC0IRpw, grRho0IR, [0., 0., 0.], surfInt )
              !surfInts(1:3, 1:3, iDatomB, iDatomA) = surfInts(1:3, 1:3, iDatomB, iDatomA) - transpose(surfInt(1:3, 1:3))

              !do idirC=1,3
              !  do idirR=1,3
              !    write(465,*) '(5.3.183), 2nd integral, (grad rho) Vxc component'
              !    write(465,*) idirR, idirC, -surfInt(idirC, idirR)
              !  end do
              !end do

              ! rho0MT grVxc0mt
              !surfInt(:, :) = cmplx(0., 0.)
              !call CalcSurfIntMTDynMat( atoms, lathar, clnu_atom, nmem_atom, mlh_atom, rho0MT(:, :, :, 1), grVxc0MT, surfInt)
              !surfInts(1:3, 1:3, iDatomB, iDatomA) = surfInts(1:3, 1:3, iDatomB, iDatomA) - transpose(surfInt(1:3, 1:3))

              ! vxc0mtlh grRho0MT MT
              !surfInt(:, :) = cmplx(0., 0.)
              !call CalcSurfIntMTDynMat( atoms, lathar, clnu_atom, nmem_atom, mlh_atom, vXC0MTlh(:, :, :, 1), grRho0MTContainer, surfInt)
              !surfInts(1:3, 1:3, iDatomB, iDatomA) = surfInts(1:3, 1:3, iDatomB, iDatomA) - transpose(surfInt(1:3, 1:3))

              ! 3/4
              ! rho0ir grExcIR
              !surfInt(:, :) = cmplx(0., 0.)
              !call CalcSurfIntIRDynMat( atoms, cell, ngdp, ngdp, gdp, gdp, rho0IRpw, grExcIR, [0., 0., 0.], surfInt )
              !surfInts(1:3, 1:3, iDatomB, iDatomA) = surfInts(1:3, 1:3, iDatomB, iDatomA) + transpose(surfInt(1:3, 1:3))

              !do idirC=1,3
              !  do idirR=1,3
              !    write(465,*) '(5.3.183), 2nd integral, rho (grad exc) component'
              !    write(465,*) idirR, idirC, surfInt(idirC, idirR)
              !  end do
              !end do

              ! 4/4
              ! grRho0IR excIRpw IR
              ! Note : The gradient of the density is to be continious because the matching coefficients found for the wavefunctions,
              ! from which the density is contructed from also ensure continuity of the gradient of the wavefunctions, of which the
              ! the gradient of the density is made from. Therefore, such integrals should vanish and they do up to 1e-6.
              !surfInt(:, :) = cmplx(0., 0.)
              !call CalcSurfIntIRDynMat( atoms, cell, ngdp, ngdp, gdp, gdp, eXCIRpw, grRho0IR, [0., 0., 0.], surfInt )
              !surfInts(1:3, 1:3, iDatomB, iDatomA) = surfInts(1:3, 1:3, iDatomB, iDatomA) + transpose(surfInt(1:3, 1:3))

              !do idirC=1,3
              !  do idirR=1,3
              !    write(465,*) '(5.3.183), 2nd integral, (grad rho) exc component'
              !    write(465,*) idirR, idirC, surfInt(idirC, idirR)
              !  end do
              !end do

              ! rho0mt grExcmt mt
              !surfInt(:, :) = cmplx(0., 0.)
              !call CalcSurfIntMTDynMat( atoms, lathar, clnu_atom, nmem_atom, mlh_atom, rho0MT(:, :, :, 1), grExcMT, surfInt)
              !surfInts(1:3, 1:3, iDatomB, iDatomA) = surfInts(1:3, 1:3, iDatomB, iDatomA) + transpose(surfInt(1:3, 1:3))

              !excmt grRho0MT
              !surfInt(:, :) = cmplx(0., 0.)
              !call CalcSurfIntMTDynMat( atoms, lathar, clnu_atom, nmem_atom, mlh_atom, eXCMTlh, grRho0MTContainer, surfInt)
              !surfInts(1:3, 1:3, iDatomB, iDatomA) = surfInts(1:3, 1:3, iDatomB, iDatomA) + transpose(surfInt(1:3, 1:3))
            end if
          end do ! iDeqatA
        end do ! iDtypeA
      end do ! iDeqatB
    end do ! iDtypeB

    iDatomA = 0
    do iDtypeA = 1, atoms%ntype
      do iDeqat = 1, atoms%neq(iDtypeA)
        iDatomA = iDatomA + 1
        iDatomB = 0
        do iDtypeB = 1, atoms%ntype
          do iDeqatB = 1, atoms%neq(iDtypeB)
            iDatomB = iDatomB + 1
            do idirR = 1, 3
              do idirC = 1, 3
                dynMatSf(idirR + (iDatomB - 1) * 3, idirC + (iDatomA - 1) * 3) = surfInts(idirR, idirC, iDatomB, iDatomA)
              end do ! idirC
            end do ! idirR
          end do ! iDeqatB
        end do ! iDtypeB
      end do ! iDeqatA
    end do ! iDtypeA

  end subroutine SetupDynMatSF


  ! Calculates the surface integral in the Sternheimer equation, where from the Hamiltonian only the kinetic energy operator is
  ! calculated here. The part with the effective potential is calculated seperately.
  !--------------------------------------------------------------------------------------------------------------------------------------
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
  !--------------------------------------------------------------------------------------------------------------------------------------
  subroutine calcSintKinEnergOvl( atoms, cell, kpts, results, iDtypeB, iDatomB, ikptBra, ikptKet, idir, nv, nobd, gBasBra, gBasKin, &
      & gBasKet, gBasKin2, zBra, zKet, sIntT, sInt )

    use m_ylm_old
    use m_sphbes
    use m_types

    implicit none

    ! Type parameter
    type(t_atoms),             intent(in)  :: atoms
    type(t_cell),              intent(in)  :: cell
    type(t_kpts),              intent(in)  :: kpts
    type(t_results),           intent(in)  :: results

    ! Scalar parameters
    integer,                   intent(in)  :: iDatomB
    integer,                   intent(in)  :: iDtypeB
    integer,                   intent(in)  :: ikptBra
    integer,                   intent(in)  :: ikptKet
    integer,                   intent(in)  :: idir

    ! Array parameters
    integer,                   intent(in)  :: nv(:)
    integer,                   intent(in)  :: nobd
    real,                      intent(in)  :: gBasBra(:, :)
    real,                      intent(in)  :: gBasKin(:, :)
    real,                      intent(in)  :: gBasKet(:, :)
    real,                      intent(in)  :: gBasKin2(:, :)
    COMPLEX,                  intent(in)  :: zBra(:, :)
    COMPLEX,                  intent(in)  :: zKet(:, :)
    complex,                   intent(out) :: sIntT(:, :)
    complex,                   intent(out) :: sInt(:, :)

    ! Local scalar variable
    ! constPreFacs : constant prefactors
    ! sumM         : sum of m-dependent factors
    ! i* : loop indices
    ! GGqNorm : holds norm of variables of arguments of exponential function being Raileigh expanded
    ! lm: index for l and m
    ! tempNoMdep : temporary complex prefactor which is not dependent on m
    complex                                :: constPreFacs
    complex                                :: sumM
    integer                                :: ibandKet
    integer                                :: ibandBra
    integer                                :: iBasK
    integer                                :: iBasB
    real                                   :: GGqNorm
    integer                                :: lm
    complex                                :: tempNoMdep

    ! basisTemp : Temparray containing the basis vectors
    ! kinEffect : Temparray containing the effect of the kinetic operator
    ! Local array variable
    complex,       allocatable             :: basisTemp(:, :)
    complex,       allocatable             :: basisTempKin(:, :)
    complex,       allocatable             :: basisTempKinzKet(:, :)
    complex,       allocatable             :: basisTempzKet(:, :)
    real,          allocatable             :: kinEffect(:)
    real,          allocatable             :: kinEffect2(:,:)
    real                                   :: GpkCart(3)
    real                                   :: GpqCart(3)
    real                                   :: GGqInt(3)
    real                                   :: GGqCart(3)
    complex                                :: ylm(4)
    real                                   :: sphBesJ(0 : 1)

    ! Reset intent(out) variable
    sIntT(:, :) = cmplx(0.0, 0.0)
    sInt(:, :) = cmplx(0.0, 0.0)

    allocate(basisTemp(nv(ikptBra), nv(ikptKet)))
    allocate(basisTempKin(nv(ikptBra), nv(ikptKet)))
    allocate(basisTempKinzKet(nv(ikptBra), nobd))
    allocate(basisTempzKet(nv(ikptBra), nobd))
    allocate(kinEffect(nv(ikptKet)))
    allocate(kinEffect2(nv(ikptBra),nv(ikptKet)))
    basisTemp(:, :) = cmplx(0.0, 0.0)
    basisTempKin(:, :) = cmplx(0.0, 0.0)
    basisTempKinzKet(:, :) = cmplx(0.0, 0.0)
    basisTempzKet(:, :) = cmplx(0.0, 0.0)
    kinEffect(:) = 0.
    kinEffect2(:,:) = 0.

    constPreFacs = -fpi_const * ImagUnit / cell%omtil * atoms%rmt(iDtypeB)**2

    do iBasK = 1, nv(ikptKet)

      ! Kinetic energy
      gpkCart(1:3) = matmul(cell%bmat(1:3, 1:3), gBasKin(1:3, iBasK) + kpts%bk(1:3, ikptKet))
      kinEffect(iBasK)= (norm2(gpkCart))**2

      ! Overlap of basis functions
      do iBasB = 1, nv(ikptBra)
        gpqCart(1:3) = matmul(cell%bmat(1:3, 1:3), gBasKin2(1:3, iBasB) + kpts%bk(1:3, ikptBra))
        kinEffect2(iBasB,iBasK) = gpqCart(1)*gpkCart(1) + gpqCart(2)*gpkCart(2)  + gpqCart(3)*gpkCart(3)
        GGqInt(1:3) = gBasKet(1:3, iBasK) - gBasBra(1:3, iBasB)
        GGqCart(1:3) = matmul(cell%bmat(1:3, 1:3), GGqInt(1:3))
        GGqNorm = norm2(GGqCart(:))
        ! Using the parity of the spherical harmonics we can put the
        ! minus into constPreFacs. We only have l = 1 because we have used the orthogonality relation of the spherical
        ! harmonic between the Y_lm from the natural coordinates and the Y_lm from the Rayleigh expansion
        call ylm4(1, GGqCart, ylm)
        call sphbes(1, GGqNorm * atoms%rmt(iDtypeB), sphBesJ)
        tempNoMdep =  sphBesJ(1) * exp(ImagUnit * tpi_const * dot_product(GGqInt(1:3), atoms%taual(1:3, iDatomB)))
        sumM = (0.0,0.0)
        ! For l = 1, we only have -1 < m < 1
        do lm = 2, 4
          ! We do not have the conjugated natural coordinate expansion coefficients because only the Y_lm were conjugated
          ! by multiplying them with a factor (-1)^m canceling away with the same factor stemming from ylm not being
          ! conjugated, although it is claimed in the Rayleigh expansion
          sumM = sumM + conjg(c_im(idir, lm - 1)) * conjg(ylm(lm))
        end do ! t
        basisTemp(iBasB, iBasK) = tempNoMdep * sumM
        !basisTempKin(iBasB, iBasK) = basisTemp(iBasB, iBasK) * 0.5 * kinEffect(iBasK)
        !!!Symmetrized Ekin:
        basisTempKin(iBasB, iBasK) = basisTemp(iBasB, iBasK) * 0.5 * kinEffect2(iBasB, iBasK)
      end do ! iBasB
    end do ! iBasK

    ! Then we multiply the band dependent quantities
    basisTempzKet(1:nv(ikptBra), 1:nobd) = &
                                       & matmul(basisTemp(1:nv(ikptBra), 1:nv(ikptKet)), zKet(1:nv(ikptKet), 1:nobd))
    basisTempKinzKet(1:nv(ikptBra), 1:nobd) = &
                                       & matmul(basisTempKin(1:nv(ikptBra), 1:nv(ikptKet)), zKet(1:nv(ikptKet), 1:nobd))
    ! TODO Think about the ikptKet in the line below when not dealing with isolators
    do ibandKet = 1, nobd
      do ibandBra = 1, nobd
        sInt(ibandBra, ibandKet)  = constPreFacs * dot_product(zBra(1:nv(ikptBra), ibandBra), basisTempzKet(1:nv(ikptBra), ibandKet))
        sIntT(ibandBra, ibandKet) = constPreFacs * dot_product(zBra(1:nv(ikptBra), ibandBra), basisTempKinzKet(1:nv(ikptBra), ibandKet))
      end do ! ibandBra
    end do ! ibandKet

  end subroutine calcSintKinEnergOvl


  subroutine calcSurfVeff( atoms, input, kpts, cell, ikptBra, ikptKet, coScale, gbas, zeta, nv, nobd, ilst, zBra, zKet, iDtypeB, iDatomB, surfInt )

    use m_gaunt, only : gaunt1
    use m_types
    use m_ylm_old
    use m_sphbes

    implicit none

    ! Type parameters
    type(t_atoms),             intent(in)  :: atoms
    type(t_input),             intent(in)  :: input
    type(t_kpts),              intent(in)  :: kpts
    type(t_cell),              intent(in)  :: cell

    ! Scalar parameters
    integer,                   intent(in)  :: ikptBra
    integer,                   intent(in)  :: ikptKet
    integer,                   intent(in)  :: coScale
    integer,                   intent(in)  :: iDtypeB
    integer,                   intent(in)  :: iDatomB
    integer,                   intent(in)  :: nobd

    ! Array parameters
    integer,                   intent(in)  :: gbas(:, :)
    complex,                   intent(in)  :: zeta(:, :)
    integer,                   intent(in)  :: nv(:, :)
    integer,                   intent(in)  :: ilst(:, :, :)
    COMPLEX,                  intent(in)  :: zBra(:, :)
    COMPLEX,                  intent(in)  :: zKet(:, :)
    complex,                   intent(out) :: surfInt(:, :)

    ! Scalar variables
    real                                   :: pref
    complex                                :: phaseFac
    integer                                :: iband
    integer                                :: oqn_l
    integer                                :: lm_pre
    complex                                :: factL
    integer                                :: mqn_m
    integer                                :: lm
    integer                                :: iG
    integer                                :: idir
    integer                                :: oqn_l2p
    integer                                :: lm2pPre
    integer                                :: lm2p
    integer                                :: oqn_l1p
    integer                                :: lm1pPre
    integer                                :: mqn_m1p
    integer                                :: lm1p
    complex                                :: gauntCoeff
    integer                                :: mqn_m2p
    integer                                :: lmaxScaled


    ! Array variables
    complex,       allocatable             :: psiKCoeff(:, :)
    complex,       allocatable             :: psiBCoeff(:, :)
    complex,       allocatable             :: basCoeff(:, :)
    complex,       allocatable             :: ylm(:)
    real,          allocatable             :: sbes(:)
    real                                   :: gpk(3)
    real                                   :: gpkCart(3)

    lmaxScaled = coScale * atoms%lmax(iDtypeB)

    allocate( psiKCoeff(nobd, (lmaxScaled + 1)**2) )
    allocate( psiBCoeff(input%neig, (lmaxScaled + 1)**2))
    allocate( basCoeff(MAXVAL(nv), (lmaxScaled + 1)**2) )
    allocate( ylm((lmaxScaled + 1)**2) )
    allocate( sbes(0: lmaxScaled) )


    psiKCoeff(:, :) = cmplx(0., 0.)
    psiBCoeff(:, :) = cmplx(0., 0.)

    surfInt(:, :) = cmplx(0., 0.)


    basCoeff(:, :) = cmplx(0., 0.)
    pref = fpi_const**2 / cell%omtil * atoms%rmt(iDtypeB)**2
    do iG = 1, nv(1, ikptKet)
      gpk(1:3) = gbas(1:3, ilst(iG, ikptKet, 1)) + kpts%bk(1:3, ikptKet)
      phaseFac = exp( ImagUnit * tpi_const * dot_product(gpk(:), atoms%taual(1:3, iDatomB)) )

      gpkCart(1:3) = matmul( cell%bmat(1:3, 1:3), gpk(1:3) )

      ylm(:) = cmplx(0., 0.)
      call ylm4( atoms%lmax(1), gpkCart, ylm )

      sbes(:) = 0.
      call sphbes( atoms%lmax(1), norm2(gpkCart(1:3)) * atoms%rmt(iDtypeB), sbes)
      do oqn_l = 0, atoms%lmax(1)!lmaxScaled
        lm_pre = oqn_l * (oqn_l + 1) + 1
        factL = ImagUnit**oqn_l * sbes(oqn_l) * phaseFac
        do mqn_m = -oqn_l, oqn_l
          lm = lm_pre + mqn_m
          basCoeff(iG, lm) = basCoeff(iG, lm) + factL * conjg(ylm(lm))
        end do ! mqn_m
      end do ! oqn_l
    end do ! iG

    do oqn_l = 0, atoms%lmax(1)!lmaxScaled
      lm_pre = oqn_l * (oqn_l + 1) + 1
      do mqn_m = -oqn_l, oqn_l
        lm = lm_pre + mqn_m
        do iband = 1, nobd
          psiKCoeff(iband, lm) = dot_product(conjg(basCoeff(1:nv(1, ikptKet), lm)), zKet(1:nv(1, ikptKet), iband))
          ! Interstitial wavefunction
!          write(2042, '(3(i8),2(f15.8))') ikpt, iband, lm, fpi_const / sqrt(cell%omtil) * psiKCoeff(iband, lm)
        end do ! iband
      end do ! mqn_m
    end do ! oqn_l

    !todo 2lmax or lmax?
      basCoeff(:, :) = cmplx(0., 0.)
      do iG = 1, nv(1, ikptBra)
        gpk(1:3) = gbas(1:3, ilst(iG, ikptBra, 1)) + kpts%bk(1:3, ikptBra)
        phaseFac = exp( -ImagUnit * tpi_const * dot_product(gpk(:), atoms%taual(1:3, iDatomB)) )

        gpkCart(1:3) = matmul( cell%bmat(1:3, 1:3), gpk(1:3) )

        ylm(:) = cmplx(0., 0.)
        call ylm4( lmaxScaled, gpkCart, ylm )

        sbes(:) = 0.
        call sphbes(lmaxScaled, norm2(gpkCart(1:3)) * atoms%rmt(iDtypeB), sbes)
        do oqn_l = 0, lmaxScaled
          lm_pre = oqn_l * (oqn_l + 1) + 1
          factL = conjg(ImagUnit)**oqn_l * sbes(oqn_l) * phaseFac
          do mqn_m = -oqn_l, oqn_l
            lm = lm_pre + mqn_m
            basCoeff(iG, lm) = basCoeff(iG, lm) + factL * ylm(lm)
          end do ! mqn_m
        end do ! oqn_l
      end do ! iG
    do oqn_l = 0, lmaxScaled
      lm_pre = oqn_l * (oqn_l + 1) + 1
      do mqn_m = -oqn_l, oqn_l
        lm = lm_pre + mqn_m
        do iband = 1, nobd
          psiBCoeff(iband, lm) = dot_product(zBra(1:nv(1, ikptBra), iband), basCoeff(1:nv(1, ikptBra), lm))
          ! Complex conjugated interstitial wavefunction
!          write(2045, '(3(i8),2(f15.8))') ikpt, iband, lm, fpi_const / sqrt(cell%omtil) * psiBCoeff(iband, lm)
        end do ! iband
      end do ! mqn_m
    end do ! oqn_l

    do oqn_l2p = 0, lmaxScaled
      lm2pPre = oqn_l2p * (oqn_l2p + 1) + 1
      do mqn_m2p = -oqn_l2p, oqn_l2p
        lm2p = lm2pPre + mqn_m2p
        do oqn_l = 0, lmaxScaled
          lm_pre = oqn_l * (oqn_l + 1) + 1
          do mqn_m = -oqn_l, oqn_l
            lm = lm_pre + mqn_m
            do oqn_l1p = abs( oqn_l2p - oqn_l ), min(oqn_l2p + oqn_l, atoms%lmax(iDtypeB))
              lm1pPre = oqn_l1p * (oqn_l1p + 1) + 1
              mqn_m1p = mqn_m2p + mqn_m
              if ( abs(mqn_m1p) > oqn_l1p ) cycle
              lm1p = lm1pPre + mqn_m1p
              gauntCoeff = gaunt1(oqn_l1p, oqn_l2p, oqn_l, mqn_m1p, mqn_m2p, mqn_m, lmaxScaled)
              do idir = 1, 3
                do iband = 1, nobd
                  surfInt(iband, idir) = surfInt(iband, idir) + pref * zeta(idir, lm2p) * gauntCoeff * psiBCoeff(iband, lm1p) * psiKCoeff(iband, lm)
                end do ! iband
              end do ! idir
            end do ! oqn_l1p
          end do ! mqn_m
        end do ! oqn_l
      end do ! mqn_m2p
    end do ! oqn_l2p

  end subroutine calcSurfVeff

!  subroutine PrepareMTSurfInt( atoms, dimens, kpts, lathar, V0Fleur, iDtype, iDatom, nRadFun, El, varphiKet1, varphiKet2, varphiBra1, varphiBra2, nmem_atom, mlh_atom,&
!      & clnu_atom, hFullNoAbcofBK, overlapNoAbcofBK, z)
  subroutine PrepareMTSurfIntDM( atoms, iDtype, chanMaxBra, chanMaxKet, lmaxBra, lmpMax, mqn_m2PrBra, mqn_m2PrKet, muKet, muBra, lambdaKet,      &
                   & lambdaBra, nRadFun, varphiKet1, varphiKet2, varphiBra1, varphiBra2, hVarphi, jacobiDet, hFullNoAbcofBK, overlapNoAbcofBK )

    use m_types
    use m_gaunt, only : gaunt1

    implicit none

    ! Type parameters
    type(t_atoms),                  intent(in)  :: atoms

    ! Scalar parameters
    integer,                        intent(in)  :: iDtype
    integer,                        intent(in)  :: chanMaxKet
    integer,                        intent(in)  :: chanMaxBra
    integer,                        intent(in)  :: lmaxBra
    integer,                        intent(in)  :: lmpMax
    integer,                        intent(in)  :: mqn_m2PrKet
    integer,                        intent(in)  :: mqn_m2PrBra
    real,                           intent(in)  :: jacobiDet

    ! Array parameters
    integer,                    intent(in)  :: muKet(-atoms%lmaxd:, -1:)
    integer,                    intent(in)  :: muBra(-atoms%lmaxd:, -1:)
    integer,                    intent(in)  :: lambdaKet(:, 0:)
    integer,                    intent(in)  :: lambdaBra(:, 0:)
    integer,                        intent(in)  :: nRadFun(0:, :)
    real,                       intent(in)  :: varphiKet1(:, :, :, -1:)
    real,                       intent(in)  :: varphiKet2(:, :, :, -1:)
    real,                       intent(in)  :: varphiBra1(:, :, :, -1:)
    real,                       intent(in)  :: varphiBra2(:, :, :, -1:)
    complex,                    intent(in)  :: hVarphi(:, :, 0:, :)
    complex,                        intent(out) :: hFullNoAbcofBK(:, :, :, :)
    complex,                        intent(out) :: overlapNoAbcofBK(:, :, :, :)

    ! Scalar variables
    integer                                     :: oqn_l
    integer                                     :: oqn_l4Pr
    integer                                     :: oqn_l5Pr
    integer                                     :: rMt
    integer                                     :: lmp
    integer                                     :: mqn_m
    integer                                     :: mqn_m4Pr
    integer                                     :: iradf
    integer                                     :: ichanKet
    integer                                     :: mqn_m3Pr
    integer                                     :: mqn_m5Pr
    integer                                     :: lm4Pr
    integer                                     :: lm5Pr
    real                                        :: gauntFactor
    integer                                     :: irel
    integer                                     :: lmp1Pr
    integer                                     :: oqn_l1Pr
    integer                                     :: mqn_m1Pr
    integer                                     :: iradf1Pr
    integer                                     :: ichanBra
    integer                                     :: idir

    ! Array variables
    complex,           allocatable              :: ovlKetuV(:, :, :, :)
    complex,           allocatable              :: hamilKetuV(:, :, :, :)

    ! Index of MT boundary on MT logarithmic mesh
    rMt = atoms%jri(iDtype)
    allocate( ovlKetuV(2, (atoms%lmaxd + 2)**2, lmpMax, 3) )
    allocate( hamilKetuV(2, (atoms%lmaxd + 2)**2, lmpMax, 3) )
    ovlKetuV(:, :, :, :) = cmplx(0., 0.)
    hamilKetuV(:, :, :, :) = cmplx(0., 0.)
    hFullNoAbcofBK(:, :, :, :) = cmplx(0., 0.)
    overlapNoAbcofBK(:, :, :, :) = cmplx(0., 0.)


    ! Sum up the Gaunt coefficients of the spherical Hamiltonian
    do idir = 1, 3
      lmp = 0
      do oqn_l = 0, atoms%lmax(iDtype)
        do mqn_m = -oqn_l, oqn_l
          ! The mqn_m4Pr for the overlap
          do iradf = 1, nRadFun(oqn_l, iDtype)
            lmp = lmp + 1
            ! Has to be here otherwise mqn_m4Pr is set wrongly from the other loop beneath
            mqn_m4Pr = muKet(mqn_m, mqn_m2PrKet)
            do ichanKet = 1, chanMaxKet
              oqn_l4Pr = lambdaKet(ichanKet, oqn_l)
              if ( (abs(mqn_m4Pr) > oqn_l4Pr) .or. (oqn_l4Pr < 0) .or. (oqn_l4PR > lmaxBra) ) cycle
              do oqn_l5Pr = 0, atoms%lmax(iDtype)
                ! do oqn_l5Pr = abs(oqn_l4Pr - 1), min(oqn_l4Pr + 1, lmaxBra), 2
                do mqn_m5Pr = -oqn_l5Pr, oqn_l5Pr
              !    mqn_m5Pr = mqn_m3Pr + mqn_m4Pr
                  do mqn_m3Pr = -1, 1
                  !if ( abs(mqn_m5Pr) > oqn_l5Pr ) cycle
                    lm5Pr = oqn_l5Pr * (oqn_l5Pr + 1) + 1 + mqn_m5Pr
                    !gauntFactor = gaunt1( oqn_l5Pr, 1, oqn_l4Pr, mqn_m5Pr, mqn_m3Pr, mqn_m4Pr, lmaxBra + 1)
                    gauntFactor = gaunt1( oqn_l5Pr, 1, oqn_l4Pr, mqn_m5Pr, mqn_m3Pr, mqn_m4Pr, lmaxBra)
                    ! todo it seems that Gaunt rules do not go well with this channel concept
                    if (abs(gauntFactor) < 1e-8) cycle
                    ovlKetUv(1, lm5Pr, lmp, idir) = ovlKetUv(1, lm5Pr, lmp, idir) + gauntFactor * c_im(idir, mqn_m3Pr + 2) &
                                                                                     & * varphiKet1(rMT, ichanKet, lmp, mqn_m2PrKet)
                    ovlKetUv(2, lm5Pr, lmp, idir) = ovlKetUv(2, lm5Pr, lmp, idir) + gauntFactor * c_im(idir, mqn_m3Pr + 2) &
                                                                                     & * varphiKet2(rMT, ichanKet, lmp, mqn_m2PrKet)
                end do ! mqn_m3Pr
               end do ! mqn_m5Pr
              end do ! oqn_l5Pr
            end do ! ichanKet
            do oqn_l4Pr = 0, atoms%lmax(iDtype)
              do mqn_m4Pr = -oqn_l4Pr, oqn_l4Pr
                lm4Pr = oqn_l4Pr * (oqn_l4Pr + 1) + mqn_m4Pr
                do oqn_l5Pr = 0, atoms%lmax(iDtype)
                  do mqn_m5Pr = -oqn_l5Pr, oqn_l5Pr
                !do oqn_l5Pr = abs(oqn_l4Pr - 1), min(oqn_l4Pr + 1, lmaxBra), 2
                    do mqn_m3Pr = -1, 1
                !    mqn_m5Pr = mqn_m3Pr + mqn_m4Pr
                    !if ( abs(mqn_m5Pr) > oqn_l5Pr ) cycle
                      lm5Pr = oqn_l5Pr * (oqn_l5Pr + 1) + 1 + mqn_m5Pr
                      gauntFactor = gaunt1( oqn_l5Pr, 1, oqn_l4Pr, mqn_m5Pr, mqn_m3Pr, mqn_m4Pr, lmaxBra)
                      if (abs(gauntFactor) < 1e-8) cycle
                      do irel = 1, 2
                        hamilKetUv(irel, lm5Pr, lmp, idir) = hamilKetUv(irel, lm5Pr, lmp, idir) + gauntFactor * &
                                                                         & c_im(idir, mqn_m3Pr + 2) * hVarphi(irel, rMT, lm4Pr, lmp)
                      end do ! irel
                    end do ! mqn_m3Pr
                  end do ! mqn_l5Pr
                end do ! oqn_l5Pr
              end do ! mqn_m4Pr
            end do ! oqn_l4Pr
          end do ! iradf
        end do ! mqn_m
      end do ! oqn_l
    end do ! idir

    do idir = 1, 3
      lmp = 0
      do oqn_l = 0, atoms%lmax(iDtype)
        do mqn_m = -oqn_l, oqn_l
          do iradf = 1, nRadFun(oqn_l, iDtype)
            lmp = lmp + 1
            lmp1Pr = 0
            do oqn_l1Pr = 0, atoms%lmax(iDtype)
              do mqn_m1Pr = -oqn_l1Pr, oqn_l1Pr
                mqn_m5Pr = muBra(mqn_m1Pr, mqn_m2PrBra)
                do iradf1Pr = 1, nRadFun(oqn_l1Pr, iDtype)
                  lmp1Pr = lmp1Pr + 1
                  do ichanBra = 1, chanMaxBra
                    oqn_l5Pr = lambdaBra(ichanBra, oqn_l1Pr)
                    if ( (abs(mqn_m5Pr) > oqn_l5Pr) .or. (oqn_l5Pr < 0)) cycle
                    lm5Pr = oqn_l5Pr * (oqn_l5Pr + 1) + 1 + mqn_m5Pr
                    overlapNoAbcofBK(lmp1Pr, lmp, idir, iDtype) = overlapNoAbcofBK(lmp1Pr, lmp, idir, iDtype) + jacobiDet &
                                    & * ( varphiBra1(rMT, ichanBra, lmp1Pr, mqn_m2PrBra) * ovlKetUv(1, lm5Pr, lmp, idir)           &
                                    &   + varphiBra2(rMt, ichanBra, lmp1Pr, mqn_m2PrBra) * ovlKetUv(2, lm5Pr, lmp, idir) )
                    hFullNoAbcofBK(lmp1Pr, lmp, idir, iDtype) = hFullNoAbcofBK(lmp1Pr, lmp, idir, iDtype) + jacobiDet &
                                    & * ( varphiBra1(rMT, ichanBra, lmp1Pr, mqn_m2PrBra) * hamilKetUv(1, lm5Pr, lmp, idir)         &
                                    &   + varphiBra2(rMt, ichanBra, lmp1Pr, mqn_m2PrBra) * hamilKetUv(2, lm5Pr, lmp, idir) )
                    ! todo there might be an inconsistency in the comparison with the benchmark surface integral at 4000 at non-spherical potential
                    !write(4001, '(6i8,2f15.8)') idir, lmp1Pr, oqn_l1Pr, mqn_m1Pr, iradf1Pr, lmp, hFullNoAbcofBK(lmp1Pr, lmp, idir, iDtype)
                  end do ! ichanBra
                end do ! iradf1Pr
              end do ! mqn_m1Pr
            end do ! oqn_l1Pr
          end do ! iradf
        end do ! mqn_m
      end do ! oqn_l
    end do ! idir

  end subroutine PrepareMTSurfIntDM

  subroutine CalcSIntMT(ikpq, ikpt, iDdir, iDatom, iDtype, lmpMax, nobd, hFullNoAbcofBK, overlapNoAbcofBK, abcofBar, abcofKet, overlap, &
                                                                                                                       & surfIntMT )
    implicit none

    ! Scalar parameters
    integer,                        intent(in)  :: ikpt
    integer,                        intent(in)  :: ikpq
    integer,                        intent(in)  :: iDdir
    integer,                        intent(in)  :: iDatom
    integer,                        intent(in)  :: iDtype
    integer,                        intent(in)  :: lmpMax

    ! Array parameters
    integer,                        intent(in)  :: nobd(:, :)
    complex,                        intent(in)  :: hFullNoAbcofBK(:, :, :, :)
    complex,                        intent(in)  :: overlapNoAbcofBK(:, :, :, :)
    complex,                        intent(in)  :: abcofBar(:, :)
    complex,                        intent(in)  :: abcofKet(:, :)
    complex,                        intent(out) :: overlap(:)
    complex,                        intent(out) :: surfIntMT(:)

    ! Scalar variables
    integer                                     :: iband

    ! Array variables
    complex,           allocatable              :: hFullNoAbcofK(:, :)
    complex,           allocatable              :: overlapNoAbcofK(:, :)
    complex,           allocatable              :: hFull(:)

    allocate( hFullNoAbcofK(lmpMax, maxval(nobd(:, :))) )
    allocate( hFull(maxval(nobd(:, :))) )
    allocate( overlapNoAbcofK( lmpMax, maxval(nobd(:, :))) )

    surfIntMT = cmplx(0., 0.)
    hFullNoAbcofK = cmplx(0., 0.)
    hFull = cmplx(0., 0.)
    overlapNoAbcofK = cmplx(0., 0.)
    overlap = cmplx(0., 0.)

    !todo ikpq redundant here? .CRG

    do iband = 1, nobd(ikpt, 1)
    ! Just from here on iDdir is relevant!
    ! delete surfIntMT for every ikpt
      hFullNoAbcofK(1:lmpMax, iband) = matmul( hFullNoAbcofBK(1:lmpMax, 1:lmpMax, iDdir, iDatom), abcofKet(1:lmpMax, iband) )
      overlapNoAbcofK(1:lmpMax, iband) = matmul( overlapNoAbcofBK(1:lmpMax, 1:lmpMax, iDdir, iDtype), abcofKet(1:lmpMax, iband) )
    end do
    do iband = 1, nobd(ikpt, 1)
      surfIntMT(iband) = dot_product( abcofBar(1:lmpMax, iband), hFullNoAbcofK(1:lmpMax, iband) )
      ! has to be iDatom because of mCofs(iDatom)
      overlap(iband) = dot_product( abcofBar(1:lmpMax, iband), overlapNoAbcofK(1:lmpMax, iband) )
    end do ! iband

    ! Note: If we only switch on the potential, this integral compares well to the interstitial version. The same is true for a
    !       constant potential. Comparing directly the wavefunction from the interstitial and the muffin-tin version, the same comes
    !       out. But if we compare the full action of the spherical Hamiltonian, there is a difference in the values which is due to
    !       the fact, that we compare the second derivative of the LAPW basis function which is not continious any more, necessarily
    !       Therefore, we belive this discrepance to come from this fact.
  end subroutine CalcSIntMT

  subroutine IRcoeffVeffUv( atoms, stars, cell, iDtype, iDatom, ngdp, coScale, gdp, vpw_effPw, veffUvIR )

    use m_types
    use m_gaunt, only : gaunt1
    use m_dfpt_init, only : convertStar2G
    use m_sphbes
    use m_ylm_old

    implicit none

    ! Type parameter
    type(t_atoms),             intent(in)  :: atoms
    type(t_stars),             intent(in)  :: stars
    type(t_cell),              intent(in)  :: cell

    ! Scalar parameter
    integer,                   intent(in)  :: iDtype
    integer,                   intent(in)  :: iDatom
    integer,                   intent(in)  :: ngdp
    integer,                   intent(in)  :: coScale

    ! Array parameter
    integer,                   intent(in)  :: gdp(:, :)
    ! should be unwarped, change this later
    complex,                   intent(in)  :: vpw_effPw(:)
    complex,                   intent(out) :: veffUvIR(:, :)

    ! Scalar variables
    integer                                :: iG
    complex                                :: phaseFac
    complex                                :: factG
    integer                                :: oqn_l2p
    integer                                :: mqn_m2p
    integer                                :: oqn_l3p
    integer                                :: mqn_m4p
    complex                                :: factL
    integer                                :: lm2p
    complex                                :: factLM
    integer                                :: mqn_m3p
    integer                                :: lm3p
    real                                   :: gauntCoeff
    integer                                :: idir
    integer                                :: lmaxScale

    ! Array variables
    real,          allocatable             :: sbes(:)
    complex,       allocatable             :: ylm(:)
    real                                   :: gCart(3)

    lmaxScale = coScale * atoms%lmax(iDtype)
    allocate( sbes(0:lmaxScale) )
    allocate( ylm((lmaxScale + 1)**2) )

    veffUvIR = cmplx(0., 0.)

!    ! Load unwarped potential from Fleur
!    allocate( vpw_effSt( stars%ng3, 1 ) )
!    vpw_effSt = 0
!    call Fopen( 1000, name='vpw_eff', status='old', form='unformatted' )
!    read( 1000 ) vpw_effSt
!    call Fclose( 1000 )

!    allocate( vpw_effPw( ngdp ) )
!    vpw_effPw = 0
!    call convertStar2G( vpw_effSt(:, 1), vpw_effPw, stars, ngdp, gdp )
!    call convertStar2G( vpw_effSt(:), vpw_effPw, stars, ngdp, gdp )
!    vpw_effPw = 0
!    do iG = 1, ngdp
!      if (all(gdp(:, iG) == 0)) then
!        vpw_effPw(iG) = 1
!      end if
!    end do ! iG

    ! Unite the unit vector with the effective potential determing a coefficient dependent on direction and lm
    do iG  = 1, ngdp

!      if (all(gdp(:, iG) == 0)) cycle
      phaseFac = exp( ImagUnit * tpi_const * dot_product(gdp(1:3, iG), atoms%taual(1:3, iDatom)) )
      factG = -fpi_const * vpw_effPw(iG) * phaseFac

      gCart(1:3) = matmul( cell%bmat(1:3, 1:3), gdp(1:3, iG) )
      sbes(:) = 0.
      call sphbes(lmaxScale, norm2(gCart(1:3)) * atoms%rmt(iDtype), sbes)

      !TODO remember to change the norm of the ylm!!!
      ylm(:) = cmplx(0., 0.)
      call ylm4( lmaxScale, gCart, ylm )

      do oqn_l2p = 0, lmaxScale
        factL = factG * ImagUnit**oqn_l2p * sbes(oqn_l2p)
        do mqn_m2p = -oqn_l2p, oqn_l2p
          lm2p = oqn_l2p * (oqn_l2p + 1) + 1 + mqn_m2p
          factLM = factL * conjg(ylm(lm2p))
          !todo discuss cutoffs!
          do oqn_l3p = abs(oqn_l2p - 1), min(oqn_l2p + 1, atoms%lmax(iDtype))
          !do oqn_l3p = abs(oqn_l2p - 1), oqn_l2p + 1
            do mqn_m4p = -1, 1
              mqn_m3p = mqn_m2p + mqn_m4p
              if ( abs(mqn_m3p) > oqn_l3p ) cycle
              lm3p = oqn_l3p * (oqn_l3p + 1) + 1 + mqn_m3p
              !gauntCoeff = gaunt1( oqn_l3p, oqn_l2p, 1, mqn_m3p, mqn_m2p, mqn_m4p, coScale * (atoms%lmaxd + 1) )
              gauntCoeff = gaunt1( oqn_l3p, oqn_l2p, 1, mqn_m3p, mqn_m2p, mqn_m4p, lmaxScale )
              do idir = 1, 3
                veffUvIR(idir, lm3p) = veffUvIR(idir, lm3p) + factLM * c_im(idir, mqn_m4p + 2) * gauntCoeff
              end do ! idir
            end do ! mqn_m4p
          end do ! oqn_l3p
        end do ! mqn_m2p
      end do ! oqn_l2p
    end do ! iG

  end subroutine  IRcoeffVeffUv

  subroutine CalcHGrVarphi( atoms, iDtype, mqn_m2PrC, lmpMax, lmaxBra, grVarphiChMout, nRadFun, grVarPhiCh1, grVarPhiCh2, grVarphiChLout, vEff0NsphGrVarphi, El, lmpT, hGrVarphi )

    use m_types, only : t_atoms

    implicit none

    ! Type parameters
    type(t_atoms),              intent(in)  :: atoms

    ! Scalar parameters
    integer,                    intent(in)  :: iDtype
    integer,                    intent(in)  :: mqn_m2PrC
    integer,                    intent(in)  :: lmpMax
    integer,                    intent(in)  :: lmaxBra

    ! Array parameters
    integer,                    intent(in)  :: grVarphiChMout(-atoms%lmaxd:, -1:)
    integer,                    intent(in)  :: nRadFun(0:, :)
    real,                       intent(in)  :: grVarphiCh1(:, :, :, -1:)
    real,                       intent(in)  :: grVarphiCh2(:, :, :, -1:)
    integer,                    intent(in)  :: grVarphiChLout(:, 0:)
    complex,                    intent(in)  :: vEff0NsphGrVarphi(:, :, :, :, -1:)
    real,                       intent(in)  :: El(:, 0:, :, :)
    integer,                    intent(in)  :: lmpT(:)
    complex,                    intent(out) :: hGrVarphi(:, :, :, :, -1:)

    ! Scalar variables
    integer                                 :: mqn_m2Pr
    integer                                 :: lmp
    integer                                 :: oqn_l
    integer                                 :: mqn_m
    integer                                 :: mqn_m3Pr
    integer                                 :: iradf
    integer                                 :: ichan
    integer                                 :: oqn_l3Pr
    integer                                 :: imesh
    integer                                 :: lm3Pr
    integer                                 :: lmMax

    ! Array variables
    real,          allocatable              :: hsphGrVarphiPart(:, :, :, :, :)

    allocate( hsphGrVarphiPart(2, atoms%jmtd, (atoms%lmaxd + 2)**2, lmpMax, -1:1))

    hsphGrVarphiPart = 0.
    hGrVarphi(:, :, :, :, :) = cmplx(0., 0.)

    do mqn_m2Pr = -1, 1
      lmp = 0
      do oqn_l = 0, atoms%lmax(iDtype)
        do mqn_m = -oqn_l, oqn_l
          mqn_m3Pr = grVarphiChMout(mqn_m, mqn_m2Pr)
          do iradf = 1, nRadFun(oqn_l, iDtype)
            lmp = lmp + 1
            ! todo Here, we have to review when introducing LOs
            if (iradf < 3) then
              do ichan = 1, 2
                oqn_l3Pr = grVarphiChLout(ichan, oqn_l)
                if ( ( abs(mqn_m3Pr) > oqn_l3Pr )  .or. ( oqn_l3Pr < 0 ) .or. ( oqn_l3Pr > lmaxBra) )  cycle
                lm3Pr = oqn_l3Pr * (oqn_l3Pr + 1) + 1 + mqn_m3Pr
                do imesh = 1, atoms%jmtd
                   hsphGrVarphiPart(1, imesh, lm3Pr, lmp, mqn_m2Pr) = hsphGrVarphiPart(1, imesh, lm3Pr, lmp, mqn_m2Pr) &
                                                                                       & + El(1, oqn_l, iDtype, 1) * grVarPhiCh1(imesh, ichan, lmp, mqn_m2Pr)
                   hsphGrVarphiPart(2, imesh, lm3Pr, lmp, mqn_m2Pr) = hsphGrVarphiPart(2, imesh, lm3Pr, lmp, mqn_m2Pr) &
                                                                                       & + El(1, oqn_l, iDtype, 1) * grVarPhiCh2(imesh, ichan, lmp, mqn_m2Pr)
                end do ! imesh
                if (iradf == 2) then
                  do imesh = 1, atoms%jri(iDtype)
                    hsphGrVarphiPart(1, imesh, lm3Pr, lmp, mqn_m2Pr) = hsphGrVarphiPart(1, imesh, lm3Pr, lmp, mqn_m2Pr)            &
                                                                                        & + grVarPhiCh1(imesh, ichan, lmp, mqn_m2Pr)
                    hsphGrVarphiPart(2, imesh, lm3Pr, lmp, mqn_m2Pr) = hsphGrVarphiPart(2, imesh, lm3Pr, lmp, mqn_m2Pr)            &
                                                                                        & + grVarPhiCh2(imesh, ichan, lmp, mqn_m2Pr)
                  end do ! imesh
                end if ! iradf == 2
              end do ! ichan
            end if ! iradf < 3
          end do ! iradf
        end do ! mqn_m
      end do ! oqn_l
    end do ! mqn_m2Pr

    !lmMax = (atoms%lmax(iDtype) + 2)**2
    lmMax = (atoms%lmax(iDtype) + 1)**2

    hGrVarphi(1:2, 1:atoms%jri(iDtype), 1:lmMax, 1:lmpT(iDtype), -1:1) = &
     &   hsphGrVarphiPart(1:2, 1:atoms%jri(iDtype), 1:lmMax, 1:lmpT(iDtype), -1:1) &
     & + vEff0NsphGrVarphi(1:2, 1:atoms%jri(iDtype), 1:lmMax, 1:lmpT(iDtype), -1:1)
   ! hGrVarphi(1:2, 1:atoms%jri(iDtype), 1:lmMax, 1:lmpT(iDtype), -1:1) = &
   !  & + vEff0NsphGrVarphi(1:2, 1:atoms%jri(iDtype), 1:lmMax, 1:lmpT(iDtype), -1:1)

  end subroutine CalcHGrVarphi

  subroutine CalcSFintIRgradPsiHepsPsi( atoms, input, stars, cell, kpts, results, Veff0, ngdp, iDtypeB, iDatomB, nobd, eig, gBas, &
                                                                                                   & mapGbas, nv, gdp, z0, surfInt )
    use m_types
    use m_dfpt_init, only : ConvertStar2G

    implicit none

    ! Type parameter
    type(t_atoms),             intent(in)  :: atoms
    type(t_input),             intent(in)  :: input
    type(t_stars),             intent(in)  :: stars
    type(t_cell),              intent(in)  :: cell
    type(t_kpts),              intent(in)  :: kpts
    type(t_results),           intent(in)  :: results
    type(t_potden),         intent(in)  :: Veff0

    ! Scalar parameter
    integer,                   intent(in)  :: ngdp
    integer,                   intent(in)  :: iDtypeB
    integer,                   intent(in)  :: iDatomB

    ! Array parameter
    integer,                   intent(in)  :: nobd(:, :)
    real,                      intent(in)  :: eig(:,:,:)
    integer,                   intent(in)  :: gBas(:, :)
    integer,                   intent(in)  :: mapGbas(:, :, :)
    integer,                   intent(in)  :: nv(:, :)
    integer,                   intent(in)  :: gdp(:, :)
    COMPLEX,                  intent(in)  :: z0(:,:,:,:)
    complex,                   intent(out) :: surfInt(:, :)

    ! Scalar variable
    integer                                :: ikpt
    integer                                :: iBas
    integer                                :: iband
    integer                                :: idirC
    integer                                :: idirR
    integer                                :: maxNobd
    integer                                :: coScale

    ! Array variables
    complex,           allocatable         :: surfIntOvl(:, :)
    complex,           allocatable         :: surfIntVeff0(:, :, :)
    complex,           allocatable         :: surfIntKinEnerg(:, :)
    complex,           allocatable         :: vpw_eff_uw(:)
    complex,           allocatable         :: veffUvIR(:, :)
    real,              allocatable         :: gBasMapped(:, :)
    complex,           allocatable         :: z0gradAction(:, :, :)
    real                                   :: gExt(3)
    real                                   :: kExt(3)

    maxNobd = maxval( nobd(:, :) )
    coScale = 1

    allocate( surfIntOvl( maxNobd, maxNobd ), surfIntVeff0( maxNobd, 3, 3 ), surfIntKinEnerg( maxNobd, maxNobd ) )
    allocate( veffUvIR(3, (coScale * atoms%lmaxd + 1)**2) )
    allocate( vpw_eff_uw( ngdp ) )
    allocate( gBasMapped(3, MAXVAL(nv)) )
    allocate( z0gradAction(MAXVAL(nv), maxNobd, 3) )

    surfInt(:, :) = cmplx(0., 0.)
    surfIntOvl(:, :) = cmplx(0., 0.)
    surfIntVeff0(:, :, :) = cmplx(0., 0.)
    surfIntKinEnerg(:, :) = cmplx(0., 0.)
    veffUvIR(:, :) = cmplx(0., 0.)
    gBasMapped(:, :) = 0.


    vpw_eff_uw = cmplx(0., 0.)
    call convertStar2G( Veff0%pw(:, 1), vpw_eff_uw, stars, ngdp, gdp )
    call IRcoeffVeffUv( atoms, stars, cell, iDtypeB, iDatomB, ngdp, coScale, gdp, vpw_eff_uw, veffUvIR )

    do ikpt = 1, kpts%nkpt

      gBasMapped(:, :) = 0.
      z0gradAction(:, :, :) = cmplx(0., 0.)
      kExt(1:3) = matmul(cell%bmat, kpts%bk(1:3, ikpt))
      gExt(:) = 0.
      ! todo interchange the order of loops
      do iBas = 1, nv(1, ikpt) !+ atoms%nlotot
        gExt(1:3) = matmul( cell%bmat(1:3, 1:3), gBas(1:3, mapGbas(iBas, ikpt, 1)))
        do iband = 1, nobd(ikpt, 1)
          do idirR = 1, 3
            z0gradAction(iBas, iband, idirR) = ImagUnit * ( kExt(idirR) + gExt(idirR) ) * z0(iBas, iband, ikpt, 1)
          end do ! idirR
        end do ! iband
        gBasMapped(1:3, iBas) = real( gBas(1:3, mapGbas(iBas, ikpt, 1)) )
      end do ! iBas

      surfIntVeff0(:, :, :) = cmplx(0., 0.)
      do idirR = 1, 3
        call calcSurfVeff( atoms, input, kpts, cell, ikpt, ikpt, coScale, gBas, veffUvIR, nv, nobd(ikpt, 1), &
          & mapGbas, z0gradAction(:, :, idirR), z0(:, :, ikpt, 1), iDtypeB, iDatomB, surfIntVeff0(:, :, idirR) )
      end do !

      do idirC  = 1, 3
        do idirR = 1, 3

          surfIntKinEnerg(:, :) = cmplx(0., 0.)
          surfIntOvl(:, :) = cmplx(0., 0.)
          call calcSintKinEnergOvl( atoms, cell, kpts, results, iDtypeB, iDatomB, ikpt, ikpt, idirC, &
            & nv(1, :), nobd(ikpt, 1), gBasMapped, gBasMapped, gBasMapped, gBasMapped, z0gradAction(:, :, idirR), z0(:, :, ikpt, 1), surfIntKinEnerg, &
            & surfIntOvl )

          do iband = 1, nobd(ikpt, 1)
            !!!latest:
            surfInt(idirR, idirC) = surfInt(idirR, idirC) + 2 * results%w_iks(iband, ikpt, 1) * ( surfIntKinEnerg(iband, iband) &
                                                   & + 0.0*surfIntVeff0(iband, idirC, idirR) - eig(iband, ikpt, 1) * surfIntOvl(iband, iband) )
!            surfInt(idirR, idirC) = surfInt(idirR, idirC) + results%w_iks(iband, ikpt, 1) * ( surfIntVeff0(iband, idirC, idirR) - eig(iband, ikpt, 1) * surfIntOvl(iband) )
          end do ! iband

        end do ! idirR
      end do ! idirC

    end do ! ikpt

  end subroutine CalcSFintIRgradPsiHepsPsi

  subroutine CalcSFintIRPsiHepsGradPsi( atoms, input, stars, cell, kpts, results, Veff0, ngdp, iDtypeB, iDatomB, nobd, eig, gBas, &
      mapGbas, nv, gdp, z0, surfInt )

    use m_types
    use m_dfpt_init, only : ConvertStar2G

    implicit none

    ! Type parameter
    type(t_atoms),             intent(in)  :: atoms
    type(t_input),             intent(in)  :: input
    type(t_stars),             intent(in)  :: stars
    type(t_cell),              intent(in)  :: cell
    type(t_kpts),              intent(in)  :: kpts
    type(t_results),           intent(in)  :: results
    type(t_potden),         intent(in)  :: Veff0

    ! Scalar parameter
    integer,                   intent(in)  :: ngdp
    integer,                   intent(in)  :: iDtypeB
    integer,                   intent(in)  :: iDatomB

    ! Array parameter
    integer,                   intent(in)  :: nobd(:, :)
    real,                      intent(in)  :: eig(:,:,:)
    integer,                   intent(in)  :: gBas(:, :)
    integer,                   intent(in)  :: mapGbas(:, :, :)
    integer,                   intent(in)  :: nv(:, :)
    integer,                   intent(in)  :: gdp(:, :)
    COMPLEX,                  intent(in)  :: z0(:,:,:,:)
    complex,                   intent(out) :: surfInt(:, :)

    ! Scalar variable
    integer                                :: ikpt
    integer                                :: iBas
    integer                                :: iband
    integer                                :: idirC
    integer                                :: idirR
    integer                                :: maxNobd
    integer                                :: coScale

    ! Array variables
    complex,           allocatable         :: surfIntOvl(:, :)
    complex,           allocatable         :: surfIntVeff0(:, :, :)
    complex,           allocatable         :: surfIntKinEnerg(:, :)
    complex,           allocatable         :: vpw_eff_uw(:)
    complex,           allocatable         :: veffUvIR(:, :)
    real,              allocatable         :: gBasMapped(:, :)
    complex,           allocatable         :: z0gradAction(:, :, :)
    real                                   :: gExt(3)
    real                                   :: kExt(3)

    maxNobd = maxval( nobd(:, :) )
    coScale = 1

    allocate( surfIntOvl( maxNobd, maxNobd ), surfIntVeff0( maxNobd, 3, 3 ), surfIntKinEnerg( maxNobd, maxNobd ) )
    allocate( veffUvIR(3, (coScale * atoms%lmaxd + 1)**2) )
    allocate( vpw_eff_uw( ngdp ) )
    allocate( gBasMapped(3, MAXVAL(nv)) )
    allocate( z0gradAction(MAXVAL(nv), maxNobd, 3) )

    surfInt(:, :) = cmplx(0., 0.)
    surfIntOvl(:, :) = cmplx(0., 0.)
    surfIntVeff0(:, :, :) = cmplx(0., 0.)
    surfIntKinEnerg(:, :) = cmplx(0., 0.)
    veffUvIR(:, :) = cmplx(0., 0.)
    gBasMapped(:, :) = 0.


    vpw_eff_uw = cmplx(0., 0.)
    call convertStar2G( Veff0%pw(:, 1), vpw_eff_uw, stars, ngdp, gdp )
    call IRcoeffVeffUv( atoms, stars, cell, iDtypeB, iDatomB, ngdp, coScale, gdp, vpw_eff_uw, veffUvIR )

    do ikpt = 1, kpts%nkpt

      gBasMapped(:, :) = 0.
      z0gradAction(:, :, :) = cmplx(0., 0.)
      kExt(1:3) = matmul(cell%bmat, kpts%bk(1:3, ikpt))
      gExt(:) = 0.
      ! todo interchange the order of loops
      do iBas = 1, nv(1, ikpt) !+ atoms%nlotot
        gExt(1:3) = matmul( cell%bmat(1:3, 1:3), gBas(1:3, mapGbas(iBas, ikpt, 1)))
        do iband = 1, nobd(ikpt, 1)
          do idirR = 1, 3
            z0gradAction(iBas, iband, idirR) = ImagUnit * ( kExt(idirR) + gExt(idirR) ) * z0(iBas, iband, ikpt, 1)
          end do ! idirR
        end do ! iband
        gBasMapped(1:3, iBas) = real( gBas(1:3, mapGbas(iBas, ikpt, 1)) )
      end do ! iBas

      surfIntVeff0(:, :, :) = cmplx(0., 0.)
      do idirR = 1, 3
        call calcSurfVeff( atoms, input, kpts, cell, ikpt, ikpt, coScale, gBas, veffUvIR, nv, nobd(ikpt, 1), &
          & mapGbas, z0(:, :, ikpt, 1), z0gradAction(:, :, idirR), iDtypeB, iDatomB, surfIntVeff0(:, :, idirR) )
      end do !

      do idirC  = 1, 3
        do idirR = 1, 3

          surfIntKinEnerg(:, :) = cmplx(0., 0.)
          surfIntOvl(:, :) = cmplx(0., 0.)
          call calcSintKinEnergOvl( atoms, cell, kpts, results, iDtypeB, iDatomB, ikpt, ikpt, idirC, &
            & nv(1, :), nobd(ikpt, 1), gBasMapped, gBasMapped, gBasMapped, gBasMapped, z0(:, :, ikpt, 1), z0gradAction(:, :, idirR), surfIntKinEnerg, &
            & surfIntOvl )

          do iband = 1, nobd(ikpt, 1)
            !!!latest:
            surfInt(idirR, idirC) = surfInt(idirR, idirC) + 2 * results%w_iks(iband, ikpt, 1) * ( surfIntKinEnerg(iband, iband) &
                                                   & + 0.0*surfIntVeff0(iband, idirC, idirR) - eig(iband, ikpt, 1) * surfIntOvl(iband, iband) )
!           surfInt(idirR, idirC) = surfInt(idirR, idirC) + results%w_iks(iband, ikpt, 1) * (surfIntVeff0(iband, idirC, idirR) &
!                                                                                       & - eig(iband, ikpt, 1) * surfIntOvl(iband) )
          end do ! iband

        end do ! idirR
      end do ! idirC

    end do ! ikpt

  end subroutine CalcSFintIRPsiHepsGradPsi

  subroutine CalcSFintIRPsigrVeffPsi( atoms, input, stars, cell, kpts, results, ngdp, iDtypeB, iDatomB, nobd, eig, gBas, &
      mapGbas, nv, gdp, z0, grVeff0IR, surfInt )

    use m_types
    use m_dfpt_init, only : ConvertStar2G

    implicit none

    ! Type parameter
    type(t_atoms),             intent(in)  :: atoms
    type(t_input),             intent(in)  :: input
    type(t_stars),             intent(in)  :: stars
    type(t_cell),              intent(in)  :: cell
    type(t_kpts),              intent(in)  :: kpts
    type(t_results),           intent(in)  :: results

    ! Scalar parameter
    integer,                   intent(in)  :: ngdp
    integer,                   intent(in)  :: iDtypeB
    integer,                   intent(in)  :: iDatomB

    ! Array parameter
    integer,                   intent(in)  :: nobd(:, :)
    real,                      intent(in)  :: eig(:,:,:)
    integer,                   intent(in)  :: gBas(:, :)
    integer,                   intent(in)  :: mapGbas(:, :, :)
    integer,                   intent(in)  :: nv(:, :)
    integer,                   intent(in)  :: gdp(:, :)
    COMPLEX,                  intent(in)  :: z0(:,:,:,:)
    complex,                   intent(in)  :: grVeff0IR(:, :)
    complex,                   intent(out) :: surfInt(:, :)

    ! Scalar variable
    integer                                :: ikpt
    integer                                :: iBas
    integer                                :: iband
    integer                                :: idirC
    integer                                :: idirR
    integer                                :: maxNobd
    integer                                :: coScale

    ! Array variables
    complex,           allocatable         :: surfIntOvl(:, :)
    complex,           allocatable         :: surfIntVeff0(:, :, :)
    complex,           allocatable         :: surfIntKinEnerg(:, :)
    complex,           allocatable         :: vpw_eff_uw(:)
    complex,           allocatable         :: veffUvIR(:, :, :)
    real,              allocatable         :: gBasMapped(:, :)
    complex,           allocatable         :: z0gradAction(:, :, :)
    real                                   :: gExt(3)
    real                                   :: kExt(3)

    maxNobd = maxval( nobd(:, :) )
    coScale = 1

    allocate( surfIntOvl( maxNobd, maxNobd ), surfIntVeff0( maxNobd, 3, 3 ), surfIntKinEnerg( maxNobd, maxNobd ) )
    allocate( veffUvIR(3, (coScale * atoms%lmaxd + 1)**2, 3) )
    allocate( gBasMapped(3, MAXVAL(nv)) )

    surfInt(:, :) = cmplx(0., 0.)
    surfIntOvl(:, :) = cmplx(0., 0.)
    surfIntVeff0(:, :, :) = cmplx(0., 0.)
    surfIntKinEnerg(:, :) = cmplx(0., 0.)
    veffUvIR(:, :, :) = cmplx(0., 0.)
    gBasMapped(:, :) = 0.


    do idirR = 1, 3
      call IRcoeffVeffUv( atoms, stars, cell, iDtypeB, iDatomB, ngdp, coScale, gdp, grVeff0IR(:, idirR), veffUvIR(:, :, idirR) )
    end do

    do ikpt = 1, kpts%nkpt

      gBasMapped(:, :) = 0.
      ! todo interchange the order of loops
      do iBas = 1, nv(1, ikpt) !+ atoms%nlotot
        gBasMapped(1:3, iBas) = real( gBas(1:3, mapGbas(iBas, ikpt, 1)) )
      end do ! iBas

      surfIntVeff0(:, :, :) = cmplx(0., 0.)
      do idirR = 1, 3
        call calcSurfVeff( atoms, input, kpts, cell, ikpt, ikpt, coScale, gBas, veffUvIR(:, :, idirR), nv, nobd(ikpt, 1), &
          & mapGbas, z0(:, :, ikpt, 1), z0(:, :, ikpt, 1), iDtypeB, iDatomB, surfIntVeff0(:, :, idirR) )
      end do !

      do idirC  = 1, 3
        do idirR = 1, 3

          surfIntKinEnerg(:, :) = cmplx(0., 0.)
          surfIntOvl(:, :) = cmplx(0., 0.)
          call calcSintKinEnergOvl( atoms, cell, kpts, results, iDtypeB, iDatomB, ikpt, ikpt, idirC, &
            & nv(1, :), nobd(ikpt, 1), gBasMapped, gBasMapped, gBasMapped, gBasMapped, z0(:, :, ikpt, 1), z0(:, :, ikpt, 1), surfIntKinEnerg, &
            & surfIntOvl )

          do iband = 1, nobd(ikpt, 1)
            surfInt(idirR, idirC) = surfInt(idirR, idirC) + results%w_iks(iband, ikpt, 1) * ( surfIntKinEnerg(iband, iband) &
                                                   & + surfIntVeff0(iband, idirC, idirR) - eig(iband, ikpt, 1) * surfIntOvl(iband, iband) )
          end do ! iband

        end do ! idirR
      end do ! idirC

    end do ! ikpt

  end subroutine CalcSFintIRPsigrVeffPsi

  subroutine CalcSFintMTgradPsiHepsPsi( fmpi, noco, nococonv, oneD, atoms, input, cell, kpts, sym, results, usdus, iDtypeB, iDatomB, nRadFun, &
      & lmpMax, nobd, nv, gBas, mapGbas, kveclo, z0, eig, hFullVarphi, iloTable, grVarphiChLout, &
      & grVarphiChMout, grVarphiCh1, grVarphiCh2, varphi1, varphi2, surfInt )

    use m_types
    use m_abcof3
    use m_constants

    implicit none

    ! Type parameters
    type(t_mpi),                  intent(in)  :: fmpi
    type(t_nococonv),                intent(in)  :: nococonv
    type(t_noco),                intent(in)  :: noco
    type(t_oneD),                intent(in)  :: oneD
    type(t_atoms),                  intent(in) :: atoms
    type(t_input),                  intent(in) :: input
    type(t_kpts),                   intent(in) :: kpts
    type(t_cell),                   intent(in) :: cell
    type(t_sym),                    intent(in) :: sym
    type(t_results),                intent(in) :: results
    type(t_usdus),                  intent(in) :: usdus

    ! Scalar parameters
    integer,                        intent(in)  :: iDtypeB
    integer,                        intent(in)  :: iDatomB
    integer,                        intent(in)  :: lmpMax

    ! Array parameters
    integer,                        intent(in)  :: nRadFun(0:, :)
    integer,                        intent(in)  :: nv(:, :)
    integer,                        intent(in)  :: nobd(:, :)
    integer,                        intent(in)  :: gBas(:, :)
    integer,                        intent(in)  :: mapGbas(:, :, :)
    integer,                        intent(in)  :: kveclo(:,:)
    COMPLEX,                       intent(in)  :: z0(:,:,:,:)
    real,                           intent(in)  :: eig(:,:,:)
    complex,                        intent(in)  :: hFullVarphi(:, :, 0:, :)
    integer,                        intent(in)  :: iloTable(:, 0:, :)
    integer,                        intent(in)  :: grVarphiChLout(:, :)
    integer,                        intent(in)  :: grVarphiChMout(:, :)
    real,                           intent(in)  :: grVarphiCh1(:, :, :, :)
    real,                           intent(in)  :: grVarphiCh2(:, :, :, :)
    real,                           intent(in)  :: varphi1(:, :, 0:)
    real,                           intent(in)  :: varphi2(:, :, 0:)
    complex,                        intent(out) :: surfInt(:, :)

    ! Type variables
    type(od_inp)                                :: odi
    type(od_sym)                                :: ods
    type(t_lapw) :: lapw

    ! Scalar variables
    integer                                     :: oqn_l
    integer                                     :: iradf
    integer                                     :: imesh
    integer                                     :: nRadFunMax
    integer                                     :: lm
    integer                                     :: mqn_m
    integer                                     :: chanMaxBra
    integer                                     :: chanMaxKet
    integer                                     :: lmaxBra
    integer                                     :: mqn_m2PrKet
    integer                                     :: lmp
    integer                                     :: mqn_m2PrR
    integer                                     :: ikpt
    integer                                     :: nmat
    integer                                     :: pMaxLocal
    integer                                     :: iband
    integer                                     :: maxNobd
    integer                                     :: idirC
    integer :: nk

    ! Array variables
    integer,           allocatable              :: muKet(:, :)
    real,              allocatable              :: varphiKet1(:, :, :, :)
    real,              allocatable              :: varphiKet2(:, :, :, :)
    integer,           allocatable              :: lambdaKet(:, :)
    complex,           allocatable              :: hFullNoAbcofBK(:, :, :, :, :)
    complex,           allocatable              :: overlapNoAbcofBK(:, :, :, :, :)
    complex,           allocatable              :: a(:, :, :)
    complex,           allocatable              :: b(:, :, :)
    complex,           allocatable              :: bascof_lo(:, :, :, :, :)
    complex,           allocatable              :: surfIntOvl(:)
    complex,           allocatable              :: surfIntHfl(:)
    integer,           allocatable              :: ngoprI(:)
    complex,           allocatable              :: ab0cof(:, :)
    complex                                     :: surfIntNat(-1:1, 3)

    ! Quantities for initialization
    nRadFunMax = maxval( nRadFun(:, iDtypeB) )
    maxNobd    = maxval( nobd(:, :) )

    allocate( muKet(-atoms%lmaxd:atoms%lmaxd, -1:1) )
    allocate( lambdaKet(2, 0:atoms%lmaxd) )
    allocate( varphiKet1(atoms%jmtd, 2, lmpMax, -1:1), varphiKet2(atoms%jmtd, 2, lmpMax, -1:1) )
    allocate( hFullNoAbcofBK(lmpMax, lmpMax, 3, atoms%nat, -1:1), overlapNoAbcofBK(lmpMax, lmpMax, 3, atoms%ntype, -1:1) )
    allocate( ngoprI(atoms%nat) )
    allocate( surfIntOvl( maxNobd ), surfIntHfl( maxNobd ) )
    allocate( ab0cof(lmpMax, maxNobd) )
    allocate( a( MAXVAL(nv), 0:atoms%lmaxd*(atoms%lmaxd+2), atoms%nat), b(MAXVAL(nv), 0:atoms%lmaxd*(atoms%lmaxd+2), atoms%nat), &
      & bascof_lo(3, -atoms%llod:atoms%llod, 4 * atoms%llod + 2, atoms%nlod, atoms%nat) )

    surfIntNat(:, :) = cmplx(0., 0.)
    surfInt(:, :)    = cmplx(0., 0.)
    ngoprI(:) = 1


    chanMaxBra = 2
    chanMaxKet = 1
    lmaxBra = atoms%lmax(iDtypeB)! + 1
    muKet(:, :) = 0
    mqn_m2PrKet = -1
    varphiKet1 = 0.
    varphiKet2 = 0.
    do mqn_m = -atoms%lmax(iDtypeB), atoms%lmax(iDtypeB)
      muKet(mqn_m, -1) = mqn_m
    end do ! mqn_m
    do oqn_l = 0, atoms%lmax(iDtypeB)
      lambdaKet(1, oqn_l) = oqn_l
    end do ! oqn_l
    lmp = 0
    do oqn_l = 0, atoms%lmax(iDtypeB)
      do mqn_m = -oqn_l, oqn_l
        do iradf = 1, nRadFun(oqn_l, iDtypeB)
          lmp = lmp + 1
          varphiKet1(atoms%jri(iDtypeB), 1, lmp, -1) = varphi1(atoms%jri(iDtypeB), iradf, oqn_l )
          varphiKet2(atoms%jri(iDtypeB), 1, lmp, -1) = varphi2(atoms%jri(iDtypeB), iradf, oqn_l )
        end do ! iradf
      end do ! mqn_m
    end do ! oqn_l

    ! todo reorganize direction loops
    hFullNoAbcofBK(:, :, :, :, :) = cmplx(0., 0.)
    overlapNoAbcofBK(:, :, :, :, :) = cmplx(0., 0.)
    do mqn_m2PrR = -1, 1
      call PrepareMTSurfIntDM( atoms, iDtypeB, chanMaxBra, chanMaxKet, lmaxBra, lmpMax, mqn_m2PrR, mqn_m2PrKet, muKet, grVarphiChMout,     &
        & lambdaKet, grVarphiChLout, nRadFun, varphiKet1, varphiKet2, grVarphiCh1, grVarphiCh2, hFullVarphi, atoms%rmt(iDtypeB)**2, hFullNoAbcofBK(:, :, :, :, mqn_m2PrR), overlapNoAbcofBK(:, :, :, :, mqn_m2PrR) )
    end do

      surfIntOvl(:) = cmplx(0., 0.)
      surfIntHfl(:) = cmplx(0., 0.)
      do ikpt = 1, kpts%nkpt
        nmat = nv(1, ikpt) + atoms%nlotot
        a(:, :, :) = cmplx(0.0, 0.0)
        b(:, :, :) = cmplx(0.0, 0.0)
        bascof_lo(:, :, :, :, :) = cmplx(0.0, 0.0)
        !call abcof3( atoms%lmaxd, atoms%ntype, atoms%nat, sym%nop, MAXVAL(nv), input%jspins, 1, atoms%lmaxd*(atoms%lmaxd+2), SIZE(z0(:,1,1,1)), &
        !  & atoms%llod, atoms%nlod, atoms%nlotot, sym%invtab, atoms%ntype, sym%mrot, ngoprI, atoms%taual, atoms%neq, atoms%lmax, &
        !  & atoms%rmt, cell%omtil, cell%bmat, cell%bbmat, kpts%bk(:, ikpt), gBas(1, mapGbas(:nv(1, ikpt), ikpt, 1)), &
        !  & gBas(2, mapGbas(:nv(1, ikpt), ikpt, 1)), gBas(3, mapGbas(:nv(1, ikpt), ikpt, 1)), nv(:, ikpt), nmat, &
        !  & usdus%us, usdus%dus, usdus%uds, usdus%duds, usdus%ddn, sym%invsat, sym%invsatnr, usdus%ulos, usdus%uulon, usdus%dulon, &
        !  & usdus%dulos, atoms%llo, atoms%nlo, atoms%l_dulo, atoms%lapw_l, kveclo(:,ikpt), odi, ods, a, b, bascof_lo )
        nk=fmpi%k_list(ikpt)
        CALL lapw%init(input, noco, nococonv, kpts, atoms, sym, nk, cell, .FALSE., fmpi)
        CALL abcof3(input, atoms, sym, 1, cell, kpts%bk(:, ikpt), lapw, &
                            usdus, oneD, 1, lapw%dim_nvd(), a, b, bascof_lo)
        ab0cof(:, :) = cmplx(0., 0.)
        do iband = 1, nobd(ikpt, 1)
          lmp = 0
          lm  = 0
          do oqn_l = 0, atoms%lmax(iDtypeB)
            do mqn_m = - oqn_l, oqn_l
              pMaxLocal = nRadFun(oqn_l, iDtypeB)
              ! p = 1
              ab0cof(lmp + 1, iband) = ImagUnit**oqn_l * dot_product( conjg(z0(:nv(1, ikpt), iband, ikpt, 1)), a(:nv(1, ikpt), lm, iDatomB) )
              ! p = 2
              ab0cof(lmp + 2, iband) = ImagUnit**oqn_l * dot_product( conjg(z0(:nv(1, ikpt), iband, ikpt, 1)), b(:nv(1, ikpt), lm, iDatomB) )
              ! Add LO contributions
              do iradf = 3, pMaxLocal
                ! p = 1
                ab0cof(lmp + 1, iband) = ab0cof(lmp + 1, iband) + ImagUnit**oqn_l * &
                  & dot_product( conjg(z0(nv(1, ikpt) + 1:nv(1, ikpt) + atoms%nlotot, iband, ikpt, 1)), &
                                                       & bascof_lo(1, mqn_m, :atoms%nlotot, iloTable(iradf, oqn_l, iDtypeB), iDatomB) )
                ! p = 2
                ab0cof(lmp + 2, iband) = ab0cof(lmp + 2, iband) + ImagUnit**oqn_l * &
                  & dot_product( conjg(z0(nv(1, ikpt) + 1:nv(1, ikpt) + atoms%nlotot, iband, ikpt, 1)), &
                                                       & bascof_lo(2, mqn_m, :atoms%nlotot, iloTable(iradf, oqn_l, iDtypeB), iDatomB) )
                ! 2 < p < LOs for that l and that atom type
                ab0cof(lmp + iradf, iband) = ImagUnit**oqn_l * dot_product( conjg(z0(nv(1, ikpt) + 1:nv(1, ikpt) + atoms%nlotot, iband, ikpt, 1)), &
                                                       & bascof_lo(3, mqn_m, :atoms%nlotot, iloTable(iradf, oqn_l, iDtypeB), iDatomB) )
              end do ! iradf


              ! This is a precalculation of the 1st and 3rd line of A.50. Actually the resulting lmp should have been primed or at
              ! least the p as lm1Pr = lm. But for sake of performance we place it here.
              ! sake of performance
              lm = lm + 1
              lmp = lmp + pMaxLocal
            end do ! mqn_m
          end do !oqn_l
        end do ! iband

        do idirC = 1, 3
          do mqn_m2PrR = -1, 1
            surfIntOvl(:) = cmplx(0., 0.)
            surfIntHfl(:) = cmplx(0., 0.)
            call CalcSintMT( ikpt, ikpt, idirC, iDatomB, iDtypeB, lmpMax, nobd, hFullNoAbcofBK(:, :, :, :, mqn_m2PrR), overlapNoAbcofBK(:, :, :, :, mqn_m2PrR), ab0cof, ab0cof,   &
                                                                                                            & surfIntOvl, surfIntHfl )
            do iband = 1, nobd(ikpt, 1)
              surfIntNat(mqn_m2PrR, idirC) = surfIntNat(mqn_m2PrR, idirC) + results%w_iks(iband, ikpt, 1) * ( &
                                                                       & surfIntHfl(iband) - eig(iband, ikpt, 1) * surfIntOvl(iband) )
            end do ! iband
          end do ! mqn_m2PrR
        end do ! idirC
      end do ! ikpt

    surfInt(1:3, 1:3) = matmul( conjg(Tmatrix0(1:3, 1:3)), surfIntNat(-1:1, 1:3) )

  end subroutine CalcSFintMTgradPsiHepsPsi

  subroutine CalcSFintMTPsiHepsGradPsi( fmpi, noco, nococonv, oneD, atoms, input, kpts, sym, cell, usdus, results, iDtypeB, iDatomB, varphi1, varphi2, nv, El, gBas, &
      & eig, lmpMax, mapGbas, nRadFun, kveclo, nobd, z0, iloTable, grVarphiChLout, grVarphiChMout, r2grVeff0SphVarphi, vEff0NsphGrVarphi, lmpT, grVarphiCh1, grVarphiCh2, surfInt )

    use m_abcof3
    use m_types
    use m_types_oneD, only : od_inp, od_sym
    USE m_constants

    implicit none

    ! Type parameter
    type(t_mpi),                  intent(in)  :: fmpi
    type(t_nococonv),                intent(in)  :: nococonv
    type(t_noco),                intent(in)  :: noco
    type(t_oneD),                intent(in)  :: oneD
    type(t_atoms),                 intent(in)  :: atoms
    type(t_input),                 intent(in)  :: input
    type(t_kpts),                  intent(in)  :: kpts
    type(t_sym),                   intent(in)  :: sym
    type(t_cell),                  intent(in)  :: cell
    type(t_usdus),                 intent(in)  :: usdus
    type(t_results),               intent(in)  :: results

    ! Scalar parameter
    integer,                       intent(in)  :: iDtypeB
    integer,                       intent(in)  :: iDatomB
    integer,                       intent(in)  :: lmpMax

    ! Array parameter
    real,                          intent(in)  :: varphi1(:, :, 0:)
    real,                          intent(in)  :: varphi2(:, :, 0:)
    integer,                       intent(in)  :: nv(:, :)
    real,                          intent(in)  :: El(:, 0:, :, :)
    COMPLEX,                      intent(in)  :: z0(:, :, :, :)
    integer,                       intent(in)  :: gbas(:, :)
    integer,                       intent(in)  :: mapGbas(:, :, :)
    integer,                       intent(in)  :: nRadFun(0:, :)
    integer,                       intent(in)  :: kveclo(:,:)
    integer,                       intent(in)  :: nobd(:, :)
    real,                          intent(in)  :: eig(:, :, :)
    integer,                       intent(in)  :: iloTable(:, 0:, :)
    integer,                       intent(in)  :: grVarphiChLout(:, 0:)
    integer,                       intent(in)  :: grVarphiChMout(-atoms%lmaxd:, -1:)
    complex,                       intent(in)  :: r2grVeff0SphVarphi(:, :, 0:, :, :)
    complex,                       intent(in)  :: vEff0NsphGrVarphi(:, :, 0:, :, -1:)
    integer,                       intent(in)  :: lmpT(:)
    real,                          intent(in)  :: grVarphiCh1(:, :, :, :)
    real,                          intent(in)  :: grVarphiCh2(:, :, :, :)
    complex,                       intent(out) :: surfInt(:, :)

    ! Type variables
    type(od_inp)                               :: odi
    type(od_sym)                               :: ods
    type(t_lapw) :: lapw

    ! Scalar variables
    integer                                    :: chanMaxBra
    integer                                    :: chanMaxKet
    integer                                    :: lmaxBra
    integer                                    :: mqn_m2PrBra
    integer                                    :: mqn_m2PrKet
    integer                                    :: mqn_m
    integer                                    :: oqn_l
    integer                                    :: lmp
    integer                                    :: iradf
    integer                                    :: mqn_m2PrR
    integer                                    :: ikpt
    integer                                    :: nmat
    integer                                    :: iband
    integer                                    :: lm
    integer                                    :: pMaxLocal
    integer                                    :: idirC
    integer                                    :: maxNobd
    integer :: nk

    ! Array variables
    integer,           allocatable             :: muKet(:, :)
    integer,           allocatable             :: muBra(:, :)
    real,              allocatable             :: varphiKet1(:, :, :, :)
    real,              allocatable             :: varphiKet2(:, :, :, :)
    real,              allocatable             :: varphiBra1(:, :, :, :)
    real,              allocatable             :: varphiBra2(:, :, :, :)
    complex,           allocatable             :: hPartNoAbcofBK(:, :, :, :, :)
    complex,           allocatable             :: varphiGrVeffVarphiNoAbcofBK(:, :, :, :, :)
    complex,           allocatable             :: overlapNoAbcofBK(:, :, :, :, :)
    complex,           allocatable             :: overlapNoAbcofBKDummy(:, :, :, :)
    integer,           allocatable             :: lambdaKet(:, :)
    integer,           allocatable             :: lambdaBra(:, :)
    complex,           allocatable             :: hGrVarphiNat(:, :, :, :, :)
    complex,           allocatable             :: overlapNat(:)
    complex,           allocatable             :: surfIntNat(:)
    complex,           allocatable             :: a(:, :, :)
    complex,           allocatable             :: b(:, :, :)
    complex,           allocatable             :: bascof_lo(:, :, :, :, :)
    integer,           allocatable             :: ngoprI(:)
    complex,           allocatable             :: ab0cof(:, :)
    complex,           allocatable             :: psiGrVeffsphPsi(:)
    complex,           allocatable             :: overlapNatDummy(:)
    complex                                    :: psiGrVeffsphPsiSum(3, 3)
    complex                                    :: surfIntNatSum(-1:1, 3)

    surfInt(:, :) = cmplx(0., 0.)

    maxNobd = maxval( nobd(:, :) )
    allocate( surfIntNat(maxNobd), overlapNat(maxNobd), overlapNatDummy(maxNobd), psiGrVeffsphPsi(maxNobd) )
    allocate( hPartNoAbcofBK(lmpMax, lmpMax, 3, atoms%nat, -1:1), overlapNoAbcofBK(lmpMax, lmpMax, 3, atoms%ntype, -1:1) )
    allocate( varphiGrVeffVarphiNoAbcofBK(lmpMax, lmpMax, 3, atoms%nat, 3), overlapNoAbcofBKDummy(lmpMax, lmpMax, 3, atoms%ntype)  )
    allocate( ab0cof(lmpMax, maxNobd) )
    allocate( a( MAXVAL(nv), 0:atoms%lmaxd*(atoms%lmaxd+2), atoms%nat), b(MAXVAL(nv), 0:atoms%lmaxd*(atoms%lmaxd+2), atoms%nat), &
      & bascof_lo(3, -atoms%llod:atoms%llod, 4 * atoms%llod + 2, atoms%nlod, atoms%nat) )
    allocate( varphiBra1(atoms%jmtd, 2, lmpMax, -1:1), varphiBra2(atoms%jmtd, 2, lmpMax, -1:1), &
            & varphiKet1(atoms%jmtd, 2, lmpMax, -1:1), varphiKet2(atoms%jmtd, 2, lmpMax, -1:1) )
    allocate( lambdaKet(2, 0:atoms%lmaxd), lambdaBra(2, 0:atoms%lmaxd) )
    allocate( muKet(-atoms%lmaxd:atoms%lmaxd, -1:1), muBra(-atoms%lmaxd:atoms%lmaxd, -1:1) )
    allocate( hGrVarphiNat(2, atoms%jmtd, (atoms%lmaxd + 2)**2, lmpMax, -1:1))
    allocate( ngoprI(atoms%nat) )

    chanMaxBra = 1
    chanMaxKet = 2
    lmaxBra = atoms%lmax(iDtypeB)
    muBra(:, :) = 0
    mqn_m2PrBra = -1
    varphiBra1 = 0.
    varphiBra2 = 0.
    do mqn_m = -atoms%lmax(iDtypeB), atoms%lmax(iDtypeB)
      muBra(mqn_m, -1) = mqn_m
    end do ! mqn_m
    do oqn_l = 0, atoms%lmax(iDtypeB)
      lambdaBra(1, oqn_l) = oqn_l
    end do ! oqn_l
    lmp = 0
    do oqn_l = 0, atoms%lmax(iDtypeB)
      do mqn_m = -oqn_l, oqn_l
        do iradf = 1, nRadFun(oqn_l, iDtypeB)
          lmp = lmp + 1
          varphiBra1(atoms%jri(iDtypeB), 1, lmp, -1) = varphi1(atoms%jri(iDtypeB), iradf, oqn_l )
          varphiBra2(atoms%jri(iDtypeB), 1, lmp, -1) = varphi2(atoms%jri(iDtypeB), iradf, oqn_l )
        end do ! iradf
      end do ! mqn_m
    end do ! oqn_l

    hGrVarphiNat(:, :, :, :, :) = cmplx(0., 0.)
    hPartNoAbcofBK(:, :, :, :, :) = cmplx(0., 0.)
    overlapNoAbcofBK(:, :, :, :, :) = cmplx(0., 0.)
    do mqn_m2PrR = -1, 1
      ! todo one does not need hGrVarphi
      ! todo check the mqn_m2PrKet
      call CalcHGrVarphi( atoms, iDtypeB, mqn_m2PrR, lmpMax, lmaxBra, grVarphiChMout, nRadFun, grVarPhiCh1, grVarPhiCh2,        &
                                                  & grVarphiChLout, vEff0NsphGrVarphi, El, lmpT, hGrVarphiNat )

    !todo adjust chanmaxket usw. parameter!!!!!!!!!!!!!!!!!!!!!NOTE
    !todo
    !todo
      call PrepareMTSurfIntDM( atoms, iDtypeB, chanMaxBra, chanMaxKet, lmaxBra, lmpMax, mqn_m2PrBra, mqn_m2PrR, grVarphiChMout,  muBra,    &
        & grVarphiChLout, lambdaBra, nRadFun, grVarphiCh1, grVarphiCh2, varphiBra1, varphiBra2, hGrVarphiNat(:, :, :, :, mqn_m2PrR), atoms%rmt(iDtypeB)**2, hPartNoAbcofBK(:, :, :, :, mqn_m2PrR), overlapNoAbcofBK(:, :, :, :, mqn_m2PrR) )
    end do ! mqn_m2PrR

    chanMaxBra = 1
    chanMaxKet = 1
    lmaxBra = atoms%lmax(iDtypeB)
    muBra(:, :) = 0
    muKet(:, :) = 0
    mqn_m2PrBra = -1
    mqn_m2PrKet = -1
    varphiBra1 = 0.
    varphiKet1 = 0.
    varphiBra2 = 0.
    varphiKet2 = 0.
    do mqn_m = -atoms%lmax(iDtypeB), atoms%lmax(iDtypeB)
      muKet(mqn_m, -1) = mqn_m
      muBra(mqn_m, -1) = mqn_m
    end do ! mqn_m
    do oqn_l = 0, atoms%lmax(iDtypeB)
      lambdaKet(1, oqn_l) = oqn_l
      lambdaBra(1, oqn_l) = oqn_l
    end do ! oqn_l
    lmp = 0
    do oqn_l = 0, atoms%lmax(iDtypeB)
      do mqn_m = -oqn_l, oqn_l
        do iradf = 1, nRadFun(oqn_l, iDtypeB)
          lmp = lmp + 1
          varphiBra1(atoms%jri(iDtypeB), 1, lmp, -1) = varphi1(atoms%jri(iDtypeB), iradf, oqn_l )
          varphiKet1(atoms%jri(iDtypeB), 1, lmp, -1) = varphi1(atoms%jri(iDtypeB), iradf, oqn_l )
          varphiBra2(atoms%jri(iDtypeB), 1, lmp, -1) = varphi2(atoms%jri(iDtypeB), iradf, oqn_l )
          varphiKet2(atoms%jri(iDtypeB), 1, lmp, -1) = varphi2(atoms%jri(iDtypeB), iradf, oqn_l )
        end do ! iradf
      end do ! mqn_m
    end do ! oqn_l
    varphiGrVeffVarphiNoAbcofBK(:, :, :, :, :) = cmplx(0., 0.)
    overlapNoAbcofBKDummy(:, :, :, :) = cmplx(0., 0.)
    do idirC = 1, 3
      call PrepareMTSurfIntDM( atoms, iDtypeB, chanMaxBra, chanMaxKet, lmaxBra, lmpMax, mqn_m2PrBra, mqn_m2PrKet, muKet, muBra,  &
        & lambdaKet, lambdaBra, nRadFun, varphiKet1, varphiKet2, varphiBra1, varphiBra2, r2grVeff0SphVarphi(:, :, :, :, idirC), 1., varphiGrVeffVarphiNoAbcofBK(:, :, :, :, idirC), overlapNoAbcofBKDummy(:, :, :, :) )
    end do ! idirC
    deallocate(overlapNoAbcofBKDummy)


      ! Note for q = 0 the matrices should be diagonal in the end for the dynamical matrix

      psiGrVeffsphPsiSum = cmplx(0., 0.)
      surfIntNatSum(:, :) = cmplx(0., 0.)
      ngoprI(:) = 1
      do ikpt = 1, kpts%nkpt

        nmat = nv(1, ikpt) + atoms%nlotot
        a(:, :, :) = cmplx(0.0, 0.0)
        b(:, :, :) = cmplx(0.0, 0.0)
        bascof_lo(:, :, :, :, :) = cmplx(0.0, 0.0)
        !call abcof3( atoms%lmaxd, atoms%ntype, atoms%nat, sym%nop, MAXVAL(nv), input%jspins, 1, atoms%lmaxd*(atoms%lmaxd+2), SIZE(z0(:,1,1,1)), &
         ! & atoms%llod, atoms%nlod, atoms%nlotot, sym%invtab, atoms%ntype, sym%mrot, ngoprI, atoms%taual, atoms%neq, atoms%lmax, &
          !& atoms%rmt, cell%omtil, cell%bmat, cell%bbmat, kpts%bk(:, ikpt), gbas(1, mapGbas(:nv(1, ikpt), ikpt, 1)), &
          !& gbas(2, mapGbas(:nv(1, ikpt), ikpt, 1)), gbas(3, mapGbas(:nv(1, ikpt), ikpt, 1)), nv(:, ikpt), nmat, &
          !& usdus%us, usdus%dus, usdus%uds, usdus%duds, usdus%ddn, sym%invsat, sym%invsatnr, usdus%ulos, usdus%uulon, usdus%dulon, &
          !& usdus%dulos, atoms%llo, atoms%nlo, atoms%l_dulo, atoms%lapw_l, kveclo(:,ikpt), odi, ods, a, b, bascof_lo )
          nk=fmpi%k_list(ikpt)
          CALL lapw%init(input, noco, nococonv, kpts, atoms, sym, nk, cell, .FALSE., fmpi)
          CALL abcof3(input, atoms, sym, 1, cell, kpts%bk(:, ikpt), lapw, &
                              usdus, oneD, 1, lapw%dim_nvd(), a, b, bascof_lo)

        ab0cof(:, :) = cmplx(0., 0.)
        do iband = 1, nobd(ikpt, 1)
          lmp = 0
          lm  = 0
          do oqn_l = 0, atoms%lmax(iDtypeB)
            do mqn_m = - oqn_l, oqn_l
              pMaxLocal = nRadFun(oqn_l, iDtypeB)
              ! p = 1
              ab0cof(lmp + 1, iband) = ImagUnit**oqn_l * dot_product( conjg(z0(:nv(1, ikpt), iband, ikpt, 1)), a(:nv(1, ikpt), lm, iDatomB) )
              ! p = 2
              ab0cof(lmp + 2, iband) = ImagUnit**oqn_l * dot_product( conjg(z0(:nv(1, ikpt), iband, ikpt, 1)), b(:nv(1, ikpt), lm, iDatomB) )
              ! Add LO contributions
              do iradf = 3, pMaxLocal
                ! p = 1
                ab0cof(lmp + 1, iband) = ab0cof(lmp + 1, iband) + ImagUnit**oqn_l * &
                  & dot_product( conjg(z0(nv(1, ikpt) + 1:nv(1, ikpt) + atoms%nlotot, iband, ikpt, 1)), &
                                                       & bascof_lo(1, mqn_m, :atoms%nlotot, iloTable(iradf, oqn_l, iDtypeB), iDatomB) )
                ! p = 2
                ab0cof(lmp + 2, iband) = ab0cof(lmp + 2, iband) + ImagUnit**oqn_l * &
                  & dot_product( conjg(z0(nv(1, ikpt) + 1:nv(1, ikpt) + atoms%nlotot, iband, ikpt, 1)), &
                                                       & bascof_lo(2, mqn_m, :atoms%nlotot, iloTable(iradf, oqn_l, iDtypeB), iDatomB) )
                ! 2 < p < LOs for that l and that atom type
                ab0cof(lmp + iradf, iband) = ImagUnit**oqn_l * dot_product( conjg(z0(nv(1, ikpt) + 1:nv(1, ikpt) + atoms%nlotot, iband, ikpt, 1)), &
                                                       & bascof_lo(3, mqn_m, :atoms%nlotot, iloTable(iradf, oqn_l, iDtypeB), iDatomB) )
              end do ! iradf


              ! This is a precalculation of the 1st and 3rd line of A.50. Actually the resulting lmp should have been primed or at
              ! least the p as lm1Pr = lm. But for sake of performance we place it here.
              ! sake of performance
              lm = lm + 1
              lmp = lmp + pMaxLocal
            end do ! mqn_m
          end do !oqn_l
        end do ! iband

        do idirC = 1, 3
          do mqn_m2PrR = -1, 1
            overlapNat(:) = cmplx(0., 0.)
            surfIntNat(:) = cmplx(0., 0.)
            call CalcSintMT( ikpt, ikpt, idirC, iDatomB, iDtypeB, lmpMax, nobd, hPartNoAbcofBK(:, :, :, :, mqn_m2PrR), &
              & overlapNoAbcofBK(:, :, :, :, mqn_m2PrR), ab0cof, ab0cof, overlapNat, surfIntNat )
            ! todo overhead here because overlap is calculated again
            ! todo array structures are not linearly run through yet
            overlapNatDummy(:) = cmplx(0., 0.)
            psiGrVeffsphPsi(:) = cmplx(0., 0.)
            call CalcSintMT( ikpt, ikpt, idirC, iDatomB, iDtypeB, lmpMax, nobd, varphiGrVeffVarphiNoAbcofBK(:, :, :, :, mqn_m2PrR + 2), &
              & overlapNoAbcofBK(:, :, :, :, mqn_m2PrR), ab0cof, ab0cof, overlapNatDummy, psiGrVeffsphPsi )
            do iband = 1, nobd(ikpt, 1)
              surfIntNatSum(mqn_m2PrR, idirC) = surfIntNatSum(mqn_m2PrR, idirC) + results%w_iks(iband, ikpt, 1) &
                                                                 & * ( surfIntNat(iband) - eig(iband, ikpt, 1) * overlapNat(iband) )
              psiGrVeffsphPsiSum(mqn_m2PrR + 2, idirC) = PsiGrVeffsphPsiSum(mqn_m2PrR + 2, idirC) &
                                                                      & + results%w_iks(iband, ikpt, 1) * ( psiGrVeffsphPsi(iband) )
            end do ! iband
          end do ! mqn_m2PrR
        end do ! idirC

      end do ! ikpt

      ! Transform the part involving gradients of basis functions into cartesian coordinates
      surfInt(1:3, 1:3) = matmul(Tmatrix0(1:3, 1:3), surfIntNatSum(-1:1, 1:3))

      ! Add term which has not been processed in natural coordinates
      surfInt(1:3, 1:3) = surfInt(1:3, 1:3) - psiGrVeffsphPsiSum(1:3, 1:3)

  end subroutine CalcSFintMTPsiHepsGradPsi

  subroutine CalcSFintIRPsi1HepsPsi( atoms, input, stars, cell, results, Veff0, kpts, qpts, ngdp, iqpt, iDtypeB, iDatomB, iDatomA, nobd,&
      & gdp, mapKpq2K, kpq2kPrVec, nv, mapGbas, gBas, z0, eig, surfInt, testGoldstein )

    use m_types
    use m_dfpt_init, only : ConvertStar2G

    implicit none

    ! Type parameters
    type(t_atoms),                 intent(in)  :: atoms
    type(t_input),                 intent(in)  :: input
    type(t_stars),                 intent(in)  :: stars
    type(t_cell),                  intent(in)  :: cell
    type(t_results),               intent(in)  :: results
    type(t_potden),             intent(in)  :: Veff0
    type(t_kpts),                  intent(in)  :: kpts
    type(t_kpts),                  intent(in)  :: qpts

    ! Scalar parameters
    integer,                       intent(in)  :: ngdp
    integer,                       intent(in)  :: iqpt
    integer,                       intent(in)  :: iDtypeB
    integer,                       intent(in)  :: iDatomB
    integer,                       intent(in)  :: iDatomA
    logical,                       intent(in)  :: testGoldstein

    ! Array parameters
    integer,                       intent(in)  :: nobd(:, :)
    integer,                       intent(in)  :: gdp(:, :)
    integer,                       intent(in)  :: mapKpq2K(:, :)
    integer,                       intent(in)  :: kpq2kPrVec(:, :, :)
    integer,                       intent(in)  :: nv(:, :)
    integer,                       intent(in)  :: mapGbas(:, :, :)
    integer,                       intent(in)  :: gBas(:, :)
    COMPLEX,                      intent(in)  :: z0(:, :, :, :)
    real,                          intent(in)  :: eig(:, :, :)
    complex,                       intent(inout) :: surfInt(:, :)

    ! Scalar variables
    integer                                    :: maxNobd
    integer                                    :: coScale
    integer                                    :: ikpt
    integer                                    :: iBas
    integer                                    :: ikpq
    integer                                    :: idirC
    integer                                    :: idirR
    integer                                    :: iband

    ! Array variables
    complex,           allocatable             :: surfIntOvl(:, :)
    complex,           allocatable             :: surfIntVeff0(:, :, :)
    complex,           allocatable             :: surfIntKinEnerg(:, :)
    complex,           allocatable             :: vpw_eff_uw(:)
    complex,           allocatable             :: veffUvIR(:, :)
    real,              allocatable             :: gBasBra(:, :)
    real,              allocatable             :: gBasKet(:, :)
    real,              allocatable             :: gBasBra2(:, :)
    complex,           allocatable             :: z1nG(:, :, :, :)
    complex,           allocatable             :: z1nGContainer(:, :, :)
    real                                       :: Gext(3)
    real                                       :: kExt(3)
    complex                                    :: overlapTest(3, 3)

    maxNobd = maxval( nobd(:, :) )
    coScale = 1

    allocate( surfIntOvl( maxNobd, maxNobd ), surfIntVeff0( maxNobd, 3, 3 ), surfIntKinEnerg( maxNobd, maxNobd ) )
    allocate( vpw_eff_uw( ngdp ) )
    allocate( veffUvIR(3, (coScale * atoms%lmaxd + 1)**2) )
    allocate( z1nG(SIZE(z0(:,1,1,1)), 3, atoms%nat, maxNobd), z1nGContainer( SIZE(z0(:,1,1,1)), maxNobd, 3) )
    allocate( gBasKet(3, MAXVAL(nv)), gBasBra(3, MAXVAL(nv)), gBasBra2(3, MAXVAL(nv))  )

    surfInt(:, :) = cmplx(0., 0.)
    surfIntOvl(:, :) = cmplx(0., 0.)
    surfIntVeff0(:, :, :) = cmplx(0., 0.)
    overlapTest(:, :) = cmplx(0., 0.)

    ! Generate plane-wave expanded potential
    vpw_eff_uw = cmplx(0., 0.)
    call convertStar2G( Veff0%pw(:, 1), vpw_eff_uw, stars, ngdp, gdp )

    ! Calculate k-independent parts
    veffUvIR(:, :) = cmplx(0., 0.)
    call IRcoeffVeffUv( atoms, stars, cell, iDtypeB, iDatomB, ngdp, coScale, gdp, vpw_eff_uw, veffUvIR )

    do ikpt = 1, kpts%nkpt
      ikpq = mapKpq2K(ikpt, iqpt)

      ! Read in z1 for k-point ikpq
      z1nGContainer = cmplx(0.0, 0.0)
      if (testGoldstein) then
        kExt(1:3) = matmul(cell%bmat, kpts%bk(1:3, ikpt))
        gExt(:) = 0.
        do iBas = 1, nv(1, ikpt) !+ atoms%nlotot
          gExt(1:3) = matmul( cell%bmat(1:3, 1:3), gBas(1:3, mapGbas(iBas, ikpt, 1)))
          do iband = 1, nobd(ikpt, 1)
            do idirR = 1, 3
              z1nGContainer(iBas, iband, idirR) = -ImagUnit * ( kExt(idirR) + gExt(idirR) ) * z0(iBas, iband, ikpt, 1)
            end do ! idirR
          end do ! iband
        end do ! iBas
      else
        call ReadInz1( atoms, ikpt, iqpt, ikpq, nobd, nv, z1nG )
      end if

      ! Generate basis set for kinetic energy action, ket wave function (the same) and bra ket function
      gBasBra(:, :) = 0.
      gBasKet(:, :) = 0.
      gBasBra2(:, :) = 0.
      ! NOTE: One has to care which k is here
      do iBas = 1, nv(1, ikpt)
        gBasKet(1:3, iBas) = real( gBas(1:3, mapGbas(iBas, ikpt, 1)) )
      end do ! iBas
      ! NOTE: One has to care which k is here
      do iBas = 1, nv(1, ikpq)
        gBasBra(1:3, iBas) = real( gBas(1:3, mapGbas(iBas, ikpq, 1)) + qpts%bk(1:3, iqpt) + kpq2kPrVec(1:3, ikpt, iqpt) )
      end do ! iBas
      do iBas = 1, nv(1, ikpq)
        gBasBra2(1:3, iBas) = real( gBas(1:3, mapGbas(iBas, ikpq, 1)) )
      end do ! iBas

      ! Rearrange z1 array
      if (.not.testGoldstein) then
        ! this is probably wrong with the bands it should be at ikpt
      z1nGContainer(:, :, :) = cmplx(0., 0.)
      do idirC = 1, 3
        do iBas = 1, nv(1, ikpq) + atoms%nlotot
          do iband = 1, nobd(ikpt, 1)
            z1nGContainer(iBas, iBand, idirC) = z1nG(iBas, idirC, iDatomA, iband)
          end do ! iband
        end do ! iBas
      end do ! idirC
    end if

      ! Calculate surface integral with effective potential
      surfIntVeff0(:, :, :) = cmplx(0., 0.)
      do idirC = 1, 3
      !todo have a look onto the directions
          call calcSurfVeff( atoms, input, kpts, cell, ikpq, ikpt, coScale, gBas, veffUvIR, nv, nobd(ikpt, 1), &
            & mapGbas, z1nGContainer(:, :, idirC), z0(:, :, ikpt, 1), iDtypeB, iDatomB, surfIntVeff0(:, :, idirC) )
      end do ! idirC

      do idirC = 1, 3
        do idirR = 1, 3
          surfIntKinEnerg(:, :) = cmplx(0., 0.)
          surfIntOvl(:, :) = cmplx(0., 0.)
            ! Calculate surface integral with kinetic energy and overlap
            call calcSintKinEnergOvl( atoms, cell, kpts, results, iDtypeB, iDatomB, ikpq, ikpt, idirR, &
              & nv(1, :), nobd(ikpt, 1), gBasBra, gBasKet, gBasKet, gBasBra2, z1nGContainer(:, :, idirC), z0(:, :, ikpt, 1), surfIntKinEnerg, &
              & surfIntOvl )
          do iband = 1, nobd(ikpt, 1)
              !!!latest:
              surfInt(idirR, idirC) = surfInt(idirR, idirC) + 2 * results%w_iks(iband, ikpt, 1) * ( surfIntKinEnerg(iband, iband) &
                                                   & + 0.0*surfIntVeff0(iband, idirR, idirC) - eig(iband, ikpt, 1) * surfIntOvl(iband, iband) )
!              surfInt(idirR, idirC) = surfInt(idirR, idirC) + 2 * results%w_iks(iband, ikpt, 1) * ( surfIntKinEnerg(iband)&
!                                                   & - eig(iband, ikpt, 1) * surfIntOvl(iband) )
!              surfInt(idirR, idirC) = surfInt(idirR, idirC) + 2 * results%w_iks(iband, ikpt, 1) * ( &
!                                                   & surfIntVeff0(iband, idirR, idirC) - eig(iband, ikpt, 1) * surfIntOvl(iband))
          end do ! iband
        end do ! idirR
      end do ! idirC
    end do ! ikpt

    if ( .false. ) then
      write(*, '(a)') 'IR Surface integral routine z1'
      write(*, '(3(2(es16.8,1x),3x))') surfInt(1, :)
      write(*, '(3(2(es16.8,1x),3x))') surfInt(2, :)
      write(*, '(3(2(es16.8,1x),3x))') surfInt(3, :)

      write(*, '(a)') 'IR Surface integral routineo'
      write(*, '(3(2(es16.8,1x),3x))') overlapTest(1, :)
      write(*, '(3(2(es16.8,1x),3x))') overlapTest(2, :)
      write(*, '(3(2(es16.8,1x),3x))') overlapTest(3, :)
    end if

  end subroutine CalcSFintIRPsi1HepsPsi

  subroutine CalcSFintIRPsiHepsPsi1( atoms, input, stars, cell, results, Veff0, kpts, qpts, ngdp, iqpt, iDtypeB, iDatomB, iDatomA,&
      & nobd, gdp, mapKpq2K, kpq2kPrVec, nv, mapGbas, gBas, z0, eig, surfInt, testGoldstein )

    use m_types
    use m_dfpt_init, only : ConvertStar2G

    implicit none

    ! Type parameters
    type(t_atoms),                 intent(in)  :: atoms
    type(t_input),                 intent(in)  :: input
    type(t_stars),                 intent(in)  :: stars
    type(t_cell),                  intent(in)  :: cell
    type(t_results),               intent(in)  :: results
    type(t_potden),             intent(in)  :: Veff0
    type(t_kpts),                  intent(in)  :: kpts
    type(t_kpts),                  intent(in)  :: qpts

    ! Scalar parameters
    integer,                       intent(in)  :: ngdp
    integer,                       intent(in)  :: iqpt
    integer,                       intent(in)  :: iDtypeB
    integer,                       intent(in)  :: iDatomB
    integer,                       intent(in)  :: iDatomA

    ! Array parameters
    integer,                       intent(in)  :: nobd(:, :)
    integer,                       intent(in)  :: gdp(:, :)
    integer,                       intent(in)  :: mapKpq2K(:, :)
    integer,                       intent(in)  :: kpq2kPrVec(:, :, :)
    integer,                       intent(in)  :: nv(:, :)
    integer,                       intent(in)  :: mapGbas(:, :, :)
    integer,                       intent(in)  :: gBas(:, :)
    COMPLEX,                      intent(in)  :: z0(:, :, :, :)
    real,                          intent(in)  :: eig(:, :, :)
    logical,                       intent(in)  :: testGoldstein
    complex,                       intent(out) :: surfInt(:, :)

    ! Scalar variables
    integer                                    :: maxNobd
    integer                                    :: coScale
    integer                                    :: ikpt
    integer                                    :: iBas
    integer                                    :: ikpq
    integer                                    :: idirC
    integer                                    :: idirR
    integer                                    :: iband

    ! Array variables
    complex,           allocatable             :: surfIntOvl(:, :)
    complex,           allocatable             :: surfIntVeff0(:, :, :)
    complex,           allocatable             :: surfIntKinEnerg(:, :)
    complex,           allocatable             :: vpw_eff_uw(:)
    complex,           allocatable             :: veffUvIR(:, :)
    real,              allocatable             :: gBasBra(:, :)
    real,              allocatable             :: gBasKet(:, :)
    real,              allocatable             :: gBasKinEnerg(:, :)
    complex,           allocatable             :: z1nG(:, :, :, :)
    complex,           allocatable             :: z1nGContainer(:, :, :)
    real                                       :: Gext(3)
    real                                       :: kExt(3)
    complex                                    :: overlapTest(3, 3)

    maxNobd = maxval( nobd(:, :) )
    coScale = 1

    allocate( surfIntOvl( maxNobd, maxNobd ), surfIntVeff0( maxNobd, 3, 3 ), surfIntKinEnerg( maxNobd, maxNobd ) )
    allocate( vpw_eff_uw( ngdp ) )
    allocate( veffUvIR(3, (coScale * atoms%lmaxd + 1)**2) )
    allocate( z1nG(SIZE(z0(:,1,1,1)), 3, atoms%nat, maxNobd), z1nGContainer( SIZE(z0(:,1,1,1)), maxNobd, 3) )
    allocate( gBasKet(3, MAXVAL(nv)), gBasBra(3, MAXVAL(nv)), gBasKinEnerg(3, MAXVAL(nv)) )

    surfInt(:, :) = cmplx(0., 0.)
    surfIntOvl(:, :) = cmplx(0., 0.)
    surfIntVeff0(:, :, :) = cmplx(0., 0.)
    overlapTest(:, :) = cmplx(0., 0.)

    ! Generate plane-wave expanded potential
    vpw_eff_uw = cmplx(0., 0.)
    call convertStar2G( Veff0%pw(:, 1), vpw_eff_uw, stars, ngdp, gdp )

    ! Calculate k-independent parts
    veffUvIR(:, :) = cmplx(0., 0.)
    call IRcoeffVeffUv( atoms, stars, cell, iDtypeB, iDatomB, ngdp, coScale, gdp, vpw_eff_uw, veffUvIR )

    do ikpt = 1, kpts%nkpt
! todo read in the z1
      ikpq = mapKpq2K(ikpt, iqpt)

      ! Read in z1 for k-point ikpq
      z1nGContainer = cmplx(0.0, 0.0)
      if ( testGoldstein) then
        kExt(1:3) = matmul(cell%bmat, kpts%bk(1:3, ikpt))
        gExt(:) = 0.
        do iBas = 1, nv(1, ikpt) !+ atoms%nlotot
          gExt(1:3) = matmul( cell%bmat(1:3, 1:3), gBas(1:3, mapGbas(iBas, ikpt, 1)))
          do iband = 1, nobd(ikpt, 1)
            do idirR = 1, 3
              z1nGContainer(iBas, iband, idirR) = -ImagUnit * ( kExt(idirR) + gExt(idirR) ) * z0(iBas, iband, ikpt, 1)
            end do ! idirR
          end do ! iband
        end do ! iBas
      else
        call ReadInz1( atoms, ikpt, iqpt, ikpq, nobd, nv, z1nG )
      end if

      ! Generate basis set for kinetic energy action, ket wave function (the same) and bra ket function
      gBasBra(:, :) = 0.
      gBasKet(:, :) = 0.
      gBasKinEnerg(:, :) = 0.
      ! NOTE: One has to care which k is here
      do iBas = 1, nv(1, ikpq)
        gBasKet(1:3, iBas) = real( gBas(1:3, mapGbas(iBas, ikpq, 1)) + qpts%bk(1:3, iqpt) + kpq2kPrVec(1:3, ikpt, iqpt) )
      end do ! iBas
      ! NOTE: One has to care which k is here
      do iBas = 1, nv(1, ikpt)
        gBasBra(1:3, iBas) = real( gBas(1:3, mapGbas(iBas, ikpt, 1)) )
      end do ! iBas
      do iBas = 1, nv(1, ikpq)
        gBasKinEnerg(1:3, iBas) = real( gBas(1:3, mapGbas(iBas, ikpq, 1)) )
      end do ! iBas

      ! Rearrange z1 array
      if ((.not.testGoldstein)) then
      z1nGContainer(:, :, :) = cmplx(0., 0.)
      do idirC = 1, 3
        do iBas = 1, nv(1, ikpq) + atoms%nlotot
          do iband = 1, nobd(ikpt, 1)
            z1nGContainer(iBas, iBand, idirC) = z1nG(iBas, idirC, iDatomA, iband)
          end do ! iband
        end do ! iBas
      end do ! idirC
      end if

      ! Calculate surface integral with effective potential
      surfIntVeff0(:, :, :) = cmplx(0., 0.)
      do idirC = 1, 3
      !todo have a look onto the directions
          call calcSurfVeff( atoms, input, kpts, cell, ikpt, ikpq, coScale, gBas, veffUvIR, nv, nobd(ikpt, 1), &
            & mapGbas, z0(:, :, ikpt, 1), z1nGContainer(:, :, idirC), iDtypeB, iDatomB, surfIntVeff0(:, :, idirC) )

      end do ! idirC

      do idirC = 1, 3
        do idirR = 1, 3
          surfIntKinEnerg(:, :) = cmplx(0., 0.)
          surfIntOvl(:, :) = cmplx(0., 0.)
          ! Calculate surface integral with kinetic energy and overlap
          call calcSintKinEnergOvl( atoms, cell, kpts, results, iDtypeB, iDatomB, ikpt, ikpq, idirR, &
            & nv(1, :), nobd(ikpt, 1), gBasBra, gBasKinEnerg, gBasKet, gBasBra, z0(:, :, ikpt, 1), z1nGContainer(:, :, idirC), surfIntKinEnerg, &
            & surfIntOvl )
          do iband = 1, nobd(ikpt, 1)
              !!!latest:
              surfInt(idirR, idirC) = surfInt(idirR, idirC) + 2 * results%w_iks(iband, ikpt, 1) * ( surfIntKinEnerg(iband, iband) &
                                                     & + 0.0*surfIntVeff0(iband, idirR, idirC) - eig(iband, ikpt, 1) * surfIntOvl(iband, iband) )
              !surfInt(idirR, idirC) = surfInt(idirR, idirC) + 2 * results%w_iks(iband, ikpt, 1) * (  surfIntKinEnerg(iband)&
              !                                       &  - eig(iband, ikpt, 1) * surfIntOvl(iband) )
            !  surfInt(idirR, idirC) = surfInt(idirR, idirC) + 2 * results%w_iks(iband, ikpt, 1) * ( &
            !                                       & surfIntVeff0(iband, idirR, idirC) - eig(iband, ikpt, 1) * surfIntOvl(iband))
          end do ! iband
        end do ! idirR
      end do ! idirC
    end do ! ikpt

    if ( .false. ) then
      write(*, '(a)') 'IR Surface integral routine'
      write(*, '(3(2(es16.8,1x),3x))') surfInt(1, :)
      write(*, '(3(2(es16.8,1x),3x))') surfInt(2, :)
      write(*, '(3(2(es16.8,1x),3x))') surfInt(3, :)

      write(*, '(a)') 'IR Surface integral routineo'
      write(*, '(3(2(es16.8,1x),3x))') overlapTest(1, :)
      write(*, '(3(2(es16.8,1x),3x))') overlapTest(2, :)
      write(*, '(3(2(es16.8,1x),3x))') overlapTest(3, :)
    end if

  end subroutine CalcSFintIRPsiHepsPsi1

  subroutine CalcSFintMTPsi1HepsPsiAndPsiHepsPsi1ExpCoeffVar( fmpi, noco, nococonv, oneD, atoms, input, sym, usdus, kpts, cell, results, iqpt, iDtypeB, iDatomB, iDatomA, lmpMax, nRadFun, eig, varphi1, varphi2, hVarphi, &
    & mapKpq2K, gBas, mapGbas, nv, kveclo, z0, lmpT, nobd, iloTable, surfInt, testGoldstein )

    use m_types
    use m_abcof3
    use m_types_oneD, only : od_inp, od_sym

    implicit none

    ! Type parameters
    type(t_mpi),                  intent(in)  :: fmpi
    type(t_nococonv),                intent(in)  :: nococonv
    type(t_noco),                intent(in)  :: noco
    type(t_oneD),                intent(in)  :: oneD
    type(t_atoms),                  intent(in)  :: atoms
    type(t_input),                  intent(in)  :: input
    type(t_sym),                    intent(in)  :: sym
    type(t_usdus),                  intent(in)  :: usdus
    type(t_kpts),                   intent(in)  :: kpts
    type(t_cell),                   intent(in)  :: cell
    type(t_results),                intent(in)  :: results

    ! Scalar parameter
    integer,                        intent(in)  :: iqpt
    integer,                        intent(in)  :: lmpMax
    integer,                        intent(in)  :: iDtypeB
    integer,                        intent(in)  :: iDatomB
    integer,                        intent(in)  :: iDatomA
    logical,                        intent(in)  :: testGoldstein

    ! Array parameters
    integer,                        intent(in)  :: nRadFun(0:, :)
    real,                           intent(in)  :: eig(:,:,:)
    complex,                        intent(in)  :: hVarphi(:, :, :, :)
    integer,                        intent(in)  :: mapKpq2K(:, :)
    integer,                        intent(in)  :: gbas(:, :)
    integer,                        intent(in)  :: mapGbas(:, :, :)
    integer,                        intent(in)  :: nv(:, :)
    integer,                        intent(in)  :: kveclo(:,:)
    COMPLEX,                       intent(in)  :: z0(:,:,:,:)
    integer,                        intent(in)  :: nobd(:, :)
    integer,                        intent(in)  :: lmpT(:)
    integer,                        intent(in)  :: iloTable(:, 0:, :)
    real,                           intent(in)  :: varphi1(:, :, 0:)
    real,                           intent(in)  :: varphi2(:, :, 0:)
    complex,                        intent(out) :: surfInt(:, :)

    ! Type variables
    type(od_inp)                                :: odi
    type(od_sym)                                :: ods
    type(t_lapw) :: lapw

    ! Scalar variables
    integer                                     :: oqn_l
    integer                                     :: iradf
    integer                                     :: imesh
    integer                                     :: mqn_m
    integer                                     :: lm
    integer                                     :: nmat
    integer                                     :: chanMaxBra
    integer                                     :: chanMaxKet
    integer                                     :: lmaxBra
    integer                                     :: mqn_m2PrBra
    integer                                     :: mqn_m2PrKet
    integer                                     :: lmp
    integer                                     :: ikpt
    integer                                     :: ikpq
    integer                                     :: maxNobd
    integer                                     :: idirC
    integer                                     :: idirR
    integer                                     :: iband
    integer                                     :: pmaxLocal
    integer :: nk

    ! Array variables
    integer,           allocatable              :: muKet(:, :)
    integer,           allocatable              :: muBra(:, :)
    complex,           allocatable              :: hFullNoAbcofBK(:, :, :, :)
    complex,           allocatable              :: overlapNoAbcofBK(:, :, :, :)
    real,              allocatable              :: varphiBra1(:, :, :, :)
    real,              allocatable              :: varphiKet1(:, :, :, :)
    real,              allocatable              :: varphiBra2(:, :, :, :)
    real,              allocatable              :: varphiKet2(:, :, :, :)
    integer,           allocatable              :: lambdaKet(:, :)
    integer,           allocatable              :: lambdaBra(:, :)
    integer,           allocatable              :: ngoprI(:)
    complex,           allocatable              :: a(:, :, :)
    complex,           allocatable              :: b(:, :, :)
    complex,           allocatable              :: bascof_lo(:, :, :, :, :)
    complex,           allocatable              :: aKpq(:, :, :)
    complex,           allocatable              :: bKpq(:, :, :)
    complex,           allocatable              :: bascof_loKpq(:, :, :, :, :)
    complex,           allocatable              :: ab1cofVec(:, :, :)
    complex,           allocatable              :: ab0cof(:, :)
    complex,           allocatable              :: z1nG(:, :, :, :)
    complex,           allocatable              :: z1Analytical(:, :, :, :)
    complex,           allocatable              :: surfIntMT(:)
    complex,           allocatable              :: overlap(:)
    real                                        :: Gext(3)
    real                                        :: kExt(3)
    integer :: iBas


    ! Quantities for initialization
    maxNobd = maxval(nobd(:, :))
    allocate( muKet(-atoms%lmaxd:atoms%lmaxd, -1:1), muBra(-atoms%lmaxd:atoms%lmaxd, -1:1) )
    allocate( hFullNoAbcofBK(lmpMax, lmpMax, 3, atoms%nat), overlapNoAbcofBK(lmpMax, lmpMax, 3, atoms%ntype) )
    allocate( varphiBra1(atoms%jmtd, 2, lmpMax, -1:1), varphiBra2(atoms%jmtd, 2, lmpMax, -1:1), &
            & varphiKet1(atoms%jmtd, 2, lmpMax, -1:1), varphiKet2(atoms%jmtd, 2, lmpMax, -1:1) )
    allocate( lambdaKet(2, 0:atoms%lmaxd), lambdaBra(2, 0:atoms%lmaxd) )
    allocate( ngoprI(atoms%nat) )
    allocate( a( MAXVAL(nv), 0:atoms%lmaxd*(atoms%lmaxd+2), atoms%nat), b(MAXVAL(nv), 0:atoms%lmaxd*(atoms%lmaxd+2), atoms%nat), &
      & bascof_lo(3, -atoms%llod:atoms%llod, 4 * atoms%llod + 2, atoms%nlod, atoms%nat) )
    allocate( aKpq( MAXVAL(nv), 0:atoms%lmaxd*(atoms%lmaxd+2), atoms%nat), bKpq(MAXVAL(nv), 0:atoms%lmaxd*(atoms%lmaxd+2), atoms%nat), &
      & bascof_loKpq(3, -atoms%llod:atoms%llod, 4 * atoms%llod + 2, atoms%nlod, atoms%nat) )
    allocate( ab1cofVec(lmpMax, maxNobd, 3), ab0cof(lmpMax, maxNobd) )
    allocate( z1nG(SIZE(z0(:,1,1,1)), 3, atoms%nat, maxNobd), z1Analytical(SIZE(z0(:,1,1,1)), 3, atoms%nat, maxNobd) )
    allocate( surfIntMT(maxNobd), overlap(maxNobd) )


    hFullNoAbcofBK(:, :, :, :) = cmplx(0., 0.)
    overlapNoAbcofBK(:, :, :, :) = cmplx(0., 0.)
    chanMaxBra = 1
    chanMaxKet = 1
    lmaxBra = atoms%lmax(iDtypeB)
    muBra(:, :) = 0
    muKet(:, :) = 0
    mqn_m2PrBra = -1
    mqn_m2PrKet = -1
    varphiBra1 = 0.
    varphiKet1 = 0.
    varphiBra2 = 0.
    varphiKet2 = 0.
    do mqn_m = -atoms%lmax(iDtypeB), atoms%lmax(iDtypeB)
      muKet(mqn_m, -1) = mqn_m
      muBra(mqn_m, -1) = mqn_m
    end do ! mqn_m
    do oqn_l = 0, atoms%lmax(iDtypeB)
      lambdaKet(1, oqn_l) = oqn_l
      lambdaBra(1, oqn_l) = oqn_l
    end do ! oqn_l
    lmp = 0
    do oqn_l = 0, atoms%lmax(iDtypeB)
      do mqn_m = -oqn_l, oqn_l
        do iradf = 1, nRadFun(oqn_l, iDtypeB)
          lmp = lmp + 1
          varphiBra1(atoms%jri(iDtypeB), 1, lmp, -1) = varphi1(atoms%jri(iDtypeB), iradf, oqn_l )
          varphiKet1(atoms%jri(iDtypeB), 1, lmp, -1) = varphi1(atoms%jri(iDtypeB), iradf, oqn_l )
          varphiBra2(atoms%jri(iDtypeB), 1, lmp, -1) = varphi2(atoms%jri(iDtypeB), iradf, oqn_l )
          varphiKet2(atoms%jri(iDtypeB), 1, lmp, -1) = varphi2(atoms%jri(iDtypeB), iradf, oqn_l )
        end do ! iradf
      end do ! mqn_m
    end do ! oqn_l

    call PrepareMTSurfIntDM( atoms, iDtypeB, chanMaxBra, chanMaxKet, lmaxBra, lmpMax, mqn_m2PrBra, mqn_m2PrKet, muKet, muBra,     &
      & lambdaKet, lambdaBra, nRadFun, varphiKet1, varphiKet2, varphiBra1, varphiBra2, hVarphi, atoms%rmt(iDtypeB)**2, hFullNoAbcofBK, overlapNoAbcofBK )

    surfInt(:, :) = cmplx(0., 0.)
    ngoprI(:) = 1
    do ikpt = 1, kpts%nkpt

      ikpq = mapKpq2K(ikpt, iqpt)

      ! Read in z1 for k-point ikpq
      if (.not.testGoldstein) then
        z1nG = cmplx(0.0, 0.0)
        call ReadInz1( atoms, ikpt, iqpt, ikpq, nobd, nv, z1nG )
      else
        z1Analytical(:, :, :, :) = cmplx(0., 0.)
        kExt(1:3) = matmul(cell%bmat, kpts%bk(1:3, ikpt))
        gExt(:) = 0.
        do iBas = 1, nv(1, ikpt) !+ atoms%nlotot
          gExt(1:3) = matmul( cell%bmat(1:3, 1:3), gBas(1:3, mapGbas(iBas, ikpt, 1)))
          do iband = 1, nobd(ikpt, 1)
            do idirR = 1, 3
              !gradz0(iBas, iband, idirR) = iu * ( kExt(idirR) + gExt(idirR) ) * z0(iBas, iband, ikpt, 1)
              z1Analytical(iBas, idirR, iDatomA, iband) = -ImagUnit * ( kExt(idirR) + gExt(idirR) ) * z0(iBas, iband, ikpt, 1)
            end do ! idirR
          end do ! iband
        end do ! iBas
      end if

      ! Calculate the basis matching coefficients at k + q = k' and the matching coefficients at k.
      nmat = nv(1, ikpq) + atoms%nlotot
      aKpq(:, :, :) = cmplx(0.0, 0.0)
      bKpq(:, :, :) = cmplx(0.0, 0.0)
      bascof_loKpq(:, :, :, :, :) = cmplx(0.0, 0.0)
      !call abcof3( atoms%lmaxd, atoms%ntype, atoms%nat, sym%nop, MAXVAL(nv), input%jspins, 1, atoms%lmaxd*(atoms%lmaxd+2), SIZE(z0(:,1,1,1)), &
        !& atoms%llod, atoms%nlod, atoms%nlotot, sym%invtab, atoms%ntype, sym%mrot, ngoprI, atoms%taual, atoms%neq, atoms%lmax, &
        !& atoms%rmt, cell%omtil, cell%bmat, cell%bbmat, kpts%bk(:, ikpq), gbas(1, mapGbas(:nv(1, ikpq), ikpq, 1)), &
        !& gbas(2, mapGbas(:nv(1, ikpq), ikpq, 1)), gbas(3, mapGbas(:nv(1, ikpq), ikpq, 1)), nv(:, ikpq), nmat, &
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
      !call abcof3( atoms%lmaxd, atoms%ntype, atoms%nat, sym%nop, MAXVAL(nv), input%jspins, 1, atoms%lmaxd*(atoms%lmaxd+2), SIZE(z0(:,1,1,1)), &
        !& atoms%llod, atoms%nlod, atoms%nlotot, sym%invtab, atoms%ntype, sym%mrot, ngoprI, atoms%taual, atoms%neq, atoms%lmax, &
        !& atoms%rmt, cell%omtil, cell%bmat, cell%bbmat, kpts%bk(:, ikpt), gbas(1, mapGbas(:nv(1, ikpt), ikpt, 1)), &
        !& gbas(2, mapGbas(:nv(1, ikpt), ikpt, 1)), gbas(3, mapGbas(:nv(1, ikpt), ikpt, 1)), nv(:, ikpt), nmat, &
        !& usdus%us, usdus%dus, usdus%uds, usdus%duds, usdus%ddn, sym%invsat, sym%invsatnr, usdus%ulos, usdus%uulon, usdus%dulon, &
        !& usdus%dulos, atoms%llo, atoms%nlo, atoms%l_dulo, atoms%lapw_l, kveclo(:,ikpt), odi, ods, a, b, bascof_lo )
        nk=fmpi%k_list(ikpt)
        CALL lapw%init(input, noco, nococonv, kpts, atoms, sym, nk, cell, .FALSE., fmpi)
        CALL abcof3(input, atoms, sym, 1, cell, kpts%bk(:, ikpt), lapw, &
                            usdus, oneD, 1, lapw%dim_nvd(), a, b, bascof_lo)

      ! Calculate the vector-like large matching coefficients using the first-order wave-function expansion coefficients gained
      ! from solving the Sternheimer equation.
      if (.not.testGoldstein) then
      ab1cofVec(:, :, :) = cmplx(0.0, 0.0)
      do idirC = 1, 3
        do iband = 1, nobd(ikpt, 1)
          lm  = 0
          lmp = 0
          do oqn_l = 0, atoms%lmax(iDtypeB)
            do mqn_m = -oqn_l, oqn_l
              pMaxLocal = nRadFun(oqn_l, iDtypeB)
              ! p = 1
              ab1cofVec(lmp + 1, iband, idirC) = ImagUnit**oqn_l * &
                             & dot_product( conjg( z1nG(:nv(1, ikpq), idirC, iDatomA, iband) ), aKpq(:nv(1, ikpq), lm, iDatomB) )
              ! p = 2
              ab1cofVec(lmp + 2, iband, idirC) = ImagUnit**oqn_l * &
                             & dot_product( conjg( z1nG(:nv(1, ikpq), idirC, iDatomA, iband) ), bKpq(:nv(1, ikpq), lm, iDatomB) )
              do iradf = 3, pMaxLocal
                ! p = 1
                ab1cofVec(lmp + 1, iband, idirC) = ab1cofVec(lmp + 1, iband, idirC) &
                  & + ImagUnit**oqn_l * dot_product( conjg( z1nG(nv(1, ikpq) + 1:nv(1, ikpq) + atoms%nlotot, idirC, iDatomA, iband) ), &
                                              & bascof_loKpq(1, mqn_m, :atoms%nlotot, iloTable(iradf, oqn_l, iDtypeB), iDatomB) )
                ! p = 2
                ab1cofVec(lmp + 2, iband, idirC) = ab1cofVec(lmp + 2, iband, idirC) &
                  & + ImagUnit**oqn_l * dot_product(conjg(z1nG(nv(1, ikpq) + 1:nv(1, ikpq) + atoms%nlotot, idirC, iDatomA, iband)), &
                                              & bascof_loKpq(2, mqn_m, :atoms%nlotot, iloTable(iradf, oqn_l, iDtypeB), iDatomB) )
                ! 2 < p < LOs for that l and that atom type
                ab1cofVec(lmp + iradf, iband, idirC) = &
                  & ImagUnit**oqn_l * dot_product( conjg( z1nG(nv(1, ikpq) + 1:nv(1, ikpq) + atoms%nlotot, idirC, iDatomA, iband) ),&
                                              & bascof_loKpq(3, mqn_m, :atoms%nlotot, iloTable(iradf, oqn_l, iDtypeB), iDatomB) )
              end do ! iradf
              lm = lm + 1
              lmp = lmp + nRadFun(oqn_l, iDtypeB)
            end do ! mqn_m
          end do ! oqn_l
        end do ! iband
      end do ! idirC
    else
      ab1cofVec(:, :, :) = cmplx(0.0, 0.0)
      do idirC = 1, 3
        do iband = 1, nobd(ikpt, 1)
          lm  = 0
          lmp = 0
          do oqn_l = 0, atoms%lmax(iDtypeB)
            do mqn_m = -oqn_l, oqn_l
              pMaxLocal = nRadFun(oqn_l, iDtypeB)
              ! p = 1
              ab1cofVec(lmp + 1, iband, idirC) = ImagUnit**oqn_l * &
                             & dot_product( conjg( z1Analytical(:nv(1, ikpq), idirC, iDatomA, iband) ), aKpq(:nv(1, ikpq), lm, iDatomB) )
              ! p = 2
              ab1cofVec(lmp + 2, iband, idirC) = ImagUnit**oqn_l * &
                             & dot_product( conjg( z1Analytical(:nv(1, ikpq), idirC, iDatomA, iband) ), bKpq(:nv(1, ikpq), lm, iDatomB) )
              do iradf = 3, pMaxLocal
                ! p = 1
                ab1cofVec(lmp + 1, iband, idirC) = ab1cofVec(lmp + 1, iband, idirC) &
                  & + ImagUnit**oqn_l * dot_product( conjg( z1Analytical(nv(1, ikpq) + 1:nv(1, ikpq) + atoms%nlotot, idirC, iDatomA, iband) ), &
                                              & bascof_loKpq(1, mqn_m, :atoms%nlotot, iloTable(iradf, oqn_l, iDtypeB), iDatomB) )
                ! p = 2
                ab1cofVec(lmp + 2, iband, idirC) = ab1cofVec(lmp + 2, iband, idirC) &
                  & + ImagUnit**oqn_l * dot_product(conjg(z1Analytical(nv(1, ikpq) + 1:nv(1, ikpq) + atoms%nlotot, idirC, iDatomA, iband)), &
                                              & bascof_loKpq(2, mqn_m, :atoms%nlotot, iloTable(iradf, oqn_l, iDtypeB), iDatomB) )
                ! 2 < p < LOs for that l and that atom type
                ab1cofVec(lmp + iradf, iband, idirC) = &
                  & ImagUnit**oqn_l * dot_product( conjg( z1Analytical(nv(1, ikpq) + 1:nv(1, ikpq) + atoms%nlotot, idirC, iDatomA, iband) ),&
                                              & bascof_loKpq(3, mqn_m, :atoms%nlotot, iloTable(iradf, oqn_l, iDtypeB), iDatomB) )
              end do ! iradf
              lm = lm + 1
              lmp = lmp + nRadFun(oqn_l, iDtypeB)
            end do ! mqn_m
          end do ! oqn_l
        end do ! iband
      end do ! idirC
      end if

      ab0cof(:, :) = cmplx(0., 0.)
      do iband = 1, nobd(ikpt, 1)
        lmp = 0
        lm  = 0
        do oqn_l = 0, atoms%lmax(iDtypeB)
          do mqn_m = - oqn_l, oqn_l
            pMaxLocal = nRadFun(oqn_l, iDtypeB)
            ! p = 1
            ab0cof(lmp + 1, iband) = ImagUnit**oqn_l * dot_product( conjg(z0(:nv(1, ikpt), iband, ikpt, 1)), a(:nv(1, ikpt), lm, iDatomB) )
            ! p = 2
            ab0cof(lmp + 2, iband) = ImagUnit**oqn_l * dot_product( conjg(z0(:nv(1, ikpt), iband, ikpt, 1)), b(:nv(1, ikpt), lm, iDatomB) )
            ! Add LO contributions
            do iradf = 3, pMaxLocal
              ! p = 1
              ab0cof(lmp + 1, iband) = ab0cof(lmp + 1, iband) + ImagUnit**oqn_l * &
                & dot_product( conjg(z0(nv(1, ikpt) + 1:nv(1, ikpt) + atoms%nlotot, iband, ikpt, 1)), &
                                                     & bascof_lo(1, mqn_m, :atoms%nlotot, iloTable(iradf, oqn_l, iDtypeB), iDatomB) )
              ! p = 2
              ab0cof(lmp + 2, iband) = ab0cof(lmp + 2, iband) + ImagUnit**oqn_l * &
                & dot_product( conjg(z0(nv(1, ikpt) + 1:nv(1, ikpt) + atoms%nlotot, iband, ikpt, 1)), &
                                                     & bascof_lo(2, mqn_m, :atoms%nlotot, iloTable(iradf, oqn_l, iDtypeB), iDatomB) )
              ! 2 < p < LOs for that l and that atom type
              ab0cof(lmp + iradf, iband) = ImagUnit**oqn_l * dot_product( conjg(z0(nv(1, ikpt) + 1:nv(1, ikpt) + atoms%nlotot, iband, ikpt, 1)), &
                                                     & bascof_lo(3, mqn_m, :atoms%nlotot, iloTable(iradf, oqn_l, iDtypeB), iDatomB) )
            end do ! iradf


            ! This is a precalculation of the 1st and 3rd line of A.50. Actually the resulting lmp should have been primed or at
            ! least the p as lm1Pr = lm. But for sake of performance we place it here.
            ! sake of performance
            lm = lm + 1
            lmp = lmp + pMaxLocal
          end do ! mqn_m
        end do !oqn_l
      end do ! iband

      do idirC = 1, 3
        do idirR = 1, 3
          ! Psi1 in Bra
          overlap(:) = cmplx(0., 0.)
          surfIntMT(:) = cmplx(0., 0.)
          !todo this is not correct lmpMax should be lmpT, also in the tests
          ! Patch.CRG
          call CalcSintMT( ikpt, ikpt, idirR, iDatomB, iDtypeB, lmpMax, nobd, hFullNoAbcofBK, overlapNoAbcofBK, ab1cofVec(:, :, idirC), ab0cof,   &
                                                                                                                  & overlap, surfIntMT )
          do iband = 1, nobd(ikpt, 1)
            surfInt(idirR, idirC) = surfInt(idirR, idirC) + 2 * results%w_iks(iband, ikpt, 1) * (surfIntMT(iband) - eig(iband, ikpt, 1) * overlap(iband))
          end do ! iband

          ! Psi1 in ket
          overlap(:) = cmplx(0., 0.)
          surfIntMT(:) = cmplx(0., 0.)
          !todo this is not correct lmpMax should be lmpT, also here
          ! Patch.CRG
          call CalcSintMT( ikpt, ikpt, idirR, iDatomB, iDtypeB, lmpMax, nobd, hFullNoAbcofBK, overlapNoAbcofBK, ab0cof, ab1cofVec(:, :, idirC), &
                                                                                                                  & overlap, surfIntMT )
          do iband = 1, nobd(ikpt, 1)
            surfInt(idirR, idirC) = surfInt(idirR, idirC) + 2 * results%w_iks(iband, ikpt, 1) * (surfIntMT(iband) - eig(iband, ikpt, 1) * overlap(iband))
          end do ! iband
        end do ! idirR
      end do ! idirC

    end do ! ikpt

  end subroutine CalcSFintMTPsi1HepsPsiAndPsiHepsPsi1ExpCoeffVar

  subroutine CalcSFintMTPsi1HepsPsiAndPsiHepsPsi1BasVarikpG( fmpi, noco, nococonv, oneD, atoms, input, sym, usdus, kpts, cell, results, lmpMax, iDtypeA,      &
      & iDatomA, nRadFun, eig, hVarphi, gBas, mapGbas, nv, kveclo, z0, nobd, lmpT, iloTable, varphi1, varphi2, surfInt )

    use m_types
    use m_abcof3

    implicit none

    ! Type parameters
    type(t_mpi),                  intent(in)  :: fmpi
    type(t_nococonv),                intent(in)  :: nococonv
    type(t_noco),                intent(in)  :: noco
    type(t_oneD),                intent(in)  :: oneD
    type(t_atoms),                  intent(in)  :: atoms
    type(t_input),                  intent(in)  :: input
    type(t_sym),                    intent(in)  :: sym
    type(t_usdus),                  intent(in)  :: usdus
    type(t_kpts),                   intent(in)  :: kpts
    type(t_cell),                   intent(in)  :: cell
    type(t_results),                intent(in)  :: results

    ! Scalar parameter
    integer,                        intent(in)  :: lmpMax
    integer,                        intent(in)  :: iDtypeA
    integer,                        intent(in)  :: iDatomA
!todo iDatomB and iDatomA richtig machen
    ! Array parameters
    integer,                        intent(in)  :: nRadFun(0:, :)
    real,                           intent(in)  :: eig(:,:,:)
    complex,                        intent(in)  :: hVarphi(:, :, :, :)
    integer,                        intent(in)  :: gbas(:, :)
    integer,                        intent(in)  :: mapGbas(:, :, :)
    integer,                        intent(in)  :: nv(:, :)
    integer,                        intent(in)  :: kveclo(:,:)
    complex,                       intent(in)  :: z0(:,:,:,:)
    integer,                        intent(in)  :: nobd(:, :)
    integer,                        intent(in)  :: lmpT(:)
    integer,                        intent(in)  :: iloTable(:, 0:, :)
    real,                           intent(in)  :: varphi1(:, :, 0:)
    real,                           intent(in)  :: varphi2(:, :, 0:)
    complex,                        intent(out) :: surfInt(:, :)

    ! Type variables
    type(od_inp)                                :: odi
    type(od_sym)                                :: ods
    type(t_lapw) :: lapw

    ! Scalar variables
    integer                                     :: oqn_l
    integer                                     :: iradf
    integer                                     :: imesh
    integer                                     :: mqn_m
    integer                                     :: lm
    integer                                     :: nmat
    integer                                     :: chanMaxBra
    integer                                     :: chanMaxKet
    integer                                     :: lmaxBra
    integer                                     :: mqn_m2PrBra
    integer                                     :: mqn_m2PrKet
    integer                                     :: lmp
    integer                                     :: ikpt
    integer                                     :: maxNobd
    integer                                     :: idirC
    integer                                     :: idirR
    integer                                     :: iband
    integer                                     :: pmaxLocal
    integer                                     :: iBas
    integer :: nk

    ! Array variables
    integer,           allocatable              :: muKet(:, :)
    integer,           allocatable              :: muBra(:, :)
    complex,           allocatable              :: hFullNoAbcofBK(:, :, :, :)
    complex,           allocatable              :: overlapNoAbcofBK(:, :, :, :)
    real,              allocatable              :: varphiBra1(:, :, :, :)
    real,              allocatable              :: varphiKet1(:, :, :, :)
    real,              allocatable              :: varphiBra2(:, :, :, :)
    real,              allocatable              :: varphiKet2(:, :, :, :)
    integer,           allocatable              :: lambdaKet(:, :)
    integer,           allocatable              :: lambdaBra(:, :)
    integer,           allocatable              :: ngoprI(:)
    complex,           allocatable              :: a(:, :, :)
    complex,           allocatable              :: b(:, :, :)
    complex,           allocatable              :: bascof_lo(:, :, :, :, :)
    complex,           allocatable              :: ab0cofSumVec(:, :, :)
    complex,           allocatable              :: ab0cofVec(:, :, :)
    complex,           allocatable              :: ab0cof(:, :)
    complex,           allocatable              :: surfIntMT(:)
    complex,           allocatable              :: overlap(:)
    complex,           allocatable              :: gBasExt(:, :)
    complex                                     :: kExt(3)

    maxNobd = maxval(nobd(:, :))
    allocate( muKet(-atoms%lmaxd:atoms%lmaxd, -1:1), muBra(-atoms%lmaxd:atoms%lmaxd, -1:1) )
    allocate( hFullNoAbcofBK(lmpMax, lmpMax, 3, atoms%nat), overlapNoAbcofBK(lmpMax, lmpMax, 3, atoms%ntype) )
    allocate( varphiBra1(atoms%jmtd, 2, lmpMax, -1:1), varphiBra2(atoms%jmtd, 2, lmpMax, -1:1), &
            & varphiKet1(atoms%jmtd, 2, lmpMax, -1:1), varphiKet2(atoms%jmtd, 2, lmpMax, -1:1) )
    allocate( lambdaKet(2, 0:atoms%lmaxd), lambdaBra(2, 0:atoms%lmaxd) )
    allocate( ngoprI(atoms%nat) )
    allocate( a( MAXVAL(nv), 0:atoms%lmaxd*(atoms%lmaxd+2), atoms%nat), b(MAXVAL(nv), 0:atoms%lmaxd*(atoms%lmaxd+2), atoms%nat), &
      & bascof_lo(3, -atoms%llod:atoms%llod, 4 * atoms%llod + 2, atoms%nlod, atoms%nat) )
    allocate( ab0cofSumVec(lmpMax, maxNobd, 3), ab0cofVec(lmpMax, maxNobd, 3), ab0cof(lmpMax, maxNobd) )
    allocate( surfIntMT(maxNobd), overlap(maxNobd) )
    allocate( gBasExt(SIZE(z0(:,1,1,1)), 3))

    hFullNoAbcofBK(:, :, :, :) = cmplx(0., 0.)
    overlapNoAbcofBK(:, :, :, :) = cmplx(0., 0.)
    chanMaxBra = 1
    chanMaxKet = 1
    lmaxBra = atoms%lmax(iDtypeA)
    muBra(:, :) = 0
    muKet(:, :) = 0
    mqn_m2PrBra = -1
    mqn_m2PrKet = -1
    varphiBra1 = 0.
    varphiKet1 = 0.
    varphiBra2 = 0.
    varphiKet2 = 0.
    do mqn_m = -atoms%lmax(iDtypeA), atoms%lmax(iDtypeA)
      muKet(mqn_m, -1) = mqn_m
      muBra(mqn_m, -1) = mqn_m
    end do ! mqn_m
    do oqn_l = 0, atoms%lmax(iDtypeA)
      lambdaKet(1, oqn_l) = oqn_l
      lambdaBra(1, oqn_l) = oqn_l
    end do ! oqn_l
    lmp = 0
    do oqn_l = 0, atoms%lmax(iDtypeA)
      do mqn_m = -oqn_l, oqn_l
        do iradf = 1, nRadFun(oqn_l, iDtypeA)
          lmp = lmp + 1
          varphiBra1(atoms%jri(iDtypeA), 1, lmp, -1) = varphi1(atoms%jri(iDtypeA), iradf, oqn_l )
          varphiKet1(atoms%jri(iDtypeA), 1, lmp, -1) = varphi1(atoms%jri(iDtypeA), iradf, oqn_l )
          varphiBra2(atoms%jri(iDtypeA), 1, lmp, -1) = varphi2(atoms%jri(iDtypeA), iradf, oqn_l )
          varphiKet2(atoms%jri(iDtypeA), 1, lmp, -1) = varphi2(atoms%jri(iDtypeA), iradf, oqn_l )
        end do ! iradf
      end do ! mqn_m
    end do ! oqn_l

    call PrepareMTSurfIntDM( atoms, iDtypeA, chanMaxBra, chanMaxKet, lmaxBra, lmpMax, mqn_m2PrBra, mqn_m2PrKet, muKet, muBra,     &
      & lambdaKet, lambdaBra, nRadFun, varphiKet1, varphiKet2, varphiBra1, varphiBra2, hVarphi, atoms%rmt(iDtypeA)**2, hFullNoAbcofBK, overlapNoAbcofBK )

    surfInt(:, :) = cmplx(0., 0.)
    ngoprI(:) = 1
    do ikpt = 1, kpts%nkpt

      nmat = nv(1, ikpt) + atoms%nlotot
      a(:, :, :) = cmplx(0.0, 0.0)
      b(:, :, :) = cmplx(0.0, 0.0)
      bascof_lo(:, :, :, :, :) = cmplx(0.0, 0.0)
      !call abcof3( atoms%lmaxd, atoms%ntype, atoms%nat, sym%nop, MAXVAL(nv), input%jspins, 1, atoms%lmaxd*(atoms%lmaxd+2), SIZE(z0(:,1,1,1)), &
        !& atoms%llod, atoms%nlod, atoms%nlotot, sym%invtab, atoms%ntype, sym%mrot, ngoprI, atoms%taual, atoms%neq, atoms%lmax, &
        !& atoms%rmt, cell%omtil, cell%bmat, cell%bbmat, kpts%bk(:, ikpt), gbas(1, mapGbas(:nv(1, ikpt), ikpt, 1)), &
        !& gbas(2, mapGbas(:nv(1, ikpt), ikpt, 1)), gbas(3, mapGbas(:nv(1, ikpt), ikpt, 1)), nv(:, ikpt), nmat, &
        !& usdus%us, usdus%dus, usdus%uds, usdus%duds, usdus%ddn, sym%invsat, sym%invsatnr, usdus%ulos, usdus%uulon, usdus%dulon, &
        !& usdus%dulos, atoms%llo, atoms%nlo, atoms%l_dulo, atoms%lapw_l, kveclo(:,ikpt), odi, ods, a, b, bascof_lo )
        nk=fmpi%k_list(ikpt)
        CALL lapw%init(input, noco, nococonv, kpts, atoms, sym, nk, cell, .FALSE., fmpi)
        CALL abcof3(input, atoms, sym, 1, cell, kpts%bk(:, ikpt), lapw, &
                            usdus, oneD, 1, lapw%dim_nvd(), a, b, bascof_lo)
      ab0cof(:, :) = cmplx(0., 0.)
      do iband = 1, nobd(ikpt, 1)
        lmp = 0
        lm  = 0
        do oqn_l = 0, atoms%lmax(iDtypeA)
          do mqn_m = - oqn_l, oqn_l
            pMaxLocal = nRadFun(oqn_l, iDtypeA)
            ! p = 1
            ab0cof(lmp + 1, iband) = ImagUnit**oqn_l * dot_product( conjg(z0(:nv(1, ikpt), iband, ikpt, 1)), a(:nv(1, ikpt), lm, iDatomA) )
            ! p = 2
            ab0cof(lmp + 2, iband) = ImagUnit**oqn_l * dot_product( conjg(z0(:nv(1, ikpt), iband, ikpt, 1)), b(:nv(1, ikpt), lm, iDatomA) )
            ! Add LO contributions
            do iradf = 3, pMaxLocal
              ! p = 1
              ab0cof(lmp + 1, iband) = ab0cof(lmp + 1, iband) + ImagUnit**oqn_l * &
                & dot_product( conjg(z0(nv(1, ikpt) + 1:nv(1, ikpt) + atoms%nlotot, iband, ikpt, 1)), &
                                                     & bascof_lo(1, mqn_m, :atoms%nlotot, iloTable(iradf, oqn_l, iDtypeA), iDatomA) )
              ! p = 2
              ab0cof(lmp + 2, iband) = ab0cof(lmp + 2, iband) + ImagUnit**oqn_l * &
                & dot_product( conjg(z0(nv(1, ikpt) + 1:nv(1, ikpt) + atoms%nlotot, iband, ikpt, 1)), &
                                                     & bascof_lo(2, mqn_m, :atoms%nlotot, iloTable(iradf, oqn_l, iDtypeA), iDatomA) )
              ! 2 < p < LOs for that l and that atom type
              ab0cof(lmp + iradf, iband) = ImagUnit**oqn_l * dot_product( conjg(z0(nv(1, ikpt) + 1:nv(1, ikpt) + atoms%nlotot, iband, ikpt, 1)), &
                                                     & bascof_lo(3, mqn_m, :atoms%nlotot, iloTable(iradf, oqn_l, iDtypeA), iDatomA) )
            end do ! iradf


            ! This is a precalculation of the 1st and 3rd line of A.50. Actually the resulting lmp should have been primed or at
            ! least the p as lm1Pr = lm. But for sake of performance we place it here.
            ! sake of performance
            lm = lm + 1
            lmp = lmp + pMaxLocal
          end do ! mqn_m
        end do !oqn_l
      end do ! iband

      kExt(1:3) = matmul(cell%bmat, kpts%bk(1:3, ikpt))

      gbasExt(:, :) = 0.
      ! We assign the G-basis vectors once so later we can just loop without jumps over iBas.
      do iBas = 1, nv(1, ikpt)
        gbasExt(iBas, 1:3) = matmul( cell%bmat(1:3, 1:3), gbas(1:3, mapGbas(iBas, ikpt, 1)))
      end do ! iBas
      do iBas = nv(1, ikpt) + 1, nv(1, ikpt) + atoms%nlotot
        gBasExt(iBas, 1:3) = matmul( cell%bmat(1:3, 1:3), gbas(1:3, mapGbas(kveclo(iBas, ikpt), ikpt, 1)) )
      end do ! iBas

      ! Calculate vector-like large matching coefficients at k-point k with the unperturbed expansion coefficients z.
      ! This vector has to be fully available so we choose idirR
      ab0cofVec(:, :, :) = cmplx(0.0, 0.0)
      ab0cofSumVec(:, :, :) = cmplx(0.0, 0.0)
      do idirR = 1, 3
        do iband = 1, nobd(ikpt, 1)
          lmp = 0
          lm  = 0
          do oqn_l = 0, atoms%lmax(iDtypeA)
            do mqn_m = -oqn_l, oqn_l
              pMaxLocal = nRadFun(oqn_l, iDtypeA)
              ! p = 1
              ab0cofVec(lmp + 1, iband, idirR) = ImagUnit**oqn_l * dot_product( conjg( z0(:nv(1, ikpt), iband, ikpt, 1)), &
                                                                  & gbasExt(:nv(1, ikpt), idirR) * a(:nv(1, ikpt), lm, iDatomA) )
              ! p = 2
              ab0cofVec(lmp + 2, iband, idirR) = ImagUnit**oqn_l * dot_product( conjg( z0(:nv(1, ikpt), iband, ikpt, 1)), &
                                                                  & gbasExt(:nv(1, ikpt), idirR) * b(:nv(1, ikpt), lm, iDatomA) )
              do iradf = 3, pMaxLocal
                ! p = 1
                ab0cofVec(lmp + 1, iband, idirR) = ab0cofVec(lmp + 1, iband, idirR) + ImagUnit**oqn_l * &
                  & dot_product( conjg( z0(nv(1, ikpt) + 1:nv(1, ikpt) + atoms%nlotot, iband, ikpt, 1) ), &
                                                   &   gbasExt(nv(1, ikpt) + 1:nv(1, ikpt) + atoms%nlotot, idirR) &
                                                   & * bascof_lo(1, mqn_m, :atoms%nlotot, iloTable(iradf, oqn_l, iDtypeA), iDatomA) )
                ! p = 2
                ab0cofVec(lmp + 2, iband, idirR) = ab0cofVec(lmp + 2, iband, idirR) + ImagUnit**oqn_l * &
                  & dot_product( conjg( z0(nv(1, ikpt) + 1:nv(1, ikpt) + atoms%nlotot, iband, ikpt, 1) ), &
                                                   &   gbasExt(nv(1, ikpt) + 1:nv(1, ikpt) + atoms%nlotot, idirR) &
                                                   & * bascof_lo(2, mqn_m, :atoms%nlotot, iloTable(iradf, oqn_l, iDtypeA), iDatomA) )
                ! 2 < p < LOs for that l and that atom type
                ab0cofVec(lmp + iradf, iband, idirR) = ImagUnit**oqn_l * dot_product( conjg( z0(nv(1, ikpt) + 1:nv(1, ikpt) + atoms%nlotot, iband, ikpt, 1) ),&
                                                   &   gbasExt(nv(1, ikpt) + 1:nv(1, ikpt) + atoms%nlotot, idirR) &
                                                   & * bascof_lo(3, mqn_m, :atoms%nlotot, iloTable(iradf, oqn_l, iDtypeA), iDatomA) )
              end do ! iradf

              ! This help array has to be calulated in a seperate loop here because earlier ab0cofVec would not be available and the
              ! next command needs ab0cofSumVec already.
              do iradf = 1, pMaxLocal
                ab0cofSumVec(lmp + iradf, iband, idirR) = ImagUnit * (ab0cof(lmp + iradf, iband) * kExt(idirR) + ab0cofVec(lmp + iradf, iband, idirR))
              end do ! iradf

              ! Precalculation of the 2nd and the 4th line in A.51. Due to performance, the lmp in the resulting quantity are not
              ! primed
              lm = lm + 1
              lmp = lmp + pMaxLocal
            end do ! mqn_m
          end do ! oqn_l
        end do ! iband
      end do ! idirR

      do idirC = 1, 3
        do idirR = 1, 3
          ! Psi1 in Bra
          overlap(:) = cmplx(0., 0.)
          surfIntMT(:) = cmplx(0., 0.)
          call CalcSintMT( ikpt, ikpt, idirR, iDatomA, iDtypeA, lmpMax, nobd, hFullNoAbcofBK, overlapNoAbcofBK, ab0cofSumVec(:, :, idirC), ab0cof,   &
                                                                                                                  & overlap, surfIntMT )
          do iband = 1, nobd(ikpt, 1)
            surfInt(idirR, idirC) = surfInt(idirR, idirC) + 2 * results%w_iks(iband, ikpt, 1) * (surfIntMT(iband) - eig(iband, ikpt, 1) * overlap(iband))
          end do ! iband

          ! Psi1 in ket
          overlap(:) = cmplx(0., 0.)
          surfIntMT(:) = cmplx(0., 0.)
          call CalcSintMT( ikpt, ikpt, idirR, iDatomA, iDtypeA, lmpMax, nobd, hFullNoAbcofBK, overlapNoAbcofBK, ab0cof, ab0cofSumVec(:, :, idirC), &
                                                                                                                  & overlap, surfIntMT )
          do iband = 1, nobd(ikpt, 1)
            surfInt(idirR, idirC) = surfInt(idirR, idirC) + 2 * results%w_iks(iband, ikpt, 1) * (surfIntMT(iband) - eig(iband, ikpt, 1) * overlap(iband))
          end do ! iband
        end do ! idirR
      end do ! idirC

    end do ! ikpt
  end subroutine CalcSFintMTPsi1HepsPsiAndPsiHepsPsi1BasVarikpG

  subroutine CalcSurfintMTPsigrVeff0Psi( fmpi, noco, nococonv, oneD, atoms, input, sym, cell, kpts, usdus, results, lmpMax, iDtypeB, iDatomB, nobd, gbas, &
      & mapGbas, nRadFun, nv, z0, kveclo, varphi1, varphi2, grVeff0MT, iloTable, surfInt )

    use m_types_oneD
    use m_types
    use m_abcof3

    implicit none

    ! Type parameter
    type(t_mpi),                  intent(in)  :: fmpi
    type(t_nococonv),                intent(in)  :: nococonv
    type(t_noco),                intent(in)  :: noco
    type(t_oneD),                intent(in)  :: oneD
    type(t_atoms),                  intent(in)  :: atoms
    type(t_input),                  intent(in)  :: input
    type(t_sym),                    intent(in)  :: sym
    type(t_cell),                   intent(in)  :: cell
    type(t_kpts),                   intent(in)  :: kpts
    type(t_usdus),                  intent(in)  :: usdus
    type(t_results),                intent(in)  :: results


    ! Scalar parameter
    integer,                        intent(in)  :: lmpMax
    integer,                        intent(in)  :: iDtypeB
    integer,                        intent(in)  :: iDatomB

    ! Array parameter
    integer,                        intent(in)  :: nobd(:, :)
    integer,                        intent(in)  :: gbas(:, :)
    integer,                        intent(in)  :: mapGbas(:, :, :)
    integer,                        intent(in)  :: nRadFun(0:, :)
    integer,                        intent(in)  :: nv(:, :)
    complex,                       intent(in)  :: z0(:,:,:,:)
    integer,                        intent(in)  :: kveclo(:,:)
    real,                           intent(in)  :: varphi1(:, :, 0:)
    real,                           intent(in)  :: varphi2(:, :, 0:)
    integer,                        intent(in)  :: iloTable(:, 0:, :)
    complex,                        intent(in)  :: grVeff0MT(:, :, :, :)
    complex,                        intent(out) :: surfInt(:, :)

    ! Type variables
    type(od_inp)                                :: odi
    type(od_sym)                                :: ods
    type(t_lapw) :: lapw

    ! Scalar variables
    integer                                     :: idirR
    integer                                     :: iband
    integer                                     :: idirC
    integer                                     :: nmat
    integer                                     :: lmp
    integer                                     :: lm
    integer                                     :: lmaxBra
    integer                                     :: chanMaxBra
    integer                                     :: chanMaxKet
    integer                                     :: mqn_m2PrBra
    integer                                     :: mqn_m2PrKet
    integer                                     :: iradf
    integer                                     :: oqn_l
    integer                                     :: mqn_m
    integer                                     :: ikpt
    integer                                    :: pMaxLocal
    integer                                     :: maxNobd
    INTEGER :: nk

    ! Array variables
    complex,           allocatable              :: ab0cofKet(:, :)
    complex,           allocatable              :: hFullNoAbcofBK(:, :, :, :)
    complex,           allocatable              :: overlapNoAbcofBK(:, :, :, :)
    complex,           allocatable              :: overlap(:)
    complex,           allocatable              :: surfIntMT(:)
    complex,           allocatable              :: a(:, :, :)
    complex,           allocatable              :: b(:, :, :)
    complex,           allocatable              :: bascof_lo(:, :, :, :, :)
    integer,           allocatable              :: ngoprI(:)
    integer,           allocatable              :: lambdaKet(:, :)
    integer,           allocatable              :: lambdaBra(:, :)
    integer,           allocatable              :: muKet(:, :)
    integer,           allocatable              :: muBra(:, :)
    real,              allocatable              :: varphiBra1(:, :, :, :)
    real,              allocatable              :: varphiKet1(:, :, :, :)
    real,              allocatable              :: varphiBra2(:, :, :, :)
    real,              allocatable              :: varphiKet2(:, :, :, :)
    complex,           allocatable              :: grVeff0Varphi(:, :, :, :)

    maxNobd = maxval( nobd(:, :) )
    allocate( ab0cofKet(lmpMax, maxNobd) )
    allocate( hFullNoAbcofBK(lmpMax, lmpMax, 3, atoms%nat), overlapNoAbcofBK(lmpMax, lmpMax, 3, atoms%ntype) )
    allocate( surfIntMT(maxNobd), overlap(maxNobd) )
    allocate( a( MAXVAL(nv), 0:atoms%lmaxd*(atoms%lmaxd+2), atoms%nat), b(MAXVAL(nv), 0:atoms%lmaxd*(atoms%lmaxd+2), atoms%nat), &
      & bascof_lo(3, -atoms%llod:atoms%llod, 4 * atoms%llod + 2, atoms%nlod, atoms%nat) )
    allocate( ngoprI(atoms%nat) )
    allocate( lambdaKet(2, 0:atoms%lmaxd), lambdaBra(2, 0:atoms%lmaxd) )
    allocate( muKet(-atoms%lmaxd:atoms%lmaxd, -1:1), muBra(-atoms%lmaxd:atoms%lmaxd, -1:1) )
    allocate( varphiBra1(atoms%jmtd, 2, lmpMax, -1:1), varphiBra2(atoms%jmtd, 2, lmpMax, -1:1), &
            & varphiKet1(atoms%jmtd, 2, lmpMax, -1:1), varphiKet2(atoms%jmtd, 2, lmpMax, -1:1) )
    allocate( grVeff0Varphi(2, atoms%jmtd, 0:(atoms%lmaxd + 3)**2 - 1, lmpMax) )
    ngoprI(:) = 1
    surfInt(:, :) = cmplx(0., 0.)

    do idirR = 1, 3

      ! Action of normal non-spherical potential onto MT basis functions without basis matching coefficients
      grVeff0Varphi = cmplx(0., 0.)
      ! todo grRho0MT is grVeff0MT change that
      call CalcFnsphVarphi( atoms, iDtypeB, 0, nRadFun, varphi1, varphi2, grVeff0MT(:, :, idirR, iDatomB), grVeff0Varphi )

      hFullNoAbcofBK(:, :, :, :) = cmplx(0., 0.)
      overlapNoAbcofBK(:, :, :, :) = cmplx(0., 0.)
      chanMaxBra = 1
      chanMaxKet = 1
      lmaxBra = atoms%lmax(iDtypeB)
      muBra(:, :) = 0
      muKet(:, :) = 0
      mqn_m2PrBra = -1
      mqn_m2PrKet = -1
      varphiBra1 = 0.
      varphiKet1 = 0.
      varphiBra2 = 0.
      varphiKet2 = 0.
      do mqn_m = -atoms%lmax(iDtypeB), atoms%lmax(iDtypeB)
        muKet(mqn_m, -1) = mqn_m
        muBra(mqn_m, -1) = mqn_m
      end do ! mqn_m
      do oqn_l = 0, atoms%lmax(iDtypeB)
        lambdaKet(1, oqn_l) = oqn_l
        lambdaBra(1, oqn_l) = oqn_l
      end do ! oqn_l
      lmp = 0
      do oqn_l = 0, atoms%lmax(iDtypeB)
        do mqn_m = -oqn_l, oqn_l
          do iradf = 1, nRadFun(oqn_l, iDtypeB)
            lmp = lmp + 1
            varphiBra1(atoms%jri(iDtypeB), 1, lmp, -1) = varphi1(atoms%jri(iDtypeB), iradf, oqn_l )
            varphiKet1(atoms%jri(iDtypeB), 1, lmp, -1) = varphi1(atoms%jri(iDtypeB), iradf, oqn_l )
            varphiBra2(atoms%jri(iDtypeB), 1, lmp, -1) = varphi2(atoms%jri(iDtypeB), iradf, oqn_l )
            varphiKet2(atoms%jri(iDtypeB), 1, lmp, -1) = varphi2(atoms%jri(iDtypeB), iradf, oqn_l )
          end do ! iradf
        end do ! mqn_m
      end do ! oqn_l

      call PrepareMTSurfIntDM( atoms, iDtypeB, chanMaxBra, chanMaxKet, lmaxBra, lmpMax, mqn_m2PrBra, mqn_m2PrKet, muKet, muBra,     &
        & lambdaKet, lambdaBra, nRadFun, varphiKet1, varphiKet2, varphiBra1, varphiBra2, grVeff0Varphi, atoms%rmt(iDtypeB)**2, hFullNoAbcofBK, overlapNoAbcofBK )

      do ikpt = 1, kpts%nkpt

        ! Calculate the basis matching coefficients at k + q = k' and the matching coefficients at k.
        nmat = nv(1, ikpt) + atoms%nlotot
        a(:, :, :) = cmplx(0.0, 0.0)
        b(:, :, :) = cmplx(0.0, 0.0)
        bascof_lo(:, :, :, :, :) = cmplx(0.0, 0.0)
        !call abcof3( atoms%lmaxd, atoms%ntype, atoms%nat, sym%nop, MAXVAL(nv), input%jspins, 1, atoms%lmaxd*(atoms%lmaxd+2), SIZE(z0(:,1,1,1)), &
        !  & atoms%llod, atoms%nlod, atoms%nlotot, sym%invtab, atoms%ntype, sym%mrot, ngoprI, atoms%taual, atoms%neq, atoms%lmax, &
        !  & atoms%rmt, cell%omtil, cell%bmat, cell%bbmat, kpts%bk(:, ikpt), gbas(1, mapGbas(:nv(1, ikpt), ikpt, 1)), &
        !  & gbas(2, mapGbas(:nv(1, ikpt), ikpt, 1)), gbas(3, mapGbas(:nv(1, ikpt), ikpt, 1)), nv(:, ikpt), nmat, &
        !  & usdus%us, usdus%dus, usdus%uds, usdus%duds, usdus%ddn, sym%invsat, sym%invsatnr, usdus%ulos, usdus%uulon, usdus%dulon, &
        !  & usdus%dulos, atoms%llo, atoms%nlo, atoms%l_dulo, atoms%lapw_l, kveclo(:,ikpt), odi, ods, a, b, bascof_lo )
        !!!!!DO jsp = 1, MERGE(1,fi%input%jspins,fi%noco%l_noco) TODO!!!!!
         !!!!!!  k_loop:DO nk_i = 1,size(fmpi%k_list)
        !nk=fmpi%k_list(nk_i)
        ! Set up lapw list
        !CALL lapw%init(input, noco, nococonv, kpts, atoms, sym, nk, cell, .FALSE., fmpi)
        nk=fmpi%k_list(ikpt)
        CALL lapw%init(input, noco, nococonv, kpts, atoms, sym, nk, cell, .FALSE., fmpi)
        CALL abcof3(input, atoms, sym, 1, cell, kpts%bk(:, ikpt), lapw, &
                            usdus, oneD, 1, lapw%dim_nvd(), a, b, bascof_lo)
        ab0cofKet(:, :) = cmplx(0., 0.)
        do iband = 1, nobd(ikpt, 1)
          lmp = 0
          lm  = 0
          do oqn_l = 0, atoms%lmax(iDtypeB)
            do mqn_m = - oqn_l, oqn_l
              pMaxLocal = nRadFun(oqn_l, iDtypeB)
              ! p = 1
              ab0cofKet(lmp + 1, iband) = ImagUnit**oqn_l * dot_product( conjg(z0(:nv(1, ikpt), iband, ikpt, 1)), a(:nv(1, ikpt), lm, iDatomB) )
              ! p = 2
              ab0cofKet(lmp + 2, iband) = ImagUnit**oqn_l * dot_product( conjg(z0(:nv(1, ikpt), iband, ikpt, 1)), b(:nv(1, ikpt), lm, iDatomB) )
              ! Add LO contributions
              do iradf = 3, pMaxLocal
                ! p = 1
                ab0cofKet(lmp + 1, iband) = ab0cofKet(lmp + 1, iband) + ImagUnit**oqn_l * &
                  & dot_product( conjg(z0(nv(1, ikpt) + 1:nv(1, ikpt) + atoms%nlotot, iband, ikpt, 1)), &
                                                       & bascof_lo(1, mqn_m, :atoms%nlotot, iloTable(iradf, oqn_l, iDtypeB), iDatomB) )
                ! p = 2
                ab0cofKet(lmp + 2, iband) = ab0cofKet(lmp + 2, iband) + ImagUnit**oqn_l * &
                  & dot_product( conjg(z0(nv(1, ikpt) + 1:nv(1, ikpt) + atoms%nlotot, iband, ikpt, 1)), &
                                                       & bascof_lo(2, mqn_m, :atoms%nlotot, iloTable(iradf, oqn_l, iDtypeB), iDatomB) )
                ! 2 < p < LOs for that l and that atom type
                ab0cofKet(lmp + iradf, iband) = ImagUnit**oqn_l * dot_product( conjg(z0(nv(1, ikpt) + 1:nv(1, ikpt) + atoms%nlotot, iband, ikpt, 1)), &
                                                       & bascof_lo(3, mqn_m, :atoms%nlotot, iloTable(iradf, oqn_l, iDtypeB), iDatomB) )
              end do ! iradf


              ! This is a precalculation of the 1st and 3rd line of A.50. Actually the resulting lmp should have been primed or at
              ! least the p as lm1Pr = lm. But for sake of performance we place it here.
              ! sake of performance
              lm = lm + 1
              lmp = lmp + pMaxLocal
            end do ! mqn_m
          end do !oqn_l
        end do ! iband

        do idirC = 1, 3

          overlap(:) = cmplx(0., 0.)
          surfIntMT(:) = cmplx(0., 0.)
          call CalcSintMT( ikpt, ikpt, idirC, iDatomB, iDtypeB, lmpMax, nobd, hFullNoAbcofBK, overlapNoAbcofBK, ab0cofKet, ab0cofKet,   &
                                                                                                          & overlap, surfIntMT )
          do iband = 1, nobd(ikpt, 1)
            surfInt(idirR, idirC) = surfInt(idirR, idirC) + results%w_iks(iband, ikpt, 1) * surfIntMT(iband)
          end do ! iband
        end do ! idirC
      end do ! ikpt
    end do ! idirR

  end subroutine CalcSurfintMTPsigrVeff0Psi

  subroutine CalcHnGrV0Varphi( atoms, sym, lathar, itype, iatom, lmpMax, El, varphi1, varphi2, nRadFun, vEff0MtSpH, vEff0MtLh, clnu_atom, &
      & nmem_atom, mlh_atom, grVarphiCh1, grVarphiCh2, grVarphiChLout, grVarphiChMout, hVarphi, vEff0NsphGrVarphi, r2grVeff0SphVarphi, r2grVeff0Varphi )

    use m_types, only : t_atoms, t_sym, t_sphhar

!    use mod_juPhonUtils, only : fopen, fclose
    implicit none

    ! Type parameter
    type(t_atoms),        intent(in)  :: atoms
    type(t_sym),        intent(in)  :: sym
    type(t_sphhar),       intent(in)  :: lathar

    ! Scalar parameter
    integer,              intent(in)  :: itype
    integer,              intent(in)  :: iatom
    integer,              intent(in)  :: lmpMax

    ! Array parameter
    real,                 intent(in)  :: El(:, 0:, :, :)
    real,                 intent(in)  :: varphi1(:,:,0:)
    real,                 intent(in)  :: varphi2(:,:,0:)
    integer,              intent(in)  :: nRadFun(0:, :)
    complex,              intent(in)  :: vEff0MtSpH(:, :)
    real,                 intent(in)  :: vEff0MtLh(:, 0:, :)
    complex,              intent(in)  :: clnu_atom(:, 0:, :)
    integer,              intent(in)  :: nmem_atom(0:, :)
    integer,              intent(in)  :: mlh_atom(:, 0:, :)
    real,                 intent(in)  :: grVarphiCh1(:, :, :, -1:)
    real,                 intent(in)  :: grVarphiCh2(:, :, :, -1:)
    integer,              intent(in)  :: grVarphiChLout(:, 0:)
    integer,              intent(in)  :: grVarphiChMout(-atoms%lmaxd:, -1:)
    complex,              intent(out) :: hVarphi(:, :, 0:, :)
    complex,              intent(out) :: vEff0NsphGrVarphi(:, :, :, :, -1:)
    complex,              intent(out) :: r2grVeff0SphVarphi(:, :, :, :, :)
    complex,              intent(out) :: r2grVeff0Varphi(:, :, :, :, :)



    ! Scalar variables
    integer                           :: lmp
    integer                           :: oqn_l
    integer                           :: mqn_m
    integer                           :: oqn_l1Pr
    integer                           :: mqn_m1Pr
    integer                           :: imesh
    integer                           :: iradf
    real                              :: rInv
    integer                           :: idir
    integer                           :: lm
    integer                           :: lm_pre
    integer                           :: lm1Pr
    integer                           :: lm1Pr_pre
    integer                           :: irel
    integer                           :: mqn_m2PrC
    integer                           :: ilh

    ! Array variables
    real,           allocatable       :: hSph(:, :, :)
    real,           allocatable       :: r2Veff0SphLh(:, :, :)
    real,           allocatable       :: r2Veff0Lh(:, :, :)
    complex,        allocatable       :: r2grVeff0SphSh( :, :, :, : )
    complex,        allocatable       :: r2grVeff0Sh( :, :, :, : )
    real,           allocatable       :: recMesh(:)
    complex,        allocatable       :: grVeff0ShRobust(:, :, :)
    complex,        allocatable       :: vnsphEff0Varphi(:, :, :, :)

!    real,           allocatable       :: vxc0mt(:, :, :, :)
!
    allocate( hSph(2, atoms%jmtd, lmpMax), vnsphEff0Varphi(2, atoms%jmtd, 0:(atoms%lmaxd + 3)**2 - 1, lmpMax) )

    hSph(:, :, :) = 0.
    vnsphEff0Varphi = cmplx(0., 0.)

    ! Spherical part only contributes until lmax and not lmax + 2 as we have only the ket here expanded until lmax and in the overlap we
    ! cut there at lmax.
    ! It should be the best compromise to close the loops here over oqn_l, mqn_m and p otherwise we have redundant code or if-clauses in
    ! the inner loops. With this solution we only count the loops of oqn_l, mqn_m and p twice.
    ! If we would have run the p loop from 1, this had led to if-clauses in the inner loops.
    lmp = 0
    do oqn_l = 0, atoms%lmax(itype)
      do mqn_m = -oqn_l, oqn_l
        lmp = lmp + 1
        do imesh = 1, atoms%jri(itype)
          hSph(1, imesh, lmp) = El(1, oqn_l, itype, 1) * varphi1(imesh, 1, oqn_l)
          hSph(2, imesh, lmp) = El(1, oqn_l, itype, 1) * varphi2(imesh, 1, oqn_l)
        end do ! imesh
        lmp = lmp + 1
        do imesh = 1, atoms%jri(itype)
          hSph(1, imesh, lmp) = varphi1(imesh, 1, oqn_l) + El(1, oqn_l, itype, 1) * varphi1(imesh, 2, oqn_l)
          hSph(2, imesh, lmp) = varphi2(imesh, 1, oqn_l) + El(1, oqn_l, itype, 1) * varphi2(imesh, 2, oqn_l)
        end do ! imesh
        do iradf = 3, nRadFun(oqn_l, itype)
          lmp = lmp + 1
          do imesh = 1, atoms%jri(itype)
            hSph(1, imesh, lmp) = El(iradf - 1, oqn_l, itype, 1) * varphi1(imesh, iradf, oqn_l)
            hSph(2, imesh, lmp) = El(iradf - 1, oqn_l, itype, 1) * varphi2(imesh, iradf, oqn_l)
          end do ! imesh
        end do ! p
      end do ! mqn_m
    end do ! oqn_l

    ! Action of normal non-spherical potential onto MT basis functions without basis matching coefficients
!    call CalcFnsphVarphi( atoms, itype, 0, nRadFun, varphi1, varphi2, vEff0MtSpH, vnsphEff0Varphi )
    call CalcFnsphVarphi( atoms, itype, 1, nRadFun, varphi1, varphi2, vEff0MtSpH, vnsphEff0Varphi )

    ! Action of complete Hamiltonian onto MT basis functions without basis matching coefficients
    lmp = 0
    do oqn_l = 0, atoms%lmax(itype)
      lm_pre = oqn_l * (oqn_l + 1)
      do mqn_m = -oqn_l, oqn_l
        lm = lm_pre + mqn_m
        do iradf = 1, nRadFun(oqn_l, itype)
          lmp = lmp + 1
          ! Add action of spherical Hamiltonian to the diagonal of lmp and l'm'p'.
          do imesh = 1, atoms%jri(itype)
            do irel = 1, 2
              hVarphi(irel, imesh, lm, lmp) = hSph(irel, imesh, lmp)
            end do ! irel
          end do ! imesh
          ! We store the bra index until lmax + 2 because the double gradient of the radial solution is expanded maximally until
          ! lmax + 2
          !do oqn_l1Pr = 0, atoms%lmax(itype) + 2
          do oqn_l1Pr = 0, atoms%lmax(itype) !+ 2
            lm1Pr_pre = oqn_l1Pr * (oqn_l1Pr + 1)
            do mqn_m1Pr = -oqn_l1Pr, oqn_l1Pr
              lm1Pr = lm1Pr_pre + mqn_m1Pr
              do imesh = 1, atoms%jri(itype)
                do irel = 1, 2
                  ! The factors i^l actually assigned to the matching coefficients in theoretical equations are excluded from abcof
                  ! and abcof3. Therefore, we multiply them here, so that, we do not need to care for them later.
                  hVarphi(irel, imesh, lm1Pr, lmp) = hVarphi(irel, imesh, lm1Pr, lmp) + vnsphEff0Varphi(irel, imesh, lm1Pr, lmp)
!                  hVarphi(irel, imesh, lm1Pr, lmp) = vnsphEff0Varphi(irel, imesh, lm1Pr, lmp)
                end do ! irel
              end do ! imesh
            end do ! mqn_m
          end do ! oqn_l
        end do ! iradf
      end do ! mqn_m1Pr
    end do ! oqn_l1Pr
    deallocate(hSph)

    ! Construction of the special non-spherical potential required for the calculation of <gradVarphi|H|gradVarphi>

    ! Remove numerically critical singularity at the core
    ! !todo test this!
    ! allocate( recMesh(atoms%jmtd) )
    ! recMesh = 0.
    ! do imesh = 1, atoms%jri(itype)
    !   rInv = 1. / atoms%rmsh(imesh, itype)
    !   recMesh(imesh) =  rInv / atoms%rmsh(imesh, itype)
    !   r2VeffSpLh(imesh, 0) = r2VeffSpLh(imesh, 0) - atoms%zatom(itype) * rInv
    ! end do ! imesh
    !!call CalcGrFLhNatAt(atoms, lathar, itype, iatom, clnu_atom, nmem_atom, mlh_atom, vEff0MtLh(:, :, itype), grVeff0ShRobust)

    ! ! The analytical gradient 1 / r of the l = 0 channel is expanded in spherical harmonics and added to the robust numerically derived
    ! ! gradient part in the l = 1 channel
    ! ! todo c_mi is the transposed vesion of c_im!
    ! do idir = -1, 1
    !   do mqn_m = -1, 1
    !     ! l * (l + 1) = 2
    !     lm = 2 + mqn_m
    !     do imesh = 1, atoms%jri(itype)
    !       grVeff0ShRobust(imesh, lm, idir) = grVeff0ShRobust(imesh, lm, idir) - recMesh(imesh) * c_mi(mqn_m + 2, idir + 2)
    !     end do ! imesh
    !   end do ! mqn_m
    ! end do ! idir
    ! deallocate(recMesh)

    allocate(r2Veff0SphLh(atoms%jmtd, 0:lathar%nlhd, atoms%ntype))
    r2Veff0SphLh(:, :, :) = cmplx(0.0, 0.0)

!    allocate( vXC0MT(atoms%jmtd, 0:lathar%nlhd, atoms%ntype, 1) )
!    vxc0mt(:, :, :, :) = cmplx(0., 0.)
!    ! XC potential in the MT read in from FLEUR
!    call fopen(1000, name='v0MTFLEUR_xc', status='old', action='read', form='unformatted')
!    read(1000) vXC0MT(:, :, :, :)
!    call fclose(1000)
!    do ilh = 0, lathar%nlhd
    do imesh = 1, atoms%jri(itype)
      r2Veff0SphLh(imesh, 0, itype) = vEff0Mtlh(imesh, 0, itype) * atoms%rmsh(imesh, itype) * atoms%rmsh(imesh, itype)
!      !r2Veff0SphLh(imesh, ilh, itype) = vEff0Mtlh(imesh, ilh, itype) * atoms%rmsh(imesh, itype) * atoms%rmsh(imesh, itype)
!      r2Veff0SphLh(imesh, ilh, itype) = vxc0mt(imesh, ilh, itype, 1) * atoms%rmsh(imesh, itype) * atoms%rmsh(imesh, itype)
!      r2Veff0SphLh(imesh, ilh, itype) = (vEff0Mtlh(imesh, ilh, itype) - vxc0mt(imesh, ilh, itype, 1)) * atoms%rmsh(imesh, itype) * atoms%rmsh(imesh, itype)
    end do
!    end do

    allocate( r2grVeff0SphSh( atoms%jmtd, (atoms%lmax(itype) + 2)**2, atoms%nat, 3 ) )
    r2grVeff0SphSh = cmplx(0.0, 0.0)

    call calcGrR2FinLH( atoms, sym, lathar, clnu_atom, nmem_atom, mlh_atom, r2Veff0SphLh, r2grVeff0SphSh )

    allocate(r2Veff0Lh(atoms%jmtd, 0:lathar%nlhd, atoms%ntype))
    r2Veff0Lh(:, :, :) = cmplx(0.0, 0.0)
    do ilh = 0, lathar%nlhd
      do imesh = 1, atoms%jri(itype)
        r2Veff0Lh(imesh, ilh, itype) = vEff0Mtlh(imesh, ilh, itype) * atoms%rmsh(imesh, itype) * atoms%rmsh(imesh, itype)
      end do ! imesh
    end do ! ilh

    allocate( r2grVeff0Sh( atoms%jmtd, (atoms%lmax(itype) + 2)**2, atoms%nat, 3 ) )
    r2grVeff0Sh = cmplx(0.0, 0.0)
    call calcGrR2FinLH( atoms, sym, lathar, clnu_atom, nmem_atom, mlh_atom, r2Veff0Lh, r2grVeff0Sh )


  !  do idir = 1, 3
  !    do oqn_l = 0, atoms%lmax(itype) + 1
  !      do mqn_m = -oqn_l, oqn_l
  !        lm = oqn_l * (oqn_l + 1) + 1 + mqn_m
  !        do imesh = 1, atoms%jmtd
  !          write(4000, '(4i8,2f15.8)') idir, oqn_l, mqn_m, imesh, r2grVeff0SphSh(imesh, lm, 1, idir)
  !        end do ! imesh
  !      end do ! mqn_m
  !    end do ! oqn_l
  !  end do ! idir


    !todo before commnit: is it critical that lmPr is running until lmax + 2
    ! NOTE: The main difference comes from the Coulomb potential, because here we have a difference in Weinert or numerical, in the
    ! xc potential, we both do a numerical gradient so do not feel a big difference between chain rule to the kernel or direct
    ! numerical gradient. In the Weinert method, for Neon we only have a contribuiton of 1e-6 from the non-spherical contributions of
    ! the potential. The main difference comes from the spherical component. Here, we only have a difference of 1e-4 between using
    ! Weinert or not.
    r2grVeff0SphVarphi(:, :, :, :, :) = cmplx(0., 0.)
    do idir = 1, 3
      call CalcFnsphVarphi( atoms, itype, 0, nRadFun, varphi1, varphi2, r2grVeff0SphSh(:, :, iatom, idir), &
                                                                                            & r2grVeff0SphVarphi(:, :, :, :, idir) )
    end do

    r2grVeff0Varphi(:, :, :, :, :) = cmplx(0., 0.)
    do idir = 1, 3
      call CalcFnsphVarphi( atoms, itype, 0, nRadFun, varphi1, varphi2, r2grVeff0Sh(:, :, iatom, idir), &
                                                                                            & r2grVeff0Varphi(:, :, :, :, idir) )
    end do

    vEff0NsphGrVarphi = cmplx(0., 0.)
    do mqn_m2PrC = -1, 1
      call calcFnsphGrVarphi( atoms, itype, 1, mqn_m2PrC, nRadFun, grVarphiCh1, grVarphiCh2, grVarphiChLout, grVarphiChMout, vEff0MtSpH, &
                                                                                        & vEff0NsphGrVarphi(:, :, :, :, mqn_m2PrC) )
!      call calcFnsphGrVarphi( atoms, itype, 0, mqn_m2PrC, nRadFun, grVarphiCh1, grVarphiCh2, grVarphiChLout, grVarphiChMout, vEff0MtSpH, &
!                                                                                        & vEff0NsphGrVarphi(:, :, :, :, mqn_m2PrC) )
    end do ! mqn_m2Pr

    !todo should be up to lmax + 1
    !todo put in the generation of the spherical potential in here
  end subroutine CalcHnGrV0Varphi

  subroutine CalcSurfIntIRDynMat( atoms, cell, ngdp1, ngdp2, gdp1, gdp2, rho0IRpw, grVext0IR, qpoint, surfInt )

    use m_types
    use m_ylm_old
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

  subroutine readInz1( atoms, ikpt, iqpt, ikpq, nobd, nv, z1nG )

    use m_types, only : t_atoms
    use m_juDFT_stop, only : juDFT_warn

    implicit none

    ! Type parameter
    type(t_atoms),                 intent(in)  :: atoms

    ! Scalar parameter
    integer,                       intent(in)  :: ikpt
    integer,                       intent(in)  :: iqpt
    integer,                       intent(in)  :: ikpq

    ! Array parameter
    integer,                       intent(in)  :: nobd(:, :)
    integer,                       intent(in)  :: nv(:, :)
    complex,                       intent(out) :: z1nG(:, :, :, :)

    ! Scalar variables
    integer                                    :: iDtype
    integer                                    :: iDeqat
    integer                                    :: iDatom
    integer                                    :: iband
    integer                                    :: idir
    integer                                    :: iBas
    logical                                    :: st

    ! Array variables
    character(len=:), allocatable              :: filename
    character(len=15)                          :: filenameTemp


    iDatom = 0
    z1nG(:, :, :, :) = cmplx(0.0, 0.0)
    do iDtype = 1, atoms%ntype
      do iDeqat = 1, atoms%neq(iDtype)
        iDatom = iDatom + 1
        if (ikpt < 10 .and. iqpt >= 10 ) then
          write(filenameTemp, '(a,i1,a,i0,a,i2.2)') 'z1a', iDatom, 'q', iqpt, 'k', ikpt
        else if( ikpt >= 10 .and. iqpt < 10) then
          write(filenameTemp, '(a,i1,a,i2.2,a,i0)') 'z1a', iDatom, 'q', iqpt, 'k', ikpt
        else if( ikpt < 10 .and. iqpt < 10) then
          write(filenameTemp, '(a,i1,a,i2.2,a,i2.2)') 'z1a', iDatom, 'q', iqpt, 'k', ikpt
        else
          write(filenameTemp, '(a,i1,a,i0,a,i0)') 'z1a', iDatom, 'q', iqpt, 'k', ikpt
        end if
        filename = trim(filenameTemp)
        inquire (file=filename, exist=st)
        if ( st ) then
          open( 1000, file=filename, status='old', action='read', form='unformatted')
          rewind(1000)
          ! ATTENTION: Do not change order of loops, otherwise wrong input!
          do iband = 1, nobd(ikpt, 1)
            do idir = 1, 3
              do iBas = 1, nv(1, ikpq) + atoms%nlotot
                read(1000) z1nG(iBas, idir, iDatom, iband)
              end do ! iBas
            end do ! idir
          end do ! iband
          close(1000)
        else
          call juDFT_warn( filename//' not found! Its z1 is set to zero.', calledby='readInz1', hint='Perform calculation to gather&
            & required z1.' )
        end if
      end do ! iDeqat
    end do ! iDtype

  end subroutine readInz1

  subroutine calcFnsphVarphi(atoms, itype, oqn_lVmin, nRadFun, varphi1, varphi2, fSh, fShNsphVarphi)

    use m_types, only : t_atoms
    use m_gaunt, only : gaunt1

    implicit none

    type(t_atoms), intent(in)  :: atoms

    integer,       intent(in)  :: itype
    integer,       intent(in)  :: oqn_lVmin

    integer,       intent(in)  :: nRadFun(0:, :)
    real,          intent(in)  :: varphi1(:, :, 0:)
    real,          intent(in)  :: varphi2(:, :, 0:)
    complex,       intent(in)  :: fSh(:, :)
    complex,       intent(out) :: fShNsphVarphi(:, :, 0:, :)

    integer                    :: lmp
    integer                    :: oqn_l
    integer                    :: mqn_m
    integer                    :: iradf
    integer                    :: oqn_lV
    integer                    :: lmV_pre
    integer                    :: mqn_mV
    integer                    :: lmV
    integer                    :: mqn_m1Pr
    integer                    :: oqn_l1Pr
    integer                    :: lm1Pr
    integer                    :: imesh
    real                       :: gauntFactor
    complex                    :: vEff0Gaunt

    ! It might be slightly faster to interchange the oqn_lV and oqn_l loops, but then then bugs are easier to produce.
    lmp = 0
    do oqn_l = 0, atoms%lmax(itype)
      do mqn_m = -oqn_l, oqn_l
        ! We have due to the p loop a factor 2 without LOs and with LOs 2 + number of LOs but otherwise we would not run linearly through
        ! the arrays.
        do iradf = 1, nRadFun(oqn_l, itype)
          lmp = lmp + 1
          do oqn_lV = oqn_lVmin, atoms%lmax(itype)
            lmV_pre = oqn_lV * (oqn_lV + 1)
            do mqn_mV = -oqn_lV, oqn_lV
              lmV = lmV_pre + mqn_mV + 1
              ! The ket is only given until lmax, therefore we need no extension of rbas1/rbas2.
              ! Gaunt selection rule for m
              mqn_m1Pr = mqn_m + mqn_mV
              ! Gaunt selection rule for the l of the bra and
              ! in the non-spherical part of this routine, we have to consider that the bra (due to twofold gradient) has contributions up
              ! to lmax + 2 and the non-spherical potential which has a cut-off of lmax and is multiplied with the ket having a cutoff of
              ! lmax. The product is a quantity expanded until 2 * lmax (+ 1). Since we multiply it with the bra the product only has to
              ! calculated until lmax + 2. Combined with the Gaunt selection rules, we calculate either until a l' from where the
              ! Gaunt coefficients are zero or if this is larger than lmax + 2 we do not need them so we NOstopNOhere. On the other
              ! hand the minimal l' should be given by the lower border of the Gaunt coefficients or by m' which is also determined
              ! according to a Gaunt selection rule. l' should always be larger or equals m'. Otherwise we get lm indices which are
              ! negative.
       !       do oqn_l1Pr = max(abs(oqn_lV - oqn_l), abs(mqn_m1Pr)), min(oqn_lV + oqn_l, atoms%lmax(itype) + 2)
              do oqn_l1Pr = max(abs(oqn_lV - oqn_l), abs(mqn_m1Pr)), min(oqn_lV + oqn_l, atoms%lmax(itype))
                lm1Pr = oqn_l1Pr * (oqn_l1Pr + 1) + mqn_m1Pr
       !         gauntFactor = gaunt1( oqn_l1Pr, oqn_lV, oqn_l, mqn_m1Pr, mqn_mV, mqn_m, atoms%lmaxd + 2)
                gauntFactor = gaunt1( oqn_l1Pr, oqn_lV, oqn_l, mqn_m1Pr, mqn_mV, mqn_m, atoms%lmaxd)
                do imesh = 1, atoms%jri(itype)
                  vEff0Gaunt = fSh(imesh, lmV) * gauntFactor
                  fShNsphVarphi(1, imesh, lm1Pr, lmp) = fShNsphVarphi(1, imesh, lm1Pr, lmp) + vEff0Gaunt * varphi1(imesh, iradf, oqn_l)
                  fShNsphVarphi(2, imesh, lm1Pr, lmp) = fShNsphVarphi(2, imesh, lm1Pr, lmp) + vEff0Gaunt * varphi2(imesh, iradf, oqn_l)
                end do ! imesh
              end do ! oqn_l
            end do ! mqn_mV
          end do ! oqn_lV
        end do ! iradf
      end do ! mqn_m
    end do ! oqn_l

  end subroutine calcFnsphVarphi

  subroutine calcGrR2FinLH(atoms, sym, lathar, clnu_atom, nmem_atom, mlh_atom, r2FlhMt, r2GrFshMt)

    use m_gaunt, only : gaunt1
    use m_types
    USE m_constants
    USE m_dfpt_init, only : derivative

    implicit none

    ! Type parameters
    ! ***************
    type(t_atoms),               intent(in)  :: atoms
    type(t_sym),               intent(in)  :: sym
    type(t_sphhar),              intent(in)  :: lathar

    ! Array parameters
    ! ****************
    complex,                     intent(in)  :: clnu_atom(:, 0:, :)
    integer,                     intent(in)  :: nmem_atom(0:, :)
    integer,                     intent(in)  :: mlh_atom(:, 0:, :)
    real,                        intent(in)  :: r2FlhMt(:, 0:, :)
    complex,        allocatable, intent(out) :: r2GrFshMt(:, :, :, :)


    ! Local Scalar Variables
    ! **********************
    ! pfac    : Prefactor
    ! tGaunt  :  Gaunt coefficient
    ! itype   : Loop index for atom types
    ! ieqat   : Loop index for equivalent atoms
    ! iatom   : Loop index for all atoms
    ! imesh   : Loop index for radial mesh point
    ! mqn_m   : Magnetic quantum number m
    ! oqn_l   : Orbital quantum number l
    ! mqn_mpp : Magnetic quantum number double primed to index the natural coordinates
    ! lm      : Collective index for orbital and magnetic quantum number
    ! symType : Index of the symmetry
    ! ilh     : Loop index for different lattice harmonics (not their members!)
    ! imem    : Loop index for members of a lattice harmonics
    real                                     :: pfac
    real                                     :: tGaunt
    integer                                  :: itype
    integer                                  :: ieqat
    integer                                  :: iatom
    integer                                  :: imesh
    integer                                  :: mqn_m
    integer                                  :: oqn_l
    integer                                  :: mqn_mpp
    integer                                  :: lm
    integer                                  :: symType
    integer                                  :: ilh
    integer                                  :: imem

    ! Local Array Variables
    ! *********************
    ! rDerFlhMt    : Radial derrivative of the incoming fuction
    ! r2GrFshMtNat : Expansion coefficients of the muffin-tin gradient applied to the incoming function. The coefficients are given
    !                in natural coordinates and multiplied by $r^2$
    real,           allocatable              :: rDerFlhMt(:)
    complex,        allocatable              :: r2GrFshMtNat(:, :, :, :)


    ! Initialization of additionaly required arrays.
    allocate( r2GrFshMt(atoms%jmtd, ( atoms%lmaxd + 1)**2, atoms%nat, 3), &
            & r2GrFshMtNat(atoms%jmtd, ( atoms%lmaxd + 1)**2, atoms%nat, 3) )
    allocate( rDerFlhMt(atoms%jmtd) )
    r2GrFshMt = cmplx(0., 0.)
    r2GrFshMtNat = cmplx(0., 0.)
    rDerFlhMt = 0.

    pfac = sqrt( fpi_const / 3. )
    do mqn_mpp = -1, 1
      iatom = 0
      do itype = 1, atoms%ntype
        do ieqat = 1, atoms%neq(itype)
          iatom = iatom + 1
          symType = sym%ntypsy(iatom)
          do ilh = 0, lathar%nlh(symType)
            oqn_l = lathar%llh(ilh, symType)
            do imem = 1, nmem_atom(ilh, iatom)
              mqn_m = mlh_atom(imem, ilh, iatom)

              ! l + 1 block
              ! oqn_l - 1 to l, so oqn_l should be < lmax not <= lmax
              if ( ( abs(mqn_m - mqn_mpp) <= oqn_l + 1 ) .and. ( abs(mqn_m) <= oqn_l ) .and. (oqn_l < atoms%lmax(itype)) ) then
                lm = ( oqn_l + 1 ) * ( oqn_l + 2 ) + 1 + mqn_m - mqn_mpp
                call Derivative( r2FlhMt(:, ilh, itype), itype, atoms, rDerFlhMt )
                tGaunt = Gaunt1( oqn_l + 1, oqn_l, 1, mqn_m - mqn_mpp, mqn_m, -mqn_mpp, atoms%lmaxd )
                do imesh = 1, atoms%jri(itype)
                  r2GrFshMtNat(imesh, lm, iatom, mqn_mpp + 2) = r2GrFshMtNat(imesh, lm, iatom, mqn_mpp + 2) + pfac * (-1)**mqn_mpp &
                    &* tGaunt * (rDerFlhMt(imesh) * clnu_atom(imem, ilh, iatom) &
                    &- ((oqn_l + 2) * r2FlhMt(imesh, ilh, itype) * clnu_atom(imem, ilh, iatom) / atoms%rmsh(imesh, itype)))
                end do ! imesh
              end if ! ( abs(mqn_m - mqn_mpp) <= oqn_l + 1 ) .and. ( abs(mqn_m) <= oqn_l )

              ! l - 1 block
              if ( ( abs(mqn_m - mqn_mpp) <= oqn_l - 1 ) .and. ( abs(mqn_m) <= oqn_l ) ) then
                if ( oqn_l - 1 == -1 ) then
                  write (*, *) 'oqn_l too low'
                end if
                lm = (oqn_l - 1) * oqn_l + 1 + mqn_m - mqn_mpp
                ! This is also a trade of between storage and performance, because derivative is called redundantly, maybe store it?
                call Derivative( r2FlhMt(:, ilh, itype), itype, atoms, rDerFlhMt )
                tGaunt = Gaunt1( oqn_l - 1, oqn_l, 1, mqn_m - mqn_mpp, mqn_m, -mqn_mpp, atoms%lmaxd )
                do imesh = 1, atoms%jri(itype)
                  r2GrFshMtNat(imesh, lm, iatom, mqn_mpp + 2) = r2GrFshMtNat(imesh, lm, iatom, mqn_mpp + 2) + pfac * (-1)**mqn_mpp &
                    & * tGaunt * (rDerFlhMt(imesh)  * clnu_atom(imem, ilh, iatom) &
                    & + ((oqn_l - 1) * r2FlhMt(imesh, ilh, itype) * clnu_atom(imem, ilh, iatom) / atoms%rmsh(imesh, itype)))
                end do ! imesh
              end if ! ( abs(mqn_m - mqn_mpp) <= oqn_l - 1 ) .and. ( abs(mqn_m) <= oqn_l )
            end do ! imem
          end do ! ilh
        end do ! ieqat
      end do ! itype
    end do ! mqn_mpp

    ! Conversion from natural to cartesian coordinates
    iatom = 0
    do itype = 1, atoms%ntype
      do ieqat = 1, atoms%neq(itype)
        iatom = iatom + 1
        do oqn_l = 0, atoms%lmax(itype)
          do mqn_m = -oqn_l, oqn_l
            lm = oqn_l * (oqn_l + 1) + 1 + mqn_m
            do imesh = 1, atoms%jri(itype)
              r2GrFshMt(imesh, lm, iatom, 1:3) = matmul( Tmatrix0(1:3, 1:3), r2GrFshMtNat(imesh, lm, iatom, 1:3) )
            end do
          end do ! mqn_m
        end do ! oqn_l
      end do ! ieqat
    end do ! itype

  end subroutine calcGrR2FinLH

  subroutine CalcChannelsGrFlpNat( atoms, itype, nRadFun, flp1, flp2, delrFlp1, delrFlp2, grFlpChLout, grFlpChMout, grFlpCh1,    &
      & grFlpCh2 )

    use m_types, only : t_atoms
    use m_gaunt, only : gaunt1

    implicit none

    ! Type parameter
    type(t_atoms),               intent(in)  :: atoms

    ! Scalar parameter
    integer,                     intent(in)  :: itype

    ! Array parameter
    ! nRadFun     : number of radial functions per l and type
    ! flp1        : 1st(large) component of function with l and p(p loops over the radial function types)
    ! flp2        : 2nd(small) component of function with l and p(p loops over the radial function types)
    ! delrFlp1    : Radial derivative of the 1st component of function with l and p(p loops over the radial function types)
    ! delrFlp2    : Radial derivative of the 2nd component of function with l and p(p loops over the radial function types)
    ! grFlpChLout : Contains the resulting outgoing l of a respective scattering channel of the gradient
    ! grFlpChMout : Contains the resulting outgoing m of a respective scattering channel of the gradient
    ! grFlpCh1    : Contains the 1st (large) component of the gradient scattering channels of the input function
    ! grFlpCh2    : Contains the 2nd (small) component of the gradient scattering channels of the input function
    integer,                     intent(in)  :: nRadFun(0:, :)
    real,                        intent(in)  :: flp1(:, :, 0:)
    real,                        intent(in)  :: flp2(:, :, 0:)
    real,                        intent(in)  :: delrFlp1(:, :, 0:)
    real,                        intent(in)  :: delrFlp2(:, :, 0:)
    integer,                     intent(out) :: grFlpChLout(:, 0:)
    integer,                     intent(out) :: grFlpChMout(-atoms%lmaxd:, -1:)
    real,                        intent(out) :: grFlpCh1(:, :, :, -1:)
    real,                        intent(out) :: grFlpCh2(:, :, :, -1:)

    ! Local Variables:
    !
    ! pfac     : contains sqrt(4 pi / 3)
    ! lm          : encodes oqn_l and mqn_m
    ! idirec      : runs over 3 directions the atom can be displaced to
    ! tempGaunt1  : auxillary variable to store a Gaunt coefficient
    ! imesh       : runs over mesh points of current grid
    ! mqn_m       : magnetic quantum number m
    ! mqn_mpp     : magnetic quantum number m", also used for indexing 3 directions the atom can be displaced to
    ! oqn_l       : orbital quantum number l

    ! Local Scalar Variables
    real                                     :: pfac
    integer                                  :: oqn_l
    integer                                  :: mqn_m
    integer                                  :: iradf
    integer                                  :: lmp
    real                                     :: gauntCoeff
    integer                                  :: imesh
    integer                                  :: mqn_m2Pr


    ! Local Array Variables

    grFlpCh1(:, :, :, :) = 0.
    grFlpCh2(:, :, :, :) = 0.


    pfac = sqrt( fpi_const / 3. )

    ! Precalculate L-output channels
    do oqn_l = 0, atoms%lmax(itype)
      grFlpChLout(1, oqn_l) = oqn_l + 1
      grFlpChLout(2, oqn_l) = oqn_l - 1
    end do ! oqn_l
    ! Set this lout = 0 - 1 to abnormal value as not allowed
    grFlpChLout(2, 0) = -9999

    ! Precalculate M-output channels
    do mqn_m2Pr = -1, 1
      do mqn_m = -atoms%lmax(itype), atoms%lmax(itype)
        grFlpChMout(mqn_m, mqn_m2Pr) = mqn_m - mqn_m2Pr
      end do ! mqn_m
    end do ! mqn_m2Pr

    do mqn_m2Pr = -1, 1
      lmp = 0
      do oqn_l = 0, atoms%lmax(itype)
        do mqn_m = -oqn_l, oqn_l
          do iradf = 1, nRadFun(oqn_l, itype)
            lmp = lmp + 1

            ! scattering channel (l + 1)
            if ( (abs(mqn_m - mqn_m2Pr) <= oqn_l + 1) .and. (oqn_l < atoms%lmax(itype))) then
              gauntCoeff = Gaunt1( oqn_l + 1, oqn_l, 1, mqn_m - mqn_m2Pr, mqn_m, -mqn_m2Pr, atoms%lmax(itype) )
              do imesh = 1, atoms%jri(itype)
                ! Consider large and small relativistic components of radial solution
                grFlpCh1(imesh, 1, lmp, mqn_m2Pr) = grFlpCh1(imesh, 1, lmp, mqn_m2Pr) + pfac * (-1)**mqn_m2Pr * gauntCoeff &
                               & * (delrFlp1(imesh, iradf, oqn_l) -  oqn_l      * flp1(imesh, iradf, oqn_l) / atoms%rmsh(imesh, itype))
                grFlpCh2(imesh, 1, lmp, mqn_m2Pr) = grFlpCh2(imesh, 1, lmp, mqn_m2Pr) + pfac * (-1)**mqn_m2Pr * gauntCoeff &
                               & * (delrFlp2(imesh, iradf, oqn_l) -  oqn_l      * flp2(imesh, iradf, oqn_l) / atoms%rmsh(imesh, itype))
              enddo ! imesh
            endif ! scattering channel (l + 1)

            ! scattering channel (l - 1)
            ! This condition ensures that oqn_l = 0 is not accepted due to the emerging false condition 0 or 1 <= -1
            if ( ( abs(mqn_m - mqn_m2Pr) <= oqn_l - 1 ) ) then
              gauntCoeff = Gaunt1( oqn_l - 1, oqn_l, 1, mqn_m - mqn_m2Pr, mqn_m, -mqn_m2Pr, atoms%lmax(itype))
              do imesh = 1, atoms%jri(itype)
                ! Consider large and small relativistic components of radial solution
                grFlpCh1(imesh, 2, lmp, mqn_m2Pr) = grFlpCh1(imesh, 2, lmp, mqn_m2Pr) + pfac * (-1)**mqn_m2Pr * gauntCoeff &
                               & * (delrFlp1(imesh, iradf, oqn_l) + (oqn_l + 1) * flp1(imesh, iradf, oqn_l) / atoms%rmsh(imesh, itype))
                grFlpCh2(imesh, 2, lmp, mqn_m2Pr) = grFlpCh2(imesh, 2, lmp, mqn_m2Pr) + pfac * (-1)**mqn_m2Pr * gauntCoeff &
                               & * (delrFlp2(imesh, iradf, oqn_l) + (oqn_l + 1) * flp2(imesh, iradf, oqn_l) / atoms%rmsh(imesh, itype))
              enddo ! imesh
            endif ! scattering channel (l - 1)

          end do ! p
        end do ! mqn_m
      end do ! oqn_l
    end do ! mqn_m2Pr

  end subroutine CalcChannelsGrFlpNat

  subroutine CalcChannelsGrGrtFlpNat( atoms, itype, lmpMax, nRadFun, grFlpChLout, grFlpChMout, grFlpCh1, grFlpCh2, delrGrFlpCh1, &
      & delrGrFlpCh2, grGrtFlpChLout, grGrtFlpChMout, grGrtFlpCh1, grGrtFlpCh2 )

    use m_types, only : t_atoms
    use m_gaunt, only : gaunt1

    implicit none

    ! Type parameter
    type(t_atoms),               intent(in)  :: atoms

    ! Scalar parameter
    integer,                     intent(in)  :: itype
    integer,                     intent(in)  :: lmpMax

    ! Array parameter
    ! nRadFun        : number of radial functions per l and type
    ! grFlpChLout    : Contains the resulting outgoing l of a respective scattering channel of the input gradient
    ! grFlpChMout    : Contains the resulting outgoing m of a respective scattering channel of the input gradient
    ! grFlpCh1       : Contains the 1st (large) component of the gradient scattering channels of the input function
    ! grFlpCh2       : Contains the 2nd (small) component of the gradient scattering channels of the input function
    ! delrgrFlpCh1   : Radial derivative of the 1st (large) component of the gradient scattering channels of the input function
    ! delrgrFlpCh2   : Radial derivative of the 2nd (small) component of the gradient scattering channels of the input function
    ! grGrtFlpChLout : Resulting output l channels to the scattering channels given in grGrtFlpCh1/2
    ! grGrtFlpChMout : Resulting output m channels to the scattering channels given in grGrtFlpCh1/2
    ! grGrtFlpCh1    : 1st(large) component of the double gradient of a function given in spherical harmonics
    ! grGrtFlpCh2    : 2nd(small) component of the double gradient of a function given in spherical harmonics
    ! Array parameter
    integer,                     intent(in)  :: nRadFun(0:, :)
    integer,                     intent(in)  :: grFlpChLout(:, 0:)
    integer,                     intent(in)  :: grFlpChMout(-atoms%lmaxd:, -1:)
    real,                        intent(in)  :: grFlpCh1(:, :, :, -1:)
    real,                        intent(in)  :: grFlpCh2(:, :, :, -1:)
    real,                        intent(in)  :: delrGrFlpCh1(:, :, :, -1:)
    real,                        intent(in)  :: delrGrFlpCh2(:, :, :, -1:)
    integer,                     intent(out) :: grGrtFlpChLout(:, 0:)
    integer,                     intent(out) :: grGrtFlpChMout(-atoms%lmaxd:, -1:, -1:)
    real,                        intent(out) :: grGrtFlpCh1(:, :, :, -1:, -1:)
    real,                        intent(out) :: grGrtFlpCh2(:, :, :, -1:, -1:)

    ! Scalar Variables
    real                                     :: pfac
    integer                                  :: mqn_m2PrC
    integer                                  :: mqn_m2PrR
    integer                                  :: lmp
    integer                                  :: oqn_l
    integer                                  :: mqn_m
    integer                                  :: iradf
    integer                                  :: iChGrGrt
    integer                                  :: iChGr
    integer                                  :: oqn_l3Pr
    integer                                  :: mqn_m3Pr
    real                                     :: gauntCoeff
    integer                                  :: imesh

    ! Array Variables

    pfac = sqrt( fpi_const / 3 )

    do mqn_m2PrC = -1, 1
      do mqn_m2PrR = -1, 1
        ! Precalculate output magnetic quantum number m
        do mqn_m = -atoms%lmax(itype), atoms%lmax(itype)
          ! We calculate grad grad^T. Hence, the incoming gradient is a column-vector
          mqn_m3Pr = grFlpChMout(mqn_m, mqn_m2PrC)
          ! The resulting m of the double gradient is calculated from the output m of the simple gradient.
          grGrtFlpChMout(mqn_m, mqn_m2PrR, mqn_m2PrC) = mqn_m3Pr - mqn_m2PrR
        end do ! mqn_m
        lmp = 0
        ! If we loop over the l, m quantum numbers of the original function and the scattering channels of its gradient, we
        ! account for all contributions.
        do oqn_l = 0, atoms%lmax(itype)
          do mqn_m = -oqn_l, oqn_l
            mqn_m3Pr = grFlpChMout(mqn_m, mqn_m2PrC)
            do iradf = 1, nRadFun(oqn_l, itype)
              lmp = lmp + 1
              iChGrGrt = 0
              do iChGr = 1, 2
                oqn_l3Pr = grFlpChLout(iChGr, oqn_l)
                if ((oqn_l3Pr < 0 ) .or. (oqn_l3Pr > atoms%lmax(itype) + 1)) cycle
                iChGrGrt = iChGrGrt + 1
                ! Still for the Gaunt coefficient conditions we need the resulting m and l quantum numbers of the simple gradient
                ! scattering channels
                if ( ( abs(mqn_m3Pr - mqn_m2PrR) <= oqn_l3Pr + 1 ) .and. ( abs(mqn_m3Pr) <= oqn_l3Pr ) ) then
                  grGrtFlpChLout(iChGrGrt, oqn_l) = oqn_l3Pr + 1
                  gauntCoeff = Gaunt1( oqn_l3Pr + 1, oqn_l3Pr, 1, mqn_m3Pr - mqn_m2PrR, mqn_m3Pr, -mqn_m2PrR, atoms%lmax(itype) + 2)
                  do imesh = 1, atoms%jri(itype)
                    ! We do not need two quantities here for large and small component but to make it consistent we use two!
                    grGrtFlpCh1(imesh, iChGrGrt, lmp, mqn_m2PrR, mqn_m2PrC) = grGrtFlpCh1(imesh, iChGrGrt, lmp, mqn_m2PrR, mqn_m2PrC) &
                                                & + pfac * (-1)**mqn_m2PrR * gauntCoeff * (delrGrFlpCh1(imesh, iChGr, lmp, mqn_m2PrC) &
                                                & - oqn_l3Pr * grFlpCh1(imesh, iChGr, lmp, mqn_m2PrC) / atoms%rmsh(imesh, itype))
                    grGrtFlpCh2(imesh, iChGrGrt, lmp, mqn_m2PrR, mqn_m2PrC) = grGrtFlpCh2(imesh, iChGrGrt, lmp, mqn_m2PrR, mqn_m2PrC) &
                                                & + pfac * (-1)**mqn_m2PrR * gauntCoeff * (delrGrFlpCh2(imesh, iChGr, lmp, mqn_m2PrC) &
                                                & - oqn_l3Pr * grFlpCh2(imesh, iChGr, lmp, mqn_m2PrC) / atoms%rmsh(imesh, itype))
                  end do ! imesh
                end if ! l + 1-scattering channel

                iChGrGrt = iChGrGrt + 1
                ! Still for the Gaunt coefficient conditions we need the resulting m and l quantum numbers of the simple gradient
                ! scattering channels
                if ( (abs(mqn_m3Pr - mqn_m2PrR) <= oqn_l3Pr - 1) .and. ( abs(mqn_m3Pr) <= oqn_l3Pr ) ) then
                  grGrtFlpChLout(iChGrGrt, oqn_l) = oqn_l3Pr - 1
                  gauntCoeff = Gaunt1( oqn_l3Pr - 1, oqn_l3Pr, 1, mqn_m3Pr - mqn_m2PrR, mqn_m3Pr, -mqn_m2PrR, atoms%lmax(itype) + 1)
                  do imesh = 1, atoms%jri(itype)
                    ! We do not need two quantities here for large and small component but to make it consistent we use two!
                    grGrtFlpCh1(imesh, iChGrGrt, lmp, mqn_m2PrR, mqn_m2PrC) = grGrtFlpCh1(imesh, iChGrGrt, lmp, mqn_m2PrR, mqn_m2PrC) &
                                                & + pfac * (-1)**mqn_m2PrR * gauntCoeff * (delrGrFlpCh1(imesh, iChGr, lmp, mqn_m2PrC) &
                                                & + (oqn_l3Pr + 1) * grFlpCh1(imesh, iChGr, lmp, mqn_m2PrC) / atoms%rmsh(imesh, itype))
                    grGrtFlpCh2(imesh, iChGrGrt, lmp, mqn_m2PrR, mqn_m2PrC) = grGrtFlpCh2(imesh, iChGrGrt, lmp, mqn_m2PrR, mqn_m2PrC) &
                                                & + pfac * (-1)**mqn_m2PrR * gauntCoeff * (delrGrFlpCh2(imesh, iChGr, lmp, mqn_m2PrC) &
                                                & + (oqn_l3Pr + 1) * grFlpCh2(imesh, iChGr, lmp, mqn_m2PrC) / atoms%rmsh(imesh, itype))
                  end do ! imesh
                end if ! l - 1-scattering channel
              end do ! iChGr
            end do ! iradf
          end do ! mqn_m
        end do ! oqn_l
      end do ! mqn_m2PrR
    end do ! mqn_m2PrC

  end subroutine CalcChannelsGrGrtFlpNat

  subroutine calcFnsphGrVarphi(atoms, itype, oqn_lVmin, mqn_m2Pr, nRadFun, grVarphiCh1, grVarphiCh2, grVarphiChLout, grVarphiChMout, fSh,&
                                                                                                                    & fShNsphVarphi)

    use m_types, only : t_atoms
    use m_gaunt, only : gaunt1

    implicit none

    type(t_atoms), intent(in)  :: atoms

    integer,       intent(in)  :: itype
    integer,       intent(in)  :: oqn_lVmin
    integer,       intent(in)  :: mqn_m2Pr

    integer,       intent(in)  :: nRadFun(0:, :)
    real,          intent(in)  :: grVarphiCh1(:, :, :, -1:)
    real,          intent(in)  :: grVarphiCh2(:, :, :, -1:)
    integer,       intent(in)  :: grVarphiChLout(:, 0:)
    integer,       intent(in)  :: grVarphiChMout(-atoms%lmaxd:, -1:)
    complex,       intent(in)  :: fSh(:, :)
    complex,       intent(out) :: fShNsphVarphi(:, :, 0:, :)

    integer                    :: lmp
    integer                    :: oqn_l
    integer                    :: mqn_m
    integer                    :: iradf
    integer                    :: oqn_lV
    integer                    :: lmV_pre
    integer                    :: mqn_mV
    integer                    :: lmV
    integer                    :: mqn_m1Pr
    integer                    :: oqn_l1Pr
    integer                    :: lm1Pr
    integer                    :: imesh
    real                       :: gauntFactor
    complex                    :: vEff0Gaunt
    integer                    :: mqn_m3Pr
    integer                    :: ichan
    integer                    :: oqn_l3Pr

    fShNsphVarphi = cmplx(0., 0.)
    ! It might be slightly faster to interchange the oqn_lV and oqn_l loops, but then then bugs are easier to produce.
    lmp = 0
    do oqn_l = 0, atoms%lmax(itype)
      do mqn_m = -oqn_l, oqn_l
        mqn_m3Pr = grVarphiChMout(mqn_m, mqn_m2Pr)
        ! We have due to the p loop a factor 2 without LOs and with LOs 2 + number of LOs but otherwise we would not run linearly through
        ! the arrays.
        do iradf = 1, nRadFun(oqn_l, itype)
          lmp = lmp + 1
          do ichan = 1, 2
            oqn_l3Pr = grVarphiChLout(ichan, oqn_l)
            if (oqn_l3Pr < 0 .or. abs(mqn_m3Pr) > oqn_l3Pr .or. oqn_l3Pr > atoms%lmax(itype)) cycle
            do oqn_lV = oqn_lVmin, atoms%lmax(itype)
              lmV_pre = oqn_lV * (oqn_lV + 1)
              do mqn_mV = -oqn_lV, oqn_lV
                lmV = lmV_pre + mqn_mV + 1
                ! The ket is only given until lmax, therefore we need no extension of rbas1/rbas2.
                ! Gaunt selection rule for m
                mqn_m1Pr = mqn_m3Pr + mqn_mV
                ! Gaunt selection rule for the l of the bra and
                ! in the non-spherical part of this routine, we have to consider that the bra (due to twofold gradient) has contributions up
                ! to lmax + 2 and the non-spherical potential which has a cut-off of lmax and is multiplied with the ket having a cutoff of
                ! lmax. The product is a quantity expanded until 2 * lmax (+ 1). Since we multiply it with the bra the product only has to
                ! calculated until lmax + 2. Combined with the Gaunt selection rules, we calculate either until a l' from where the
                ! Gaunt coefficients are zero or if this is larger than lmax + 2 we do not need them so we NOstopNOhere. On the other
                ! hand the minimal l' should be given by the lower border of the Gaunt coefficients or by m' which is also determined
                ! according to a Gaunt selection rule. l' should always be larger or equals m'. Otherwise we get lm indices which are
                ! negative.
!                do oqn_l1Pr = max(abs(oqn_lV - oqn_l3Pr), abs(mqn_m1Pr)), min(oqn_lV + oqn_l3Pr, atoms%lmax(itype) + 2)
                do oqn_l1Pr = max(abs(oqn_lV - oqn_l3Pr), abs(mqn_m1Pr)), min(oqn_lV + oqn_l3Pr, atoms%lmax(itype))
                  lm1Pr = oqn_l1Pr * (oqn_l1Pr + 1) + mqn_m1Pr
!                  gauntFactor = gaunt1( oqn_l1Pr, oqn_lV, oqn_l3Pr, mqn_m1Pr, mqn_mV, mqn_m3Pr, atoms%lmaxd + 2)
                  gauntFactor = gaunt1( oqn_l1Pr, oqn_lV, oqn_l3Pr, mqn_m1Pr, mqn_mV, mqn_m3Pr, atoms%lmaxd)
                  do imesh = 1, atoms%jri(itype)
                    vEff0Gaunt = fSh(imesh, lmV) * gauntFactor
                    fShNsphVarphi(1, imesh, lm1Pr, lmp) = fShNsphVarphi(1, imesh, lm1Pr, lmp) + vEff0Gaunt &
                                                                                        & * grVarphiCh1(imesh, ichan, lmp, mqn_m2Pr)
                    fShNsphVarphi(2, imesh, lm1Pr, lmp) = fShNsphVarphi(2, imesh, lm1Pr, lmp) + vEff0Gaunt &
                                                                                        & * grVarphiCh2(imesh, ichan, lmp, mqn_m2Pr)
                  end do ! imesh
                end do ! oqn_l
              end do ! mqn_mV
            end do ! oqn_lV
          end do ! ichan
        end do ! iradf
      end do ! mqn_m
    end do ! oqn_l

  end subroutine calcFnsphGrVarphi

end module m_jpSetupDynMatSF
