!--------------------------------------------------------------------------------
! Copyright (c) 2021 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------
MODULE m_dfpt_init

    USE m_juDFT
    USE m_types
    USE m_constants
    USE m_radfun
    USE m_radflo
    USE m_eig66_io, ONLY : read_eig
    USE m_intgrf,     ONLY : intgrf_init

    IMPLICIT NONE

CONTAINS
    SUBROUTINE dfpt_init(juPhon, sym, input, atoms, sphhar, stars, cell, noco, nococonv, &
                       & kpts, fmpi, results, enpara, rho, vTot, eig_id, nvfull, usdus, rho0, &
                       & grRho0, vTot0, grVTot0, ngdp, El, recG, ngdp2km, gdp2Ind, gdp2iLim, GbasVec, ilst, nRadFun, iloTable, ilo2p, &
                       & uuilonout, duilonout, ulouilopnout, kveclo, rbas1, rbas2, gridf, z0, grVxcIRKern, dKernMTGPts, &
                       & gausWts, ylm, qpwcG, rho1MTCoreDispAt, grVeff0MT_init, grVeff0MT_main, grVext0IR_DM, grVext0MT_DM, &
                       & grVCoul0IR_DM_SF, grVCoul0MT_DM_SF, grVeff0IR_DM, grVeff0MT_DM, tdHS0, loosetdout, nocc, rhoclean, oldmode, xcpot, grRho)

        USE m_jpGrVeff0, ONLY : GenGrVeff0
        USE m_npy

        TYPE(t_juPhon),   INTENT(IN)  :: juPhon
        TYPE(t_sym),      INTENT(IN)  :: sym
        TYPE(t_input),    INTENT(IN)  :: input
        TYPE(t_atoms),    INTENT(IN)  :: atoms
        TYPE(t_sphhar),   INTENT(IN)  :: sphhar
        TYPE(t_stars),    INTENT(IN)  :: stars
        TYPE(t_cell),     INTENT(IN)  :: cell
        TYPE(t_noco),     INTENT(IN)  :: noco
        TYPE(t_nococonv), INTENT(IN)  :: nococonv
        TYPE(t_kpts),     INTENT(IN)  :: kpts
        TYPE(t_mpi),      INTENT(IN)  :: fmpi
        TYPE(t_results),  INTENT(IN)  :: results
        TYPE(t_enpara),   INTENT(IN)  :: enpara
        TYPE(t_potden),   INTENT(IN)  :: rho, vTot, rhoclean
        INTEGER,          INTENT(IN)  :: eig_id
        INTEGER,          INTENT(IN)  :: nvfull(:, :)

        TYPE(t_usdus),    INTENT(OUT) :: usdus
        TYPE(t_jpPotden), INTENT(OUT) :: rho0, grRho0, vTot0, grVTot0
        TYPE(t_potden),   INTENT(OUT)  :: grRho

        INTEGER,          INTENT(OUT) :: ngdp

        ! New input:
        LOGICAL, INTENT(IN)        :: oldmode
        CLASS(t_xcpot), INTENT(IN) :: xcpot

        REAL,             ALLOCATABLE, INTENT(OUT) :: El(:, :, :, :)
        INTEGER,          ALLOCATABLE, INTENT(OUT) :: recG(:, :)
        INTEGER,                       INTENT(OUT) :: ngdp2km
        INTEGER,          ALLOCATABLE, INTENT(OUT) :: gdp2Ind(:, :, :)
        INTEGER,                       INTENT(OUT) :: gdp2iLim(2, 3)
        INTEGER,          ALLOCATABLE, INTENT(OUT) :: GbasVec(:, :)
        INTEGER,          ALLOCATABLE, INTENT(OUT) :: ilst(:, :, :)
        INTEGER,          ALLOCATABLE, INTENT(OUT) :: nRadFun(:,:)
        INTEGER,          ALLOCATABLE, INTENT(OUT) :: iloTable(:, :, :)
        INTEGER,          ALLOCATABLE, INTENT(OUT) :: ilo2p(:, :)
        REAL,             ALLOCATABLE, INTENT(OUT) :: uuilonout(:, :)
        REAL,             ALLOCATABLE, INTENT(OUT) :: duilonout(:, :)
        REAL,             ALLOCATABLE, INTENT(OUT) :: ulouilopnout(:, :, :)
        INTEGER,          ALLOCATABLE, INTENT(OUT) :: kveclo(:,:)
        REAL,             ALLOCATABLE, INTENT(OUT) :: rbas1(:, :, :, :, :)
        REAL,             ALLOCATABLE, INTENT(OUT) :: rbas2(:, :, :, :, :)
        REAL,             ALLOCATABLE, INTENT(OUT) :: gridf(:, :)
        COMPLEX,          ALLOCATABLE, INTENT(OUT) :: z0(:, :, :, :)
        complex,           allocatable, intent(out) :: grVxcIRKern(:)
        real,              allocatable, intent(out) :: dKernMTGPts(:, :, :)
        real,              allocatable, intent(out) :: gausWts(:)
        complex,           allocatable, intent(out) :: ylm(:, :)
        complex,           allocatable, intent(out) :: qpwcG(:, :)
        complex,           allocatable, intent(out) :: rho1MTCoreDispAt(:, :, :, :)
        complex,           allocatable, intent(out) :: grVeff0MT_init(:, :, :, :)
        complex,           allocatable, intent(out) :: grVeff0MT_main(:, :, :, :)
        complex,           allocatable, intent(out) :: grVext0IR_DM(:, :)
        complex,           allocatable, intent(out) :: grVext0MT_DM(:, :, :, :)
        complex,           allocatable, intent(out) :: grVCoul0IR_DM_SF(:, :)
        complex,           allocatable, intent(out) :: grVCoul0MT_DM_SF(:, :, :, :)
        complex,           allocatable, intent(out) :: grVeff0IR_DM(:, :)
        complex,           allocatable, intent(out) :: grVeff0MT_DM(:, :, :, :)
        type(t_tlmplm),                 intent(out) :: tdHS0
        COMPLEX, ALLOCATABLE,           INTENT(OUT) :: loosetdout(:, :, :, :)
        INTEGER, ALLOCATABLE,           INTENT(OUT) :: nocc(:, :)

        TYPE(t_lapw) :: lapw
        TYPE(t_mat)  :: zMat

        TYPE(t_cdnvalJob) :: cdnvaljob

        INTEGER      :: nG, nSpins, nbasfcn, iSpin, ik, iG, indexG, ifind, nk, iType, ilo, l_lo, oqn_l, iGrid, iOrd
        INTEGER      :: nodeu, noded, xInd, yInd, zInd, iStar
        LOGICAL      :: l_real, addG, harSw, extSw, xcSw, testGoldstein, grRhoTermSw
        REAL         :: bkpt(3), wronk

        INTEGER, ALLOCATABLE :: GbasVec_temp(:, :), ev_list(:)
        REAL,    ALLOCATABLE :: we(:), eig(:)
        REAL,    ALLOCATABLE :: f(:, :, :)
        REAL,    ALLOCATABLE :: g(:, :, :)
        REAL,    ALLOCATABLE :: flo(:, :, :)
        real,              allocatable              :: acoff(:)
        real,              allocatable              :: alpha(:)
        complex,           allocatable              :: grVeff0IRDummy(:, :)
        complex,           allocatable              :: grVeff0MTAdd_init(:, :, :, :)
        complex,           allocatable              :: grVeff0IR_DMhxc(:, :)
        INTEGER,           ALLOCATABLE              :: GbasVec_eig_loc(:, :, :, :)

!#ifdef CPP_MPI
!        INTEGER ierr
!#endif

        call uniteEnergyParameters( atoms, enpara, input, El )

        ! TODO: This is a hack:
        ALLOCATE(kveclo(atoms%nlotot, kpts%nkpt))
        kveclo = 0

        IF (noco%l_noco) THEN
           nSpins = 1
        ELSE
           nSpins = input%jspins
        ENDIF

        ALLOCATE(ilst(MAXVAL(nvfull), kpts%nkpt, input%jspins), GbasVec_temp(3, 0), GbasVec(3, 0) )

        !nbasfcn = MERGE(MAXVAL(nvfull(1,:))+MAXVAL(nvfull(2,:))+2*atoms%nlotot,MAXVAL(nvfull(1,:))+atoms%nlotot,noco%l_noco)

        nbasfcn = MAXVAL(nvfull(1,:))+atoms%nlotot

        ALLOCATE(z0(nbasfcn, input%neig, kpts%nkpt, nSpins))
        ALLOCATE(nocc(kpts%nkpt, input%jspins))

        ALLOCATE(GbasVec_eig_loc(3, input%neig, kpts%nkpt, MERGE(1,input%jspins,noco%l_noco)))
        GbasVec_eig_loc = 0

        z0 = CMPLX(0.0,0.0)
        nocc = 0

        l_real = sym%invs.AND.(.NOT.noco%l_soc).AND.(.NOT.noco%l_noco).AND.atoms%n_hia==0

        !CALL zMat%init(l_real, nbasfcn, input%neig)

        DO iSpin = 1, nspins
            CALL cdnvalJob%init(fmpi, input, kpts, noco, results, iSpin)
            DO ik = 1,kpts%nkpt

                !nk=cdnvalJob%k_list(ik)

                !bkpt=kpts%bk(:, nk)
                bkpt=kpts%bk(:, ik)

                !CALL lapw%init(input, noco, nococonv, kpts, atoms, sym, nk, cell, .FALSE., fmpi)
                CALL lapw%init(input, noco, nococonv, kpts, atoms, sym, ik, cell, fmpi)
                GbasVec_eig_loc(:, :lapw%nv(iSpin) + atoms%nlotot, ik, iSpin) = lapw%gvec(:, :, iSpin)

                ! Kinda like this for lapw%kvec and nkvec. Irrelevant for now.
                !DO lo = 1,nlo(n)
                !    kveclo(nkvec_sv+1:nkvec_sv+nkvec(lo,1)) =
                !    +                                            kvec(1:nkvec(lo,1),lo)
                !    nkvec_sv = nkvec_sv+nkvec(lo,1)
                !    nkvec(lo,:) = 0
                !END DO

                ev_list = cdnvaljob%compact_ev_list(ik, l_empty = .FALSE.)
                nocc(ik, iSpin) = SIZE(ev_list)
                ev_list = cdnvaljob%compact_ev_list(ik, l_empty = .TRUE.)
                !nocc(ik, iSpin) = cdnvaljob%noccbd(nk)
                !we  = cdnvalJob%weights(ev_list, nk)
                !eig = results%eig(ev_list, nk, iSpin)

                !IF (fmpi%irank == 0) THEN
                nbasfcn = MERGE(lapw%nv(1)+lapw%nv(2)+2*atoms%nlotot,lapw%nv(1)+atoms%nlotot,noco%l_noco)
                CALL zMat%init(l_real, nbasfcn, SIZE(ev_list))
                CALL read_eig(eig_id, ik, iSpin, list = ev_list, zmat=zMat)
                    !CALL read_eig(eig_id, ik, iSpin, neig = results%neig(ik, iSpin), &
                    !                                & eig = results%eig(:, ik, iSpin))
                !END IF

                IF (l_real) THEN
                    !z0(:nvfull(iSpin, ik), :results%neig(ik, iSpin), ik, iSpin) = CMPLX(1.0,0.0) * zMat%data_r(:nvfull(iSpin, ik), :results%neig(ik, iSpin))
                    z0(:lapw%nv(iSpin), :lapw%nv(iSpin), ik, iSpin) = CMPLX(1.0,0.0) * zMat%data_r(:lapw%nv(iSpin), :lapw%nv(iSpin))
                ELSE
                    !z0(:nvfull(iSpin, ik), :results%neig(ik, iSpin), ik, iSpin) = zMat%data_c(:nvfull(iSpin, ik), :results%neig(ik, iSpin))
                    z0(:lapw%nv(iSpin), :lapw%nv(iSpin), ik, iSpin) = zMat%data_c(:lapw%nv(iSpin), :lapw%nv(iSpin))
                END IF

!#ifdef CPP_MPI
!                CALL MPI_BARRIER(fmpi%mpi_comm,ierr)
!#endif
            END DO
        END DO

        ! Store basis vectors G in a compressed way, as they occur multiple times
        ! for different k-points. The array ilst stores the vector's index if it
        ! has occured for a different k-point or spin already.
        indexG = 0
        ilst = -1
        DO iSpin = 1, input%jspins
            DO ik = 1, kpts%nkpt
                DO iG = 1, nvfull(iSpin, ik)
                    IF (uBound(GbasVec, 2).EQ.0) then
                        addG = .true.
                    ELSE
                        DO ifind =  1, uBound(GbasVec, 2)
                            IF (.NOT.ALL(ABS(GbasVec(:, ifind) - GbasVec_eig_loc(:, iG, ik, iSpin)).LE.1E-8)) THEN
                                addG = .TRUE.
                            ELSE
                                ilst(iG, ik, iSpin) = ifind
                                addG = .FALSE.
                                EXIT
                            END IF
                        END DO ! ifind
                    END IF

                    IF (addG) THEN
                        IF (uBound(GbasVec, 2).EQ.0) THEN
                            DEALLOCATE(GbasVec, GbasVec_temp)
                            ALLOCATE(GbasVec(3, 1), GbasVec_temp(3, 1))
                        ELSE
                            GbasVec_temp = GbasVec
                            DEALLOCATE(GbasVec)
                            ALLOCATE(GbasVec(3, indexG + 1))
                            GbasVec(:, :indexG) = GbasVec_temp
                        END IF

                        indexG = indexG + 1
                        ilst(iG, ik, iSpin) = indexG
                        GbasVec(1, indexG) = GbasVec_eig_loc(1, iG, ik, iSpin)
                        GbasVec(2, indexG) = GbasVec_eig_loc(2, iG, ik, iSpin)
                        GbasVec(3, indexG) = GbasVec_eig_loc(3, iG, ik, iSpin)

                        IF (uBound(GbasVec, 2).NE.0) THEN
                            DEALLOCATE(GbasVec_temp)
                            ALLOCATE(GbasVec_temp(3, indexG))
                        END IF
                    END IF
                END DO ! iG
            END DO ! ik
        END DO ! iSpin

        DEALLOCATE(GbasVec_temp)

        allocate( nRadFun(0:atoms%lmaxd, atoms%ntype) )
        ! The number of radial functions p for a given atom type and orbital quantum number l is at least 2 ( u and udot ).
        nRadFun = 2

        allocate( iloTable(2 + atoms%nlod, 0: atoms%lmaxd, atoms%ntype) )
        ! Mapping array gives zero if for a given atom type and orbital quantum number l there is no local orbital radial solution u_LO
        iloTable = 0

        ! Generate Arrays which count the number of radial solutions for a given atom type and orbital quantum number l and give the
        ! number of the local orbital if the index which counts the radial solutions for a given atom type and orbital quantum number l
        ! is given as well as the orbital quantum number l and the atom type themselves.
        ALLOCATE(ilo2p(MAXVAL(atoms%nlo), atoms%ntype))
        ilo2p = 0
        DO iSpin = 1, input%jspins
            DO iType = 1, atoms%ntype
                DO ilo = 1, atoms%nlo(iType)
                    l_lo = atoms%llo(ilo, iType)
                    nRadFun(l_lo, iType) = nRadFun(l_lo, iType) + 1
                    ilo2p(ilo, iType) = nRadFun(l_lo, iType)
                    iloTable(nRadFun(l_lo, iType), l_lo, iType) = ilo
                END DO ! ilo
            END DO ! itype
        END DO ! isp

        ALLOCATE(rbas1(atoms%jmtd, 2 + atoms%nlod, 0:atoms%lmaxd, atoms%ntype, input%jspins)) !maybe samller dimension see Markus
        ALLOCATE(rbas2(atoms%jmtd, 2 + atoms%nlod, 0:atoms%lmaxd, atoms%ntype, input%jspins))

        ! This is initialized to the neutral element of multiplication to avoid division by zero
        rbas1 = 1
        rbas2 = 1

        ALLOCATE(f(atoms%jmtd, 2, 0:atoms%lmaxd), &
               & g(atoms%jmtd, 2, 0:atoms%lmaxd) ) ! second index runs over vector-component
        ALLOCATE(flo(atoms%jmtd, 2, atoms%nlod) )  ! second index runs over vector-component

        CALL usdus%init(atoms,input%jspins)
        ALLOCATE(uuilonout(atoms%nlod, atoms%ntype))
        ALLOCATE(duilonout(atoms%nlod, atoms%ntype))
        ALLOCATE(ulouilopnout(atoms%nlod, atoms%nlod, atoms%ntype))
        uuilonout    = usdus%uuilon(:, :, 1)
        duilonout    = usdus%duilon(:, :, 1)
        ulouilopnout = usdus%ulouilopn(:, :, :, 1)

        ! Generate radial solutions with from Fleur recycled routines radfun and radflo. Rearrange them into rbas1 and rbas2 to not give
        ! LOs a special treatment but to only loop over all radial functions and treat LOs implicitely
        DO iSpin = 1, input%jspins
          DO iType = 1, atoms%ntype
            DO oqn_l = 0, atoms%lmax(itype)

              ! The El from oqn_l = 3 upwards are the same, so we just use El(oqn_l = 4) here.

              ! Juphon
              !call radfun( oqn_l, enpara%el0(oqn_l, itype, isp), V0sFleur%vr0(1, itype, 1), atoms%jri(itype), atoms%rmsh(1, itype), &
              !  & atoms%dx(itype), atoms%jmtd, f(1, 1, oqn_l), g(1, 1, oqn_l), us(oqn_l, itype), dus(oqn_l, itype), uds(oqn_l, itype), &
              !  & duds(oqn_l, itype), ddn(oqn_l, itype), nodeu, noded, wronk )

              ! Max_FLEUR:
              !       CALL radfun(l,iType,jspin,enpara%el0(l,iType,jspin),vrTmp,atoms,&
                !          f(1,1,l),g(1,1,l),usdus,nodeu,noded,wronk) !!! Some Hubbard shenanigans here.

              ! new:
              CALL radfun(oqn_l, iType, iSpin, enpara%el0(oqn_l ,iType, iSpin), vTot%mt(:, 0, iType, iSpin), atoms, &
                        & f(1, 1, oqn_l), g(1, 1, oqn_l), usdus, nodeu, noded, wronk)

              ! Juphon
              !call radflo( atoms%ntype, atoms%nlod, dimens%jspd, atoms%jmtd, atoms%lmaxd, itype, isp, enpara%ello0(1, 1, isp), &
              !  & V0sFleur%vr0(1, itype, 1), atoms%jri(itype), atoms%rmsh(1, itype), atoms%dx(itype), f, g, atoms%llo, atoms%nlo, &
              !  & atoms%l_dulo(1, itype), mpi%irank, atoms%ulo_der, ulos, dulos, uulon, dulon, uloulopn, uuilon, duilon, ulouilopn, flo)

              ! Max_FLEUR:
              !        CALL radflo(atoms,iType,jspin,enpara%ello0(1,1,jspin),vTot%mt(:,0,iType,jspin),f,g,fmpi,&
                !          usdus,usdus%uuilon(1,1,jspin),usdus%duilon(1,1,jspin),usdus%ulouilopn(1,1,1,jspin),flo)

              ! new:
              IF (atoms%nlo(iType).GE.1) THEN
                  CALL radflo(atoms, iType, iSpin, enpara%ello0(1, 1, iSpin), vTot%mt(:, 0, iType, iSpin), &
                            & f, g, fmpi, usdus, usdus%uuilon(:, :, iSpin), usdus%duilon(:, :, iSpin), usdus%ulouilopn(:, :, :, iSpin), flo)
              END IF

              DO igrid = 1, atoms%jri(itype) ! was jmtd in former times, which is not necessary

                rbas1(iGrid, 1, oqn_l, iType, iSpin) = f(iGrid, 1, oqn_l)
                rbas2(iGrid, 1, oqn_l, iType, iSpin) = f(iGrid, 2, oqn_l)
                rbas1(iGrid, 2, oqn_l, iType, iSpin) = g(iGrid, 1, oqn_l)
                rbas2(iGrid, 2, oqn_l, iType, iSpin) = g(iGrid, 2, oqn_l)

              END DO
              DO iOrd = 3, nRadFun(oqn_l, iType)
                  ! p > 2 : LO indices
                DO iGrid = 1, atoms%jri(iType)
                 rbas1(iGrid, iOrd, oqn_l, iType, iSpin) = flo(iGrid, 1, iloTable(iOrd, oqn_l, iType))
                 rbas2(iGrid, iOrd, oqn_l, iType, iSpin) = flo(iGrid, 2, iloTable(iOrd, oqn_l, iType))
                END DO
              END DO
            END DO
          END DO
        END DO

        CALL Intgrf_init(atoms%ntype, atoms%jmtd, atoms%jri, atoms%dx, atoms%rmsh, gridf)

        ! Free rho and V of their additional factors.
        !do isp = 1,input%jspins
        !  do itype = 1,atoms%ntype
        !    do imesh = 1,atoms%jri(itype)
        !      V0Fleur%vr(imesh,0,itype,isp) = V0Fleur%vr(imesh,0,itype,isp)*sqrt(fpi)/atoms%rmsh(imesh,itype)
        !    end do
        !  end do
        !end do

        ! Initialize unsymmetrized potdens and G vectors.

        call genPotDensGvecs(stars, cell, input, ngdp, ngdp2km, recG, gdp2Ind, gdp2iLim, .false.)

        call calcIRdVxcKern(stars, recG, ngdp, rho%pw(:, 1), grVxcIRKern, oldmode, xcpot)

        call calcMTdVxcKern(atoms, sphhar, sym, rhoclean%mt(:, :, :, 1), sphhar%nmem, sphhar%clnu, sphhar%mlh, gausWts, ylm, dKernMTGPts, oldmode, xcpot)

        ! Calculate parameters of Gauss curve (pseudo-core density in MT) and Fourier transform of pseudo-core density for IR
        call calcPsDensMT( fmpi, atoms, cell, sym, stars, input, ngdp, acoff, alpha, qpwcG, recG )

        ! Calculate the gradient of the Gauss function representing the pseudo-core density. This term is used to correct the core-tail
        ! correction in the displaced atom itself, as there we know the real core-density that can be derived; see the calculation of
        ! the linear charge density variation.
        allocate( rho1MTCoreDispAt(atoms%jmtd, 1 : 4, atoms%nat, 3) )
        rho1MTCoreDispAt = cmplx(0.0, 0.0)
        call calcPDCinAlph( atoms, alpha, acoff, rho1MTCoreDispAt )

        !nG=COUNT(stars%ig>0)
        nG = ngdp
        CALL rho0%init_jpPotden(1, 0, nG, atoms%jmtd, atoms%lmaxd, atoms%nat, input%jspins, .FALSE.)
        CALL grRho0%init_jpPotden(3, 0, nG, atoms%jmtd, atoms%lmaxd + juPhon%jplPlus, atoms%nat, input%jspins, .FALSE.)
        CALL vTot0%init_jpPotden(1, 0, nG, atoms%jmtd, atoms%lmaxd, atoms%nat, input%jspins, .TRUE.)
        CALL grVTot0%init_jpPotden(3, 0, nG, atoms%jmtd, atoms%lmaxd + juPhon%jplPlus, atoms%nat, input%jspins, .TRUE.)

        ! Unpack the star coefficients onto a G-representation.
        !CALL stars_to_pw(sym, stars, input%jspins, nG, rho%pw, rho0%pw)
        !CALL stars_to_pw(sym, stars, input%jspins, nG, vTot%pw, vTot0%pw)
        !CALL stars_to_pw(sym, stars, input%jspins, nG, vTot%pw_w, vTot0%pw_w)
        CALL convertStar2G(rho%pw(:, 1), rho0%pw, stars, nG, recG, cell)
        CALL convertStar2G(vTot%pw(:, 1), vTot0%pw, stars, nG, recG, cell)
        CALL convertStar2G(vTot%pw_w(:, 1), vTot0%pw_w, stars, nG, recG, cell)

        ! Unpack the lattice harmonics onto spherical harmonics.
        ! Radial factors: rho has r^2 globally in l [2], vTot has r/sqrt(4pi) in l=0 [0].
        ! The factors are removed in the transformation.
        CALL lh_to_sh(sym, atoms, sphhar, input%jspins, 2, rho%mt, rho0%mt)
        CALL lh_to_sh(sym, atoms, sphhar, input%jspins, 0, vTot%mt, vTot0%mt)

        ! Construct the interstitial gradients.
        CALL pw_gradient(input%jspins, nG, recG, cell%bmat, rho0%pw, grRho0%pw)
        !CALL pw_gradient(input%jspins, nG, recG, cell%bmat, vTot0%pw, grVTot0%pw)

        ! Construct the muffin tin gradients.
        !CALL mt_gradient(input%jspins, atoms, juPhon%jplPlus, rho0%mt, grRho0%mt)
        !CALL mt_gradient(input%jspins, atoms, juPhon%jplPlus, vTot0%mt, grVTot0%mt)
        call mt_gradient_old(atoms, sphhar, sym, sphhar%clnu, sphhar%nmem, sphhar%mlh, rho%mt(:, :, :, 1), grRho0%mt(:, :, :, 1, 1, :) )

        CALL grRho%copyPotDen(rho)
        CALL grRho%resetPotDen()

        DO zInd = -stars%mx3, stars%mx3
           DO yInd = -stars%mx2, stars%mx2
              DO xInd = -stars%mx1, stars%mx1
                 iStar = stars%ig(xInd, yInd, zInd)
                 IF (iStar.EQ.0) CYCLE
                 !IF (stars%sk3(iStar).GT.stars%gmax) CYCLE
                 grRho%pw(iStar,:) = rho%pw(iStar,:) * cmplx(0.0,dot_product([1.0,0.0,0.0],matmul(real([xInd,yInd,zInd]),cell%bmat)))
              END DO
           END DO
        END DO

        ! Gradient of external unperturbed potential without terms canceling in the Sternheimer SCC with the linear variation of the
        ! external potential; see documentation of GenGrVeff0 for details. Required for 1st Sternheimer SCC iteration.
        harSw = .false.
        extSw = .true.
        xcSw = .false.
        testGoldstein = .false.
        grRhoTermSw = .TRUE.
        call GenGrVeff0( atoms, cell, stars, ngdp, harSw, extSw, xcSw, recG, rho0%pw(:, :, 1 ,1), rho0%mt(:, :, :, :, 1 ,1), &
                       & grRho0%pw(:, 1, 1 ,:), grRho0%mt(:, :, :, 1, 1, :), gausWts, ylm, dKernMTGPts, grVxcIRKern, &
                       & testGoldstein, grRhoTermSw, grVeff0IRdummy, grVeff0MT_init )
        ! Only need MT part for Sternheimer SCC
        deallocate ( grVeff0IRdummy )

        ! Add full muffin-tin gradient of xc and Hartree potential to muffin-tin gradient contribution of the external potential with
        ! volume term that should cancel in the Sternheimer SCC with the linear variation of the Hartree and the xc potential; see
        ! documentation of GenGrVeff0 for details. Required for 1st Sternheimer SCC iteration.
        harSw = .true.
        extSw = .false.
        xcSw = .true.
        testGoldstein = .false.
        grRhoTermSw = .true.
        call GenGrVeff0( atoms, cell, stars, ngdp, harSw, extSw, xcSw, recG, rho0%pw(:, :, 1 ,1), rho0%mt(:, :, :, :, 1 ,1), &
                       & grRho0%pw(:, 1, 1 ,:), grRho0%mt(:, :, :, 1, 1, :), gausWts, ylm, dKernMTGPts, grVxcIRKern, &
                       & testGoldstein, grRhoTermSw, grVeff0IRdummy, grVeff0MTAdd_init )
        ! Only need MT part for Sternheimer SCC
        deallocate ( grVeff0IRdummy )

        !if (compPhon.or.anfix) then
        !  iatom = 0
        !  do itype = 1, atoms%ntype
        !    do ieqat = 1, atoms%neq(itype)
        !      iatom = iatom + 1
        !      do idir = 1, 3
        !        do oqn_l = 0, atoms%lmax(iatom)
        !          do mqn_m = -oqn_l, oqn_l
        !            do imesh = 1, atoms%jri(iatom)
        !              grVeff0MT_init(imesh, lm, idir, iatom) = grVeff0MT_init(imesh, lm, idir, iatom) &
        !                                                                                   & + grVeff0MTAdd_init(imesh, lm, idir, iatom)
        !            end do ! imesh
        !          end do ! mqn_m
        !        end do ! oqn_l
        !      end do ! idir
        !    end do ! iDeqat
        !  end do ! iDtype
        !end if
        deallocate(grVeff0MTAdd_init)

        ! Full gradient of external unperturbed potential required for the Hellmann-Feynman contributions of the dynamical matrix.
        harSw = .false.
        extSw = .true.
        xcSw = .false.
        testGoldstein = .false.
        grRhoTermSw = .true.
        call GenGrVeff0( atoms, cell, stars, ngdp, harSw, extSw, xcSw, recG, rho0%pw(:, :, 1 ,1), rho0%mt(:, :, :, :, 1 ,1), &
                       & grRho0%pw(:, 1, 1 ,:), grRho0%mt(:, :, :, 1, 1, :), gausWts, ylm, dKernMTGPts, grVxcIRKern, &
                       & testGoldstein, grRhoTermSw, grVext0IR_DM, grVext0MT_DM )

        ! Gradient of effective unperturbed potential without terms canceling in the SternheimerSCC with the linear variation of the
        ! effective potential; see documentation of GenGrVeff0 for details. Required for regular iterations and the final Sternheimer
        ! SCC iteration.
        harSw = .true.
        extSw = .true.
        xcSw = .true.
        grRhoTermSw = .false.
        testGoldstein = .false.
        call GenGrVeff0( atoms, cell, stars, ngdp, harSw, extSw, xcSw, recG, rho0%pw(:, :, 1 ,1), rho0%mt(:, :, :, :, 1 ,1), &
                       & grRho0%pw(:, 1, 1 ,:), grRho0%mt(:, :, :, 1, 1, :), gausWts, ylm, dKernMTGPts, grVxcIRKern, &
                       & testGoldstein, grRhoTermSw, grVeff0IRdummy, grVeff0MT_main )
        ! Only need MT part for Sternheimer SCC
        deallocate ( grVeff0IRdummy )

        ! Gradient of Coulomb unperturbed potential without terms canceling in the SternheimerSCC with the linear variation of the
        ! Coulomb potential. Required for dynamical matrix surface term
        harSw = .true.
        extSw = .true.
        xcSw = .false.
        grRhoTermSw = .false.
        testGoldstein = .false.
        call GenGrVeff0( atoms, cell, stars, ngdp, harSw, extSw, xcSw, recG, rho0%pw(:, :, 1 ,1), rho0%mt(:, :, :, :, 1 ,1), &
                       & grRho0%pw(:, 1, 1 ,:), grRho0%mt(:, :, :, 1, 1, :), gausWts, ylm, dKernMTGPts, grVxcIRKern, &
                       & testGoldstein, grRhoTermSw, grVCoul0IR_DM_SF, grVCoul0MT_DM_SF )
        ! Only need MT part for Sternheimer SCC

        ! Full gradient of effective unperturbed potential required for the Hellmann-Feynman contributions of the dynamical matrix.
        harSw = .true.
        extSw = .true.
        xcSw = .true.
        grRhoTermSw = .true.
        testGoldstein = .false.

        call GenGrVeff0( atoms, cell, stars, ngdp, harSw, extSw, xcSw, recG, rho0%pw(:, :, 1 ,1), rho0%mt(:, :, :, :, 1 ,1), &
                       & grRho0%pw(:, 1, 1 ,:), grRho0%mt(:, :, :, 1, 1, :), gausWts, ylm, dKernMTGPts, grVxcIRKern, &
                       & testGoldstein, grRhoTermSw, grVeff0IR_DM, grVeff0MT_DM )

        call tlmplm4H0( atoms, enpara, usdus, input, tdHS0, 1, rbas1, rbas2, usdus%uuilon(:, :, 1), usdus%duilon(:, :, 1), usdus%ulouilopn(:, :, :, 1), ilo2p, vTot0%mt(:, :, :, 1, 1, 1), loosetdout )

    END SUBROUTINE dfpt_init

    SUBROUTINE stars_to_pw(sym, stars, jspins, nG, rhostar, rhopw)

        TYPE(t_sym),   INTENT(IN)  :: sym
        TYPE(t_stars), INTENT(IN)  :: stars
        INTEGER,       INTENT(IN)  :: jspins
        INTEGER,       INTENT(IN)  :: nG

        COMPLEX,       INTENT(IN)  :: rhostar(:, :)
        COMPLEX,       INTENT(OUT) :: rhopw(:, :, :, :)

        INTEGER                    :: iStar, iG, iGx, iGy, iGz, iSpin

        rhopw = CMPLX(0.0,0.0)

        DO iSpin = 1, jspins
            iG = 0
            DO iGz= -stars%mx3, stars%mx3
                DO iGy = -stars%mx2, stars%mx2
                    DO iGx = -stars%mx1, stars%mx1
                        iStar = stars%ig(iGx, iGy, iGz)
                        IF (iStar.NE.0) THEN
                            iG = iG + 1
                            rhopw(iG, iSpin, 1, 1) = rhostar(iStar, iSpin) &
                                    ! & * stars%rgphs(iGx, iGy, iGz) & TODO
                                    & * stars%nstr(iStar) / sym%nop
                        END IF
                    END DO
                END DO
            END DO
        END DO

    END SUBROUTINE stars_to_pw

    !---------------------------------------------------------------------------------------------------------------------------------
    !> @author
    !> Christian-Roman Gerhorst, Forschungszentrum Jülich: IAS1 / PGI1
    !>
    !> @brief
    !> Converts expansion coefficient of star expansion into expansion coefficients of plane-wave coefficients for interstitial
    !> unperturbed density or potential.
    !>
    !> @details
    !> See als Equation 7.4 (dissertation CRG)
    !>
    !> @param[in]  stars       : Stars type, see types.f90
    !> @param[in]  ngdp        : Number of G-vectors for potentials and densities
    !> @param[in]  stExpandQ   : Expansion coefficients of quantity expanded in stars
    !> @param[in]  gdp         : G-vectors of potentials and densities
    !> @param[out] gVecExpandQ : Expansion coefficients of quantity expanded in plane waves
    !---------------------------------------------------------------------------------------------------------------------------------
    subroutine convertStar2G(stExpandQ, gVecExpandQ, stars, ngdp, gdp, cell)

      ! Scalar Arguments
      type(t_stars),  intent(in)   :: stars
      integer,        intent(in)   :: ngdp

      ! Array Arguments
      complex,        intent(in)   :: stExpandQ(:)
      integer,        intent(in)   :: gdp(:, :)
      complex,        intent(out)  :: gVecExpandQ(ngdp)
      type(t_cell), optional, intent(in) :: cell

      ! Local Scalar Variables
      integer                      :: iGvec
      real                         :: Gext(3)
      gVecExpandQ = 0
      do iGvec = 1, ngdp
        gVecExpandQ(iGvec) = stExpandQ(stars%ig(gdp(1, iGvec), gdp(2, iGvec), gdp(3, iGvec))) * stars%rgphs(gdp(1, iGvec), gdp(2, iGvec), gdp(3, iGvec))
      end do
    end subroutine convertStar2G

    SUBROUTINE lh_to_sh(sym, atoms, sphhar, jspins, radfact, rholh, rhosh)

        TYPE(t_sym),    INTENT(IN)  :: sym
        TYPE(t_atoms),  INTENT(IN)  :: atoms
        TYPE(t_sphhar), INTENT(IN)  :: sphhar
        INTEGER,        INTENT(IN)  :: jspins, radfact
        REAL,           INTENT(IN)  :: rholh(:, 0:, :, :)
        COMPLEX,        INTENT(OUT) :: rhosh(:, :, :, :, :, :)

        INTEGER :: iSpin, iType, iEqat, iAtom, ilh, iMem, ilm, iR
        INTEGER :: ptsym, l, m
        REAL    :: factor

        rhosh = CMPLX(0.0,0.0)

        DO iSpin = 1, jspins
            iAtom = 0
            DO iType = 1, atoms%ntype
                DO iEqat = 1, atoms%neq(iType)
                    iAtom = iAtom + 1
                    ptsym = sym%ntypsy(iAtom)
                    DO ilh = 0, sphhar%nlh(ptsym)
                        l = sphhar%llh(iLH, ptsym)
                        DO iMem = 1, sphhar%nmem(ilh, ptsym)
                            m = sphhar%mlh(iMem, ilh, ptsym)
                            ilm = l * (l+1) + m + 1
                            DO iR = 1, atoms%jri(iType)
                                IF ((radfact.EQ.0).AND.(l.EQ.0)) THEN
                                    factor = atoms%rmsh(iR, iType) / sfp_const
                                ELSE IF (radfact.EQ.2) THEN
                                    factor = atoms%rmsh(iR, iType)**2
                                ELSE
                                    factor = 1.0
                                END IF
                                rhosh(iR, ilm, iAtom, iSpin, 1, 1) = &
                              & rhosh(iR, ilm, iAtom, iSpin, 1, 1) + &
                              & rholh(iR, ilh, iType, iSpin) * sphhar%clnu(iMem, ilh, ptsym) / factor
                            END DO
                        END DO
                    END DO
                END DO
            END DO
        END DO

    END SUBROUTINE lh_to_sh

    SUBROUTINE sh_to_lh(sym, atoms, sphhar, jspins, radfact, rhosh, rholhreal, rholhimag)

        ! WARNING: This routine will not fold back correctly for activated sym-
        !          metry and gradients (rho in l=0 and lattice harmonics do not
        !          allow l=1 --> gradient in l=1 is lost)

        TYPE(t_sym),    INTENT(IN)  :: sym
        TYPE(t_atoms),  INTENT(IN)  :: atoms
        TYPE(t_sphhar), INTENT(IN)  :: sphhar
        INTEGER,        INTENT(IN)  :: jspins, radfact
        COMPLEX,        INTENT(IN)  :: rhosh(:, :, :, :)
        REAL,           INTENT(OUT) :: rholhreal(:, 0:, :, :), rholhimag(:, 0:, :, :)

        INTEGER :: iSpin, iType, iEqat, iAtom, ilh, iMem, ilm, iR
        INTEGER :: ptsym, l, m
        REAL    :: factor

        rholhreal = 0.0
        rholhimag = 0.0

        DO iSpin = 1, jspins
            DO iType = 1, atoms%ntype
                iAtom = SUM(atoms%neq(:iType-1)) + 1
                ptsym = sym%ntypsy(iAtom)
                DO ilh = 0, sphhar%nlh(ptsym)
                    l = sphhar%llh(iLH, ptsym)
                    DO iMem = 1, sphhar%nmem(ilh, ptsym)
                        m = sphhar%mlh(iMem, ilh, ptsym)
                        ilm = l * (l+1) + m + 1
                        DO iR = 1, atoms%jri(iType)
                           IF ((radfact.EQ.0).AND.(l.EQ.0)) THEN
                               factor = atoms%rmsh(iR, iType) / sfp_const
                           ELSE IF (radfact.EQ.2) THEN
                               factor = atoms%rmsh(iR, iType)**2
                           ELSE
                               factor = 1.0
                           END IF
                            rholhreal(iR, ilh, iType, iSpin) = &
                          & rholhreal(iR, ilh, iType, iSpin) + &
                          &  real(rhosh(iR, ilm, iatom, iSpin) * conjg(sphhar%clnu(iMem, ilh, ptsym))) * factor
                            rholhimag(iR, ilh, iType, iSpin) = &
                          & rholhimag(iR, ilh, iType, iSpin) + &
                          & aimag(rhosh(iR, ilm, iatom, iSpin) * conjg(sphhar%clnu(iMem, ilh, ptsym))) * factor
                        END DO
                    END DO
                END DO
            END DO
        END DO

    END SUBROUTINE sh_to_lh

    SUBROUTINE pw_gradient(jspins, nG, recG, bmat, rhopw, grRhopw)

        INTEGER,       INTENT(IN)  :: jspins
        INTEGER,       INTENT(IN)  :: nG
        INTEGER,       INTENT(IN)  :: recG(:, :)
        REAL,          INTENT(IN)  :: bmat(:, :)

        COMPLEX,       INTENT(IN)  :: rhopw(:, :, :, :)
        COMPLEX,       INTENT(OUT) :: grRhopw(:, :, :, :)

        INTEGER :: iG, iSpin
        REAL    :: gExt(3)

        grRhopw = CMPLX(0.0,0.0)

        DO iSpin = 1, jspins
            DO iG = 1, nG
                gExt = matmul( bmat, recG(:, iG) )
                grRhopw(iG, iSpin, 1, :) = ImagUnit * gExt * rhopw(iG, iSpin, 1, 1)
            END DO
        END DO

    END SUBROUTINE pw_gradient

    SUBROUTINE mt_gradient(jspins, atoms, lplus, rhosh, grRhosh)

        TYPE(t_atoms), INTENT(IN)  :: atoms
        INTEGER,       INTENT(IN)  :: jspins, lplus
        COMPLEX,       INTENT(IN)  :: rhosh(:, :, :, :, :, :)
        COMPLEX,       INTENT(OUT) :: grRhosh(:, :, :, :, :, :)

        INTEGER :: lmaxgrad(atoms%ntype)

        lmaxgrad = atoms%lmax + lplus

    END SUBROUTINE

    !---------------------------------------------------------------------------------------------------------------------------------
    !> @author
    !> Christian-Roman Gerhorst, Forschungszentrum Jülich: IAS1 / PGI1
    !>
    !> @brief
    !> Calculates the spherical harmonic expansion coefficients of the muffin-tin gradient applied to an arbitrary function multiplied
    !> by $r^2$. The resulting gradient expansion coefficients are multiplied by a factor $r^2$.
    !>
    !> @note
    !> The ingoing function is assumed to be multiplied with $r$. The outgoing resulting function is also multiplied by $r$.
    !>
    !> @param[in]  atoms     : Contains atoms-related quantities; definition of its members in types.F90 file.
    !> @param[in]  lathar    : Contains entities concerning the lattice harmonics; more precise definition in type.F90 file.
    !> @param[in]  clnu_atom : Coefficients to transform from lattice harmonics to spherical harmonics without any symmetry.
    !> @param[in]  nmem_atom : Number of member per lattice harmonic in a system without any symmetry.
    !> @param[in]  mlh_atom  : Magnetic quantum numbers of the members of any lattice harmonic in a system without any symmetry.
    !> @param[in]  r2FlhMt   : Lattice harmonic coefficients of muffin-tin quantity multiplied by a factor of r**2.
    !> @param[out] r2GrFshMt : Spherical harmonic coefficients of muffin-tin quantity's gradient multiplied by a factor of r**2
    !---------------------------------------------------------------------------------------------------------------------------------
    subroutine mt_gradient_old(atoms, lathar, sym, clnu_atom, nmem_atom, mlh_atom, r2FlhMt, GrFshMt)

      use m_gaunt, only : gaunt1

      ! Type parameters
      ! ***************
      type(t_atoms),               intent(in)  :: atoms
      type(t_sphhar),              intent(in)  :: lathar
      type(t_sym),                 intent(in)  :: sym

      ! Array parameters
      ! ****************
      complex,                     intent(in)  :: clnu_atom(:, 0:, :)
      integer,                     intent(in)  :: nmem_atom(0:, :)
      integer,                     intent(in)  :: mlh_atom(:, 0:, :)
      real,                        intent(in)  :: r2FlhMt(:, 0:, :)
      complex,                     intent(out) :: GrFshMt(:, :, :, :)


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
      allocate( r2GrFshMtNat(atoms%jmtd, ( atoms%lmaxd + 1)**2, atoms%nat, 3) )
      allocate( rDerFlhMt(atoms%jmtd) )
      GrFshMt = cmplx(0., 0.)
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
                grFshMt(imesh, lm, iatom, 1:3) = matmul( Tmatrix0(1:3, 1:3), r2GrFshMtNat(imesh, lm, iatom, 1:3) ) / atoms%rmsh(imesh, itype)**2
              end do
            end do ! mqn_m
          end do ! oqn_l
        end do ! ieqat
      end do ! itype

    end subroutine mt_gradient_old

    subroutine derivative(f, itype, atoms, df)

       integer,       intent(in)  :: itype
       type(t_atoms), intent(in)  :: atoms
       real,          intent(in)  :: f(atoms%jri(itype))
       real,          intent(out) :: df(atoms%jri(itype))
       real                       :: h, r, d21, d32, d43, d31, d42, d41, df1, df2, s
       real                       :: y0, y1, y2
       integer                    :: i, n

       n = atoms%jri(itype)
       h = atoms%dx(itype)
       r = atoms%rmsh(1, itype)

       ! use Lagrange interpolation of 3rd order (and averaging) for points 3 to n
       d21 = r * (exp(h)-1) ; d32 = d21 * exp(h) ; d43 = d32 * exp(h)
       d31 = d21 + d32      ; d42 = d32 + d43
       d41 = d31 + d43
       df(1) =   d31*d41 / (d21*d32*d42) * f(2) + ( -1d0/d21 - 1d0/d31 - 1d0/d41) * f(1)&
      &        - d21*d41 / (d31*d32*d43) * f(3) + d21*d31 / (d41*d42*d43) * f(4)
       df(2) = - d32*d42 / (d21*d31*d41) * f(1) + (  1d0/d21 - 1d0/d32 - 1d0/d42) * f(2)&
      &        + d21*d42 / (d31*d32*d43) * f(3) - d21*d32 / (d41*d42*d43) * f(4)
       df1   =   d32*d43 / (d21*d31*d41) * f(1) - d31*d43 / (d21*d32*d42) * f(2) +&
      &  ( 1d0/d31 + 1d0/d32 - 1d0/d43 ) * f(3) + d31*d32 / (d41*d42*d43) * f(4)
       do i = 3, n - 2
          d21 = d32 ; d32 = d43 ; d43 = d43 * exp(h)
          d31 = d42 ; d42 = d42 * exp(h)
          d41 = d41 * exp(h)
          df2   = - d32*d42 / (d21*d31*d41) * f(i-1) + ( 1d0/d21 - 1d0/d32 - 1d0/d42) * f(i) + &
      &             d21*d42 / (d31*d32*d43) * f(i+1) - d21*d32 / (d41*d42*d43) * f(i+2)
          df(i) = ( df1 + df2 ) / 2
          df1   = d32*d43 / (d21*d31*d41) * f(i-1) - d31*d43 / (d21*d32*d42) * f(i) +&
      &    ( 1d0/d31 + 1d0/d32 - 1d0/d43 ) * f(i+1) + d31*d32 / (d41*d42*d43) * f(i+2)
       enddo
       df(n-1) = df1
       df(n)   = - d42*d43 / (d21*d31*d41) * f(n-3) + d41*d43 / (d21*d32*d42) * f(n-2) -&
      &            d41*d42 / (d31*d32*d43) * f(n-1) + ( 1d0/d41 + 1d0/d42 + 1d0/d43 ) * f(n)
       ! for first two points use Lagrange interpolation of second order for log(f(i))
       ! or, as a fall-back, Lagrange interpolation with the conditions f(1), f(2), f(3), f'(3).
       s = sign(1d0,f(1))
       if(sign(1d0,f(2)) /= s .or. sign(1d0,f(3))  /= s .or. any(abs(f(:3)) < 1e0)) then
          d21   = r * (exp(h)-1)
          d32   = d21 * exp(h)
          d31   = d21 + d32
          s     = df(3) / (d31*d32) - f(1) / (d21*d31**2) + f(2) / (d21*d32**2) - f(3) / (d31**2*d32) - f(3) / (d31*d32**2)
          df(1) = - (d21+d31) / (d21*d31) * f(1) + d31 / (d21*d32) * f(2) - d21 / (d31*d32) * f(3) + d21*d31 * s

          df(2) = - (d21-d32) / (d21*d32) * f(2) - d32 / (d21*d31) * f(1) + d21 / (d31*d32) * f(3) - d21*d32 * s
       else
          y0    = log(abs(f(1)))
          y1    = log(abs(f(2)))
          y2    = log(abs(f(3)))
          df(1) = ( - 3*y0/2 + 2*y1 - y2/2 ) * f(1) / (h*r)
          df(2) = (y2-y0)/2                  * f(2) / (h*r*exp(h))
       endif
    end subroutine derivative

    !---------------------------------------------------------------------------------------------------------------------------------
    !> @author
    !> Christian-Roman Gerhorst, Forschungszentrum Jülich: IAS1 / PGI1
    !>
    !> @brief
    !> This subroutine generates the G-vectors which are used for constructing the potential and the density.
    !>
    !> @details
    !> They are fulfilling |G| < Gmax
    !>
    !> @param[in]  starsT   : Contains stars-related quantities; definition of its members in types.F90 file.
    !> @param[in]  cellT    : Contains unit-cell related quantities; definition of its members in types.F90 file.
    !> @param[out] input    : Input type, see types.f90.
    !> @param[out] ngpd     : Number of G-vectors used for density and potential.
    !> @param[out] ngdp2km  : Number of G-vectors for potentials and densities which are smaller than 2 kmax.
    !> @param[out] gdp      : Contains G-vectors in internal coordinates used for density and potential, see for dimensions.
    !> @param[out] gdp2Ind  : Stores the index of a potential and density G-Vector. The dimensions are the G-vector components.
    !> @param[out] gdp2iLim : Stores the min and maxvals of gdp2Ind
    !> @param[in]  testMode : Indexarray stores all G-vectors until Gmax instead of 2kmax if set .true.
    !---------------------------------------------------------------------------------------------------------------------------------
    subroutine genPotDensGvecs(starsT, cellT, inputT, ngdp, ngdp2km, gdp, gdp2Ind, gdp2iLim, testMode )

      use m_juDFT

      ! Variables:
      ! Gx      : loop variable to run over all possible x-components of an internal G-vector
      ! Gy      : loop variable to run over all possible y-components of an internal G-vector
      ! Gz      : loop variable to run over all possible z-components of an internal G-vector
      ! gdptemp : auxiliary array for intermediate storage of accepted G-vectors
      ! Gint    : temporary variable to store current candidate for a G-vector in internal coordinates
      ! Gext    : temporary variable to store current candidate for a G-vector in external coordinates

      ! Scalar parameters
      type(t_stars),          intent(in)   :: starsT
      type(t_cell),           intent(in)   :: cellT
      type(t_input),          intent(in)   :: inputT

      integer,                intent(out)  :: ngdp
      integer,                intent(out)  :: ngdp2km
      logical,                intent(in)   :: testMode

      ! Array parameter
      integer,  allocatable,  intent(out)  :: gdp(:, :)
      integer,  allocatable,  intent(out)  :: gdp2Ind(:, :, :)
      integer,                intent(out)  :: gdp2iLim(2, 3)

      ! Local scalar variables
      integer                              :: Gx, Gy, Gz, iG
      integer                              :: ngrest

      ! Local array variables
      integer                              :: gdptemp2kmax(3, (2 * starsT%mx1 + 1) * (2 * starsT%mx2 + 1) * (2 * starsT%mx3 +  1))
      integer                              :: gdptemprest(3, (2 * starsT%mx1 + 1) * (2 * starsT%mx2 + 1) * (2 * starsT%mx3 +  1))
      integer                              :: Gint(3)
      real                                 :: Gext(3)

      ! From all possible G-vectors in a box, only these are accepted which are element of a sphere with radius gmax. As benchmark,
      ! it is checked, whether the current G-vector candidate is contained within any star
      ngdp = 0
      gdptemp2kmax = 0
      gdptemprest = 0
      ngdp2km = 0
      ngrest = 0
      ! The idea to sort the G-vectors until 2kmax before all other Gvectors stems from M. Betzinger
      do Gx = -starsT%mx1, starsT%mx1
        do Gy = -starsT%mx2, starsT%mx2
          do Gz = -starsT%mx3, starsT%mx3
            Gint = [Gx, Gy, Gz]
            Gext =  matmul(cellT%bmat, Gint) !transform from internal to external coordinates
            if (norm2(Gext) <= inputT%gmax) then
!  #ifdef DEBUG_MODE
!              if (starsT%ig(Gx, Gy, Gz) <= 0) then
!                call juDFT_error('Inconsistency in determination of G-vectors for potential or density', calledby='genPotDensGvecs', &
!                   & hint='Check whether selection methods correct.', file='jpPotDens.F90', line=52)
!              endif
!  #endif
              ngdp = ngdp + 1
              ! Sort G-vectors
              if ( norm2(Gext) <= 2 * inputT%rkmax ) then
                ngdp2km = ngdp2km + 1
                gdptemp2kmax(:, ngdp2km) = Gint(:)
              else
                ngrest = ngrest + 1
                gdptemprest(:, ngrest) = Gint(:)
              end if
            endif
          enddo !Gz
        enddo !Gy
      enddo !Gx
      allocate(gdp(3, ngdp))
      gdp(:, :ngdp2km) = gdptemp2kmax(:, :ngdp2km)
      gdp(:, ngdp2km + 1 : ngdp) = gdptemprest(:, :ngrest)


      ! Create mapping array from G-vector to G-vector index up to Gmax
      gdp2iLim(1, 1) = minval(gdp(1, :))
      gdp2iLim(2, 1) = maxval(gdp(1, :))
      gdp2iLim(1, 2) = minval(gdp(2, :))
      gdp2iLim(2, 2) = maxval(gdp(2, :))
      gdp2iLim(1, 3) = minval(gdp(3, :))
      gdp2iLim(2, 3) = maxval(gdp(3, :))
      allocate(gdp2Ind(gdp2iLim( 1, 1) : gdp2iLim( 2, 1), gdp2iLim( 1, 2) : gdp2iLim( 2, 2), gdp2iLim( 1, 3) : gdp2iLim( 2, 3)))
      gdp2Ind = 0
      if (testMode) then
        do iG = 1, ngdp2km
          gdp2Ind(gdptemp2kmax(1, iG), gdptemp2kmax(2, iG), gdptemp2kmax(3, iG)) = iG
        end do
      else
        do iG = 1, ngdp
          gdp2Ind(gdptemp2kmax(1, iG), gdptemp2kmax(2, iG), gdptemp2kmax(3, iG)) = iG
        end do
      end if

    end subroutine genPotDensGvecs

    !---------------------------------------------------------------------------------------------------------------------------------
    !> @author
    !> Christian-Roman Gerhorst, Forschungszentrum Jülich: IAS1 / PGI1
    !>
    !> @brief
    !> Calculating the plane-wave expansion coefficients of the functional derivative of any interstitial kernel with respect to the
    !> unperturbed density.
    !>
    !> @note
    !> This routine is tested as we compared Fleur and juPhon potentials
    !> @note
    !> At the moment, we use the x-alpha potential, but in principle in this routine the call to libxc can be integrated
    !> @note
    !> This routine is similiar to the routine that calculates the interstitial x-alpha kernel in Fleur.
    !>
    !> @param[in] stars : Stars type, see types.f90
    !> @param[in] gdp   : G-vectors of potentials and densities
    !> @param[in] ngdp  : Number of G-vectors for potentials and densities
    !> @param[in] qpw   : Star coefficients of the unperturbed and converged interstitial density parsed from Fleur
    !> @param[in] fg3G  : Plane wave interstitial coefficients of the functional derivative of the xc-kernel with respect to the
    !>                    density
    !---------------------------------------------------------------------------------------------------------------------------------
    subroutine calcIRdVxcKern(stars, gdp, ngdp, qpw, fg3G, oldmode, xcpot)

      use m_fft3d

      ! Type parameter
      type(t_stars),              intent(in) :: stars

      ! Scalar parameter
      integer,                    intent(in) :: ngdp

      ! Array parameter
      integer,                    intent(in) :: gdp(:, :)
      complex,                    intent(in) :: qpw(:)
      complex,       allocatable, intent(out):: fg3G(:)

      ! New input:
      LOGICAL, INTENT(IN)        :: oldmode
      CLASS(t_xcpot), INTENT(IN) :: xcpot

      ! Scalar variables
      integer                                :: ifftd
      integer                                :: iG

      ! Array variables
      real,          allocatable             :: af3(:), bf3(:)
      real,          allocatable             :: VxcIRKern(:)
      complex,       allocatable             :: fg3(:)

      ! Init and allocate output variable
      allocate(fg3G(ngdp))
      fg3G(:) = cmplx(0., 0.)

      ! Size of FFT mesh
      ifftd = 27 * stars%mx1 * stars%mx2 * stars%mx3
      allocate( af3(0:ifftd - 1), bf3(0:ifftd - 1), VxcIRKern(ifftd), fg3(stars%ng3)) ! TODO: Was n3d originally! Correct?
      af3(:) = 0.
      bf3(:) = 0.
      VxcIRKern(:) = 0.
      fg3(:) = cmplx(0., 0.)

      ! FFT of qpw from reciprocal to direct space
      ! We can use the fft which is desiged for star expanded quantities because here we have the qpw which is given in stars.
      call fft3d(af3, bf3, qpw, stars, +1)

      if ( any( abs(bf3) >= 1e-8 ) ) then
        write(*, *) 'Warning: FFT in calcIRdVxcKern has complex contributions!'
      end if
      bf3 = 0

      ! Calculate functional derivative of kernel based on qpw in direct space because we know the representation in direct space.
      call calcKernDerOnGrid(ifftd, 1., af3, VxcIRKern, oldmode, xcpot)

      ! Back-FFT to direct space. We only have made a functional derivative which does not break the symmetry so we can still use this
      ! kind of fft. We use the fft of the new fleur version, if we have no optional parameter the function is the same as in old
      ! fleur.
      CALL fft3d(VxcIRKern, bf3, fg3, stars, -1)

      ! Transform from star representation to plane wave representation
      do iG = 1, ngdp
        fg3G(iG) = fg3G(iG) + fg3(stars%ig(gdp(1, iG), gdp(2, iG), gdp(3, iG))) * stars%rgphs(gdp(1, iG), gdp(2, iG), gdp(3, iG))
      end do

    end subroutine calcIRdVxcKern

    !---------------------------------------------------------------------------------------------------------------------------------
    !> @author
    !> Christian-Roman Gerhorst, Forschungszentrum Jülich: IAS1 / PGI1
    !>
    !> @brief
    !> Calculating the functional derivative of any MT kernel with respect to the unperturbed density
    !>
    !> @param[in] atoms       : Atoms type, see types.f90
    !> @param[in] dimens      : Dimension type, see types.f90
    !> @param[in] sphhar      : Lattice harmonics type, see types.f90.
    !> @param[in] rho0MT      : Lattice harmonic coefficients of unperturbed and converged muffin-tin density parsed from Fleur
    !> @param[in] nmem_atom   : Number of lattice harmonic members for every atom.
    !> @param[in] clnu_atom   : Phase mediating between stars and plane waves.
    !> @param[in] mlh_atom    : Magnetic quantum number m of lattice harmonic members for every atom.
    !> @param[in] gWghts      : Gauss weights for Gauss quadrature
    !> @param[in] ylm         : Set of spherical harmonics whose arguments are the unit vectors of the Gauss mesh points up to lmax
    !> @param[in] dKernMTGPts : Spherical harmonic muffin-tin coefficients of the functional derivative of the xc-kernel with
    !>                         respect to the density
    !---------------------------------------------------------------------------------------------------------------------------------
    subroutine calcMTdVxcKern(atoms, sphhar, sym, rho0MT, nmem_atom, clnu_atom, mlh_atom, gWghts, ylm, dKernMTGPts, oldmode, xcpot)

      use m_gaussp
      use m_ylm_old

      ! Type parameters
      type(t_atoms),                  intent(in)  :: atoms
      type(t_sphhar),                 intent(in)  :: sphhar
      type(t_sym),                    intent(in)  :: sym

      ! Array parameters
      real,                           intent(in)  :: rho0MT(:, 0:, :)
      integer,                        intent(in)  :: nmem_atom(0:, :)
      complex,                        intent(in)  :: clnu_atom(:, 0:, :)
      integer,                        intent(in)  :: mlh_atom(:, 0:, :)
      real,              allocatable, intent(out) :: gWghts(:) ! gaussian weights belonging to gausPts
      complex,           allocatable, intent(out) :: ylm(:, :)
      real,              allocatable, intent(out) :: dKernMTGPts(:, :, :)
      ! New input:
      LOGICAL, INTENT(IN)        :: oldmode
      CLASS(t_xcpot), INTENT(IN) :: xcpot

      ! Scalar local variables
      integer                                     :: igmesh
      integer                                     :: iatom
      integer                                     :: itype
      integer                                     :: ieqat
      integer                                     :: ptSym
      integer                                     :: irmesh
      integer                                     :: ilh
      integer                                     :: oqn_l
      integer                                     :: imem
      integer                                     :: mqn_m
      integer                                     :: lm

      ! Local Array Variables
      real,              allocatable              :: gPts(:, :) ! gaussian points to exactly integrate spherial harmonics
      complex,           allocatable              :: rhoMTGpts(:, :), ylmtemp(:)

      ! Replaced dimension%nspd by atoms%nsp()
      allocate( gWghts(atoms%nsp()), ylmtemp( ( atoms%lmaxd + 1 )**2), ylm(atoms%nsp(), ( atoms%lmaxd + 1 )**2), &
        & dKernMTGPts(atoms%nsp(), atoms%jmtd, atoms%nat), gPts(3, atoms%nsp()), rhoMTGpts(atoms%nsp(), atoms%jmtd) )

      gWghts(:) = 0.
      ylm(:, :) = cmplx(0., 0.)
      ylmtemp = cmplx(0., 0.)
      dKernMTGPts(:, :, :) = 0.
      gPts(:, :) = 0.
      rhoMTGpts(:, :) = cmplx(0., 0.)

      ! generates dimension%nspd points on a spherical shell with radious 1.0. The angular mesh is equidistant in phi, theta are the
      ! zeros of the legendre polynomials.?
      ! Gaussian points to exactly integrate product of spherical harmonics up to lmax
      ! todo understand the gaussp routine ! Gauss quadrature on wikipedia!
      call gaussp( atoms%lmaxd, gPts, gWghts )

      ! The ylms are note filled linearily as a compromise to perform a linear run through all the other arrays
      call ylmnorm_init(atoms%lmaxd)
      do igmesh = 1, atoms%nsp()
        call ylm4( atoms%lmaxd, gPts(:, igmesh), ylmtemp )
        ylm(igmesh, :) = ylmtemp
      end do

      dKernMTGPts(:, :, :) = 0
      iatom = 0
      do itype = 1, atoms%ntype
        do ieqat = 1, atoms%neq(itype)
          iatom = iatom + 1
          ! Evaluate the unperturbed density on the Gauss mesh then calculate the functional derivative of the density
          ptSym = sym%ntypsy(iatom)
          rhoMTGpts(:, :) = cmplx(0., 0.)
          do ilh = 0, sphhar%nlh(ptsym)
            oqn_l = sphhar%llh(ilh, ptsym)
            do imem = 1, nmem_atom(ilh, iatom)
              mqn_m = mlh_atom(imem, ilh, iatom)
              lm = oqn_l * ( oqn_l + 1 ) + mqn_m + 1
              do irmesh = 1, atoms%jri(itype)
                do igmesh = 1, atoms%nsp()
                  rhoMTGpts(igmesh, irmesh) = rhoMTGpts(igmesh, irmesh) + rho0MT(irmesh, ilh, itype) * clnu_atom(imem, ilh, iatom) &
                                            & * ylm(igmesh, lm)
                end do ! igmesh
              end do ! irmesh
            end do ! imem
          end do ! ilh

          if ( any( abs(aimag( rhoMTGpts )) > 1e-7 ) ) then
            write(*, *) 'Warning rhoMTGpts has imaginary components.'
          end if
          ! Calculate functional derivative of the density on the Gauss mesh for every atom
          do irmesh = 1, atoms%jri(itype)
            call calcKernDerOnGrid(atoms%nsp(), 1., real(rhoMTGpts(:, irmesh)), dKernMTGPts(:, irmesh, iatom), oldmode, xcpot)
          end do ! irmesh
        end do ! ieqat
      end do ! itype

    end subroutine calcMTdVxcKern

    !---------------------------------------------------------------------------------------------------------------------------------
    !> @author
    !> Christian-Roman Gerhorst, Forschungszentrum Jülich: IAS1 / PGI1
    !>
    !> @brief
    !> Calculates the x-alpha functional derivative of the xc-kernel with respect to the unperturbed density in real space.
    !>
    !> @param[in] nGridPts        : Number of grid points of the real space mesh of rhoMTGpts
    !> @param[in] alpha           : Parameter alpha in the x-alpha kernel.
    !> @param[in] rhoMTGpts       : Unperturbed density on a real-space mesh
    !> @param[in] grVxcMTKernGPts : Functional derivative of the x-alpha kernel with respect to the density.
    !---------------------------------------------------------------------------------------------------------------------------------
    subroutine calcKernDerOnGrid(nGridPts, alpha, rhoMTGpts, grVxcMTKernGPts, oldmode, xcpot)

      USE m_types
      USE m_types_xcpot_libxc
      USE m_npy

      integer,           intent(in)  :: nGridPts
      real,              intent(in)  :: alpha
      real,              intent(in)  :: rhoMTGpts(:)
      real,              intent(out) :: grVxcMTKernGPts(:)
      ! New input needed for the expansion to generalized dVxc:
      LOGICAL, INTENT(IN)        :: oldmode
      CLASS(t_xcpot), INTENT(IN) :: xcpot

      real                           :: prfac
      integer                        :: imesh
      real                           :: rhoCapped
      LOGICAL :: l_libxc

      ! Libxc:
      real :: rhoMTGptsdummy(SIZE(rhoMTGpts), 1)
      real :: grVxcMTKernGPtsdummy(SIZE(grVxcMTKernGPts), 1)

      ! Oldmode for
      grVxcMTKernGPts(:) = 0.
      prfac = 1. / 9. / pi_const
      do imesh = 1, nGridPts
        rhoCapped = max( 1e-15, rhoMTGpts(imesh) )
        rhoMTGptsdummy(imesh, 1) = rhoCapped
        ! We merged all factors from FLEUR together converted it in hartree units (division by 2) and multiplied it with the 1 / 3
        ! from the derivative which is the prefactor used here so we have 1 / 3 * (3 / pi)**(1/3) * rho**(-2 / 3)
        grVxcMTKernGPts(imesh) = - prfac**(1. / 3.) * alpha * rhoCapped**(-2. / 3.)
      end do

      !CALL save_npy("fxc_inbuild.npy", grVxcMTKernGPts)

      SELECT TYPE(xcpot)
      TYPE IS (t_xcpot_libxc)
         l_libxc=.TRUE.
      END SELECT

#ifdef CPP_LIBXC
      IF (l_libxc) THEN
          CALL xcpot%get_fxc(1, rhoMTGptsdummy, grVxcMTKernGPtsdummy)
          !CALL save_npy("fxc_libxc.npy", grVxcMTKernGPtsdummy(:, 1))
          !STOP
      ELSE
          CALL judft_error("You tried to use Fleur DFPT with a non-libxc functional. Please fix that.")
      END IF
#else
      CALL judft_error("You compiled Fleur without libxc but want to use DFPT. Please fix that.")
#endif

    ! The .npy comparison showed this is ok.
    IF (.NOT.oldmode) THEN
        grVxcMTKernGPts = grVxcMTKernGPtsdummy(:, 1)
    END IF

    end subroutine calcKernDerOnGrid

    !---------------------------------------------------------------------------------------------------------------------------------
    !> @author
    !> Christian-Roman Gerhorst, Forschungszentrum Jülich: IAS1 / PGI1
    !>
    !> @brief
    !> Calculates the parameters for the gradient of the Gauss function pseudo-core density and the resulting interstitial core-tail
    !> correction for q = 0 that still have to be modulated with the respective q.
    !>
    !> @details
    !> See also 7.26b (dissertation CRG)
    !>
    !> @note
    !> This routine is based on a routine from FLEUR
    !>
    !>  @param[in]  atoms   : Atoms type, see types.f90
    !>  @param[in]  cell    : Unit cell type, see types.f90.
    !>  @param[in]  sym     : Symmetries type, see types.f90.
    !>  @param[in]  stars   : Stars type, see types.f90.
    !>  @param[in]  dimens  : Dimension type, see types.f90.
    !>  @param[in]  input   : Input type, see types.f90.
    !>  @param[in]  logUnit : Unit number for juPhon.log.
    !>  @param[in]  ngdp2km : Number of G-vectors for potentials and densities which are smaller than 2 kmax.
    !>  @param[in]  gdp     : G-vectors of potentials and densities
    !>  @param[out] acoff   : Uppercase A of Gauss curve pseudo core-density, see 7.26b (dissertation CRG)
    !>  @param[out] alpha   : Lowercase a of Gauss curve pseudo core-density, see 7.26b (dissertation CRG)
    !>  @param[out] qpwcG   : Plane-wave expansion coefficients of the FFT of the gradient of the pseudo core-density ( Gauss curve )
    !>
    !---------------------------------------------------------------------------------------------------------------------------------
    subroutine calcPsDensMT( fmpi, atoms, cell, sym, stars, input, ngdp2km, acoff, alpha, qpwcG, gdp )

      use m_constants
      use m_diflgr
      use m_types
      use m_juDFT_stop, only : juDFT_error, juDFT_warn
      USE m_cdn_io
      !
      !     .. Parameters ..
      TYPE(t_mpi),                    intent(in)            :: fmpi
      TYPE(t_atoms),                  intent(in)            :: atoms
      TYPE(t_cell),                   intent(in)            :: cell
      TYPE(t_sym),                    intent(in)            :: sym
      TYPE(t_stars),                  intent(in)            :: stars
      TYPE(t_input),                  intent(in)            :: input

      ! Scalar Parameters
      integer,                        intent(in)            :: ngdp2km

      ! Array Parameters
      integer,                        intent(in)            :: gdp(:, :)
      real,              allocatable, intent(out)           :: acoff(:)
      real,              allocatable, intent(out)           :: alpha(:)
      complex,           allocatable, intent(out)           :: qpwcG(:, :)

      ! Local Scalars
      real                                                  :: dif
      real                                                  :: dxx
      real                                                  :: tl_14
      integer                                               :: j
      integer                                               :: j1
      integer                                               :: ncmsh
      integer                                               :: ispin
      integer                                               :: itype
      integer                                               :: imesh
      integer,                                    parameter :: method1 = 2
      integer,                                    parameter :: method2 = 1
      LOGICAL :: l_CoreDenPresent

      ! Local Arrays
      integer,           allocatable                        :: mshc(:)
      real,              allocatable                        :: rat(:, :)
      !real,              allocatable                        :: rh(:, :)
      REAL                             :: rh(atoms%msh,atoms%ntype,input%jspins)
      REAL                             :: qint(atoms%ntype,input%jspins)
      REAL                             :: tec(atoms%ntype,input%jspins)

      !
      !     OLD VERSION:
      !     The previous version to calculate the overlap of the
      !     core tails was done in real space:
      !     A three dimensional real space grid was set up. On each
      !     of these grid points the charge density of all atoms inside
      !     the unit cell and neigboring unit cells was calculated.
      !     This calculation required a lagrange fit since the core
      !     charge is on a radial mesh. The same was done in the vacuum
      !     region, except on each z-plane we worked on a two dimension
      !     grid. After the charge density was generated on a equidistant
      !     grid the charge density was FFTed into G-space. The set up
      !     of the charge density on the real space grid is rather time
      !     consuming. 3-loops are required for the 3D grid
      !                1-loop over all atoms in the unit cells
      !                Larange interpolation
      !                3D and 2D FFTs
      !     In order to save time the non-spherical contributions inside
      !     the sphere had been ignored. It turns out that the later
      !     approximation is pure in the context of force calculations.
      !
      !     PRESENT VERSION:
      !     The present version is written from scratch. It is based on the
      !     idea that the FFT of an overlap of spherically symmetric
      !     charges can be expressed by the product of
      !
      !     sum_natype{ F(G,ntype) * sum_atom(atype) {S(\vec{G},atom)}}
      !
      !     of form factor F and structure factor S. The Form factor
      !     depends only G while the structure factor depends on \vec{G}
      !     and can build up in G-space. F of a gaussian chargedensity can
      !     be calculated analytically.

      !     The core-tails to the vacuum region are described by an
      !     exponentially decaying function into the vacuum:

      !     rho(r_||,z,ivac)= sum_n{ rho(n,ivac) * exp(-kappa*(z-z_v))
      !                                          * exp(iG(n)r_||) }

      !     And the plane waves are expanded into lattice harmonics
      !     up to a l_cutoff. Tests of the accuracy inside the sphere
      !     have shown that a reduction of the l_cutoff inside the
      !     in order to save time leads to sizable errors and should
      !     be omitted.

      !     rho_L(r) =  4 \pi i^l \sum_{g =|= 0}  \rho_int(g) r_i^{2} \times
      !                              j_{l} (gr_i) \exp{ig\xi_i} Y^*_{lm} (g)

      !     Tests have shown that the present version is about 10 times
      !     faster than the previous one also all nonspherical terms are
      !     included up to l=8 and the previous one included only l=0.
      !     For l=16 it is still faster by factor 2.

      !     coded                  Stefan Bl"ugel, IFF Nov. 1997
      !     tested                 RObert Abt    , IFF Dez. 1997
      !*****************************************************************
      !
      !     ..
      !----> Abbreviation
      !
      !     mshc     : maximal radial meshpoint for which the radial coretail
      !                density is larger than tol_14
      !     method1  : two different ways to calculate the derivative of the
      !                charge density at the sphere boundary.
      !                (1) use subroutine diflgr based on lagrange interpol.
      !                (2) use two point formular in real space,
      !                    see notes of SB.
      !                Tests have shown that (2) is more accurate.
      !     method2  : two different integration routines to calculate form
      !                factor of coretails outside the sphere.
      !                (1) use subroutine intgrz to integrate the tails from
      !                    outside to inside.
      !                (2) use subroutine intgr3 to integrate the tails from
      !                    muffin-tin radius to outside and include correction
      !                    for start up.
      !                Tests have shown that (1) is more accurate.
      !
      !

      allocate( acoff(atoms%ntype), alpha(atoms%ntype) )
      allocate( mshc(atoms%ntype), rat(atoms%msh, atoms%ntype) )
      !allocate( rh(dimens%msh, atoms%ntype))
      allocate( qpwcG(ngdp2km, atoms%nat) )

      mshc(:) = 0
      rat(:, :) = 0.
      acoff(:) = 0.
      alpha(:) = 0.
      rh = 0.
      tl_14 = 1e-14
      qpwcG = cmplx(0., 0.)

      ! Old juPhon read
      ! Read in core density from cdnc
      !call fopen(1000, name='cdnc', status='old', action='read', form='unformatted')
      !do ispin = 1, input%jspins
    !    if ( ispin == 1 ) then
    !      rewind( 1000 )
    !    end if
    !    do itype = 1, atoms%ntype

    !      ncmsh = nint( log( (atoms%rmt(itype) + 10.0) / atoms%rmsh(1, itype) ) / atoms%dx(itype) + 1 )
    !      ncmsh = min( ncmsh, dimens%msh )
          ! read in core density
          ! rh should be dependent on spin SPIN
    !      read( 1000 ) (rh( imesh, itype ), imesh=1, ncmsh)

          ! skip kinetic enrgy of the core
    !      read( 1000 )
    !    end do
        ! only important if more than one spin later
    !    read( 1000 )
     ! end do
      !call fclose(1000)

      !!! READ FROM Max_FLEUR
      l_CoreDenPresent = .FALSE.
      IF (input%kcrel==0) THEN
         ! Generate input file ecore for subsequent GW calculation
         ! 11.2.2004 Arno Schindlmayr
!         IF ((input%gw==1 .or. input%gw==3).AND.(fmpi%irank==0)) THEN
!            OPEN (15,file='ecore',status='unknown', action='write',form='unformatted')
!         END IF

         rh = 0.0
         tec = 0.0
         qint = 0.0
         IF (input%frcor) THEN
            IF (fmpi%irank==0) THEN
               IF(isCoreDensityPresent()) THEN
                  CALL readCoreDensity(input,atoms,rh,tec,qint)
                  l_coreDenPresent = .TRUE.
               END IF
            END IF
 !  #ifdef CPP_MPI
!            CALL MPI_BCAST(l_CoreDenPresent,1,MPI_LOGICAL,0,fmpi%mpi_comm,ierr)
!            CALL mpi_bc_coreDen(fmpi,atoms,input,rh,tec,qint)
 !  #endif
         END IF
      END IF

      ! (1) set up radial mesh beyond muffin-tin radius
      ! (2) cut_off core tails from noise

      nloop : do itype = 1, atoms%ntype
        if ( atoms%econf(itype)%num_core_states > 0 ) THEN
          do j = 1, atoms%jri(itype)
            rat(j, itype) = atoms%rmsh(j, itype)
          end do
          dxx = exp(atoms%dx(itype))
          do j = atoms%jri(itype) + 1, atoms%msh
            rat(j, itype) = rat(j - 1, itype) * dxx
          end do
          do j = atoms%jri(itype) - 1, atoms%msh
            rh(j, itype, 1) = rh(j, itype, 1) / (fpi_const * rat(j, itype) * rat(j, itype))
          end do
          do j = atoms%msh, atoms%jri(itype), -1
            if ( rh(j,itype, 1) > tl_14 ) THEN
              mshc(itype) = j
              cycle nloop
            end if
          end do
          mshc(itype) = atoms%jri(itype)
        end if
      end do nloop

      ! the core density inside the spheres is replaced by a
      ! gaussian-like pseudo density : n(r) = acoff*exp(-alpha*r*r)
      ! acoff and alpha determined to obtain a continous and
      ! differentiable density at the sphere boundary.
      ! IF mshc = jri  either core tail too small or no core (i.e. H)
      !
      do itype = 1,atoms%ntype
        if ( ( mshc(itype) > atoms%jri(itype) ) .and. ( ( atoms%econf(itype)%num_core_states > 0 ) ) ) then
          j1 = atoms%jri(itype) - 1
          if ( method1 == 1) then
             dif = diflgr( rat(j1, itype), rh(j1, itype, 1) )
             !write (logUnit, fmt=8000) itype, rh(atoms%jri(itype),itype, 1), dif
             alpha(itype) = -0.5 * dif / ( rh(atoms%jri(itype), itype, 1) * atoms%rmt(itype) )
          else if ( method1 == 2) then
             alpha(itype) = log( rh(j1,itype, 1) / rh(atoms%jri(itype),itype, 1) )
             alpha(itype) = alpha(itype) / ( atoms%rmt(itype) * atoms%rmt(itype) * ( 1.0 - EXP( -2.0 * atoms%dx(itype) ) ) )
          else
             !write ( logUnit, '('' error in choice of method1 in cdnovlp '')' )
             call juDFT_error("error in choice of method1 in cdnovlp", calledby ="cdnovlp")
          end if
          acoff(itype) = rh(atoms%jri(itype), itype, 1) * exp( alpha(itype) * atoms%rmt(itype) * atoms%rmt(itype) )
          !write (logUnit, FMT=8010) alpha(itype), acoff(itype)
          do j = 1, atoms%jri(itype) - 1
             rh(j, itype, 1) = acoff(itype) * exp( -alpha(itype) * rat(j, itype)**2 )
          end do
        else
          alpha(itype) = 0.0
        end if
      end do

  !8000 FORMAT (/, 10x, 'core density and its first derivative at sph. bound. for atom type', i2, ' is', 3x, 2e15.7 )
  !8010 FORMAT (/, 10x, 'alpha=', f10.5,5x, 'acoff=', f10.5)

      ! Tolerance when alpha is assumed to be zero.
      tl_14 = 1.0e-10!-14
      ! calculate the fourier transform of the core-pseudocharge
      call ft_of_CorePseudocharge( atoms, mshc, alpha, tl_14, rh, acoff, stars, method2, rat, cell, sym, ngdp2km, gdp, qpwcG )

    end subroutine calcPsDensMT

    !---------------------------------------------------------------------------------------------------------------------------------
    !> @author
    !> Christian-Roman Gerhorst, Forschungszentrum Jülich: IAS1 / PGI1
    !>
    !> @brief
    !> Calculates the gradient of the pseudo-core Gauss function density in the displaced muffin-tin
    !>
    !> @details
    !> See also 7.26b (dissertation CRG)
    !>
    !> @param[in]   atoms             : Atoms type, see types.f90
    !> @param[in]   lcApsDens         : Lowercase A of Gauss curve pseudo core-density, see 7.26b (dissertation CRG)
    !> @param[in]   ucApsDens         : Uppercase A of Gauss curve pseudo core-density, see 7.26b (dissertation CRG)
    !> @param[out]  rho1MTAlphPSCcorr : Spherical-harmonic expansion coefficients of the gradient of the pseudo core-density ( Gauss
    !>                                  curve ) at the displaced atom.
    !---------------------------------------------------------------------------------------------------------------------------------
    subroutine calcPDCinAlph( atoms, lcApsDens, ucApsDens, rho1MTAlphPSCcorr )

      use m_types, only : t_atoms

      implicit none

      ! Type parameter
      type(t_atoms),              intent(in)  :: atoms

      ! Array parameter
      real,                       intent(in)  :: lcApsDens(:) ! lowercase A of pseudoDensity
      real,                       intent(in)  :: ucApsDens(:) ! uppercase A of pseudoDensity
      complex,                    intent(out) :: rho1MTAlphPSCcorr(:, : , :, :)

      ! Scalar variables
      integer                                 :: itype
      integer                                 :: imesh
      integer                                 :: lm
      integer                                 :: iatom
      integer                                 :: ieqat
      integer                                 :: idir
      integer                                 :: ikpt

      ! Array variables
      real,          allocatable              :: prFac(:)
      real,          allocatable              :: beforeSum(:, :)

      allocate( prFac(atoms%ntype), beforeSum(atoms%jmtd, atoms%ntype) )
      prFac(:) = 0.
      beforeSum(:, :) = 0.

      rho1MTAlphPSCcorr(:, :, :, :) = cmplx(0., 0.)

      ! todo is this really a itype story or iatom for the lc and uc Aps
      prFac(:) = 0.
      do itype = 1, atoms%ntype
        prFac(itype) = 2. * lcApsDens(itype) * ucApsDens(itype)
        beforeSum(:, :) = 0.
        do imesh = 1, atoms%jri(itype)
          beforeSum(imesh, itype) = prFac(itype) * atoms%rmsh(imesh, itype) * exp( -lcApsDens(itype) * atoms%rmsh(imesh, itype)**2 )
        end do
      end do

      ! todo with testing!!!! transpose c_im
      do idir = 1, 3
        iatom = 0
        do itype = 1, atoms%ntype
          do ieqat = 1, atoms%neq(itype)
            iatom = iatom + 1
            do lm = 2, 4
              do imesh = 1, atoms%jri(itype)
                rho1MTAlphPSCcorr(imesh, lm, iatom, idir) = rho1MTAlphPSCcorr(imesh, lm, iatom, idir) + beforeSum(imesh, itype) &
                                                                                                                & * c_im(idir, lm - 1)
              end do
            end do
          end do
        end do
      end do

    end subroutine calcPDCinAlph

    subroutine ft_of_CorePseudocharge(atoms,mshc,alpha,&
          tol_14,rh,acoff,stars,method2,rat,cell,sym, ngpqdp, gpqdp, qpwcG)

    !=====> calculate the fourier transform of the core-pseudocharge
    !
    !     qpw(\vec{G}) = Sum_{n} [ F(G,n) * Sum_{atm{n}} S(\vec{G},atm) ]
    !                  n = atom_type
    !                  F = Formfactor = F_in_sphere + F_outsphere
    !                  S = Structure factor

    USE m_types

!      type(t_mpi)      ,intent(in) :: mpi
    type(t_atoms)    ,intent(in) :: atoms
    integer          ,intent(in) :: mshc(atoms%ntype)
    real             ,intent(in) :: alpha(atoms%ntype), tol_14
    real             ,intent(in) :: rh(atoms%msh,atoms%ntype)
    real             ,intent(in) :: acoff(atoms%ntype)
    type(t_stars)    ,intent(in) :: stars
    integer          ,intent(in) :: method2
    real             ,intent(in) :: rat(atoms%msh,atoms%ntype)
    type(t_cell)     ,intent(in) :: cell
    ! 
     
    type(t_sym)      ,intent(in) :: sym
    integer, intent(in) :: ngpqdp
    integer, intent(in) :: gpqdp(:, :)
    !complex         ,intent(out) :: qpwc(stars%ng3, atoms%ntype)
    !complex         ,intent(out) :: qpwc(stars%ng3, atoms%nat)
    complex         ,intent(out) :: qpwcG(:, :)

    integer :: ieqat, iatom

!     ..Local variables
    integer nat1, n, n_out_p, k
    complex czero

!     ..Local arrays
!      real :: qf(stars%ng3)
    real, allocatable :: qfG(:)
    !complex qpwc_at(stars%ng3)
    complex, allocatable :: qpwc_atG(:)
!#ifdef CPP_MPI
!      external mpi_bcast
!      complex :: qpwc_loc(stars%ng3)
!      integer :: ierr
!      include "mpif.h"
!#endif

    allocate(qfG(ngpqdp), qpwc_atG(ngpqdp))
    qfG(:) = 0.
    qpwc_atG(:) = cmplx(0., 0.)
    czero = (0.0,0.0)
!#ifdef CPP_MPI
!      DO k = 1 , stars%ng3
!         qpwc_loc(k) = czero
!      ENDDO
!#endif
    !DO  n = 1 + mpi%irank, atoms%ntype, mpi%isize
    !  DO k = 1 , stars%ng3
    !      qpwc(k, n) = czero
    !  ENDDO
    !end do
    qpwcG(:, :) = cmplx(0., 0.)

    !
    !*****> start loop over the atom type
    !
    iatom = 0
!      DO  n = 1 + mpi%irank, atoms%ntype, mpi%isize
    DO  n = 1 , atoms%ntype
        IF ( ( mshc(n) .GT. atoms%jri(n) ).AND.&
            &        ( alpha(n) .GT. tol_14 ) )    THEN

            n_out_p = mshc(n)-atoms%jri(n)+1

            ! (1) Form factor for each atom type

            CALL FormFactor_forAtomType(atoms,method2,n_out_p,&
                               atoms%rmt(n),atoms%jri(n),atoms%dx(n),mshc(n),rat(:,n), &
                               rh(:,n),alpha(n),stars,cell,acoff(n), ngpqdp, gpqdp, qfG)

            ! (2) structure constant for each atom of atom type

!              nat1 = 1
!              IF (n>1) THEN
!                 DO k = 1, n-1
!                    nat1 = nat1 + atoms%neq(k)
!                 END DO
!              END IF
            do ieqat = 1, atoms%neq(n)
              iatom = iatom + 1
       !     CALL StructureConst_forAtom(nat1,stars ,sym,&
            CALL StructureConst_forAtom(iatom,stars ,sym,&
                               atoms%neq(n),atoms%nat,atoms%taual, ngpqdp, gpqdp, &
!                                 cell,qf,qpwc_at)
                               cell,qfG,qpwc_atG)
!#ifdef CPP_MPI
!              DO k = 1, stars%ng3
!                 qpwc_loc(k) = qpwc_loc(k)  + qpwc_at(k)
!              END DO
!#else
!              DO k = 1 , stars%ng3
            DO k = 1 , ngpqdp
            !qpwc(k, n) = qpwc(k, n) + qpwc_at(k) !this is changed because we don't need the sum c.g.
            qpwcG(k, iatom) = qpwcG(k, iatom) + qpwc_atG(k) !this is changed because we don't need the sum c.g.
            END DO
          end do
!#endif

        END IF
     ENDDO
!#ifdef CPP_MPI
!       CALL mpi_allreduce(qpwc_loc,qpwc,stars%ng3,CPP_MPI_COMPLEX,mpi_sum, &
!               mpi%mpi_comm,ierr)
!#endif

    end subroutine ft_of_CorePseudocharge

    subroutine StructureConst_forAtom(nat1,stars , sym,&
                        neq,natd,taual,ngpqdp, gpqdp, cell,qfG,qpwc_atG)
!       calculates structure constant for each atom of atom type

    USE m_types
    USE m_spgrot
    USE m_constants
     

     integer          ,intent(in) :: nat1
     type(t_stars)    ,intent(in) :: stars
      
     type(t_sym)      ,intent(in) :: sym
     integer          ,intent(in) :: neq,natd
     real             ,intent(in) :: taual(3,natd)
     type(t_cell)     ,intent(in) :: cell
     !real             ,intent(in) :: qf(stars%ng3)
     real             ,intent(in) :: qfG(:)
     !complex         ,intent(out) :: qpwc_at(:)
     integer, intent(in) :: ngpqdp
     integer, intent(in) :: gpqdp(:, :)
     complex         ,intent(out) :: qpwc_atG(:)
     complex :: iu

!      ..Local variables
    integer k, nat2, nat, j
    real x
    complex sf, czero

!      ..Local arrays
     integer kr(3,sym%nop)
!       real    kro(3 %ods%nop)
     complex phas(sym%nop)
    
    czero = (0.0,0.0)
    !DO k = 1 , stars%ng3
    !    qpwc_at(k) = czero
    !ENDDO
    qpwc_atG(:) = cmplx(0., 0.)
    iu = cmplx(0., 1.)

    !
    !    first G=0
    !
    !k=1
    !qpwc_at(k)      = qpwc_at(k)      + neq * qf(k)
    !qpwc_at(k)      = qpwc_at(k)      + qf(k)

    !
    !    then  G>0
    !
!!$OMP PARALLEL DO DEFAULT(none) &
!!$OMP SHARED(stars ,sym,neq,natd,nat1,taual,cell,qf,qpwc_at) &
!!$OMP FIRSTPRIVATE(czero) &
! !$OMP PRIVATE(k,kr,phas,nat2,nat,sf,j,x,kro)
    DO  k = 1,ngpqdp
      ! G=0
      !if (all(gpqdp(:, k) == 0)) then
      !  qpwc_atG(k)      = qpwc_atG(k)      + qf(k)
      !  cycle
      !end if
    !DO  k = 2,stars%ng3
!          
        !    CALL spgrot(&
        !         &                       sym%nop,sym%symor,sym%mrot,sym%tau,sym%invtab,&
        !         &                       stars%kv3(:,k),&
        !         &                       kr,phas)
            !
            ! ----> start loop over equivalent atoms
            !
            ! nat2 = nat1 + neq - 1
            ! DO  nat = nat1,nat2
               !  sf = czero
        !         DO  j = 1,sym%nop
                 !    x = -tpi_const* ( kr(1,j) * taual(1,nat)&
                 !        &          + kr(2,j) * taual(2,nat)&
                 !        &          + kr(3,j) * taual(3,nat))
                     !x = -tpi_const* ( kr(1,j) * taual(1,nat1)&
                     !    &          + kr(2,j) * taual(2,nat1)&
                     !    &          + kr(3,j) * taual(3,nat1))
                     !gb      sf = sf + CMPLX(COS(x),SIN(x))*phas(j)
                     !sf = sf + CMPLX(COS(x),SIN(x))*conjg(phas(j))
        !         ENDDO
                 !sf = sf / REAL( sym%nop )
                 qpwc_atG(k)      = qpwc_atG(k)      + exp(-tpi_const * iu * dot_product(gpqdp(1:3, k), taual(1:3, nat1))) * qfG(k)
             !ENDDO
!          
    ENDDO
!!$OMP END PARALLEL DO
    end subroutine StructureConst_forAtom

!----------------------------------------------------------------------
    subroutine FormFactor_forAtomType(atoms,method2,n_out_p,&
                     rmt,jri,dx,mshc,rat,&
                     rh,alpha,stars,cell,acoff, ngpqdp, gpqdp, qfG)

    USE m_types
    USE m_constants
    USE m_rcerf
    USE m_intgr, ONLY : intgr3,intgz0

    type(t_atoms)    ,intent(in) :: atoms
    integer          ,intent(in) :: method2, n_out_p
    real             ,intent(in) :: rmt
    integer          ,intent(in) :: jri
    real             ,intent(in) :: dx
    integer          ,intent(in) :: mshc
    real             ,intent(in) :: rat(atoms%msh)
    real             ,intent(in) :: rh(atoms%msh)
    real             ,intent(in) :: alpha
    type(t_stars)    ,intent(in) :: stars
    type(t_cell)     ,intent(in) :: cell
    real             ,intent(in) :: acoff
    integer,         intent(in) :: ngpqdp
    integer,         intent(in) :: gpqdp(:, :)
!      real            ,intent(out) :: qf(stars%ng3)
    real            ,intent(out) :: qfG(:)


!     ..Local variables
    real f11, f12, ar, g, ai, qfin, qfout, gr, a4, alpha3, zero
    integer k, ir, j
    logical tail
    real :: gExt(3)

!     ..Local arrays
    real rhohelp(atoms%msh)

    zero = 0.0
!    DO k = 1,stars%ng3
!      qf(k) = 0.0
!    END DO
     qfG(:) = 0.


    tail = .FALSE.
    f11 = tpi_const * rmt * rh(jri) / alpha
    f12 = acoff * ( pi_const/alpha ) *SQRT(pi_const/alpha)
    ar  = SQRT( alpha ) * rmt

!!$OMP PARALLEL DO DEFAULT(none) &
!!$OMP SHARED(stars,f11,f12,ar,method2,n_out_p,jri,rat,rh,dx,tail) &
!!$OMP SHARED(alpha,cell,mshc,rmt,qfG) &
!!$OMP FIRSTPRIVATE(zero) &
!!$OMP PRIVATE(k,g,ai,qfin,ir,j,rhohelp,qfout,gr,a4,alpha3)
    !DO  k = 1,stars%ng3
    DO  k = 1,ngpqdp
        !g = stars%sk3(k)
        gExt(1:3) = matmul( cell%bmat(1:3, 1:3), gpqdp(1:3, k) )
        g = norm2(gExt(1:3))
        !    first G=0
        !IF ( k.EQ.1 ) THEN
        IF ( g < 1e-7 ) THEN
            ai = zero
            !
            ! ---->     calculate form factor inside the mt-sphere
            !           (use analytic integration of gaussian)
            !
            qfin = - f11 + f12 * rcerf(ar,ai)
            !
            ! ---->     calculate form factor outside the mt-sphere
            !           (do numerical integration of tails)
            !
            IF ( method2 .EQ. 1) THEN

                DO ir = 1 , n_out_p
                   j  = jri+ir-1
                   rhohelp(mshc+1-j) =  rat(j) * rat(j) &
                         &                       * rat(j) *  rh(j)
                END DO
                CALL intgz0(rhohelp,dx,n_out_p,qfout,tail)

            ELSE

                DO ir = 1 , n_out_p
                   j  = jri+ir-1
                   rhohelp(ir) = rat(j) * rat(j) * rh(j)
                END DO
                CALL intgr3(rhohelp,rat(jri),dx,&
                         &                        n_out_p,qfout)
                !--->             have to remove the small r-correction from intgr3
                qfout=qfout-rmt*rhohelp(1)

            END IF

            qfout = fpi_const * qfout
            !
        ELSE
            !    then  G>0
            ai = 0.5*g/SQRT(alpha)
            gr = g*rmt
            a4 = 0.25/alpha
            !
            ! ---->     calculate form factor inside the mt-sphere
            !           (use analytic integration of gaussian)
            !
            qfin = - f11 * SIN(gr)/gr &
                         &                + f12 * rcerf(ar,ai) * EXP(-a4*g*g)
            !
            ! ---->     calculate form factor outside the mt-sphere
            !           (do numerical integration of tails)

            IF ( method2 .EQ. 1) THEN

                DO ir = 1 , n_out_p
                    j  = jri+ir-1
                    rhohelp(mshc-jri+2-ir) =  rat(j)*rat(j) &
                                  &                                     * rh(j) * SIN( g*rat(j) )
                END DO
                !
                !--->       note we use here the integration routine for vacuum. Because
                !           the vacuum integration is made for an inwards integration
                !           from outside to inside. Outside the starting value will be
                !           nearly zero since the core density is small. if the tail
                !           correction (tail=.true.) is included for the integrals, the
                !           integrand is from infinity inward. This might not be
                !           necessary. Further the integration routine is made for
                !           equidistant meshpoints, therefore the term r(i) of
                !           dr/di = dx*r(i) is included in rhohelp


                CALL intgz0(rhohelp,dx,n_out_p,qfout,tail)

            ELSE

                DO ir = 1 , n_out_p
                    j  = jri+ir-1
                    rhohelp(ir) = rat(j) * rh(j) * SIN(g*rat(j))
                END DO
                CALL intgr3(rhohelp,rat(jri),dx,&
                            &                        n_out_p,qfout)
                !--->             have to remove the small r-correction from intgr3
                !roa...correction.from.intgr3.......................
                IF (rhohelp(1)*rhohelp(2).GT.zero) THEN
                    alpha3 = 1.0 + LOG(rhohelp(2)/rhohelp(1))/dx
                    IF (alpha3.GT.zero)&
                            &                 qfout = qfout - rat(jri)*rhohelp(1)/alpha3
                    ENDIF
                !roa...end.correction...............................


                END IF

                qfout = fpi_const * qfout / g
                !
        END IF
        !
        qfG(k)    = (qfin + qfout)/cell%omtil
    ENDDO
!!$OMP END PARALLEL DO
    end subroutine FormFactor_forAtomType

    subroutine tlmplm4H0( atoms, enpara, usdus, input, td, jsp, rbas1, rbas2, uuilon, duilon, ulouilopn, ilo2p, vr0Sph, loosetdout )

      use m_types, only : t_atoms, t_enpara, t_usdus, t_input, t_tlmplm
      !use m_intgr, only : intgr3LinIntp
      use m_intgr, only : intgr3LinIntp
      use m_gaunt, only : gaunt1
      use m_juDFT_stop, only : juDFT_error
      use m_jpSternhPulaySurface, only : tlo4HS0

      implicit none

      ! Type Parameters
      type(t_atoms),               intent(in)  :: atoms
      type(t_enpara),              intent(in)  :: enpara
      type(t_usdus),               intent(in)  :: usdus
      type(t_input),               intent(in)  :: input
      type(t_tlmplm),              intent(out) :: td

      ! Scalar Parameters
      integer,                     intent(in)  :: jsp

      ! Array Parameters
      real,                        intent(in)  :: rbas1(:,:,0:,:,:)
      real,                        intent(in)  :: rbas2(:,:,0:,:,:)
      real,                        intent(in)  :: uuilon(:,:)
      real,                        intent(in)  :: duilon(:,:)
      real,                        intent(in)  :: ulouilopn(:,:,:)
      integer,                     intent(in)  :: ilo2p(:, :)
      complex,                     intent(in)  :: vr0Sph(:, :, :)
      COMPLEX, ALLOCATABLE,        INTENT(OUT) :: loosetdout(:, :, :, :) ! TODO: This is a hacky solution since the type was modded.

      ! Scalar Variables
      complex                                  :: cil
      real                                     :: tempReal
      real                                     :: tempImag
      integer                                  :: i
      integer                                  :: l
      integer                                  :: l2
      integer                                  :: lamda
      integer                                  :: lm
      integer                                  :: lmin
      integer                                  :: lmin0
      integer                                  :: lmp
      integer                                  :: lmpl
      integer                                  :: lmplm
      integer                                  :: lmx
      integer                                  :: lmxx
      integer                                  :: lp
      integer                                  :: lp1
      integer                                  :: lpl
      integer                                  :: mp
      integer                                  :: mu
      integer                                  :: n
      integer                                  :: na
      integer                                  :: m
      integer                                  :: err
      integer                                  :: mlotot
      integer                                  :: mlolotot
      integer                                  :: mlot_d
      integer                                  :: mlolot_d
      integer                                  :: ieqat
      integer                                  :: lmsph
      real                                     :: tempRbas
      integer :: lmd, lmplmd

      ! Array Variables
      integer,        allocatable              :: indt(:)
      complex,        allocatable              :: dvd(:, :)
      complex,        allocatable              :: dvu(:, :)
      complex,        allocatable              :: uvd(:, :)
      complex,        allocatable              :: uvu(:, : )
      real,           allocatable              :: xReal(:)
      real,           allocatable              :: xImag(:)

      ! Initialization which is taken from eigen.F90 in Fleur
      mlotot = 0
      mlolotot = 0
      DO n = 1, atoms%ntype
        do na = 1, atoms%neq(n)
         mlotot = mlotot + atoms%nlo(n)
         mlolotot = mlolotot + atoms%nlo(n) * (atoms%nlo(n) + 1) / 2
        end do
      ENDDO
      mlot_d = max(mlotot,1)
      mlolot_d = max(mlolotot,1)

      err = 0
      lmd=atoms%lmaxd*(atoms%lmaxd+2)
      lmplmd=(lmd*(lmd+3))/2

      allocate( indt(0:lmplmd) )
      allocate( dvd(0:atoms%lmaxd * (atoms%lmaxd + 3) / 2,(atoms%lmaxd + 1)**2 ) )
      allocate( dvu(0:atoms%lmaxd * (atoms%lmaxd + 3) / 2,(atoms%lmaxd + 1)**2 ) )
      allocate( uvd(0:atoms%lmaxd * (atoms%lmaxd + 3) / 2,(atoms%lmaxd + 1)**2 ) )
      allocate( uvu(0:atoms%lmaxd * (atoms%lmaxd + 3) / 2,(atoms%lmaxd + 1)**2 ) )
      allocate( xReal(atoms%jmtd), xImag(atoms%jmtd) )
      !allocate(td%tuu(0:lmplmd,atoms%nat,1),stat=err)
      !allocate(td%tud(0:lmplmd,atoms%nat,1),stat=err)
      !allocate(td%tdd(0:lmplmd,atoms%nat,1),stat=err)
      !allocate(td%tdu(0:lmplmd,atoms%nat,1),stat=err)
      ALLOCATE(loosetdout(0:lmplmd,atoms%nat,1,4),stat=err)
      allocate(td%tdulo(0:lmd,-atoms%llod:atoms%llod,mlot_d,1,1),stat=err) ! TODO: This needed a scond spin.
      allocate(td%tuulo(0:lmd,-atoms%llod:atoms%llod,mlot_d,1,1),stat=err) ! TODO: Same.
      allocate(td%tuloulo(-atoms%llod:atoms%llod,-atoms%llod:atoms%llod,mlolot_d,1,1), stat=err) ! TODO: Same.
      allocate(td%ind(0:lmd,0:lmd,atoms%nat,1),stat=err )
      if (err.ne.0) then
         !write (logUnit,'(a)') 'eigen: an error occured during allocation of'
         !write (logUnit,'(a)') 'the tlmplm%tuu, tlmplm%tdd etc.: ',err,'  size: ',mlotot
         CALL juDFT_error("eigen: Error during allocation of tlmplm, tdd  etc.",calledby ="eigen")
      end if

      indt(:) = 0
      dvd(:, :) = cmplx(0., 0.)
      dvu(:, :) = cmplx(0., 0.)
      uvd(:, :) = cmplx(0., 0.)
      uvu(:, :) = cmplx(0., 0.)
      xReal(:) = 0.
      xImag(:) = 0.
      !td%tuu(:, :, :) = cmplx(0., 0.)
      !td%tdd(:, :, :) = cmplx(0., 0.)
      !td%tud(:, :, :) = cmplx(0., 0.)
      !td%tdu(:, :, :) = cmplx(0., 0.)
      loosetdout(:, :, :, :) = cmplx(0., 0.)
      td%tdulo(:, :, :, :, :) = cmplx(0., 0.)
      td%tuulo(:, :, :, :, :) = cmplx(0., 0.)
      td%tuloulo(:, :, :, :, :) = cmplx(0., 0.)
      td%ind(:, :, :, :) = -9999

      na = 0
      do n = 1, atoms%ntype
        do ieqat = 1, atoms%neq(n)
        na = na + 1
         !generate the irreducible integrals (u(l'):v(lamda,nu):u(l)) for l' >= l
         do lp = 0, atoms%lmax(n)
            ! NOTE: lp1 is integer
            ! Generation of triangular numbers. NOTE: lp1 is integer
            lp1 = (lp * (lp + 1)) / 2
            do l = 0, lp
               lpl = lp1 + l
               ! loop over non-spherical components of the potential the spherical part (l = 0) will be mixed with the kinetic energy.
              do lamda = 1, atoms%lmax(n)
                lmin = lp - l
                lmx = lp + l
                ! We only have a contribution according to the Gaunt selection rules. The triangular condition can be derived from an
                ! inequalities system so that we can set up the triangular conditions for every index of the Gaunt coefficients.
                ! Furthermore l'+l+lamda should be even
                if ((mod(lamda + lmx, 2) .eq. 1) .or. (lamda.lt.lmin) .or. (lamda.gt.lmx)) then
                  ! No contribution
                  do m = -lamda, lamda ! changed
                    lm = lamda * (lamda + 1) + m + 1 ! changed
                    uvu(lpl, lm) = 0.0
                    dvd(lpl, lm) = 0.0
                    uvd(lpl, lm) = 0.0
                    dvu(lpl, lm) = 0.0
                  end do
                else
                  ! Gaunt coefficient is not zero for this branch. As we have complex expansion coefficients, we have to calculate
                  ! 2 real integrals. This was done for security reasons and to have a similiar method to tlmplm4V which can be
                  ! tested.
                  ! todo In a second step we can use the lattice harmonic coefficients to speed up the routine. In principle we still
                  ! have the Fleur symmetry for the unperturbed quantities but should benchmark the alternative routine against this
                  ! algorithm.
                  xReal = 0.
                  xImag = 0.
                  do m = -lamda, lamda
                    lm = lamda * (lamda + 1) + m + 1
                    ! Calculate the integral <u|V|u>
                    do i = 1, atoms%jri(n)
                      tempRbas = ( rbas1(i, 1, lp, n, 1) * rbas1(i, 1, l, n, 1) + rbas2(i, 1, lp, n, 1) * rbas2(i, 1, l, n, 1) )
                      xReal(i) =  tempRbas * real(vr0SpH(i, lm, na))
                      xImag(i) =  tempRbas * aimag(vr0SpH(i, lm, na))
                    end do
                    call intgr3LinIntp(xReal, atoms%rmsh(1,n), atoms%dx(n), atoms%jri(n), tempReal, 1)
                    call intgr3LinIntp(xImag, atoms%rmsh(1,n), atoms%dx(n), atoms%jri(n), tempImag, 1)
                    uvu(lpl, lm) = cmplx(tempReal, tempImag)
                    ! Calculate the integral <uDot|V|u>
                    do i = 1, atoms%jri(n)
                      tempRbas = ( rbas1(i, 2, lp, n, 1) * rbas1(i, 1, l, n, 1) + rbas2(i, 2, lp, n, 1) * rbas2(i, 1, l, n, 1) )
                      xReal(i) = tempRbas * real(vr0SpH(i, lm, na))
                      xImag(i) = tempRbas * aimag(vr0SpH(i, lm, na))
                    end do
                    call intgr3LinIntp(xReal, atoms%rmsh(1,n), atoms%dx(n), atoms%jri(n), tempReal, 1)
                    call intgr3LinIntp(xImag, atoms%rmsh(1,n), atoms%dx(n), atoms%jri(n), tempImag, 1)
                    dvu(lpl, lm) = cmplx(tempReal, tempImag)
                    ! Calculate the integral <u|V|uDot>
                    do i = 1,atoms%jri(n)
                      tempRbas = (rbas1(i, 1, lp, n, 1) * rbas1(i, 2, l, n, 1) + rbas2(i, 1, lp, n, 1) * rbas2(i, 2, l, n, 1) )
                      xReal(i) = tempRbas * real(vr0SpH(i, lm, na))
                      xImag(i) = tempRbas * aimag(vr0SpH(i, lm, na))
                    end do
                    call intgr3LinIntp(xReal, atoms%rmsh(1,n), atoms%dx(n), atoms%jri(n), tempReal, 1)
                    call intgr3LinIntp(xImag, atoms%rmsh(1,n), atoms%dx(n), atoms%jri(n), tempImag, 1)
                    uvd(lpl, lm) = cmplx(tempReal, tempImag)
                    ! Calculte the integral <uDot|V|uDot>
                    do i = 1,atoms%jri(n)
                      tempRbas = (rbas1(i, 2, lp, n, 1) * rbas1(i, 2, l, n, 1) + rbas2(i, 2, lp, n, 1) * rbas2(i, 2, l, n, 1) )
                      xReal(i) = tempRbas * real(vr0SpH(i, lm, na))
                      xImag(i) = tempRbas * aimag(vr0SpH(i, lm, na))
                    end do
                    call intgr3LinIntp(xReal, atoms%rmsh(1,n), atoms%dx(n), atoms%jri(n), tempReal, 1)
                    call intgr3LinIntp(xImag, atoms%rmsh(1,n), atoms%dx(n), atoms%jri(n), tempImag, 1)
                    dvd(lpl, lm) = cmplx(tempReal, tempImag)
                  end do ! m
                end if ! Gaunt selection rules
              end do ! lh
            end do ! l
          end do ! lp
          ! generate the various t(l'm',lm) matrices for l'm'.ge.lm
          indt(:) = 0
          ! loop over l'm'
          do lp = 0, atoms%lmax(n)
            ! NOTE: lp1 is integer
            lp1 = (lp * (lp + 1)) / 2
            do mp = -lp, lp
              lmp = lp * (lp + 1) + mp
              ! NOTE: lmpl is integer
              lmpl = (lmp * (lmp + 1)) / 2
              ! Loop over the non-spherical (l \= 0) components of the potential
              do lamda = 1, atoms%lmax(n)
                lmin0 = abs(lp - lamda)
                ! l' - lambda > l' and -l' + lambda > l' leads to 0 <= lamda <= 2 l'
                IF (lmin0.GT.lp) cycle
                ! lmxx is in principle lp. To ensure, that the oddness of lambda is always compensated (here strictly speaking only at
                ! lmxx) by l, the modulo operation is subtracted. It is subtracted to not exceed lmxx. The lmxx boundary was given by
                ! Gaunt selection rules, actually reading l' + lambda but as we have l <= l' lmax is only l'. We group lambda and the
                ! modulo operation in lmxx + l' + lambda, so that we only have to ensure l + l' is even.
                lmxx = lp - mod(lamda, 2)
                do mu = -lamda, lamda
                  !collection index of lamda and mu!
                  lmsph = lamda * (lamda + 1) + 1 + mu
                  ! selection rule mp = m + mu!
                  m = mp - mu
                  ! lmin = max(|l' - lamda|, |m' - mu|) to ensure that |m = m' - mu| <= l while fulfilling m = m' - mu.
                  ! Imagine for example l' = 1, lamda = 1, m' = -1,mu = 1. Independent of the selection rule, l >= |m|
                  lmin = max(lmin0, abs(m))
                  ! Serves only for ensuring l + l' + lambda = even, if either lmxx is odd or lmin odd, with l2, lmxx gives its
                  ! oddness or eveness to lmin, so as lmxx was l' and l is interated in steps of size 2 l + l' is always even, because
                  ! either l and l' are even or l and l' are odd, therefore l + l' + lambda is even because the eveness is ensured by
                  ! the modulo operation mod(lambda, 2).
                  l2 = abs(lmxx - lmin)
                  ! Corrects lmin so that mod(lmin, 2) = mod(lmxx, 2), lmin is corrected upwards to not exceed the limits given by
                  ! either the Gaunt selection rules (soft) or l >= |m| (hard) condition.
                  lmin = lmin + mod(l2, 2)
                  do l = lmin, lmxx, 2
                    ! collection index of l and m
                    lm = l * (l+1) + m
                    ! lm is always <= lmp!
                    IF (lm.GT.lmp) CYCLE
                    ! index for uvu, dvd, uvd, dvu to access only those relevant after the Gaunt selection rules
                    lpl = lp1 + l
                    ! index similiar to lpl but with ms in it (after the selection Gaunt selection rules)
                    lmplm = lmpl + lm
                    ! After having found out that the Gaunt coefficient is not zero we multiply it.
                    cil = gaunt1(lp, lamda, l, mp, mu, m, atoms%lmaxd)

                    !td%tuu(lmplm, na, jsp) = td%tuu(lmplm, na, jsp) + cil * uvu(lpl, lmsph)
                    !td%tdd(lmplm, na, jsp) = td%tdd(lmplm, na, jsp) + cil * dvd(lpl, lmsph)
                    !td%tud(lmplm, na, jsp) = td%tud(lmplm, na, jsp) + cil * uvd(lpl, lmsph)
                    !td%tdu(lmplm, na, jsp) = td%tdu(lmplm, na, jsp) + cil * dvu(lpl, lmsph)
                    loosetdout(lmplm, na, jsp, 1) = loosetdout(lmplm, na, jsp, 1) + cil * uvu(lpl, lmsph) ! uu
                    loosetdout(lmplm, na, jsp, 2) = loosetdout(lmplm, na, jsp, 2) + cil * uvd(lpl, lmsph) ! ud
                    loosetdout(lmplm, na, jsp, 3) = loosetdout(lmplm, na, jsp, 3) + cil * dvu(lpl, lmsph) ! du
                    loosetdout(lmplm, na, jsp, 4) = loosetdout(lmplm, na, jsp, 4) + cil * dvd(lpl, lmsph) ! dd
                    ! Logical matrix where there are non-vanishing matrix entries in the tlmplm matrices
                    indt(lmplm) = 1
                  end do ! l
                end do ! mem
              end do ! lh
            end do ! mp
          end do ! lp

          ! set up mapping array
          do lp = 0, atoms%lmax(n)
            do mp = -lp, lp
              lmp = lp * (lp + 1) + mp
              do l = 0, atoms%lmax(n)
                do m = -l, l
                  lm = l * (l + 1) + m
                  if (lmp.ge.lm) then
                    lmplm = (lmp * (lmp + 1)) / 2 + lm
                    if (indt(lmplm).NE.0) THEN
                      td%ind(lmp, lm, na, jsp) = lmplm
                    else
                      td%ind(lmp, lm, na, jsp) = -9999
                    end if
                  else
                    ! As we have lm > lmp here (this was not calculated within the routine), we have to transpose, i.e. interchanging
                    ! lmp and lm for the packed storage index.
                    lmplm = (lm* (lm + 1)) / 2 + lmp
                    if (indt(lmplm).NE.0) THEN
                      td%ind(lmp, lm, na, jsp) = -lmplm
                    else
                      td%ind(lmp, lm, na, jsp) = -9999
                    end if
                  end if
                end do ! m
              end do ! l
            end do ! mp
          end do ! lp

          ! set up the t-matrices for the local orbitals, if there are any
          if (.false.) then
            if (atoms%nlo(n).ge.1) then
              call tlo4HS0(atoms, enpara, usdus, input, td, input%jspins, jsp, n, na, vr0SpH(:, :, na), rbas1, rbas2, uuilon, &
                                                                                                           & duilon, ulouilopn, ilo2p)
            end if
          end if

        end do !ieqat
      end do ! n

    end subroutine tlmplm4H0

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
    subroutine uniteEnergyParameters( atoms, enpara, input, El )

      use m_types

      implicit none

      type(t_atoms),               intent(in)  :: atoms
      type(t_enpara),              intent(in)  :: enpara
      type(t_input),               intent(in)  :: input
      real,           allocatable, intent(out) :: El(:, :, :, :)

      integer                                  :: oqn_l !loop variable
      integer                                  :: itype !loop variable
      integer                                  :: isp   !loop variable
      integer                                  :: ilo   !loop variable

      ! Fill in regular energy parameter to index p = 1 similiar to radial solutions
      allocate(El(1 + atoms%nlod, 0:atoms%lmaxd, atoms%ntype, input%jspins))
      do isp = 1, input%jspins
        do itype = 1, atoms%ntype
          do oqn_l = 0, atoms%lmaxd
            El(1, oqn_l, itype, isp) = enpara%el0(oqn_l, itype, isp)
          end do
        end do
      end do

      ! Fill in LO energy parameter to index p (indices LOs for p > 1)
      do isp = 1, input%jspins
        do itype = 1, atoms%ntype
          do ilo = 1, atoms%nlo(itype)
            El(1 + ilo, atoms%llo(ilo, itype), itype, isp) = enpara%ello0(ilo, itype, isp)
          end do
        end do
      end do

    end subroutine uniteEnergyParameters

END MODULE m_dfpt_init
