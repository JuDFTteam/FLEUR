!--------------------------------------------------------------------------------
! Copyright (c) 2021 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------
MODULE m_dfpt_init

    USE m_types
    USE m_constants
    USE m_radfun
    USE m_radflo
    USE m_eig66_io, ONLY : read_eig
    USE m_intgrf,     ONLY : intgrf_init

    IMPLICIT NONE

CONTAINS
    SUBROUTINE dfpt_init(juPhon, sym, input, atoms, sphhar, stars, cell, noco, nococonv, &
                       & kpts, fmpi, results, enpara, rho, vTot, eig_id, nvfull, GbasVec_eig, usdus, rho0, &
                       & grRho0, ngdp, recG, GbasVec, ilst, nRadFun, iloTable, ilo2p, &
                       & uuilonout, duilonout, ulouilopnout, rbas1, rbas2, gridf, z0)

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
        TYPE(t_potden),   INTENT(IN)  :: rho, vTot
        INTEGER,          INTENT(IN)  :: eig_id
        INTEGER,          INTENT(IN)  :: nvfull(:, :), GbasVec_eig(:, :, :, :)

        TYPE(t_usdus),    INTENT(OUT) :: usdus
        TYPE(t_jpPotden), INTENT(OUT) :: rho0, grRho0

        INTEGER,          INTENT(OUT) :: ngdp

        INTEGER,          ALLOCATABLE, INTENT(OUT) :: recG(:, :)
        INTEGER,          ALLOCATABLE, INTENT(OUT) :: GbasVec(:, :)
        INTEGER,          ALLOCATABLE, INTENT(OUT) :: ilst(:, :, :)
        INTEGER,          ALLOCATABLE, INTENT(OUT) :: nRadFun(:,:)
        INTEGER,          ALLOCATABLE, INTENT(OUT) :: iloTable(:, :, :)
        INTEGER,          ALLOCATABLE, INTENT(OUT) :: ilo2p(:, :)
        REAL,             ALLOCATABLE, INTENT(OUT) :: uuilonout(:, :)
        REAL,             ALLOCATABLE, INTENT(OUT) :: duilonout(:, :)
        REAL,             ALLOCATABLE, INTENT(OUT) :: ulouilopnout(:, :, :)
        REAL,             ALLOCATABLE, INTENT(OUT) :: rbas1(:, :, :, :, :)
        REAL,             ALLOCATABLE, INTENT(OUT) :: rbas2(:, :, :, :, :)
        REAL,             ALLOCATABLE, INTENT(OUT) :: gridf(:, :)
        COMPLEX,          ALLOCATABLE, INTENT(OUT) :: z0(:, :, :, :)

        TYPE(t_lapw) :: lapw
        TYPE(t_mat)  :: zMat

        TYPE(t_cdnvalJob) :: cdnvaljob

        INTEGER      :: nG, nSpins, nbasfcn, iSpin, ik, iG, indexG, ifind, nk, iType, ilo, l_lo, oqn_l, iGrid, iOrd
        INTEGER      :: nodeu, noded
        LOGICAL      :: l_real, addG
        REAL         :: bkpt(3), wronk

        INTEGER, ALLOCATABLE :: GbasVec_temp(:, :), ev_list(:), nocc(:, :)
        REAL,    ALLOCATABLE :: we(:), eig(:)
        REAL,    ALLOCATABLE :: f(:, :, :)
        REAL,    ALLOCATABLE :: g(:, :, :)
        REAL,    ALLOCATABLE :: flo(:, :, :)

!#ifdef CPP_MPI
!        INTEGER ierr
!#endif

        IF (noco%l_noco) THEN
           nSpins = 1
        ELSE
           nSpins = input%jspins
        ENDIF

        ALLOCATE(ilst(MAXVAL(nvfull), kpts%nkpt, input%jspins), GbasVec_temp(3, 0), GbasVec(3, 0) )

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
                            IF (.NOT.ALL(ABS(GbasVec(:, ifind) - GbasVec_eig(:, iG, ik, iSpin)).LE.1E-8)) THEN
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
                        GbasVec(1, indexG) = GbasVec_eig(1, iG, ik, iSpin)
                        GbasVec(2, indexG) = GbasVec_eig(2, iG, ik, iSpin)
                        GbasVec(3, indexG) = GbasVec_eig(3, iG, ik, iSpin)

                        IF (uBound(GbasVec, 2).NE.0) THEN
                            DEALLOCATE(GbasVec_temp)
                            ALLOCATE(GbasVec_temp(3, indexG))
                        END IF
                    END IF
                END DO ! iG
            END DO ! ik
        END DO ! iSpin

        DEALLOCATE(GbasVec_temp)

        !nbasfcn = MERGE(MAXVAL(nvfull(1,:))+MAXVAL(nvfull(2,:))+2*atoms%nlotot,MAXVAL(nvfull(1,:))+atoms%nlotot,noco%l_noco)

        nbasfcn = MAXVAL(nvfull(1,:))+atoms%nlotot

        ALLOCATE(z0(nbasfcn, input%neig, kpts%nkpt, nSpins))
        ALLOCATE(nocc(kpts%nkpt, input%jspins))

        z0 = CMPLX(0.0,0.0)
        nocc = 0

        l_real = sym%invs.AND.(.NOT.noco%l_soc).AND.(.NOT.noco%l_noco).AND.atoms%n_hia==0

        CALL zMat%init(l_real, nbasfcn, input%neig)

        DO iSpin = 1, nspins
            CALL cdnvalJob%init(fmpi, input, kpts, noco, results, iSpin)
            DO ik = 1,kpts%nkpt

                nk=cdnvalJob%k_list(ik)

                bkpt=kpts%bk(:, nk)

                CALL lapw%init(input, noco, nococonv, kpts, atoms, sym, nk, cell, .FALSE., fmpi)
                ev_list = cdnvaljob%compact_ev_list(nk, l_empty = .TRUE.)
                nocc(ik, iSpin) = cdnvaljob%noccbd(nk)
                we  = cdnvalJob%weights(ev_list, nk)
                eig = results%eig(ev_list, nk, iSpin)

                IF (fmpi%irank == 0) THEN
                    CALL read_eig(eig_id, nk, iSpin, list = ev_list,  zmat=zMat)
                    !CALL read_eig(eig_id, ik, iSpin, neig = results%neig(ik, iSpin), &
                    !                                & eig = results%eig(:, ik, iSpin))
                END IF

                IF (l_real) THEN
                    z0(:nvfull(iSpin, ik), :results%neig(ik, iSpin), ik, iSpin) = CMPLX(1.0,0.0) * zMat%data_r(:nvfull(iSpin, ik), :results%neig(ik, iSpin))
                ELSE
                    z0(:nvfull(iSpin, ik), :results%neig(ik, iSpin), ik, iSpin) = zMat%data_c(:nvfull(iSpin, ik), :results%neig(ik, iSpin))
                END IF

!#ifdef CPP_MPI
!                CALL MPI_BARRIER(fmpi%mpi_comm,ierr)
!#endif
            END DO
        END DO

        ! The number of radial functions p for a given atom type and orbital quantum number l is at least 2 ( u and udot ).
        nRadFun = 2

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
              CALL radflo(atoms, iType, iSpin, enpara%ello0(1, 1, iSpin), vTot%mt(:, 0, iType, iSpin), &
                        & f, g, fmpi, usdus, usdus%uuilon(:, :, iSpin), usdus%duilon(:, :, iSpin), usdus%ulouilopn(:, :, :, iSpin), flo)

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

        ! TODO: Meditate about the correct normalizations on all rho/V MT/pw.

        ! Initialize unsymmetrized potdens and G vectors.
        nG=COUNT(stars%ig>0)
        ngdp=nG
        CALL rho0%init_jpPotden(1, 0, nG, atoms%jmtd, atoms%lmaxd, atoms%nat, input%jspins)
        CALL grRho0%init_jpPotden(3, 0, nG, atoms%jmtd, atoms%lmaxd + juPhon%jplPlus, atoms%nat, input%jspins)

        ALLOCATE(recG(nG, 3))
        !v0%init_potden_simple(1, 0, nG, atoms%jmtd, atoms%lmaxd, atoms%nat, input%jspins)

        ! Unpack the star coefficients onto a G-representation.
        CALL stars_to_pw(sym, stars, input%jspins, nG, rho%pw, rho0%pw, recG)

        ! Unpack the lattice harmonics onto spherical harmonics.
        CALL lh_to_sh(sym, atoms, sphhar, input%jspins, rho%mt, rho0%mt)

        ! Construct the interstitial gradients.
        CALL pw_gradient(input%jspins, nG, recG, cell%bmat, rho0%pw, grRho0%pw)

        ! Construct the muffin tin gradients.
        CALL mt_gradient(input%jspins, atoms, juPhon%jplPlus, rho0%mt, grRho0%mt)

    END SUBROUTINE dfpt_init

    SUBROUTINE stars_to_pw(sym, stars, jspins, nG, rhostar, rhopw, recG)

        TYPE(t_sym),   INTENT(IN)  :: sym
        TYPE(t_stars), INTENT(IN)  :: stars
        INTEGER,       INTENT(IN)  :: jspins
        INTEGER,       INTENT(IN)  :: nG

        COMPLEX,       INTENT(IN)  :: rhostar(:, :)
        COMPLEX,       INTENT(OUT) :: rhopw(:, :, :, :)
        INTEGER,       INTENT(OUT) :: recG(:, :)

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

    SUBROUTINE lh_to_sh(sym, atoms, sphhar, jspins, rholh, rhosh )

        TYPE(t_sym),    INTENT(IN)  :: sym
        TYPE(t_atoms),  INTENT(IN)  :: atoms
        TYPE(t_sphhar), INTENT(IN)  :: sphhar
        INTEGER,        INTENT(IN)  :: jspins
        REAL,           INTENT(IN)  :: rholh(:, 0:, :, :)
        COMPLEX,        INTENT(OUT) :: rhosh(:, :, :, :, :, :)

        INTEGER :: iSpin, iType, iEqat, iAtom, ilh, iMem, ilm, iR
        INTEGER :: ptsym, l, m

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
                                rhosh(iR, ilm, iAtom, iSpin, 1, 1) = &
                              & rhosh(iR, ilm, iAtom, iSpin, 1, 1) + &
                              & rholh(iR, ilh, iType, iSpin) * sphhar%clnu(iMem, ilh, ptsym)
                            END DO
                        END DO
                    END DO
                END DO
            END DO
        END DO

    END SUBROUTINE lh_to_sh

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
                gExt = matmul( bmat, recG(iG, :) )
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

END MODULE m_dfpt_init
