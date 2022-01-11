!--------------------------------------------------------------------------------
! Copyright (c) 2021 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------
MODULE m_dfpt_init

    USE m_types
    USE m_constants
    USE m_eig66_io, ONLY : read_eig

    IMPLICIT NONE

CONTAINS
    SUBROUTINE dfpt_init(juPhon, sym, input, atoms, sphhar, stars, cell, noco, nococonv, &
                       & kpts, fmpi, results, rho, eig_id, nvfull, GbasVec_eig, rho0, &
                       & grRho0, ngdp, recG, GbasVec, ilst, z0)

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
        TYPE(t_potden),   INTENT(IN)  :: rho
        INTEGER,          INTENT(IN)  :: eig_id
        INTEGER,          INTENT(IN)  :: nvfull(:, :), GbasVec_eig(:, :, :, :)

        TYPE(t_jpPotden), INTENT(OUT) :: rho0, grRho0

        INTEGER,          INTENT(OUT) :: ngdp

        INTEGER,          ALLOCATABLE, INTENT(OUT) :: recG(:, :)
        INTEGER,          ALLOCATABLE, INTENT(OUT) :: GbasVec(:, :)
        INTEGER,          ALLOCATABLE, INTENT(OUT) :: ilst(:, :, :)
        COMPLEX,          ALLOCATABLE, INTENT(OUT) :: z0(:, :, :, :)

        TYPE(t_lapw) :: lapw
        TYPE(t_mat)  :: zMat

        TYPE(t_cdnvalJob) :: cdnvaljob

        INTEGER      :: nG, nSpins, nbasfcn, iSpin, ik, iG, indexG, ifind, nk
        LOGICAL      :: l_real, addG
        REAL         :: bkpt(3)

        INTEGER, ALLOCATABLE :: GbasVec_temp(:, :), ev_list(:), nocc(:, :)
        REAL,    ALLOCATABLE :: we(:), eig(:)

!#ifdef CPP_MPI
!        INTEGER ierr
!#endif

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
