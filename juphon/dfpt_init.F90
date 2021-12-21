!--------------------------------------------------------------------------------
! Copyright (c) 2021 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------
MODULE m_dfpt_init

    USE m_types

    IMPLICIT NONE

CONTAINS
    SUBROUTINE dfpt_init(juPhon, sym, input, atoms, sphhar, stars, rho, rho0, recG)

        TYPE(t_juPhon),   INTENT(IN)  :: juPhon
        TYPE(t_sym),      INTENT(IN)  :: sym
        TYPE(t_input),    INTENT(IN)  :: input
        TYPE(t_atoms),    INTENT(IN)  :: atoms
        TYPE(t_sphhar),   INTENT(IN)  :: sphhar
        TYPE(t_stars),    INTENT(IN)  :: stars
        TYPE(t_potden),   INTENT(IN)  :: rho

        TYPE(t_jpPotden), INTENT(OUT) :: rho0

        INTEGER,          ALLOCATABLE, INTENT(OUT) :: recG(:, :)

        INTEGER :: nG

        ! Initialize unsymmetrized potdens and G vectors.
        nG=COUNT(stars%ig>0)
        CALL rho0%init_jpPotden(1, 0, nG, atoms%jmtd, atoms%lmaxd, atoms%nat, input%jspins)
        ALLOCATE(recG(nG, 3))
        !v0%init_potden_simple(1, 0, nG, atoms%jmtd, atoms%lmaxd, atoms%nat, input%jspins)

        ! Unpack the star coefficients onto a G-representation.
        CALL stars_to_pw(sym, stars, input%jspins, nG, rho%pw, rho0%pw, recG)

        ! Unpack the lattice harmonics onto spherical harmonics.
        CALL lh_to_sh(sym, atoms, sphhar, input%jspins, rho%mt, rho0%mt)

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

        iG    = 0
        rhopw = CMPLX(0.0,0.0)

        DO iSpin = 1, jspins
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

END MODULE m_dfpt_init
