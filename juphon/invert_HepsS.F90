!--------------------------------------------------------------------------------
! Copyright (c) 2022 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------

MODULE m_invert_HepsS

    IMPLICIT NONE

CONTAINS
    SUBROUTINE invert_HepsS(atoms, noco, juPhon, lapwkpr, zMatkpr, eignukpr, eignuk, nekpr, nocck, l_real, invHepsS)
        ! Subroutine to calculate (H-\epsilon_{\nu k} S)^(-1) as
        ! z_{k'}(\epsilon_{k'}-\epsilon_{\nu k})^(-1)z_{k'}^H, i.e.
        ! in the spectral representation.

        USE m_types

        TYPE(t_atoms),   INTENT(IN) :: atoms
        TYPE(t_noco),    INTENT(IN) :: noco
        TYPE(t_juPhon),  INTENT(IN) :: juPhon
        TYPE(t_lapw),    INTENT(IN) :: lapwkpr
        CLASS(t_mat),    INTENT(IN) :: zMatkpr

        INTEGER,         INTENT(IN) :: nekpr, nocck
        REAL,            INTENT(IN) :: eignukpr(:), eignuk(:)
        LOGICAL,         INTENT(IN) :: l_real

        CLASS(t_mat), ALLOCATABLE, INTENT(OUT) :: invHepsS(:)

        INTEGER :: nbasfcn, nu, iGpr, iG, nupr
        REAL    :: deps, invdeps

        ALLOCATE(invHepsS(nocck))

        nbasfcn = MERGE(lapwkpr%nv(1)+lapwkpr%nv(2)+2*atoms%nlotot,lapwkpr%nv(1)+atoms%nlotot,noco%l_noco)
        DO nu = 1, nocck
            DO iGpr = 1, nbasfcn
                DO iG = 1, nbasfcn
                    DO nupr = 1, nekpr
                        CALL invHepsS(nu)%init(l_real, nbasfcn, nbasfcn)
                        deps = eignukpr(nupr)-eignuk(nocck)
                        invdeps = 0.0
                        IF (deps.GT.juPhon%eDiffcut) invdeps = 1.0 / deps
                        IF (l_real) THEN
                            invHepsS(nu)%data_r(iGpr,iG) = zMatkpr%data_r(iGpr, nupr) * &
                                                         & zMatkpr%data_r(iG, nupr) * invdeps
                        ELSE
                            invHepsS(nu)%data_c(iGpr,iG) = zMatkpr%data_c(iGpr, nupr) * &
                                                   & CONJG(zMatkpr%data_c(iG, nupr)) * invdeps
                        END IF
                    END DO
                END DO
            END DO
        END DO

    END SUBROUTINE invert_HepsS
 END MODULE m_invert_HepsS
