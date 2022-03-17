!--------------------------------------------------------------------------------
! Copyright (c) 2022 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------

MODULE m_dfpt_hs_int
CONTAINS
    ! Constructs the interstitial perturbed Hamiltonian and overlap matrix
    SUBROUTINE dfpt_hs_int(noco, starsq, lapwq, lapw, fmpi, bbmat, isp, vpw, smat, hmat)

        USE m_types
        USE m_hs_int_direct

        IMPLICIT NONE

        TYPE(t_noco),INTENT(IN)       :: noco
        TYPE(t_stars),INTENT(IN)      :: starsq
        REAL, INTENT(IN)              :: bbmat(3, 3)
        TYPE(t_lapw),INTENT(IN)       :: lapwq, lapw
        TYPE(t_mpi),INTENT(IN)        :: fmpi
        INTEGER,INTENT(IN)            :: isp
        COMPLEX,INTENT(IN)            :: vpw(:, :)
        CLASS(t_mat),INTENT(INOUT)     :: smat(:,:),hmat(:,:)

        INTEGER :: iSpinPr,iSpin, iMatPr, iMat, iTkin
        LOGICAL :: l_smat
        COMPLEX, ALLOCATABLE :: vpw_temp(:)

        IF (noco%l_noco.AND.isp.EQ.2) RETURN !was done already

        ALLOCATE(vpw_temp(SIZE(vpw, 1)))

        DO iSpinPr = MERGE(1, isp, noco%l_noco), MERGE(2, isp, noco%l_noco)
            ! co:
            ! iSpinPr = isp, SIZE(smat, 1) = 1 (?)
            ! noco:
            ! iSpinPr = 1...2, SIZE(smat, 1) = 2 (?)
            !iispin = MIN(iSpinPr, SIZE(smat, 1))
            ! co:
            ! iispin = 1
            ! noco:
            ! iispin = 1...2
            ! --> alternative: iispin = MERGE(iSpinPr, 1, noco%l_noco) ?
            iMatPr = MERGE(iSpinPr, 1, noco%l_noco)
            DO iSpin=MERGE(1,isp,noco%l_noco),MERGE(2,isp,noco%l_noco)
                iMat = MERGE(iSpin, 1, noco%l_noco)
                iTkin = 0
                ! 1, 2, 3, 4 == 11, 22, 21, 12:
                IF ((iSpinPr.EQ.1).AND.(iSpin.EQ.1)) vpw_temp = vpw(:, 1)
                IF ((iSpinPr.EQ.2).AND.(iSpin.EQ.2)) vpw_temp = vpw(:, 2)
                IF ((iSpinPr.EQ.2).AND.(iSpin.EQ.1)) vpw_temp = vpw(:, 3)
                IF ((iSpinPr.EQ.1).AND.(iSpin.EQ.2)) vpw_temp = vpw(:, 4)

                l_smat = iSpinPr.EQ.iSpin

                IF (iSpinPr.EQ.iSpin) iTkin = 1

                CALL hs_int_direct(fmpi, lapwq%gvec(:, :, iSpinPr), lapw%gvec(:,:,iSpin), lapwq%bkpt, lapw%bkpt, lapw%nv(iSpinPr), lapw%nv(iSpin), &
                                 & starsq, bbmat, vpw_temp, hmat(iMatPr, iMat), smat(iMatPr, iMat), l_smat, .TRUE., iTkin, 1)
            END DO
        END DO
    END SUBROUTINE dfpt_hs_int
END MODULE m_dfpt_hs_int
