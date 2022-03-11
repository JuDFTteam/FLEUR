!--------------------------------------------------------------------------------
! Copyright (c) 2022 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------

MODULE m_hs_int_direct
CONTAINS
    SUBROUTINE hs_int_direct(fmpi, gvecPr, gvec, kvecPr, kvec, nvPr, nv, stars, cell, vpw, hmat, smat, l_smat, l_fullj, iTkin)
        ! Calculates matrix elements of the form
        ! <\phi_{k'G'}|M|\phi_{kG}>
        ! for different use cases in the DFT/DFPT scf loop and operators M.
        !
        ! The spin loop and distinction is found on a higher level and translates
        ! into new switches.
        !
        ! DFT:
        ! M = \Theta_{IR} * (T + V) goes into hmat, M = \Theta_{IR} into smat.
        ! [vpw = \Theta_{IR} * V]
        ! [stars%ustep = \Theta_{IR}]
        ! [l_smat = F for offdiags, l_fullj = F]
        ! [iTkin = 0 for offdiags, 1 for l_useapw, 2 else]
        !
        ! DFPT:
        ! M = \Theta_{IR}^{(1)} * (T + V) + \Theta_{IR} * V^{(1)} goes into hmat,
        ! M = \Theta_{IR}^{(1)} into smat.
        ! [vpw = \Theta_{IR}^{(1)} * V + \Theta_{IR} * V^{(1)}]
        ! [stars%ustep = \Theta_{IR}^{(1)}]
        ! [l_smat = ?, l_fullj = T]
        ! [iTkin = 0 for offdiags, 1 else]

        USE m_types

        IMPLICIT NONE

        TYPE(t_stars),INTENT(IN)      :: stars
        TYPE(t_cell),INTENT(IN)       :: cell
        INTEGER ,INTENT(IN)           :: gvecPr(:, :), gvec(:, :)
        INTEGER, INTENT(IN)           :: nvPr, nv, iTkin
        REAL, INTENT(IN)              :: kvecPr(3), kvec(3)
        TYPE(t_mpi),INTENT(IN)        :: fmpi

        COMPLEX,INTENT(IN)            :: vpw(:)
        CLASS(t_mat),INTENT(INOUT)     :: hmat, smat

        LOGICAL, INTENT(IN) :: l_smat, l_fullj

        INTEGER :: i, j, i0, jmax, gPrG(3)
        INTEGER :: gInd
        COMPLEX :: th, ts, phase
        REAL    :: bvecPr(3), bvec(3), r2

        !$OMP PARALLEL DO SCHEDULE(dynamic) DEFAULT(none) &
        !$OMP SHARED(fmpi, stars, cell, vpw, gvecPr, gvec, kvecPr, kvec) &
        !$OMP SHARED(nvPr, nv, l_smat, l_fullj, iTkin)&
        !$OMP SHARED(hmat, smat)&
        !$OMP PRIVATE(gPrG, i0, i, j, jmax, gInd, phase, bvecPr, bvec, r2, th, ts)
        DO  i = fmpi%n_rank + 1, nvPr, fmpi%n_size
            i0 = (i-1) / fmpi%n_size + 1
            jmax = MERGE(nv, MIN(i, nv), l_fullj)
            DO  j = 1, jmax
                gPrG = gvecPr(:,i) - gvec(:,j)

                gInd = stars%ig(gPrG(1), gPrG(2), gPrG(3))
                IF (gInd.EQ.0) CYCLE

                phase = stars%rgphs(gPrG(1), gPrG(2), gPrG(3))

                th = phase * vpw(gInd)

                IF (iTkin.GT.0) THEN
                    bvecPr = kvecPr + gvecPr(:,i)
                    bvec = kvec + gvec(:,j)

                    IF (iTkin.EQ.1) THEN ! Symmetric Dirac form
                        r2 = 0.5 * DOT_PRODUCT(MATMUL(bvecPr,cell%bbmat),bvec)
                    ELSE IF (iTkin.EQ.2) THEN ! Symmetrized Laplace form
                        r2 = 0.25 * DOT_PRODUCT(MATMUL(bvecPr,cell%bbmat),bvecPr)
                        r2 = r2 + 0.25 * DOT_PRODUCT(MATMUL(bvec,cell%bbmat),bvec)
                    ELSE ! Pure Laplace form
                        r2 = 0.5 * DOT_PRODUCT(MATMUL(bvec,cell%bbmat),bvec)
                    END IF
                    th = th + phase * r2 * stars%ustep(gInd)
                END IF

                IF (l_smat) THEN
                    ts = phase*stars%ustep(gInd)
                ELSE
                    ts = CMPLX(0.0, 0.0)
                END IF

                IF (hmat%l_real) THEN
                    hmat%data_r(j,i0) = REAL(th)
                    smat%data_r(j,i0) = REAL(ts)
                ELSE
                    hmat%data_c(j,i0) = th
                    smat%data_c(j,i0) = ts
                END IF
            END DO
        END DO
        !$OMP END PARALLEL DO
    END SUBROUTINE hs_int_direct
END MODULE m_hs_int_direct
