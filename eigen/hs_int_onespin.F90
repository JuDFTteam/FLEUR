!--------------------------------------------------------------------------------
! Copyright (c) 2022 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------

MODULE m_hs_int_onespin
CONTAINS
    SUBROUTINE hs_int_onespin(input, fmpi, lapw1, lapw2, stars, ispin, jspin, cell, vpw, hmat, smat, l_dfpths)

        USE m_types

        IMPLICIT NONE

        TYPE(t_input),INTENT(IN)      :: input
        TYPE(t_stars),INTENT(IN)      :: stars
        TYPE(t_cell),INTENT(IN)       :: cell
        TYPE(t_lapw),INTENT(IN)       :: lapw1, lapw2
        TYPE(t_mpi),INTENT(IN)        :: fmpi

        INTEGER,INTENT(IN)            :: ispin, jspin
        COMPLEX,INTENT(IN)            :: vpw(:,:)
        CLASS(t_mat),INTENT(INOUT)     :: hmat, smat

        LOGICAL, INTENT(IN) :: l_dfpths

        INTEGER :: i, j, i0, jmax, ii(3)
        INTEGER :: in
        COMPLEX :: th,ts,phase
        REAL    :: b1(3),b2(3),r2

        DO  i = fmpi%n_rank+1,lapw1%nv(ispin),fmpi%n_size
           i0=(i-1)/fmpi%n_size+1
           !--->    loop over (k+g)
           jmax = MERGE(lapw2%nv(jspin),MIN(i,lapw2%nv(jspin)),l_dfpths)
           DO  j = 1, jmax
              ii = lapw1%gvec(:,i,ispin) - lapw2%gvec(:,j,jspin)
              IF (ispin==1.AND.jspin==2.AND..NOT.l_dfpths) THEN
                 ii=-1*ii
                 in = stars%ig(ii(1),ii(2),ii(3))
                 IF (in.EQ.0) CYCLE
                 th = stars%rgphs(ii(1),ii(2),ii(3))*conjg(vpw(in,3))
                 ts=0.0
              ELSEIF(ispin==2.and.jspin==1.AND..NOT.l_dfpths) THEN
              !   ii = -1*ii
                 in = stars%ig(ii(1),ii(2),ii(3))
                 IF (in.EQ.0) CYCLE
                 th = stars%rgphs(ii(1),ii(2),ii(3))*vpw(in,3)
                 ts=0.0
              ELSE
                 !-->     determine index and phase factor
                 in = stars%ig(ii(1),ii(2),ii(3))
                 IF (in.EQ.0) CYCLE
                 phase = stars%rgphs(ii(1),ii(2),ii(3))
                 ts = phase*stars%ustep(in)
                 IF (input%l_useapw.OR.l_dfpths) THEN
                    b1=lapw1%bkpt+lapw1%gvec(:,i,ispin)
                    b2=lapw2%bkpt+lapw2%gvec(:,j,jspin)
                    r2 = DOT_PRODUCT(MATMUL(b2,cell%bbmat),b1)
                    th = phase*(0.5*r2*stars%ustep(in)+vpw(in,ispin))
                 ELSE
                    th = phase* (0.25* (lapw1%rk(i,ispin)**2+lapw2%rk(j,jspin)**2)*stars%ustep(in) + vpw(in,ispin))
                 ENDIF
              ENDIF
              !--->    determine matrix element and store
              IF (hmat%l_real) THEN
                 hmat%data_r(j,i0) = REAL(th)
                 smat%data_r(j,i0) = REAL(ts)
              else
                 hmat%data_c(j,i0) = th
                 smat%data_c(j,i0) = ts
              endif
           ENDDO
        ENDDO
    END SUBROUTINE hs_int_onespin
END MODULE m_hs_int_onespin
