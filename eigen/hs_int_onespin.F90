!--------------------------------------------------------------------------------
! Copyright (c) 2022 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------

MODULE m_hs_int_onespin
CONTAINS
    SUBROUTINE hs_int_onespin(input, fmpi, gvec1, gvec2, bkpt1, bkpt2, nv1, nv2, rk1, rk2, stars, ispin, jspin, cell, vpw, hmat, smat, l_dfpths)

        USE m_types

        IMPLICIT NONE

        TYPE(t_input),INTENT(IN)      :: input
        TYPE(t_stars),INTENT(IN)      :: stars
        TYPE(t_cell),INTENT(IN)       :: cell
        INTEGER ,INTENT(IN)           :: gvec1(:, :), gvec2(:, :)
        INTEGER, INTENT(IN)           :: nv1, nv2
        REAL, INTENT(IN)              :: bkpt1(3), bkpt2(3), rk1(:), rk2(:)
        TYPE(t_mpi),INTENT(IN)        :: fmpi

        INTEGER,INTENT(IN)            :: ispin, jspin
        COMPLEX,INTENT(IN)            :: vpw(:,:)
        CLASS(t_mat),INTENT(INOUT)     :: hmat, smat

        LOGICAL, INTENT(IN) :: l_dfpths

        INTEGER :: i, j, i0, jmax, ii(3)
        INTEGER :: in
        COMPLEX :: th,ts,phase
        REAL    :: b1(3),b2(3),r2

        !$OMP PARALLEL DO SCHEDULE(dynamic) DEFAULT(none) &
        !$OMP SHARED(fmpi,stars,input,cell,vpw,gvec1,gvec2,bkpt1,bkpt2,rk1,rk2) &
        !$OMP SHARED(nv1,nv2,ispin,jspin,l_dfpths)&
        !$OMP SHARED(hmat,smat)&
        !$OMP PRIVATE(ii,i0,i,j,jmax,in,phase,b1,b2,r2,th,ts)
        DO  i = fmpi%n_rank+1,nv1,fmpi%n_size
           i0=(i-1)/fmpi%n_size+1
           !--->    loop over (k+g)
           jmax = MERGE(nv2,MIN(i,nv2),l_dfpths)
           DO  j = 1, jmax
              ii = gvec1(:,i) - gvec2(:,j)
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
                    b1=bkpt1+gvec1(:,i)
                    b2=bkpt2+gvec2(:,j)
                    r2 = DOT_PRODUCT(MATMUL(b2,cell%bbmat),b1)
                    th = phase*(0.5*r2*stars%ustep(in)+vpw(in,ispin))
                 ELSE
                    th = phase* (0.25* (rk1(i)**2+rk2(j)**2)*stars%ustep(in) + vpw(in,ispin))
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
        !$OMP END PARALLEL DO
    END SUBROUTINE hs_int_onespin
END MODULE m_hs_int_onespin
