!--------------------------------------------------------------------------------
! Copyright (c) 2016 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------

MODULE m_hs_int
CONTAINS
  !Subroutine to construct the interstitial Hamiltonian and overlap matrix
  SUBROUTINE hs_int(input,noco,stars,lapw,mpi,cell,isp,bkpt,vpw,&
       smat,hmat)
    USE m_types
    IMPLICIT NONE
    TYPE(t_input),INTENT(IN)      :: input
    TYPE(t_noco),INTENT(IN)       :: noco
    TYPE(t_stars),INTENT(IN)      :: stars
    TYPE(t_cell),INTENT(IN)       :: cell
    TYPE(t_lapw),INTENT(IN)       :: lapw
    TYPE(t_mpi),INTENT(IN)        :: mpi
    INTEGER,INTENT(IN)            :: isp
    REAL,INTENT(IN)               :: bkpt(3)
    COMPLEX,INTENT(IN)            :: vpw(:,:)
    TYPE(t_lapwmat),INTENT(INOUT) :: smat(:,:),hmat(:,:)


    INTEGER :: ispin,jspin,vpw_spin !spin indices
    INTEGER :: i,j,ii(3),iispin,jjspin
    INTEGER :: in
    COMPLEX :: th,ts,phase
    REAL    :: b1(3),b2(3),r2

    IF (noco%l_noco.AND.isp==2) RETURN !was done already
    DO ispin=MERGE(1,isp,noco%l_noco),MERGE(2,isp,noco%l_noco)
       iispin=MAX(ispin,SIZE(smat,1))
       DO jspin=MERGE(1,isp,noco%l_noco),MERGE(2,isp,noco%l_noco)
          jjspin=MAX(jspin,SIZE(smat,1))
          IF (jspin==ispin) THEN
             vpw_spin=ispin
          ELSE
             vpw_spin=3
          ENDIF
          !$OMP PARALLEL DO SCHEDULE(dynamic) DEFAULT(none) &
          !$OMP SHARED(mpi,lapw,stars,input,bkpt,cell,vpw) &
          !$OMP SHARED(jjspin,iispin,ispin,jspin,vpw_spin)&
          !$OMP SHARED(hmat,smat)&
          !$OMP PRIVATE(ii,i,j,in,phase,b1,b2,r2,th,ts)
          DO  i = mpi%n_rank+1,lapw%nv(ispin),mpi%n_size
             !--->    loop over (k+g)
             DO  j = 1,i  
                !-->     determine index and phase factor
                ii = lapw%gvec(:,i,ispin) - lapw%gvec(:,j,jspin)
                in = stars%ig(ii(1),ii(2),ii(3))
                IF (in.EQ.0) CYCLE
                phase = stars%rgphs(ii(1),ii(2),ii(3))
                !+APW_LO
                IF (input%l_useapw) THEN
                   b1=bkpt+lapw%gvec(:,i,ispin)
                   b2=bkpt+lapw%gvec(:,j,jspin)
                   r2 = DOT_PRODUCT(MATMUL(b2,cell%bbmat),b1)   
                   th = phase*(0.5*r2*stars%ustep(in)+vpw(in,vpw_spin))
                ELSE
                   IF (vpw_spin==3.AND.jspin==2) THEN
                      th = phase* (0.25* (lapw%rk(i,ispin)**2+lapw%rk(j,jspin)**2)*stars%ustep(in) + CONJG(vpw(in,vpw_spin)))
                   ELSE
                      th = phase* (0.25* (lapw%rk(i,ispin)**2+lapw%rk(j,jspin)**2)*stars%ustep(in) + vpw(in,vpw_spin))
                   ENDIF
                ENDIF
                !-APW_LO
                !--->    determine matrix element and store
                ts = phase*stars%ustep(in)
                IF (hmat(1,1)%l_real) THEN
                   hmat(jjspin,iispin)%data_r(j,i) = REAL(th)
                   smat(jjspin,iispin)%data_r(j,i) = REAL(ts)
                else
                   hmat(jjspin,iispin)%data_c(j,i) = th
                   smat(jjspin,iispin)%data_c(j,i) = ts
                endif
             ENDDO
          ENDDO
          !$OMP END PARALLEL DO
       ENDDO
    ENDDO
  END SUBROUTINE hs_int
END MODULE m_hs_int
