!--------------------------------------------------------------------------------
! Copyright (c) 2016 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------

MODULE m_hs_int
CONTAINS
  !Subroutine to construct the interstitial Hamiltonian and overlap matrix
  SUBROUTINE hs_int(input,noco,stars,lapw,mpi,cell,isp,vpw,&
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
    COMPLEX,INTENT(IN)            :: vpw(:,:)
    CLASS(t_mat),INTENT(INOUT)     :: smat(:,:),hmat(:,:)


    INTEGER :: ispin,jspin !spin indices
    INTEGER :: i,j,ii(3),iispin,jjspin,i0
    INTEGER :: in
    COMPLEX :: th,ts,phase
    REAL    :: b1(3),b2(3),r2

    IF (noco%l_noco.AND.isp==2) RETURN !was done already
    DO ispin=MERGE(1,isp,noco%l_noco),MERGE(2,isp,noco%l_noco)
       iispin=MIN(ispin,SIZE(smat,1))
       DO jspin=MERGE(1,isp,noco%l_noco),MERGE(2,isp,noco%l_noco)
          jjspin=MIN(jspin,SIZE(smat,1))
        
          !$OMP PARALLEL DO SCHEDULE(dynamic) DEFAULT(none) &
          !$OMP SHARED(mpi,lapw,stars,input,cell,vpw) &
          !$OMP SHARED(jjspin,iispin,ispin,jspin)&
          !$OMP SHARED(hmat,smat)&
          !$OMP PRIVATE(ii,i0,i,j,in,phase,b1,b2,r2,th,ts)
          DO  i = mpi%n_rank+1,lapw%nv(ispin),mpi%n_size
             i0=(i-1)/mpi%n_size+1
             !--->    loop over (k+g)
             DO  j = 1,MIN(i,lapw%nv(jspin))  
                ii = lapw%gvec(:,i,ispin) - lapw%gvec(:,j,jspin)
                IF (ispin==1.AND.jspin==2) THEN
                   ii=-1*ii
                   in = stars%ig(ii(1),ii(2),ii(3))
                   IF (in.EQ.0) CYCLE
                   th = stars%rgphs(ii(1),ii(2),ii(3))*conjg(vpw(in,3))           
                   ts=0.0
                ELSEIF(ispin==2.and.jspin==1) THEN
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
                   IF (input%l_useapw) THEN
                      b1=lapw%bkpt+lapw%gvec(:,i,ispin)
                      b2=lapw%bkpt+lapw%gvec(:,j,jspin)
                      r2 = DOT_PRODUCT(MATMUL(b2,cell%bbmat),b1)   
                      th = phase*(0.5*r2*stars%ustep(in)+vpw(in,ispin))
                   ELSE
                      th = phase* (0.25* (lapw%rk(i,ispin)**2+lapw%rk(j,jspin)**2)*stars%ustep(in) + vpw(in,ispin))
                   ENDIF
                ENDIF
                !--->    determine matrix element and store
                IF (hmat(1,1)%l_real) THEN
                   hmat(jjspin,iispin)%data_r(j,i0) = REAL(th)
                   smat(jjspin,iispin)%data_r(j,i0) = REAL(ts)
                else
                   hmat(jjspin,iispin)%data_c(j,i0) = th
                   smat(jjspin,iispin)%data_c(j,i0) = ts
                endif
             ENDDO
          ENDDO
          !$OMP END PARALLEL DO
       ENDDO
    ENDDO
  END SUBROUTINE hs_int
END MODULE m_hs_int
