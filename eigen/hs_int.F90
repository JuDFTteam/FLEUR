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
    TYPE(t_lapwmat),INTENT(INOUT) :: smat,hmat


    INTEGER :: sb,ispin,jspin,vpw_spin !spin indices
    INTEGER :: i,j,ii,jj               !loop over k+g, offset in noco case
    INTEGER :: i1,i2,i3,in,i0
    COMPLEX :: th,ts,phase
    REAL    :: b1(3),b2(3),r2

    DO sb=1,MERGE(3,1,noco%l_noco) !loop over three possible spin-blocks (UpUp,DownDown,DownUp)
       CALL lapw%spinblock(noco%l_noco,sb,isp,ispin,jspin,ii,jj,vpw_spin) !determine offsets and spin indices in noco case

       !$OMP PARALLEL DO SCHEDULE(dynamic) DEFAULT(none) &
       !$OMP SHARED(mpi,lapw,stars,input,bkpt,cell,vpw) &
       !$OMP SHARED(ii,jj,ispin,jspin,vpw_spin)&
       !$OMP SHARED(hmat,smat,sb)&
       !$OMP PRIVATE(i,j,i1,i2,i3,in,phase,b1,b2,r2,th,ts)
       DO  i0 = 1,smat%matsize2
          i=smat%local_rk_map(i0,ispin)
          if (i<0) CYCLE !this is an LO
          !--->    loop over (k+g)
          DO  j = 1,MERGE(lapw%nv(1),i,sb==3)  
             !-->     determine index and phase factor
             i1 = lapw%k1(i,ispin) - lapw%k1(j,jspin)
             i2 = lapw%k2(i,ispin) - lapw%k2(j,jspin)
             i3 = lapw%k3(i,ispin) - lapw%k3(j,jspin)
             in = stars%ig(i1,i2,i3)
             IF (in.EQ.0) CYCLE
             phase = stars%rgphs(i1,i2,i3)
             !+APW_LO
             IF (input%l_useapw) THEN
                b1(1) = bkpt(1)+lapw%k1(i,ispin) ; b2(1) = bkpt(1)+lapw%k1(j,jspin)
                b1(2) = bkpt(2)+lapw%k2(i,ispin) ; b2(2) = bkpt(2)+lapw%k2(j,jspin)
                b1(3) = bkpt(3)+lapw%k3(i,ispin) ; b2(3) = bkpt(3)+lapw%k3(j,jspin)
                r2 = DOT_PRODUCT(MATMUL(b2,cell%bbmat),b1)   

                th = phase*(0.5*r2*stars%ustep(in)+vpw(in,vpw_spin))
             ELSE
                th = phase* (0.25* (lapw%rk(i,ispin)**2+lapw%rk(j,jspin)**2)*stars%ustep(in) + vpw(in,vpw_spin))
             ENDIF
             !-APW_LO
             !--->    determine matrix element and store
             ts = phase*stars%ustep(in)
             if (hmat%l_real) THEN
                hmat%data_r(jj+j,ii+i0) = REAL(th)
                smat%data_r(jj+j,ii+i0) = REAL(ts)
             else
                hmat%data_c(jj+j,ii+i0) = th
                smat%data_c(jj+j,ii+i0) = ts
             endif
          ENDDO
       ENDDO
       !$OMP END PARALLEL DO
    ENDDO !Spinblock in noco case
  END SUBROUTINE hs_int
END MODULE m_hs_int
