MODULE m_abcof_soc
CONTAINS
  SUBROUTINE abcof_soc(input,atoms,sym, cell,lapw,ne,usdus,&
                   noco,jspin,oneD,nat_start,nat_stop,nat_l,&
                   acof,bcof,ccof,zMat)
    !     ************************************************************
    !     subroutine constructs the a,b coefficients of the linearized
    !     m.t. wavefunctions for each band and atom.       c.l. fu
    !     ************************************************************
#include "cpp_double.h"

    USE m_constants, ONLY : tpi_const
    USE m_setabc1lo
    USE m_sphbes
    USE m_dsphbs
    USE m_ylm
    USE m_types
    USE m_juDFT
    IMPLICIT NONE
    TYPE(t_input),INTENT(IN)  :: input
    TYPE(t_usdus),INTENT(IN)  :: usdus
    TYPE(t_lapw),INTENT(IN)   :: lapw
    TYPE(t_oneD),INTENT(IN)   :: oneD
    TYPE(t_noco),INTENT(IN)   :: noco
    TYPE(t_sym),INTENT(IN)    :: sym
    TYPE(t_cell),INTENT(IN)   :: cell
    TYPE(t_atoms),INTENT(IN)  :: atoms
    TYPE(t_mat),INTENT(IN)    :: zMat
    !     ..
    !     .. Scalar Arguments ..
    INTEGER, INTENT (IN) :: ne,nat_start,nat_stop,nat_l
    INTEGER, INTENT (IN) :: jspin
    !     ..
    !     .. Array Arguments ..
    COMPLEX, INTENT (OUT) :: acof(:,0:,:)!(nobd,0:atoms%lmaxd*(atoms%lmaxd+2),nat_l)
    COMPLEX, INTENT (OUT) :: bcof(:,0:,:)!(nobd,0:atoms%lmaxd*(atoms%lmaxd+2),nat_l)
    COMPLEX, INTENT (OUT) :: ccof(-atoms%llod:,:,:,:)!(-llod:llod,nobd,atoms%nlod,nat_l)
    !     ..
    !     .. Local Scalars ..
    COMPLEX cexp,phase,c_0,c_1,c_2,ctmp,term1
    REAL const,df,r1,s,tmk,wronk
    INTEGER i,j,k,l,ll1,lm,n,natom,nn,iatom,jatom,lmp,m,nkvec,nbasf
    INTEGER inv_f,ie,ilo,iintsp,nintsp,nvmax,lo,natom_l,na2
    !     ..
    !     .. Local Arrays ..
    INTEGER nbasf0(atoms%nlod,atoms%nat)
    REAL dfj(0:atoms%lmaxd),fj(0:atoms%lmaxd),fg(3),fgp(3),fk(3),fkp(3)
    REAL alo1(atoms%nlod),blo1(atoms%nlod),clo1(atoms%nlod)
    COMPLEX ylm( (atoms%lmaxd+1)**2 )
    LOGICAL apw(0:atoms%lmaxd,atoms%ntype)
    REAL,    ALLOCATABLE :: work_r(:)
    COMPLEX, ALLOCATABLE :: work_c(:)
 !$ COMPLEX, ALLOCATABLE :: acof_l(:,:),bcof_l(:,:),ccof_l(:,:,:)

    CALL timestart("abcof")

    IF (zmat%l_real) THEN
       IF (noco%l_soc.AND.sym%invs) CALL judft_error("BUG in abcof, SOC&INVS but real?")
       IF (noco%l_noco) CALL judft_error("BUG in abcof, l_noco but real?")
    ENDIF

    const = 2 * tpi_const/SQRT(cell%omtil)

    acof(:,:,:)   = CMPLX(0.0,0.0)
    bcof(:,:,:)   = CMPLX(0.0,0.0)
    ccof(:,:,:,:) = CMPLX(0.0,0.0)

    !+APW_LO
    DO n = 1, atoms%ntype
       DO l = 0,atoms%lmax(n)
          apw(l,n) = .FALSE.
          DO lo = 1,atoms%nlo(n)
             IF (atoms%l_dulo(lo,n)) apw(l,n) = .TRUE.
          ENDDO
          IF ((input%l_useapw).AND.(atoms%lapw_l(n).GE.l)) apw(l,n) = .FALSE.
       ENDDO
       DO lo = 1,atoms%nlo(n)
          IF (atoms%l_dulo(lo,n)) apw(atoms%llo(lo,n),n) = .TRUE.
       ENDDO
    ENDDO
    !+APW_LO
    !
    nintsp = 1
    !---> loop over the interstitial spin
    DO iintsp = 1,nintsp
       !
       nvmax=lapw%nv(jspin)
       IF (noco%l_ss) nvmax=lapw%nv(iintsp)
       natom_l = 0

       !---> loop over atom types
       DO n = 1,atoms%ntype
          CALL setabc1lo(atoms,n,usdus,jspin,alo1,blo1,clo1)
          !  ----> loop over equivalent atoms
          DO nn = 1,atoms%neq(n)
             natom = SUM(atoms%neq(:n-1)) + nn
             IF ((natom.GE.nat_start).AND.(natom.LE.nat_stop)) THEN
                natom_l = natom_l + 1
                IF (sym%invsat(natom).EQ.2) THEN
                   jatom = sym%invsatnr(natom)
                ELSE
                   jatom = natom
                ENDIF
                IF (zmat%l_real) THEN
                   ALLOCATE ( work_r(ne) )
                ELSE
                   ALLOCATE ( work_c(ne) )
                ENDIF
                !--->        loop over lapws
#ifndef CPP_OLDINTEL
                !$OMP PARALLEL &
                !$OMP& DEFAULT(none)&
                !$OMP& PRIVATE(k,i,work_r,work_c,fg,fk,s,r1,fj,dfj,l,df,wronk,tmk,phase,lo,nkvec,na2,nbasf,&
                !$OMP& j,fkp,fgp,ylm,ll1,m,c_0,c_1,c_2,lmp,inv_f,lm,term1,ctmp,acof_l,bcof_l,ccof_l)&
                !$OMP& SHARED(n,nn,natom,natom_l,noco,atoms,sym,cell,oneD,lapw,nvmax,ne,zMat,usdus,iintsp,&
                !$OMP& alo1,blo1,clo1,jatom,jspin,apw,const,nbasf0,acof,bcof,ccof,nat_start,nat_stop)
                !$   ALLOCATE(acof_l(size(acof,1),0:size(acof,2)-1),bcof_l(size(bcof,1),0:size(bcof,2)-1))
                !$   ALLOCATE(ccof_l(-atoms%llod:atoms%llod,size(ccof,2),size(ccof,3)))
                !$   acof_l=0 ; bcof_l=0 ; ccof_l=0 
                !$OMP DO
#endif
                DO k = 1,nvmax
                   IF (zmat%l_real) THEN
                      work_r(:ne)=zMat%data_r(k,:ne)
                   ELSE
                      work_c(:ne)=zMat%data_c(k,:ne)
                   END IF

                   fg = lapw%gvec(:,k,jspin) 
                   fk = lapw%bkpt + fg
                   s =  DOT_PRODUCT(fk,MATMUL(cell%bbmat,fk))
                   s = SQRT(s)
                   r1 = atoms%rmt(n)*s
                   CALL sphbes(atoms%lmax(n),r1,fj)
                   CALL dsphbs(atoms%lmax(n),r1,fj,dfj)
                   !  ----> construct a and b coefficients
                   DO l = 0,atoms%lmax(n)
                      df = s*dfj(l)
                      wronk = usdus%uds(l,n,jspin)*usdus%dus(l,n,jspin) - usdus%us(l,n,jspin)*usdus%duds(l,n,jspin)
                      IF (apw(l,n)) THEN
                         fj(l) = 1.0*const * fj(l)/usdus%us(l,n,jspin)
                         dfj(l) = 0.0
                      ELSE
                         dfj(l) = const* (usdus%dus(l,n,jspin)*fj(l)-df*usdus%us(l,n,jspin))/wronk
                         fj(l) = const* (df*usdus%uds(l,n,jspin)-fj(l)*usdus%duds(l,n,jspin))/wronk
                      ENDIF
                   ENDDO ! loop over l
                   tmk = tpi_const* DOT_PRODUCT(fk(:),atoms%taual(:,jatom))
                   phase = CMPLX(COS(tmk),SIN(tmk))
                   fkp = MATMUL(fk(:),cell%bmat) ! fkr
                   fgp = MATMUL(fg(:),cell%bmat) ! fgr
                   !     ----> generate spherical harmonics
                   CALL ylm4(atoms%lmax(n),fkp,ylm)
                   !     ----> loop over l
                   DO l = 0,atoms%lmax(n)
                      ll1 = l* (l+1)
                      !     ----> loop over m
                      DO m = -l,l
                         lm = ll1 + m
                         c_0 = CONJG(ylm(lm+1))*phase
                         c_1 = c_0 *  fj(l)
                         c_2 = c_0 * dfj(l)
                         IF (sym%invsat(natom).EQ.2) THEN
                           lmp = ll1 - m
                           inv_f = (-1)**(l-m)
                           c_1 =  conjg(c_1) * inv_f
                           c_2 =  conjg(c_2) * inv_f
                           IF (zmat%l_real) THEN
                              !$ IF (.false.) THEN
                              acof(:ne,lmp,natom_l) = acof(:ne,lmp,natom_l) +  c_1 * work_r(:ne)
                              bcof(:ne,lmp,natom_l) = bcof(:ne,lmp,natom_l) +  c_2 * work_r(:ne)
                              !$ ENDIF
                              !$ acof_l(:ne,lmp) = acof_l(:ne,lmp) +  c_1 * work_r(:ne)
                              !$ bcof_l(:ne,lmp) = bcof_l(:ne,lmp) +  c_2 * work_r(:ne)
                           ELSE
                              !$ IF (.false.) THEN
                              acof(:ne,lmp,natom_l) = acof(:ne,lmp,natom_l) +  c_1 * work_c(:ne)
                              bcof(:ne,lmp,natom_l) = bcof(:ne,lmp,natom_l) +  c_2 * work_c(:ne)
                              !$ ENDIF
                              !$ acof_l(:ne,lmp) = acof_l(:ne,lmp) +  c_1 * work_c(:ne)
                              !$ bcof_l(:ne,lmp) = bcof_l(:ne,lmp) +  c_2 * work_c(:ne)
                           END IF
                         ELSE
                           !     ----> loop over bands
                           IF (zmat%l_real) THEN
                              !$ IF (.false.) THEN
                              acof(:ne,lm,natom_l) = acof(:ne,lm,natom_l) +  c_1 * work_r(:ne)
                              bcof(:ne,lm,natom_l) = bcof(:ne,lm,natom_l) +  c_2 * work_r(:ne)
                              !$ ENDIF
                              !$ acof_l(:ne,lm) = acof_l(:ne,lm) +  c_1 * work_r(:ne)
                              !$ bcof_l(:ne,lm) = bcof_l(:ne,lm) +  c_2 * work_r(:ne)
                           ELSE
                              !$ IF (.false.) THEN
                              acof(:ne,lm,natom_l) = acof(:ne,lm,natom_l) +  c_1 * work_c(:ne)
                              bcof(:ne,lm,natom_l) = bcof(:ne,lm,natom_l) +  c_2 * work_c(:ne)
                              !$ ENDIF
                              !$ acof_l(:ne,lm) = acof_l(:ne,lm) +  c_1 * work_c(:ne)
                              !$ bcof_l(:ne,lm) = bcof_l(:ne,lm) +  c_2 * work_c(:ne)
                           END IF
                         ENDIF
                      ENDDO ! loop over m
                   ENDDO ! loop over l
                   DO lo=1,atoms%nlo(n)
                      DO nkvec=1,lapw%nkvec(lo,jatom) 
                         IF (k==lapw%kvec(nkvec,lo,jatom)) THEN !check if this k-vector has LO attached
                            term1 = 2 * tpi_const/SQRT(cell%omtil) * ((atoms%rmt(n)**2)/2) * phase
                            IF ((sym%invsat(natom)==0).OR.(sym%invsat(natom)==1)) THEN
                               na2=natom
                            ELSE
                               na2 = sym%invsatnr(natom)
                            ENDIF
                            nbasf=lapw%nv(iintsp)+lapw%index_lo(lo,na2)+nkvec
                            l = atoms%llo(lo,n)
                            ll1 = l* (l+1)
                            DO i = 1,ne
                               DO m = -l,l
                                  lm = ll1 + m
                                  !+gu_con
                                  IF ((sym%invsat(natom)==0).OR.(sym%invsat(natom)==1)) THEN
                                     IF (zMat%l_real) THEN
                                        ctmp = zMat%data_r(nbasf,i)*term1*CONJG(ylm(ll1+m+1))
                                     ELSE
                                        ctmp = zMat%data_c(nbasf,i)*term1*CONJG(ylm(ll1+m+1))
                                     ENDIF
                                     !$ IF (.false.) THEN
                                     acof(i,lm,natom_l) = acof(i,lm,natom_l) + ctmp*alo1(lo)
                                     bcof(i,lm,natom_l) = bcof(i,lm,natom_l) + ctmp*blo1(lo)
                                     ccof(m,i,lo,natom_l) = ccof(m,i,lo,natom_l) + ctmp*clo1(lo)
                                     !$ ENDIF
                                     !$ acof_l(i,lm)   = acof_l(i,lm)   +  ctmp*alo1(lo)
                                     !$ bcof_l(i,lm)   = bcof_l(i,lm)   +  ctmp*blo1(lo)
                                     !$ ccof_l(m,i,lo) = ccof_l(m,i,lo) +  ctmp*clo1(lo)
                                  ELSE
                                     ctmp = zMat%data_c(nbasf,i)*CONJG(term1)*ylm(ll1+m+1)*(-1)**(l-m)
                                     lmp = ll1 - m
                                     !$ IF (.false.) THEN
                                     acof(i,lmp,natom_l) = acof(i,lmp,natom_l) +ctmp*alo1(lo)
                                     bcof(i,lmp,natom_l) = bcof(i,lmp,natom_l) +ctmp*blo1(lo)
                                     ccof(-m,i,lo,natom_l) = ccof(-m,i,lo,natom_l) +ctmp*clo1(lo)
                                     !$ ENDIF
                                     !$ acof_l(i,lmp)   = acof_l(i,lmp)   +  ctmp*alo1(lo)
                                     !$ bcof_l(i,lmp)   = bcof_l(i,lmp)   +  ctmp*blo1(lo)
                                     !$ ccof_l(-m,i,lo) = ccof_l(-m,i,lo) +  ctmp*clo1(lo)
                                  ENDIF
                               END DO
                            END DO
                         ENDIF
                      ENDDO
                   END DO
                ENDDO ! loop over LAPWs (k)
#ifndef CPP_OLDINTEL
                !$OMP END DO
                !$OMP CRITICAL
                !$ acof(:,:,natom_l)   = acof(:,:,natom_l)   + acof_l(:,:)
                !$ bcof(:,:,natom_l)   = bcof(:,:,natom_l)   + bcof_l(:,:)
                !$ ccof(:,:,:,natom_l) = ccof(:,:,:,natom_l) + ccof_l(:,:,:)
                !$OMP END CRITICAL
                !$ DEALLOCATE (acof_l,bcof_l,ccof_l)
                !$OMP END PARALLEL
#endif
                IF (zmat%l_real) THEN
                   DEALLOCATE(work_r)
                ELSE               
                   DEALLOCATE(work_c)
                ENDIF
             ENDIF ! nat_start <= natom <= nat_stop
          ENDDO    ! loop over equivalent atoms
       ENDDO       ! loop over atom types
    ENDDO       ! loop over interstitial spin

    CALL timestop("abcof")

  END SUBROUTINE abcof_soc
END MODULE m_abcof_soc
