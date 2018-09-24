MODULE m_abcof_soc
CONTAINS
  SUBROUTINE abcof_soc(input,atoms,sym, cell,lapw,ne,usdus,&
                   noco,jspin,oneD,nat_start,nat_stop,nat_l,&
                   acof,bcof,ccof,zMat,eig,force)
    !     ************************************************************
    !     subroutine constructs the a,b coefficients of the linearized
    !     m.t. wavefunctions for each band and atom.       c.l. fu
    !     ************************************************************
#include "cpp_double.h"

    USE m_constants, ONLY : tpi_const
    USE m_setabc1lo
    USE m_sphbes
    USE m_dsphbs
    USE m_abclocdn_soc
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
    TYPE(t_force),OPTIONAL,INTENT(INOUT) :: force
    !     ..
    !     .. Scalar Arguments ..
    INTEGER, INTENT (IN) :: ne,nat_start,nat_stop,nat_l
    INTEGER, INTENT (IN) :: jspin
    !     ..
    !     .. Array Arguments ..
    COMPLEX, INTENT (OUT) :: acof(:,0:,:)!(nobd,0:dimension%lmd,nat_l)
    COMPLEX, INTENT (OUT) :: bcof(:,0:,:)!(nobd,0:dimension%lmd,nat_l)
    COMPLEX, INTENT (OUT) :: ccof(-atoms%llod:,:,:,:)!(-llod:llod,nobd,atoms%nlod,nat_l)
    REAL,    OPTIONAL, INTENT (IN) :: eig(:)!(dimension%neigd)
    !     ..
    !     .. Local Scalars ..
    COMPLEX cexp,phase,c_0,c_1,c_2,ci
    REAL const,df,r1,s,tmk,wronk
    REAL s2h, s2h_e(ne)
    INTEGER i,j,k,l,ll1,lm,n,nap,natom,nn,iatom,jatom,lmp,m,nkvec
    INTEGER inv_f,ie,ilo,kspin,iintsp,nintsp,nvmax,lo,inap,natom_l
    LOGICAL l_force
    !     ..
    !     .. Local Arrays ..
    INTEGER nbasf0(atoms%nlod,atoms%nat)
    REAL dfj(0:atoms%lmaxd),fj(0:atoms%lmaxd),fg(3),fgp(3),fgr(3),fk(3),fkp(3),fkr(3)
    REAL alo1(atoms%nlod),blo1(atoms%nlod),clo1(atoms%nlod)
    COMPLEX ylm( (atoms%lmaxd+1)**2 )
    COMPLEX ccchi(2,2)
    LOGICAL apw(0:atoms%lmaxd,atoms%ntype)
    REAL,    ALLOCATABLE :: work_r(:)
    COMPLEX, ALLOCATABLE :: work_c(:)

    CALL timestart("abcof")

    IF (zmat%l_real) THEN
       IF (noco%l_soc.AND.sym%invs) CALL judft_error("BUG in abcof, SOC&INVS but real?")
       IF (noco%l_noco) CALL judft_error("BUG in abcof, l_noco but real?")
    ENDIF

    ci = CMPLX(0.0,1.0)
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
                IF (atoms%invsat(natom).EQ.2) THEN
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
                !$OMP PARALLEL DO &
                !$OMP& DEFAULT(none)&
                !$OMP& PRIVATE(k,i,work_r,work_c,ccchi,kspin,fg,fk,s,r1,fj,dfj,l,df,wronk,tmk,phase,lo,nkvec,&
                !$OMP& inap,nap,j,fgr,fgp,s2h,s2h_e,fkr,fkp,ylm,ll1,m,c_0,c_1,c_2,lmp,inv_f,lm)&
                !$OMP& SHARED(n,nn,natom,natom_l,noco,atoms,sym,cell,oneD,lapw,nvmax,ne,zMat,usdus,ci,iintsp,eig,l_force,&
                !$OMP& alo1,blo1,clo1,jatom,jspin,apw,const,nbasf0,acof,bcof,ccof,force,nat_start,nat_stop)
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
                         dfj(l) = 0.0d0
                      ELSE
                         dfj(l) = const* (usdus%dus(l,n,jspin)*fj(l)-df*usdus%us(l,n,jspin))/wronk
                         fj(l) = const* (df*usdus%uds(l,n,jspin)-fj(l)*usdus%duds(l,n,jspin))/wronk
                      ENDIF
                   ENDDO ! loop over l
                   tmk = tpi_const* (fk(1)*atoms%taual(1,jatom)+&
                                     fk(2)*atoms%taual(2,jatom)+&
                                     fk(3)*atoms%taual(3,jatom))
                   phase = CMPLX(COS(tmk),SIN(tmk))
                   IF (oneD%odi%d1) THEN
                      inap = oneD%ods%ngopr(jatom)
                   ELSE
                      nap = atoms%ngopr(jatom)
                      inap = sym%invtab(nap)
                   END IF
                   DO j = 1,3
                      fkr(j) = 0.0
                      fgr(j) = 0.0
                      DO i = 1,3
                         IF (oneD%odi%d1) THEN
                            fkr(j) = fkr(j) + fk(i)*oneD%ods%mrot(i,j,inap)
                            fgr(j) = fgr(j) + fg(i)*oneD%ods%mrot(i,j,inap)
                         ELSE
                            fkr(j) = fkr(j) + fk(i)*sym%mrot(i,j,inap)
                            fgr(j) = fgr(j) + fg(i)*sym%mrot(i,j,inap)
                         END IF
                      END DO
                   END DO
                   fkp = MATMUL(fkr,cell%bmat)
                   fgp = MATMUL(fgr,cell%bmat)
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
                         IF (atoms%invsat(natom).EQ.2) THEN
                           lmp = ll1 - m
                           inv_f = (-1)**(l-m)
                           c_1 =  conjg(c_1) * inv_f
                           c_2 =  conjg(c_2) * inv_f
                           IF (zmat%l_real) THEN
                              acof(:ne,lmp,natom_l) = acof(:ne,lmp,natom_l) +  c_1 * work_r(:ne)
                              bcof(:ne,lmp,natom_l) = bcof(:ne,lmp,natom_l) +  c_2 * work_r(:ne)
                           ELSE
                              acof(:ne,lmp,natom_l) = acof(:ne,lmp,natom_l) +  c_1 * work_c(:ne)
                              bcof(:ne,lmp,natom_l) = bcof(:ne,lmp,natom_l) +  c_2 * work_c(:ne)
                           END IF
                         ELSE
                            !     ----> loop over bands
                            IF (zmat%l_real) THEN
                               acof(:ne,lm,natom_l) = acof(:ne,lm,natom_l) +  c_1 * work_r(:ne)
                               bcof(:ne,lm,natom_l) = bcof(:ne,lm,natom_l) +  c_2 * work_r(:ne)
                            ELSE
                               acof(:ne,lm,natom_l) = acof(:ne,lm,natom_l) +  c_1 * work_c(:ne)
                               bcof(:ne,lm,natom_l) = bcof(:ne,lm,natom_l) +  c_2 * work_c(:ne)
                            END IF
                         ENDIF
                      ENDDO ! loop over m
                   ENDDO ! loop over l
                   DO lo=1,atoms%nlo(n)
                      DO nkvec=1,lapw%nkvec(lo,jatom) 
                         IF (k==lapw%kvec(nkvec,lo,jatom)) THEN !check if this k-vector has LO attached
                            CALL abclocdn_soc(atoms,sym,noco,lapw,cell,ccchi(:,jspin),iintsp,phase,ylm,&
                                          n,natom,natom_l,k,nkvec,lo,ne,alo1,blo1,clo1,acof,bcof,ccof,zMat,l_force,fgp,force)
                         ENDIF
                      ENDDO
                   END DO
                ENDDO ! loop over LAPWs (k)
#ifndef CPP_OLDINTEL
                !$OMP END PARALLEL DO
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
