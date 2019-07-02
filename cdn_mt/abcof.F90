MODULE m_abcof
CONTAINS
  SUBROUTINE abcof(input,atoms,sym, cell,lapw,ne,usdus,&
                   noco,jspin,oneD, acof,bcof,ccof,zMat,eig,force)
    !     ************************************************************
    !     subroutine constructs the a,b coefficients of the linearized
    !     m.t. wavefunctions for each band and atom.       c.l. fu
    !     ************************************************************
#include "cpp_double.h"

    USE m_constants, ONLY : tpi_const, ImagUnit
    USE m_setabc1lo
    USE m_sphbes
    USE m_dsphbs
    USE m_abclocdn
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
    INTEGER, INTENT (IN) :: ne
    INTEGER, INTENT (IN) :: jspin
    !     ..
    !     .. Array Arguments ..
    COMPLEX, INTENT (OUT) :: acof(:,0:,:)!(nobd,0:dimension%lmd,atoms%nat)
    COMPLEX, INTENT (OUT) :: bcof(:,0:,:)!(nobd,0:dimension%lmd,atoms%nat)
    COMPLEX, INTENT (OUT) :: ccof(-atoms%llod:,:,:,:)!(-llod:llod,nobd,atoms%nlod,atoms%nat)
    REAL,    OPTIONAL, INTENT (IN) :: eig(:)!(dimension%neigd)
    !     ..
    !     .. Local Scalars ..
    COMPLEX cexp,phase,c_0,c_1,c_2
    REAL const,df,r1,s,tmk,wronk,qss(3)
    REAL s2h, s2h_e(ne)
    INTEGER i,j,k,l,ll1,lm ,n,nap,natom,nn,iatom,jatom,lmp,m,nkvec
    INTEGER inv_f,ie,ilo,kspin,iintsp,nintsp,nvmax,lo,inap
    LOGICAL l_force
    !     ..
    !     .. Local Arrays ..
    INTEGER nbasf0(atoms%nlod,atoms%nat)
    REAL dfj(0:atoms%lmaxd),fj(0:atoms%lmaxd),fg(3),fgp(3),fgr(3),fk(3),fkp(3),fkr(3)
    REAL alo1(atoms%nlod),blo1(atoms%nlod),clo1(atoms%nlod)
    COMPLEX ylm( (atoms%lmaxd+1)**2 )
    COMPLEX ccchi(2,2)
    LOGICAL enough(atoms%nat),apw(0:atoms%lmaxd,atoms%ntype)
    REAL,    ALLOCATABLE :: work_r(:)
    COMPLEX, ALLOCATABLE :: work_c(:)

    CALL timestart("abcof")

    IF (zmat%l_real) THEN
       IF (noco%l_soc.AND.sym%invs) CALL judft_error("BUG in abcof, SOC&INVS but real?")
       IF (noco%l_noco) CALL judft_error("BUG in abcof, l_noco but real?")
    ENDIF

    const = 2 * tpi_const/SQRT(cell%omtil)

    acof(:,:,:)   = CMPLX(0.0,0.0)
    bcof(:,:,:)   = CMPLX(0.0,0.0)
    ccof(:,:,:,:) = CMPLX(0.0,0.0)
    l_force = .FALSE.
    IF(PRESENT(eig).AND.input%l_f) THEN
       l_force = .TRUE.
    END IF
    IF(l_force) THEN
       force%acoflo  = CMPLX(0.0,0.0)
       force%bcoflo  = CMPLX(0.0,0.0)
       force%e1cof   = CMPLX(0.0,0.0)
       force%e2cof   = CMPLX(0.0,0.0)
       force%aveccof = CMPLX(0.0,0.0)
       force%bveccof = CMPLX(0.0,0.0)
       force%cveccof = CMPLX(0.0,0.0)
    END IF

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
    IF (noco%l_ss) nintsp = 2
    !---> loop over the interstitial spin
    DO iintsp = 1,nintsp
       !
       nvmax=lapw%nv(jspin)
       IF (noco%l_ss) nvmax=lapw%nv(iintsp)
       IF (iintsp .EQ. 1) THEN
          qss= - noco%qss/2
       ELSE
          qss= + noco%qss/2
       ENDIF

       !---> loop over atom types
       !$OMP PARALLEL DO &
       !$OMP& DEFAULT(none)&
       !$OMP& PRIVATE(n,nn,natom,k,i,work_r,work_c,ccchi,kspin,fg,fk,s,r1,fj,dfj,l,df,wronk,tmk,phase,lo,nkvec,&
       !$OMP& alo1,blo1,clo1,inap,nap,j,fgr,fgp,s2h,s2h_e,fkr,fkp,ylm,ll1,m,c_0,c_1,c_2,jatom,lmp,inv_f,lm)&
       !$OMP& SHARED(noco,atoms,sym,cell,oneD,lapw,nvmax,ne,zMat,usdus,iintsp,eig,l_force,&
       !$OMP& jspin,qss,apw,const,nbasf0,enough,acof,bcof,ccof,force)
       DO n = 1,atoms%ntype
          CALL setabc1lo(atoms,n,usdus,jspin,alo1,blo1,clo1)
          
          !  ----> loop over equivalent atoms
          DO nn = 1,atoms%neq(n)
             natom = 0
             DO i = 1, n-1
                natom = natom + atoms%neq(i)
             ENDDO
             natom = natom + nn
             IF ((sym%invsat(natom).EQ.0) .OR. (sym%invsat(natom).EQ.1)) THEN
                !--->        loop over lapws
                IF (zmat%l_real) THEN
                   ALLOCATE ( work_r(ne) )
                ELSE
                   ALLOCATE ( work_c(ne) )
                ENDIF
                DO k = 1,nvmax
                   IF (.NOT.noco%l_noco) THEN
                      IF (zmat%l_real) THEN
                         work_r(:ne)=zMat%data_r(k,:ne)
                      ELSE
                         work_c(:ne)=zMat%data_c(k,:ne)
                      END IF
                   ENDIF

                   IF (noco%l_noco) THEN
                      !--->            generate the complex conjgates of the spinors (chi)
                      ccchi(1,1) = CONJG( EXP(-ImagUnit*noco%alph(n)/2)*COS(noco%beta(n)/2))
                      ccchi(1,2) = CONJG(-EXP(-ImagUnit*noco%alph(n)/2)*SIN(noco%beta(n)/2))
                      ccchi(2,1) = CONJG( EXP( ImagUnit*noco%alph(n)/2)*SIN(noco%beta(n)/2))
                      ccchi(2,2) = CONJG( EXP( ImagUnit*noco%alph(n)/2)*COS(noco%beta(n)/2))
                      IF (noco%l_ss) THEN
                         !--->              the coefficients of the spin-down basis functions are
                         !--->              stored in the second half of the eigenvector
                         kspin = (iintsp-1)*(lapw%nv(1)+atoms%nlotot)
                         DO i = 1,ne
                            work_c(i) = ccchi(iintsp,jspin)*zMat%data_c(kspin+k,i)
                         ENDDO
                      ELSE
                         !--->              perform sum over the two interstitial spin directions
                         !--->              and take into account the spin boundary conditions
                         !--->              (jspin counts the local spin directions inside each MT)
                         kspin = lapw%nv(1)+atoms%nlotot
                         DO i = 1,ne
                            work_c(i) = ccchi(1,jspin)*zMat%data_c(k,i)&
                                 &                        + ccchi(2,jspin)*zMat%data_c(kspin+k,i)
                         ENDDO
                      ENDIF
                   ENDIF ! (noco%l_noco)
                   IF (noco%l_ss) THEN
                      fg = lapw%gvec(:,k,iintsp) + qss
                   ELSE
                      fg = lapw%gvec(:,k,jspin) + qss
                   ENDIF ! (noco%l_ss)
                   fk = lapw%bkpt + fg
                   s =  DOT_PRODUCT(fk,MATMUL(cell%bbmat,fk))
                   IF(l_force) THEN
                      s2h = 0.5*s
                      s2h_e(:ne) = s2h-eig(:ne)
                   END IF
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
                   tmk = tpi_const* (fk(1)*atoms%taual(1,natom)+&
                                     fk(2)*atoms%taual(2,natom)+&
                                     fk(3)*atoms%taual(3,natom))
                   phase = CMPLX(COS(tmk),SIN(tmk))
                   IF (oneD%odi%d1) THEN
                      inap = oneD%ods%ngopr(natom)
                   ELSE
                      nap = sym%ngopr(natom)
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
                         !     ----> loop over bands
                         IF (zmat%l_real) THEN
                            acof(:ne,lm,natom) = acof(:ne,lm,natom) +  c_1 * work_r(:ne)
                            bcof(:ne,lm,natom) = bcof(:ne,lm,natom) +  c_2 * work_r(:ne)
                         ELSE
                            acof(:ne,lm,natom) = acof(:ne,lm,natom) +  c_1 * work_c(:ne)
                            bcof(:ne,lm,natom) = bcof(:ne,lm,natom) +  c_2 * work_c(:ne)
                         END IF

                         IF (atoms%l_geo(n).AND.l_force) THEN
                            IF (zmat%l_real) THEN
                               force%e1cof(:ne,lm,natom) = force%e1cof(:ne,lm,natom) + c_1 * work_r(:ne) * s2h_e(:ne)
                               force%e2cof(:ne,lm,natom) = force%e2cof(:ne,lm,natom) + c_2 * work_r(:ne) * s2h_e(:ne)
                               DO i = 1,3
                                  force%aveccof(i,:ne,lm,natom) = force%aveccof(i,:ne,lm,natom) + c_1 * work_r(:ne) * fgp(i)
                                  force%bveccof(i,:ne,lm,natom) = force%bveccof(i,:ne,lm,natom) + c_2 * work_r(:ne) * fgp(i)
                               END DO
                            ELSE
                               force%e1cof(:ne,lm,natom) = force%e1cof(:ne,lm,natom) + c_1 * work_c(:ne) * s2h_e(:ne)
                               force%e2cof(:ne,lm,natom) = force%e2cof(:ne,lm,natom) + c_2 * work_c(:ne) * s2h_e(:ne)
                               DO i = 1,3
                                  force%aveccof(i,:ne,lm,natom) = force%aveccof(i,:ne,lm,natom) + c_1 * work_c(:ne) * fgp(i)
                                  force%bveccof(i,:ne,lm,natom) = force%bveccof(i,:ne,lm,natom) + c_2 * work_c(:ne) * fgp(i)
                               END DO
                            END IF
                         ENDIF

                         IF (noco%l_soc.AND.sym%invs) THEN
                            IF (sym%invsat(natom).EQ.1) THEN
                               jatom = sym%invsatnr(natom)
                               lmp = ll1 - m
                               inv_f = (-1)**(l-m)
                               c_1 =  CONJG(c_1) * inv_f
                               c_2 =  CONJG(c_2) * inv_f
                               CALL CPP_BLAS_caxpy(ne,c_1,work_c,1, acof(1,lmp,jatom),1)
                               CALL CPP_BLAS_caxpy(ne,c_2,work_c,1, bcof(1,lmp,jatom),1)
                               IF (atoms%l_geo(n).AND.l_force) THEN
                                  CALL CPP_BLAS_caxpy(ne,c_1,work_c*s2h_e,1, force%e1cof(1,lmp,jatom),1)
                                  CALL CPP_BLAS_caxpy(ne,c_2,work_c*s2h_e,1, force%e2cof(1,lmp,jatom),1)
                                  DO i = 1,3
                                     CALL CPP_BLAS_caxpy(ne,c_1,work_c*fgp(i),1, force%aveccof(i,1,lmp,jatom),3)
                                     CALL CPP_BLAS_caxpy(ne,c_2,work_c*fgp(i),1, force%bveccof(i,1,lmp,jatom),3)
                                  END DO
                               END IF
                            ENDIF
                         ENDIF
                      ENDDO ! loop over m
                   ENDDO ! loop over l
                   DO lo=1,atoms%nlo(n)
                      DO nkvec=1,lapw%nkvec(lo,natom)
                         IF (k==lapw%kvec(nkvec,lo,natom)) THEN !check if this k-vector has LO attached
                            CALL abclocdn(atoms,sym,noco,lapw,cell,ccchi(:,jspin),iintsp,phase,ylm,&
                                          n,natom,k,nkvec,lo,ne,alo1,blo1,clo1,acof,bcof,ccof,zMat,l_force,fgp,force)
                         ENDIF
                      ENDDO
                   END DO
                ENDDO ! loop over LAPWs
                IF (zmat%l_real) THEN
                   DEALLOCATE(work_r)
                ELSE               
                   DEALLOCATE(work_c)
                ENDIF
             ENDIF  ! invsatom == ( 0 v 1 )
          ENDDO    ! loop over equivalent atoms
       ENDDO       ! loop over atom types
       !$OMP END PARALLEL DO
    ENDDO       ! loop over interstitial spin
    IF (noco%l_soc.AND.sym%invs) THEN

       !
       !                           -p,n       (l+m)   p,n  *
       ! Usually, we exploit that A     = (-1)      (A    )  if p and -p are the positions
       !                           l,m                l,-m
       ! of two atoms related by inversion symmetry and the coefficients are considered to
       ! be in the local frame of the representative atom. This is possible, if z is real.
       ! After SOC, however, the eigenvectors z are complex and this is no longer possible
       ! so the z has to enter, not z*. This is done within the k-loop.
       !                                    -p,n       m   p,n  *
       ! When called from hsohelp, we need A     = (-1)  (A    ) because we don't have to
       !                                     l,m           l,-m                    l
       ! rotate, but in the sums in hsoham only products A*  A   enter and the (-1) cancels.
       !                                                  lm  lm
    ELSE
       iatom = 0
       DO n = 1,atoms%ntype
          DO nn = 1,atoms%neq(n)
             iatom = iatom + 1
             IF (sym%invsat(iatom).EQ.1) THEN
                jatom = sym%invsatnr(iatom)
                cexp = EXP(tpi_const*ImagUnit*DOT_PRODUCT(atoms%taual(:,jatom)&
                     &             + atoms%taual(:,iatom),lapw%bkpt))
                DO ilo = 1,atoms%nlo(n)
                   l = atoms%llo(ilo,n)
                   DO m = -l,l
                      inv_f = (-1.0)**(m+l)
                      DO ie = 1,ne
                         ccof(m,ie,ilo,jatom) = inv_f * cexp * CONJG(  ccof(-m,ie,ilo,iatom))
                         IF(l_force) THEN
                            force%acoflo(m,ie,ilo,jatom) = inv_f * cexp * CONJG(force%acoflo(-m,ie,ilo,iatom))
                            force%bcoflo(m,ie,ilo,jatom) = inv_f * cexp * CONJG(force%bcoflo(-m,ie,ilo,iatom))
                            DO i = 1,3
                               force%cveccof(i,m,ie,ilo,jatom) = -inv_f * cexp * CONJG(force%cveccof(i,-m,ie,ilo,iatom))
                            END DO
                         END IF
                      ENDDO
                   ENDDO
                ENDDO
                DO l = 0,atoms%lmax(n)
                   ll1 = l* (l+1)
                   DO m =-l,l
                      lm  = ll1 + m
                      lmp = ll1 - m
                      inv_f = (-1.0)**(m+l)
                      DO ie = 1,ne
                         acof(ie,lm,jatom) = inv_f * cexp * CONJG(acof(ie,lmp,iatom))
                         bcof(ie,lm,jatom) = inv_f * cexp * CONJG(bcof(ie,lmp,iatom))
                      END DO
                      IF (atoms%l_geo(n).AND.l_force) THEN
                         DO ie = 1,ne
                            force%e1cof(ie,lm,jatom) = inv_f * cexp * CONJG(force%e1cof(ie,lmp,iatom))
                            force%e2cof(ie,lm,jatom) = inv_f * cexp * CONJG(force%e2cof(ie,lmp,iatom))
                            DO i = 1,3
                               force%aveccof(i,ie,lm,jatom) = -inv_f * cexp * CONJG(force%aveccof(i,ie,lmp,iatom))
                               force%bveccof(i,ie,lm,jatom) = -inv_f * cexp * CONJG(force%bveccof(i,ie,lmp,iatom))
                            END DO
                         END DO
                      END IF
                   ENDDO
                ENDDO
             ENDIF
          ENDDO
       ENDDO
    ENDIF

    CALL timestop("abcof")

  END SUBROUTINE abcof
END MODULE m_abcof
