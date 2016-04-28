MODULE m_abcof
CONTAINS
  SUBROUTINE abcof(atoms,nobd,sym, cell, bkpt,lapw,ne,z,usdus,&
       noco,jspin,kveclo,oneD, acof,bcof,ccof)
    !     ************************************************************
    !     subroutine constructs the a,b coefficients of the linearized
    !     m.t. wavefunctions for each band and atom.       c.l. fu
    !     ************************************************************
#include "cpp_double.h"

    USE m_constants, ONLY : tpi_const
    USE m_setabc1locdn
    USE m_sphbes
    USE m_dsphbs
    USE m_abclocdn
    USE m_ylm
    USE m_types
    IMPLICIT NONE
    TYPE(t_usdus),INTENT(IN)   :: usdus
    TYPE(t_lapw),INTENT(IN)   :: lapw
    TYPE(t_oneD),INTENT(IN)   :: oneD
    TYPE(t_noco),INTENT(IN)   :: noco
    TYPE(t_sym),INTENT(IN)    :: sym
    TYPE(t_cell),INTENT(IN)   :: cell
    TYPE(t_atoms),INTENT(IN)  :: atoms
    !     ..
    !     .. Scalar Arguments ..
    INTEGER, INTENT (IN) :: nobd
    INTEGER, INTENT (IN) :: ne
    INTEGER, INTENT (IN) :: jspin
    !     ..
    !     .. Array Arguments ..
    INTEGER, INTENT (IN) :: kveclo(atoms%nlotot)
    REAL,    INTENT (IN) :: bkpt(3)
#if ( !defined(CPP_INVERSION) || defined(CPP_SOC) )
    COMPLEX, INTENT (IN) :: z(:,:)!(dimension%nbasfcn,dimension%neigd)
#else
    REAL,    INTENT (IN) :: z(:,:)!(dimension%nbasfcn,dimension%neigd)
#endif
    COMPLEX, INTENT (OUT):: acof(:,0:,:)!(nobd,0:dimension%lmd,atoms%natd)
    COMPLEX, INTENT (OUT):: bcof(:,0:,:)!(nobd,0:dimension%lmd,atoms%natd)
    COMPLEX, INTENT (OUT):: ccof(-atoms%llod:,:,:,:)!(-llod:llod,nobd,atoms%nlod,atoms%natd)
    !     ..
    !     .. Local Scalars ..
    COMPLEX cexp,phase,c_0,c_1,c_2,ci
    REAL const,df,r1,s,tmk,wronk,qss1,qss2,qss3
    INTEGER i,j,k,l,ll1,lm ,n,nap,natom,nn,iatom,jatom,lmp,m
    INTEGER inv_f,ie,ilo,kspin,iintsp,nintsp,nvmax,lo,inap
    !     ..
    !     .. Local Arrays ..
    INTEGER kvec(2*(2*atoms%llod+1),atoms%nlod,atoms%natd  )
    INTEGER nbasf0(atoms%nlod,atoms%natd),nkvec(atoms%nlod,atoms%natd)
    REAL dfj(0:atoms%lmaxd),fj(0:atoms%lmaxd),fk(3),fkp(3),fkr(3)
    REAL alo1(atoms%nlod,atoms%ntypd),blo1(atoms%nlod,atoms%ntypd),clo1(atoms%nlod,atoms%ntypd)
    COMPLEX ylm( (atoms%lmaxd+1)**2 )
    COMPLEX ccchi(2,2)
    !$    COMPLEX, ALLOCATABLE :: acof_loc(:,:), bcof_loc(:,:)
    !$    COMPLEX, ALLOCATABLE :: acof_inv(:,:), bcof_inv(:,:)
    LOGICAL enough(atoms%natd),apw(0:atoms%lmaxd,atoms%ntypd)
#if ( !defined(CPP_INVERSION) || defined(CPP_SOC) )
    COMPLEX, ALLOCATABLE :: work(:)
#else
    REAL,    ALLOCATABLE :: work(:)
#endif
    !     ..
    !     .. Intrinsic Functions ..
    INTRINSIC cmplx,conjg,exp,sqrt
    !     ..
    ci = cmplx(0.0,1.0)
    const = 2 * tpi_const/sqrt(cell%omtil)
    !
    acof(:,:,:) = cmplx(0.0,0.0)
    bcof(:,:,:) = cmplx(0.0,0.0)
    !     ..
    !+APW_LO
    DO n = 1, atoms%ntype
       DO l = 0,atoms%lmax(n)
          apw(l,n) = .false.
          DO lo = 1,atoms%nlo(n)
             IF (atoms%l_dulo(lo,n)) apw(l,n) = .true.
          ENDDO
#ifdef CPP_APW
          IF (atoms%lapw_l(n).GE.l) apw(l,n) = .false.
#endif
       ENDDO
       DO lo = 1,atoms%nlo(n)
          IF (atoms%l_dulo(lo,n)) apw(atoms%llo(lo,n),n) = .true.
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
       CALL setabc1locdn(jspin, atoms,lapw,ne,noco,iintsp, sym,usdus,&
            kveclo, enough,nkvec,kvec,nbasf0,ccof, alo1,blo1,clo1)
       !
       IF (iintsp .EQ. 1) THEN
          qss1= - noco%qss(1)/2
          qss2= - noco%qss(2)/2
          qss3= - noco%qss(3)/2
       ELSE
          qss1= + noco%qss(1)/2
          qss2= + noco%qss(2)/2
          qss3= + noco%qss(3)/2
       ENDIF

       !---> loop over atom types
       natom = 0
       DO n = 1,atoms%ntype
          !  ----> loop over equivalent atoms
          DO nn = 1,atoms%neq(n)
             natom = natom + 1
             IF ((atoms%invsat(natom).EQ.0) .OR. (atoms%invsat(natom).EQ.1)) THEN
                !--->        loop over lapws
                !$OMP PARALLEL IF(enough(natom)) &
                !$OMP& DEFAULT(none)&
                !$OMP& PRIVATE(k,i,work,ccchi,kspin,fk,s,r1,fj,dfj,l,df,wronk,tmk,phase,&
                !$OMP& inap,nap,j,fkr,fkp,ylm,ll1,m,c_0,c_1,c_2,jatom,lmp,inv_f,lm,&
                !$OMP& acof_loc,bcof_loc,acof_inv,bcof_inv)&
                !$OMP& SHARED(noco,atoms,sym,cell,oneD,lapw,nvmax,ne,z,usdus,n,ci,iintsp,&
                !$OMP& jspin,bkpt,qss1,qss2,qss3,&
                !$OMP& apw,const,natom,&
                !$OMP& nobd,&
                !$OMP& alo1,blo1,clo1,kvec,nbasf0,nkvec,enough,acof,bcof)&
                !$OMP& REDUCTION(+:ccof)
                ALLOCATE ( work(nobd) )

                !$    ALLOCATE(acof_loc(nobd,0:size(acof,2)-1),bcof_loc(nobd,0:size(acof,2)-1))
                !$    acof_loc(:,:) = cmplx(0.0,0.0)
                !$    bcof_loc(:,:) = cmplx(0.0,0.0)

#if ( defined(CPP_SOC) && defined(CPP_INVERSION) )
                !$ IF (invsat(natom).EQ.1) THEN
                !$    ALLOCATE(acof_inv(nobd,0:lmd),bcof_inv(nobd,0:size(acof,2)-1))
                !$    acof_inv(:,:) = cmplx(0.0,0.0)
                !$    bcof_inv(:,:) = cmplx(0.0,0.0)
                !$ ENDIF
#endif


!!!!

                !$OMP  DO
                DO k = 1,nvmax
                   IF (.NOT.noco%l_noco) THEN
                      DO i = 1,ne
                         work(i) = z(k,i)
                      ENDDO
                   ENDIF

                   IF (noco%l_noco) THEN
                      !--->            generate the complex conjgates of the spinors (chi)
                      ccchi(1,1) = conjg( exp(-ci*noco%alph(n)/2)*cos(noco%beta(n)/2))
                      ccchi(1,2) = conjg(-exp(-ci*noco%alph(n)/2)*sin(noco%beta(n)/2))
                      ccchi(2,1) = conjg( exp( ci*noco%alph(n)/2)*sin(noco%beta(n)/2))
                      ccchi(2,2) = conjg( exp( ci*noco%alph(n)/2)*cos(noco%beta(n)/2))
                      IF (noco%l_ss) THEN
                         !--->              the coefficients of the spin-down basis functions are
                         !--->              stored in the second half of the eigenvector
                         kspin = (iintsp-1)*(lapw%nv(1)+atoms%nlotot)
                         DO i = 1,ne
                            work(i) = ccchi(iintsp,jspin)*z(kspin+k,i)
                         ENDDO
                      ELSE
                         !--->              perform sum over the two interstitial spin directions
                         !--->              and take into account the spin boundary conditions
                         !--->              (jspin counts the local spin directions inside each MT)
                         kspin = lapw%nv(1)+atoms%nlotot
                         DO i = 1,ne
                            work(i) = ccchi(1,jspin)*z(k,i)&
                                 &                        + ccchi(2,jspin)*z(kspin+k,i)
                         ENDDO
                      ENDIF
                   ENDIF ! (noco%l_noco)
                   IF (noco%l_ss) THEN
                      fk(1) = bkpt(1) + lapw%k1(k,iintsp) + qss1
                      fk(2) = bkpt(2) + lapw%k2(k,iintsp) + qss2
                      fk(3) = bkpt(3) + lapw%k3(k,iintsp) + qss3
                   ELSE
                      fk(1) = bkpt(1) + lapw%k1(k,jspin) + qss1
                      fk(2) = bkpt(2) + lapw%k2(k,jspin) + qss2
                      fk(3) = bkpt(3) + lapw%k3(k,jspin) + qss3
                   ENDIF ! (noco%l_ss)
                   s=  dot_product(fk,matmul(cell%bbmat,fk))
                   s = sqrt(s)
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
                        &                     fk(2)*atoms%taual(2,natom)+&
                        &                     fk(3)*atoms%taual(3,natom))
                   phase = cmplx(cos(tmk),sin(tmk))
                   IF (oneD%odi%d1) THEN
                      inap = oneD%ods%ngopr(natom)
                   ELSE
                      nap = atoms%ngopr(natom)
                      inap = sym%invtab(nap)
                   END IF
                   DO j = 1,3
                      fkr(j) = 0.
                      DO i = 1,3
                         IF (oneD%odi%d1) THEN
                            fkr(j) = fkr(j) + fk(i)*oneD%ods%mrot(i,j,inap)
                         ELSE
                            fkr(j) = fkr(j) + fk(i)*sym%mrot(i,j,inap)
                         END IF
                      ENDDO
                   ENDDO
                   fkp=matmul(fkr,cell%bmat)
                   !     ----> generate spherical harmonics
                   CALL ylm4(atoms%lmax(n),fkp,ylm)
                   !     ----> loop over l
                   DO l = 0,atoms%lmax(n)
                      ll1 = l* (l+1)
                      !     ----> loop over m
                      DO m = -l,l
                         lm = ll1 + m
                         c_0 = conjg(ylm(lm+1))*phase
                         c_1 = c_0 *  fj(l)
                         c_2 = c_0 * dfj(l)
                         !     ----> loop over bands
                         !$                 if (.false.) THEN
                         DO i = 1,ne
                            acof(i,lm,natom) = acof(i,lm,natom) + &
                                 &                                  c_1 * work(i)
                         ENDDO
                         DO i = 1,ne
                            bcof(i,lm,natom) = bcof(i,lm,natom) +&
                                 &                                  c_2 * work(i)
                         ENDDO
                         !$                 endif
                         !$                 DO i = 1,ne
                         !$                   acof_loc(i,lm) = acof_loc(i,lm) + c_1 * work(i)
                         !$                 ENDDO
                         !$                 DO i = 1,ne
                         !$                   bcof_loc(i,lm) = bcof_loc(i,lm) + c_2 * work(i)
                         !$                 ENDDO
#if ( defined(CPP_SOC) && defined(CPP_INVERSION) )
                         IF (atoms%invsat(natom).EQ.1) THEN
                            jatom = sym%invsatnr(natom)
                            lmp = ll1 - m
                            inv_f = (-1)**(l-m)
                            c_1 =  conjg(c_1) * inv_f
                            c_2 =  conjg(c_2) * inv_f
                            !$                   if (.false.) THEN
                            CALL CPP_BLAS_caxpy(ne,c_1,work,1,&
                                 &                                   acof(1,lmp,jatom),1)
                            CALL CPP_BLAS_caxpy(ne,c_2,work,1,&
                                 &                                   bcof(1,lmp,jatom),1)
                            !$                   endif
                            !$                   CALL CPP_BLAS_caxpy(ne,c_1,work,1,&
                            !$                                       acof_inv(1,lmp),1)
                            !$                   CALL CPP_BLAS_caxpy(ne,c_2,work,1,&
                            !$                                       bcof_inv(1,lmp),1)
                         ENDIF
#endif
                      ENDDO ! loop over m
                   ENDDO ! loop over l
                   IF (.NOT.enough(natom)) THEN
                      CALL abclocdn(atoms,sym, noco,ccchi(1,jspin),kspin,iintsp,const,phase,ylm,n,natom,k,&
                           s,nvmax,ne,z,nbasf0,alo1,blo1,clo1,kvec(1,1,natom),nkvec,enough,acof,bcof,ccof)
                   ENDIF
                ENDDO ! loop over LAPWs
                !$OMP END DO
                !$OMP CRITICAL
                !$      acof(:,:,natom) = acof(:,:,natom) + acof_loc(:,:)
                !$      bcof(:,:,natom) = bcof(:,:,natom) + bcof_loc(:,:)
#if ( defined(CPP_SOC) && defined(CPP_INVERSION) )
                !$      IF (invsat(natom).EQ.1) THEN
                !$        jatom = invsatnr(natom)
                !$        acof(:,:,jatom) = acof(:,:,jatom) + acof_inv(:,:)
                !$        bcof(:,:,jatom) = bcof(:,:,jatom) + bcof_inv(:,:)
                !$      ENDIF
#endif
                !$OMP END CRITICAL
                !$    DEALLOCATE(acof_loc,bcof_loc)
#if ( defined(CPP_SOC) && defined(CPP_INVERSION) )
                !$    DEALLOCATE(acof_inv,bcof_inv)
#endif
                DEALLOCATE(work)
                !$OMP END PARALLEL
             ENDIF  ! invsatom == ( 0 v 1 )
          ENDDO    ! loop over equivalent atoms
       ENDDO       ! loop over atom types
    ENDDO       ! loop over interstitial spin

#if ( defined(CPP_SOC) && defined(CPP_INVERSION) )
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
#else
    iatom = 0
    DO n = 1,atoms%ntype
       DO nn = 1,atoms%neq(n)
          iatom = iatom + 1
          IF (atoms%invsat(iatom).EQ.1) THEN
             jatom = sym%invsatnr(iatom)
             cexp = exp(tpi_const*ci*dot_product(atoms%taual(:,jatom)&
                  &             + atoms%taual(:,iatom),(/bkpt(1),bkpt(2),bkpt(3)/)))
             DO ilo = 1,atoms%nlo(n)
                l = atoms%llo(ilo,n)
                DO m = -l,l
                   inv_f = (-1.0)**(m+l)
                   DO ie = 1,ne
                      ccof(m,ie,ilo,jatom) = inv_f * cexp *conjg(  ccof(-m,ie,ilo,iatom))
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
                      acof(ie,lm,jatom) = inv_f * cexp * conjg(acof(ie,lmp,iatom))
                   ENDDO
                   DO ie = 1,ne
                      bcof(ie,lm,jatom) = inv_f * cexp * conjg(bcof(ie,lmp,iatom))
                   ENDDO
                ENDDO
             ENDDO
          ENDIF
       ENDDO
    ENDDO
#endif

  END SUBROUTINE abcof
END MODULE m_abcof
