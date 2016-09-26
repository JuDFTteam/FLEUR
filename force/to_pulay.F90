MODULE m_topulay
!     ************************************************************
!     put together all lm quantities for pulay forces
!     al la Yu et al equa A12,A17,A20
!     ************************************************************
CONTAINS
  SUBROUTINE to_pulay(&
       input,atoms,nobd,sym,lapw,noco,cell,bkpt,z,ne,eig,&
       usdus,kveclo,jspin,oneD,&
       acof,bcof,e1cof,e2cof,aveccof,bveccof,&
       ccof,acoflo,bcoflo,cveccof)
    !
    USE m_constants, ONLY : tpi_const
    USE m_setabc1locdn
    USE m_sphbes
    USE m_dsphbs
    USE m_abclocdnpulay
    USE m_ylm
    USE m_types
    IMPLICIT NONE
    TYPE(t_input),INTENT(IN)  :: input
    TYPE(t_oneD),INTENT(IN)   :: oneD
    TYPE(t_noco),INTENT(IN)   :: noco
    TYPE(t_sym),INTENT(IN)    :: sym
    TYPE(t_cell),INTENT(IN)   :: cell
    TYPE(t_atoms),INTENT(IN)  :: atoms
    TYPE(t_usdus),INTENT(IN)  :: usdus
    TYPE(t_lapw),INTENT(IN)   :: lapw
!     ..
!     .. Scalar Arguments ..
      INTEGER, INTENT (IN) :: nobd   
      INTEGER, INTENT (IN) :: ne      
      INTEGER, INTENT (IN) :: jspin 
!     ..
!     .. Array Arguments ..
      INTEGER, INTENT (IN) :: kveclo(atoms%nlotot)
      REAL,    INTENT (IN) :: bkpt(3)  
      REAL,    INTENT (IN) :: eig(:)!(dimension%neigd)  
#if ( defined(CPP_INVERSION) && !defined(CPP_SOC) )
      REAL,    INTENT (IN) :: z(:,:)!(dimension%nbasfcn,dimension%neigd)
#else
      COMPLEX, INTENT (IN) :: z(:,:)!(dimension%nbasfcn,dimension%neigd)
#endif
      COMPLEX, INTENT (OUT)::      acof(:,0:,:)!(nobd,0:dimension%lmd,atoms%natd)
      COMPLEX, INTENT (OUT)::      bcof(:,0:,:)!(nobd,0:dimension%lmd,atoms%natd)
      COMPLEX, INTENT (OUT)::      ccof(-atoms%llod:atoms%llod,nobd,atoms%nlod,atoms%natd)
      COMPLEX, INTENT (OUT)::    acoflo(-atoms%llod:atoms%llod,nobd,atoms%nlod,atoms%natd)
      COMPLEX, INTENT (OUT)::    bcoflo(-atoms%llod:atoms%llod,nobd,atoms%nlod,atoms%natd)
      COMPLEX, INTENT (OUT)::     e1cof(:,0:,:)!(nobd,0:dimension%lmd,atoms%natd)
      COMPLEX, INTENT (OUT)::     e2cof(:,0:,:)!(nobd,0:dimension%lmd,atoms%natd)
      COMPLEX, INTENT (OUT):: aveccof(:,:,0:,:)!(3,nobd,0:dimension%lmd,atoms%natd)
      COMPLEX, INTENT (OUT):: bveccof(:,:,0:,:)!(3,nobd,0:dimension%lmd,atoms%natd)
      COMPLEX, INTENT (OUT):: cveccof(3,-atoms%llod:atoms%llod,nobd,atoms%nlod,atoms%natd)
!-odim
!+odim
!     ..
!     .. Local Scalars ..
      COMPLEX phase,c_0,c_1,c_2
      REAL const,df,r1,s,tmk,wronk,s2h,s2h_e,qss1,qss2,qss3
      REAL t1,t2,t3
      INTEGER i,ie,j,k,l,lm ,n,natom,nn,ll1,iatom,jatom,lmp,ilo,m
      INTEGER inv_f,lo,nintsp,iintsp,nvmax,kspin,nap,inap
!     ..
!     .. Local Arrays ..
      INTEGER kvec(2*(2*atoms%llod+1) ,atoms%nlod,atoms%natd )
      INTEGER nbasf0(atoms%nlod,atoms%natd),nkvec(atoms%nlod,atoms%natd)
      REAL alo1(atoms%nlod,atoms%ntypd),blo1(atoms%nlod,atoms%ntypd),clo1(atoms%nlod,atoms%ntypd)
      REAL dfj(0:atoms%lmaxd),fg(3),fgp(3),fgr(3),fj(0:atoms%lmaxd),fk(3), fkp(3),fkr(3)
      COMPLEX ylm( (atoms%lmaxd+1)**2 ),ccchi(2,2)
      LOGICAL enough(atoms%natd),apw(0:atoms%lmaxd,atoms%ntypd)
      COMPLEX, ALLOCATABLE :: aaux(:),baux(:)
#if ( defined(CPP_INVERSION) && !defined(CPP_SOC) )
      REAL,    ALLOCATABLE :: work(:)
#else
      COMPLEX, ALLOCATABLE :: work(:)
#endif
!     ..
!     .. Data statements ..
      COMPLEX,PARAMETER:: czero=cmplx(.0,0.0)
      COMPLEX,PARAMETER:: ci = cmplx(0.0,1.0)
!     ..
      ALLOCATE ( aaux(nobd),baux(nobd),work(nobd) )
      const = 2 * tpi_const/sqrt(cell%omtil)
!
!     preset lm quantities of Pulay forces a la yu et al
!
      acof(:,:,:) = czero ;  acoflo(:,:,:,:) = czero
      bcof(:,:,:) = czero ;  bcoflo(:,:,:,:) = czero
      e1cof(:,:,:) = czero ; aveccof(:,:,:,:) = czero
      e2cof(:,:,:) = czero ; bveccof(:,:,:,:) = czero
      cveccof(:,:,:,:,:) = czero
!
      DO n = 1, atoms%ntype
         DO l = 0,atoms%lmax(n)
           apw(l,n) = .false.
           DO lo = 1,atoms%nlo(n)
             IF (atoms%l_dulo(lo,n)) apw(l,n) = .true.
           ENDDO
           IF ((input%l_useapw).AND.(atoms%lapw_l(n).GE.l)) apw(l,n) = .false.

         ENDDO
         DO lo = 1,atoms%nlo(n)
           IF (atoms%l_dulo(lo,n)) apw(atoms%llo(lo,n),n) = .true.
         ENDDO
      ENDDO
!
      nintsp = 1
      IF (noco%l_ss) nintsp = 2
!---> loop over the interstitial spin
      DO iintsp = 1,nintsp
      nvmax = lapw%nv(jspin)
      IF (noco%l_ss) nvmax = lapw%nv(iintsp)
!
      CALL setabc1locdn(jspin,atoms,lapw,ne,noco,iintsp,&
           sym,usdus, kveclo, enough,nkvec,kvec,nbasf0,ccof,&
           alo1,blo1,clo1)
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
!
!     ----> loop over lapws
!
      nvmax = lapw%nv(jspin)
      IF (noco%l_ss) nvmax = lapw%nv(iintsp)
!
      DO  k = 1,nvmax
         IF (.NOT.noco%l_noco) THEN
            DO ie = 1,ne
               work(ie) = z(k,ie)
            ENDDO
         ENDIF
         IF (noco%l_ss) THEN
            fg(1) = lapw%k1(k,iintsp) + qss1
            fg(2) = lapw%k2(k,iintsp) + qss2
            fg(3) = lapw%k3(k,iintsp) + qss3
         ELSE
            fg(1) = lapw%k1(k,jspin) + qss1
            fg(2) = lapw%k2(k,jspin) + qss2
            fg(3) = lapw%k3(k,jspin) + qss3
         ENDIF
         DO i = 1,3
            fk(i) = bkpt(i) + fg(i)
         END DO
!-gu
         s=dot_product(fk,matmul(cell%bbmat,fk))
         s2h = 0.5*s
         s = sqrt(s)
!     ----> loop over atom types
         natom = 0
         DO  n = 1,atoms%ntype
            IF (noco%l_noco) THEN
!--->          generate the complex conjgates of the spinors (chi)
               ccchi(1,1) = conjg( exp(-ci*noco%alph(n)/2)*cos(noco%beta(n)/2))
               ccchi(1,2) = conjg(-exp(-ci*noco%alph(n)/2)*sin(noco%beta(n)/2))
               ccchi(2,1) = conjg( exp(ci*noco%alph(n)/2)*sin(noco%beta(n)/2))
               ccchi(2,2) = conjg( exp(ci*noco%alph(n)/2)*cos(noco%beta(n)/2))
               IF (noco%l_ss) THEN
!--->             the coefficients of the spin-down basis functions are
!--->             stored in the second half of the eigenvector
                  kspin = (iintsp-1)*(lapw%nv(1)+atoms%nlotot)
                  DO i = 1,ne
                     work(i) = ccchi(iintsp,jspin)*z(kspin+k,i)
                  ENDDO
               ELSE
!--->             perform sum over the two interstitial spin directions
!--->             and take into account the spin boundary conditions
!--->             (jspin counts the local spin directions inside each MT)
                  kspin = lapw%nv(1)+atoms%nlotot
                  DO ie = 1,ne
                     work(ie) = ccchi(1,jspin)*z(k,ie)&
     &                        + ccchi(2,jspin)*z(kspin+k,ie)
                  ENDDO
               ENDIF
            ENDIF
            r1 = atoms%rmt(n)*s
            CALL sphbes(atoms%lmax(n),r1, fj)
            CALL dsphbs(atoms%lmax(n),r1,fj, dfj)
!     ----> construct a and b coefficients

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
            END DO

!     ----> loop over equivalent atoms
            DO  nn = 1,atoms%neq(n)
               natom = natom + 1
!-inv
               IF ((atoms%invsat(natom).EQ.0) .OR. (atoms%invsat(natom).EQ.1)) THEN
!
! taual is the actual position (in internal ccordinates) of
! the considered atom. It can be rotated into the representative
! by operation R_equiv
               tmk = tpi_const* (fk(1)*atoms%taual(1,natom)+&
     &                     fk(2)*atoms%taual(2,natom)+&
     &                     fk(3)*atoms%taual(3,natom))
               phase = cmplx(cos(tmk),sin(tmk))

!  ROTATION:  kr = kaux*R_equiv
!      R_equiv rotates equivalent atom into representative
!      Note: For forces it would be only necessary to
!      calculate force for representative and rotate the force itself
               IF (oneD%odi%d1) THEN
                  inap = oneD%ods%ngopr(natom)
!                 nap = ods%ngopr(natom)
!                 inap = ods%invtab(nap)
               ELSE
                  nap = atoms%ngopr(natom)
                  inap = sym%invtab(nap)
               END IF
               DO j = 1,3
                  fkr(j) = 0.
                  DO i = 1,3
                     IF (.NOT.oneD%odi%d1) THEN
                        fkr(j) = fkr(j) + fk(i)*sym%mrot(i,j,inap)
                     ELSE
                        fkr(j) = fkr(j) + fk(i)*oneD%ods%mrot(i,j,inap)
                     END IF             
                  END DO
               END DO
!
!       TRANSFORM (k+g) TO CART. SYSTEM
               fkp=matmul(fkr,cell%bmat)
  !     ----> generate spherical harmonics
               CALL ylm4(atoms%lmax(n),fkp, ylm)
!
!  needed vor aveccof,bveccof (equ. A18)
!  ROTATION:  kr = kaux*R_equiv
               DO j = 1,3
                  fgr(j) = 0.
                  DO i = 1,3
                     IF (.NOT.oneD%odi%d1) THEN
                        fgr(j) = fgr(j) + fg(i)*sym%mrot(i,j,inap)
                     ELSE
                        fgr(j) = fgr(j) + fg(i)*oneD%ods%mrot(i,j,inap)
                     END IF
                  END DO
               END DO
!       TRANSFORM RECIPROCAL LATTICE VECTOR g TO CART. SYSTEM
fgp=matmul(fgr,cell%bmat)
!
!     ----> loop over l,m
               DO l = 0,atoms%lmax(n)
                  ll1 = l* (l+1)
                  DO m = -l,l
                     lm = ll1 + m
                     c_0 = conjg(ylm(lm+1)) * phase
                     c_1 = c_0 *  fj(l)
                     c_2 = c_0 * dfj(l)
!     ----> loop over bands
                     DO ie = 1,ne
                       aaux(ie) = c_1 * work(ie)
                       baux(ie) = c_2 * work(ie)
                     ENDDO
                     DO ie = 1,ne
                       acof(ie,lm,natom) = acof(ie,lm,natom) + aaux(ie)
                       bcof(ie,lm,natom) = bcof(ie,lm,natom) + baux(ie)
                     ENDDO

                     IF ( atoms%l_geo(n) ) THEN
                       DO ie = 1,ne
!
                          s2h_e = s2h-eig(ie)
                          e1cof(ie,lm,natom) = e1cof(ie,lm,natom) + aaux(ie)* s2h_e
                          e2cof(ie,lm,natom) = e2cof(ie,lm,natom) + baux(ie)* s2h_e
                          DO i = 1,3
                             aveccof(i,ie,lm,natom) = aveccof(i,ie,lm,natom) + aaux(ie)*fgp(i)
                          END DO
                          DO i = 1,3
                             bveccof(i,ie,lm,natom) = bveccof(i,ie,lm,natom) + baux(ie)*fgp(i)
                          END DO
!
                       END DO
                     ENDIF
#if ( defined(CPP_SOC) && defined(CPP_INVERSION) )
                     IF (atoms%invsat(natom).EQ.1) THEN
                        jatom = sym%invsatnr(natom)
                        lmp = ll1 - m
                        inv_f = (-1)**(l-m)
                        c_1 =  conjg(c_1) * inv_f
                        c_2 =  conjg(c_2) * inv_f
                        DO ie = 1,ne
                          aaux(ie) = c_1 * work(ie)
                          baux(ie) = c_2 * work(ie)
                        ENDDO
                        DO ie = 1,ne
                          acof(ie,lmp,jatom) = acof(ie,lmp,jatom) + aaux(ie)
                          bcof(ie,lmp,jatom) = bcof(ie,lmp,jatom) + baux(ie)
                        ENDDO
                        IF ( atoms%l_geo(n) ) THEN
                          DO ie = 1,ne
!
                             s2h_e = s2h-eig(ie)
                             e1cof(ie,lmp,jatom) = e1cof(ie,lmp,jatom) + aaux(ie)* s2h_e
                             e2cof(ie,lmp,jatom) = e2cof(ie,lmp,jatom) + baux(ie)* s2h_e
                             DO i = 1,3
                               aveccof(i,ie,lmp,jatom) = aveccof(i,ie,lmp,jatom) - aaux(ie)*fgp(i)
                             END DO
                             DO i = 1,3
                               bveccof(i,ie,lmp,jatom) = bveccof(i,ie,lmp,jatom) - baux(ie)*fgp(i)
                             END DO
!
                          END DO
                        ENDIF
                      ENDIF   ! atoms%invsat(na) = 1
#endif
!
! END m loop
                  END DO
!
! END l loop
               END DO
               IF (.NOT.enough(natom)) THEN
                  CALL abclocdn_pulay(atoms, sym, noco,ccchi(1,jspin),kspin, iintsp,const,phase,ylm,n,natom,&
                      k,fgp,s,nvmax,ne,z,nbasf0, alo1,blo1,clo1,kvec(1,1,natom), nkvec,&
                      enough,acof,bcof,ccof, acoflo,bcoflo,aveccof,bveccof,cveccof)
               END IF
!-inv
            END IF
   enddo
   enddo
!
!     ---> end loop plane waves
   enddo
!---> end loop over interstitial spin
      ENDDO
!
! -> calculate inversion-symmetric atoms
!
#if ( defined(CPP_SOC) && defined(CPP_INVERSION) )
!
! See the correpsonding remarks in abcof() at this place. Note, in addition that
!
!      -p,n          (l+m)      p,n  *
!  Avec      = - (-1)      (Avec    )  i.e. an additional minus sign enters.
!       l,m                    l,-m
!
#else
      iatom = 0
      DO n = 1,atoms%ntype
         DO nn = 1,atoms%neq(n)
            iatom = iatom + 1
            IF (atoms%invsat(iatom).EQ.1) THEN
               jatom = sym%invsatnr(iatom)
               DO ilo = 1,atoms%nlo(n)
                  l = atoms%llo(ilo,n)
                  DO m = -l,l
                     inv_f = (-1.0)**(m+l)
                     DO ie = 1,ne
                        acoflo(m,ie,ilo,jatom) = inv_f * conjg(acoflo(-m,ie,ilo,iatom))
                        bcoflo(m,ie,ilo,jatom) = inv_f * conjg(bcoflo(-m,ie,ilo,iatom))
                          ccof(m,ie,ilo,jatom) = inv_f * conjg(  ccof(-m,ie,ilo,iatom))
                        DO i = 1,3
                           cveccof(i,m,ie,ilo,jatom) = -inv_f * conjg(cveccof(i,-m,ie,ilo,iatom))
                        ENDDO
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
                        acof(ie,lm,jatom) = inv_f * conjg(acof(ie,lmp,iatom))
                        bcof(ie,lm,jatom) = inv_f * conjg(bcof(ie,lmp,iatom))
                     ENDDO
!+fo
                     IF ( atoms%l_geo(n) ) THEN
                       DO ie = 1,ne
                          e1cof(ie,lm,jatom) = inv_f * conjg(e1cof(ie,lmp,iatom))
                          e2cof(ie,lm,jatom) = inv_f * conjg(e2cof(ie,lmp,iatom))
                          DO i = 1,3
                             aveccof(i,ie,lm,jatom) = -inv_f * conjg(aveccof(i,ie,lmp,iatom))
                             bveccof(i,ie,lm,jatom) = -inv_f * conjg(bveccof(i,ie,lmp,iatom))
                          END DO
                       ENDDO
                     ENDIF
!-fo
                  ENDDO
               ENDDO
            ENDIF
         ENDDO
      ENDDO
#endif
      DEALLOCATE ( aaux,baux,work )

      END SUBROUTINE to_pulay
      END MODULE m_topulay
