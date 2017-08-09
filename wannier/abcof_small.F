      MODULE m_abcof_small
      CONTAINS
      SUBROUTINE abcof_small(
     >                 sym,atoms,usdus,lapw,noco,zMat,
     >                 lmaxd,ntypd,neigd,nobd,natd,nop,nvd,jspd,
     >                 lmd,nbasfcn,llod,nlod,nlotot,invtab,
     >                 ntype,mrot,ngopr,taual,neq,lmax,rmt,omtil,
     >                 bmat,bbmat,bkpt,k1,k2,k3,nv,nmat,ne,z,
     >                 us,dus,uds,duds,ddn,invsat,invsatnr,
     >                 ulos,uulon,dulon,dulos,llo,nlo,l_dulo,lapw_l,
     >                 l_noco,l_ss,jspin,alph,beta,qss,kveclo,odi,ods,
     <                 acof,bcof,ccof,nig)
c     ************************************************************
c     subroutine constructs the a,b coefficients of the linearized
c     m.t. wavefunctions for each band and atom.       c.l. fu
c     ************************************************************
#include "cpp_double.h"

      USE m_constants
      USE m_setabc1locdn
      USE m_sphbes
      USE m_dsphbs
      USE m_abclocdn
      USE m_ylm
      USE m_types

      IMPLICIT NONE

      TYPE(t_sym),INTENT(IN)     :: sym
      TYPE(t_atoms),INTENT(IN)   :: atoms
      TYPE(t_usdus),INTENT(IN)   :: usdus
      TYPE(t_lapw),INTENT(IN)    :: lapw
      TYPE(t_noco),INTENT(IN)    :: noco
      TYPE(t_zMat),INTENT(IN)    :: zMat

C     ..
C     .. Scalar Arguments ..
      INTEGER, INTENT (IN) :: lmaxd,ntypd,neigd,nobd,natd,nop,nvd,jspd
      INTEGER, INTENT (IN) :: ne,ntype,nmat,nbasfcn,llod,nlod,lmd
      REAL,    INTENT (IN) :: omtil
      INTEGER, INTENT (IN) :: jspin,nlotot,nig
      LOGICAL, INTENT (IN) :: l_noco,l_ss
C     ..
C     .. Array Arguments ..
      INTEGER, INTENT (IN) :: mrot(3,3,nop),ngopr(natd),lmax(ntypd)
      INTEGER, INTENT (IN) :: k1(nvd,jspd),k2(nvd,jspd),k3(nvd,jspd)
      INTEGER, INTENT (IN) :: nv(jspd),lapw_l(ntypd),invtab(nop)
      INTEGER, INTENT (IN) :: neq(ntypd),invsat(natd),invsatnr(natd)
      INTEGER, INTENT (IN) :: nlo(ntypd),llo(nlod,ntypd),kveclo(nlotot)
      REAL,    INTENT (IN) :: bmat(3,3),bbmat(3,3),bkpt(3)
      REAL,    INTENT (IN) :: taual(3,natd),rmt(ntypd)
      REAL,    INTENT (IN) :: dus(0:lmaxd,ntypd),duds(0:lmaxd,ntypd)
      REAL,    INTENT (IN) ::  us(0:lmaxd,ntypd), uds(0:lmaxd,ntypd)
      REAL,    INTENT (IN) :: ulos(nlod,ntypd),uulon(nlod,ntypd)
      REAL,    INTENT (IN) :: dulon(nlod,ntypd),dulos(nlod,ntypd)
      REAL,    INTENT (IN) :: ddn(0:lmaxd,ntypd)
      REAL,    INTENT (IN) :: alph(ntypd),beta(ntypd),qss(3)
#if ( !defined(CPP_INVERSION) || defined(CPP_SOC) )
      COMPLEX, INTENT (IN) :: z(nbasfcn,neigd)
#else
      REAL,    INTENT (IN) :: z(nbasfcn,neigd)
#endif
      LOGICAL, INTENT (IN) :: l_dulo(nlod,ntypd)
      COMPLEX, INTENT (OUT):: acof(nobd,0:lmd,natd)
      COMPLEX, INTENT (OUT):: bcof(nobd,0:lmd,natd)
      COMPLEX, INTENT (OUT):: ccof(-llod:llod,nobd,nlod,natd)
c-odim
      TYPE (od_inp), INTENT (IN) :: odi
      TYPE (od_sym), INTENT (IN) :: ods
c+odim
C     ..
C     .. Local Scalars ..
      COMPLEX phase,cexp,c_0,c_1,c_2,ci
      REAL const,df,r1,s,tmk,wronk,tpi,qss1,qss2,qss3
      INTEGER i,j,k,l,ll1,lm,m,n,nap,natom,nn,iatom,jatom,lmp
      INTEGER inv_f,ie,ilo,kspin,iintsp,nintsp,nvmax,lo,inap
C     ..
C     .. Local Arrays ..
      INTEGER kvec(2*(2*llod+1),nlod,natd)
      INTEGER nbasf0(nlod,natd),nkvec(nlod,natd)
      REAL dfj(0:lmaxd),fj(0:lmaxd),fk(3),fkp(3),fkr(3)
      REAL alo1(nlod,ntypd),blo1(nlod,ntypd),clo1(nlod,ntypd)
      COMPLEX ylm( (lmaxd+1)**2 )
      COMPLEX ccchi(2,2)
      LOGICAL enough(natd),apw(0:lmaxd,ntypd)
#if ( !defined(CPP_INVERSION) || defined(CPP_SOC) )
      COMPLEX, ALLOCATABLE :: work(:)
#else
      REAL,    ALLOCATABLE :: work(:)
#endif
C     ..
      ci = cmplx(0.0,1.0)
      tpi = 2 * pimach()
      const = 2 * tpi/sqrt(omtil)
c
      acof(:,:,:) = cmplx(0.0,0.0)
      bcof(:,:,:) = cmplx(0.0,0.0)
      ALLOCATE ( work(nobd) )
C     ..
c+APW_LO
      DO n = 1, ntype
         DO l = 0,lmax(n)
           apw(l,n) = .false.
           DO lo = 1,nlo(n)
             IF (l_dulo(lo,n)) apw(l,n) = .true.
           ENDDO
#ifdef CPP_APW
           IF (lapw_l(n).GE.l) apw(l,n) = .false.
#endif
         ENDDO
         DO lo = 1,nlo(n)
           IF (l_dulo(lo,n)) apw(llo(lo,n),n) = .true.
         ENDDO
      ENDDO
c+APW_LO
c
      nintsp = 1
      IF (l_ss) nintsp = 2
c---> loop over the interstitial spin
      DO iintsp = 1,nintsp
      nvmax = nv(jspin)
      IF (l_ss) nvmax = nv(iintsp)
c
      CALL setabc1locdn(
     >                  jspin,atoms,lapw,ne,noco,iintsp,sym,usdus,
     >                  kveclo,
     <                  enough,nkvec,kvec,nbasf0,ccof,
     <                  alo1,blo1,clo1)
c
      IF (iintsp .EQ. 1) THEN
         qss1= - qss(1)/2
         qss2= - qss(2)/2
         qss3= - qss(3)/2
      ELSE
         qss1= + qss(1)/2
         qss2= + qss(2)/2
         qss3= + qss(3)/2
      ENDIF

      k = nig

c---> loop over atom types
      natom = 0
      DO n = 1,ntype
c  ----> loop over equivalent atoms
         DO nn = 1,neq(n)
            natom = natom + 1
           IF ((invsat(natom).EQ.0) .OR. (invsat(natom).EQ.1)) THEN
c--->        loop over lapws
             !DO k = 1,nvmax
               IF (.NOT.l_noco) THEN
                 DO i = 1,ne
                   work(i) = z(k,i)
                 ENDDO
               ENDIF

               IF (l_noco) THEN
c--->            generate the complex conjgates of the spinors (chi)
                 ccchi(1,1) = conjg( exp(-ci*alph(n)/2)*cos(beta(n)/2))
                 ccchi(1,2) = conjg(-exp(-ci*alph(n)/2)*sin(beta(n)/2))
                 ccchi(2,1) = conjg( exp( ci*alph(n)/2)*sin(beta(n)/2))
                 ccchi(2,2) = conjg( exp( ci*alph(n)/2)*cos(beta(n)/2))
                 IF (l_ss) THEN
c--->              the coefficients of the spin-down basis functions are
c--->              stored in the second half of the eigenvector
                   kspin = (iintsp-1)*(nv(1)+nlotot)
                   DO i = 1,ne
                      work(i) = ccchi(iintsp,jspin)*z(kspin+k,i)
                   ENDDO
                 ELSE
c--->              perform sum over the two interstitial spin directions
c--->              and take into account the spin boundary conditions
c--->              (jspin counts the local spin directions inside each MT)
                   kspin = nv(1)+nlotot
                   DO i = 1,ne
                      work(i) = ccchi(1,jspin)*z(k,i)
     +                        + ccchi(2,jspin)*z(kspin+k,i)
                   ENDDO
                 ENDIF
               ENDIF ! (l_noco)
               IF (l_ss) THEN
                 fk(1) = bkpt(1) + k1(k,iintsp) + qss1
                 fk(2) = bkpt(2) + k2(k,iintsp) + qss2
                 fk(3) = bkpt(3) + k3(k,iintsp) + qss3
               ELSE
                 fk(1) = bkpt(1) + k1(k,jspin) + qss1
                 fk(2) = bkpt(2) + k2(k,jspin) + qss2
                 fk(3) = bkpt(3) + k3(k,jspin) + qss3
               ENDIF ! (l_ss)
               s = dot_product(fk,matmul(bbmat,fk))
!               s = dotirp(fk,fk,bbmat)
               s = sqrt(s)
               r1 = rmt(n)*s
               CALL sphbes(
     >                     lmax(n),r1,
     <                     fj)
               CALL dsphbs(
     >                     lmax(n),r1,fj,
     <                     dfj)
c  ----> construct a and b coefficients
               DO l = 0,lmax(n)
                 df = s*dfj(l)
                 wronk = uds(l,n)*dus(l,n) - us(l,n)*duds(l,n)
                 IF (apw(l,n)) THEN
                   fj(l) = 1.0*const * fj(l)/us(l,n)
                   dfj(l) = 0.0d0
                 ELSE
                   dfj(l) = const* (dus(l,n)*fj(l)-df*us(l,n))/wronk
                   fj(l) = const* (df*uds(l,n)-fj(l)*duds(l,n))/wronk
                 ENDIF
               ENDDO ! loop over l
               tmk = tpi* (fk(1)*taual(1,natom)+
     +                     fk(2)*taual(2,natom)+
     +                     fk(3)*taual(3,natom))
               phase = cmplx(cos(tmk),sin(tmk))
               IF (odi%d1) THEN
                  inap = ods%ngopr(natom)
               ELSE
                  nap = ngopr(natom)
                  inap = invtab(nap)
               END IF
               DO j = 1,3
                 fkr(j) = 0.
                 DO i = 1,3
                   IF (odi%d1) THEN
                     fkr(j) = fkr(j) + fk(i)*ods%mrot(i,j,inap)
                   ELSE
                     fkr(j) = fkr(j) + fk(i)*mrot(i,j,inap)
                   END IF
                 ENDDO
               ENDDO
               fkp=MATMUL(fkr,bmat)
!               CALL cotra3(fkr,fkp,bmat)
c     ----> generate spherical harmonics
               CALL ylm4(
     >                   lmax(n),fkp,
     <                   ylm)
c     ----> loop over l
               DO l = 0,lmax(n)
                 ll1 = l* (l+1)
c     ----> loop over m
                 DO m = -l,l
                   lm = ll1 + m
                   c_0 = conjg(ylm(lm+1))*phase
                   c_1 = c_0 *  fj(l)
                   c_2 = c_0 * dfj(l)
c     ----> loop over bands
                   DO i = 1,ne
                     acof(i,lm,natom) = acof(i,lm,natom) + 
     +                                  c_1 * work(i)
                   ENDDO
                   DO i = 1,ne
                     bcof(i,lm,natom) = bcof(i,lm,natom) +
     +                                  c_2 * work(i)
                   ENDDO
#if ( defined(CPP_SOC) && defined(CPP_INVERSION) )
                   IF (invsat(natom).EQ.1) THEN
                     jatom = invsatnr(natom)
                     lmp = ll1 - m
                     inv_f = (-1)**(l-m)
                     c_1 =  conjg(c_1) * inv_f
                     c_2 =  conjg(c_2) * inv_f
                     CALL CPP_BLAS_caxpy(ne,c_1,work,1,
     X                                   acof(1,lmp,jatom),1)
                     CALL CPP_BLAS_caxpy(ne,c_2,work,1,
     X                                   bcof(1,lmp,jatom),1)
                   ENDIF
#endif
                 ENDDO ! loop over m
               ENDDO ! loop over l
               IF (.NOT.enough(natom)) THEN
                 write(*,*)'.not.enough(natom)'
                 CALL abclocdn(
     >                 atoms,sym,noco,ccchi(1,jspin),kspin,iintsp,
     >                 const,phase,ylm,n,natom,k,s,nvmax,
     >                 ne,nbasf0,alo1,blo1,clo1,kvec(1,1,natom),
     <                 nkvec,enough(natom),acof,bcof,ccof,zMat)
               ENDIF
             !ENDDO ! loop over LAPWs
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
      DO n = 1,ntype
         DO nn = 1,neq(n)
            iatom = iatom + 1
            IF (invsat(iatom).EQ.1) THEN
              jatom = invsatnr(iatom)
              cexp = exp(tpi*ci*dot_product(taual(:,jatom)
     +             + taual(:,iatom),(/bkpt(1),bkpt(2),bkpt(3)/)))
               DO ilo = 1,nlo(n)
                  l = llo(ilo,n)
                  DO m = -l,l
                     inv_f = (-1.0)**(m+l)
                     DO ie = 1,ne
                        ccof(m,ie,ilo,jatom) = inv_f * cexp *
     +                               conjg(  ccof(-m,ie,ilo,iatom))
                     ENDDO
                  ENDDO
               ENDDO
               DO l = 0,lmax(n)
                  ll1 = l* (l+1)
                  DO m =-l,l
                     lm  = ll1 + m
                     lmp = ll1 - m
                     inv_f = (-1.0)**(m+l)
                     DO ie = 1,ne
                        acof(ie,lm,jatom) = inv_f * cexp * 
     *                                      conjg(acof(ie,lmp,iatom))
                     ENDDO
                     DO ie = 1,ne
                        bcof(ie,lm,jatom) = inv_f * cexp * 
     *                                      conjg(bcof(ie,lmp,iatom))
                     ENDDO
                  ENDDO
               ENDDO
            ENDIF
         ENDDO
      ENDDO
#endif
      DEALLOCATE ( work )
c
      END SUBROUTINE abcof_small
      END MODULE m_abcof_small
