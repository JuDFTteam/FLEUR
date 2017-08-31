!--------------------------------------------------------------------------------
! Copyright (c) 2016 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------

      MODULE m_wann_real
c     ********************************************************
c     calculates the value of the periodic part of the
c     wavefunction at the given real-grid point p(:)
c                          Y.Mokrousov 16.8.6
c     ********************************************************
      CONTAINS
      SUBROUTINE wann_real(
     >                  p,n,na,iv,iflag,bkpt,iband,
     >                  n3d,nmzxyd,n2d,ntypsd,lmaxd,jmtd,
     >                  natd,ntypd,nmzd,nop,nop2,mrot,tau,invtab,
     >                  nq3,nvac,invs,z1,delz,nmz,nmzxy,nq2,
     >                  lmax,rmsh,jri,pos,ngopr,ntypsy,nvd,
     >                  omtil,amat,bmat,odi,ods,nlod,llod,nlo,llo,
     >                  ff,gg,flo,acof,bcof,ccof,zMat,
     >                  nv,k1,k2,k3,lmd,nbasfcn,l_ss,qss,jspin,addnoco,
     <                  xdnout)

      USE m_types
      USE m_ylm
      USE m_constants
      IMPLICIT NONE

      TYPE(t_zMat),INTENT(IN)        :: zMat

C     .. Scalar Arguments ..
      INTEGER, INTENT (IN) :: n3d,nmzxyd,n2d,ntypsd,llod,nlod,iband
      INTEGER, INTENT (IN) :: lmaxd,jmtd,ntypd,natd,nmzd
      INTEGER, INTENT (IN) :: iflag,n,na,iv,lmd,nv,nvd,nbasfcn
      INTEGER, INTENT (IN) :: nq3,nvac,nmz,nmzxy,nq2,nop,nop2
      LOGICAL, INTENT (IN) :: invs,l_ss
      REAL,    INTENT (IN) :: z1,delz,omtil,bkpt(3),qss(3)
      INTEGER, INTENT (IN) :: jspin,addnoco
      COMPLEX, INTENT (OUT):: xdnout
c-odim
      TYPE (od_inp), INTENT (IN) :: odi
      TYPE (od_sym), INTENT (IN) :: ods
c+odim
C     ..
C     .. Array Arguments ..
      INTEGER, INTENT (IN) :: ngopr(natd),ntypsy(natd),jri(ntypd)
      INTEGER, INTENT (IN) :: lmax(ntypd),mrot(3,3,nop),invtab(nop)
      INTEGER, INTENT (IN) :: nlo(ntypd),llo(nlod,ntypd)
      REAL,    INTENT (IN) :: amat(3,3),bmat(3,3),pos(3,natd)
      REAL,    INTENT (IN) :: rmsh(jmtd,ntypd),tau(3,nop)
      INTEGER, INTENT (IN) :: k1(nvd),k2(nvd),k3(nvd) 
      COMPLEX, INTENT (IN) :: ccof(-llod:llod,nlod,natd)
      COMPLEX, INTENT (IN) :: acof(0:lmd,natd)
      COMPLEX, INTENT (IN) :: bcof(0:lmd,natd)
      REAL,    INTENT (IN) :: ff(ntypd,jmtd,2,0:lmaxd)
      REAL,    INTENT (IN) :: gg(ntypd,jmtd,2,0:lmaxd)
      REAL,    INTENT (IN) :: flo(ntypd,jmtd,2,nlod)
      REAL,    INTENT (INOUT) :: p(3)

C     ..
C     .. Local Scalars ..
      REAL delta,sx,xx1,xx2,rrr,phi,const,arg,tpi,arg1
      INTEGER i,j,jp3,jr,k,l,nd,nopa,ivac,lm,m,gzi,kk
      INTEGER kk1,kk2,kk3
      COMPLEX ci,const2,s,xd1,xd2,const3
C     ..
C     .. Local Arrays ..
      COMPLEX sf2(n2d),sf3(n3d),ylm((lmaxd+1)**2)
      REAL rcc(3),x(3),rcc2(3)
      REAL bqpt(3)

C     ..
      ci = cmplx(0.,1.)
      tpi   = 2 * pimach()
      const = 1./(sqrt(omtil))

c..define the factor e^{-ikr}
      rcc2=matmul(bmat,p)/tpi_const

      bqpt = 0.0
!      if(l_ss.and.jspin.eq.1) then 
!         bqpt = -qss/2.0
!      elseif(l_ss.and.jspin.eq.2) then
!         bqpt = +qss/2.0
!      endif

      arg = -tpi*(   (bkpt(1)+bqpt(1))*rcc2(1) 
     >             + (bkpt(2)+bqpt(2))*rcc2(2)
     >             + (bkpt(3)+bqpt(3))*rcc2(3)  )

      arg1 = tpi*( bkpt(1)*rcc2(1) + bkpt(2)*rcc2(2) + bkpt(3)*rcc2(3) )
      const2 = cmplx(cos(arg),sin(arg))
      const3 = cmplx(cos(arg1),sin(arg1))
c     write (6,*) 'bkpt,const2,const3=',bkpt(:),const2,const3

      ivac=iv

      IF (iflag.EQ.0) GO TO 20
      IF (iflag.EQ.1) GO TO 40
c     ---> interstitial part
      rcc=matmul(bmat,p)/tpi_const
      xdnout = cmplx(0.,0.)
c     write (6,*) 'nv,nvd=',nv,nvd
      IF (zMat%l_real) THEN
         DO k = 1,nv
c           write (6,*) 'k1,k2,k3=',k1(k),k2(k),k3(k)
c           write (6,*) 'z(k,iband)=', z(k,iband)
            arg = tpi * ((k1(k))*rcc(1)+(k2(k))*rcc(2)+(k3(k))*rcc(3))
            xdnout = xdnout + zMat%z_r(k+addnoco,iband)*
     +                        cmplx(cos(arg),sin(arg))*const
            IF (((abs(p(1)-2.2).le.0.0001).and.(abs(p(2)).le.0.0001))
     &    .or.((abs(p(2)-2.2).le.0.0001).and.(abs(p(1)).le.0.0001)))then
c              write (6,*) 'p(i)=',p(1:2)
c              write (6,*) 'G=',k1(k),k2(k),k3(k)
c              write (6,*) 'z(k,iband)=',z(k,iband)
c              write (6,*) 'val=',z(k,iband)*cmplx(cos(arg),sin(arg))
            ENDIF
         END DO
      ELSE
         DO k = 1,nv
c           write (6,*) 'k1,k2,k3=',k1(k),k2(k),k3(k)
c           write (6,*) 'z(k,iband)=', z(k,iband)
            arg = tpi * ((k1(k))*rcc(1)+(k2(k))*rcc(2)+(k3(k))*rcc(3))
            xdnout = xdnout + zMat%z_c(k+addnoco,iband)*
     +                        cmplx(cos(arg),sin(arg))*const
            IF (((abs(p(1)-2.2).le.0.0001).and.(abs(p(2)).le.0.0001))
     &    .or.((abs(p(2)-2.2).le.0.0001).and.(abs(p(1)).le.0.0001)))then
c              write (6,*) 'p(i)=',p(1:2)
c              write (6,*) 'G=',k1(k),k2(k),k3(k)
c              write (6,*) 'z(k,iband)=',z(k,iband)
c              write (6,*) 'val=',z(k,iband)*cmplx(cos(arg),sin(arg))
            ENDIF
         END DO
      END IF
c     write (6,*) 'ir:p(i)',p(:)
      RETURN
c     ---> vacuum part
   20 CONTINUE
      xdnout = cmplx(0.,0.)
      return
c     ----> m.t. part
   40 CONTINUE
      

      nd = ntypsy(na)
      nopa = ngopr(na)
      nopa=1
     
      IF (odi%d1) nopa = ods%ngopr(na)
      sx = 0.0
      DO 50 i = 1,3
         x(i) = p(i) - pos(i,na)
         sx = sx + x(i)*x(i)
   50 CONTINUE
      sx = sqrt(sx)
      IF (nopa.NE.1) THEN
c... switch to internal units
         rcc=matmul(bmat,p)/tpi_const
c... rotate into representative
         DO 70 i = 1,3
            p(i) = 0.
            DO 60 j = 1,3
              IF (.NOT.odi%d1) THEN
               p(i) = p(i) + mrot(i,j,nopa)*rcc(j)
              ELSE
               p(i) = p(i) + ods%mrot(i,j,nopa)*rcc(j)
              END IF
   60       CONTINUE
   70    CONTINUE
c... switch back to cartesian units
         x=matmul(amat,p)/tpi_const
      END IF
      DO 80 j = jri(n),2,-1
         IF (sx.GE.rmsh(j,n)) GO TO 90
   80 CONTINUE
   90 jr = j
      CALL ylm4(
     >          lmax(n),x,
     <          ylm)
      xd1 = cmplx(0.,0.)
      xd2 = cmplx(0.,0.)
      DO l = 0,lmax(n)
c        if (p(1).eq.0. .and. p(2).eq.0. .and. p(3).eq.0)then
c               write (6,*) 'ff(l,300)=',ff(1,300,1,l)
c               write (6,*) 'ff(l,300)=',ff(1,300,2,l)
c               write (6,*) 'gg(l,300)=',gg(1,300,1,l)
c               write (6,*) 'gg(l,300)=',gg(1,300,2,l)
c        endif
       DO 110 m = -l,l
        lm = l*(l+1)+m
        s = ylm(lm+1)*(ci)**l
c       if (p(1).eq.0. .and. p(2).eq.0. .and. p(3).eq.0)then
c              write (6,*) 'acof=',acof(lm,1)
c              write (6,*) 'bcof=',bcof(lm,1)
c       endif
        xd1 = xd1 + (acof(lm,na)*cmplx(ff(n,jr,1,l),0.)+
     +               bcof(lm,na)*cmplx(gg(n,jr,1,l),0.))*s/
     /               (rmsh(jr,n)) 
c    /               (rmsh(jr,n)*rmsh(jr,n))
        IF (jr.EQ.1) GO TO 110
        xd2 = xd2 + (acof(lm,na)*cmplx(ff(n,jr+1,1,l),0.)+
     +               bcof(lm,na)*cmplx(gg(n,jr+1,1,l),0.))*s/  
     /               (rmsh(jr+1,n))
c    /               (rmsh(jr+1,n)*rmsh(jr+1,n))
  110  CONTINUE
      ENDDO
c..contributions from the local orbitals
      IF (nlo(n).GE.1) THEN
       DO l = 1,nlo(n)
        DO 111 m = -llo(l,n),llo(l,n)
         lm = llo(l,n)*(llo(l,n)+1)+m
         s = ylm(lm+1)*(ci)**llo(l,n) 
         xd1 = xd1 + ccof(m,l,na)*flo(n,jr,1,l)*s/
     /               (rmsh(jr,n))         
         IF (jr.EQ.1) GO TO 111
         xd2 = xd2 + ccof(m,l,na)*flo(n,jr+1,1,l)*s/
     /               (rmsh(jr+1,n))         
  111   CONTINUE
       ENDDO
      ENDIF    
      IF (jr.EQ.1) THEN
         xdnout = xd1
      ELSE
         xdnout = xd1 + (xd2-xd1) *
     +                  (sx-rmsh(jr,n)) / (rmsh(jr+1,n)-rmsh(jr,n))
         
      END IF
      xdnout = xdnout*const2
c     write (6,*) 'mt:p(i)',p(:)
 8000 FORMAT (2f10.6)
c
      RETURN
      END SUBROUTINE wann_real
      END MODULE m_wann_real
