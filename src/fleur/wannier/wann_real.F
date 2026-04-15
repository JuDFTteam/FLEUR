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
     >                  atoms,cell,sym,
     >                  ff,gg,flo,acof,bcof,ccof,zMat,
     >                  nv,k1,k2,k3,lmd,nbasfcn,
     >                  l_ss,qss,jspin,addnoco,
     <                  xdnout)

      USE m_types
      USE m_ylm
      USE m_constants
      use m_juDFT

      IMPLICIT NONE

      TYPE(t_atoms), INTENT(IN) :: atoms
      TYPE(t_cell),  INTENT(IN) :: cell
      TYPE(t_sym),   INTENT(IN) :: sym
      TYPE(t_mat),INTENT(IN)    :: zMat

C     .. Scalar Arguments ..
      INTEGER, INTENT (IN) :: iband
      INTEGER, INTENT (IN) :: iflag,n,na,iv,lmd,nv,nbasfcn
      LOGICAL, INTENT (IN) :: l_ss
      REAL,    INTENT (IN) :: bkpt(3),qss(3)
      INTEGER, INTENT (IN) :: jspin,addnoco
      COMPLEX, INTENT (OUT):: xdnout

C     ..
C     .. Array Arguments ..
      INTEGER, INTENT (IN) :: k1(:),k2(:),k3(:)
      COMPLEX, INTENT (IN) :: ccof(-atoms%llod:,:,:)
      COMPLEX, INTENT (IN) :: acof(0:,:)
      COMPLEX, INTENT (IN) :: bcof(0:,:)
      REAL,    INTENT (IN) :: ff(:,:,:,0:)
      REAL,    INTENT (IN) :: gg(:,:,:,0:)
      REAL,    INTENT (IN) :: flo(:,:,:,:)
      REAL,    INTENT (INOUT) :: p(3)

C     ..
C     .. Local Scalars ..
      REAL delta,sx,xx1,xx2,rrr,phi,const,arg,tpi,arg1
      INTEGER i,j,jp3,jr,k,l,nd,nopa,ivac,lm,m,gzi,kk
      INTEGER kk1,kk2,kk3
      COMPLEX const2,s,xd1,xd2,const3
C     ..
C     .. Local Arrays ..
      COMPLEX ylm((atoms%lmaxd+1)**2)
      REAL rcc(3),x(3),rcc2(3)
      REAL bqpt(3)

      call timestart("wann_real")
      tpi   = 2 * pimach()
      const = 1./(sqrt(cell%omtil))

c..define the factor e^{-ikr}
      rcc2=matmul(cell%bmat,p)/tpi_const

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
c     write (oUnit,*) 'bkpt,const2,const3=',bkpt(:),const2,const3

      ivac=iv

      IF (iflag.EQ.0) GO TO 20
      IF (iflag.EQ.1) GO TO 40
c     ---> interstitial part
      rcc=matmul(cell%bmat,p)/tpi_const
      xdnout = cmplx(0.,0.)
c     write (oUnit,*) 'nv,nvd=',nv,nvd
      IF (zMat%l_real) THEN
         DO k = 1,nv
c           write (oUnit,*) 'k1,k2,k3=',k1(k),k2(k),k3(k)
c           write (oUnit,*) 'z(k,iband)=', z(k,iband)
            arg = tpi * ((k1(k))*rcc(1)+(k2(k))*rcc(2)+(k3(k))*rcc(3))
            xdnout = xdnout + zMat%data_r(k+addnoco,iband)*
     +                        cmplx(cos(arg),sin(arg))*const
            IF (((abs(p(1)-2.2).le.0.0001).and.(abs(p(2)).le.0.0001))
     &    .or.((abs(p(2)-2.2).le.0.0001).and.(abs(p(1)).le.0.0001)))then
c              write (oUnit,*) 'p(i)=',p(1:2)
c              write (oUnit,*) 'G=',k1(k),k2(k),k3(k)
c              write (oUnit,*) 'z(k,iband)=',z(k,iband)
c              write (oUnit,*) 'val=',z(k,iband)*cmplx(cos(arg),sin(arg))
            ENDIF
         END DO
      ELSE
         DO k = 1,nv
c           write (oUnit,*) 'k1,k2,k3=',k1(k),k2(k),k3(k)
c           write (oUnit,*) 'z(k,iband)=', z(k,iband)
            arg = tpi * ((k1(k))*rcc(1)+(k2(k))*rcc(2)+(k3(k))*rcc(3))
            xdnout = xdnout + zMat%data_c(k+addnoco,iband)*
     +                        cmplx(cos(arg),sin(arg))*const
            IF (((abs(p(1)-2.2).le.0.0001).and.(abs(p(2)).le.0.0001))
     &    .or.((abs(p(2)-2.2).le.0.0001).and.(abs(p(1)).le.0.0001)))then
c              write (oUnit,*) 'p(i)=',p(1:2)
c              write (oUnit,*) 'G=',k1(k),k2(k),k3(k)
c              write (oUnit,*) 'z(k,iband)=',z(k,iband)
c              write (oUnit,*) 'val=',z(k,iband)*cmplx(cos(arg),sin(arg))
            ENDIF
         END DO
      END IF
c     write (oUnit,*) 'ir:p(i)',p(:)
      call timestop("wann_real")
      RETURN
c     ---> vacuum part
   20 CONTINUE
      xdnout = cmplx(0.,0.)
      call timestop("wann_real")
      return
c     ----> m.t. part
   40 CONTINUE
      

      nd = sym%ntypsy(na)
      nopa = sym%ngopr(na)
      nopa=1
     
   
      sx = 0.0
      DO 50 i = 1,3
         x(i) = p(i) - atoms%pos(i,na)
         sx = sx + x(i)*x(i)
   50 CONTINUE
      sx = sqrt(sx)
      IF (nopa.NE.1) THEN
c... switch to internal units
         rcc=matmul(cell%bmat,p)/tpi_const
c... rotate into representative
         DO 70 i = 1,3
            p(i) = 0.
            DO 60 j = 1,3
               p(i) = p(i) + sym%mrot(i,j,nopa)*rcc(j)
            
   60       CONTINUE
   70    CONTINUE
c... switch back to cartesian units
         x=matmul(cell%amat,p)/tpi_const
      END IF
      DO 80 j = atoms%jri(n),2,-1
         IF (sx.GE.atoms%rmsh(j,n)) GO TO 90
   80 CONTINUE
   90 jr = j
      CALL ylm4(
     >          atoms%lmax(n),x,
     <          ylm)
      xd1 = cmplx(0.,0.)
      xd2 = cmplx(0.,0.)
      DO l = 0,atoms%lmax(n)
c        if (p(1).eq.0. .and. p(2).eq.0. .and. p(3).eq.0)then
c               write (oUnit,*) 'ff(l,300)=',ff(1,300,1,l)
c               write (oUnit,*) 'ff(l,300)=',ff(1,300,2,l)
c               write (oUnit,*) 'gg(l,300)=',gg(1,300,1,l)
c               write (oUnit,*) 'gg(l,300)=',gg(1,300,2,l)
c        endif
       DO 110 m = -l,l
        lm = l*(l+1)+m
        s = ylm(lm+1)*(ImagUnit)**l
c       if (p(1).eq.0. .and. p(2).eq.0. .and. p(3).eq.0)then
c              write (oUnit,*) 'acof=',acof(lm,1)
c              write (oUnit,*) 'bcof=',bcof(lm,1)
c       endif
        xd1 = xd1 + (acof(lm,na)*cmplx(ff(n,jr,1,l),0.)+
     +               bcof(lm,na)*cmplx(gg(n,jr,1,l),0.))*s/
     /               (atoms%rmsh(jr,n)) 
c    /               (atoms%rmsh(jr,n)*atoms%rmsh(jr,n))
        IF (jr.EQ.1) GO TO 110
        xd2 = xd2 + (acof(lm,na)*cmplx(ff(n,jr+1,1,l),0.)+
     +               bcof(lm,na)*cmplx(gg(n,jr+1,1,l),0.))*s/  
     /               (atoms%rmsh(jr+1,n))
c    /               (atoms%rmsh(jr+1,n)*atoms%rmsh(jr+1,n))
  110  CONTINUE
      ENDDO
c..contributions from the local orbitals
      IF (atoms%nlo(n).GE.1) THEN
       DO l = 1,atoms%nlo(n)
        DO 111 m = -atoms%llo(l,n),atoms%llo(l,n)
         lm = atoms%llo(l,n)*(atoms%llo(l,n)+1)+m
         s = ylm(lm+1)*(ImagUnit)**atoms%llo(l,n) 
         xd1 = xd1 + ccof(m,l,na)*flo(n,jr,1,l)*s/
     /               (atoms%rmsh(jr,n))         
         IF (jr.EQ.1) GO TO 111
         xd2 = xd2 + ccof(m,l,na)*flo(n,jr+1,1,l)*s/
     /               (atoms%rmsh(jr+1,n))         
  111   CONTINUE
       ENDDO
      ENDIF    
      IF (jr.EQ.1) THEN
         xdnout = xd1
      ELSE
         xdnout = xd1 + (xd2-xd1) *
     +                  (sx-atoms%rmsh(jr,n)) /
     /                  (atoms%rmsh(jr+1,n)-atoms%rmsh(jr,n))
         
      END IF
      xdnout = xdnout*const2
c     write (oUnit,*) 'mt:p(i)',p(:)
 8000 FORMAT (2f10.6)
c
      call timestop("wann_real")
      RETURN
      END SUBROUTINE wann_real
      END MODULE m_wann_real
