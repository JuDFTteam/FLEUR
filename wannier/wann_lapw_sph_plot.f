!--------------------------------------------------------------------------------
! Copyright (c) 2016 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------

      MODULE m_wann_lapw_sph_plot

c*****************************************************************
c            plot wannierfunction in muffin tins
c            Frank Freimuth, November 2006
c*****************************************************************
      CONTAINS
      SUBROUTINE wann_lapw_sph_plot(ff,gg,flo,acof,bcof,ccof,x,
     >         nlo,jmtd,lmaxd,nlod,llod,lmd,rmsh,lmaxn,llo,jri,
     <         xdnout)

      USE m_ylm

      implicit none
      integer,intent(in)::jmtd,lmaxd,nlod,llod,lmd,lmaxn,nlo
      integer,intent(in)::llo(nlod),jri
      real,intent(in)::rmsh(jmtd)
      real,intent(in)::x(3)
      real,intent(in)::ff(jmtd,0:lmaxd)
      real,intent(in)::gg(jmtd,0:lmaxd)
      real,intent(in)::flo(jmtd,nlod)
      complex,intent(in)::acof(0:lmd)
      complex,intent(in)::bcof(0:lmd)
      complex,intent(in)::ccof(-llod:llod,nlod)
      complex,intent(out)::xdnout

      real sx
      integer i,j,jr,l,m,lm
      complex ylm((lmaxd+1)**2),xd1,xd2,ci,s

      ci=cmplx(0.,1.)
      sx = 0.0
      DO 50 i = 1,3
         sx = sx + x(i)*x(i)
   50 CONTINUE
      sx = sqrt(sx)
      DO 80 j = jri-1,2,-1
         IF (sx.GE.rmsh(j)) GO TO 90
   80 CONTINUE
   90 jr = j
      CALL ylm4(
     >          lmaxn,x,
     <          ylm)
      xd1 = cmplx(0.,0.)
      xd2 = cmplx(0.,0.)
      DO l = 0,lmaxn
       DO 110 m = -l,l
        lm = l*(l+1)+m
        s = ylm(lm+1)*(ci)**l
        xd1 = xd1 + (acof(lm)*cmplx(ff(jr,l),0.)+
     +               bcof(lm)*cmplx(gg(jr,l),0.))*s/
     /               (rmsh(jr)) 
c        print*,"xd1=",xd1
        IF (jr.EQ.1) GO TO 110
        xd2 = xd2 + (acof(lm)*cmplx(ff(jr+1,l),0.)+
     +               bcof(lm)*cmplx(gg(jr+1,l),0.))*s/  
     /               (rmsh(jr+1))

  110  CONTINUE
      ENDDO
c..contributions from the local orbitals
      IF (nlo.GE.1) THEN
       DO l = 1,nlo
        DO 111 m = -llo(l),llo(l)
         lm = llo(l)*(llo(l)+1)+m

         s = ylm(lm+1)*(ci)**l
         xd1 = xd1 + ccof(m,l)*flo(jr,l)*s/
     /               (rmsh(jr))         
 
         IF (jr.EQ.1) GO TO 111
         xd2 = xd2 + ccof(m,l)*flo(jr+1,l)*s/
     /               (rmsh(jr+1))         

 

  111   CONTINUE
       ENDDO
      ENDIF    
      IF (jr.EQ.1) THEN
         xdnout = xd1
      ELSE
         xdnout = xd1 + (xd2-xd1) *
     +                  (sx-rmsh(jr)) / (rmsh(jr+1)-rmsh(jr))
         
      END IF

      END SUBROUTINE wann_lapw_sph_plot
      END MODULE m_wann_lapw_sph_plot
