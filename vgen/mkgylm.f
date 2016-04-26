      MODULE m_mkgylm
      CONTAINS
      SUBROUTINE mkgylm(
     >                  jspins,rv,thet,nsp,nspd,jspd,
     >                  rh,rhdr,rhdt,rhdf,rhdrr,rhdtt,rhdff,
     >                  rhdtf,rhdrt,rhdrf,
     <                  agr,agru,agrd,g2r,g2ru,g2rd,gggr,
     <                  gggru,gggrd,grgru,grgrd,gzgr)
c.....------------------------------------------------------------------
c     by use of charge density and its polar coord. gradient components
cc    calculate agr and others used to evaluate gradient
cc    contributions to potential and energy. t.a. 1996.
c.....------------------------------------------------------------------
c     ro=sum(ro*ylh), rdr=sum(drr*ylh), drdr=sum(ddrr*ylh),
cc    rdt=sum(ro*ylht1), rdtt=sum(ro*ylht2), ...
cc    rdf=sum(ro*ylhf1), rdff=sum(ro*ylhf2), ...
cc    rdtf=sum(ro*ylhtf), rdff=sum(ro*ylhf2), ...
cc    rdrt=sum(drr*ylht1),rdrf=sum(drr*ylhf1),

c     agr: abs(grad(ro)), g2r: laplacian(ro),
cc    gggr: grad(ro)*grad(agr),
cc    grgru,d: grad(ro)*grad(rou),for rod., gzgr: grad(zeta)*grad(ro).

c     dagrr,-t,-f: d(agr)/dr, d(agr)/dth/r, d(agr)/dfi/r/sint.
c.....------------------------------------------------------------------
      IMPLICIT NONE
C  ..
C  .. Arguments ..
      INTEGER, INTENT (IN) :: nspd,jspd
      REAL,    INTENT (IN) :: rv
      REAL,    INTENT (IN) :: thet(nspd)
      REAL,    INTENT (IN) :: rh(nspd,jspd),rhdr(nspd,jspd)
      REAL,    INTENT (IN) :: rhdf(nspd,jspd),rhdrr(nspd,jspd)
      REAL,    INTENT (IN) :: rhdtt(nspd,jspd),rhdff(nspd,jspd)
      REAL,    INTENT (IN) :: rhdtf(nspd,jspd),rhdrt(nspd,jspd)
      REAL,    INTENT (IN) :: rhdrf(nspd,jspd),rhdt(nspd,jspd)
      REAL,    INTENT (OUT) :: agr(nspd),agru(nspd),agrd(nspd)
      REAL,    INTENT (OUT) :: g2ru(nspd),g2rd(nspd),gggr(nspd)
      REAL,    INTENT (OUT) :: gggru(nspd),gzgr(nspd),g2r(nspd)
      REAL,    INTENT (OUT) :: gggrd(nspd),grgru(nspd),grgrd(nspd)
C  ..
C  .. Locals ..
      INTEGER i,jspins,nsp
      REAL    chsml, dagrf,dagrfd,dagrfu,dagrr,dagrrd,dagrru,dagrt,
     +        dagrtd,dagrtu,drdr,dzdfs,dzdr,dzdtr,grf,grfd,grfu,grr,
     +        grrd,grru,grt,grtd,grtu,rdf,rdfd,rdff,rdffd,rdffu,rdfu,
     +        rdr,rdrd,rdrf,rdrfd,rdrfu,rdrrd,rdrru,rdrt,rdrtd,rdrtu,
     +        rdru,rdt,rdtd,rdtf,rdtfd,rdtfu,rdtt,rdttd,rdttu,rdtu,
     +        ro,ro2,rod,rou,rv1,rv2,rv3,rvsin1,sint1,sint2,tant1

      chsml = 1.e-10

      rv1 = rv
      rv2 = rv1**2
      rv3 = rv1**3

      IF (jspins.EQ.1) THEN

        points_1 : DO i = 1,nsp

          agr(i)   = 0.0
          agru(i)  = 0.0
          agrd(i)  = 0.0
          g2r(i)   = 0.0
          g2ru(i)  = 0.0
          g2rd(i)  = 0.0
          gggr(i)  = 0.0
          gggru(i) = 0.0
          gggrd(i) = 0.0
          grgru(i) = 0.0
          grgrd(i) = 0.0
          gzgr(i)  = 0.0

          ro = rh(i,1)

          IF (ro.LT.chsml) CYCLE points_1

          sint1 = sin(thet(i))
          sint2 = sint1**2
          tant1 = tan(thet(i))
          rvsin1 = rv1*sint1


          rou   = ro/2
          rdru  = rhdr(i,1)/2
          rdtu  = rhdt(i,1)/2
          rdfu  = rhdf(i,1)/2
          rdrru = rhdrr(i,1)/2
          rdttu = rhdtt(i,1)/2
          rdffu = rhdff(i,1)/2
          rdtfu = rhdtf(i,1)/2
          rdrtu = rhdrt(i,1)/2
          rdrfu = rhdrf(i,1)/2

          rod = rou
          rdrd = rdru
          rdtd = rdtu
          rdfd = rdfu
          rdrrd = rdrru
          rdttd = rdttu
          rdffd = rdffu
          rdtfd = rdtfu
          rdrtd = rdrtu
          rdrfd = rdrfu

          rdr = rdru + rdrd
          rdt = rdtu + rdtd
          rdf = rdfu + rdfd
          drdr = rdrru + rdrrd
          rdtt = rdttu + rdttd
          rdff = rdffu + rdffd
          rdtf = rdtfu + rdtfd
          rdrt = rdrtu + rdrtd
          rdrf = rdrfu + rdrfd

          ro2 = ro**2

          grr = rdr
          grt = rdt/rv1
          grf = rdf/rvsin1

          agr(i) = sqrt(grr**2+grt**2+grf**2)

          IF (agr(i).LT.chsml) CYCLE points_1

          dagrr = (rdr*drdr*rv3+rdt* (rdrt*rv1-rdt)+
     +            rdf* (rdrf*rv1-rdf)/sint2)/agr(i)/rv3

          dagrt = (rdr*rdrt*rv2+rdt*rdtt+rdf* (-rdf/tant1+rdtf)/sint2)/
     +            (agr(i)*rv3)

          dagrf = (rdr*rdrf*rv2+rdt*rdtf+rdf*rdff/sint2)/
     +            (agr(i)*rv3*sint1)

          g2r(i) = drdr + 2.e0*rdr/rv1 + (rdtt+rdt/tant1+rdff/sint2)/rv2


          dzdr = ((rdru-rdrd)*ro- (rou-rod)*rdr)/ro2

!-->      dzdtr,dzdfs vanish by definition.

          dzdtr = 0.0
          dzdfs = 0.0

          gggr(i) = grr*dagrr + grt*dagrt + grf*dagrf

          gzgr(i) = dzdr*grr + dzdtr*grt + dzdfs*grf

          grru = rdru
          grtu = rdtu/rv1
          grfu = rdfu/rvsin1

          agru(i) = sqrt(grru**2+grtu**2+grfu**2)

          dagrru = (rdru*rdrru*rv3+rdtu* (rdrtu*rv1-rdtu)+
     +             rdfu* (rdrfu*rv1-rdfu)/sint2)/agru(i)/rv3

          dagrtu = (rdru*rdrtu*rv2+rdtu*rdttu+
     +             rdfu* (-rdfu/tant1+rdtfu)/sint2)/ (agru(i)*rv3)

          dagrfu = (rdru*rdrfu*rv2+rdtu*rdtfu+rdfu*rdffu/sint2)/
     +             (agru(i)*rv3*sint1)

          g2ru(i) = rdrru + 2.e0*rdru/rv1 +
     +              (rdttu+rdtu/tant1+rdffu/sint2)/rv2

          gggru(i) = grru*dagrru + grtu*dagrtu + grfu*dagrfu

          grgru(i) = grr*grru + grt*grtu + grf*grfu

          grrd = rdrd
          grtd = rdtd/rv1
          grfd = rdfd/rvsin1

          agrd(i) = sqrt(grrd**2+grtd**2+grfd**2)

          dagrrd = (rdrd*rdrrd*rv3+rdtd* (rdrtd*rv1-rdtd)+
     +             rdfd* (rdrfd*rv1-rdfd)/sint2)/agrd(i)/rv3

          dagrtd = (rdrd*rdrtd*rv2+rdtd*rdttd+
     +             rdfd* (-rdfd/tant1+rdtfd)/sint2)/ (agrd(i)*rv3)

          dagrfd = (rdrd*rdrfd*rv2+rdtd*rdtfd+rdfd*rdffd/sint2)/
     +             (agrd(i)*rv3*sint1)

          g2rd(i) = rdrrd + 2*rdrd/rv1 +
     +              (rdttd+rdtd/tant1+rdffd/sint2)/rv2

          gggrd(i) = grrd*dagrrd + grtd*dagrtd + grfd*dagrfd

          grgrd(i) = grr*grrd + grt*grtd + grf*grfd


        ENDDO points_1

      ELSE

        points_2 : DO i = 1,nsp

          agr(i) = 0.0
          agru(i) = 0.0
          agrd(i) = 0.0
          g2r(i) = 0.0
          g2ru(i) = 0.0
          g2rd(i) = 0.0
          gggr(i) = 0.0
          gggru(i) = 0.0
          gggrd(i) = 0.0
          grgru(i) = 0.0
          grgrd(i) = 0.0
          gzgr(i) = 0.0

          ro = rh(i,1) + rh(i,jspins)

          IF (ro.LT.chsml) CYCLE points_2

          sint1 = sin(thet(i))
          sint2 = sint1**2
          tant1 = tan(thet(i))
          rvsin1 = rv1*sint1

          rou = rh(i,1)
          rdru = rhdr(i,1)
          rdtu = rhdt(i,1)
          rdfu = rhdf(i,1)
          rdrru = rhdrr(i,1)
          rdttu = rhdtt(i,1)
          rdffu = rhdff(i,1)
          rdtfu = rhdtf(i,1)
          rdrtu = rhdrt(i,1)
          rdrfu = rhdrf(i,1)

          rod = rh(i,jspins)
          rdrd = rhdr(i,jspins)
          rdtd = rhdt(i,jspins)
          rdfd = rhdf(i,jspins)
          rdrrd = rhdrr(i,jspins)
          rdttd = rhdtt(i,jspins)
          rdffd = rhdff(i,jspins)
          rdtfd = rhdtf(i,jspins)
          rdrtd = rhdrt(i,jspins)
          rdrfd = rhdrf(i,jspins)

          rdr = rdru + rdrd
          rdt = rdtu + rdtd
          rdf = rdfu + rdfd
          drdr = rdrru + rdrrd
          rdtt = rdttu + rdttd
          rdff = rdffu + rdffd
          rdtf = rdtfu + rdtfd
          rdrt = rdrtu + rdrtd
          rdrf = rdrfu + rdrfd

          ro2 = ro**2

          grr = rdr
          grt = rdt/rv1
          grf = rdf/rvsin1

          agr(i) = sqrt(grr**2+grt**2+grf**2)

          IF (agr(i).LT.chsml) CYCLE points_2

          dagrr = (rdr*drdr*rv3+rdt* (rdrt*rv1-rdt)+
     +            rdf* (rdrf*rv1-rdf)/sint2)/agr(i)/rv3

          dagrt = (rdr*rdrt*rv2+rdt*rdtt+rdf* (-rdf/tant1+rdtf)/sint2)/
     +            (agr(i)*rv3)

          dagrf = (rdr*rdrf*rv2+rdt*rdtf+rdf*rdff/sint2)/
     +            (agr(i)*rv3*sint1)

          g2r(i) = drdr + 2.e0*rdr/rv1 + (rdtt+rdt/tant1+rdff/sint2)/rv2

          dzdr = ((rdru-rdrd)*ro- (rou-rod)*rdr)/ro2

c         dzdtr,dzdfs vanish by definition.
          dzdtr = 0.0
          dzdfs = 0.0

          gggr(i) = grr*dagrr + grt*dagrt + grf*dagrf

          gzgr(i) = dzdr*grr + dzdtr*grt + dzdfs*grf

          grru = rdru
          grtu = rdtu/rv1
          grfu = rdfu/rvsin1

          agru(i) = sqrt(grru**2+grtu**2+grfu**2)

          dagrru = (rdru*rdrru*rv3+rdtu* (rdrtu*rv1-rdtu)+
     +             rdfu* (rdrfu*rv1-rdfu)/sint2)/agru(i)/rv3

          dagrtu = (rdru*rdrtu*rv2+rdtu*rdttu+
     +             rdfu* (-rdfu/tant1+rdtfu)/sint2)/ (agru(i)*rv3)

          dagrfu = (rdru*rdrfu*rv2+rdtu*rdtfu+rdfu*rdffu/sint2)/
     +             (agru(i)*rv3*sint1)

          g2ru(i) = rdrru + 2.e0*rdru/rv1 +
     +              (rdttu+rdtu/tant1+rdffu/sint2)/rv2

          gggru(i) = grru*dagrru + grtu*dagrtu + grfu*dagrfu

          grgru(i) = grr*grru + grt*grtu + grf*grfu


          grrd = rdrd
          grtd = rdtd/rv1
          grfd = rdfd/rvsin1

          agrd(i) = sqrt(grrd**2+grtd**2+grfd**2)

          dagrrd = (rdrd*rdrrd*rv3+rdtd* (rdrtd*rv1-rdtd)+
     +             rdfd* (rdrfd*rv1-rdfd)/sint2)/agrd(i)/rv3

          dagrtd = (rdrd*rdrtd*rv2+rdtd*rdttd+
     +             rdfd* (-rdfd/tant1+rdtfd)/sint2)/ (agrd(i)*rv3)

          dagrfd = (rdrd*rdrfd*rv2+rdtd*rdtfd+rdfd*rdffd/sint2)/
     +             (agrd(i)*rv3*sint1)



          g2rd(i) = rdrrd + 2*rdrd/rv1 +
     +              (rdttd+rdtd/tant1+rdffd/sint2)/rv2

          gggrd(i) = grrd*dagrrd + grtd*dagrtd + grfd*dagrfd

          grgrd(i) = grr*grrd + grt*grtd + grf*grfd


        ENDDO points_2

      ENDIF

      RETURN
      END SUBROUTINE mkgylm
      END MODULE m_mkgylm
