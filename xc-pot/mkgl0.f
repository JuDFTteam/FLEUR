      MODULE m_mkgl0
c.....------------------------------------------------------------------
c     make quantities for vxcallg. for paramag. case
c       dens,drr,ddrr are charge density and its gradient for nonmag.
cc      one spin.

c     agr: abs(grad(ro)), g2r: laplacian(ro),
cc    gggr: grad(ro)*grad(agr),
cc    grgru,d: grad(ro)*grad(rou),for rod., gzgr: grad(zeta)*grad(ro).
c.....------------------------------------------------------------------
      CONTAINS
      SUBROUTINE mkgl0(
     >                 mshd,msh,jspd,jspins,rad, densi,drri,ddrri, 
     <                 agr,agru,agrd, g2r,g2ru,g2rd,
     <                 gggr,gggru,gggrd, grgru,grgrd, gzgr)

      IMPLICIT NONE

      INTEGER, INTENT (IN) :: mshd,msh,jspd,jspins
      REAL,    INTENT (IN) :: rad(mshd),densi(mshd,jspd)
      REAL,    INTENT (IN) :: drri(mshd,jspd),ddrri(mshd,jspd)
      REAL,    INTENT (OUT) :: agr(mshd),agru(mshd),agrd(mshd)
      REAL,    INTENT (OUT) :: g2r(mshd),g2ru(mshd),g2rd(mshd)
      REAL,    INTENT (OUT) :: gggr(mshd),gggru(mshd),gggrd(mshd)
      REAL,    INTENT (OUT) :: grgru(mshd),grgrd(mshd),gzgr(mshd)
c
      REAL dagrr,dagrru,ddrr,ddrrd,ddrru,drr,drrd,drru,dzdr,ro,
     +     rod,rou,rv,spnf
      INTEGER i
c.....------------------------------------------------------------------
c     ..
      spnf = 1./(3-jspins)
      DO i = 1,msh

          rv = rad(i)
          ro =  spnf * (densi(i,1) + densi(i,jspins)) 
          rou = spnf * densi(i,1)
          rod = spnf * densi(i,jspins) 

          drr =   spnf * (drri(i,1) + drri(i,jspins))
          drru =  spnf * drri(i,1)
          drrd =  spnf * drri(i,jspins)
          ddrr =  spnf * (ddrri(i,1) + ddrri(i,jspins))
          ddrru = spnf * ddrri(i,1)
          ddrrd = spnf * ddrri(i,jspins)

          agr(i) = abs(drr)
          agru(i) = abs(drru)
          agrd(i) = agru(i)

          dagrr = drr*ddrr/agr(i)
          dagrru = drru*ddrru/agru(i)
 
          gggr(i) = drr*dagrr
          gggru(i) = drru*dagrru
          gggrd(i) = gggru(i)

          dzdr = ((drru-drrd)*ro- (rou-rod)*drr)/ro**2

          gzgr(i) = dzdr*drr

          g2r(i) = ddrr + 2*drr/rv
          g2ru(i) = ddrru + 2*drru/rv
          g2rd(i) = ddrrd + 2*drrd/rv

          grgru(i) = drr*drru
          grgrd(i) = drr*drrd

      ENDDO

      END SUBROUTINE mkgl0
      END MODULE m_mkgl0
