      MODULE m_mkgz
c.....------------------------------------------------------------------
c     by use of cartesian x,y,z components of charge density gradients,
c     make the quantities
cc      agrt,agru,agrd,g2rt,g2ru,g2rd,gggrt,gggru,gggrd,
cc      gzgr
cc    used to calculate gradient contribution to xc potential and
cc    energy.
c.....------------------------------------------------------------------
      CONTAINS 
      SUBROUTINE mkgz(
     >                nmzdf,jspins,rh1,rh2,rhdz1,rhdz2,rhdzz1,rhdzz2,
     <          agrt,agru,agrd,g2rt,g2ru,g2rd,gggrt,gggru,gggrd,gzgr)

      IMPLICIT NONE
      INTEGER, INTENT (IN) :: nmzdf,jspins
      REAL,    INTENT (IN) :: rh1(nmzdf),rhdz1(nmzdf),rhdzz1(nmzdf)
      REAL,    INTENT (IN) :: rh2(nmzdf),rhdz2(nmzdf),rhdzz2(nmzdf)
      REAL,    INTENT (OUT) :: agrt(nmzdf),agru(nmzdf),agrd(nmzdf) 
      REAL,    INTENT (OUT) :: g2rt(nmzdf),g2ru(nmzdf),g2rd(nmzdf)
      REAL,    INTENT (OUT) :: gggrt(nmzdf),gggru(nmzdf),gggrd(nmzdf)
      REAL,    INTENT (OUT) :: gzgr(nmzdf)

      INTEGER i
      REAL    vlt,dvzt,dvzzt,vlu,dvzu,dvzzu,vld,dvzd,dvzzd
      REAL    dagrzt,dagrzu,dagrzd,dztadz,sml

      sml = 1.e-14

      IF (jspins == 1) THEN

        DO i = 1,nmzdf

          vlu = max(rh1(i)/2,sml)
          dvzu = rhdz1(i)/2
          dvzzu = rhdzz1(i)/2
          vld = vlu
          dvzd = dvzu
          dvzzd = dvzzu

          vlt = vlu+vld
          dvzt = dvzu+dvzd
          dvzzt = dvzzu+dvzzd

c         agrt: abs(grad(ro)), u,d for up and down.

          agrt(i) = max(abs(dvzt),sml)
          agru(i) = max(abs(dvzu),sml)
          agrd(i) = max(abs(dvzd),sml)

          dagrzt= dvzt*dvzzt/agrt(i)
          dagrzu= dvzu*dvzzu/agru(i)
          dagrzd= dvzd*dvzzd/agrd(i)

          gggrt(i) = dvzt*dagrzt
          gggru(i) = dvzu*dagrzu
          gggrd(i) = dvzd*dagrzd

c         dztadz=d(zeta)/dz,..

          dztadz = (dvzu-dvzd)/vlt - (vlu-vld)*dvzt/vlt**2

c         gzgr=grad(zeta)*grad(ro).

          gzgr(i) = dztadz*dvzt

c         g2rt: grad(grad(ro))

          g2rt(i)  = dvzzt
          g2ru(i)  = dvzzu
          g2rd(i)  = dvzzd

        ENDDO

      ELSE

        DO i = 1,nmzdf

          vlu = max(rh1(i),sml)
          dvzu = rhdz1(i)
          dvzzu = rhdzz1(i)
          vld = max(rh2(i),sml)
          dvzd = rhdz2(i)
          dvzzd = rhdzz2(i)

          vlt = vlu+vld
          dvzt = dvzu+dvzd
          dvzzt = dvzzu+dvzzd

c         agrt: abs(grad(ro)), u,d for up and down.

          agrt(i) = max(abs(dvzt),sml)
          agru(i) = max(abs(dvzu),sml)
          agrd(i) = max(abs(dvzd),sml)

          dagrzt= dvzt*dvzzt/agrt(i)
          dagrzu= dvzu*dvzzu/agru(i)
          dagrzd= dvzd*dvzzd/agrd(i)

          gggrt(i) = dvzt*dagrzt
          gggru(i) = dvzu*dagrzu
          gggrd(i) = dvzd*dagrzd

c         dztadz=d(zeta)/dz,..

          dztadz = (dvzu-dvzd)/vlt - (vlu-vld)*dvzt/vlt**2

c         gzgr=grad(zeta)*grad(ro).

          gzgr(i) = dztadz*dvzt

c         g2rt: grad(grad(ro))

          g2rt(i)  = dvzzt
          g2ru(i)  = dvzzu
          g2rd(i)  = dvzzd


        ENDDO
      ENDIF

      END SUBROUTINE mkgz
      END MODULE m_mkgz
