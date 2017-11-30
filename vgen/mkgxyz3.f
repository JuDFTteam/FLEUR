      MODULE m_mkgxyz3
c.....------------------------------------------------------------------
c     by use of cartesian x,y,z components of charge density gradients,
c     make the quantities
cc      agrt,agru,agrd,g2rt,g2ru,g2rd,gggrt,gggru,gggrd,gzgr
cc    used to calculate gradient contribution to xc potential and
cc    energy.
c.....------------------------------------------------------------------
      CONTAINS
      SUBROUTINE mkgxyz3(
     >                   ndm,jsdm,ng3,jspins,vl,
     >                   dvx,dvy,dvz,dvxx,dvyy,dvzz,dvyz,dvzx,dvxy,
     <                   agrt,agru,agrd, g2rt,g2ru,g2rd, 
     <                   gggrt,gggru,gggrd,gzgr)

      IMPLICIT NONE
      INTEGER, INTENT (IN) :: ndm,ng3,jsdm,jspins
      REAL,    INTENT (IN) :: vl(ndm,jsdm)
      REAL,    INTENT (IN) :: dvx(ndm,jsdm),dvy(ndm,jsdm),dvz(ndm,jsdm)
      REAL, INTENT (IN) :: dvxx(ndm,jsdm),dvyy(ndm,jsdm),dvzz(ndm,jsdm)
      REAL, INTENT (IN) :: dvyz(ndm,jsdm),dvzx(ndm,jsdm),dvxy(ndm,jsdm)

      REAL, INTENT (OUT) :: agrt(ndm),agru(ndm),agrd(ndm)
      REAL, INTENT (OUT) :: g2rt(ndm),g2ru(ndm),g2rd(ndm)
      REAL, INTENT (OUT) :: gggrt(ndm),gggru(ndm),gggrd(ndm),gzgr(ndm)

      REAL vlt,dvxt,dvyt,dvzt,dvxxt,dvyyt,dvzzt,dvyzt,dvzxt,dvxyt,
     &     vlu,dvxu,dvyu,dvzu,dvxxu,dvyyu,dvzzu,dvyzu,dvzxu,dvxyu,
     &     vld,dvxd,dvyd,dvzd,dvxxd,dvyyd,dvzzd,dvyzd,dvzxd,dvxyd,
     &     dagrxt,dagrxd,dagrxu,dagryt,dagryd,dagryu,dagrzt,dagrzd,
     +     dagrzu,dzdx,dzdy,dzdz,
     +     sml
      INTEGER i

      sml = 1.e-14

      DO i = 1,ndm
          agrt(i) = 0.0
          agru(i) = 0.0
          agrd(i) = 0.0
          gggrt(i) = 0.0
          gggru(i) = 0.0
          gggrd(i) = 0.0
          gzgr(i) = 0.0
          g2rt(i) = 0.0
          g2ru(i) = 0.0
          g2rd(i) = 0.0
      ENDDO

      IF (jspins.eq.1) THEN

        DO 10 i = 1,ng3

          vlu=max(vl(i,1)/2,sml)
          dvxu=dvx(i,1)/2
          dvyu=dvy(i,1)/2
          dvzu=dvz(i,1)/2
          dvxxu=dvxx(i,1)/2
          dvyyu=dvyy(i,1)/2
          dvzzu=dvzz(i,1)/2
          dvyzu=dvyz(i,1)/2
          dvzxu=dvzx(i,1)/2
          dvxyu=dvxy(i,1)/2

          vld=vlu
          dvxd=dvxu
          dvyd=dvyu
          dvzd=dvzu
          dvxxd=dvxxu
          dvyyd=dvyyu
          dvzzd=dvzzu
          dvyzd=dvyzu
          dvzxd=dvzxu
          dvxyd=dvxyu


          vlt = vlu + vld

          dvxt = dvxu + dvxd
          dvyt = dvyu + dvyd
          dvzt = dvzu + dvzd
          dvxxt = dvxxu + dvxxd
          dvyyt = dvyyu + dvyyd
          dvzzt = dvzzu + dvzzd
          dvyzt = dvyzu + dvyzd
          dvzxt = dvzxu + dvzxd
          dvxyt = dvxyu + dvxyd

c         agr: abs(grad(ro)), t,u,d for total, up and down.

          agrt(i) = max(sqrt(dvxt**2+dvyt**2+dvzt**2),sml)
          agru(i) = max(sqrt(dvxu**2+dvyu**2+dvzu**2),sml)
          agrd(i) = max(sqrt(dvxd**2+dvyd**2+dvzd**2),sml)

          dagrxt = (dvxt*dvxxt+dvyt*dvxyt+dvzt*dvzxt)/agrt(i)
          dagrxu = (dvxu*dvxxu+dvyu*dvxyu+dvzu*dvzxu)/agru(i)
          dagrxd = (dvxd*dvxxd+dvyd*dvxyd+dvzd*dvzxd)/agrd(i)

          dagryt = (dvxt*dvxyt+dvyt*dvyyt+dvzt*dvyzt)/agrt(i)
          dagryu = (dvxu*dvxyu+dvyu*dvyyu+dvzu*dvyzu)/agru(i)
          dagryd = (dvxd*dvxyd+dvyd*dvyyd+dvzd*dvyzd)/agrd(i)

          dagrzt = (dvxt*dvzxt+dvyt*dvyzt+dvzt*dvzzt)/agrt(i)
          dagrzu = (dvxu*dvzxu+dvyu*dvyzu+dvzu*dvzzu)/agru(i)
          dagrzd = (dvxd*dvzxd+dvyd*dvyzd+dvzd*dvzzd)/agrd(i)

          gggrt(i) = dvxt*dagrxt + dvyt*dagryt + dvzt*dagrzt
          gggru(i) = dvxu*dagrxu + dvyu*dagryu + dvzu*dagrzu
          gggrd(i) = dvxd*dagrxd + dvyd*dagryd + dvzd*dagrzd

c         dzdx=d(zeta)/dx,..

          dzdx = (dvxu-dvxd)/vlt - (vlu-vld)*dvxt/vlt**2
          dzdy = (dvyu-dvyd)/vlt - (vlu-vld)*dvyt/vlt**2
          dzdz = (dvzu-dvzd)/vlt - (vlu-vld)*dvzt/vlt**2

c         gzgr=grad(zeta)*grad(ro).

          gzgr(i) = dzdx*dvxt + dzdy*dvyt + dzdz*dvzt

c         g2r: grad(grad(ro))

          g2rt(i) = dvxxt + dvyyt + dvzzt
          g2ru(i) = dvxxu + dvyyu + dvzzu
          g2rd(i) = dvxxd + dvyyd + dvzzd


  10    ENDDO

      ELSE

        DO 20 i = 1,ng3

          vlu = max(vl(i,1),sml)
          dvxu=dvx(i,1)
          dvyu=dvy(i,1)
          dvzu=dvz(i,1)
          dvxxu=dvxx(i,1)
          dvyyu=dvyy(i,1)
          dvzzu=dvzz(i,1)
          dvyzu=dvyz(i,1)
          dvzxu=dvzx(i,1)
          dvxyu=dvxy(i,1)

          vld = max(vl(i,jspins),sml)
          dvxd=dvx(i,jspins)
          dvyd=dvy(i,jspins)
          dvzd=dvz(i,jspins)
          dvxxd=dvxx(i,jspins)
          dvyyd=dvyy(i,jspins)
          dvzzd=dvzz(i,jspins)
          dvyzd=dvyz(i,jspins)
          dvzxd=dvzx(i,jspins)
          dvxyd=dvxy(i,jspins)

          vlt = vlu + vld

          dvxt = dvxu + dvxd
          dvyt = dvyu + dvyd
          dvzt = dvzu + dvzd
          dvxxt = dvxxu + dvxxd
          dvyyt = dvyyu + dvyyd
          dvzzt = dvzzu + dvzzd
          dvyzt = dvyzu + dvyzd
          dvzxt = dvzxu + dvzxd
          dvxyt = dvxyu + dvxyd

c         agr: abs(grad(ro)), t,u,d for total, up and down.

          agrt(i) = max(sqrt(dvxt**2+dvyt**2+dvzt**2),sml)
          agru(i) = max(sqrt(dvxu**2+dvyu**2+dvzu**2),sml)
          agrd(i) = max(sqrt(dvxd**2+dvyd**2+dvzd**2),sml)

          dagrxt = (dvxt*dvxxt+dvyt*dvxyt+dvzt*dvzxt)/agrt(i)
          dagrxu = (dvxu*dvxxu+dvyu*dvxyu+dvzu*dvzxu)/agru(i)
          dagrxd = (dvxd*dvxxd+dvyd*dvxyd+dvzd*dvzxd)/agrd(i)

          dagryt = (dvxt*dvxyt+dvyt*dvyyt+dvzt*dvyzt)/agrt(i)
          dagryu = (dvxu*dvxyu+dvyu*dvyyu+dvzu*dvyzu)/agru(i)
          dagryd = (dvxd*dvxyd+dvyd*dvyyd+dvzd*dvyzd)/agrd(i)

          dagrzt = (dvxt*dvzxt+dvyt*dvyzt+dvzt*dvzzt)/agrt(i)
          dagrzu = (dvxu*dvzxu+dvyu*dvyzu+dvzu*dvzzu)/agru(i)
          dagrzd = (dvxd*dvzxd+dvyd*dvyzd+dvzd*dvzzd)/agrd(i)

          gggrt(i) = dvxt*dagrxt + dvyt*dagryt + dvzt*dagrzt
          gggru(i) = dvxu*dagrxu + dvyu*dagryu + dvzu*dagrzu
          gggrd(i) = dvxd*dagrxd + dvyd*dagryd + dvzd*dagrzd

c         dzdx=d(zeta)/dx,..

          dzdx = (dvxu-dvxd)/vlt -  (vlu-vld)*dvxt/vlt**2
          dzdy = (dvyu-dvyd)/vlt -  (vlu-vld)*dvyt/vlt**2
          dzdz = (dvzu-dvzd)/vlt -  (vlu-vld)*dvzt/vlt**2

c         gzgr=grad(zeta)*grad(ro).

          gzgr(i) = dzdx*dvxt + dzdy*dvyt + dzdz*dvzt

c         g2r: grad(grad(ro))

          g2rt(i) = dvxxt + dvyyt + dvzzt
          g2ru(i) = dvxxu + dvyyu + dvzzu
          g2rd(i) = dvxxd + dvyyd + dvzzd

  20    ENDDO

      ENDIF
    
      END SUBROUTINE mkgxyz3
      END MODULE m_mkgxyz3
