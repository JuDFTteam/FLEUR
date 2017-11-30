!--------------------------------------------------------------------------------
! Copyright (c) 2016 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------

      MODULE m_od_mkgxyz3
      CONTAINS
      SUBROUTINE od_mkgxyz3(
     >     ndm,jsdm,ng3,jspins,vl,r,
     >     dvz,dvp,dvr,dvzz,dvpp,dvrr,dvpr,dvrz,dvzp,
     <     agrt,agru,agrd, g2rt,g2ru,g2rd, gggrt,gggru,gggrd,
     <     gzgr)

c     by use of cylindrical z,phi,r components of charge density gradients,
c     make the quantities
c     agrt,agru,agrd,g2rt,g2ru,g2rd,gggrt,gggru,gggrd,gzgr
c     used to calculate gradient contribution to xc potential and
c     energy.
c     Have no time so writing after this crazy subroutines of t.asada
c     Y.Mokrousov, 15-16 october 2002

      IMPLICIT NONE

      integer, intent (in) ::  ndm,ng3,jsdm,jspins
      real,    intent (in) ::  vl(ndm,jsdm),r
      real,    intent (in) ::  dvz(ndm,jsdm),dvp(ndm,jsdm)
      real,    intent (in) ::  dvr(ndm,jsdm),dvzz(ndm,jsdm)
      real,    intent (in) ::  dvpp(ndm,jsdm),dvrr(ndm,jsdm)
      real,    intent (in) ::  dvpr(ndm,jsdm),dvrz(ndm,jsdm)
      real,    intent (in) ::  dvzp(ndm,jsdm)

      real,    intent (out) :: agrt(ndm),agru(ndm),agrd(ndm),g2rt(ndm)
      real,    intent (out) :: g2ru(ndm),g2rd(ndm),gggrt(ndm)
      real,    intent (out) :: gggru(ndm),gggrd(ndm)
      real,    intent (out) :: gzgr(ndm)

c     local

      real    :: vlt,dvzt,dvpt,dvrt,dvzzt,dvppt,dvrrt,dvprt,dvrzt,dvzpt
      real    :: vlu,dvzu,dvpu,dvru,dvzzu,dvppu,dvrru,dvpru,dvrzu,dvzpu
      real    :: vld,dvzd,dvpd,dvrd,dvzzd,dvppd,dvrrd,dvprd,dvrzd,dvzpd
      real    :: dagrzt,dagrzd,dagrzu,dagrpt,dagrpd,dagrpu,dagrrt,dagrrd
      real    :: dagrru,dzdz,dzdp,dzdr
      real    :: sml
      integer :: i

c     so what do all these hieroglyphs mean?
c     for the description of input parameters go to vvacxcg_1.F 
c     as far as in a cylindrical case r and phi variables are somehow 
c     connected, all the staff will depend on the r-coordinate of the point
c     at which it is calculated, which is noted by r
c     I use three types of densities: rho(up),rho(down),rho(total)
c     so that all the corresponding variables are noted by 'u','d', and 't'
c     in the end, special case is variable zeta which is an expression:
c     zeta = (rho(up)-rho(down))/rho(total)
c     different types of differential expressions are calculated here:
c     agrt,agru,agrd:
c                  --->  abs(grad(rho))
c     g2rt,g2ru,g2rd:
c                  --->  laplasian(rho)
c     gggrt,gggru,gggrd :
c                  --->  grad(rho)*grad(abs(grad(rho)))
c     gzgr:
c                  --->  grad(zeta)*grad(rho(total))  

      sml = 1.e-14

      do i = 1,ndm
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
      end do

         if (jspins.eq.1) then
            
            
            do 10 i = 1,ng3
               vlu=max(vl(i,1)/2.,sml)
               dvzu=dvz(i,1)/2.
               dvpu=dvp(i,1)/2.
               dvru=dvr(i,1)/2.
               dvzzu=dvzz(i,1)/2.
               dvppu=dvpp(i,1)/2.
               dvrru=dvrr(i,1)/2.
               dvpru=dvpr(i,1)/2.
               dvrzu=dvrz(i,1)/2.
               dvzpu=dvzp(i,1)/2.
               
               vld=vlu
               dvzd=dvzu
               dvpd=dvpu
               dvrd=dvru
               dvzzd=dvzzu
               dvppd=dvppu
               dvrrd=dvrru
               dvprd=dvpru
               dvrzd=dvrzu
               dvzpd=dvzpu
               
               
               vlt = vlu + vld
               
               dvzt = dvzu + dvzd
               dvpt = dvpu + dvpd
               dvrt = dvru + dvrd
               dvzzt = dvzzu + dvzzd
               dvppt = dvppu + dvppd
               dvrrt = dvrru + dvrrd
               dvprt = dvpru + dvprd
               dvrzt = dvrzu + dvrzd
               dvzpt = dvzpu + dvzpd
               
c     agr: abs(grad(ro)), t,u,d for total, up and down.
           agrt(i) = max(sqrt(dvzt**2+(dvpt**2)/(r**2)+dvrt**2),sml)
           agru(i) = max(sqrt(dvzu**2+(dvpu**2)/(r**2)+dvru**2),sml)
           agrd(i) = max(sqrt(dvzd**2+(dvpd**2)/(r**2)+dvrd**2),sml)
c     d(abs(grad(rho)))/dz
               dagrzt = (dvzt*dvzzt+(dvpt*dvzpt)/(r**2)+dvrt*dvrzt)/
     /              agrt(i)
               dagrzu = (dvzu*dvzzu+(dvpu*dvzpu)/(r**2)+dvru*dvrzu)/
     /              agru(i)
               dagrzd = (dvzd*dvzzd+(dvpd*dvzpd)/(r**2)+dvrd*dvrzd)/
     /              agrd(i)
c     d(abs(grad(rho)))/d(phi)
               dagrpt = (dvzt*dvzpt+(dvpt*dvppt)/(r**2)+dvrt*dvprt)/
     /              agrt(i)
               dagrpu = (dvzu*dvzpu+(dvpu*dvppu)/(r**2)+dvru*dvpru)/
     /              agru(i)
               dagrpd = (dvzd*dvzpd+(dvpd*dvppd)/(r**2)+dvrd*dvprd)/
     /              agrd(i)
c     d(abs(grad(rho)))/dr
               dagrrt = (dvzt*dvrzt+dvrt*dvrrt+
     +              (dvpt*dvprt)/(r**2)-(dvpt**2)/(r**3) )/agrt(i)
               dagrru = (dvzu*dvrzu+dvru*dvrru+
     +              (dvpu*dvpru)/(r**2)-(dvpu**2)/(r**3) )/agru(i)
               dagrrd = (dvzd*dvrzd+dvrd*dvrrd+
     +              (dvpd*dvprd)/(r**2)-(dvpd**2)/(r**3) )/agrd(i)
c     grad(rho)*grad(abs(grad(rho)))
               gggrt(i) = dvzt*dagrzt + (dvpt*dagrpt)/(r**2) + dvrt*
     *              dagrrt
               gggru(i) = dvzu*dagrzu + (dvpu*dagrpu)/(r**2) + dvru*
     *              dagrru
               gggrd(i) = dvzd*dagrzd + (dvpd*dagrpd)/(r**2) + dvrd*
     *              dagrrd
c     dzdz=d(zeta)/dz,dzdp=d(zeta)/dp,dzdr=d(zeta)/dr
               dzdz = (dvzu-dvzd)/vlt - (vlu-vld)*dvzt/(vlt**2)
               dzdp = (dvpu-dvpd)/vlt - (vlu-vld)*dvpt/(vlt**2)
               dzdr = (dvru-dvrd)/vlt - (vlu-vld)*dvrt/(vlt**2)
c     gzgr=grad(zeta)*grad(ro)
               gzgr(i) = dzdz*dvzt + (dzdp*dvpt)/(r**2) + dzdr*dvrt
c     g2r: laplasian(rho)
               g2rt(i) = dvzzt + (dvppt)/(r**2) + dvrrt + dvrt/r
               g2ru(i) = dvzzu + (dvppu)/(r**2) + dvrru + dvru/r
               g2rd(i) = dvzzd + (dvppd)/(r**2) + dvrrd + dvrd/r
 10         end do
            
         else
            
            do 20 i = 1,ng3
               vlu=max(vl(i,1),sml)
               dvzu=dvz(i,1)
               dvpu=dvp(i,1)
               dvru=dvr(i,1)
               dvzzu=dvzz(i,1)
               dvppu=dvpp(i,1)
               dvrru=dvrr(i,1)
               dvpru=dvpr(i,1)
               dvrzu=dvrz(i,1)
               dvzpu=dvzp(i,1)
               
               vld=max(vl(i,jspins),sml)
               dvzd=dvz(i,jspins)
               dvpd=dvp(i,jspins)
               dvrd=dvr(i,jspins)
               dvzzd=dvzz(i,jspins)
               dvppd=dvpp(i,jspins)
               dvrrd=dvrr(i,jspins)
               dvprd=dvpr(i,jspins)
               dvrzd=dvrz(i,jspins)
               dvzpd=dvzp(i,jspins)
               
               vlt = vlu + vld
               
               dvzt = dvzu + dvzd
               dvpt = dvpu + dvpd
               dvrt = dvru + dvrd
               dvzzt = dvzzu + dvzzd
               dvppt = dvppu + dvppd
               dvrrt = dvrru + dvrrd
               dvprt = dvpru + dvprd
               dvrzt = dvrzu + dvrzd
               dvzpt = dvzpu + dvzpd
c     agr: abs(grad(ro)), t,u,d for total, up and down.
               agrt(i) = max(sqrt(dvzt**2+(dvpt**2)/(r**2)+dvrt**2),sml)
               agru(i) = max(sqrt(dvzu**2+(dvpu**2)/(r**2)+dvru**2),sml)
               agrd(i) = max(sqrt(dvzd**2+(dvpd**2)/(r**2)+dvrd**2),sml)
c     d(abs(grad(rho)))/dz
               dagrzt = (dvzt*dvzzt+(dvpt*dvzpt)/(r**2)+dvrt*dvrzt)/
     /              agrt(i)
               dagrzu = (dvzu*dvzzu+(dvpu*dvzpu)/(r**2)+dvru*dvrzu)/
     /              agru(i)
               dagrzd = (dvzd*dvzzd+(dvpd*dvzpd)/(r**2)+dvrd*dvrzd)/
     /              agrd(i)
c     d(abs(grad(rho)))/d(phi)
               dagrpt = (dvzt*dvzpt+(dvpt*dvppt)/(r**2)+dvrt*dvprt)/
     /              agrt(i)
               dagrpu = (dvzu*dvzpu+(dvpu*dvppu)/(r**2)+dvru*dvpru)/
     /              agru(i)
               dagrpd = (dvzd*dvzpd+(dvpd*dvppd)/(r**2)+dvrd*dvprd)/
     /              agrd(i)
c     d(abs(grad(rho)))/dr
               dagrrt = (dvzt*dvrzt+dvrt*dvrrt+
     +              (dvpt*dvprt)/(r**2)-(dvpt**2)/(r**3) )/agrt(i)
               dagrru = (dvzu*dvrzu+dvru*dvrru+
     +              (dvpu*dvpru)/(r**2)-(dvpu**2)/(r**3) )/agru(i)
               dagrrd = (dvzd*dvrzd+dvrd*dvrrd+
     +              (dvpd*dvprd)/(r**2)-(dvpd**2)/(r**3) )/agrd(i)
c     grad(rho)*grad(abs(grad(rho)))
               gggrt(i) = dvzt*dagrzt + (dvpt*dagrpt)/(r**2) + dvrt*
     *              dagrrt
               gggru(i) = dvzu*dagrzu + (dvpu*dagrpu)/(r**2) + dvru*
     *              dagrru
               gggrd(i) = dvzd*dagrzd + (dvpd*dagrpd)/(r**2) + dvrd*
     *              dagrrd
c     dzdz=d(zeta)/dz,dzdp=d(zeta)/dp,dzdr=d(zeta)/dr
               dzdz = (dvzu-dvzd)/vlt - (vlu-vld)*dvzt/(vlt**2)
               dzdp = (dvpu-dvpd)/vlt - (vlu-vld)*dvpt/(vlt**2)
               dzdr = (dvru-dvrd)/vlt - (vlu-vld)*dvrt/(vlt**2)
c     gzgr=grad(zeta)*grad(ro)
               gzgr(i) = dzdz*dvzt + (dzdp*dvpt)/(r**2) + dzdr*dvrt
c     g2r: laplasian(rho)
               g2rt(i) = dvzzt + (dvppt)/(r**2) + dvrrt + dvrt/r
               g2ru(i) = dvzzu + (dvppu)/(r**2) + dvrru + dvru/r
               g2rd(i) = dvzzd + (dvppd)/(r**2) + dvrrd + dvrd/r
               
 20         end do
            
         end if
      
      END SUBROUTINE od_mkgxyz3
      END MODULE m_od_mkgxyz3
