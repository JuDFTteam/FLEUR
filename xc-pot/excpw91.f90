MODULE m_excpw91
!.....-----------------------------------------------------------------
!..... pw91 exchange-correlation energy density in hartree.
!.....------------------------------------------------------------------
CONTAINS
   SUBROUTINE excpw91( &
      jspins,mirm,irmx,rh,agr,agru,agrd, &
      g2r,g2ru,g2rd,gggr,gggru,gggrd,gzgr, &
      exc, &
      idsprs,isprsv,sprsv)
      USE m_corl91
      USE m_corg91
      USE m_xch91
      USE m_constants

      IMPLICIT NONE

! .. Arguments ..
      INTEGER, INTENT (IN) :: jspins,irmx,mirm,isprsv,idsprs
      REAL,    INTENT (IN) :: sprsv

      REAL,    INTENT (IN) :: rh(mirm,jspins)
      REAL,    INTENT (IN) :: agr(mirm),agru(mirm),agrd(mirm)
      REAL,    INTENT (IN) :: g2r(mirm),g2ru(mirm),g2rd(mirm)
      REAL,    INTENT (IN) :: gggr(mirm),gggru(mirm)
      REAL,    INTENT (IN) :: gggrd(mirm),gzgr(mirm)
      REAL,    INTENT (OUT) :: exc(mirm)

! .. local variables
      REAL :: ro,zta,alf,alfc,c13,c23,c43,c53, &
              cedg,cedl,dbrod,dbrou,dsprs,ec,ecrs,eczta,fk,gz, &
              ro13,ro2,rod,rod3,rod43,rod53,rou,rou3,rou43,rou53, &
              rs,sd,sk,su,tc,td,tksg,tu,uc,ud,uu,vc, &
              vcgd,vcgu,vcld,vclu,vxgd,vxgu,vxld, &
              vxlu,wc,xced,xedg,xedgd,xedgu,xedl,xedld,xedlu

      INTEGER :: i
      REAL, PARAMETER :: sml = 1.e-14
      REAL, PARAMETER :: smlc = 2.01e-14

      DO i = 1,irmx

         IF (jspins == 1) THEN
            rou=rh(i,1)/2
            rou=max(rou,sml)
            rod=rou
         ELSE
            rou=rh(i,1)
            rod=rh(i,jspins)
            rou=max(rou,sml)
            rod=max(rod,sml)
         ENDIF

         ro=rou+rod

         !.....
         !       xedl,xedg: exchange energy density (local,grad.exp.) in ry.
         !       cedl,cedg: exchange energy density (local,grad.expnd.) in ry.
         !.....
         xedl = 0.0e0
         cedl = 0.0e0
         xedg = 0.0e0
         cedg = 0.0e0

         if(ro <= smlc) go to 200

         zta=(rou-rod)/ro
         if(zta > 1.e0-sml) zta = 1.e0 - sml
         if(zta < -1.e0-sml) zta = -1.e0 + sml

         !.....
         c13 = 1.e0/3.e0
         c23 = 2.e0/3.e0
         c43 = 4.e0/3.e0
         c53 = 5.e0/3.e0
         !.....  alf=-3*(3/4*pai)**(1/3).
         alf = -1.861051473e0
         !.....
         ro2 = ro*ro
         ro13 = ro**c13
         !.....
         rou3 = rou**3
         rou43 = rou**c43
         !.....
         rod3 = rod**3
         rod43 = rod**c43
         !.....
         !       gr2=drr*drr
         !       gr2u=drru**2
         !       drrd=drr-drru
         !       gr2d=drrd**2
         !       ddrrd=ddrr-ddrru
         !.....
         !.....  gz,gz2,gz3: for wang-perdew ssf.
         gz = ((1.e0+zta)**c23+ (1.e0-zta)**c23)/2.e0
         !.....
         !.....
         rs = 0.620350491e0/ro13
         !.....
         !.....  xedl: exchange-energy-density in ry.
         xedl = alf* (rou43+rod43)
         !.....
         !....  .gradient correction.
         !.....
         !         if(abs(agr(i)).lt.sml) go to 200
         !.....

         !.....
         dsprs = 1.e0
         if(idsprs == 1) dsprs = 1.e-19
         if(isprsv == 1) dsprs = dsprs*sprsv

         !.....
         !      agr,agru,agrd: abs(grad(rho)), for all, up, and down.
         !c     gr2,gr2u,gr2d: grad(rho_all)**2, grad(rho_up)**2, grad(rho_d)**2.
         !      g2r,g2ru,g2rd: laplacian rho_all, _up and _down.
         !      gggru,-d: grad(rho)*grad(abs(grad(rho))) for all,up and down.
         !      grgru,-d: grad(rho_all)*grad(rhor_up) and for down.

         !         g2r=ddrr+2*drr/rv
         !.....
         rou53 = rou**c53
         !.....
         !.....  edrru: d(abs(d(rou)/dr))/dr, edrrd for down.
         !         edrru=ddrru
         !         if(drru.lt.0.) edrru=-ddrru
         !.....
         !         agr,agbru,-d: abs(grad(rho)),for rou, rod.
         !         gggru,-d: grad(rho)*grad(abs(grad(rho))) for up and down.
         !.....  su:at ro=2*rou. 1/(2(3*pai**2)**(1/3))*|grad(rou)|/rou**(4/3).
         su = 0.128278244e0*agru(i)/rou43
         !         if(su.gt.huges) go to 200
         !         g2ru=ddrru+2*drru/rv
         tu = .016455307e0*g2ru(i)/rou53
         uu = 0.002110857e0*gggru(i)/rou3

         dbrou = rou*2

         call xch91(dbrou,su,uu,tu,xedlu,xedgu,vxlu,vxgu)

         xedl = xedlu/2
         xedgu = dsprs*xedgu

         !.....
         rod53 = rod**c53
         !         edrrd=ddrrd
         !         if(drrd.lt.0.) edrrd=-ddrrd

         sd = 0.128278244e0*agrd(i)/rod43
         !         if(sd.gt.huges) go to 200

         !         g2rd=ddrrd+2*drrd/rv

         td = .016455307e0*g2rd(i)/rod53
         ud = 0.002110857e0*gggrd(i)/rod3

         dbrod = rod*2

         call xch91(dbrod,sd,ud,td,xedld,xedgd,vxld,vxgd)

         xedl = xedl + xedld/2
         xedgd = dsprs*xedgd

         xedg = (xedgu+xedgd)/2

         !....   cro: c(n) of (6),phys.rev..b33,8822('86). in ry.
         !....   dcdr: d(cro)/d(ro).
         !.....  0.001625816=1.745*f(=0.11)*cro(rs=0).

         call corl91(rs,zta,ec,vclu,vcld,ecrs,eczta,alfc)

         cedl = ec*2.e0*ro

         fk = 1.91915829e0/rs
         sk = sqrt(4.e0*fk/pi_const)
         tksg = 2.e0*sk*gz
         tc = agr(i)/ (ro*tksg)
         uc = gggr(i)/ (ro2*tksg**3)
         vc = g2r(i)/ (ro*tksg**2)
         wc = gzgr(i)/ (ro*tksg**2)

         call corg91(fk,sk,gz,ec,ecrs,eczta,rs,zta,tc,uc,vc,wc,cedg, &
                     vcgu,vcgd)

         cedg = dsprs*cedg*ro*2.e0

200      continue

         xced = (xedl+cedl+xedg+cedg)
         if(ro > smlc) xced = xced/ro

         !     in Ry.

         exc(i) = xced

      END DO

   END SUBROUTINE excpw91
END MODULE m_excpw91
