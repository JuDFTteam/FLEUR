      MODULE m_easypbe
c----------------------------------------------------------------------
c     easypbe---exchpbe
c             --corpbe  --- pbecor2
c----------------------------------------------------------------------
c easypbe is a driver for the pbe subroutines, using simple inputs
c k. burke, may 14, 1996.
c----------------------------------------------------------------------
      CONTAINS
      SUBROUTINE easypbe(icorr,
     >                   up,agrup,delgrup,uplap,
     >                   dn,agrdn,delgrdn,dnlap,
     >                   agr,delgr,lcor,lpot,
     <                   exlsd,vxuplsd,vxdnlsd,
     <                   eclsd,vcuplsd,vcdnlsd,
     <                   expbe,vxuppbe,vxdnpbe,
     <                   ecpbe,vcuppbe,vcdnpbe,
     <                   vxupsr,vxdnsr)

      USE m_exchpbe
      USE m_corpbe
      IMPLICIT NONE

      ! .. Arguments ..
      INTEGER, INTENT (IN) :: icorr    ! selects different X-enhancement in exchpbe
      INTEGER, INTENT (IN) :: lcor     ! flag to do correlation(=0=>don't)
      INTEGER, INTENT (IN) :: lpot     ! flag to do potential  (=0=>don't)
      REAL,    INTENT (IN) :: up,dn    ! density (spin up & down)
      REAL,    INTENT (IN) :: agrup,agrdn     ! |grad up|, |grad down|
      REAL,    INTENT (IN) :: delgrup,delgrdn ! delgrup=(grad up).(grad |grad up|)  (in pw91: gggru)
      REAL,    INTENT (IN) :: uplap,dnlap     ! grad^2 up=laplacian of up           (in pw91: g2ru)
      REAL,    INTENT (IN) :: agr,delgr       ! agr=|grad rho|, 
                                              ! delgr=(grad rho).(grad |grad rho|)  (in pw91: gggr)
      REAL,    INTENT (OUT) :: vxuplsd,vxdnlsd ! up/down lsd exchange potential
      REAL,    INTENT (OUT) :: vcuplsd,vcdnlsd ! up/down lsd correlation potential
      REAL,    INTENT (OUT) :: exlsd,eclsd     ! lsd exchange / correlation energy density
      REAL,    INTENT (OUT) :: vxuppbe,vxdnpbe ! as above, but pbe quantities
      REAL,    INTENT (OUT) :: vcuppbe,vcdnpbe
      REAL,    INTENT (OUT) :: expbe,ecpbe     ! note that : exlsd=int d^3r rho(r) exlsd(r)
      REAL,    INTENT (OUT) :: vxupsr,vxdnsr

      ! .. local variables ..
      ! for exchange:
      REAL fk ! local fermi wavevector for 2*up=(3 pi^2 (2up))^(1/3)
      REAL s  ! dimensionless density gradient=|grad rho|/ (2*fk*rho)_(rho=2*up)
      REAL u  ! delgrad/(rho^2*(2*fk)**3)_(rho=2*up)
      REAL v  ! laplacian/(rho*(2*fk)**2)_(rho=2*up)
      ! for correlation:
      REAL zet    ! (up-dn)/rho
      REAL g      ! phi(zeta)
      REAL rs     ! (3/(4pi*rho))^(1/3)=local seitz radius=alpha/fk
      REAL sk     ! ks=thomas-fermi screening wavevector=sqrt(4fk/pi)
      REAL twoksg ! 2*ks*phi
      REAL t      ! correlation dimensionless gradient=|grad rho|/(2*ks*phi*rho)
      REAL uu     ! delgrad/(rho^2*twoksg^3)
      REAL rholap ! laplacian
      REAL vv     ! laplacian/(rho*twoksg^2)
      REAL ww     ! (|grad up|^2-|grad dn|^2-zet*|grad rho|^2)/(rho*twoksg)^2
      REAL ec     ! lsd correlation energy
      REAL vcup   ! lsd up correlation potential
      REAL vcdn   ! lsd down correlation potential
      REAL h      ! gradient correction to correlation energy
      REAL dvcup  ! gradient correction to up correlation potential
      REAL dvcdn  ! gradient correction to down correlation potential

      REAL rdum,eta,exdnlsd,exdnpbe,exuplsd,exuppbe,rho,rho2
c     ..
      REAL, PARAMETER :: thrd = 1.e0/3.e0
      REAL, PARAMETER :: thrd2 = 2.e0*thrd
      REAL, PARAMETER :: pi = 3.1415926535897932384626433832795e0
      REAL, PARAMETER :: pi32 = 29.608813203268075856503472999628e0 ! 3 pi**2
      REAL, PARAMETER :: alpha=1.91915829267751300662482032624669e0 ! (9pi/4)**thrd

      rho2 = 2.e0*up

c----------------------------------------------------------------------
c pbe exchange
c use  ex[up,dn]=0.5*(ex[2*up]+ex[2*dn]) (i.e., exact spin-scaling)
c do up exchange
c----------------------------------------------------------------------

      IF (rho2.GT.1e-18) THEN
          fk = (pi32*rho2)**thrd
          s = 2.e0*agrup/ (2.e0*fk*rho2)
          u = 4.e0*delgrup/ (rho2*rho2* (2.e0*fk)**3)
          v = 2.e0*uplap/ (rho2* (2.e0*fk)**2)

          CALL exchpbe(icorr,rho2,s,u,v,0,lpot,exuplsd,vxuplsd,rdum)
          CALL exchpbe(icorr,rho2,s,u,v,1,lpot,exuppbe,vxuppbe,vxupsr)

      ELSE

          exuplsd = 0.e0
          vxuplsd = 0.e0
          exuppbe = 0.e0
          vxuppbe = 0.e0
          vxupsr  = 0.e0

      ENDIF

c repeat for down
      rho2 = 2.e0*dn

      IF (rho2.gt.1e-18) THEN

          fk = (pi32*rho2)**thrd
          s = 2.e0*agrdn/ (2.e0*fk*rho2)
          u = 4.e0*delgrdn/ (rho2*rho2* (2.e0*fk)**3)
          v = 2.e0*dnlap/ (rho2* (2.e0*fk)**2)

          CALL exchpbe(icorr,rho2,s,u,v,0,lpot,exdnlsd,vxdnlsd,rdum)
          CALL exchpbe(icorr,rho2,s,u,v,1,lpot,exdnpbe,vxdnpbe,vxdnsr)

      ELSE

          exdnlsd = 0.e0
          vxdnlsd = 0.e0
          exdnpbe = 0.e0
          vxdnpbe = 0.e0
          vxdnsr  = 0.e0

      ENDIF

c construct total density and contribution to ex
      rho = up + dn
      exlsd = (exuplsd*up+exdnlsd*dn)/rho
      expbe = (exuppbe*up+exdnpbe*dn)/rho

      IF (lcor.EQ.0) RETURN
c----------------------------------------------------------------------
c now do correlation
c----------------------------------------------------------------------

      IF (rho.lt.1.e-18) RETURN

      zet = (up-dn)/rho

c9999+
c     eta: eta should not be smaller than 1.e-16.
c          otherwise will cause floating invalid.
c          if bigger, the last digit may differ
c          from the run by aix without this zet-guard.
      eta = 1.e-16
      zet = min(zet,1.0-eta)
      zet = max(zet,-1.0+eta)
c9999-

      g = ((1.e0+zet)**thrd2+ (1.e0-zet)**thrd2)/2.e0
      fk = (pi32*rho)**thrd
      rs = alpha/fk
      sk = sqrt(4.e0*fk/pi)
      twoksg = 2.e0*sk*g
      t = agr/ (twoksg*rho)
      uu = delgr/ (rho*rho*twoksg**3)
      rholap = uplap + dnlap
      vv = rholap/ (rho*twoksg**2)
      ww = (agrup**2-agrdn**2-zet*agr**2)/ (rho*rho*twoksg**2)

      CALL corpbe(
     >            icorr,rs,zet,t,uu,vv,ww,1,lpot,
     <            ec,vcup,vcdn,h,dvcup,dvcdn)

      eclsd = ec
      ecpbe = ec + h
      vcuplsd = vcup
      vcdnlsd = vcdn
      vcuppbe = vcup + dvcup
      vcdnpbe = vcdn + dvcdn

      END SUBROUTINE easypbe
      END MODULE m_easypbe
