      MODULE m_excepbe
c.....-----------------------------------------------------------------
c.....epbe(easy_pbe) exchange-correlation energy density  in hartree.
c     excepbe - easypbe
c.....------------------------------------------------------------------
      CONTAINS
      SUBROUTINE excepbe(
     >                   xcpot,jspins,mirm,irmx,
     >                   rh,agr,agru,agrd,
     +                   g2ru,g2rd,gggr,gggru,gggrd,
     <                   exc)

      USE m_easypbe
      USE m_types_xcpot_data

      IMPLICIT NONE

! .. Arguments ..
      TYPE(t_xcpot_data),INTENT(IN)::xcpot
      INTEGER, INTENT (IN) :: irmx,jspins,mirm
      REAL,    INTENT (IN) :: rh(mirm,jspins)
      REAL,    INTENT (IN) :: agr(mirm),agru(mirm),agrd(mirm)
      REAL,    INTENT (IN) :: g2ru(mirm),g2rd(mirm),gggr(mirm)
      REAL,    INTENT (IN) :: gggru(mirm),gggrd(mirm)
      REAL,    INTENT (OUT) :: exc(mirm)


      ! .. local variables ..
      INTEGER lcor,lpot,i
      REAL ro,rou,rod,xedl,cedl,xedg,cedg,xced
      REAL up,agrup,delgrup,uplap,dn,agrdn,delgrdn,dnlap,
     +     agrt,delgrt,
     +     exlsd,vxuplsd,vxdnlsd,eclsd,vcuplsd,vcdnlsd,
     +     expbe,vxuppbe,vxdnpbe,ecpbe,vcuppbe,vcdnpbe,
     +     vxupsr,vxdnsr

      REAL, PARAMETER :: sml = 1.e-14
      REAL, PARAMETER :: smlc = 2.01e-14

!$OMP parallel do default(private)
!$OMP& SHARED(xcpot,jspins,mirm,irmx)
!$OMP& SHARED(rh,agr,agru,agrd)
!$OMP& SHARED(g2ru,g2rd,gggr,gggru,gggrd)
!$OMP& SHARED(exc)
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

c.....
c       xedl,xedg: exchange energy density (local,grad.exp.) in ry.
c       cedl,cedg: exchange energy density (local,grad.expnd.) in ry.
c.....
          xedl = 0.0e0
          cedl = 0.0e0
          xedg = 0.0e0
          cedg = 0.0e0

        IF (ro > smlc) THEN

          lcor=1
          lpot=1
          up=rou
          agrup=agru(i)
          delgrup=gggru(i)
          uplap=g2ru(i)
          dn=rod
          agrdn=agrd(i)
          delgrdn=gggrd(i)
          dnlap=g2rd(i)
          agrt=agr(i)
          delgrt=gggr(i)

          CALL easypbe (xcpot,
     &           up,agrup,delgrup,uplap,dn,agrdn,delgrdn,dnlap,
     1           agrt,delgrt,lcor,lpot,
     1           exlsd,vxuplsd,vxdnlsd,eclsd,vcuplsd,vcdnlsd,
     1           expbe,vxuppbe,vxdnpbe,ecpbe,vcuppbe,vcdnpbe,
     1           vxupsr,vxdnsr)

          xedl=exlsd
          cedl=eclsd
          xedg=expbe-exlsd
          cedg=ecpbe-eclsd

        ENDIF ! ro > smlc

        xced = (xedl+cedl+xedg+cedg)

        exc(i) = xced*2  ! in ry

      ENDDO
!$OMP end parallel do
      END SUBROUTINE excepbe
      END MODULE m_excepbe
