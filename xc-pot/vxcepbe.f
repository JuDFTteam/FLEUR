      MODULE m_vxcepbe
c.....-----------------------------------------------------------------
c.....epbe(easy_pbe) exchange-correlation potential in hartree.
c     vxcepbe - easypbe 
c.....------------------------------------------------------------------
      CONTAINS
      SUBROUTINE vxcepbe(
     >                   icorr,jspins,mirm,irmx,
     >                   rh,agr,agru,agrd,
     +                   g2ru,g2rd,gggr,gggru,gggrd,
     <                   vx,vxc)
      
      Use m_easypbe
      USE m_icorrkeys

      IMPLICIT NONE

      ! .. Arguments ..
      INTEGER, INTENT (IN) :: icorr,irmx,jspins,mirm
      REAL,    INTENT (IN) :: rh(mirm,jspins)
      REAL,    INTENT (IN) :: agr(mirm),agru(mirm),agrd(mirm)
      REAL,    INTENT (IN) :: g2ru(mirm),g2rd(mirm),gggr(mirm)
      REAL,    INTENT (IN) :: gggru(mirm),gggrd(mirm)
      REAL,    INTENT (OUT):: vx(mirm,jspins),vxc(mirm,jspins)

      ! .. local variables ..
      INTEGER lcor,lpot,i
      REAL ro,rou,rod,xcptu,xcptd,
     +     vxlu,vclu,vxld,vcld,vxgu,vcgu,vxgd,vcgd
      REAL up,agrup,delgrup,uplap,dn,agrdn,delgrdn,dnlap,
     +     agrt,delgrt,
     +     exlsd,vxuplsd,vxdnlsd,eclsd,vcuplsd,vcdnlsd,
     +     expbe,vxuppbe,vxdnpbe,ecpbe,vcuppbe,vcdnpbe,
     +     vxupsr,vxdnsr

      REAL, PARAMETER :: sml = 1.e-14
      REAL, PARAMETER :: smlc = 2.01e-14

      DO i = 1,irmx

        IF (jspins.eq.1) THEN
          rou=rh(i,1)/2
          rou=max(rou,sml)
          rod=rou
        ELSE
          rou=rh(i,1)
          rod=rh(i,jspins)
          rou=max(rou,sml)
          rod=max(rod,sml)
        ENDIF

c.....
c       vxlu,vxld,vxgu,vxgd: exchange potential in ry.(local,grad),
cc        (up,dw).
c       vclu,vcld,vcgu,vcgd: correl. potential in ry.(local,grad),
cc        (up,dw).
c       all later in hartree.
c.....
        vxlu   = 0.0e0
        vclu   = 0.0e0
        vxld   = 0.0e0
        vcld   = 0.0e0
        vxgu   = 0.0e0
        vcgu   = 0.0e0
        vxgd   = 0.0e0
        vcgd   = 0.0e0
        vxupsr = 0.0e0
        vxdnsr = 0.0e0

c.....
        ro=rou+rod

        IF (ro > smlc) THEN

          lcor    = 1
          lpot    = 1
          up      = rou
          agrup   = agru(i)
          delgrup = gggru(i)
          uplap   = g2ru(i)
          dn      = rod
          agrdn   = agrd(i)
          delgrdn = gggrd(i)
          dnlap   = g2rd(i)
          agrt    = agr(i)
          delgrt  = gggr(i)

          CALL easypbe(icorr,
     &           up,agrup,delgrup,uplap,dn,agrdn,delgrdn,dnlap,
     1           agrt,delgrt,lcor,lpot,
     1           exlsd,vxuplsd,vxdnlsd,eclsd,vcuplsd,vcdnlsd,
     1           expbe,vxuppbe,vxdnpbe,ecpbe,vcuppbe,vcdnpbe,
     1           vxupsr,vxdnsr)

          vxlu=vxuplsd
          vclu=vcuplsd
          vxgu=vxuppbe-vxuplsd
          vcgu=vcuppbe-vcuplsd

          vxld=vxdnlsd
          vcld=vcdnlsd
          vxgd=vxdnpbe-vxdnlsd
          vcgd=vcdnpbe-vcdnlsd

        END IF ! ro > smlc

        xcptu = vxlu + vclu + vxgu + vcgu
        xcptd = vxld + vcld + vxgd + vcgd

        IF ( icorr.EQ.icorr_hse .OR. icorr.EQ.icorr_hseloc .OR.
     +       icorr.EQ.icorr_vhse ) THEN
          vx(i,1)       = vxupsr * 2
          vx(i,jspins)  = vxdnsr * 2
        ELSE
          vx(i,1)       = (vxlu + vxgu)*2
          vx(i,jspins)  = (vxld + vxgd)*2
        END IF

        vxc(i,1)      = xcptu*2    ! transform to Ry, will be converted
        vxc(i,jspins) = xcptd*2    ! back to htr in calling routine

      END DO

      END SUBROUTINE vxcepbe
      END MODULE m_vxcepbe
