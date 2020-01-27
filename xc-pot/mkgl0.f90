MODULE m_mkgl0
!      ------------------------------------------------------------------
!      make quantities for vxcallg. for paramag. case
!        dens,drr,ddrr are charge density and its gradient for nonmag.
!        one spin.

!      agr: abs(grad(ro)), g2r: laplacian(ro),
!      gggr: grad(ro)*grad(agr),
!      grgru,d: grad(ro)*grad(rou),for rod., gzgr: grad(zeta)*grad(ro).
!      ------------------------------------------------------------------
CONTAINS
   SUBROUTINE mkgl0(jspins,rad, densi,drri,ddrri, grad)

      USE m_types
      IMPLICIT NONE

      INTEGER, INTENT (IN) :: jspins
      REAL,    INTENT (IN) :: rad(:),densi(:,:)
      REAL,    INTENT (IN) :: drri(:,:),ddrri(:,:)
      TYPE(t_gradients)::grad

      REAL dagrr,dagrru,ddrr,ddrrd,ddrru,drr,drrd,drru,dzdr,ro,rod,rou,rv,spnf
      INTEGER i
!     ------------------------------------------------------------------

      spnf = 1./(3-jspins)
      IF (allocated(grad%sigma)) THEN
         DO i=1,size(rad)
            grad%sigma(1,i)=spnf * drri(i,1)*spnf * drri(i,1)
            IF (jspins>1) THEN
               grad%sigma(2,i)=spnf * drri(i,1)*spnf * drri(i,2)
               grad%sigma(3,i)=spnf * drri(i,2)*spnf * drri(i,2)
            ENDIF
         ENDDO
         RETURN
      ENDIF
      DO i = 1,size(rad)

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

         grad%agrt(i) = abs(drr)
         grad%agru(i) = abs(drru)
         grad%agrd(i) = grad%agru(i)

         dagrr = drr*ddrr/grad%agrt(i)
         dagrru = drru*ddrru/grad%agru(i)

         grad%gggrt(i) = drr*dagrr
         grad%gggru(i) = drru*dagrru
         grad%gggrd(i) = grad%gggru(i)

         dzdr = ((drru-drrd)*ro- (rou-rod)*drr)/ro**2

         grad%gzgr(i) = dzdr*drr

         grad%g2rt(i) = ddrr + 2*drr/rv
         grad%g2ru(i) = ddrru + 2*drru/rv
         grad%g2rd(i) = ddrrd + 2*drrd/rv

         grad%grgru(i) = drr*drru
         grad%grgrd(i) = drr*drrd

      ENDDO

   END SUBROUTINE mkgl0
END MODULE m_mkgl0
