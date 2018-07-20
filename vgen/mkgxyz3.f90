!--------------------------------------------------------------------------------
! Copyright (c) 2016 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------
MODULE m_mkgxyz3
  !.....------------------------------------------------------------------
  !c     by use of cartesian x,y,z components of charge density gradients,
  !c     make the quantities
  !cc      agrt,agru,agrd,g2rt,g2ru,g2rd,gggrt,gggru,gggrd,gzgr
  !cc    used to calculate gradient contribution to xc potential and
  !cc    energy.
  !c.....------------------------------------------------------------------
CONTAINS
  SUBROUTINE mkgxyz3(vl,dvx,dvy,dvz,dvxx,dvyy,dvzz,dvyz,dvzx,dvxy,grad)
    USE m_types
    IMPLICIT NONE
    REAL,    INTENT (IN) :: vl(:,:)
    REAL,    INTENT (IN) :: dvx(:,:),dvy(:,:),dvz(:,:)
    REAL, INTENT (IN) :: dvxx(:,:),dvyy(:,:),dvzz(:,:)
    REAL, INTENT (IN) :: dvyz(:,:),dvzx(:,:),dvxy(:,:)

    TYPE(t_gradients),INTENT(INOUT)::grad

    REAL vlt,dvxt,dvyt,dvzt,dvxxt,dvyyt,dvzzt,dvyzt,dvzxt,dvxyt,&
         vlu,dvxu,dvyu,dvzu,dvxxu,dvyyu,dvzzu,dvyzu,dvzxu,dvxyu,&
         vld,dvxd,dvyd,dvzd,dvxxd,dvyyd,dvzzd,dvyzd,dvzxd,dvxyd,&
         dagrxt,dagrxd,dagrxu,dagryt,dagryd,dagryu,dagrzt,dagrzd,&
         dagrzu,dzdx,dzdy,dzdz,sml
    INTEGER i,js,jspins,nsp

    nsp=SIZE(dvx,1)
    jspins=SIZE(dvx,2)
    sml = 1.e-14

    IF (ALLOCATED(grad%gr)) THEN
       !      Gradients for libxc
       DO js=1,jspins
          DO i=1,nsp
             grad%gr(:,i,js)=(/dvx(i,js),dvy(i,js),dvz(i,js)/)
          ENDDO
       END DO
       IF(ALLOCATED(grad%sigma)) THEN
          !Use only contracted gradients for libxc
          IF (jspins==1) THEN
             DO i=1,nsp
                grad%sigma(1,i) = dvx(i,1)*dvx(i,1)+dvy(i,1)*dvy(i,1)+dvz(i,1)*dvz(i,1)
             ENDDO
          ELSE
             DO i=1,nsp
                grad%sigma(1,i) = dvx(i,1)*dvx(i,1) + dvy(i,1)*dvy(i,1) + dvz(i,1)*dvz(i,1)
                grad%sigma(2,i) = dvx(i,1)*dvx(i,2) + dvy(i,1)*dvy(i,2) + dvz(i,1)*dvz(i,2)
                grad%sigma(3,i) = dvx(i,2)*dvx(i,2) + dvy(i,2)*dvy(i,2) + dvz(i,2)*dvz(i,2)
             ENDDO
          ENDIF
       END IF
       IF(ALLOCATED(grad%laplace)) THEN
          DO js=1,jspins
             DO i=1,nsp
                grad%laplace(i,js)= dvxx(i,js)+dvyy(i,js)+dvzz(i,js)
             ENDDO
          ENDDO
       ENDIF
       RETURN
    ENDIF

    IF (ANY(SHAPE(vl).NE.SHAPE(dvx))) CALL judft_error("Gradients for internal GGA called with inconsistent sizes",hint="This is a bug")
    
    DO i = 1,size(grad%agrt)
       grad%agrt(i)  = 0.0
       grad%agru(i)  = 0.0
       grad%agrd(i)  = 0.0
       grad%gggrt(i) = 0.0
       grad%gggru(i) = 0.0
       grad%gggrd(i) = 0.0
       grad%gzgr(i)  = 0.0
       grad%g2rt(i)  = 0.0
       grad%g2ru(i)  = 0.0
       grad%g2rd(i)  = 0.0
    ENDDO

    IF (jspins.eq.1) THEN

       DO i = 1,nsp

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

          !         agr: abs(grad(ro)), t,u,d for total, up and down.

          grad%agrt(i) = max(sqrt(dvxt**2+dvyt**2+dvzt**2),sml)
          grad%agru(i) = max(sqrt(dvxu**2+dvyu**2+dvzu**2),sml)
          grad%agrd(i) = max(sqrt(dvxd**2+dvyd**2+dvzd**2),sml)

          dagrxt = (dvxt*dvxxt+dvyt*dvxyt+dvzt*dvzxt)/grad%agrt(i)
          dagrxu = (dvxu*dvxxu+dvyu*dvxyu+dvzu*dvzxu)/grad%agru(i)
          dagrxd = (dvxd*dvxxd+dvyd*dvxyd+dvzd*dvzxd)/grad%agrd(i)

          dagryt = (dvxt*dvxyt+dvyt*dvyyt+dvzt*dvyzt)/grad%agrt(i)
          dagryu = (dvxu*dvxyu+dvyu*dvyyu+dvzu*dvyzu)/grad%agru(i)
          dagryd = (dvxd*dvxyd+dvyd*dvyyd+dvzd*dvyzd)/grad%agrd(i)

          dagrzt = (dvxt*dvzxt+dvyt*dvyzt+dvzt*dvzzt)/grad%agrt(i)
          dagrzu = (dvxu*dvzxu+dvyu*dvyzu+dvzu*dvzzu)/grad%agru(i)
          dagrzd = (dvxd*dvzxd+dvyd*dvyzd+dvzd*dvzzd)/grad%agrd(i)

          grad%gggrt(i) = dvxt*dagrxt + dvyt*dagryt + dvzt*dagrzt
          grad%gggru(i) = dvxu*dagrxu + dvyu*dagryu + dvzu*dagrzu
          grad%gggrd(i) = dvxd*dagrxd + dvyd*dagryd + dvzd*dagrzd

          !         dzdx=d(zeta)/dx,..

          dzdx = (dvxu-dvxd)/vlt - (vlu-vld)*dvxt/vlt**2
          dzdy = (dvyu-dvyd)/vlt - (vlu-vld)*dvyt/vlt**2
          dzdz = (dvzu-dvzd)/vlt - (vlu-vld)*dvzt/vlt**2

          !         gzgr=grad(zeta)*grad(ro).

          grad%gzgr(i) = dzdx*dvxt + dzdy*dvyt + dzdz*dvzt

          !         g2r: grad(grad(ro))

          grad%g2rt(i) = dvxxt + dvyyt + dvzzt
          grad%g2ru(i) = dvxxu + dvyyu + dvzzu
          grad%g2rd(i) = dvxxd + dvyyd + dvzzd


       ENDDO

    ELSE

       DO i = 1,nsp

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

          !c         agr: abs(grad(ro)), t,u,d for total, up and down.

          grad%agrt(i) = max(sqrt(dvxt**2+dvyt**2+dvzt**2),sml)
          grad%agru(i) = max(sqrt(dvxu**2+dvyu**2+dvzu**2),sml)
          grad%agrd(i) = max(sqrt(dvxd**2+dvyd**2+dvzd**2),sml)

          dagrxt = (dvxt*dvxxt+dvyt*dvxyt+dvzt*dvzxt)/grad%agrt(i)
          dagrxu = (dvxu*dvxxu+dvyu*dvxyu+dvzu*dvzxu)/grad%agru(i)
          dagrxd = (dvxd*dvxxd+dvyd*dvxyd+dvzd*dvzxd)/grad%agrd(i)

          dagryt = (dvxt*dvxyt+dvyt*dvyyt+dvzt*dvyzt)/grad%agrt(i)
          dagryu = (dvxu*dvxyu+dvyu*dvyyu+dvzu*dvyzu)/grad%agru(i)
          dagryd = (dvxd*dvxyd+dvyd*dvyyd+dvzd*dvyzd)/grad%agrd(i)

          dagrzt = (dvxt*dvzxt+dvyt*dvyzt+dvzt*dvzzt)/grad%agrt(i)
          dagrzu = (dvxu*dvzxu+dvyu*dvyzu+dvzu*dvzzu)/grad%agru(i)
          dagrzd = (dvxd*dvzxd+dvyd*dvyzd+dvzd*dvzzd)/grad%agrd(i)

          grad%gggrt(i) = dvxt*dagrxt + dvyt*dagryt + dvzt*dagrzt
          grad%gggru(i) = dvxu*dagrxu + dvyu*dagryu + dvzu*dagrzu
          grad%gggrd(i) = dvxd*dagrxd + dvyd*dagryd + dvzd*dagrzd

          !c         dzdx=d(zeta)/dx,..

          dzdx = (dvxu-dvxd)/vlt -  (vlu-vld)*dvxt/vlt**2
          dzdy = (dvyu-dvyd)/vlt -  (vlu-vld)*dvyt/vlt**2
          dzdz = (dvzu-dvzd)/vlt -  (vlu-vld)*dvzt/vlt**2

          !c         gzgr=grad(zeta)*grad(ro).

          grad%gzgr(i) = dzdx*dvxt + dzdy*dvyt + dzdz*dvzt

          !c         g2r: grad(grad(ro))

          grad%g2rt(i) = dvxxt + dvyyt + dvzzt
          grad%g2ru(i) = dvxxu + dvyyu + dvzzu
          grad%g2rd(i) = dvxxd + dvyyd + dvzzd

       ENDDO

    ENDIF

  END SUBROUTINE mkgxyz3
END MODULE m_mkgxyz3
