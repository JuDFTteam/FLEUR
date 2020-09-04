!--------------------------------------------------------------------------------
! Copyright (c) 2016 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------
MODULE m_mkgxyz3
   USE m_judft
   !-----------------------------------------------------------------------------
   ! Using the cartesian components and derivatives of a charge density rho on
   ! the real space grid, make the following quantaties:
   !
   ! gr(js):      grad(rho_js)
   ! sigma:       |grad(rho)|^2 for jspins==1, otherwise three components, namely
   !              |grad(rho_up)|^2,grad(rho_up)*grad(rho_down),|grad(rho_down)|^2
   ! laplace(js): laplace(rho_js)
   !
   ! and these older components:
   ! agrt/u/d:    |grad(rho)| for total density/spin-up/spin-down
   ! g2rt/u/d:    laplace(rho)
   ! gggrt/u/d:   (grad(rho))*(grad(|grad(rho)|)) [scalar product]
   ! gzgr:        (grad(zeta))*(grad(rho)) for zeta=(rho_up-rho_down)/rho
   !
   ! which are used to calculate gradient contributions to the xc potential and
   ! energy.
   !
   ! vl is rho, dv[i][j] are the partial derivatives along one/two directions.
   !
   ! Modified so only allocated old quantities require calculations. A.N. 2019
   !
   ! Quantities fo libxc are calculated as well. D.W./M.R. 2018
   !
   ! Original script by T.A. 1996
   !-----------------------------------------------------------------------------
CONTAINS
   SUBROUTINE mkgxyz3(vl,dvx,dvy,dvz,dvxx,dvyy,dvzz,dvyz,dvxz,dvxy,idx,grad)
      USE m_types
      IMPLICIT NONE
      REAL, INTENT (IN)                :: vl(:,:)
      REAL, INTENT (IN)                :: dvx(:,:),dvy(:,:),dvz(:,:)
      REAL, INTENT (IN)                :: dvxx(:,:),dvyy(:,:),dvzz(:,:)
      REAL, INTENT (IN)                :: dvyz(:,:),dvxz(:,:),dvxy(:,:)
      INTEGER ,intent(in)              :: idx
      TYPE(t_gradients), INTENT(INOUT) :: grad

      REAL vlt,dvxt,dvyt,dvzt,dvxxt,dvyyt,dvzzt,dvyzt,dvxzt,dvxyt, &
           vlu,dvxu,dvyu,dvzu,dvxxu,dvyyu,dvzzu,dvyzu,dvxzu,dvxyu, &
           vld,dvxd,dvyd,dvzd,dvxxd,dvyyd,dvzzd,dvyzd,dvxzd,dvxyd, &
           dagrxt,dagrxd,dagrxu,dagryt,dagryd,dagryu,dagrzt,dagrzd,&
           dagrzu,dzdx,dzdy,dzdz
      REAL, PARAMETER  :: sml = 1.e-14
      INTEGER i,js,jspins,nsp

      nsp=SIZE(dvx,1)
      jspins=SIZE(dvx,2)

      ! Gradients for libxc and sourcefree calculations.
      IF (ALLOCATED(grad%gr)) THEN
         DO js=1,jspins
            DO i=1,nsp
               grad%gr(:,i+idx,js)=(/dvx(i,js),dvy(i,js),dvz(i,js)/)
            END DO
         END DO
         ! Use contracted gradients only for libxc.
         IF(ALLOCATED(grad%sigma)) THEN
            IF (jspins==1) THEN
               DO i=1,nsp
                  grad%sigma(1,i+idx) = dvx(i,1)*dvx(i,1) + dvy(i,1)*dvy(i,1) + dvz(i,1)*dvz(i,1)
               END DO
            ELSE
               DO i=1,nsp
                  grad%sigma(1,i+idx) = dvx(i,1)*dvx(i,1) + dvy(i,1)*dvy(i,1) + dvz(i,1)*dvz(i,1)
                  grad%sigma(2,i+idx) = dvx(i,1)*dvx(i,2) + dvy(i,1)*dvy(i,2) + dvz(i,1)*dvz(i,2)
                  grad%sigma(3,i+idx) = dvx(i,2)*dvx(i,2) + dvy(i,2)*dvy(i,2) + dvz(i,2)*dvz(i,2)
               END DO
            END IF
         END IF
         IF(ALLOCATED(grad%laplace)) THEN
            DO js=1,jspins
               DO i=1,nsp
                  grad%laplace(i+idx,js)= dvxx(i,js)+dvyy(i,js)+dvzz(i,js)
               END DO
            END DO
         END IF
         RETURN ! Do not calculate arrays for in-build GGA.
      END IF

      IF (ANY(SHAPE(vl).NE.SHAPE(dvx))) CALL judft_error("Gradients for internal GGA called with inconsistent sizes",hint="This is a bug")

 	   IF(ALLOCATED(grad%agrt)) THEN
         DO i = 1,nsp
            grad%agrt(idx+i)  = 0.0
            grad%agru(idx+i)  = 0.0
            grad%agrd(idx+i)  = 0.0
            grad%gggrt(idx+i) = 0.0
            grad%gggru(idx+i) = 0.0
            grad%gggrd(idx+i) = 0.0
            grad%gzgr(idx+i)  = 0.0
            grad%g2rt(idx+i)  = 0.0
            grad%g2ru(idx+i)  = 0.0
            grad%g2rd(idx+i)  = 0.0
         END DO

         IF (jspins.eq.1) THEN

            DO i = 1,nsp

               vlu   = max(vl(i,1)/2,sml)
               dvxu  = dvx(i,1)/2
               dvyu  = dvy(i,1)/2
               dvzu  = dvz(i,1)/2
               dvxxu = dvxx(i,1)/2
               dvyyu = dvyy(i,1)/2
               dvzzu = dvzz(i,1)/2
               dvyzu = dvyz(i,1)/2
               dvxzu = dvxz(i,1)/2
               dvxyu = dvxy(i,1)/2

               vld   = vlu
               dvxd  = dvxu
               dvyd  = dvyu
               dvzd  = dvzu
               dvxxd = dvxxu
               dvyyd = dvyyu
               dvzzd = dvzzu
               dvyzd = dvyzu
               dvxzd = dvxzu
               dvxyd = dvxyu

               vlt = vlu + vld
               dvxt  = dvxu  + dvxd
               dvyt  = dvyu  + dvyd
               dvzt  = dvzu  + dvzd
               dvxxt = dvxxu + dvxxd
               dvyyt = dvyyu + dvyyd
               dvzzt = dvzzu + dvzzd
               dvyzt = dvyzu + dvyzd
               dvxzt = dvxzu + dvxzd
               dvxyt = dvxyu + dvxyd

               grad%agrt(idx+i) = max(sqrt(dvxt**2 + dvyt**2 + dvzt**2),sml)
               grad%agru(idx+i) = max(sqrt(dvxu**2 + dvyu**2 + dvzu**2),sml)
               grad%agrd(idx+i) = max(sqrt(dvxd**2 + dvyd**2 + dvzd**2),sml)

               dagrxt = (dvxt*dvxxt + dvyt*dvxyt + dvzt*dvxzt) / grad%agrt(idx+i)
               dagrxu = (dvxu*dvxxu + dvyu*dvxyu + dvzu*dvxzu) / grad%agru(idx+i)
               dagrxd = (dvxd*dvxxd + dvyd*dvxyd + dvzd*dvxzd) / grad%agrd(idx+i)

               dagryt = (dvxt*dvxyt + dvyt*dvyyt + dvzt*dvyzt) / grad%agrt(idx+i)
               dagryu = (dvxu*dvxyu + dvyu*dvyyu + dvzu*dvyzu) / grad%agru(idx+i)
               dagryd = (dvxd*dvxyd + dvyd*dvyyd + dvzd*dvyzd) / grad%agrd(idx+i)

               dagrzt = (dvxt*dvxzt + dvyt*dvyzt + dvzt*dvzzt) / grad%agrt(idx+i)
               dagrzu = (dvxu*dvxzu + dvyu*dvyzu + dvzu*dvzzu) / grad%agru(idx+i)
               dagrzd = (dvxd*dvxzd + dvyd*dvyzd + dvzd*dvzzd) / grad%agrd(idx+i)

               grad%gggrt(idx+i) = dvxt*dagrxt + dvyt*dagryt + dvzt*dagrzt
               grad%gggru(idx+i) = dvxu*dagrxu + dvyu*dagryu + dvzu*dagrzu
               grad%gggrd(idx+i) = dvxd*dagrxd + dvyd*dagryd + dvzd*dagrzd

               dzdx = (dvxu-dvxd)/vlt - (vlu-vld)*dvxt/vlt**2
               dzdy = (dvyu-dvyd)/vlt - (vlu-vld)*dvyt/vlt**2
               dzdz = (dvzu-dvzd)/vlt - (vlu-vld)*dvzt/vlt**2

               grad%gzgr(idx+i) = dzdx*dvxt + dzdy*dvyt + dzdz*dvzt

               grad%g2rt(idx+i) = dvxxt + dvyyt + dvzzt
               grad%g2ru(idx+i) = dvxxu + dvyyu + dvzzu
               grad%g2rd(idx+i) = dvxxd + dvyyd + dvzzd
            ENDDO
         ELSE
            DO i = 1,nsp

               vlu   = max(vl(i,1),sml)
               dvxu  = dvx(i,1)
               dvyu  = dvy(i,1)
               dvzu  = dvz(i,1)
               dvxxu = dvxx(i,1)
               dvyyu = dvyy(i,1)
               dvzzu = dvzz(i,1)
               dvyzu = dvyz(i,1)
               dvxzu = dvxz(i,1)
               dvxyu = dvxy(i,1)

               vld   = max(vl(i,jspins),sml)
               dvxd  = dvx(i,jspins)
               dvyd  = dvy(i,jspins)
               dvzd  = dvz(i,jspins)
               dvxxd = dvxx(i,jspins)
               dvyyd = dvyy(i,jspins)
               dvzzd = dvzz(i,jspins)
               dvyzd = dvyz(i,jspins)
               dvxzd = dvxz(i,jspins)
               dvxyd = dvxy(i,jspins)

               vlt = vlu + vld

               dvxt  = dvxu  + dvxd
               dvyt  = dvyu  + dvyd
               dvzt  = dvzu  + dvzd
               dvxxt = dvxxu + dvxxd
               dvyyt = dvyyu + dvyyd
               dvzzt = dvzzu + dvzzd
               dvyzt = dvyzu + dvyzd
               dvxzt = dvxzu + dvxzd
               dvxyt = dvxyu + dvxyd

               grad%agrt(idx+i) = max(sqrt(dvxt**2 + dvyt**2 + dvzt**2),sml)
               grad%agru(idx+i) = max(sqrt(dvxu**2 + dvyu**2 + dvzu**2),sml)
               grad%agrd(idx+i) = max(sqrt(dvxd**2 + dvyd**2 + dvzd**2),sml)

               dagrxt = (dvxt*dvxxt + dvyt*dvxyt + dvzt*dvxzt) / grad%agrt(idx+i)
               dagrxu = (dvxu*dvxxu + dvyu*dvxyu + dvzu*dvxzu) / grad%agru(idx+i)
               dagrxd = (dvxd*dvxxd + dvyd*dvxyd + dvzd*dvxzd) / grad%agrd(idx+i)

               dagryt = (dvxt*dvxyt + dvyt*dvyyt + dvzt*dvyzt) / grad%agrt(idx+i)
               dagryu = (dvxu*dvxyu + dvyu*dvyyu + dvzu*dvyzu) / grad%agru(idx+i)
               dagryd = (dvxd*dvxyd + dvyd*dvyyd + dvzd*dvyzd) / grad%agrd(idx+i)

               dagrzt = (dvxt*dvxzt + dvyt*dvyzt + dvzt*dvzzt) / grad%agrt(idx+i)
               dagrzu = (dvxu*dvxzu + dvyu*dvyzu + dvzu*dvzzu) / grad%agru(idx+i)
               dagrzd = (dvxd*dvxzd + dvyd*dvyzd + dvzd*dvzzd) / grad%agrd(idx+i)

               grad%gggrt(idx+i) = dvxt*dagrxt + dvyt*dagryt + dvzt*dagrzt
               grad%gggru(idx+i) = dvxu*dagrxu + dvyu*dagryu + dvzu*dagrzu
               grad%gggrd(idx+i) = dvxd*dagrxd + dvyd*dagryd + dvzd*dagrzd

               dzdx = (dvxu-dvxd)/vlt -  (vlu-vld)*dvxt/vlt**2
               dzdy = (dvyu-dvyd)/vlt -  (vlu-vld)*dvyt/vlt**2
               dzdz = (dvzu-dvzd)/vlt -  (vlu-vld)*dvzt/vlt**2

               grad%gzgr(idx+i) = dzdx*dvxt + dzdy*dvyt + dzdz*dvzt

               grad%g2rt(idx+i) = dvxxt + dvyyt + dvzzt
               grad%g2ru(idx+i) = dvxxu + dvyyu + dvzzu
               grad%g2rd(idx+i) = dvxxd + dvyyd + dvzzd

            END DO
         END IF
      END IF

   END SUBROUTINE mkgxyz3
END MODULE m_mkgxyz3
