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
   SUBROUTINE mkgxyz3(vl,dvx,dvy,dvz,dvxx,dvyy,dvzz,dvyz,dvxz,dvxy,grad)
      USE m_types
      IMPLICIT NONE
      REAL, INTENT (IN)                :: vl(:,:)
      REAL, INTENT (IN)                :: dvx(:,:),dvy(:,:),dvz(:,:)
      REAL, INTENT (IN)                :: dvxx(:,:),dvyy(:,:),dvzz(:,:)
      REAL, INTENT (IN)                :: dvyz(:,:),dvxz(:,:),dvxy(:,:)

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
               grad%gr(:,i,js)=(/dvx(i,js),dvy(i,js),dvz(i,js)/)
            END DO
         END DO
         ! Use contracted gradients only for libxc.
         IF(ALLOCATED(grad%sigma)) THEN
            IF (jspins==1) THEN
               DO i=1,nsp
                  grad%sigma(1,i) = dvx(i,1)*dvx(i,1) + dvy(i,1)*dvy(i,1) + dvz(i,1)*dvz(i,1)
               END DO
            ELSE
               DO i=1,nsp
                  grad%sigma(1,i) = dvx(i,1)*dvx(i,1) + dvy(i,1)*dvy(i,1) + dvz(i,1)*dvz(i,1)
                  grad%sigma(2,i) = dvx(i,1)*dvx(i,2) + dvy(i,1)*dvy(i,2) + dvz(i,1)*dvz(i,2)
                  grad%sigma(3,i) = dvx(i,2)*dvx(i,2) + dvy(i,2)*dvy(i,2) + dvz(i,2)*dvz(i,2)
               END DO
            END IF
         END IF
         IF(ALLOCATED(grad%laplace)) THEN
            DO js=1,jspins
               DO i=1,nsp
                  grad%laplace(i,js)= dvxx(i,js)+dvyy(i,js)+dvzz(i,js)
               END DO
            END DO
         END IF
         RETURN ! Do not calculate arrays for in-build GGA.
      END IF

      IF (ANY(SHAPE(vl).NE.SHAPE(dvx))) CALL judft_error("Gradients for internal GGA called with inconsistent sizes",hint="This is a bug")

 	   IF(ALLOCATED(grad%agrt)) THEN  
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

               grad%agrt(i) = max(sqrt(dvxt**2 + dvyt**2 + dvzt**2),sml)
               grad%agru(i) = max(sqrt(dvxu**2 + dvyu**2 + dvzu**2),sml)
               grad%agrd(i) = max(sqrt(dvxd**2 + dvyd**2 + dvzd**2),sml)

               dagrxt = (dvxt*dvxxt + dvyt*dvxyt + dvzt*dvxzt) / grad%agrt(i)
               dagrxu = (dvxu*dvxxu + dvyu*dvxyu + dvzu*dvxzu) / grad%agru(i)
               dagrxd = (dvxd*dvxxd + dvyd*dvxyd + dvzd*dvxzd) / grad%agrd(i)

               dagryt = (dvxt*dvxyt + dvyt*dvyyt + dvzt*dvyzt) / grad%agrt(i)
               dagryu = (dvxu*dvxyu + dvyu*dvyyu + dvzu*dvyzu) / grad%agru(i)
               dagryd = (dvxd*dvxyd + dvyd*dvyyd + dvzd*dvyzd) / grad%agrd(i)
   
               dagrzt = (dvxt*dvxzt + dvyt*dvyzt + dvzt*dvzzt) / grad%agrt(i)
               dagrzu = (dvxu*dvxzu + dvyu*dvyzu + dvzu*dvzzu) / grad%agru(i)
               dagrzd = (dvxd*dvxzd + dvyd*dvyzd + dvzd*dvzzd) / grad%agrd(i)

               grad%gggrt(i) = dvxt*dagrxt + dvyt*dagryt + dvzt*dagrzt
               grad%gggru(i) = dvxu*dagrxu + dvyu*dagryu + dvzu*dagrzu
               grad%gggrd(i) = dvxd*dagrxd + dvyd*dagryd + dvzd*dagrzd

               dzdx = (dvxu-dvxd)/vlt - (vlu-vld)*dvxt/vlt**2
               dzdy = (dvyu-dvyd)/vlt - (vlu-vld)*dvyt/vlt**2
               dzdz = (dvzu-dvzd)/vlt - (vlu-vld)*dvzt/vlt**2

               grad%gzgr(i) = dzdx*dvxt + dzdy*dvyt + dzdz*dvzt
   
               grad%g2rt(i) = dvxxt + dvyyt + dvzzt
               grad%g2ru(i) = dvxxu + dvyyu + dvzzu
               grad%g2rd(i) = dvxxd + dvyyd + dvzzd
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
   
               grad%agrt(i) = max(sqrt(dvxt**2 + dvyt**2 + dvzt**2),sml)
               grad%agru(i) = max(sqrt(dvxu**2 + dvyu**2 + dvzu**2),sml)
               grad%agrd(i) = max(sqrt(dvxd**2 + dvyd**2 + dvzd**2),sml)

               dagrxt = (dvxt*dvxxt + dvyt*dvxyt + dvzt*dvxzt) / grad%agrt(i)
               dagrxu = (dvxu*dvxxu + dvyu*dvxyu + dvzu*dvxzu) / grad%agru(i)
               dagrxd = (dvxd*dvxxd + dvyd*dvxyd + dvzd*dvxzd) / grad%agrd(i)

               dagryt = (dvxt*dvxyt + dvyt*dvyyt + dvzt*dvyzt) / grad%agrt(i)
               dagryu = (dvxu*dvxyu + dvyu*dvyyu + dvzu*dvyzu) / grad%agru(i)
               dagryd = (dvxd*dvxyd + dvyd*dvyyd + dvzd*dvyzd) / grad%agrd(i)

               dagrzt = (dvxt*dvxzt + dvyt*dvyzt + dvzt*dvzzt) / grad%agrt(i)
               dagrzu = (dvxu*dvxzu + dvyu*dvyzu + dvzu*dvzzu) / grad%agru(i)
               dagrzd = (dvxd*dvxzd + dvyd*dvyzd + dvzd*dvzzd) / grad%agrd(i)

               grad%gggrt(i) = dvxt*dagrxt + dvyt*dagryt + dvzt*dagrzt
               grad%gggru(i) = dvxu*dagrxu + dvyu*dagryu + dvzu*dagrzu
               grad%gggrd(i) = dvxd*dagrxd + dvyd*dagryd + dvzd*dagrzd

               dzdx = (dvxu-dvxd)/vlt -  (vlu-vld)*dvxt/vlt**2
               dzdy = (dvyu-dvyd)/vlt -  (vlu-vld)*dvyt/vlt**2
               dzdz = (dvzu-dvzd)/vlt -  (vlu-vld)*dvzt/vlt**2

               grad%gzgr(i) = dzdx*dvxt + dzdy*dvyt + dzdz*dvzt

               grad%g2rt(i) = dvxxt + dvyyt + dvzzt
               grad%g2ru(i) = dvxxu + dvyyu + dvzzu
               grad%g2rd(i) = dvxxd + dvyyd + dvzzd

            END DO
         END IF
      END IF

   END SUBROUTINE mkgxyz3
END MODULE m_mkgxyz3
