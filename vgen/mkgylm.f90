!--------------------------------------------------------------------------------
! Copyright (c) 2016 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------
MODULE m_mkgylm
   !-----------------------------------------------------------------------------
   ! Using the components and derivatives of a charge density rho in spherical
   ! coordinates on the real space grid, make the following quantaties:
   !
   ! gr(js):      [partial_r (rho),partial_theta (rho),partial_phi (rho)]
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
   ! rh is rho, rhd[i][j] are the partial derivatives along one/two directions.
   ! rv and thet are the radius and polar angles and t/f in derivatives stand for
   ! theta/phi.
   !
   ! Modified so only allocated old quantities require calculations. A.N. 2019
   !
   ! Quantities fo libxc are calculated as well. D.W./M.R. 2018
   !
   ! Original script by T.A. 1996
   !-----------------------------------------------------------------------------
CONTAINS
   SUBROUTINE mkgylm(jspins,rv,thet,nsp, rh,rhdr,rhdt,rhdf,rhdrr,rhdtt,rhdff, rhdtf,rhdrt,rhdrf,grad,kt)
      USE m_types
      IMPLICIT NONE
      REAL, INTENT (IN)                :: rv
      REAL, INTENT (IN)                :: thet(:)
      REAL, INTENT (IN)                :: rh(:,:)
      REAL, INTENT (IN)                :: rhdr(:,:),rhdt(:,:),rhdf(:,:)
      REAL, INTENT (IN)                :: rhdrr(:,:),rhdtt(:,:),rhdff(:,:)
      REAL, INTENT (IN)                :: rhdtf(:,:),rhdrf(:,:),rhdrt(:,:)
      TYPE(t_gradients), INTENT(INOUT) :: grad
      INTEGER, INTENT(IN)              :: jspins,nsp
      INTEGER, INTENT(IN)              :: kt !index of first point to use in gradients

      REAL sc,sint1,tant1,dagrru,dagrtu,dagrfu,dagrtd,dagrrd,dagrfd
      REAL, PARAMETER  :: chsml = 1.e-10
      INTEGER i,js,fac

      ! Gradients for libxc calculations.
      IF (ALLOCATED(grad%gr)) THEN
         DO js=1,jspins
            DO i=1,nsp
               grad%gr(:,kt+i,js)=[rhdr(i,js),rhdt(i,js)/rv,rhdf(i,js)/(rv*sin(thet(i)))]
            END DO
         END DO
         ! Use contracted gradients only for libxc.
         IF (ALLOCATED(grad%sigma)) THEN
            IF (jspins==1) THEN
               DO i=1,nsp
                  grad%sigma(1,kt+i)= dot_PRODUCT(grad%gr(:,kt+i,1),grad%gr(:,kt+i,1))
               END DO
            ELSE
               DO i=1,nsp
                  grad%sigma(1,kt+i)= dot_PRODUCT(grad%gr(:,kt+i,1),grad%gr(:,kt+i,1))
                  grad%sigma(2,kt+i)= dot_PRODUCT(grad%gr(:,kt+i,1),grad%gr(:,kt+i,2))
                  grad%sigma(3,kt+i)= dot_PRODUCT(grad%gr(:,kt+i,2),grad%gr(:,kt+i,2))
               END DO
            END IF
         END IF
         IF (ALLOCATED(grad%laplace)) THEN
            !Lapacian of density
            !fac=MERGE(2,1,jspins==1)
            fac=1
            DO js=1,jspins
               DO i=1,nsp
                  grad%laplace(kt+i,js) = (rhdrr(i,js) + 2.e0*rhdr(i,js)/rv +&
                                           (rhdtt(i,js)+rhdt(i,js)/TAN(thet(i))+rhdff(i,js)/SIN(thet(i))**2)/rv**2)/fac
               ENDDO
            ENDDO
         ENDIF
         ! Cartesian components of the gradient for sourcefree calculations. TODO: Need phi here as well.
         !IF (ALLOCATED(grad%grxyz)) THEN
         !   grad%grxyz=SIN(th)*COS(ph)*grad%gr(1,kt+k,1) + COS(th)*COS(ph)*grad%gr(2,kt+k,1)/r - SIN(ph)*grad%gr(3,kt+k,1)/(r*SIN(th))
         !END IF
         RETURN ! Do not calculate arrays for in-build GGA.
      END IF

      !     Old code for in-build xcpots
      IF(allocated(grad%agrt)) THEN
        sc=merge(0.5,1.0,jspins==1)
        DO i = 1,nsp
          IF ((rh(i,1) + rh(i,jspins))*sc<chsml) THEN
            grad%agrt(kt+i)  = 0.0
            grad%agru(kt+i)  = 0.0
            grad%agrd(kt+i)  = 0.0
            grad%g2rt(kt+i)  = 0.0
            grad%g2ru(kt+i)  = 0.0
            grad%g2rd(kt+i)  = 0.0
            grad%gggrt(kt+i) = 0.0
            grad%gggru(kt+i) = 0.0
            grad%gggrd(kt+i) = 0.0
            grad%grgru(kt+i) = 0.0
            grad%grgrd(kt+i) = 0.0
            grad%gzgr(kt+i)  = 0.0
          ELSE
            sint1 = sin(thet(i))
            tant1 = tan(thet(i))

            grad%agrt(kt+i) = sqrt((rhdr(i,1)*sc  + rhdr(i,jspins)*sc)**2+&
                              ((rhdt(i,1)*sc  + rhdt(i,jspins)*sc)/rv)**2+&
                      ((rhdf(i,1)*sc  + rhdf(i,jspins)*sc)/(rv*sint1))**2)

            IF (grad%agrt(kt+i)<chsml) THEN
              grad%agru(kt+i)  = 0.0
              grad%agrd(kt+i)  = 0.0
              grad%g2rt(kt+i)  = 0.0
              grad%g2ru(kt+i)  = 0.0
              grad%g2rd(kt+i)  = 0.0
              grad%gggrt(kt+i) = 0.0
              grad%gggru(kt+i) = 0.0
              grad%gggrd(kt+i) = 0.0
              grad%grgru(kt+i) = 0.0
              grad%grgrd(kt+i) = 0.0
              grad%gzgr(kt+i)  = 0.0
            ELSE
              grad%g2rt(kt+i)= (rhdrr(i,1)*sc  + rhdrr(i,jspins)*sc)+2.0*(rhdr(i,1)*sc  &
              + rhdr(i,jspins)*sc)/rv+((rhdtt(i,1)*sc  + rhdtt(i,jspins)*sc)+(rhdt(i,1)*sc &
              + rhdt(i,jspins)*sc)/tant1+(rhdff(i,1)*sc  + rhdff(i,jspins)*sc)/sint1**2)/rv**2

              grad%gggrt(kt+i) = (rhdr(i,1)*sc  + rhdr(i,jspins)*sc)&
              *(((rhdr(i,1)*sc  + rhdr(i,jspins)*sc)*(rhdrr(i,1)*sc &
              +rhdrr(i,jspins)*sc)*rv**3+(rhdt(i,1)*sc  + rhdt(i,jspins)*sc) &
              * ((rhdrt(i,1)*sc  + rhdrt(i,jspins)*sc)*rv-(rhdt(i,1)*sc  + rhdt(i,jspins)*sc))&
              + (rhdf(i,1)*sc  + rhdf(i,jspins)*sc)* ((rhdrf(i,1)*sc  + rhdrf(i,jspins)*sc)*rv&
              -(rhdf(i,1)*sc  + rhdf(i,jspins)*sc))/sint1**2)/grad%agrt(kt+i)/rv**3) &
              + ((rhdt(i,1)*sc  + rhdt(i,jspins)*sc)/rv)*(((rhdr(i,1)*sc  + rhdr(i,jspins)*sc)&
              *(rhdrt(i,1)*sc  + rhdrt(i,jspins)*sc)*rv**2+(rhdt(i,1)*sc  + rhdt(i,jspins)*sc)*(rhdtt(i,1)*sc  &
              + rhdtt(i,jspins)*sc)+(rhdf(i,1)*sc  + rhdf(i,jspins)*sc)&
              * (-1.*(rhdf(i,1)*sc  + rhdf(i,jspins)*sc)/tant1+(rhdtf(i,1)*sc  &
              + rhdtf(i,jspins)*sc))/sint1**2)/ (grad%agrt(kt+i)*rv**3)) &
              + ((rhdf(i,1)*sc  + rhdf(i,jspins)*sc)/(rv*sint1))*(((rhdr(i,1)*sc  + rhdr(i,jspins)*sc)*(rhdrf(i,1)*sc  &
              + rhdrf(i,jspins)*sc)*rv**2+(rhdt(i,1)*sc  + rhdt(i,jspins)*sc)*(rhdtf(i,1)*sc  &
              + rhdtf(i,jspins)*sc)+(rhdf(i,1)*sc  + rhdf(i,jspins)*sc)*(rhdff(i,1)*sc  &
              + rhdff(i,jspins)*sc)/sint1**2)/(grad%agrt(kt+i)*rv**3*sint1))

               grad%gzgr(kt+i) = ((rhdr(i,1)*sc-rhdr(i,jspins)*sc)*(rh(i,1) + rh(i,jspins))*sc &
               - (rh(i,1)*sc-rh(i,jspins)*sc)*(rhdr(i,1)*sc  + rhdr(i,jspins)*sc))/((rh(i,1) &
               + rh(i,jspins))*sc)**2*(rhdr(i,1)*sc  + rhdr(i,jspins)*sc)


               grad%agru(kt+i) = sqrt((rhdr(i,1)*sc)**2+(rhdt(i,1)*sc/rv)**2+(rhdf(i,1)*sc/(rv*sint1))**2)

               dagrru = (rhdr(i,1)*sc*rhdrr(i,1)*sc*rv**3+rhdt(i,1)*sc* (rhdrt(i,1)*sc*rv-rhdt(i,1)*sc)+&
                         rhdf(i,1)*sc* (rhdrf(i,1)*sc*rv-rhdf(i,1)*sc)/sint1**2)/grad%agru(kt+i)/rv**3

               dagrtu = (rhdr(i,1)*sc*rhdrt(i,1)*sc*rv**2+rhdt(i,1)*sc*rhdtt(i,1)*sc&
               + rhdf(i,1)*sc* (-rhdf(i,1)*sc/tant1+rhdtf(i,1)*sc)/sint1**2)/ (grad%agru(kt+i)*rv**3)

               dagrfu = (rhdr(i,1)*sc*rhdrf(i,1)*sc*rv**2+rhdt(i,1)*sc*rhdtf(i,1)*sc&
               +rhdf(i,1)*sc*rhdff(i,jspins)*sc/sint1**2)/ (grad%agru(kt+i)*rv**3*sint1)

               grad%g2ru(kt+i) = rhdrr(i,1)*sc + 2.e0*rhdr(i,1)*sc/rv + (rhdtt(i,1)*sc&
               +rhdt(i,1)*sc/tant1+rhdff(i,jspins)*sc/sint1**2)/rv**2

               grad%gggru(kt+i) = rhdr(i,1)*sc*dagrru + rhdt(i,1)*sc/rv*dagrtu + rhdf(i,1)*sc/(rv*sint1)*dagrfu

               grad%grgru(kt+i) = (rhdr(i,1)*sc  + rhdr(i,jspins)*sc)*rhdr(i,1)*sc &
               + ((rhdt(i,1)*sc  + rhdt(i,jspins)*sc)/rv)*rhdt(i,1)*sc/rv &
               + ((rhdf(i,1)*sc  + rhdf(i,jspins)*sc)/(rv*sint1))*rhdf(i,1)*sc/(rv*sint1)

               grad%agrd(kt+i) = sqrt((rhdr(i,jspins)*sc)**2+(rhdt(i,jspins)*sc/rv)**2+(rhdf(i,jspins)*sc/(rv*sint1))**2)

               dagrrd = (rhdr(i,jspins)*sc*rhdrr(i,jspins)*sc*rv**3+rhdt(i,jspins)*sc* (rhdrt(i,jspins)*sc*rv-rhdt(i,jspins)*sc)&
               + rhdf(i,jspins)*sc* (rhdrf(i,jspins)*sc*rv-rhdf(i,jspins)*sc)/sint1**2)/grad%agrd(kt+i)/rv**3

               dagrtd = (rhdr(i,jspins)*sc*rhdrt(i,jspins)*sc*rv**2+rhdt(i,jspins)*sc*rhdtt(i,jspins)*sc&
               + rhdf(i,jspins)*sc* (-rhdf(i,jspins)*sc/tant1+rhdtf(i,jspins)*sc)/sint1**2)/ (grad%agrd(kt+i)*rv**3)

               dagrfd = (rhdr(i,jspins)*sc*rhdrf(i,jspins)*sc*rv**2+rhdt(i,jspins)*sc*rhdtf(i,jspins)*sc&
               +rhdf(i,jspins)*sc*rhdff(i,jspins)*sc/sint1**2)/ (grad%agrd(kt+i)*rv**3*sint1)

               grad%g2rd(kt+i) = rhdrr(i,jspins)*sc + 2*rhdr(i,jspins)*sc/rv + (rhdtt(i,jspins)*sc+rhdt(i,jspins)*sc/tant1&
               +rhdff(i,jspins)*sc/sint1**2)/rv**2

               grad%gggrd(kt+i) = rhdr(i,jspins)*sc*dagrrd + rhdt(i,jspins)*sc/rv*dagrtd + rhdf(i,jspins)*sc/(rv*sint1)*dagrfd

               grad%grgrd(kt+i) = (rhdr(i,1)*sc  + rhdr(i,jspins)*sc)*rhdr(i,jspins)*sc &
               + ((rhdt(i,1)*sc  + rhdt(i,jspins)*sc)/rv)*rhdt(i,jspins)*sc/rv &
               + ((rhdf(i,1)*sc  + rhdf(i,jspins)*sc)/(rv*sint1))*rhdf(i,jspins)*sc/(rv*sint1)
             ENDIF
           ENDIF
         ENDDO
      ENDIF
   END SUBROUTINE mkgylm
END MODULE m_mkgylm
