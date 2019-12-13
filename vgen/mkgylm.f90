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

      REAL dagrf,dagrfd,dagrfu,dagrr,dagrrd,dagrru,dagrt, &
           dagrtd,dagrtu,drdr,dzdfs,dzdr,dzdtr,grf,grfd,grfu,grr, &
           grrd,grru,grt,grtd,grtu,rdf,rdfd,rdff,rdffd,rdffu,rdfu, &
           rdr,rdrd,rdrf,rdrfd,rdrfu,rdrrd,rdrru,rdrt,rdrtd,rdrtu, &
           rdru,rdt,rdtd,rdtf,rdtfd,rdtfu,rdtt,rdttd,rdttu,rdtu, &
           ro,ro2,rod,rou,rv1,rv2,rv3,rvsin1,sint1,sint2,tant1
      REAL, PARAMETER  :: chsml = 1.e-10
      INTEGER i,js,fac

      rv1 = rv
      rv2 = rv1**2
      rv3 = rv1**3

      ! Gradients for libxc calculations.
      IF (ALLOCATED(grad%gr)) THEN
         DO js=1,jspins
            DO i=1,nsp
               grad%gr(:,kt+i,js)=[rhdr(i,js),rhdt(i,js),rhdf(i,js)]
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
                  grad%laplace(kt+i,js) = (rhdrr(i,js) + 2.e0*rhdr(i,js)/rv1 +&
                                           (rhdtt(i,js)+rhdt(i,js)/TAN(thet(i))+rhdff(i,js)/SIN(thet(i))**2)/rv2)/fac
               ENDDO
            ENDDO
         ENDIF
         ! Cartesian components of the gradient for sourcefree calculations. TODO: Need phi here as well.
         !IF (ALLOCATED(grad%grxyz)) THEN
         !   grad%grxyz=SIN(th)*COS(ph)*grad%gr(1,kt+k,1) + COS(th)*COS(ph)*grad%gr(2,kt+k,1)/r - SIN(ph)*grad%gr(3,kt+k,1)/(r*SIN(th))
         !END IF
         RETURN ! Do not calculate arrays for in-build GGA.
      END IF

      IF (ALLOCATED(grad%gr)) THEN
         
      END IF
      !     Old code for in-build xcpots
      IF(allocated(grad%agrt)) THEN
         IF (jspins==1) THEN

            points_1 : DO i = 1,nsp

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

               ro = rh(i,1)

               IF (ro<chsml) CYCLE points_1
               sint1 = sin(thet(i))
               sint2 = sint1**2
               tant1 = tan(thet(i))
               rvsin1 = rv1*sint1

               rou   = ro/2
               rdru  = rhdr(i,1)/2
               rdtu  = rhdt(i,1)/2
               rdfu  = rhdf(i,1)/2
               rdrru = rhdrr(i,1)/2
               rdttu = rhdtt(i,1)/2
               rdffu = rhdff(i,1)/2
               rdtfu = rhdtf(i,1)/2
               rdrtu = rhdrt(i,1)/2
               rdrfu = rhdrf(i,1)/2

               rod   = rou
               rdrd  = rdru
               rdtd  = rdtu
               rdfd  = rdfu
               rdrrd = rdrru
               rdttd = rdttu
               rdffd = rdffu
               rdtfd = rdtfu
               rdrtd = rdrtu
               rdrfd = rdrfu

               rdr  = rdru  + rdrd
               rdt  = rdtu  + rdtd
               rdf  = rdfu  + rdfd
               drdr = rdrru + rdrrd
               rdtt = rdttu + rdttd
               rdff = rdffu + rdffd
               rdtf = rdtfu + rdtfd
               rdrt = rdrtu + rdrtd
               rdrf = rdrfu + rdrfd

               ro2 = ro**2

               grr = rdr
               grt = rdt/rv1
               grf = rdf/rvsin1

               grad%agrt(kt+i) = sqrt(grr**2+grt**2+grf**2)

               IF (grad%agrt(kt+i)<chsml) CYCLE points_1

               dagrr = (rdr*drdr*rv3+rdt* (rdrt*rv1-rdt)+ &
                        rdf* (rdrf*rv1-rdf)/sint2)/grad%agrt(kt+i)/rv3

               dagrt =(rdr*rdrt*rv2+rdt*rdtt+rdf* (-rdf/tant1+rdtf)/sint2)/&
                       (grad%agrt(kt+i)*rv3)

               dagrf = (rdr*rdrf*rv2+rdt*rdtf+rdf*rdff/sint2)/&
                       (grad%agrt(kt+i)*rv3*sint1)

               grad%g2rt(kt+i)=drdr+2.0*rdr/rv1+&
                                (rdtt+rdt/tant1+rdff/sint2)/rv2

               dzdr = ((rdru-rdrd)*ro- (rou-rod)*rdr)/ro2

               !--   >      dzdtr,dzdfs vanish by definition.

               dzdtr = 0.0
               dzdfs = 0.0

               grad%gggrt(kt+i) = grr*dagrr + grt*dagrt + grf*dagrf

               grad%gzgr(kt+i) = dzdr*grr + dzdtr*grt + dzdfs*grf

               grru = rdru
               grtu = rdtu/rv1
               grfu = rdfu/rvsin1

               grad%agru(kt+i) = sqrt(grru**2+grtu**2+grfu**2)

               dagrru = (rdru*rdrru*rv3+rdtu* (rdrtu*rv1-rdtu)+&
                         rdfu* (rdrfu*rv1-rdfu)/sint2)/grad%agru(kt+i)/rv3

               dagrtu = (rdru*rdrtu*rv2+rdtu*rdttu+&
                         rdfu* (-rdfu/tant1+rdtfu)/sint2)/ (grad%agru(kt+i)*rv3)

               dagrfu = (rdru*rdrfu*rv2+rdtu*rdtfu+rdfu*rdffu/sint2)/&
                        (grad%agru(kt+i)*rv3*sint1)

               grad%g2ru(kt+i) = rdrru + 2.e0*rdru/rv1 +&
                                 (rdttu+rdtu/tant1+rdffu/sint2)/rv2

               grad%gggru(kt+i) = grru*dagrru + grtu*dagrtu + grfu*dagrfu

               grad%grgru(kt+i) = grr*grru + grt*grtu + grf*grfu

               grrd = rdrd
               grtd = rdtd/rv1
               grfd = rdfd/rvsin1

               grad%agrd(kt+i) = sqrt(grrd**2+grtd**2+grfd**2)

               dagrrd = (rdrd*rdrrd*rv3+rdtd* (rdrtd*rv1-rdtd)+&
                         rdfd* (rdrfd*rv1-rdfd)/sint2)/grad%agrd(kt+i)/rv3

               dagrtd = (rdrd*rdrtd*rv2+rdtd*rdttd+&
                         rdfd* (-rdfd/tant1+rdtfd)/sint2)/ (grad%agrd(kt+i)*rv3)

               dagrfd = (rdrd*rdrfd*rv2+rdtd*rdtfd+rdfd*rdffd/sint2)/&
                        (grad%agrd(kt+i)*rv3*sint1)

               grad%g2rd(kt+i) = rdrrd + 2*rdrd/rv1 +&
                                 (rdttd+rdtd/tant1+rdffd/sint2)/rv2

               grad%gggrd(kt+i) = grrd*dagrrd + grtd*dagrtd + grfd*dagrfd

               grad%grgrd(kt+i) = grr*grrd + grt*grtd + grf*grfd

            ENDDO points_1

         ELSE

            points_2 : DO i = 1,nsp

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

               ro = rh(i,1) + rh(i,jspins)

               IF (ro<chsml) CYCLE points_2

               sint1 = sin(thet(i))
               sint2 = sint1**2
               tant1 = tan(thet(i))
               rvsin1 = rv1*sint1

               rou   = rh(i,1)
               rdru  = rhdr(i,1)
               rdtu  = rhdt(i,1)
               rdfu  = rhdf(i,1)
               rdrru = rhdrr(i,1)
               rdttu = rhdtt(i,1)
               rdffu = rhdff(i,1)
               rdtfu = rhdtf(i,1)
               rdrtu = rhdrt(i,1)
               rdrfu = rhdrf(i,1)

               rod   = rh(i,jspins)
               rdrd  = rhdr(i,jspins)
               rdtd  = rhdt(i,jspins)
               rdfd  = rhdf(i,jspins)
               rdrrd = rhdrr(i,jspins)
               rdttd = rhdtt(i,jspins)
               rdffd = rhdff(i,jspins)
               rdtfd = rhdtf(i,jspins)
               rdrtd = rhdrt(i,jspins)
               rdrfd = rhdrf(i,jspins)

               rdr   = rdru  + rdrd
               rdt   = rdtu  + rdtd
               rdf   = rdfu  + rdfd
               drdr  = rdrru + rdrrd
               rdtt  = rdttu + rdttd
               rdff  = rdffu + rdffd
               rdtf  = rdtfu + rdtfd
               rdrt  = rdrtu + rdrtd
               rdrf  = rdrfu + rdrfd

               ro2 = ro**2

               grr = rdr
               grt = rdt/rv1
               grf = rdf/rvsin1

               grad%agrt(kt+i) = sqrt(grr**2+grt**2+grf**2)

               IF (grad%agrt(kt+i)<chsml) CYCLE points_2

               dagrr = (rdr*drdr*rv3+rdt* (rdrt*rv1-rdt)+ rdf* (rdrf*rv1-rdf)/sint2)/grad%agrt(kt+i)/rv3

               dagrt =(rdr*rdrt*rv2+rdt*rdtt+rdf* (-rdf/tant1+rdtf)/sint2)/ (grad%agrt(kt+i)*rv3)

               dagrf = (rdr*rdrf*rv2+rdt*rdtf+rdf*rdff/sint2)/(grad%agrt(kt+i)*rv3*sint1)

               grad%g2rt(kt+i)= drdr+2.0*rdr/rv1+(rdtt+rdt/tant1+rdff/sint2)/rv2

               dzdr = ((rdru-rdrd)*ro- (rou-rod)*rdr)/ro2

               !     dzdtr,dzdfs vanish by definition.
               dzdtr = 0.0
               dzdfs = 0.0

               grad%gggrt(kt+i) = grr*dagrr + grt*dagrt + grf*dagrf

               grad%gzgr(kt+i) = dzdr*grr + dzdtr*grt + dzdfs*grf

               grru = rdru
               grtu = rdtu/rv1
               grfu = rdfu/rvsin1

               grad%agru(kt+i) = sqrt(grru**2+grtu**2+grfu**2)

               dagrru = (rdru*rdrru*rv3+rdtu* (rdrtu*rv1-rdtu)+&
                         rdfu* (rdrfu*rv1-rdfu)/sint2)/grad%agru(kt+i)/rv3

               dagrtu = (rdru*rdrtu*rv2+rdtu*rdttu+ rdfu* (-rdfu/tant1+rdtfu)/sint2)/ (grad%agru(kt+i)*rv3)

               dagrfu = (rdru*rdrfu*rv2+rdtu*rdtfu+rdfu*rdffu/sint2)/ (grad%agru(kt+i)*rv3*sint1)

               grad%g2ru(kt+i) = rdrru + 2.e0*rdru/rv1 + (rdttu+rdtu/tant1+rdffu/sint2)/rv2

               grad%gggru(kt+i) = grru*dagrru + grtu*dagrtu + grfu*dagrfu

               grad%grgru(kt+i) = grr*grru + grt*grtu + grf*grfu

               grrd = rdrd
               grtd = rdtd/rv1
               grfd = rdfd/rvsin1

               grad%agrd(kt+i) = sqrt(grrd**2+grtd**2+grfd**2)

               dagrrd = (rdrd*rdrrd*rv3+rdtd* (rdrtd*rv1-rdtd)+ rdfd* (rdrfd*rv1-rdfd)/sint2)/grad%agrd(kt+i)/rv3

               dagrtd = (rdrd*rdrtd*rv2+rdtd*rdttd+ rdfd* (-rdfd/tant1+rdtfd)/sint2)/ (grad%agrd(kt+i)*rv3)

               dagrfd = (rdrd*rdrfd*rv2+rdtd*rdtfd+rdfd*rdffd/sint2)/ (grad%agrd(kt+i)*rv3*sint1)

               grad%g2rd(kt+i) = rdrrd + 2*rdrd/rv1 + (rdttd+rdtd/tant1+rdffd/sint2)/rv2

               grad%gggrd(kt+i) = grrd*dagrrd + grtd*dagrtd + grfd*dagrfd

               grad%grgrd(kt+i) = grr*grrd + grt*grtd + grf*grfd

            ENDDO points_2

         ENDIF
      ENDIF
   END SUBROUTINE mkgylm
END MODULE m_mkgylm
