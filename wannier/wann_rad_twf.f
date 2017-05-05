!--------------------------------------------------------------------------------
! Copyright (c) 2016 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------

      module m_wann_rad_twf
      use m_juDFT
      contains
      subroutine wann_rad_twf(
     >               nwfs,jmtd,natd,ind,rwf,zona,regio,
     >               us,dus,uds,duds,ff,gg,lmaxd,ikpt,
     >               ntypd,ntp,jri,rmsh,dx,
     >               nlod,flo,llo,nlo,
     <               rads)
c*********************************************************************
c..this subroutine calculates on the radial grid the radial function
c..corresponding to the twf, its mt sphere parameters, diffusivity
c..and region
c
c..Y. Mokrousov
c
c..tidied-up && extension of rwf-parameter
c..Frank Freimuth
c*********************************************************************

      use m_intgr, only : intgr3

      implicit none

      integer, intent (in)  :: nwfs,natd,jmtd,ntypd,lmaxd,ikpt,nlod
      integer, intent (in)  :: rwf(nwfs),ind(nwfs),ntp(natd),jri(ntypd)
      real,    intent (in)  :: zona(nwfs),regio(nwfs),rmsh(jmtd,ntypd)
      real,    intent (in)  :: us(0:lmaxd,ntypd),dus(0:lmaxd,ntypd)
      real,    intent (in)  :: uds(0:lmaxd,ntypd),duds(0:lmaxd,ntypd)
      real,    intent (in)  :: ff(ntypd,jmtd,2,0:lmaxd) 
      real,    intent (in)  :: gg(ntypd,jmtd,2,0:lmaxd),dx(ntypd) 
      real,    intent (in)  :: flo(ntypd,jmtd,2,nlod)
      integer, intent (in)  :: llo(nlod,ntypd)
      integer, intent (in)  :: nlo(ntypd)
      real,    intent (out) :: rads(nwfs,0:3,jmtd,2)

      integer :: nwf,nat,ntyp,j,n,l,lo
      real    :: rr,alpha,const,fact,wronk,radi
      real    :: acft,bcft,aa,bb,a1,b1,rho,radf(jmtd)

      do nwf = 1,nwfs

       if (abs(regio(nwf)-1.00000).ge.0.00001) then
          CALL juDFT_error("no projections outside muffin tins",calledby
     +         ="wann_rad_twf")
       endif

       nat = ind(nwf)
       ntyp = ntp(nat)
       alpha = zona(nwf)

       rads (nwf,:,:,:) = 0.
 
       if (rwf(nwf).eq.-5) then 
c*********************************************************************
c..For rwf.eq.-5 the following way of implementing is realized.
c..Suppose, the twf has its center on some atom. The 
c..self-consistent potential V_l(r) and the energy parameters
c..E_l are known. A linear combination of the u_l and \dot{u}_l
c..is constructed, with the coefficients A and B such that 
c..on the mt boundary A*u_l + B*\dot{u}_l is smoothly continuous
c..with the 1s function 2*((zona)**(1.5))*exp(-r*zona)
c*********************************************************************
         aa = 2.*(zona(nwf)**(1.5))
         bb = zona(nwf)
        
         a1 = aa*exp(-bb*rmsh(jri(ntyp),ntyp))
         b1 = -aa*bb*exp(-bb*rmsh(jri(ntyp),ntyp))

         do l = 0,3

            wronk = us(l,ntyp)*duds(l,ntyp)-uds(l,ntyp)*dus(l,ntyp)

            acft = (a1*duds(l,ntyp) - b1*uds(l,ntyp))/wronk
            bcft = (b1*us(l,ntyp) - a1*dus(l,ntyp))/wronk
 
            do j = 1,jri(ntyp)
               rads(nwf,l,j,:) = acft*ff(ntyp,j,:,l) +
     +                           bcft*gg(ntyp,j,:,l)
            enddo 

         enddo
       elseif(rwf(nwf).eq.0) then
c***********************************************************
c..For rwf.eq.0 the zona parameter mixes the radial function
c..and its energy derivative linearly.
c***********************************************************
         do l = 0,3
            do j = 1,jri(ntyp)
               rads(nwf,l,j,:) = ff(ntyp,j,:,l)*(1.-abs(zona(nwf))) +
     +                         zona(nwf)*gg(ntyp,j,:,l)
            enddo 
         enddo   
c**************************************
c..project onto local orbitals
c**************************************
       elseif(rwf(nwf).eq.-6)then
         IF(nlod<1) CALL juDFT_error("nlod<1",calledby="wann_rad_twf"
     +         )
          do l=0,3
           do j=1,jri(ntyp)
            rads(nwf,l,j,:)=flo(ntyp,j,:,1)
           enddo!j
          enddo!l
       elseif(rwf(nwf).eq.-7)then
         IF(nlod<2) CALL juDFT_error("nlod<2",calledby="wann_rad_twf"
     +         )
          do l=0,3
           do j=1,jri(ntyp)
            rads(nwf,l,j,:)=flo(ntyp,j,:,2)
           enddo!j
          enddo!l
       elseif(rwf(nwf).eq.-8)then
          IF(nlod<3) CALL juDFT_error("nlod<3",calledby ="wann_rad_twf"
     +         )
          do l=0,3
           do j=1,jri(ntyp)
            rads(nwf,l,j,:)=flo(ntyp,j,:,3)
           enddo!j
          enddo!l
       elseif(rwf(nwf).eq.-9)then
          rads(nwf,:,:,:)=0
          do lo=1,nlo(ntyp)
           l=llo(lo,ntyp)
           do j=1,jri(ntyp)
            rads(nwf,l,j,:)=flo(ntyp,j,:,lo)     
           enddo!j
          enddo !lo
        
       elseif((-5.lt.rwf(nwf)).and.(rwf(nwf).lt.0)) then
c************************************************************
c..use the radial with l=abs(rwf)-1 for all components of the
c..trial wave function
c************************************************************
           do l = 0,3
            do j = 1,jri(ntyp)
               rads(nwf,l,j,:) = 
     =             ff(ntyp,j,:,abs(rwf(nwf))-1)*(1.-abs(zona(nwf))) +
     +             zona(nwf)*gg(ntyp,j,:,abs(rwf(nwf))-1)
            enddo 
           enddo   

c**************************************************************
c..Hydrogen-like test orbitals.
c..here we have tabulated the radial functions from orbitron.
c..for 1 <= rwf <= 6 different components of the test orbital
c..are constructed from different radials.
c  For Hybrides (e.g. sp3) it might be better to use the same
c  radials for all components. See below (7 <= rwf <= 12)
c**************************************************************
       elseif (rwf(nwf).eq.1) then
c..n=1
          do j = 1,jri(ntyp)
             rho = 2.*alpha*rmsh(j,ntyp)
             rads(nwf,0,j,1) = 2.*((alpha)**(1.5))*
     *            exp(-rho/2.)            
          enddo  

       elseif (rwf(nwf).eq.2) then
c..n=2
          do j = 1,jri(ntyp)
             rho = alpha*(rmsh(j,ntyp))
             rads(nwf,0,j,1) = (1./(2.*sqrt(2.)))*
     *              (2.-rho)*
     *              ((alpha)**(1.5))*exp(-rho/2.)
             rads(nwf,1,j,1) = (1./(2.*sqrt(6.)))*
     *              rho*
     *              ((alpha)**(1.5))*exp(-rho/2.)
          enddo          

       elseif (rwf(nwf).eq.3) then
c..n=3
          do j = 1,jri(ntyp)
             rho = 2*alpha*rmsh(j,ntyp)/3.
             rads(nwf,0,j,1) = (1./(9.*sqrt(3.)))*
     *              (6. - 6.*rho + rho**2)*
     *              ((alpha)**(1.5))*exp(-rho/2.)
             rads(nwf,1,j,1) = (1./(9.*sqrt(6.)))*
     *              (4. - rho)*rho*
     *              ((alpha)**(1.5))*exp(-rho/2.)
             rads(nwf,2,j,1) = (1./(9.*sqrt(30.)))*
     *               rho*rho*
     *              ((alpha)**(1.5))*exp(-rho/2.)
          enddo

       elseif (rwf(nwf).eq.4) then
c..n=4
          do j = 1,jri(ntyp)
             rho = 2*alpha*rmsh(j,ntyp)/4.
             rads(nwf,0,j,1) = (1./96.)*
     *              (24. - 36.*rho + 12.*(rho**2) - (rho**3))*
     *              ((alpha)**(1.5))*exp(-rho/2.)
             rads(nwf,1,j,1) = (1./(32.*sqrt(15.)))*
     *              rho*(20. - 10.*rho + rho**2)*
     *              ((alpha)**(1.5))*exp(-rho/2.)
             rads(nwf,2,j,1) = (1./(96.*sqrt(5.)))*
     *              rho*rho*(6. - rho)*
     *              ((alpha)**(1.5))*exp(-rho/2.)
             rads(nwf,3,j,1) = (1./(96.*sqrt(35.)))*
     *              rho*rho*rho*
     *              ((alpha)**(1.5))*exp(-rho/2.)
          enddo

       elseif (rwf(nwf).eq.5) then
          do j = 1,jri(ntyp)
             rho = 2*alpha*rmsh(j,ntyp)/5.
             rads(nwf,0,j,1) = (1./(300.*sqrt(5.)))*
     *              (120. - 240.*rho + 120.*(rho**2) - 20.*(rho**3)
     +                            +(rho**4))*
     *              ((alpha)**(1.5))*exp(-rho/2.)
             rads(nwf,1,j,1) = (1./(150.*sqrt(30.)))*
     *              rho*(120. - 90.*rho + 18.*(rho**2) - (rho**3))*
     *              ((alpha)**(1.5))*exp(-rho/2.)
             rads(nwf,2,j,1) = (1./(150.*sqrt(70.)))*
     *              rho*rho*(42. - 14.*rho + (rho*rho))*
     *              ((alpha)**(1.5))*exp(-rho/2.)
             rads(nwf,3,j,1) = (1./(300.*sqrt(70.)))*
     *              rho*rho*rho*(8. - rho)*
     *              ((alpha)**(1.5))*exp(-rho/2.)
          enddo
        elseif (rwf(nwf).eq.6) then

          do j = 1,jri(ntyp)
             rho = 2*alpha*rmsh(j,ntyp)/6.
             rads(nwf,0,j,1) = (1./2160.*(sqrt(6.)))*
     *              (720. - 1800.*rho + 1200.*rho*rho-
     -               300.*rho*rho*rho + 30.*(rho**4) - (rho**5))*
     *              ((alpha)**(1.5))*exp(-rho/2.)
             rads(nwf,1,j,1) = (1./432.*(sqrt(210.)))*
     *              rho*(840. - 840.*rho + 252.*rho*rho-
     -               28.*rho*rho*rho + (rho**4))*
     *              ((alpha)**(1.5))*exp(-rho/2.)
          enddo

        elseif (rwf(nwf).eq.7) then
c..n=1
          do j = 1,jri(ntyp)
             rho = 2.*alpha*rmsh(j,ntyp)
             do l=0,3
             rads(nwf,l,j,1) = 2.*((alpha)**(1.5))*
     *            exp(-rho/2.)            
             enddo
          enddo  

        elseif (rwf(nwf).eq.8) then
c..n=2
          do j = 1,jri(ntyp)
             rho = alpha*(rmsh(j,ntyp))
             do l=0,3
             rads(nwf,l,j,1) = (1./(2.*sqrt(2.)))*
     *              (2.-rho)*
     *              ((alpha)**(1.5))*exp(-rho/2.)
             enddo
          enddo          

        elseif (rwf(nwf).eq.9) then
c..n=3
          do j = 1,jri(ntyp)
             rho = 2*alpha*rmsh(j,ntyp)/3.
             do l=0,3
             rads(nwf,l,j,1) = (1./(9.*sqrt(3.)))*
     *              (6. - 6.*rho + rho**2)*
     *              ((alpha)**(1.5))*exp(-rho/2.)
             enddo
          enddo

        elseif (rwf(nwf).eq.10) then
c..n=4
          do j = 1,jri(ntyp)
             rho = 2*alpha*rmsh(j,ntyp)/4.
             do l=0,3
             rads(nwf,l,j,1) = (1./96.)*
     *              (24. - 36.*rho + 12.*(rho**2) - (rho**3))*
     *              ((alpha)**(1.5))*exp(-rho/2.)
             enddo
          enddo
        elseif (rwf(nwf).eq.11) then

          do j = 1,jri(ntyp)
             rho = 2*alpha*rmsh(j,ntyp)/5.
             do l=0,3
             rads(nwf,l,j,1) = (1./(300.*sqrt(5.)))*
     *              (120. - 240.*rho + 120.*(rho**2) - 20.*(rho**3)
     +                            +(rho**4))*
     *              ((alpha)**(1.5))*exp(-rho/2.)
             enddo
          enddo

        elseif (rwf(nwf).eq.12) then

          do j = 1,jri(ntyp)
             rho = 2*alpha*rmsh(j,ntyp)/6.
             do l=0,3
             rads(nwf,0,j,1) = 0*(1./2160.*(sqrt(6.)))*
     *              (720. - 1800.*rho + 1200.*rho*rho-
     -               300.*rho*rho*rho + 30.*(rho**4) - (rho**5))*
     *              ((alpha)**(1.5))*exp(-rho/2.)
             enddo
          enddo

        else 

           CALL juDFT_error("radial function is not tabulated" ,calledby
     +          ="wann_rad_twf")


       endif

       if (ikpt.eq.1) then
        do l = 0,3
         do j = 1,jri(ntyp)
          radf(j) = rmsh(j,ntyp)*rmsh(j,ntyp)*rads(nwf,l,j,1)*
     *                rads(nwf,l,j,1)
c         radf(j) = rads(nwf,l,j,1)*rads(nwf,l,j,1)
         enddo
         call intgr3(radf,rmsh(1,ntyp),dx(ntyp),jri(ntyp),radi)
         write (6,*)
         write (6,*) 'Wannier Function N:',nwf
         write (6,*) 'angular momentum',l
         write (6,*) 'radial function at the MT boundary:',
     &                rads(nwf,l,jri(ntyp),1)
         write (6,*) 'norma =',radi
        enddo
       endif 
 
      enddo ! by twfs

      return

      end subroutine wann_rad_twf
      end module m_wann_rad_twf
