!--------------------------------------------------------------------------------
! Copyright (c) 2016 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------

      MODULE m_modcylk
      use m_juDFT
c     generates Makdonald's function K_m(x) for a given
c     order m and point x
c                         Y.Mokrousov
      CONTAINS
      SUBROUTINE modcylk(
     >                   m,x,
     <                   kJ) 
      

      implicit none
                             
      integer, intent  (in) :: m
      real,    intent  (in) :: x
      real,    intent (out) :: kJ
                               
      real,    parameter :: zero = 0.0
      real,    parameter :: gamma = 
     =     0.5772156649015328608
                               
      integer :: mm,i,mass,fact,j,m1
      real :: quot,t,k0,k1,kJ1,kJ2,psi,a
c     real, allocatable :: aux(:)
                                
    
      if (x.le.zero)  CALL juDFT_error("x.le.zero",calledby="modcylk")

c     K_0 and K_1 part
      
      if (x.gt.1.0) then
         t = 1./x
         k0 = exp(-x)*sqrt(t)*(1.2533141373 -
     -    0.1566641816*t      + 0.0881112782*(t**2)  -
     -    0.0913909546*(t**3) + 0.1344569228*(t**4)  -
     -    0.2299850328*(t**5) + 0.3792409730*(t**6)  -
     -    0.5247277331*(t**7) + 0.5575368367*(t**8)  -
     -    0.4262632912*(t**9) + 0.2184518096*(t**10) - 
     -    0.0668097672*(t**11)+ 0.0091893830*(t**12))
         k1 = exp(-x)*sqrt(t)*(1.2533141373 +
     +    0.4699927013*t      - 0.1468582957*(t**2)  +
     +    0.1280426636*(t**3) - 0.1736431637*(t**4)  +
     +    0.2847618149*(t**5) - 0.4594342117*(t**6)  +
     +    0.6283380681*(t**7) - 0.6632295430*(t**8)  +
     +    0.5050238576*(t**9) - 0.2581303765*(t**10) + 
     +    0.0788000118*(t**11)- 0.0108241775*(t**12))
      end if
      if (x.le.1.0) then
         t = 1./x
         k0 = - gamma - log(x/2.)
         k1 = t
         do i = 1,8
            psi = -gamma
            fact = 1
            do j = 1,i
               psi = psi + 1./j
               fact = fact*j
            end do
            if (i.le.6) then
            k0 = k0 + ((x/2.)**(2*i))*(psi-log(x/2.))/
     /           (fact*fact)
            end if
            k1 = k1 + ((x/2.)**(2*i-1))*(0.5 - i*
     *           (psi - log(x/2.)))/(fact*fact)
         end do
      end if
      if (m.eq.0) then
         kJ = k0
         return
      end if
      if (m.eq.1 .or. m.eq.-1) then 
         kJ = k1
         return
      end if

c     forward recursion

      kJ1 = k0
      kJ2 = k1
      m1 = int(abs(m))
      do mm = 2,m1
         a = kJ2
         kJ2 = 2*(mm-1)*t*kJ2 + kJ1
         kJ1 = a
      end do
      kJ = kJ2
     
      END SUBROUTINE modcylk
      END MODULE m_modcylk
