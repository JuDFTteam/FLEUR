!--------------------------------------------------------------------------------
! Copyright (c) 2016 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------

      MODULE m_modcyli
      use m_juDFT

c     generates modified bessel functions I_m for a given x and 
c     order m
c                Y.Mokrousov, 21 oct 2002
      CONTAINS
      SUBROUTINE modcyli(
     >                   m1,x,
     <                   iJ)


      USE m_constants, ONLY : pimach
      IMPLICIT NONE

      integer, intent  (in) :: m1
      real,    intent  (in) :: x
      real,    intent (out) :: iJ

      real, parameter :: zero = 0.0
      real, parameter :: twelve = 12.0

      integer :: mm,i,mass,m,j
      real    :: quot,sm,qm,iJ0,fact1,fact2,iiJ,nk,fact
      real, allocatable :: aux(:)

      intrinsic abs,int,sqrt

      m = int(abs(m1))
     
      iJ = 0.
      if (x.lt.zero)  CALL juDFT_error("x.lt.zero",calledby="modcyli")
      if (x.eq.zero) then
        if (m.eq.0) then
           iJ = 1.0
           return
        else
           iJ = 0.0
           return
        end if
      else
         sm = sqrt(real(m))
         qm = m*m
         if (x.le.twelve .and. x.le.4*sm) then
c            write (*,*) '<12,<4sqrt(n)'
            fact1 = 1.
            fact = 1.
            do i = 1,m
               fact = fact*i
            end do
            iJ = ((x/2.)**m)/fact
            fact2 = 1.
            do i=1,31
               fact1 = fact1*i
               fact2 = fact
               do j=1,i
                  fact2 = fact2*(m+j)
               end do
               iJ = iJ + ((x/2.)**m)*(((x*x)/4.)**i)/
     /              (fact1*fact2)
            end do
            return
         elseif (x.gt.twelve .and. x.ge.6*qm) then
c           write (*,*) '>12,>6n2'
            iJ = exp(x)/sqrt(2*pimach()*x)
            fact1 = 1.
            iiJ = 1.
            nk = 1.
            do i = 1,31
               fact1 = fact1*i
               nk = nk*(4*(m*m)-(2*i-1)*(2*i-1))
               iiJ = iiJ + ((-1)**i)*nk/(fact1*((8*x)**i))
            end do
            iJ = iJ*iiJ
            return
         else
            mass = int( m + 10 + x )
            allocate ( aux(0:mass) )       
            aux(mass) = 0.0
            aux(mass-1) = 1.0e-22  
            do i=mass-2,0,-1
               aux(i) = 2*(i+1)*aux(i+1)/x + aux(i+2)
c               if (i.lt.m .and. x.gt.6*qm .and. i.ne.0) then
c                  quot = aux(i)
c                  iJ = aux(m)/quot
c                  return
c               end if
            end do
            quot = aux(0)
            if (x.le.twelve) then
c               write (*,*) 'fucking here'
               fact1 = 1.                                             
               iJ0 = 1.
               do i=1,31
                  fact1 = fact1*i
                  iJ0 = iJ0 + (((x*x)/4.)**i)/(fact1*fact1)
               end do
               iJ = aux(m)*iJ0/quot
               deallocate (aux)
               return
            elseif (x.gt.twelve) then
               iJ0 = exp(x)/sqrt(2*pimach()*x)
               fact1 = 1.
               iiJ = 1.
               nk = 1.
               do i = 1,31
                  fact1 = fact1*i
                  nk = nk*(-(2*i-1)*(2*i-1))
                  iiJ = iiJ + ((-1)**i)*nk/(fact1*((8*x)**i))
               end do
               iJ0 = iJ0*iiJ
               iJ = aux(m)*iJ0/quot
               deallocate (aux)
               return
            end if
         end if
      end if

      END SUBROUTINE modcyli
      END MODULE m_modcyli
