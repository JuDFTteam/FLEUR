!--------------------------------------------------------------------------------
! Copyright (c) 2016 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------

      MODULE m_cylbes     
      use m_juDFT
c********************************************************************
c generates cylindrical Bessel functions for a given x and for the orders 
c from mmax to -mmax
c                                  Y. Mokrousov
c********************************************************************
      CONTAINS
      SUBROUTINE cylbes(
     >                  mmax,x,
     <                  fJ) 

      IMPLICIT NONE
!     ..
!     ..Arguments ..
      INTEGER, INTENT  (IN) :: mmax
      REAL,    INTENT  (IN) :: x
      REAL,    INTENT (OUT) :: fJ(-mmax:mmax)
!
!     .. Parameters ..
      REAL,    PARAMETER :: zero = 0.0
!     ..Locals ..
      INTEGER :: m,i,mass
      REAL :: quot
      REAL, ALLOCATABLE :: aux(:)
!     ..

      IF (x.LT.zero)  CALL juDFT_error("cylbes2",calledby="cylbes")

      IF (x.EQ.zero) THEN
        fJ(0) = 1.

        DO m=1,mmax
         fJ(m) = 0.
         fJ(-m) = 0.
        END DO
        RETURN
      END IF 

      mass = INT( mmax + 50 + x )
      ALLOCATE ( aux(0:mass) )       
      aux(mass) = 0.0
      aux(mass-1) = 1.0e-22   
     
      DO i=mass-2,0,-1
        aux(i) = 2*(i+1)*aux(i+1)/x - aux(i+2)
      END DO

      quot = aux(0)

      DO i=1,INT( mass/2. )
        quot = quot + 2*aux(2*i)
      END DO  

      fJ(0) = aux(0)/quot

      DO m=1,mmax
        fJ(m) = aux(m)/quot
        fJ(-m) = ((-1)**m)*fJ(m)
      END DO

      DEALLOCATE ( aux )

      RETURN
      END SUBROUTINE cylbes
      END MODULE m_cylbes 
 





