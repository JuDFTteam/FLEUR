      MODULE m_dcylbs
!********************************************************************
!     calculates the derivative of the cylindrical bessel functions
!     for the given argument x amd m=-mmax,mmax
!     for l=0,lmax and argument x
!     note that the cylindrical bessel functions fj(m), m=-mmax,mmax are
!     needed (call cylbes to generate those)
!                                            Y.Mokrousov
!********************************************************************
      CONTAINS
      SUBROUTINE dcylbs(
     >                  mmax,x,fJ,
     <                  dfJ)

      IMPLICIT NONE
!     ..
!     .. Arguments ..
      INTEGER, INTENT  (IN) :: mmax
      REAL,    INTENT  (IN) :: x
      REAL,    INTENT  (IN) ::  fJ(-mmax:mmax)
      REAL,    INTENT (OUT) :: dfJ(-mmax:mmax)
!
!     .. Parameters ..
      REAL,    PARAMETER :: xlim = 1.0e-04, zero = 0.0
!     ..
!     .. Local Scalars ..
      INTEGER m
      REAL a,b
     

      dfJ(0) = -fJ(1)

      IF (x .EQ. zero) THEN
         dfJ(1) = 0.5
         dfJ(-1) = -0.5        
         DO m = 2,mmax
           dfJ(m) = zero
           dfJ(-m) = zero
         END DO
      ELSE 
c---> possible calculation of the derivatives of cylindrical Bessel 
c---> functions using first two terms of series representation of
c---> them for small x, error is ~xlim**2          (Y.Mokrousov)
c         IF (x .LT. xlim) THEN
c           a = 0.5
c           b = 0.75
c           dfJ(1) = a - b*((x/2.)**2)
c           dfJ(-1) = -dfJ(1) 
c           DO m = 2,mmax
c              a = (1./(m-1))*a
c              b = b*(2+m)/((m+1)**2)
c              dfJ(m) = a*((x/2.)**(m-1)) - b*((x/2.)**(m+1))
c              dfJ(-m) = ((-1)**m)*dfJ(m)
c           END DO
c         ELSE
c------------------------------------------------------------------

         DO m = 1,mmax
            dfJ(m) = fJ(m-1) - m*fJ(m)/x
            dfJ(-m) = ((-1)**m)*dfJ(m)
         ENDDO

c     ENDIF

      END IF

      RETURN

      END SUBROUTINE dcylbs

      END MODULE m_dcylbs




