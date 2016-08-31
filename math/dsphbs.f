      MODULE m_dsphbs
!********************************************************************
!     calculates the derivative of the spherical bessel functions
!     dfj(l) = d jl(x)/dx
!     for l=0,lmax and argument x
!     note that the spherical bessel functions fj(l), l=0,lmax are
!     needed (call sphbes to generate these)
!                                            m. weinert
!********************************************************************
      CONTAINS
      SUBROUTINE dsphbs(
     >                  lmax,x,fj,
     <                  dfj)

      IMPLICIT NONE
!     ..
!     .. Arguments ..
      INTEGER, INTENT  (IN) :: lmax
      REAL,    INTENT  (IN) :: x
      REAL,    INTENT  (IN) ::  fj(0:lmax)
      REAL,    INTENT (OUT) :: dfj(0:lmax)
!
!     .. Parameters ..
      REAL,    PARAMETER :: xlim = 1.0e-04
!     ..
!     .. Local Scalars ..
      REAL fac,x2
      INTEGER l

      dfj(0) = -fj(1)
!--->    small x limit
      IF (x.LT.xlim) THEN
         x2 = 0.5*x*x
         fac = 1./3.
         DO l = 1,lmax
            dfj(l) = fac* (l-x2* (l+2)/ (2*l+3))
            fac = x*fac/ (2*l+3)
         ENDDO
      ELSE
!--->    obtain dfj using recurrence relationship
         DO l = 1,lmax
            dfj(l) = fj(l-1) - (l+1)*fj(l)/x
         ENDDO
      END IF

      RETURN
      END SUBROUTINE dsphbs
      END MODULE m_dsphbs
