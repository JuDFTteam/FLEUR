      MODULE m_sphbes
      use m_juDFT
c********************************************************************
c calculate spherical Bessel functions and derivatives of sqrt(e)*r
c                                   P. Marksteiner and E. Badralexe
c********************************************************************
      CONTAINS
      SUBROUTINE sphbes(
     >                  lmax,x,
     <                  fj)

      IMPLICIT NONE
!     ..
!     .. Arguments ..
      INTEGER, INTENT  (IN) :: lmax
      REAL,    INTENT  (IN) :: x
      REAL,    INTENT (OUT) :: fj(0:lmax)
!
!     .. Parameters ..
      REAL,    PARAMETER :: small = 1.0e-03 , zero = 0.0
!     ..
!     .. Locals ..
      INTEGER i,l,n
      REAL fac,quot,xinv,xx
      REAL aux0, aux1, aux2
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC abs,cos,sin

      IF (lmax==0) THEN
         IF (x.GE.small) THEN
            fj(0) = sin(x)/x
         ELSE
            fj(0) = 1 - x*x/6.* (1.-x*x/20.* (1.-x*x/42.))
         END IF
         RETURN
      ENDIF
      IF (x.LT.zero)  CALL juDFT_error("sphbes2",calledby="sphbes")
      xx = x*x
      IF (x.GE.small) THEN
         xinv = 1./x
         fj(0) = sin(x)*xinv
         fj(1) = (fj(0)-cos(x))*xinv
      ELSE
         fj(0) = 1 - xx/6.* (1.-xx/20.* (1.-xx/42.))
         fj(1) = (x/3.)* (1.-xx/10.* (1.-xx/28.))
         fac = xx/15.
         DO l=2,lmax
           fj(l) = fac * ( 1. - xx / (4*l+6) )
           fac = x * fac / (2*l+3)
         ENDDO
         RETURN
      END IF
      IF (lmax.LT.x) THEN

         DO l = 2,lmax
            fj(l) = (2*l-1)*xinv*fj(l-1) - fj(l-2)
         ENDDO

      ELSE IF (lmax.GE.2) THEN
         n = INT( lmax + 10 + x )
!
! downward recursion from arbitrary starting values
!
         ! aux(n) = 0.
         ! aux(n-1) = 1.
         ! aux(i) = (2*i+3)*xinv*aux(i+1) - aux(i+2)
         aux1 = 0.
         aux0 = 1.
         DO i = n-2,lmax+1,-1
            aux2 = aux1
            aux1 = aux0
            aux0 = (2*i+3)*xinv*aux1 - aux2
         ENDDO
         DO i = lmax,2,-1
            aux2 = aux1
            aux1 = aux0
            aux0 = (2*i+3)*xinv*aux1 - aux2
            fj(i) = aux0
         ENDDO
         aux2 = aux1
         aux1 = aux0
         aux0 = (2*1+3)*xinv*aux1 - aux2
         aux2 = aux1
         aux1 = aux0
         aux0 = (2*0+3)*xinv*aux1 - aux2
!
! normalize with j0 or j1, whichever is larger
!
         IF (abs(fj(0)).LT.abs(fj(1))) THEN
            quot = fj(1)/aux1
         ELSE
            quot = fj(0)/aux0
         ENDIF
         DO l = 2,lmax
            fj(l) = fj(l)*quot
         ENDDO
      END IF

      RETURN
      END SUBROUTINE sphbes
      END MODULE m_sphbes
