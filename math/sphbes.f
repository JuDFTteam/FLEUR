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
      INTEGER i,l,min,n
      REAL fac,quot,xinv,xx
      REAL :: aux(0:int(lmax+10+x))
      !REAL, ALLOCATABLE :: aux(:)
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC abs,cos,sin

      IF (lmax.LT.1)  CALL juDFT_error("sphbes1",calledby="sphbes")
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
       !  ALLOCATE( aux(0:n) )
!
! downward recursion from arbitrary starting values
!
         aux(n) = 0.
         aux(n-1) = 1.
         DO i = n - 1,1,-1
            aux(i-1) = (2*i+1)*xinv*aux(i) - aux(i+1)
         ENDDO
!
! normalize with j0 or j1, whichever is larger
!
         min = 0
         IF (abs(fj(0)).LT.abs(fj(1))) min = 1
         quot = fj(min)/aux(min)
         DO l = 2,lmax
            fj(l) = aux(l)*quot
         ENDDO
        ! DEALLOCATE( aux )
      END IF

      RETURN
      END SUBROUTINE sphbes
      END MODULE m_sphbes
