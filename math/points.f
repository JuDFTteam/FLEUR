      MODULE m_points
      CONTAINS
      SUBROUTINE points(x,n)
c     *********************************************************
c     generate random points, in internal coordinates,
c     within the unit cell omega-tilda
c     *********************************************************
      USE m_qranf
      IMPLICIT NONE
C     .. Scalar Arguments ..
      INTEGER n
C     ..
C     .. Array Arguments ..
      REAL x(3,n)
C     ..
C     .. Local Scalars ..
      REAL r
      INTEGER i,j
C     ..

C     ..
C     .. Intrinsic Functions ..
      INTRINSIC sqrt
C     ..
      r = sqrt(13.)
      j = 1
      DO  i = 1,n
         x(1,i) = qranf(r,j)
         x(2,i) = qranf(r,j)
         x(3,i) = qranf(r,j) - 0.5
      ENDDO
      RETURN
      END SUBROUTINE
      END
