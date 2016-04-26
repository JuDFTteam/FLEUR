      MODULE m_qranf
      CONTAINS
      REAL FUNCTION qranf(x,j)
c     **********************************************************
c     quasi random generator in the interval (0.,1.)
c     **********************************************************

      IMPLICIT NONE
C     .. Scalar Arguments ..
      REAL x
      INTEGER j
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC aint
C     ..
      j = j + 1
      qranf = j*x
      qranf = qranf - aint(qranf)
      RETURN
      END FUNCTION
      END
