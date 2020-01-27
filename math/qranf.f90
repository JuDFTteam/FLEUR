MODULE m_qranf
CONTAINS
   REAL FUNCTION qranf(x,j)
!     **********************************************************
!     quasi random generator in the interval (0.,1.)
!     **********************************************************

      IMPLICIT NONE
!     .. Scalar Arguments ..
      REAL x
      INTEGER j
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC aint
!     ..
      j = j + 1
      qranf = j*x
      qranf = qranf - aint(qranf)
      RETURN
   END FUNCTION
END
