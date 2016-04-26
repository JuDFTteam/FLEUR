c................................................................coreerr
      SUBROUTINE coreerr(err,var,s,nsol,pow,qow,piw,qiw)
C   ********************************************************************
C   *                                                                  *
C   *   CALCULATE THE MISMATCH OF THE RADIAL WAVE FUNCTIONS AT THE     *
C   *   POINT  NMATCH  FOR OUT- AND INWARD INTEGRATION                 *
C   *                                                                  *
C   ********************************************************************
C
c
C     .. Scalar Arguments ..

      IMPLICIT NONE
      INTEGER nsol,s
C     ..
C     .. Array Arguments ..
      REAL err(4),piw(2,2),pow(2,2),qiw(2,2),qow(2,2),var(4)
C     ..
C     .. Local Scalars ..
      INTEGER t
C     ..
      err(1) = pow(s,s) - piw(s,s)*var(2)
      err(2) = qow(s,s) - qiw(s,s)*var(2)
c
      IF (nsol.EQ.1) RETURN
c
      t = 3 - s
c
      err(1) = err(1) + pow(s,t)*var(3) - piw(s,t)*var(2)*var(4)
      err(2) = err(2) + qow(s,t)*var(3) - qiw(s,t)*var(2)*var(4)
      err(3) = pow(t,s) - piw(t,s)*var(2) + pow(t,t)*var(3) -
     +         piw(t,t)*var(2)*var(4)
      err(4) = qow(t,s) - qiw(t,s)*var(2) + qow(t,t)*var(3) -
     +         qiw(t,t)*var(2)*var(4)
c
      RETURN
      END
