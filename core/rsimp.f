      REAL FUNCTION rsimp(mrad,f,r,irn,dx)
c.....................................................rsimp

      IMPLICIT NONE
c===================
c     radial integration via simpson
C     ..
C     .. Scalar Arguments ..
      INTEGER,INTENT (IN) :: mrad
      REAL dx
      INTEGER irn
C     ..
C     .. Array Arguments ..
      REAL f(mrad),r(mrad)
C     ..
C     .. Local Scalars ..
      REAL s
      INTEGER i,isw,nl,np
C     ..
      isw = 0
      rsimp = 0.0
      IF (irn.LE.2) RETURN
      IF (irn/2*2.EQ.irn) isw = 1
      np = irn - isw
      s = f(1)*r(1) + f(np)*r(np)
      nl = np - 1
      DO 10 i = 2,nl,2
         s = s + 4.0*f(i)*r(i)
   10 CONTINUE
      nl = nl - 1
      IF (nl.LT.3) GO TO 30
      DO 20 i = 3,nl,2
         s = s + 2.0*f(i)*r(i)
   20 CONTINUE
   30 s = s*dx/3.0
      IF (isw.EQ.1) GO TO 40
      rsimp = s
      RETURN
   40 rsimp = s + (f(irn)*r(irn)+f(irn-1)*r(irn-1))*0.50*dx
      RETURN
      END
