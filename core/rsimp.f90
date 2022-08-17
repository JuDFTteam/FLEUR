MODULE m_rsimp

   CONTAINS

   REAL FUNCTION rsimp(mrad,f,r,irn,dx)

      IMPLICIT NONE

!     radial integration via simpson
      INTEGER,INTENT (IN) :: mrad
      REAL dx
      INTEGER irn

      REAL f(mrad),r(mrad)

      REAL s
      INTEGER i,isw,nl,np

      isw = 0
      rsimp = 0.0
      IF (irn.LE.2) RETURN
      IF (irn/2*2.EQ.irn) isw = 1
      np = irn - isw
      s = f(1)*r(1) + f(np)*r(np)
      nl = np - 1
      DO i = 2,nl,2
         s = s + 4.0*f(i)*r(i)
      END DO
      nl = nl - 1
      IF (.NOT.(nl.LT.3)) THEN
         DO i = 3,nl,2
            s = s + 2.0*f(i)*r(i)
         END DO
      END IF
      s = s*dx/3.0
      IF (.NOT.(isw.EQ.1)) THEN
         rsimp = s
         RETURN
      END IF
      rsimp = s + (f(irn)*r(irn)+f(irn-1)*r(irn-1))*0.50*dx
      RETURN
   END FUNCTION rsimp
END MODULE m_rsimp
