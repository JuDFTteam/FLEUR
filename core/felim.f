      MODULE m_felim
c.........................................................felim
c energy limits setup and checking
c
      CONTAINS
      SUBROUTINE felim(
     >                 mrad,lll,zz,nqn,vv,rc,
     <                 elim)
C     -----------------
C--->  FIND  E-LIMIT
C     -----------------

      IMPLICIT NONE
C     ..
C     .. Scalar Arguments ..
      INTEGER, INTENT (IN) :: mrad,lll,nqn
      REAL   , INTENT (IN) :: zz
      REAL   , INTENT (OUT):: elim
C     ..
C     .. Array Arguments ..
      REAL   , INTENT (IN) :: rc(mrad),vv(mrad)
C     ..
C     .. Local Scalars ..
      REAL val
      INTEGER n
C     ..
      IF (lll.EQ.0) THEN
         elim = -2.0*zz*zz/ (1.50*nqn*nqn)
      ELSE
         elim = vv(1) + lll/rc(1)**2
         DO 10 n = 2,mrad
            val = vv(n) + lll/rc(n)**2
            IF (val.LE.elim) elim = val
   10    CONTINUE
      END IF
      END SUBROUTINE felim
      END MODULE m_felim
