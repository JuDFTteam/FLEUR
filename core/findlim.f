      MODULE m_findlim
c......................................................findlim
c finds turning point and practical "infinity"
c
      CONTAINS
      SUBROUTINE findlim(
     >                   mrad,lll,ec,vv,rc,
     <                   nmatch,nzero)

      IMPLICIT NONE
C     ..
C     .. Scalar Arguments ..
      REAL,    INTENT (IN) :: ec
      INTEGER, INTENT (IN) :: mrad,lll
      INTEGER, INTENT (OUT):: nmatch,nzero
C     ..
C     .. Array Arguments ..
      REAL   , INTENT (IN) :: rc(mrad),vv(mrad)
C     ..
C     .. Local Scalars ..
      REAL unend
      INTEGER n,nn
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC mod
C     ..
C     .. Data statements ..
      DATA unend/150.0/
C     ..
C                   --------------------
C--->                  FIND    NZERO
C                   --------------------
      DO 10 n = 1, (mrad-1)
         IF ((vv(n)-ec)*rc(n)**2.GT.unend) THEN
            IF (mod(n,2).EQ.0) THEN
               nzero = n + 1
            ELSE
               nzero = n
            END IF
            GO TO 20
         END IF
   10 CONTINUE
      nzero = mrad - 1
      WRITE (6,FMT=
     +'('' NRC='',I4,'' L='',I2,
     +  ''  NZERO SET TO  (NRC-1) ='',I4)') mrad,lll,(mrad-1)
   20 CONTINUE
C                     --------------------
C--->                   FIND    NMATCH
C                     --------------------
      n = nzero + 1
      DO nn = 1,nzero
         n = n - 1
!        IF ( (vv(n) + lll/rc(n)**2 - ec) < 0.0 ) THEN
         IF ((vv(n)-ec).LT.0.0) THEN
            nmatch = n
            RETURN
         END IF
      ENDDO
      WRITE (6,FMT=
     +'(//,''  STOP IN <<CORE>>'',/,
     +     '' NRC='',I2,'' L='',I2,/,
     +     ''  NO MATCHING-RADIUS FOUND FOR  EC='',F10.3)') mrad,lll,ec

      

      END SUBROUTINE findlim
      END MODULE m_findlim
