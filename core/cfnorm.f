      SUBROUTINE cfnorm(mrad,is,it,nsol,nmatch,jtop,var,gc,fc,rc,rc2,
     +                  dx,gck,fck)
c...........................................................cfnorm
c wavefunctions normalization
c
      IMPLICIT NONE
c
C     .. Scalar Arguments ..
      INTEGER, INTENT (IN) :: mrad
      REAL dx
      INTEGER is,it,jtop,nmatch,nsol
C     ..
C     .. Array Arguments ..
      REAL fc(2,2,mrad),fck(2,2,mrad),gc(2,2,mrad),gck(2,2,mrad),
     +     rc(mrad),rc2(mrad),var(4)
C     ..
C     .. Local Scalars ..
      REAL xnorm
      INTEGER i,j,k,n
C     ..
C     .. Local Arrays ..
      REAL rint(mrad)
C     ..
C     .. External Functions ..
      REAL rsimp
      EXTERNAL rsimp
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC sqrt
C     ..
      DO n = 1,mrad
         DO k = 1,nsol
            gck(k,is,n) = 0.0
            fck(k,is,n) = 0.0
         END DO
      END DO
C                         ---------------------------------
C--->                     NORMALIZE WAVEFUNCTIONS ACCORDING
C                               TO MATCHING CONDITIONS
C                         ---------------------------------
C                                    INWARD - SOLUTION
      DO 30 n = nmatch,jtop
         DO 20 j = 1,nsol
            DO 10 i = 1,nsol
               gc(i,j,n) = gc(i,j,n)*var(2)
               fc(i,j,n) = fc(i,j,n)*var(2)
   10       CONTINUE
   20    CONTINUE
   30 CONTINUE
c
      IF (nsol.EQ.2) THEN
C                                   OUTWARD - SOLUTION
         DO 50 n = 1, (nmatch-1)
            DO 40 i = 1,nsol
               gc(i,it,n) = gc(i,it,n)*var(3)
               fc(i,it,n) = fc(i,it,n)*var(3)
   40       CONTINUE
   50    CONTINUE
C                                    INWARD - SOLUTION
         DO 70 n = nmatch,jtop
            DO 60 i = 1,nsol
               gc(i,it,n) = gc(i,it,n)*var(4)
               fc(i,it,n) = fc(i,it,n)*var(4)
   60       CONTINUE
   70    CONTINUE
      END IF
C                                    SUM FOR EACH KAPPA
      DO n = 1,jtop
         DO k = 1,nsol
            gck(k,is,n) = 0.00
            fck(k,is,n) = 0.00
         END DO
      END DO
      DO n = 1,jtop
         DO k = 1,nsol
            DO j = 1,nsol
               gck(k,is,n) = gck(k,is,n) + gc(k,j,n)
               fck(k,is,n) = fck(k,is,n) + fc(k,j,n)
            END DO
         END DO
      END DO
C                       -----------------------------------
C                       CALCULATE  NORM  AND NORMALIZE TO 1
C                       -----------------------------------
      DO n = 1,jtop
         rint(n) = 0.00
         DO k = 1,nsol
            rint(n) = rint(n) + rc2(n)* (gck(k,is,n)**2+fck(k,is,n)**2)
         END DO
      END DO
      xnorm = rsimp(mrad,rint,rc,jtop,dx)
      xnorm = 1.00/sqrt(xnorm)
      DO n = 1,jtop
         DO k = 1,nsol
            gck(k,is,n) = gck(k,is,n)*xnorm
            fck(k,is,n) = fck(k,is,n)*xnorm
         END DO
      END DO
      RETURN
      END
