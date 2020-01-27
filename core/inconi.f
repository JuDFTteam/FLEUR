c.........................................................inconi
      SUBROUTINE inconi(l,xmj,e,csq,rr,xx1,xx2)
c initial point for income regular solution of dirac eq.
c order kap1=-L-1, kap2=L
c
C----->infinity point solutions construction
C     .. Scalar Arguments ..

      IMPLICIT NONE
      REAL csq,e,rr,xmj
      INTEGER l
C     ..
C     .. Array Arguments ..
      REAL xx1(4),xx2(4)
C     ..
C     .. Local Scalars ..
      REAL beta,bova
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC abs,exp,sqrt
C     ..
      beta = sqrt(-e-e*e/csq)
      bova = -beta/ (1.+e/csq)
      IF (abs(xmj).LT.l) THEN
         xx1(1) = exp(-beta*rr)
         xx1(2) = bova*xx1(1)
         xx1(3) = 0.
         xx1(4) = 0.
c
         xx2(1) = 0.
         xx2(2) = 0.
         xx2(3) = exp(-beta*rr)
         xx2(4) = bova*xx2(3)
      ELSE
c      XX1(1) =1.D-05
         xx1(1) = exp(-beta*rr)
         xx1(2) = bova*xx1(1)
c not needed
         xx1(3) = 0.0
         xx1(4) = 0.0
         xx2(1) = 0.0
         xx2(2) = 0.0
         xx2(3) = 0.0
         xx2(4) = 0.0
      END IF
      RETURN
      END
