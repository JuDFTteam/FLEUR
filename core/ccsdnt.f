      MODULE m_ccsdnt
c...........................................................ccdnt
c charge and spin density calculations
c
      CONTAINS
      SUBROUTINE ccsdnt(
     >                  mrad,is,jtop,nsol,
     >                  l,xmj,kap1,kap2,gck,fck,rc2,
     <                  rhochr,rhospn)
      IMPLICIT NONE
C     ..
C     .. Scalar Arguments ..
      INTEGER, INTENT (IN) :: mrad
      REAL xmj
      INTEGER is,jtop,kap1,kap2,l,nsol
C     ..
C     .. Array Arguments ..
      REAL fck(2,2,mrad),gck(2,2,mrad),rc2(mrad)
      REAL, INTENT (OUT) :: rhochr(mrad),rhospn(mrad)
C     ..
C     .. Local Scalars ..
      REAL cg1,cg2,cg4,cg5,cg8,cgo
      INTEGER ir,k,n
C     ..
C     .. Local Arrays ..
      REAL cgd(2),cgmd(2)
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC abs,sqrt
C     ..
      DO ir = 1,mrad
         rhochr(ir) = 0.0
         rhospn(ir) = 0.0
      END DO
c                       -----------------------------------
c                       coeffisients for spin-density
c                       -----------------------------------
c      KAP1 = - L - 1
c      KAP2 = + L
      cg1 = -xmj/ (kap1+0.50)
      cg5 = -xmj/ (-kap1+0.50)
      cgd(1) = cg1
      cgmd(1) = cg5
      IF (abs(xmj).GT.l) THEN
         cg2 = 0.00
         cg4 = 0.00
         cg8 = 0.00
         cgd(2) = 0.00
         cgo = 0.00
         cgmd(2) = 0.00
      ELSE
         cg2 = -sqrt(1.0- (xmj/ (kap1+0.50))**2)
         cg4 = -xmj/ (kap2+0.50)
         cg8 = -xmj/ (-kap2+0.50)
         cgd(2) = cg4
         cgo = cg2
         cgmd(2) = cg8
      END IF
C
      DO n = 1,jtop
         DO k = 1,nsol
            rhochr(n) = rhochr(n) + rc2(n)*
     +                  (gck(k,is,n)**2+fck(k,is,n)**2)
            rhospn(n) = rhospn(n) + rc2(n)*
     +                  (gck(k,is,n)*gck(k,is,n)*cgd(k)-
     +                  fck(k,is,n)*fck(k,is,n)*cgmd(k))
         END DO
      END DO
c
      IF (nsol.GT.1) THEN
         DO n = 1,jtop
            rhospn(n) = rhospn(n) + rc2(n)*
     +                  (gck(1,is,n)*gck(2,is,n)*cgo*2.)
         END DO
      END IF
c
      END SUBROUTINE ccsdnt
      END MODULE m_ccsdnt
