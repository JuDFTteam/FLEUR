      MODULE m_crtail
c
c initial point for income regular solution of dirac eq.
c order kap1=-L-1, kap2=L
c
      CONTAINS
      SUBROUTINE crtail(
     >                  mrad,e,rc,nsol,nzero,csq,
     X                  gc,fc)

      IMPLICIT NONE
C     ..
C     .. Scalar Arguments ..
      INTEGER, INTENT (IN) :: mrad,nsol,nzero
      REAL   , INTENT (IN) :: e,csq 
C     ..
C     .. Array Arguments ..
      REAL   , INTENT (IN) :: rc(mrad)
      REAL   , INTENT (INOUT) :: fc(2,2,mrad),gc(2,2,mrad)
C     ..
C     .. Local Scalars ..
      REAL beta,rr
      INTEGER i,ir,j
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC exp,sqrt
c
C----->infinity point solutions construction
c
      beta = sqrt(-e-e*e/csq)
      DO j = 1,nsol
         DO i = 1,nsol
            DO ir = nzero + 1,mrad
               rr = rc(ir) - rc(nzero)
               gc(i,j,ir) = gc(i,j,nzero)*exp(-beta*rr)
               fc(i,j,ir) = fc(i,j,nzero)*exp(-beta*rr)
            END DO
         END DO
      END DO

      END SUBROUTINE crtail
      END MODULE m_crtail
