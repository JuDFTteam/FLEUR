      MODULE m_diff
C   ******************************************************************  WAB07620
C   *                                                                *  WAB07630
C   *    <SUBROUTINE DIFF> DIFFERENTIATES THE POTENTIAL <V> IN THE   *  WAB07640
C   *    FIRST <N> NET POINTS.                                       *  WAB07650
C   *    <DX> IS THE INCREMENT IN THE LOGARITHMIC NET.               *  WAB07660
C   *    THE RESULT IS STORED IN <VM>.                               *  WAB07670
C   *                                                                *  WAB07680
C   ******************************************************************  WAB07690
      CONTAINS
      SUBROUTINE diff(
     >                mrad,v,dx,n,
     <                vm)

      IMPLICIT NONE

C     .. Scalar Arguments ..
      INTEGER, INTENT (IN) :: mrad,n
      REAL   , INTENT (IN) :: dx
C     ..
C     .. Array Arguments ..
      REAL   , INTENT (IN) ::   v(mrad)
      REAL   , INTENT (OUT)::  vm(mrad)
C     ..
C     .. Local Scalars ..
      INTEGER i,nm2
C     ..
      nm2 = n - 2
      vm(1) = v(1)
      vm(2) = ((6.0*v(3)+6.6666666667*v(5)+1.20*v(7))-
     +        (2.450*v(2)+7.50*v(4)+3.750*v(6)+1.0/6.0*v(8)))/dx

      DO i = 3,nm2
         vm(i) = ((v(i-2)+8.0*v(i+1))- (8.0*v(i-1)+v(i+2)))/12.0/dx
      ENDDO

      vm(n-1) = (v(n)-v(n-2))/2.0/dx
      vm(n) = (v(n-2)/2.0-2.0*v(n-1)+1.50*v(n))/dx

      END SUBROUTINE diff
      END MODULE m_diff
