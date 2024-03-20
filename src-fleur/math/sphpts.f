      MODULE m_sphpts
      CONTAINS
      SUBROUTINE sphpts(p,n,r,pos)
c     *******************************************************
c     generates points on sphere at pos with radius r
c     e. wimmer     feb. 1980
c     modified to give a better distribution of points
c     m. weinert    jan. 1982
c     *******************************************************
      USE m_qranf
      USE m_constants, ONLY : tpi_const
      IMPLICIT NONE
C     .. Scalar Arguments ..
      REAL r
      INTEGER n
C     ..
C     .. Array Arguments ..
      REAL p(3,n),pos(3)
C     ..
C     .. Local Scalars ..
      REAL phi,t,tc,x,xr,y,yr,z
      INTEGER i,j
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC cos,sin,sqrt
C     ..
      j = 0
      xr = sqrt(13.e0)
      yr = sqrt(7.e0)
      DO  i = 1,n
         tc = 2.e0*qranf(xr,j) - 1.e0
         phi = tpi_const*qranf(yr,j)
         t = sqrt(1.e0-tc*tc)
         x = t*cos(phi)
         y = t*sin(phi)
         z = tc
         p(1,i) = r*x + pos(1)
         p(2,i) = r*y + pos(2)
         p(3,i) = r*z + pos(3)
      ENDDO   
      RETURN
      END SUBROUTINE
      END
