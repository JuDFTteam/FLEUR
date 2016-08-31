      MODULE m_diflgr
      CONTAINS
      REAL FUNCTION diflgr(x,f)
c     **********************************************************
c     differentiate the function f,
c     given at the points x0,x1,x2,
c     at the point x1 by lagrange interpolation
c     e.wimmer
c     ***********************************************************
C     .. Array Arguments ..

      IMPLICIT NONE
      REAL f(0:2),x(0:2)
C     ..
      diflgr = (x(1)-x(2))*f(0)/ ((x(0)-x(1))* (x(0)-x(2))) +
     +         (2.0e0*x(1)-x(2)-x(0))*f(1)/ ((x(1)-x(0))* (x(1)-x(2))) +
     +         (x(1)-x(0))*f(2)/ ((x(2)-x(0))* (x(2)-x(1)))
c     print *,'x=',x
c     print *,'f=',f
c     print *,'diflgr=',diflgr
      RETURN
      END FUNCTION
      END
