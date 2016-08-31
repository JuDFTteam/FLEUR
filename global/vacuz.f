      MODULE m_vacuz
      CONTAINS
      SUBROUTINE vacuz(
     >     e,vz,vz0,nmz,dz,
     <     uz,duz,u)
      use m_juDFT
c*********************************************************************
c     integrates the vacuum wavefunction fo energy e<0 inward from the
c     last mesh point using the schrodinger equation and assuming that
c     this point is past the classical turning point.
!     (vz(i),i = 1,nmz) contains the potential on a linear mesh.
c     vz0 is the vacuum zero.
c     dz  is the spacing of the mesh points.
c     uz,duz are the value and outward normal derivative (i.e. into
c     the bulk) evaluated at the first mesh point.
c     (u(i),i=1,nmz) contains the normalized wavefunction
c     based on code by m. weinert 
c*********************************************************************

      USE m_intgr, ONLY : intgz0
      IMPLICIT NONE
C     ..
C     .. Scalar Arguments ..
      INTEGER, INTENT (IN) :: nmz
      REAL,    INTENT (IN) :: dz,e,vz0
      REAL,    INTENT (OUT):: duz,uz
C     ..
C     .. Array Arguments ..
      REAL, INTENT (IN) :: vz(nmz)
      REAL, INTENT (OUT)::  u(nmz)
C     ..
C     .. Local Scalars ..
      REAL a,eps,eru,erw,h,sk1,sk2,sk3,sk4,sl1,sl2,sl3,sl4,tnorm,u1,u10,
     +     u1p,vme,w1,w10,w1p,yn,zn,znm
      INTEGER i,it,n,n1
      LOGICAL tail
C     ..
C     .. Local Arrays ..
      REAL w(nmz),wp(nmz)
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC abs,exp,sqrt
C     ..
C     .. Data statements ..
      DATA eps/1.e-06/
C     ..
      IF (e.GE.vz0) THEN
         WRITE (6,*) 'e>vz0; e=',e,'vz0=',vz0
         CALL juDFT_error("e>vz0",calledby ="vacuz",hint
     +        ="Vacuum energy parameter too high")
      ENDIF
c---  >    set up initial conditions
      znm = (nmz-1)*dz
      a = sqrt(2.* (vz0-e))
      u(nmz) = exp(-a*znm)
      w(nmz) = a*u(nmz)
      wp(nmz) = a*a*u(nmz)
c---  >    use 4nd order runge-kutta to first few mesh points
      h = dz
      DO i = 1,4
         n = nmz + 1 - i
         yn = u(n)
         zn = w(n)
         sk1 = h*zn
         sl1 = h*wp(n)
         sk2 = h* (zn+0.5*sl1)
         sl2 = h* ((vz(n)+vz(n-1)-e)* (yn+0.5*sk1))
         sk3 = h* (zn+0.5*sl2)
         sl3 = h* ((vz(n)+vz(n-1)-e)* (yn+0.5*sk2))
         sk4 = h* (zn+sl3)
         sl4 = h* (2.* (vz(n-1)-e)* (yn+sk3))
         u(n-1) = yn + (sk1+2.*sk2+2.*sk3+sk4)/6.
         w(n-1) = zn + (sl1+2.*sl2+2.*sl3+sl4)/6.
         wp(n-1) = 2.* (vz(n-1)-e)*u(n-1)
      ENDDO
c---  >    use adams-bashforth-moulton predictor-corrector for rest
      DO  i = 5,nmz - 1
c---  >    adams-bashforth predictor (order h**5)
         n = nmz + 1 - i
         n1 = n - 1
         u10 = u(n) + h* (55.*w(n)-59.*w(n+1)+37.*w(n+2)-9.*w(n+3))/24.
         w10 = w(n) + h* (55.*wp(n)-59.*wp(n+1)+37.*wp(n+2)-9.*wp(n+3))/
     +        24.
c---  >    evaluate derivative at next point( n-1 corresponds to m+1)
         vme = 2.* (vz(n-1)-e)
         u1p = w10
         w1p = vme*u10
c---  >    adams-moulton correctors (order h**6)
         DO it = 1,10
            u1 = u(n) + h* (251.*u1p+646.*w(n)-264.*w(n+1)+106.*w(n+2)-
     +           19.*w(n+3))/720.
            w1 = w(n) + h* (251.*w1p+646.*wp(n)-264.*wp(n+1)+
     +           106.*wp(n+2)-19.*wp(n+3))/720.
c---  >    final evaluation
            u1p = w1
            w1p = vme*u1
c---  >    test quality of corrector and iterate if necessary
            eru = abs(u1-u10)/ (abs(u1)+abs(h*u1p))
            erw = abs(w1-w10)/ (abs(w1)+abs(h*w1p))
            IF (eru.LT.eps .AND. erw.LT.eps) EXIT
            u10 = u1
            w10 = w1
         ENDDO
         IF (it>10) THEN
            WRITE (6,FMT=8000) n1,eru,erw
 8000       FORMAT (' ***vacuz - step may be too big - mesh point',i5,
     +           ', eru,erw=',1p,2e16.7)
         ENDIF
c---  >    store values
         u(n1) = u1
         w(n1) = w1
         wp(n1) = w1p
      ENDDO
c---  >    normalize function
      tail = .true.
      DO i = 1,nmz
         wp(i) = u(nmz+1-i)*u(nmz+1-i)
      ENDDO
      CALL intgz0(wp,dz,nmz,tnorm,tail)
      tnorm = 1./sqrt(tnorm)
      DO  i = 1,nmz
         u(i) = tnorm*u(i)
      ENDDO
      uz = u(1)
      duz = tnorm*w(1)
      RETURN
      END SUBROUTINE
      END 
