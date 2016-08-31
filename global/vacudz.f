      MODULE m_vacudz
      CONTAINS
      SUBROUTINE vacudz(
     >     e,vz,vz0,nmz,dz,
     <     udz,dudz,ddnv,ue,
     >     duz,u)
c*********************************************************************
c     integrates the vacuum energy derivative wavefunction ue for the
c     energy e>0 using the schrodinger equation. the wavefunction u
c     is necessary (obtained from a previous call to vacuz). the
c     conventions are the same as in vacuz. ue will be orthogonal to
c     u. udz and dudz are the value and normal derivative at the
c     vacuum boundary.
c     based on code by m. weinert
c*********************************************************************
      use m_juDFT
      USE m_intgr, ONLY : intgz0
      IMPLICIT NONE
C     ..
C     .. Scalar Arguments ..
      INTEGER, INTENT (IN) :: nmz
      REAL,    INTENT (IN) :: duz,dz,e,vz0
      REAL,    INTENT (OUT):: ddnv,dudz,udz
C     ..
C     .. Array Arguments ..
      REAL, INTENT (IN) ::  u(nmz),vz(nmz)
      REAL, INTENT (OUT):: ue(nmz)
C     ..
C     .. Local Scalars ..
      REAL a,eru,erw,h,sk1,sk2,sk3,sk4,sl1,sl2,sl3,sl4,tnorm,u1,u10,
     +     u1p,vme,w1,w10,w1p,yn,zn,znm
      INTEGER i,it,n,n1
      LOGICAL tail
C     ..
C     .. Local Arrays ..
      REAL we(nmz),wp(nmz)
C     ..
C     .. Data statements ..
      REAL,PARAMETER:: eps=1.e-06
C     ..
      IF (e.GE.vz0)  CALL
     +     juDFT_error("e >= vz0",calledby ="vacudz",hint =
     +     'Vacuum energy-parameter too high')
c---  >    set up initial conditions
      znm = (nmz-1)*dz
      a = sqrt(2.* (vz0-e))
      ue(nmz) = (a*znm-0.5)*u(nmz)/ (a*a)
      we(nmz) = (a*znm-1.5)*u(nmz)/a
      wp(nmz) = 2.* (vz(nmz)-e)*ue(nmz) - 2.*u(nmz)
c---  >    use 4nd order runge-kutta to first few mesh points
      h = dz
      DO  i = 1,4
         n = nmz + 1 - i
         yn = ue(n)
         zn = we(n)
         sk1 = h*zn
         sl1 = h*wp(n)
         sk2 = h* (zn+0.5*sl1)
         sl2 = h* ((vz(n)+vz(n-1)-e)* (yn+0.5*sk1)- (u(n)+u(n-1)))
         sk3 = h* (zn+0.5*sl2)
         sl3 = h* ((vz(n)+vz(n-1)-e)* (yn+0.5*sk2)- (u(n)+u(n-1)))
         sk4 = h* (zn+sl3)
         sl4 = h* (2.* (vz(n-1)-e)* (yn+sk3)-2.*u(n-1))
         ue(n-1) = yn + (sk1+2.*sk2+2.*sk3+sk4)/6.
         we(n-1) = zn + (sl1+2.*sl2+2.*sl3+sl4)/6.
         wp(n-1) = 2.* (vz(n-1)-e)*ue(n-1) - 2.*u(n-1)
      ENDDO
c---  >    use adams-bashforth-moulton predictor-corrector for rest
      DO i = 5,nmz - 1
c---  >    adams-bashforth predictor (order h**5)
         n = nmz + 1 - i
         n1 = n - 1
         u10 = ue(n) + h* (55.*we(n)-59.*we(n+1)+37.*we(n+2)-
     +        9.*we(n+3))/24.
         w10 = we(n) + h* (55.*wp(n)-59.*wp(n+1)+37.*wp(n+2)-
     +        9.*wp(n+3))/24.
c---  >    evaluate derivative at next point( n-1 corresponds to m+1)
         vme = 2.* (vz(n-1)-e)
         u1p = w10
         w1p = vme*u10 - 2.*u(n-1)
c---  >    adams-moulton correctors: (order h**6)
         DO  it = 1,10
            u1 = ue(n) + h* (251.*u1p+646.*we(n)-264.*we(n+1)+
     +           106.*we(n+2)-19.*we(n+3))/720.
            w1 = we(n) + h* (251.*w1p+646.*wp(n)-264.*wp(n+1)+
     +           106.*wp(n+2)-19.*wp(n+3))/720.
c---  >    final evaluation
            u1p = w1
            w1p = vme*u1 - 2.*u(n-1)
c---  >    test quality of corrector and iterate if necessary
            eru = abs(u1-u10)/ (abs(u1)+abs(h*u1p))
            erw = abs(w1-w10)/ (abs(w1)+abs(h*w1p))
            IF (eru.LT.eps .AND. erw.LT.eps) EXIT
            u10 = u1
            w10 = w1
         ENDDO
         IF (it>10) THEN
            WRITE (6,FMT=8000) n1,eru,erw
 8000       FORMAT (' ***vacudz - step may be too big - mesh point',i5,
     +           ', eru,erw=',1p,2e16.7)
         ENDIF
c---  >    store values
         ue(n1) = u1
         we(n1) = w1
         wp(n1) = w1p
      ENDDO
c---  >    ensure orthogonality
      tail = .true.
      DO  i = 1,nmz
         wp(i) = ue(nmz+1-i)*u(nmz+1-i)
      ENDDO
      CALL intgz0(wp,dz,nmz,tnorm,tail)
      ue(:) = ue(:) - tnorm*u(:)
      udz = ue(1)
      dudz = we(1) - tnorm*duz
c---  >    determine normalization
      DO  i = 1,nmz
         wp(i) = ue(nmz+1-i)*ue(nmz+1-i)
      enddo
      CALL intgz0(wp,dz,nmz,ddnv,tail)
      RETURN
      END subroutine
      END
