      MODULE m_radsra
      use m_juDFT
c*********************************************************************
c     calculates the scalar relativistic wavefuction for energy e and
c     angular momentum l for the potential vr by integrating outward.
c     the large (small) component is returned in p (q). for non-
c     relativistic case, set cin=0.        based on code by m. weinert 
c*********************************************************************
      CONTAINS
      SUBROUTINE radsra(
     >                  e,l,vr,r0,h,jri,c,
     <                  us,dus,nodes,p,q)

      USE m_intgr, ONLY : intgr0
      IMPLICIT NONE
C     ..
C     .. Scalar Arguments ..
      INTEGER, INTENT (IN) :: l
      INTEGER, INTENT (IN) :: jri
      INTEGER, INTENT (OUT):: nodes
      REAL,    INTENT (IN) :: e,h,r0
      REAL,    INTENT (IN) :: c
      REAL,    INTENT (OUT):: dus,us
C     ..
C     .. Array Arguments ..
      REAL,    INTENT (IN) :: vr(:)
      REAL,    INTENT (OUT):: p(:),q(:)
C     ..
C     .. Local Scalars ..
      REAL dr,drh,erp,erq,fl1,p0,p1,p1p,q0,q1,q1p,r,rh,rh2,rm,rve,s,
     +     sk1,sk2,sk3,sk4,sl1,sl2,sl3,sl4,t,yn,zn,cin,cin2
      INTEGER i,it
C     ..
C     .. Local Arrays ..
      REAL pp(jri),qp(jri)
C     ..
C     .. Data statements ..
      REAL,PARAMETER  :: eps=1.e-06
C     ..
      cin = 1.0/c
      cin2 = cin*cin

      IF (jri>SIZE(vr).OR.jri>SIZE(p))  CALL juDFT_error
     +    ("BUG: data dimension in radsra too small",calledby ="radsra")
   

c--->    set up initial conditions
      fl1 = l* (l+1)
      s = sqrt(fl1+1-cin2*vr(1)*vr(1)) - 1.
      rm = 2.*r0 + cin2* (r0*e-vr(1))
      p(1) = r0** (l+1)
      q(1) = s*p(1)/rm
      pp(1) = rm*q(1) + p(1)
      qp(1) = (fl1/rm+vr(1)-r0*e)*p(1) - q(1)
c--->    use 4th order runge-kutta to get first few mesh points
      dr = exp(h)
      drh = sqrt(dr)
      r = r0
      DO i = 1,5
         rh2 = drh*r
         rh = dr*r
         sk1 = h*pp(i)
         sl1 = h*qp(i)
         rve = 0.5* (vr(i)+vr(i+1)) - rh2*e
         rm = 2.*rh2 - cin2*rve
         yn = p(i) + 0.5*sk1
         zn = q(i) + 0.5*sl1
         sk2 = h* (rm*zn+yn)
         sl2 = h* ((fl1/rm+rve)*yn-zn)
         yn = p(i) + 0.5*sk2
         zn = q(i) + 0.5*sl2
         sk3 = h* (rm*zn+yn)
         sl3 = h* ((fl1/rm+rve)*yn-zn)
         rve = vr(i+1) - rh*e
         rm = 2.*rh - cin2*rve
         yn = p(i) + sk3
         zn = q(i) + sl3
         sk4 = h* (rm*zn+yn)
         sl4 = h* ((fl1/rm+rve)*yn-zn)
         p(i+1) = p(i) + (sk1+2.*sk2+2.*sk3+sk4)/6.
         q(i+1) = q(i) + (sl1+2.*sl2+2.*sl3+sl4)/6.
         pp(i+1) = rm*q(i+1) + p(i+1)
         qp(i+1) = (fl1/rm+rve)*p(i+1) - q(i+1)
         r = rh
      ENDDO
      nodes = 0
c--->    adams-bashforth-moulton predictor-corrector
      predictor: DO i = 6,jri - 1
         r = r*dr
c--->    predictor
         p0 = p(i) + h* (4277.*pp(i)-7923.*pp(i-1)+9982.*pp(i-2)-
     +        7298.*pp(i-3)+2877.*pp(i-4)-475.*pp(i-5))/1440.
         q0 = q(i) + h* (4277.*qp(i)-7923.*qp(i-1)+9982.*qp(i-2)-
     +        7298.*qp(i-3)+2877.*qp(i-4)-475.*qp(i-5))/1440.
c--->    evaluate derivatives at next point
         rve = vr(i+1) - r*e
         rm = 2.*r - cin2*rve
         p1p = rm*q0 + p0
         q1p = (fl1/rm+rve)*p0 - q0
c--->    corrector
         corrector: DO it = 1,5
            p1 = p(i) + h* (475.*p1p+1427.*pp(i)-798.*pp(i-1)+
     +           482.*pp(i-2)-173.*pp(i-3)+27.*pp(i-4))/1440.
            q1 = q(i) + h* (475.*q1p+1427.*qp(i)-798.*qp(i-1)+
     +           482.*qp(i-2)-173.*qp(i-3)+27.*qp(i-4))/1440.
c--->       final evaluation
            p1p = rm*q1 + p1
            q1p = (fl1/rm+rve)*p1 - q1
c--->       test quality of corrector and iterate if necessary
            erp = abs(p1-p0)/ (abs(p1)+abs(h*p1p))
            erq = abs(q1-q0)/ (abs(q1)+abs(h*p1p))
            IF (erp.LT.eps .AND. erq.LT.eps) EXIT corrector
            p0 = p1
            q0 = q1
         ENDDO corrector
c--->    store values
         p(i+1) = p1
         q(i+1) = q1
         pp(i+1) = p1p
         qp(i+1) = q1p
         nodes = nodes + 0.501*abs(sign(1.0,p(i+1))-sign(1.0,p(i)))
      ENDDO predictor

c--->    normalize function
      DO i = 1,jri
         q(i) = cin*q(i)
      ENDDO
      DO i = 1,jri
         qp(i) = p(i)*p(i) + q(i)*q(i)
      ENDDO
      CALL intgr0(qp,r0,h,jri,t)
      t = 1.0/sqrt(t)
      us = t*p(jri)/r
      dus = t* (pp(jri)-p(jri))/ (r*r)
      DO i = 1,jri
         p(i) = t*p(i)
         q(i) = t*q(i)
      ENDDO

      END SUBROUTINE radsra
      END MODULE m_radsra
