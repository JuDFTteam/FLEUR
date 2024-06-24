      MODULE m_radsrd
      use m_juDFT

c*********************************************************************
c     calculates the energy derivative of the scalar relativistic
c     wavefuction for energy e and angular momentum l.
c     the large (small) component is returned in pe (qe). for non-
c     relativistic case, set cin=0.
c     the functions p,q and the radial derivative dus are required.
c                 based on code by m. weinert  
c*********************************************************************

c*********************************************************************
c
c Uncomment the following, if you want to calculate an "exact" udot
c including relativistic corrections in the inhomogeneous part of the
c Dirac-DGL. Note that these corrections are inconsistent with the
c construction of the Hamiltonian in hssphn, where the inhomogeneous
c part is assumed to be non-relativistic.
c
c C. Friedrich June 2009
c*********************************************************************
# define RELATIVISTIC_CORRECTIONS

# ifdef RELATIVISTIC_CORRECTIONS
#   define CORR *(1.+cin2*fl1/rm**2)
# else
#   define CORR
# endif

      CONTAINS
      SUBROUTINE radsrd(
     >                  e,l,vr,r0,h,jri,c,
     <                  ud,dud,ddn,nodes,pe,qe,
     >                  p,q,dus)

      USE m_intgr, ONLY : intgr0
      IMPLICIT NONE
C     ..
C     .. Scalar Arguments ..
      INTEGER, INTENT (IN) :: l
      INTEGER, INTENT (IN) :: jri
      INTEGER, INTENT (OUT):: nodes
      REAL,    INTENT (IN) :: c
      REAL,    INTENT (IN) :: dus,e,h,r0
      REAL,    INTENT (OUT):: ddn,dud,ud
C     ..
C     .. Array Arguments ..
      REAL,    INTENT (IN) :: p(:),q(:),vr(:)
      REAL,    INTENT (OUT):: pe(:),qe(:)
C     ..
C     .. Local Scalars ..
      REAL dr,drh,erp,erq,fl1,p0,p1,p1p,q0,q1,q1p,r,rh,rh2,rm,rve,
     +     sk1,sk2,sk3,sk4,sl1,sl2,sl3,sl4,t,t1,t2,yn,zn,cin,cin2
      INTEGER i,it
C     ..
C     .. Local Arrays ..
      REAL pp(size(p)),qp(size(p))
C     ..
C     .. Data statements ..
      REAL,PARAMETER :: eps=1.e-06
C     ..
      cin = 1.0/c
      cin2 = cin*cin
c
      IF(jri>size(p))  CALL juDFT_error("dimension too small",
     +     calledby="radsrd")
c--->    set up initial conditions
      fl1 = l* (l+1)
      rm = 2.*r0 + cin2* (r0*e-vr(1))
      pe(1) = 0
      qe(1) = 0
      pp(1) =   r0*cin*q(1)
      qp(1) = - r0*p(1) CORR
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
         yn = pe(i) + 0.5*sk1
         zn = qe(i) + 0.5*sl1
         t1 = 0.5*rh2*cin* (q(i)+q(i+1))
         t2 = 0.5*rh2* (p(i)+p(i+1)) CORR
         sk2 = h* (rm*zn+yn+t1)
         sl2 = h* ((fl1/rm+rve)*yn-zn-t2)
         yn = pe(i) + 0.5*sk2
         zn = qe(i) + 0.5*sl2
         sk3 = h* (rm*zn+yn+t1)
         sl3 = h* ((fl1/rm+rve)*yn-zn-t2)
         rve = vr(i+1) - rh*e
         rm = 2.*rh - cin2*rve
         yn = pe(i) + sk3
         zn = qe(i) + sl3
         sk4 = h* (rm*zn+yn+rh*cin*q(i+1))
         sl4 = h* ((fl1/rm+rve)*yn-zn-rh*p(i+1) CORR )
         pe(i+1) = pe(i) + (sk1+2.*sk2+2.*sk3+sk4)/6.
         qe(i+1) = qe(i) + (sl1+2.*sl2+2.*sl3+sl4)/6.
         pp(i+1) = rm*qe(i+1) + pe(i+1) + rh*cin*q(i+1)
         qp(i+1) = (fl1/rm+rve)*pe(i+1) - qe(i+1) -
     +             rh*p(i+1) CORR
         r = rh
      ENDDO
      nodes = 0

c--->    adams-bashforth-moulton predictor-corrector
      predictor: DO i = 6,jri - 1
         r = r*dr
c--->    predictor
         p0 = pe(i) + h* (4277.*pp(i)-7923.*pp(i-1)+9982.*pp(i-2)-
     +        7298.*pp(i-3)+2877.*pp(i-4)-475.*pp(i-5))/1440.
         q0 = qe(i) + h* (4277.*qp(i)-7923.*qp(i-1)+9982.*qp(i-2)-
     +        7298.*qp(i-3)+2877.*qp(i-4)-475.*qp(i-5))/1440.
c--->    evaluate derivatives at next point
         rve = vr(i+1) - r*e
         rm = 2.*r - cin2*rve
         t1 = cin*r*q(i+1)
         t2 = r*p(i+1) CORR
         p1p = rm*q0 + p0 + t1
         q1p = (fl1/rm+rve)*p0 - q0 - t2
c--->    corrector
         corrector: DO it = 1,5
            p1 = pe(i) + h* (475.*p1p+1427.*pp(i)-798.*pp(i-1)+
     +           482.*pp(i-2)-173.*pp(i-3)+27.*pp(i-4))/1440.
            q1 = qe(i) + h* (475.*q1p+1427.*qp(i)-798.*qp(i-1)+
     +           482.*qp(i-2)-173.*qp(i-3)+27.*qp(i-4))/1440.
c--->    final evaluation
            p1p = rm*q1 + p1 + t1
            q1p = (fl1/rm+rve)*p1 - q1 - t2
c--->    test quality of corrector and iterate if necessary
            erp = abs(p1-p0)/ (abs(p1)+abs(h*p1p))
            erq = abs(q1-q0)/ (abs(q1)+abs(h*p1p))
            IF (erp.LT.eps .AND. erq.LT.eps) EXIT corrector
            p0 = p1
            q0 = q1
         ENDDO corrector
c--->    store values
         pe(i+1) = p1
         qe(i+1) = q1
         pp(i+1) = p1p
         qp(i+1) = q1p
         nodes = nodes + 0.501*abs(sign(1.0,pe(i+1))-sign(1.0,pe(i)))
      ENDDO predictor
c--->    ensure orthogonality
      DO i = 1,jri
         qe(i) = cin*qe(i)
      ENDDO
      DO i = 1,jri
         qp(i) = p(i)*pe(i) + q(i)*qe(i)
      ENDDO
      CALL intgr0(qp,r0,h,jri,t)
      dud = (pp(jri)-pe(jri))/ (r*r)
      DO i = 1,jri
         pe(i) = pe(i) - t*p(i)
         qe(i) = qe(i) - t*q(i)
      ENDDO
      ud = pe(jri)/r
      dud = dud - t*dus
      DO i = 1,jri
         qp(i) = pe(i)*pe(i) + qe(i)*qe(i)
      ENDDO
      CALL intgr0(qp,r0,h,jri,ddn)

      END SUBROUTINE radsrd
      END MODULE m_radsrd
